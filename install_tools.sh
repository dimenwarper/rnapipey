#!/usr/bin/env bash
set -euo pipefail

# ── rnapipey tool installer ──────────────────────────────────────────────────
# Downloads data and clones git-based tools that can't be installed via pixi.
# Conda/pip packages (Infernal, ViennaRNA, PyMOL, Protenix, RNAdvisor) are
# managed by pixi — just run `pixi install`.
#
# Usage: ./install_tools.sh --all
#        ./install_tools.sh --rfam --rhofold
#        ./install_tools.sh --help
# ─────────────────────────────────────────────────────────────────────────────

RNAPIPEY_HOME="${RNAPIPEY_HOME:-$HOME/.rnapipey}"
DATA_DIR="$RNAPIPEY_HOME/data"
TOOLS_DIR="$RNAPIPEY_HOME/tools"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Colors ───────────────────────────────────────────────────────────────────

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

ok()   { echo -e "${GREEN}[OK]${NC}   $*"; }
skip() { echo -e "${YELLOW}[SKIP]${NC} $*"; }
fail() { echo -e "${RED}[FAIL]${NC} $*"; }
info() { echo -e "${CYAN}[INFO]${NC} $*"; }
header() { echo -e "\n${BOLD}── $* ──${NC}"; }

# ── Helpers ──────────────────────────────────────────────────────────────────

mkdir_p() {
    mkdir -p "$1"
}

# ── Data & tool installers ──────────────────────────────────────────────────

install_rfam() {
    header "Rfam database"
    local rfam_dir="$DATA_DIR/rfam"
    local rfam_cm="$rfam_dir/Rfam.cm"
    local rfam_clanin="$rfam_dir/Rfam.clanin"

    mkdir_p "$rfam_dir"

    if [[ -f "$rfam_cm" && -f "$rfam_clanin" ]]; then
        skip "Rfam database already downloaded."
    else
        info "Downloading Rfam.cm.gz from EBI FTP..."
        curl -fSL --progress-bar -o "$rfam_dir/Rfam.cm.gz" \
            "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
        info "Extracting Rfam.cm.gz..."
        gunzip -f "$rfam_dir/Rfam.cm.gz"

        info "Downloading Rfam.clanin..."
        curl -fSL --progress-bar -o "$rfam_clanin" \
            "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin"

        ok "Rfam database downloaded."
    fi

    # Press the CM database for fast searching
    if [[ -f "$rfam_cm.i1m" ]]; then
        skip "Rfam.cm already pressed."
    else
        info "Running cmpress on Rfam.cm (this may take a few minutes)..."
        pixi run cmpress "$rfam_cm"
        ok "Rfam.cm pressed."
    fi
}

install_rhofold() {
    header "RhoFold+"
    local rhofold_dir="$TOOLS_DIR/RhoFold"
    local weights_dir="$DATA_DIR/rhofold_weights"

    mkdir_p "$TOOLS_DIR"
    mkdir_p "$DATA_DIR"

    # Clone repo
    if [[ -d "$rhofold_dir" ]]; then
        skip "RhoFold+ repo already cloned."
    else
        info "Cloning RhoFold+..."
        git clone https://github.com/ml4bio/RhoFold.git "$rhofold_dir"
        ok "RhoFold+ cloned."
    fi

    # Install into pixi env
    info "Installing RhoFold+ into pixi environment..."
    pixi run pip install -e "$rhofold_dir" 2>/dev/null || {
        (cd "$rhofold_dir" && pixi run python setup.py install) 2>/dev/null || {
            fail "Could not install RhoFold+ automatically. You may need to install it manually."
            info "  cd $rhofold_dir && pixi run pip install -e ."
        }
    }

    # Download pretrained weights
    if [[ -d "$weights_dir" && -n "$(ls -A "$weights_dir" 2>/dev/null)" ]]; then
        skip "RhoFold+ weights already downloaded."
    else
        info "Downloading RhoFold+ pretrained weights from HuggingFace..."
        mkdir_p "$weights_dir"
        if command -v git-lfs &>/dev/null || git lfs version &>/dev/null 2>&1; then
            git clone https://huggingface.co/yangjianfeng/RhoFold "$weights_dir/hf_repo" 2>/dev/null && \
                mv "$weights_dir/hf_repo"/* "$weights_dir/" 2>/dev/null && \
                rm -rf "$weights_dir/hf_repo" && \
                ok "RhoFold+ weights downloaded." || {
                    fail "Could not download RhoFold+ weights. Install git-lfs and retry."
                    info "  git lfs install && git clone https://huggingface.co/yangjianfeng/RhoFold $weights_dir"
                }
        else
            fail "git-lfs not found. Install it to download RhoFold+ weights."
            info "  pixi add git-lfs"
            info "  git lfs install"
            info "  git clone https://huggingface.co/yangjianfeng/RhoFold $weights_dir"
        fi
    fi
}

install_spotrna() {
    header "SPOT-RNA"
    local spotrna_dir="$TOOLS_DIR/SPOT-RNA"

    mkdir_p "$TOOLS_DIR"

    if [[ -d "$spotrna_dir" ]]; then
        skip "SPOT-RNA repo already cloned."
    else
        info "Cloning SPOT-RNA..."
        git clone https://github.com/jaswindersingh2/SPOT-RNA.git "$spotrna_dir"
        ok "SPOT-RNA cloned."
    fi

    # Download models
    if [[ -d "$spotrna_dir/SPOT-RNA-models" ]]; then
        skip "SPOT-RNA models already present."
    else
        info "Downloading SPOT-RNA models..."
        (cd "$spotrna_dir" && bash download_models.sh 2>/dev/null) || {
            mkdir -p "$spotrna_dir/SPOT-RNA-models"
            fail "Could not auto-download SPOT-RNA models. Run manually:"
            info "  cd $spotrna_dir && bash download_models.sh"
        }
    fi

    echo ""
    info "${YELLOW}NOTE:${NC} SPOT-RNA requires TensorFlow 1.14 / Python 3.6."
    info "It needs a separate conda environment to run:"
    info "  conda create -n spotrna python=3.6"
    info "  conda activate spotrna"
    info "  cd $spotrna_dir && pip install -r requirements.txt"
}

install_protenix() {
    header "Protenix"
    local protenix_data_dir="$DATA_DIR/protenix"

    if pixi run pip show protenix 2>/dev/null | grep -q "Location" && \
       pixi run python -c "from protenix.predictor import predict" 2>/dev/null; then
        skip "Protenix already installed (real package)."
    else
        info "Installing Protenix from GitHub..."
        pixi run pip install protenix@git+https://github.com/bytedance/Protenix.git
        if pixi run protenix --help &>/dev/null; then
            ok "Protenix installed."
        else
            fail "Protenix install may have failed. Try manually:"
            info "  pixi run pip install protenix@git+https://github.com/bytedance/Protenix.git"
        fi
    fi

    # Download CCD data (components.cif) required for inference
    if [[ -f "$protenix_data_dir/common/components.cif" ]]; then
        skip "Protenix CCD data already downloaded."
    else
        info "Downloading Protenix CCD data for inference..."
        mkdir_p "$protenix_data_dir"
        local download_script
        download_script="$(pixi run python -c "import protenix; from pathlib import Path; print(Path(protenix.__file__).parent.parent)" 2>/dev/null)"
        if [[ -f "$download_script/scripts/database/download_protenix_data.sh" ]]; then
            PROTENIX_ROOT_DIR="$protenix_data_dir" bash "$download_script/scripts/database/download_protenix_data.sh" --inference_only
        else
            # Fallback: download script from GitHub
            PROTENIX_ROOT_DIR="$protenix_data_dir" bash <(curl -fsSL \
                "https://raw.githubusercontent.com/bytedance/Protenix/main/scripts/database/download_protenix_data.sh") \
                --inference_only
        fi

        if [[ -f "$protenix_data_dir/common/components.cif" ]]; then
            ok "Protenix CCD data downloaded."
        else
            fail "CCD data download may have failed. Try manually:"
            info "  export PROTENIX_ROOT_DIR=$protenix_data_dir"
            info "  curl -sL https://raw.githubusercontent.com/bytedance/Protenix/main/scripts/database/download_protenix_data.sh | bash -s -- --inference_only"
        fi
    fi
}

install_simrna() {
    header "SimRNA"
    local simrna_dir="$TOOLS_DIR/SimRNA"

    mkdir_p "$TOOLS_DIR"

    if [[ -f "$simrna_dir/SimRNA" ]]; then
        skip "SimRNA already installed."
    else
        local simrna_url
        case "$(uname -s)-$(uname -m)" in
            Darwin-*)
                simrna_url="https://genesilico.pl/software/simrna/version_3.20/SimRNA_64bitIntel_MacOSX_staticLibs.tgz"
                ;;
            Linux-x86_64)
                simrna_url="https://genesilico.pl/software/simrna/version_3.20/SimRNA_64bitIntel_Linux.tgz"
                ;;
            *)
                simrna_url="https://genesilico.pl/software/simrna/version_3.20/SimRNA_32bitIntel_Linux.tgz"
                ;;
        esac

        info "Downloading SimRNA from genesilico.pl ($(uname -s) $(uname -m))..."
        local tmp_tgz="$TOOLS_DIR/SimRNA.tgz"
        curl -fSL --progress-bar -o "$tmp_tgz" "$simrna_url"

        info "Extracting SimRNA..."
        mkdir_p "$simrna_dir"
        tar xzf "$tmp_tgz" -C "$simrna_dir" --strip-components=1
        rm -f "$tmp_tgz"

        if [[ -f "$simrna_dir/SimRNA" ]]; then
            chmod +x "$simrna_dir/SimRNA"
            ok "SimRNA installed."
        else
            fail "SimRNA binary not found after extraction. Check archive structure."
        fi
    fi
}

install_docker() {
    header "Docker"

    if command -v docker &>/dev/null; then
        skip "Docker CLI already installed."
        # Verify daemon is running
        if docker info &>/dev/null 2>&1; then
            ok "Docker daemon is running."
        else
            info "Docker CLI found but daemon not running."
            info "  Linux:  sudo systemctl start docker"
            info "  macOS:  colima start"
        fi
        return 0
    fi

    case "$(uname -s)" in
        Linux)
            info "Installing Docker on Linux..."
            if command -v apt-get &>/dev/null; then
                info "Using apt (Debian/Ubuntu)..."
                sudo apt-get update -qq
                sudo apt-get install -y -qq docker.io docker-compose-plugin
                sudo systemctl enable --now docker
                sudo usermod -aG docker "$USER" 2>/dev/null || true
                ok "Docker installed. You may need to log out/in for group changes."
            elif command -v dnf &>/dev/null; then
                info "Using dnf (Fedora/RHEL)..."
                sudo dnf install -y -q docker docker-compose-plugin
                sudo systemctl enable --now docker
                sudo usermod -aG docker "$USER" 2>/dev/null || true
                ok "Docker installed. You may need to log out/in for group changes."
            elif command -v yum &>/dev/null; then
                info "Using yum (CentOS/RHEL)..."
                sudo yum install -y -q docker docker-compose-plugin
                sudo systemctl enable --now docker
                sudo usermod -aG docker "$USER" 2>/dev/null || true
                ok "Docker installed. You may need to log out/in for group changes."
            else
                fail "No supported package manager found. Install Docker manually:"
                info "  https://docs.docker.com/engine/install/"
                return 1
            fi
            ;;
        Darwin)
            if command -v brew &>/dev/null; then
                info "Installing Docker CLI and Colima via Homebrew..."
                brew install docker docker-compose colima
                ok "Docker CLI + Colima installed."

                if ! colima status &>/dev/null 2>&1; then
                    info "Starting Colima..."
                    colima start --memory 4 --cpu 2
                    ok "Colima started."
                fi
            else
                fail "Homebrew not found. Install Docker manually:"
                info "  https://docs.docker.com/get-docker/"
                return 1
            fi
            ;;
        *)
            fail "Unsupported OS: $(uname -s). Install Docker manually:"
            info "  https://docs.docker.com/get-docker/"
            return 1
            ;;
    esac
}

# ── Config generation ────────────────────────────────────────────────────────

generate_config() {
    header "Generating configs/local.yaml"

    local config_file="$SCRIPT_DIR/configs/local.yaml"
    mkdir_p "$SCRIPT_DIR/configs"

    # Discover paths
    local cmscan_path="cmscan"
    local cmfetch_path="cmfetch"
    local cmalign_path="cmalign"
    local rfam_cm_path=""
    local rfam_clanin_path=""
    local rnafold_path="RNAfold"
    local spotrna_path=""
    local rhofold_script=""
    local rhofold_model_dir=""
    local simrna_binary=""
    local simrna_data=""

    [[ -f "$DATA_DIR/rfam/Rfam.cm" ]] && rfam_cm_path="$DATA_DIR/rfam/Rfam.cm"
    [[ -f "$DATA_DIR/rfam/Rfam.clanin" ]] && rfam_clanin_path="$DATA_DIR/rfam/Rfam.clanin"
    [[ -f "$TOOLS_DIR/SPOT-RNA/SPOT-RNA.py" ]] && spotrna_path="$TOOLS_DIR/SPOT-RNA/SPOT-RNA.py"

    # RhoFold+ inference script
    if [[ -d "$TOOLS_DIR/RhoFold" ]]; then
        local candidate
        for candidate in "$TOOLS_DIR/RhoFold/inference.py" "$TOOLS_DIR/RhoFold/rhofold/inference.py"; do
            if [[ -f "$candidate" ]]; then
                rhofold_script="$candidate"
                break
            fi
        done
    fi
    [[ -d "$DATA_DIR/rhofold_weights" && -n "$(ls -A "$DATA_DIR/rhofold_weights" 2>/dev/null)" ]] && \
        rhofold_model_dir="$DATA_DIR/rhofold_weights"

    # SimRNA (if manually installed)
    [[ -f "$TOOLS_DIR/SimRNA/SimRNA" ]] && simrna_binary="$TOOLS_DIR/SimRNA/SimRNA"
    [[ -d "$TOOLS_DIR/SimRNA/data" ]] && simrna_data="$TOOLS_DIR/SimRNA/data"

    # Protenix CCD data
    local protenix_data_dir=""
    [[ -f "$DATA_DIR/protenix/common/components.cif" ]] && protenix_data_dir="$DATA_DIR/protenix"

    cat > "$config_file" <<YAML
# rnapipey local configuration — generated by install_tools.sh
# Paths point to tools installed in $RNAPIPEY_HOME

tools:
  # Infernal / Rfam
  cmscan: ${cmscan_path}
  cmfetch: ${cmfetch_path}
  cmalign: ${cmalign_path}
  rfam_cm: "${rfam_cm_path}"
  rfam_clanin: "${rfam_clanin_path}"

  # ViennaRNA
  rnafold: ${rnafold_path}

  # SPOT-RNA (optional — needs separate Python 3.6 env)
  spotrna: "${spotrna_path}"

  # RhoFold+
  rhofold:
    script: "${rhofold_script}"
    model_dir: "${rhofold_model_dir}"
    device: "cpu"

  # SimRNA (manual install required — academic license)
  simrna:
    binary: "${simrna_binary}"
    data_dir: "${simrna_data}"
    replicas: 10
    steps: 10000000
    clustering_top_n: 5

  # Protenix
  protenix:
    binary: protenix
    model: ""
    data_dir: "${protenix_data_dir}"

  # RNAdvisor
  rnadvisor:
    docker: false
    metrics:
      - rsRNASP
      - DFIRE
      - RASP
      - MCQ
YAML

    ok "Generated $config_file"
}

# ── CLI ──────────────────────────────────────────────────────────────────────

usage() {
    cat <<EOF
${BOLD}rnapipey tool installer${NC}

Conda/pip packages (Infernal, ViennaRNA, PyMOL, Protenix, RNAdvisor) are
managed by pixi. Run ${BOLD}pixi install${NC} first.

This script handles data downloads and git-based tools that pixi can't manage.

${BOLD}USAGE${NC}
    ./install_tools.sh [OPTIONS]

${BOLD}OPTIONS${NC}
    --all           Download/clone everything below
    --rfam          Download Rfam database + cmpress
    --rhofold       Clone RhoFold+, install into pixi env, download weights
    --spotrna       Clone SPOT-RNA & download models
    --simrna        Download SimRNA binary
    --protenix      Install Protenix from GitHub
    --docker        Install Docker CLI + Colima (macOS) for RNAdvisor
    --config-only   Only (re)generate configs/local.yaml
    --help          Show this help

${BOLD}EXAMPLES${NC}
    pixi install && ./install_tools.sh --all
    ./install_tools.sh --rfam --rhofold
    ./install_tools.sh --config-only

${BOLD}ENVIRONMENT${NC}
    RNAPIPEY_HOME   Base directory for data/tools (default: \$HOME/.rnapipey)
EOF
}

main() {
    if [[ $# -eq 0 ]]; then
        usage
        exit 0
    fi

    local do_all=false
    local do_rfam=false
    local do_rhofold=false
    local do_spotrna=false
    local do_simrna=false
    local do_protenix=false
    local do_docker=false
    local do_config_only=false

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --all)         do_all=true ;;
            --rfam)        do_rfam=true ;;
            --rhofold)     do_rhofold=true ;;
            --spotrna)     do_spotrna=true ;;
            --simrna)      do_simrna=true ;;
            --protenix)    do_protenix=true ;;
            --docker)      do_docker=true ;;
            --config-only) do_config_only=true ;;
            --help|-h)     usage; exit 0 ;;
            *)
                fail "Unknown option: $1"
                echo "Run ./install_tools.sh --help for usage."
                exit 1
                ;;
        esac
        shift
    done

    echo -e "${BOLD}rnapipey tool installer${NC}"
    echo -e "Data/tools directory: ${CYAN}$RNAPIPEY_HOME${NC}"
    echo ""

    if $do_config_only; then
        generate_config
        exit 0
    fi

    if $do_all || $do_rfam;    then install_rfam;    fi
    if $do_all || $do_rhofold; then install_rhofold; fi
    if $do_all || $do_spotrna; then install_spotrna; fi
    if $do_all || $do_simrna;    then install_simrna;    fi
    if $do_all || $do_protenix; then install_protenix; fi
    if $do_all || $do_docker;   then install_docker;   fi

    generate_config

    header "Done"
    echo ""
    ok "Installation complete."
    info "Run your pipeline:  ${BOLD}pixi run rnapipey check -c configs/local.yaml${NC}"
    echo ""
}

main "$@"
