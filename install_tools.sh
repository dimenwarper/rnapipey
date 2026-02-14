#!/usr/bin/env bash
set -euo pipefail

# ── rnapipey tool installer ──────────────────────────────────────────────────
# Installs external bioinformatics tools used by the rnapipey pipeline.
# Usage: ./install_tools.sh --all    (install everything)
#        ./install_tools.sh --infernal --viennarna   (selective)
#        ./install_tools.sh --help
# ─────────────────────────────────────────────────────────────────────────────

RNAPIPEY_HOME="${RNAPIPEY_HOME:-$HOME/.rnapipey}"
DATA_DIR="$RNAPIPEY_HOME/data"
TOOLS_DIR="$RNAPIPEY_HOME/tools"
CONDA_ENV="rnapipey"
PYTHON_VER="3.10"
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

ensure_conda() {
    if ! command -v conda &>/dev/null; then
        fail "conda not found. Install Miniconda/Mambaforge first."
        exit 1
    fi
}

ensure_env() {
    if ! conda env list | grep -qw "$CONDA_ENV"; then
        info "Creating conda environment '$CONDA_ENV' (Python $PYTHON_VER)..."
        conda create -y -n "$CONDA_ENV" "python>=$PYTHON_VER" pip
        ok "Environment '$CONDA_ENV' created."
    else
        skip "Conda environment '$CONDA_ENV' already exists."
    fi
}

# Run a command inside the rnapipey conda env.
run_in_env() {
    conda run --no-banner -n "$CONDA_ENV" "$@"
}

mkdir_p() {
    mkdir -p "$1"
}

# ── Tool installers ─────────────────────────────────────────────────────────

install_infernal() {
    header "Infernal"
    if run_in_env command -v cmscan &>/dev/null; then
        skip "Infernal already installed in '$CONDA_ENV'."
        return 0
    fi
    info "Installing Infernal..."
    conda install -y -n "$CONDA_ENV" -c bioconda infernal
    ok "Infernal installed."
}

install_viennarna() {
    header "ViennaRNA"
    if run_in_env command -v RNAfold &>/dev/null; then
        skip "ViennaRNA already installed in '$CONDA_ENV'."
        return 0
    fi
    info "Installing ViennaRNA..."
    conda install -y -n "$CONDA_ENV" -c bioconda viennarna
    ok "ViennaRNA installed."
}

install_pymol() {
    header "PyMOL"
    if run_in_env command -v pymol &>/dev/null; then
        skip "PyMOL already installed in '$CONDA_ENV'."
        return 0
    fi
    info "Installing PyMOL (open-source)..."
    conda install -y -n "$CONDA_ENV" -c conda-forge pymol-open-source
    ok "PyMOL installed."
}

install_protenix() {
    header "Protenix"
    if run_in_env pip show protenix &>/dev/null; then
        skip "Protenix already installed in '$CONDA_ENV'."
        return 0
    fi
    info "Installing Protenix..."
    run_in_env pip install protenix
    ok "Protenix installed."
}

install_rnadvisor() {
    header "RNAdvisor"
    if run_in_env pip show rnadvisor &>/dev/null; then
        skip "RNAdvisor already installed in '$CONDA_ENV'."
        return 0
    fi
    info "Installing RNAdvisor..."
    run_in_env pip install rnadvisor
    ok "RNAdvisor installed."
}

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
        wget -q --show-progress -O "$rfam_dir/Rfam.cm.gz" \
            "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"
        info "Extracting Rfam.cm.gz..."
        gunzip -f "$rfam_dir/Rfam.cm.gz"

        info "Downloading Rfam.clanin..."
        wget -q --show-progress -O "$rfam_clanin" \
            "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.clanin"

        ok "Rfam database downloaded."
    fi

    # Press the CM database for fast searching
    if [[ -f "$rfam_cm.i1m" ]]; then
        skip "Rfam.cm already pressed."
    else
        info "Running cmpress on Rfam.cm (this may take a few minutes)..."
        run_in_env cmpress "$rfam_cm"
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

    # Install into conda env
    info "Installing RhoFold+ into '$CONDA_ENV' env..."
    run_in_env pip install -e "$rhofold_dir" 2>/dev/null || {
        # Fallback: some versions use setup.py
        (cd "$rhofold_dir" && run_in_env python setup.py install) 2>/dev/null || {
            fail "Could not install RhoFold+ automatically. You may need to install it manually."
            info "  cd $rhofold_dir && pip install -e ."
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
            info "  conda install -c conda-forge git-lfs"
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
            # Fallback: try wget
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

install_simrna() {
    header "SimRNA"
    echo ""
    info "${YELLOW}SimRNA requires a manual download.${NC}"
    info "It is distributed under an academic license from genesilico.pl."
    info ""
    info "Steps:"
    info "  1. Visit https://genesilico.pl/SimRNAweb/ and request a download"
    info "  2. Extract to $TOOLS_DIR/SimRNA/"
    info "  3. Update configs/local.yaml with the binary and data_dir paths"
    echo ""
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
    device: "cuda:0"

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

${BOLD}USAGE${NC}
    ./install_tools.sh [OPTIONS]

${BOLD}OPTIONS${NC}
    --all           Install everything (recommended)
    --infernal      Install Infernal (conda)
    --viennarna     Install ViennaRNA (conda)
    --pymol         Install PyMOL (conda)
    --protenix      Install Protenix (pip)
    --rnadvisor     Install RNAdvisor (pip)
    --rfam          Download Rfam database
    --rhofold       Clone & install RhoFold+, download weights
    --spotrna       Clone SPOT-RNA & download models
    --simrna        Print SimRNA manual install instructions
    --config-only   Only (re)generate configs/local.yaml
    --help          Show this help

${BOLD}EXAMPLES${NC}
    ./install_tools.sh --all
    ./install_tools.sh --infernal --viennarna --rfam
    ./install_tools.sh --config-only

${BOLD}ENVIRONMENT${NC}
    RNAPIPEY_HOME   Base directory for data/tools (default: \$HOME/.rnapipey)

Tools are installed into the conda environment '${CONDA_ENV}'.
Data and cloned repos go under \$RNAPIPEY_HOME.
EOF
}

main() {
    if [[ $# -eq 0 ]]; then
        usage
        exit 0
    fi

    local do_all=false
    local do_infernal=false
    local do_viennarna=false
    local do_pymol=false
    local do_protenix=false
    local do_rnadvisor=false
    local do_rfam=false
    local do_rhofold=false
    local do_spotrna=false
    local do_simrna=false
    local do_config_only=false

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --all)        do_all=true ;;
            --infernal)   do_infernal=true ;;
            --viennarna)  do_viennarna=true ;;
            --pymol)      do_pymol=true ;;
            --protenix)   do_protenix=true ;;
            --rnadvisor)  do_rnadvisor=true ;;
            --rfam)       do_rfam=true ;;
            --rhofold)    do_rhofold=true ;;
            --spotrna)    do_spotrna=true ;;
            --simrna)     do_simrna=true ;;
            --config-only) do_config_only=true ;;
            --help|-h)    usage; exit 0 ;;
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

    ensure_conda
    ensure_env

    if $do_all || $do_infernal;  then install_infernal;  fi
    if $do_all || $do_viennarna; then install_viennarna; fi
    if $do_all || $do_pymol;     then install_pymol;     fi
    if $do_all || $do_protenix;  then install_protenix;  fi
    if $do_all || $do_rnadvisor; then install_rnadvisor; fi
    if $do_all || $do_rfam;      then install_rfam;      fi
    if $do_all || $do_rhofold;   then install_rhofold;   fi
    if $do_all || $do_spotrna;   then install_spotrna;   fi
    if $do_all || $do_simrna;    then install_simrna;    fi

    generate_config

    header "Done"
    echo ""
    ok "Installation complete."
    info "Activate the environment:  ${BOLD}conda activate $CONDA_ENV${NC}"
    info "Check tool availability:   ${BOLD}rnapipey check -c configs/local.yaml${NC}"
    echo ""
}

main "$@"
