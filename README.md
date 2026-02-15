# rnapipey

RNA 3D structure prediction pipeline. Wraps multiple prediction methods behind a unified CLI, runs them on the same input, scores the results, and picks the best model.

## Pipeline

```
Sequence (FASTA) → Rfam/Infernal (MSA) → RNAfold (2D) → 3D Predictors → RNAdvisor (scoring) → Report
```

**Supported 3D predictors:**
- **RhoFold+** — deep learning, end-to-end (best for shallow MSAs)
- **SimRNA** — physics-based Monte Carlo with 2D restraints
- **Protenix** — open-source AlphaFold3 reproduction

## Install

Requires [pixi](https://pixi.sh).

### 1. Python packages & conda tools

This single command installs Python, all conda packages (Infernal, ViennaRNA, PyMOL, PyTorch), and all pip packages (Protenix, RNAdvisor, huggingface-hub):

```bash
pixi install
```

### 2. Data, models & git-based tools

For tools that need git clones, binary downloads, or large data files, use the install script:

```bash
# Install everything (Rfam, RhoFold+, SPOT-RNA, SimRNA, Docker)
./install_tools.sh --all

# Or pick specific ones
./install_tools.sh --rfam --rhofold
./install_tools.sh --docker        # Docker + Colima for RNAdvisor scoring
./install_tools.sh --simrna        # SimRNA binary (auto-detects macOS/Linux)

# Regenerate configs/local.yaml from already-installed tools
./install_tools.sh --config-only
```

This downloads databases to `~/.rnapipey/data/`, clones git-based tools to `~/.rnapipey/tools/`, and generates `configs/local.yaml` with all paths filled in. Run `./install_tools.sh --help` for the full option list.

### 3. RhoFold+ Python package

After cloning RhoFold+ with the install script, install it into the pixi environment:

```bash
pixi run pip install -e ~/.rnapipey/tools/RhoFold
```

### Notes

- **SPOT-RNA** needs a separate Python 3.6 environment (TensorFlow 1.14). The script clones the repo and prints setup instructions.
- **Docker** is required for RNAdvisor scoring. On macOS, `--docker` installs the Docker CLI + [Colima](https://github.com/abiosoft/colima) (lightweight Docker VM) via Homebrew.

## Usage

```bash
# Check which external tools are available
pixi run rnapipey check -c configs/local.yaml

# Run with specific predictors
pixi run rnapipey run input.fasta -o results/ -c configs/local.yaml --rhofold

# Run multiple predictors
pixi run rnapipey run input.fasta -o results/ -c configs/local.yaml --rhofold --simrna

# Run all predictors
pixi run rnapipey run input.fasta -o results/ -c configs/local.yaml --all

# Verbose output
pixi run rnapipey run input.fasta -o results/ -c configs/local.yaml --rhofold -v

# Regenerate report from existing results
pixi run rnapipey report results/
```

### Ensemble prediction

Generate multiple structures with different seeds for ensemble analysis and clustering:

```bash
# 20 structures with RhoFold+
pixi run rnapipey run input.fasta -o results/ -c configs/local.yaml --rhofold --nstruct 20

# 20 structures with SimRNA (physics-based MC sampling)
pixi run rnapipey run input.fasta -o results/ -c configs/local.yaml --simrna --nstruct 20
```

RhoFold+ uses **batch inference** — the model loads once and runs all seeds in a single process, so 20 structures takes ~10 min instead of ~3 hours.

When `--nstruct > 1`, the pipeline automatically clusters output structures by RMSD and scores cluster representatives.

#### Structural diversity with MC Dropout and input noise

RhoFold+ is fully deterministic — different seeds alone produce identical structures. To get meaningful ensemble diversity, use `--mc-dropout` and/or `--noise-scale`:

```bash
# MC Dropout only (re-enables dropout layers at inference)
pixi run rnapipey run input.fasta -o results/ -c configs/local.yaml \
  --rhofold --nstruct 10 --mc-dropout

# Input noise only (adds Gaussian noise to post-embedding features)
pixi run rnapipey run input.fasta -o results/ -c configs/local.yaml \
  --rhofold --nstruct 10 --noise-scale 0.1

# Both (recommended for maximum diversity)
pixi run rnapipey run input.fasta -o results/ -c configs/local.yaml \
  --rhofold --nstruct 10 --mc-dropout --noise-scale 0.1
```

The first seed (seed 0) always runs **vanilla** (no dropout, no noise) to preserve the deterministic baseline. Stochastic methods are applied only to the remaining seeds. Without these flags, behavior is unchanged from before.

### Multi-GPU

Distribute ensemble runs across multiple GPUs with `--device`:

```bash
# 4 GPUs: each loads the model once, runs 5 seeds
pixi run rnapipey run input.fasta -o results/ -c configs/local.yaml --rhofold --nstruct 20 \
  --device cuda:0,cuda:1,cuda:2,cuda:3

# Single GPU
pixi run rnapipey run input.fasta -o results/ -c configs/local.yaml --rhofold --device cuda:0
```

Multi-GPU support works for RhoFold+ and Protenix. Seeds are distributed round-robin across devices.

Or activate the environment and run directly:

```bash
pixi shell
rnapipey run input.fasta -o results/ -c configs/local.yaml --all
```

## Configuration

The install script generates `configs/local.yaml` with paths to all detected tools. You can also copy and edit `configs/default.yaml`:

```bash
pixi run rnapipey run input.fasta -c my_config.yaml --all
```

Key config fields: tool binary paths, RhoFold+ inference script & weights, SimRNA binary & data directory, RNAdvisor scoring metrics, and device (cpu/cuda).

## External tools

| Tool | Purpose | Install |
|------|---------|---------|
| [Infernal](http://eddylab.org/infernal/) + [Rfam](https://rfam.org/) | Sequence analysis, MSA | `pixi install` + `./install_tools.sh --rfam` |
| [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) | Secondary structure (RNAfold) | `pixi install` |
| [SPOT-RNA](https://github.com/jaswindersingh2/SPOT-RNA) | Pseudoknot detection (optional) | `./install_tools.sh --spotrna` |
| [RhoFold+](https://github.com/ml4bio/RhoFold) | DL 3D prediction | `./install_tools.sh --rhofold` + `pixi run pip install -e` |
| [SimRNA](https://genesilico.pl/SimRNAweb/) | Physics-based 3D prediction | `./install_tools.sh --simrna` |
| [Protenix](https://github.com/bytedance/Protenix) | AF3-class 3D prediction | `pixi install` |
| [RNAdvisor](https://github.com/EvryRNA/rnadvisor) | Model scoring (requires Docker) | `pixi install` + `./install_tools.sh --docker` |
| [PyMOL](https://pymol.org/) | Visualization (optional) | `pixi install` |
| [Docker](https://www.docker.com/) + [Colima](https://github.com/abiosoft/colima) | Container runtime for RNAdvisor | `./install_tools.sh --docker` |

## Output

```
output/
├── input/query.fasta
├── 01_sequence_analysis/    # Rfam hits, MSA
├── 02_secondary_structure/  # dot-bracket, base pair probs
├── 03_3d_prediction/        # per-predictor PDB/CIF files
├── 04_scoring/              # RNAdvisor scores + ranking
├── 05_visualization/        # summary.md + PyMOL scripts
├── logs/                    # per-tool stdout/stderr + pipeline log
└── pipeline_state.json      # tracks progress for resume
```

The pipeline saves state after each stage, so re-running the same command resumes from where it left off.

## Quick test

```bash
pixi run rnapipey run test_input.fasta -o test_results/ -c configs/local.yaml --rhofold -v
```
