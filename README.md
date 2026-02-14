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

Requires [pixi](https://pixi.sh). This single command installs Python, all conda packages (Infernal, ViennaRNA, PyMOL), all pip packages (Protenix, RNAdvisor), and rnapipey itself:

```bash
pixi install
```

### Data & git-based tools

For tools that need git clones or large data downloads (Rfam, RhoFold+, SPOT-RNA), use the install script:

```bash
# Download/clone everything
./install_tools.sh --all

# Or pick specific ones
./install_tools.sh --rfam --rhofold

# Regenerate configs/local.yaml from already-installed tools
./install_tools.sh --config-only
```

This downloads databases to `~/.rnapipey/data/`, clones git-based tools to `~/.rnapipey/tools/`, and generates `configs/local.yaml` with all paths filled in. Run `./install_tools.sh --help` for the full option list.

**Note:** SimRNA requires a manual download (academic license from [genesilico.pl](https://genesilico.pl/SimRNAweb/)), and SPOT-RNA needs a separate Python 3.6 environment — the script will clone it and print setup instructions.

## Usage

```bash
# Check which external tools are available
pixi run rnapipey check

# Run with specific predictors
pixi run rnapipey run input.fasta -o results/ --rhofold --simrna

# Run all predictors
pixi run rnapipey run input.fasta -o results/ --all

# With options
pixi run rnapipey run input.fasta -o results/ --all --skip-infernal --spotrna -v

# Regenerate report from existing results
pixi run rnapipey report results/
```

Or activate the environment and run directly:

```bash
pixi shell
rnapipey run input.fasta -o results/ --all
```

## Configuration

Copy and edit `configs/default.yaml` to point to your tool installations:

```bash
pixi run rnapipey run input.fasta -c my_config.yaml --all
```

Key config fields are paths to external tool binaries/scripts (RhoFold+ inference script, SimRNA binary, Rfam database, etc).

## External tools

These are installed by pixi or the install script — rnapipey wraps them via subprocess:

| Tool | Purpose | Install |
|------|---------|---------|
| [Infernal](http://eddylab.org/infernal/) + [Rfam](https://rfam.org/) | Sequence analysis, MSA | `pixi install` (conda) + `./install_tools.sh --rfam` (data) |
| [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) | Secondary structure (RNAfold) | `pixi install` (conda) |
| [SPOT-RNA](https://github.com/jaswindersingh2/SPOT-RNA) | Pseudoknot detection (optional) | `./install_tools.sh --spotrna` |
| [RhoFold+](https://github.com/ml4bio/RhoFold) | DL 3D prediction | `./install_tools.sh --rhofold` |
| [SimRNA](https://genesilico.pl/SimRNAweb/) | Physics-based 3D prediction | Manual (academic license) |
| [Protenix](https://github.com/bytedance/Protenix) | AF3-class 3D prediction | `pixi install` (pip) |
| [RNAdvisor](https://github.com/EvryRNA/rnadvisor) | Model scoring | `pixi install` (pip) |
| [PyMOL](https://pymol.org/) | Visualization (optional) | `pixi install` (conda) |

## Output

```
output/
├── input/query.fasta
├── 01_sequence_analysis/    # Rfam hits, MSA
├── 02_secondary_structure/  # dot-bracket, base pair probs
├── 03_3d_prediction/        # per-predictor PDB/CIF files
├── 04_scoring/              # RNAdvisor scores + ranking
├── 05_visualization/        # summary.md + PyMOL scripts
├── logs/
└── pipeline_state.json      # tracks progress for resume
```

The pipeline saves state after each stage, so re-running the same command resumes from where it left off.
