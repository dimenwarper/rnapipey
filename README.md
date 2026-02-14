# rnapipe

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

```bash
uv sync
```

## Usage

```bash
# Check which external tools are available
rnapipe check

# Run with specific predictors
rnapipe run input.fasta -o results/ --rhofold --simrna

# Run all predictors
rnapipe run input.fasta -o results/ --all

# With options
rnapipe run input.fasta -o results/ --all --skip-infernal --spotrna -v

# Regenerate report from existing results
rnapipe report results/
```

## Configuration

Copy and edit `configs/default.yaml` to point to your tool installations:

```bash
rnapipe run input.fasta -c my_config.yaml --all
```

Key config fields are paths to external tool binaries/scripts (RhoFold+ inference script, SimRNA binary, Rfam database, etc).

## External tools

These must be installed separately — rnapipe wraps them via subprocess:

| Tool | Purpose | Install |
|------|---------|---------|
| [Infernal](http://eddylab.org/infernal/) + [Rfam](https://rfam.org/) | Sequence analysis, MSA | `conda install -c bioconda infernal` |
| [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) | Secondary structure (RNAfold) | `conda install -c bioconda viennarna` |
| [SPOT-RNA](https://github.com/jaswindersingh2/SPOT-RNA) | Pseudoknot detection (optional) | See repo |
| [RhoFold+](https://github.com/ml4bio/RhoFold) | DL 3D prediction | See repo |
| [SimRNA](https://genesilico.pl/SimRNAweb/) | Physics-based 3D prediction | See website |
| [Protenix](https://github.com/bytedance/Protenix) | AF3-class 3D prediction | See repo |
| [RNAdvisor](https://github.com/EvryRNA/rnadvisor) | Model scoring | Docker or pip |
| [PyMOL](https://pymol.org/) | Visualization (optional) | `conda install -c conda-forge pymol-open-source` |

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
