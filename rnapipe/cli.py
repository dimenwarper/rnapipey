"""CLI entry point for rnapipe."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

from rnapipe import __version__
from rnapipe.config import load_config
from rnapipe.pipeline import Pipeline
from rnapipe.utils import setup_logging

app = typer.Typer(
    name="rnapipe",
    help="RNA 3D structure prediction pipeline",
    add_completion=False,
)
console = Console()


def version_callback(value: bool) -> None:
    if value:
        console.print(f"rnapipe {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        False, "--version", "-V", callback=version_callback, is_eager=True,
        help="Show version and exit.",
    ),
) -> None:
    pass


@app.command()
def run(
    input_fasta: Path = typer.Argument(
        ..., exists=True, readable=True, help="Input FASTA file (single RNA sequence)."
    ),
    output_dir: Path = typer.Option(
        "./rnapipe_output", "-o", "--output-dir", help="Output directory."
    ),
    config: Optional[Path] = typer.Option(
        None, "-c", "--config", help="YAML config file."
    ),
    # Predictor selection
    rhofold: bool = typer.Option(False, "--rhofold", help="Run RhoFold+."),
    simrna: bool = typer.Option(False, "--simrna", help="Run SimRNA."),
    protenix: bool = typer.Option(False, "--protenix", help="Run Protenix."),
    all_predictors: bool = typer.Option(
        False, "--all", help="Run all available predictors."
    ),
    # Optional stages
    skip_infernal: bool = typer.Option(
        False, "--skip-infernal", help="Skip Infernal/Rfam search."
    ),
    spotrna: bool = typer.Option(
        False, "--spotrna", help="Also run SPOT-RNA for pseudoknots."
    ),
    skip_scoring: bool = typer.Option(
        False, "--skip-scoring", help="Skip model scoring with RNAdvisor."
    ),
    # Misc
    verbose: bool = typer.Option(False, "-v", "--verbose", help="Verbose logging."),
) -> None:
    """Run the full RNA 3D structure prediction pipeline."""
    output_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(log_file=output_dir / "logs" / "rnapipe.log", verbose=verbose)

    cfg = load_config(config)

    # Determine which predictors to run
    predictors: list[str] = []
    if all_predictors:
        predictors = ["rhofold", "simrna", "protenix"]
    else:
        if rhofold:
            predictors.append("rhofold")
        if simrna:
            predictors.append("simrna")
        if protenix:
            predictors.append("protenix")

    if not predictors:
        console.print(
            "[yellow]Warning:[/yellow] No 3D predictors selected. "
            "Use --rhofold, --simrna, --protenix, or --all."
        )
        console.print("Running stages 1-2 only (sequence analysis + secondary structure).")

    console.print(f"\n[bold]rnapipe v{__version__}[/bold]")
    console.print(f"  Input:      {input_fasta}")
    console.print(f"  Output:     {output_dir}")
    console.print(f"  Predictors: {', '.join(predictors) if predictors else 'none'}")
    console.print()

    pipeline = Pipeline(cfg, output_dir.resolve())
    pipeline.run(
        fasta_path=input_fasta.resolve(),
        predictors=predictors,
        skip_infernal=skip_infernal,
        run_spotrna=spotrna,
        skip_scoring=skip_scoring,
    )
    console.print(f"\n[green]Done![/green] Results in: {output_dir}")


@app.command()
def check(
    config: Optional[Path] = typer.Option(
        None, "-c", "--config", help="YAML config file."
    ),
) -> None:
    """Check that required external tools are installed and accessible."""
    from rnapipe.tools.infernal import InfernalTool
    from rnapipe.tools.viennarna import ViennaRNATool
    from rnapipe.tools.spotrna import SPOTRNATool
    from rnapipe.tools.rhofold import RhoFoldTool
    from rnapipe.tools.simrna import SimRNATool
    from rnapipe.tools.protenix import ProtenixTool
    from rnapipe.tools.rnadvisor import RNAdvisorTool
    from rnapipe.utils import which

    cfg = load_config(config)
    tmp = Path("/tmp/rnapipe_check")
    tmp.mkdir(exist_ok=True)

    table = Table(title="rnapipe — Tool Availability")
    table.add_column("Tool", style="bold")
    table.add_column("Status")
    table.add_column("Notes")

    checks = [
        ("Infernal (cmscan)", InfernalTool(cfg.tools, tmp), "Rfam search + MSA"),
        ("ViennaRNA (RNAfold)", ViennaRNATool(cfg.tools, tmp), "Secondary structure"),
        ("SPOT-RNA", SPOTRNATool(cfg.tools, tmp), "Pseudoknot detection (optional)"),
        ("RhoFold+", RhoFoldTool(cfg.tools.rhofold, tmp), "DL 3D prediction"),
        ("SimRNA", SimRNATool(cfg.tools.simrna, tmp), "Physics-based 3D prediction"),
        ("Protenix", ProtenixTool(cfg.tools.protenix, tmp), "AF3-class 3D prediction"),
        ("RNAdvisor", RNAdvisorTool(cfg.tools.rnadvisor, tmp), "Model scoring"),
    ]

    for name, tool, notes in checks:
        available = tool.check()
        status = "[green]OK[/green]" if available else "[red]NOT FOUND[/red]"
        table.add_row(name, status, notes)

    # Also check for PyMOL
    pymol_ok = which("pymol") is not None
    table.add_row(
        "PyMOL",
        "[green]OK[/green]" if pymol_ok else "[yellow]NOT FOUND[/yellow]",
        "Visualization (optional)",
    )

    console.print()
    console.print(table)
    console.print()


@app.command()
def report(
    output_dir: Path = typer.Argument(
        ..., exists=True, help="Existing pipeline output directory."
    ),
) -> None:
    """Re-generate the summary report and PyMOL scripts from existing results."""
    import json
    from rnapipe.report import generate_report, generate_pymol_scripts
    from rnapipe.tools.base import ToolResult
    from rnapipe.utils import ensure_dir

    vis_dir = ensure_dir(output_dir / "05_visualization")
    fasta_path = output_dir / "input" / "query.fasta"

    # Reconstruct predictor results from existing PDB files
    pred_dir = output_dir / "03_3d_prediction"
    predictor_results: dict[str, ToolResult] = {}
    if pred_dir.exists():
        for subdir in pred_dir.iterdir():
            if subdir.is_dir():
                pdbs = list(subdir.glob("*.pdb")) + list(subdir.rglob("*.cif"))
                if pdbs:
                    predictor_results[subdir.name] = ToolResult(
                        success=True,
                        output_files={"pdb": pdbs[0]},
                    )

    # Reconstruct scoring result
    scoring_result = None
    scores_file = output_dir / "04_scoring" / "rnadvisor_scores.json"
    if scores_file.exists():
        scores = json.loads(scores_file.read_text())
        scoring_result = ToolResult(
            success=True,
            metrics={"scores": scores},
        )

    # Reconstruct SS result
    ss_result = None
    dot_file = output_dir / "02_secondary_structure" / "rnafold.dot"
    if dot_file.exists():
        content = dot_file.read_text().strip().splitlines()
        for line in reversed(content):
            if "(" in line or "." in line:
                parts = line.rsplit(" (", 1)
                if len(parts) == 2:
                    db = parts[0].strip()
                    try:
                        mfe = float(parts[1].rstrip(")"))
                    except ValueError:
                        mfe = 0.0
                    ss_result = ToolResult(
                        success=True,
                        metrics={"dot_bracket": db, "mfe": mfe},
                    )
                    break

    generate_report(vis_dir, fasta_path, predictor_results, scoring_result, ss_result)
    generate_pymol_scripts(vis_dir, predictor_results, scoring_result)
    console.print(f"[green]Report regenerated:[/green] {vis_dir / 'summary.md'}")
