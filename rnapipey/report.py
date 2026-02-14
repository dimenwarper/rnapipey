"""Report generation: markdown summary and PyMOL visualization scripts."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Any

from rnapipey.tools.base import ToolResult
from rnapipey.utils import read_fasta


def generate_report(
    vis_dir: Path,
    fasta_path: Path,
    predictor_results: dict[str, ToolResult],
    scoring_result: ToolResult | None,
    ss_result: ToolResult | None,
    ensemble_result: Any = None,
) -> Path:
    """Write a markdown summary report."""
    report_path = vis_dir / "summary.md"
    lines: list[str] = []

    # Header
    lines.append("# rnapipey — Structure Prediction Report")
    lines.append(f"\nGenerated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    # Input info
    records = read_fasta(fasta_path) if fasta_path.exists() else []
    if records:
        seq = records[0].sequence
        gc = (seq.count("G") + seq.count("C")) / len(seq) * 100
        lines.append("## Input")
        lines.append(f"- **Sequence ID**: {records[0].id}")
        lines.append(f"- **Length**: {len(seq)} nt")
        lines.append(f"- **GC content**: {gc:.1f}%")
        lines.append(f"- **Sequence**: `{seq}`")
        lines.append("")

    # Secondary structure
    if ss_result and ss_result.success:
        db = ss_result.metrics.get("dot_bracket", "")
        mfe = ss_result.metrics.get("mfe", "N/A")
        lines.append("## Secondary Structure (RNAfold)")
        lines.append(f"- **Dot-bracket**: `{db}`")
        lines.append(f"- **MFE**: {mfe} kcal/mol")
        lines.append("")

    # 3D Predictions
    lines.append("## 3D Structure Predictions")
    lines.append("")
    if predictor_results:
        lines.append("| Predictor | Status | Output | Structures |")
        lines.append("|-----------|--------|--------|------------|")
        for name, result in predictor_results.items():
            status = "OK" if result.success else "FAILED"
            pdb = result.output_files.get("pdb", "—")
            if pdb and isinstance(pdb, Path):
                pdb = pdb.name
            all_pdbs = result.output_files.get("all_pdbs", [])
            n_structs = len(all_pdbs) if all_pdbs else (1 if result.success else 0)
            lines.append(f"| {name} | {status} | {pdb} | {n_structs} |")
        lines.append("")
    else:
        lines.append("No 3D predictions were run.\n")

    # Confidence Metrics (per-predictor pLDDT/pTM)
    confidence_rows = []
    for name, result in predictor_results.items():
        if not result.success:
            continue
        metrics = result.metrics
        row: dict[str, str] = {"Predictor": name}
        if "plddt_mean" in metrics:
            row["pLDDT (mean)"] = f"{metrics['plddt_mean']:.1f}"
        if "ptm" in metrics:
            row["pTM"] = f"{metrics['ptm']:.3f}"
        if "iptm" in metrics:
            row["ipTM"] = f"{metrics['iptm']:.3f}"
        if "ranking_score" in metrics:
            row["Ranking Score"] = f"{metrics['ranking_score']:.3f}"
        if len(row) > 1:
            confidence_rows.append(row)

    if confidence_rows:
        lines.append("## Confidence Metrics")
        lines.append("")
        all_cols = []
        for row in confidence_rows:
            for k in row:
                if k not in all_cols:
                    all_cols.append(k)
        lines.append("| " + " | ".join(all_cols) + " |")
        lines.append("|" + "|".join(["-------"] * len(all_cols)) + "|")
        for row in confidence_rows:
            vals = [row.get(c, "—") for c in all_cols]
            lines.append("| " + " | ".join(vals) + " |")
        lines.append("")

    # Ensemble Analysis
    if ensemble_result is not None:
        _append_ensemble_section(lines, ensemble_result)

    # Scoring
    if scoring_result and scoring_result.success:
        lines.append("## Model Scoring (RNAdvisor)")
        lines.append("")
        scores = scoring_result.metrics.get("scores", {})
        ranking = scoring_result.metrics.get("ranking", [])
        best = scoring_result.metrics.get("best_model")

        if best:
            lines.append(f"**Best model**: {best}\n")

        if ranking:
            lines.append("**Ranking**:")
            for i, name in enumerate(ranking):
                lines.append(f"{i+1}. {name}")
            lines.append("")

        if scores:
            # Build score table
            all_metrics = set()
            for s in scores.values():
                all_metrics.update(s.keys())
            metrics_list = sorted(all_metrics)

            header = "| Model | " + " | ".join(metrics_list) + " |"
            sep = "|-------|" + "|".join(["-------"] * len(metrics_list)) + "|"
            lines.append(header)
            lines.append(sep)
            for model, model_scores in scores.items():
                vals = [f"{model_scores.get(m, 'N/A')}" for m in metrics_list]
                lines.append(f"| {model} | " + " | ".join(vals) + " |")
            lines.append("")

    # Footer
    lines.append("---")
    lines.append("*Generated by [rnapipey](https://github.com/your-repo/rnapipey)*")

    report_path.write_text("\n".join(lines) + "\n")
    return report_path


def _append_ensemble_section(lines: list[str], ensemble_result: Any) -> None:
    """Append ensemble analysis sections to the report."""
    clusters = ensemble_result.clusters
    n_total = len(ensemble_result.pdb_files)
    n_clusters = len(clusters)
    n_consensus = sum(1 for c in clusters if c.is_consensus)

    lines.append("## Ensemble Analysis")
    lines.append("")
    lines.append(f"- **Total structures**: {n_total}")
    lines.append(f"- **Number of clusters**: {n_clusters}")
    lines.append(f"- **Consensus clusters**: {n_consensus} (structures from 2+ predictors)")
    lines.append("")

    # Cluster summary table
    lines.append("### Cluster Summary")
    lines.append("")
    lines.append("| Cluster | Size | Predictors | Representative | Mean RMSD (A) | Consensus |")
    lines.append("|---------|------|------------|----------------|---------------|-----------|")
    for c in clusters:
        preds = ", ".join(sorted(set(c.member_predictors)))
        rep_name = c.representative.name
        consensus = "Yes" if c.is_consensus else "No"
        lines.append(
            f"| {c.cluster_id} | {len(c.members)} | {preds} | {rep_name} | {c.mean_rmsd:.2f} | {consensus} |"
        )
    lines.append("")

    # Consensus structures section
    consensus_clusters = [c for c in clusters if c.is_consensus]
    if consensus_clusters:
        lines.append("### Consensus Structures")
        lines.append("")
        lines.append(
            "These clusters contain structures from multiple predictors, "
            "indicating higher confidence in the fold:"
        )
        lines.append("")
        for c in consensus_clusters:
            preds = ", ".join(sorted(set(c.member_predictors)))
            lines.append(
                f"- **Cluster {c.cluster_id}**: {len(c.members)} structures "
                f"from {preds} (mean RMSD: {c.mean_rmsd:.2f} A, "
                f"representative: `{c.representative.name}`)"
            )
        lines.append("")


def generate_pymol_scripts(
    vis_dir: Path,
    predictor_results: dict[str, ToolResult],
    scoring_result: ToolResult | None,
    ensemble_result: Any = None,
) -> None:
    """Generate PyMOL .pml scripts for visualization."""
    # Collect all PDB paths
    pdbs: list[tuple[str, Path]] = []
    for name, result in predictor_results.items():
        if result.success:
            pdb = result.output_files.get("pdb")
            if pdb and isinstance(pdb, Path) and pdb.exists():
                pdbs.append((name, pdb))

    if not pdbs:
        return

    # Colors for different predictors
    colors = ["marine", "red", "forest", "orange", "purple", "cyan"]

    # view_all.pml: load all models, align, color by predictor
    all_lines = [
        "# rnapipey: Load and compare all predicted models",
        "bg_color white",
        "set cartoon_ring_mode, 3",
        "set cartoon_nucleic_acid_mode, 4",
        "",
    ]
    for i, (name, pdb) in enumerate(pdbs):
        color = colors[i % len(colors)]
        all_lines.append(f"load {pdb}, {name}")
        all_lines.append(f"color {color}, {name}")
        all_lines.append(f"show cartoon, {name}")
        all_lines.append("")

    if len(pdbs) > 1:
        ref = pdbs[0][0]
        for name, _ in pdbs[1:]:
            all_lines.append(f"align {name}, {ref}")
        all_lines.append("")

    all_lines.extend([
        "zoom",
        "set ray_opaque_background, 0",
        "# Use: ray 1200, 900 to render",
    ])
    (vis_dir / "view_all.pml").write_text("\n".join(all_lines) + "\n")

    # view_best.pml: load best model only
    best_name = None
    if scoring_result and scoring_result.success:
        best_name = scoring_result.metrics.get("best_model")

    # Fall back to first predictor
    if not best_name:
        best_name = pdbs[0][0]

    best_pdb = None
    for name, pdb in pdbs:
        if name == best_name or pdb.name == best_name:
            best_pdb = pdb
            break
    if not best_pdb:
        best_pdb = pdbs[0][1]

    best_lines = [
        "# rnapipey: View best predicted model",
        "bg_color white",
        "set cartoon_ring_mode, 3",
        "set cartoon_nucleic_acid_mode, 4",
        "",
        f"load {best_pdb}, best_model",
        "show cartoon, best_model",
        "color marine, best_model",
        "",
        "# Color by element for bases",
        "util.cnc best_model",
        "",
        "# Show nucleotide bases as sticks",
        "select bases, best_model and (name C2+C4+C5+C6+C8+N1+N2+N3+N4+N6+N7+N9+O2+O4+O6)",
        "show sticks, bases",
        "",
        "zoom",
        "set ray_opaque_background, 0",
        "# Use: ray 1200, 900 to render",
    ]
    (vis_dir / "view_best.pml").write_text("\n".join(best_lines) + "\n")

    # view_clusters.pml: ensemble visualization (when clustering was performed)
    if ensemble_result is not None:
        _generate_cluster_pymol(vis_dir, ensemble_result)


def _generate_cluster_pymol(vis_dir: Path, ensemble_result: Any) -> None:
    """Generate PyMOL script for cluster visualization."""
    cluster_colors = [
        "marine", "red", "forest", "orange", "purple", "cyan",
        "yellow", "salmon", "lime", "slate",
    ]
    rep_color = "white"

    lines = [
        "# rnapipey: Ensemble cluster visualization",
        "# Structures colored by cluster, representatives highlighted",
        "bg_color white",
        "set cartoon_ring_mode, 3",
        "set cartoon_nucleic_acid_mode, 4",
        "",
    ]

    loaded: list[str] = []
    ref_name = None

    for cluster in ensemble_result.clusters:
        color = cluster_colors[(cluster.cluster_id - 1) % len(cluster_colors)]

        for member in cluster.members:
            if not member.exists():
                continue
            obj_name = f"c{cluster.cluster_id}_{member.stem}"
            lines.append(f"load {member}, {obj_name}")
            lines.append(f"color {color}, {obj_name}")
            lines.append(f"show cartoon, {obj_name}")

            # Highlight representative with thicker cartoon
            if member == cluster.representative:
                lines.append(f"set cartoon_tube_radius, 0.4, {obj_name}")
                lines.append(f"# Representative of cluster {cluster.cluster_id}")

            loaded.append(obj_name)
            if ref_name is None:
                ref_name = obj_name
            lines.append("")

    # Align all to first
    if ref_name and len(loaded) > 1:
        lines.append("# Align all structures")
        for obj in loaded[1:]:
            lines.append(f"align {obj}, {ref_name}")
        lines.append("")

    lines.extend([
        "zoom",
        "set ray_opaque_background, 0",
        "# Use: ray 1200, 900 to render",
    ])
    (vis_dir / "view_clusters.pml").write_text("\n".join(lines) + "\n")
