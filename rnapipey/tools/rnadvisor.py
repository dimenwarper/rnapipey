"""RNAdvisor 2 wrapper for RNA 3D model scoring and evaluation."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

from rnapipey.config import RNAdvisorConfig
from rnapipey.tools.base import BaseTool, ToolResult
from rnapipey.utils import which

logger = logging.getLogger("rnapipey")


class RNAdvisorTool(BaseTool):
    config: RNAdvisorConfig

    @property
    def name(self) -> str:
        return "rnadvisor"

    def check(self) -> bool:
        if self.config.docker:
            return which("docker") is not None
        return which("rnadvisor") is not None

    def run(self, **kwargs: Any) -> ToolResult:
        pdb_files: list[Path] = kwargs["pdb_files"]
        if not pdb_files:
            return ToolResult(success=False, error_message="No PDB files to score")

        scores: dict[str, dict[str, float]] = {}
        for pdb in pdb_files:
            if pdb is None or not pdb.exists():
                continue
            model_scores = self._score_single(pdb)
            if model_scores:
                scores[pdb.name] = model_scores

        if not scores:
            return ToolResult(
                success=False,
                error_message="No models could be scored",
            )

        # Consensus ranking: average normalized rank across metrics
        ranking = self._consensus_rank(scores)

        # Write results
        scores_file = self.work_dir / "rnadvisor_scores.json"
        scores_file.write_text(json.dumps(scores, indent=2))

        ranking_file = self.work_dir / "ranking.txt"
        ranking_lines = [f"{i+1}. {name} (avg_rank: {rank:.2f})"
                         for i, (name, rank) in enumerate(ranking)]
        ranking_file.write_text("\n".join(ranking_lines) + "\n")

        return ToolResult(
            success=True,
            output_files={
                "scores_json": scores_file,
                "ranking": ranking_file,
            },
            metrics={
                "scores": scores,
                "ranking": [name for name, _ in ranking],
                "best_model": ranking[0][0] if ranking else None,
            },
        )

    def _score_single(self, pdb: Path) -> dict[str, float]:
        """Score a single PDB file with RNAdvisor."""
        import csv
        import shutil

        scores_str = ",".join(self.config.metrics)

        # RNAdvisor expects --pred_dir (a directory of PDBs)
        # Stage each PDB into its own temp dir
        staging_dir = self.work_dir / f"_stage_{pdb.stem}"
        staging_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(pdb, staging_dir / pdb.name)

        out_csv = self.work_dir / f"scores_{pdb.stem}.csv"

        cmd = [
            "rnadvisor",
            "--pred_dir", str(staging_dir),
            "--scores", scores_str,
            "--out_path", str(out_csv),
        ]

        result = self._run_cmd(cmd)
        if result.returncode != 0:
            logger.warning("RNAdvisor failed for %s: %s", pdb.name, result.stderr[:200])
            return {}

        # Parse CSV output
        try:
            if out_csv.exists():
                with open(out_csv) as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        return {k: float(v) for k, v in row.items()
                                if k not in ("", "name", "pdb", "file") and v}
        except (ValueError, OSError) as e:
            logger.warning("Could not parse RNAdvisor output for %s: %s", pdb.name, e)
        return {}

    def _consensus_rank(
        self, scores: dict[str, dict[str, float]]
    ) -> list[tuple[str, float]]:
        """Rank models by average rank across all metrics.

        Lower scores are better for energy-based metrics (rsRNASP, DFIRE, RASP).
        Higher scores are better for MCQ. We handle this by ranking direction.
        """
        # Metrics where lower = better
        lower_is_better = {"rsRNASP", "DFIRE", "RASP", "DFIRE-RNA"}
        models = list(scores.keys())
        if not models:
            return []

        all_metrics = set()
        for s in scores.values():
            all_metrics.update(s.keys())

        rank_sums: dict[str, float] = {m: 0.0 for m in models}
        n_metrics = 0

        for metric in all_metrics:
            # Get values for this metric across models
            vals = [(model, scores[model].get(metric)) for model in models]
            vals = [(m, v) for m, v in vals if v is not None]
            if not vals:
                continue

            reverse = metric not in lower_is_better
            sorted_vals = sorted(vals, key=lambda x: x[1], reverse=reverse)
            for rank, (model, _) in enumerate(sorted_vals):
                rank_sums[model] += rank + 1
            n_metrics += 1

        if n_metrics == 0:
            return [(m, 0.0) for m in models]

        avg_ranks = [(m, rank_sums[m] / n_metrics) for m in models]
        return sorted(avg_ranks, key=lambda x: x[1])
