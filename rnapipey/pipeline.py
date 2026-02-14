"""Pipeline orchestrator: runs stages in sequence, tracks state for resume."""

from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path
from typing import Any

from rnapipey.config import PipelineConfig
from rnapipey.tools.base import ToolResult
from rnapipey.tools.infernal import InfernalTool
from rnapipey.tools.protenix import ProtenixTool
from rnapipey.tools.rhofold import RhoFoldTool
from rnapipey.tools.rnadvisor import RNAdvisorTool
from rnapipey.tools.simrna import SimRNATool
from rnapipey.tools.spotrna import SPOTRNATool
from rnapipey.tools.viennarna import ViennaRNATool
from rnapipey.report import generate_report, generate_pymol_scripts
from rnapipey.utils import ensure_dir, read_fasta, write_fasta

logger = logging.getLogger("rnapipey")

PREDICTOR_CLASSES = {
    "rhofold": RhoFoldTool,
    "simrna": SimRNATool,
    "protenix": ProtenixTool,
}

PREDICTOR_CONFIG_ATTR = {
    "rhofold": "rhofold",
    "simrna": "simrna",
    "protenix": "protenix",
}


class Pipeline:
    def __init__(self, config: PipelineConfig, output_dir: Path):
        self.config = config
        self.output_dir = output_dir
        self.state_file = output_dir / "pipeline_state.json"
        self.state = self._load_state()
        self.logs_dir = ensure_dir(output_dir / "logs")

    def run(
        self,
        fasta_path: Path,
        predictors: list[str],
        skip_infernal: bool = False,
        run_spotrna: bool = False,
        skip_scoring: bool = False,
    ) -> None:
        """Run the full pipeline."""
        # Copy input
        input_dir = ensure_dir(self.output_dir / "input")
        query_fasta = input_dir / "query.fasta"
        if not query_fasta.exists():
            shutil.copy2(fasta_path, query_fasta)

        records = read_fasta(query_fasta)
        if not records:
            raise ValueError(f"No sequences found in {fasta_path}")
        logger.info(
            "Input: %s (%d nt)", records[0].id, len(records[0].sequence)
        )

        # Stage 1: Sequence Analysis
        msa_path: Path | None = None
        if not skip_infernal:
            msa_path = self._run_stage1(query_fasta)
        else:
            logger.info("Skipping Infernal (--skip-infernal)")

        # Stage 2: Secondary Structure
        ss_result = self._run_stage2(query_fasta, run_spotrna)
        dot_bracket = ss_result.metrics.get("dot_bracket", "") if ss_result else ""

        # Stage 3: 3D Prediction
        predictor_results = self._run_stage3(
            query_fasta, predictors, msa_path, dot_bracket
        )

        # Stage 4: Scoring
        scoring_result = None
        if not skip_scoring:
            scoring_result = self._run_stage4(predictor_results)
        else:
            logger.info("Skipping scoring (--skip-scoring)")

        # Stage 5: Report
        self._run_stage5(predictor_results, scoring_result, ss_result)

        logger.info("Pipeline complete. Results in: %s", self.output_dir)

    # -- Stage implementations --

    def _run_stage1(self, fasta_path: Path) -> Path | None:
        """Stage 1: Infernal / Rfam search."""
        if self._completed("stage1"):
            logger.info("Stage 1 (Infernal) already completed, skipping")
            # Try to recover MSA path
            msa = self.output_dir / "01_sequence_analysis" / "alignment.sto"
            return msa if msa.exists() else None

        work_dir = ensure_dir(self.output_dir / "01_sequence_analysis")
        tool = InfernalTool(self.config.tools, work_dir, self.logs_dir)

        if not tool.check():
            logger.warning(
                "Infernal/Rfam not available. Skipping sequence analysis. "
                "Set tools.rfam_cm in config to enable."
            )
            self._mark("stage1", "skipped")
            return None

        logger.info("Stage 1: Running Infernal (cmscan against Rfam)...")
        result = tool.run(fasta_path=fasta_path)
        self._mark("stage1", "completed")

        if result.success:
            family = result.metrics.get("rfam_family")
            if family:
                logger.info("Rfam family: %s", family)
            else:
                logger.info("No Rfam family match found")
            return result.output_files.get("msa")
        else:
            logger.warning("Infernal failed: %s", result.error_message)
            return None

    def _run_stage2(
        self, fasta_path: Path, run_spotrna: bool
    ) -> ToolResult | None:
        """Stage 2: Secondary structure prediction."""
        if self._completed("stage2"):
            logger.info("Stage 2 (Secondary structure) already completed, skipping")
            dot_file = self.output_dir / "02_secondary_structure" / "rnafold.dot"
            if dot_file.exists():
                content = dot_file.read_text().strip().splitlines()
                # Extract dot-bracket from last line
                for line in reversed(content):
                    if "(" in line or "." in line:
                        db = line.split()[0]
                        return ToolResult(
                            success=True,
                            metrics={"dot_bracket": db},
                        )
            return None

        work_dir = ensure_dir(self.output_dir / "02_secondary_structure")
        tool = ViennaRNATool(self.config.tools, work_dir, self.logs_dir)

        if not tool.check():
            logger.warning("RNAfold not found. Skipping secondary structure prediction.")
            self._mark("stage2", "skipped")
            return None

        logger.info("Stage 2: Running RNAfold...")
        result = tool.run(fasta_path=fasta_path)

        if result.success:
            db = result.metrics.get("dot_bracket", "")
            mfe = result.metrics.get("mfe", 0.0)
            logger.info("Secondary structure: %s (MFE: %.2f kcal/mol)", db, mfe)

        # Optional: SPOT-RNA for pseudoknots
        if run_spotrna:
            spot_tool = SPOTRNATool(self.config.tools, work_dir, self.logs_dir)
            if spot_tool.check():
                logger.info("Running SPOT-RNA for pseudoknot detection...")
                spot_result = spot_tool.run(fasta_path=fasta_path)
                if spot_result.success:
                    logger.info(
                        "SPOT-RNA structure: %s",
                        spot_result.metrics.get("dot_bracket", ""),
                    )

        self._mark("stage2", "completed")
        return result

    def _run_stage3(
        self,
        fasta_path: Path,
        predictors: list[str],
        msa_path: Path | None,
        secondary_structure: str,
    ) -> dict[str, ToolResult]:
        """Stage 3: 3D structure prediction with selected methods."""
        results: dict[str, ToolResult] = {}
        pred_dir = ensure_dir(self.output_dir / "03_3d_prediction")

        for pred_name in predictors:
            stage_key = f"stage3_{pred_name}"
            if self._completed(stage_key):
                logger.info("Stage 3 (%s) already completed, skipping", pred_name)
                # Try to recover PDB
                work_dir = pred_dir / pred_name
                pdbs = list(work_dir.glob("*.pdb")) + list(work_dir.rglob("*.cif"))
                if pdbs:
                    results[pred_name] = ToolResult(
                        success=True,
                        output_files={"pdb": pdbs[0]},
                    )
                continue

            if pred_name not in PREDICTOR_CLASSES:
                logger.error("Unknown predictor: %s", pred_name)
                continue

            work_dir = ensure_dir(pred_dir / pred_name)
            config_attr = PREDICTOR_CONFIG_ATTR[pred_name]
            tool_config = getattr(self.config.tools, config_attr)
            tool = PREDICTOR_CLASSES[pred_name](tool_config, work_dir, self.logs_dir)

            if not tool.check():
                logger.warning("%s not available, skipping", pred_name)
                self._mark(stage_key, "skipped")
                continue

            logger.info("Stage 3: Running %s...", pred_name)
            kwargs: dict[str, Any] = {"fasta_path": fasta_path}
            if pred_name in ("rhofold", "protenix") and msa_path:
                kwargs["msa_path"] = msa_path
            if pred_name == "simrna":
                kwargs["secondary_structure"] = secondary_structure

            result = tool.run(**kwargs)
            results[pred_name] = result
            self._mark(stage_key, "completed" if result.success else "failed")

            if result.success:
                pdb = result.output_files.get("pdb")
                logger.info("%s completed: %s", pred_name, pdb)
            else:
                logger.warning("%s failed: %s", pred_name, result.error_message)

        return results

    def _run_stage4(
        self, predictor_results: dict[str, ToolResult]
    ) -> ToolResult | None:
        """Stage 4: Model evaluation and scoring."""
        if self._completed("stage4"):
            logger.info("Stage 4 (Scoring) already completed, skipping")
            scores_file = self.output_dir / "04_scoring" / "rnadvisor_scores.json"
            if scores_file.exists():
                data = json.loads(scores_file.read_text())
                return ToolResult(success=True, metrics={"scores": data})
            return None

        # Collect all PDB files from predictors
        pdb_files: list[Path] = []
        for name, result in predictor_results.items():
            if result.success:
                pdb = result.output_files.get("pdb")
                if pdb and pdb.exists():
                    pdb_files.append(pdb)

        if not pdb_files:
            logger.warning("No PDB files to score")
            self._mark("stage4", "skipped")
            return None

        work_dir = ensure_dir(self.output_dir / "04_scoring")
        tool = RNAdvisorTool(self.config.tools.rnadvisor, work_dir, self.logs_dir)

        if not tool.check():
            logger.warning("RNAdvisor not available, skipping scoring")
            self._mark("stage4", "skipped")
            return None

        logger.info("Stage 4: Scoring %d models with RNAdvisor...", len(pdb_files))
        result = tool.run(pdb_files=pdb_files)
        self._mark("stage4", "completed" if result.success else "failed")

        if result.success:
            best = result.metrics.get("best_model")
            logger.info("Best model: %s", best)

        return result

    def _run_stage5(
        self,
        predictor_results: dict[str, ToolResult],
        scoring_result: ToolResult | None,
        ss_result: ToolResult | None,
    ) -> None:
        """Stage 5: Generate report and visualization scripts."""
        vis_dir = ensure_dir(self.output_dir / "05_visualization")
        fasta_path = self.output_dir / "input" / "query.fasta"

        logger.info("Stage 5: Generating report and visualization scripts...")
        generate_report(
            vis_dir, fasta_path, predictor_results, scoring_result, ss_result
        )
        generate_pymol_scripts(vis_dir, predictor_results, scoring_result)
        self._mark("stage5", "completed")

    # -- State management --

    def _load_state(self) -> dict[str, str]:
        if self.state_file.exists():
            return json.loads(self.state_file.read_text())
        return {}

    def _save_state(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.state_file.write_text(json.dumps(self.state, indent=2))

    def _completed(self, stage: str) -> bool:
        return self.state.get(stage) == "completed"

    def _mark(self, stage: str, status: str) -> None:
        self.state[stage] = status
        self._save_state()
