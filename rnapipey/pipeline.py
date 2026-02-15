"""Pipeline orchestrator: runs stages in sequence, tracks state for resume."""

from __future__ import annotations

import json
import logging
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
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
        nstruct: int = 1,
        devices: list[str] | None = None,
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
            query_fasta, predictors, msa_path, dot_bracket,
            nstruct=nstruct, devices=devices or [],
        )

        # Stage 3b: Ensemble clustering (when nstruct > 1)
        ensemble_result = None
        if nstruct > 1 and predictor_results:
            ensemble_result = self._run_stage3b_clustering(
                predictor_results, self.config.ensemble.cluster_cutoff
            )

        # Stage 4: Scoring
        scoring_result = None
        if not skip_scoring:
            scoring_result = self._run_stage4(
                predictor_results, ensemble_result
            )
        else:
            logger.info("Skipping scoring (--skip-scoring)")

        # Stage 5: Report
        self._run_stage5(
            predictor_results, scoring_result, ss_result, ensemble_result
        )

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
        nstruct: int = 1,
        devices: list[str] | None = None,
    ) -> dict[str, ToolResult]:
        """Stage 3: 3D structure prediction with selected methods."""
        devices = devices or []
        results: dict[str, ToolResult] = {}
        pred_dir = ensure_dir(self.output_dir / "03_3d_prediction")

        for pred_name in predictors:
            stage_key = f"stage3_{pred_name}"
            if self._completed(stage_key):
                logger.info("Stage 3 (%s) already completed, skipping", pred_name)
                # Try to recover PDB files
                work_dir = pred_dir / pred_name
                pdbs = sorted(work_dir.rglob("*.pdb")) + sorted(work_dir.rglob("*.cif"))
                if pdbs:
                    results[pred_name] = ToolResult(
                        success=True,
                        output_files={"pdb": pdbs[0], "all_pdbs": pdbs},
                    )
                continue

            if pred_name not in PREDICTOR_CLASSES:
                logger.error("Unknown predictor: %s", pred_name)
                continue

            config_attr = PREDICTOR_CONFIG_ATTR[pred_name]
            tool_config = getattr(self.config.tools, config_attr)

            if pred_name == "rhofold" and nstruct > 1:
                result = self._run_rhofold_ensemble(
                    fasta_path, msa_path, tool_config, pred_dir,
                    nstruct=nstruct, devices=devices,
                )
            elif pred_name == "protenix":
                result = self._run_protenix(
                    fasta_path, msa_path, tool_config, pred_dir,
                    nstruct=nstruct, devices=devices,
                )
            elif pred_name == "simrna":
                result = self._run_simrna(
                    fasta_path, secondary_structure, tool_config, pred_dir,
                    nstruct=nstruct,
                )
            else:
                # Single run (rhofold nstruct=1 or unknown)
                work_dir = ensure_dir(pred_dir / pred_name)
                tool = PREDICTOR_CLASSES[pred_name](tool_config, work_dir, self.logs_dir)
                if not tool.check():
                    logger.warning("%s not available, skipping", pred_name)
                    self._mark(stage_key, "skipped")
                    continue

                logger.info("Stage 3: Running %s...", pred_name)
                kwargs: dict[str, Any] = {"fasta_path": fasta_path}
                if pred_name in ("rhofold", "protenix") and msa_path:
                    kwargs["msa_path"] = msa_path
                if devices:
                    kwargs["device"] = devices[0]
                result = tool.run(**kwargs)

            results[pred_name] = result
            self._mark(stage_key, "completed" if result.success else "failed")

            if result.success:
                pdb = result.output_files.get("pdb")
                all_pdbs = result.output_files.get("all_pdbs", [])
                logger.info(
                    "%s completed: %s (%d structures)",
                    pred_name, pdb, len(all_pdbs) if all_pdbs else 1,
                )
            else:
                logger.warning("%s failed: %s", pred_name, result.error_message)

        return results

    def _run_rhofold_ensemble(
        self,
        fasta_path: Path,
        msa_path: Path | None,
        tool_config: Any,
        pred_dir: Path,
        nstruct: int,
        devices: list[str] | None = None,
    ) -> ToolResult:
        """Run RhoFold+ ensemble using batch inference (load model once per GPU)."""
        devices = devices or []
        rhofold_dir = ensure_dir(pred_dir / "rhofold")
        all_seeds = list(range(nstruct))

        # Availability check
        check_tool = RhoFoldTool(tool_config, rhofold_dir, self.logs_dir)
        if not check_tool.check():
            return ToolResult(
                success=False, error_message="RhoFold+ not available",
            )

        n_gpus = max(len(devices), 1)

        if n_gpus <= 1:
            # Single GPU (or CPU): one batch call with all seeds
            logger.info(
                "Stage 3: Running RhoFold+ batch ensemble (%d structures, 1 model load)...",
                nstruct,
            )
            tool = RhoFoldTool(tool_config, rhofold_dir, self.logs_dir)
            kwargs: dict[str, Any] = {
                "fasta_path": fasta_path,
                "seeds": all_seeds,
                "output_base_dir": rhofold_dir,
            }
            if msa_path:
                kwargs["msa_path"] = msa_path
            if devices:
                kwargs["device"] = devices[0]
            return tool.run_batch(**kwargs)

        # Multi-GPU: split seeds round-robin, one batch per GPU in parallel
        seed_groups: list[list[int]] = [[] for _ in range(n_gpus)]
        for i, seed in enumerate(all_seeds):
            seed_groups[i % n_gpus].append(seed)

        logger.info(
            "Stage 3: Running RhoFold+ batch ensemble (%d structures across %d GPUs, %d model loads)...",
            nstruct, n_gpus, n_gpus,
        )

        def _run_gpu_batch(gpu_idx: int, seeds: list[int]) -> ToolResult:
            device = devices[gpu_idx]
            gpu_dir = ensure_dir(rhofold_dir / f"gpu_{gpu_idx}")
            tool = RhoFoldTool(tool_config, gpu_dir, self.logs_dir)
            logger.info(
                "Running RhoFold+ batch (seeds %s) on %s...", seeds, device,
            )
            kw: dict[str, Any] = {
                "fasta_path": fasta_path,
                "seeds": seeds,
                "output_base_dir": gpu_dir,
                "device": device,
            }
            if msa_path:
                kw["msa_path"] = msa_path
            return tool.run_batch(**kw)

        gpu_results: list[ToolResult] = []
        with ThreadPoolExecutor(max_workers=n_gpus) as executor:
            futures = {
                executor.submit(_run_gpu_batch, idx, grp): idx
                for idx, grp in enumerate(seed_groups)
                if grp
            }
            for future in as_completed(futures):
                gpu_results.append(future.result())

        # Merge results from all GPUs
        all_pdbs: list[Path] = []
        total_time = 0.0
        any_success = False
        all_plddt: list[float] = []

        for r in gpu_results:
            total_time += r.runtime_seconds
            if r.success:
                any_success = True
                pdbs = r.output_files.get("all_pdbs", [])
                if pdbs:
                    all_pdbs.extend(pdbs)
                elif r.output_files.get("pdb"):
                    all_pdbs.append(r.output_files["pdb"])
                per_run = r.metrics.get("plddt_per_run", [])
                all_plddt.extend(per_run)
            else:
                logger.warning("RhoFold+ GPU batch failed: %s", r.error_message)

        if not any_success:
            return ToolResult(
                success=False,
                error_message="All RhoFold+ GPU batches failed",
                runtime_seconds=total_time,
            )

        merged_metrics: dict[str, Any] = {}
        if all_plddt:
            merged_metrics["plddt_mean"] = float(sum(all_plddt) / len(all_plddt))
            merged_metrics["plddt_per_run"] = all_plddt

        return ToolResult(
            success=True,
            output_files={"pdb": all_pdbs[0], "all_pdbs": all_pdbs},
            metrics=merged_metrics,
            runtime_seconds=total_time,
        )

    def _run_protenix(
        self,
        fasta_path: Path,
        msa_path: Path | None,
        tool_config: Any,
        pred_dir: Path,
        nstruct: int,
        devices: list[str] | None = None,
    ) -> ToolResult:
        """Run Protenix with multiple seeds, parallelised across GPUs."""
        devices = devices or []
        seeds = list(range(42, 42 + nstruct))

        # Single device (or no device): one invocation with all seeds
        if len(devices) <= 1:
            work_dir = ensure_dir(pred_dir / "protenix")
            tool = ProtenixTool(tool_config, work_dir, self.logs_dir)

            if not tool.check():
                logger.warning("Protenix not available, skipping")
                return ToolResult(success=False, error_message="Protenix not available")

            logger.info("Stage 3: Running Protenix (seeds %s)...", seeds)
            kwargs: dict[str, Any] = {
                "fasta_path": fasta_path,
                "seeds": seeds,
            }
            if msa_path:
                kwargs["msa_path"] = msa_path
            if devices:
                kwargs["device"] = devices[0]
            return tool.run(**kwargs)

        # Multiple devices: split seeds across GPUs, run in parallel
        n_gpus = len(devices)
        seed_groups: list[list[int]] = [[] for _ in range(n_gpus)]
        for i, seed in enumerate(seeds):
            seed_groups[i % n_gpus].append(seed)

        logger.info(
            "Stage 3: Running Protenix (%d seeds across %d GPUs)...",
            len(seeds), n_gpus,
        )

        def _run_group(gpu_idx: int, group_seeds: list[int]) -> ToolResult:
            device = devices[gpu_idx]
            work_dir = ensure_dir(pred_dir / "protenix" / f"gpu_{gpu_idx}")
            tool = ProtenixTool(tool_config, work_dir, self.logs_dir)
            if not tool.check():
                return ToolResult(
                    success=False,
                    error_message="Protenix not available",
                )
            logger.info(
                "Running Protenix (seeds %s) on %s...", group_seeds, device,
            )
            kw: dict[str, Any] = {
                "fasta_path": fasta_path,
                "seeds": group_seeds,
                "device": device,
            }
            if msa_path:
                kw["msa_path"] = msa_path
            return tool.run(**kw)

        group_results: list[ToolResult] = []
        with ThreadPoolExecutor(max_workers=n_gpus) as executor:
            futures = {
                executor.submit(_run_group, idx, grp): idx
                for idx, grp in enumerate(seed_groups)
                if grp  # skip empty groups
            }
            for future in as_completed(futures):
                group_results.append(future.result())

        # Merge results from all GPU groups
        all_structures: list[Path] = []
        merged_metrics: dict[str, Any] = {}
        total_time = 0.0
        any_success = False

        for r in group_results:
            total_time += r.runtime_seconds
            if r.success:
                any_success = True
                pdbs = r.output_files.get("all_pdbs", [])
                if pdbs:
                    all_structures.extend(pdbs)
                elif r.output_files.get("pdb"):
                    all_structures.append(r.output_files["pdb"])
                merged_metrics.update(r.metrics)

        if not any_success:
            return ToolResult(
                success=False,
                error_message="All Protenix GPU runs failed",
                runtime_seconds=total_time,
            )

        return ToolResult(
            success=True,
            output_files={
                "pdb": all_structures[0] if all_structures else None,
                "all_pdbs": all_structures,
            },
            metrics=merged_metrics,
            runtime_seconds=total_time,
        )

    def _run_simrna(
        self,
        fasta_path: Path,
        secondary_structure: str,
        tool_config: Any,
        pred_dir: Path,
        nstruct: int,
    ) -> ToolResult:
        """Run SimRNA with nstruct override for replicas/clustering."""
        work_dir = ensure_dir(pred_dir / "simrna")
        tool = SimRNATool(tool_config, work_dir, self.logs_dir)

        if not tool.check():
            logger.warning("SimRNA not available, skipping")
            return ToolResult(success=False, error_message="SimRNA not available")

        logger.info("Stage 3: Running SimRNA (nstruct=%d)...", nstruct)
        return tool.run(
            fasta_path=fasta_path,
            secondary_structure=secondary_structure,
            nstruct=nstruct if nstruct > 1 else None,
        )

    def _run_stage3b_clustering(
        self,
        predictor_results: dict[str, ToolResult],
        cutoff: float,
    ) -> Any:
        """Stage 3b: Cluster all structures across predictors by RMSD."""
        from rnapipey.ensemble import cluster_structures, EnsembleResult

        stage_key = "stage3b_clustering"
        if self._completed(stage_key):
            logger.info("Stage 3b (Clustering) already completed, skipping")
            return None

        all_pdbs: list[Path] = []
        all_labels: list[str] = []

        for pred_name, result in predictor_results.items():
            if not result.success:
                continue
            pdbs = result.output_files.get("all_pdbs", [])
            if not pdbs:
                pdb = result.output_files.get("pdb")
                if pdb and isinstance(pdb, Path) and pdb.exists():
                    pdbs = [pdb]
            for p in pdbs:
                if isinstance(p, Path) and p.exists():
                    all_pdbs.append(p)
                    all_labels.append(pred_name)

        if len(all_pdbs) < 2:
            logger.info("Only %d structure(s), skipping clustering", len(all_pdbs))
            self._mark(stage_key, "skipped")
            return None

        logger.info(
            "Stage 3b: Clustering %d structures (cutoff=%.1f A)...",
            len(all_pdbs), cutoff,
        )

        try:
            ensemble_result = cluster_structures(all_pdbs, all_labels, cutoff)
            self._mark(stage_key, "completed")
            return ensemble_result
        except Exception as e:
            logger.warning("Clustering failed: %s", e)
            self._mark(stage_key, "failed")
            return None

    def _run_stage4(
        self,
        predictor_results: dict[str, ToolResult],
        ensemble_result: Any = None,
    ) -> ToolResult | None:
        """Stage 4: Model evaluation and scoring."""
        if self._completed("stage4"):
            logger.info("Stage 4 (Scoring) already completed, skipping")
            scores_file = self.output_dir / "04_scoring" / "rnadvisor_scores.json"
            if scores_file.exists():
                data = json.loads(scores_file.read_text())
                return ToolResult(success=True, metrics={"scores": data})
            return None

        # Collect PDB files to score
        pdb_files: list[Path] = []

        if ensemble_result is not None:
            # Score only cluster representatives to save time
            for cluster in ensemble_result.clusters:
                rep = cluster.representative
                if rep.exists():
                    pdb_files.append(rep)
            logger.info(
                "Ensemble mode: scoring %d cluster representatives (not all %d structures)",
                len(pdb_files), len(ensemble_result.pdb_files),
            )
        else:
            # Score all PDB files from predictors
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
        ensemble_result: Any = None,
    ) -> None:
        """Stage 5: Generate report and visualization scripts."""
        vis_dir = ensure_dir(self.output_dir / "05_visualization")
        fasta_path = self.output_dir / "input" / "query.fasta"

        logger.info("Stage 5: Generating report and visualization scripts...")
        generate_report(
            vis_dir, fasta_path, predictor_results, scoring_result, ss_result,
            ensemble_result=ensemble_result,
        )
        generate_pymol_scripts(
            vis_dir, predictor_results, scoring_result,
            ensemble_result=ensemble_result,
        )
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
