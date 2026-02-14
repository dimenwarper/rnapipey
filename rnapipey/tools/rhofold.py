"""RhoFold+ wrapper for deep learning-based RNA 3D structure prediction."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from rnapipey.config import RhoFoldConfig
from rnapipey.tools.base import BaseTool, ToolResult

logger = logging.getLogger("rnapipey")


class RhoFoldTool(BaseTool):
    config: RhoFoldConfig

    @property
    def name(self) -> str:
        return "rhofold"

    def check(self) -> bool:
        script = self.config.script
        if not script:
            return False
        return Path(script).exists()

    def run(self, **kwargs: Any) -> ToolResult:
        fasta_path: Path = kwargs["fasta_path"]
        msa_path: Path | None = kwargs.get("msa_path")
        seed: int | None = kwargs.get("seed")
        device: str | None = kwargs.get("device")

        cmd = [
            "python", self.config.script,
            "--input_fas", str(fasta_path),
            "--output_dir", str(self.work_dir),
            "--single_seq_pred", "True",  # always enable single-seq mode as fallback
        ]
        if self.config.model_dir:
            cmd.extend(["--ckpt", str(self.config.model_dir)])

        # Device: kwarg overrides config
        effective_device = device if device else self.config.device
        if effective_device:
            cmd.extend(["--device", effective_device])

        if msa_path and msa_path.exists():
            cmd.extend(["--input_a3m", str(msa_path)])

        # Set PYTHONHASHSEED for reproducibility when seed is specified
        env: dict[str, str] | None = None
        if seed is not None:
            env = {"PYTHONHASHSEED": str(seed)}

        result = self._run_cmd(cmd, env=env)
        if result.returncode != 0:
            return ToolResult(
                success=False,
                error_message=f"RhoFold+ failed: {result.stderr[:300]}",
                runtime_seconds=result.runtime_seconds,
            )

        # Find output PDB
        pdb_files = list(self.work_dir.glob("*.pdb"))
        pdb = pdb_files[0] if pdb_files else None

        # RhoFold+ also outputs confidence scores and secondary structure
        ss_files = list(self.work_dir.glob("*.ct"))
        npz_files = list(self.work_dir.glob("*.npz"))

        # Extract pLDDT from .npz output
        metrics: dict[str, Any] = {}
        if npz_files:
            try:
                import numpy as np
                data = np.load(npz_files[0])
                if "plddt" in data:
                    metrics["plddt_mean"] = float(np.mean(data["plddt"]))
                    metrics["plddt_per_residue"] = data["plddt"].tolist()
            except Exception:
                logger.debug("Could not extract pLDDT from %s", npz_files[0])

        return ToolResult(
            success=pdb is not None,
            output_files={
                "pdb": pdb,
                "all_pdbs": pdb_files,
                "secondary_structure": ss_files[0] if ss_files else None,
                "distogram": npz_files[0] if npz_files else None,
            },
            metrics=metrics,
            runtime_seconds=result.runtime_seconds,
        )
