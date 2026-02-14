"""RhoFold+ wrapper for deep learning-based RNA 3D structure prediction."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from rnapipe.config import RhoFoldConfig
from rnapipe.tools.base import BaseTool, ToolResult

logger = logging.getLogger("rnapipe")


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

        cmd = [
            "python", self.config.script,
            "--input_fas", str(fasta_path),
            "--output_dir", str(self.work_dir),
            "--single_seq_pred",  # always enable single-seq mode as fallback
        ]
        if self.config.model_dir:
            cmd.extend(["--ckpt", str(self.config.model_dir)])
        if self.config.device:
            cmd.extend(["--device", self.config.device])
        if msa_path and msa_path.exists():
            cmd.extend(["--input_a3m", str(msa_path)])

        result = self._run_cmd(cmd)
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

        return ToolResult(
            success=pdb is not None,
            output_files={
                "pdb": pdb,
                "secondary_structure": ss_files[0] if ss_files else None,
                "distogram": npz_files[0] if npz_files else None,
            },
            metrics={},
            runtime_seconds=result.runtime_seconds,
        )
