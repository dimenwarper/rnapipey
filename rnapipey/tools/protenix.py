"""Protenix (open-source AlphaFold3) wrapper for RNA 3D structure prediction."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

from rnapipey.config import ProtenixConfig
from rnapipey.tools.base import BaseTool, ToolResult
from rnapipey.utils import read_fasta, which

logger = logging.getLogger("rnapipey")


class ProtenixTool(BaseTool):
    config: ProtenixConfig

    @property
    def name(self) -> str:
        return "protenix"

    def check(self) -> bool:
        try:
            import protenix  # noqa: F401
            return True
        except ImportError:
            return False

    def run(self, **kwargs: Any) -> ToolResult:
        fasta_path: Path = kwargs["fasta_path"]
        msa_path: Path | None = kwargs.get("msa_path")

        records = read_fasta(fasta_path)
        if not records:
            return ToolResult(success=False, error_message="Empty FASTA file")

        # Build Protenix input JSON
        input_json = self._build_input_json(records[0].sequence, records[0].id)
        json_path = self.work_dir / "input.json"
        json_path.write_text(json.dumps(input_json, indent=2))

        import sys

        cmd = [
            sys.executable, "-m", "protenix.predict",
            "--input", str(json_path),
            "--output_dir", str(self.work_dir),
        ]
        if self.config.model:
            cmd.extend(["--model_dir", str(self.config.model)])

        result = self._run_cmd(cmd, timeout=86400)
        if result.returncode != 0:
            return ToolResult(
                success=False,
                error_message=f"Protenix failed: {result.stderr[:300]}",
                runtime_seconds=result.runtime_seconds,
            )

        # Find output structure files (CIF or PDB)
        cif_files = list(self.work_dir.rglob("*.cif"))
        pdb_files = list(self.work_dir.rglob("*.pdb"))
        structure = cif_files[0] if cif_files else (pdb_files[0] if pdb_files else None)

        # Find confidence scores
        json_files = list(self.work_dir.rglob("*confidence*.json"))
        confidence = {}
        if json_files:
            try:
                confidence = json.loads(json_files[0].read_text())
            except (json.JSONDecodeError, OSError):
                pass

        return ToolResult(
            success=structure is not None,
            output_files={
                "pdb": structure,
                "input_json": json_path,
            },
            metrics=confidence,
            runtime_seconds=result.runtime_seconds,
        )

    def _build_input_json(self, sequence: str, name: str) -> dict:
        """Build Protenix inference input JSON for a single RNA chain."""
        return {
            "name": f"rnapipey_{name}",
            "modelSeeds": [42],
            "sequences": [
                {
                    "rnaSequence": {
                        "sequence": sequence,
                        "count": 1,
                    }
                }
            ],
        }
