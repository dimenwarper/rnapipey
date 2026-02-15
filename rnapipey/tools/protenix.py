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
        return which(self.config.binary) is not None

    def run(self, **kwargs: Any) -> ToolResult:
        fasta_path: Path = kwargs["fasta_path"]
        msa_path: Path | None = kwargs.get("msa_path")
        seeds: list[int] = kwargs.get("seeds", [42])
        device: str | None = kwargs.get("device")

        records = read_fasta(fasta_path)
        if not records:
            return ToolResult(success=False, error_message="Empty FASTA file")

        # Build Protenix input JSON with multiple seeds
        input_json = self._build_input_json(
            records[0].sequence, records[0].id, seeds=seeds
        )
        json_path = self.work_dir / "input.json"
        json_path.write_text(json.dumps(input_json, indent=2))

        cmd = [
            self.config.binary, "pred",
            "-i", str(json_path),
            "-o", str(self.work_dir),
        ]
        if self.config.model:
            cmd.extend(["-n", str(self.config.model)])

        # Set CUDA_VISIBLE_DEVICES to pin this invocation to a specific GPU
        env: dict[str, str] | None = None
        if device:
            # Extract GPU index from device string (e.g. "cuda:2" -> "2")
            gpu_id = device.split(":")[-1] if ":" in device else device
            env = {"CUDA_VISIBLE_DEVICES": gpu_id}

        result = self._run_cmd(cmd, timeout=86400, env=env)
        if result.returncode != 0:
            return ToolResult(
                success=False,
                error_message=f"Protenix failed: {result.stderr[:300]}",
                runtime_seconds=result.runtime_seconds,
            )

        # Find output structure files (CIF or PDB)
        cif_files = sorted(self.work_dir.rglob("*.cif"))
        pdb_files = sorted(self.work_dir.rglob("*.pdb"))
        all_structures = cif_files + pdb_files
        structure = all_structures[0] if all_structures else None

        # Find confidence scores
        json_files = list(self.work_dir.rglob("*confidence*.json"))
        confidence: dict[str, Any] = {}
        if json_files:
            try:
                confidence = json.loads(json_files[0].read_text())
            except (json.JSONDecodeError, OSError):
                pass

        return ToolResult(
            success=structure is not None,
            output_files={
                "pdb": structure,
                "all_pdbs": all_structures,
                "input_json": json_path,
            },
            metrics=confidence,
            runtime_seconds=result.runtime_seconds,
        )

    def _build_input_json(
        self, sequence: str, name: str, seeds: list[int] | None = None
    ) -> list[dict]:
        """Build Protenix inference input JSON for a single RNA chain."""
        return [
            {
                "name": f"rnapipey_{name}",
                "modelSeeds": seeds or [42],
                "sequences": [
                    {
                        "rnaSequence": {
                            "sequence": sequence,
                            "count": 1,
                        }
                    }
                ],
            }
        ]
