"""SPOT-RNA wrapper for pseudoknot-aware secondary structure prediction."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from rnapipe.config import ToolsConfig
from rnapipe.tools.base import BaseTool, ToolResult
from rnapipe.utils import which

logger = logging.getLogger("rnapipe")


class SPOTRNATool(BaseTool):
    config: ToolsConfig

    @property
    def name(self) -> str:
        return "spotrna"

    def check(self) -> bool:
        spotrna = self.config.spotrna
        if not spotrna:
            return False
        return Path(spotrna).exists() or which(spotrna) is not None

    def run(self, **kwargs: Any) -> ToolResult:
        fasta_path: Path = kwargs["fasta_path"]

        cmd = [
            "python", self.config.spotrna,
            "--inputs", str(fasta_path),
            "--outputs", str(self.work_dir),
        ]
        result = self._run_cmd(cmd)
        if result.returncode != 0:
            return ToolResult(
                success=False,
                error_message=f"SPOT-RNA failed: {result.stderr[:300]}",
                runtime_seconds=result.runtime_seconds,
            )

        # SPOT-RNA outputs .ct and .bpseq files
        ct_files = list(self.work_dir.glob("*.ct"))
        bpseq_files = list(self.work_dir.glob("*.bpseq"))

        dot_bracket = ""
        if bpseq_files:
            dot_bracket = self._bpseq_to_dotbracket(bpseq_files[0])

        dot_file = self.work_dir / "spotrna.dot"
        if dot_bracket:
            dot_file.write_text(dot_bracket + "\n")

        return ToolResult(
            success=True,
            output_files={
                "dot": dot_file if dot_bracket else None,
                "ct": ct_files[0] if ct_files else None,
                "bpseq": bpseq_files[0] if bpseq_files else None,
            },
            metrics={"dot_bracket": dot_bracket},
            runtime_seconds=result.runtime_seconds,
        )

    def _bpseq_to_dotbracket(self, bpseq_path: Path) -> str:
        """Convert .bpseq format to dot-bracket notation (handles pseudoknots)."""
        pairs: dict[int, int] = {}
        length = 0
        for line in bpseq_path.read_text().splitlines():
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            idx, _base, partner = int(parts[0]), parts[1], int(parts[2])
            length = max(length, idx)
            if partner > 0 and idx < partner:
                pairs[idx] = partner

        # Simple assignment: first layer (), pseudoknots []
        structure = list("." * length)
        used = set()
        # Sort pairs, assign non-crossing ones to () first
        sorted_pairs = sorted(pairs.items())
        layer1 = []
        layer2 = []
        for i, j in sorted_pairs:
            crossing = any(i < k < j < l or k < i < l < j for k, l in layer1)
            if not crossing:
                layer1.append((i, j))
            else:
                layer2.append((i, j))

        for i, j in layer1:
            structure[i - 1] = "("
            structure[j - 1] = ")"
        for i, j in layer2:
            structure[i - 1] = "["
            structure[j - 1] = "]"

        return "".join(structure)
