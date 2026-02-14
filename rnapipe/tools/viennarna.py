"""ViennaRNA (RNAfold) wrapper for secondary structure prediction."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from rnapipe.config import ToolsConfig
from rnapipe.tools.base import BaseTool, ToolResult
from rnapipe.utils import read_fasta, which

logger = logging.getLogger("rnapipe")


class ViennaRNATool(BaseTool):
    config: ToolsConfig

    @property
    def name(self) -> str:
        return "rnafold"

    def check(self) -> bool:
        return which(self.config.rnafold) is not None

    def run(self, **kwargs: Any) -> ToolResult:
        fasta_path: Path = kwargs["fasta_path"]
        records = read_fasta(fasta_path)
        if not records:
            return ToolResult(success=False, error_message="Empty FASTA file")

        seq = records[0].sequence
        dot_file = self.work_dir / "rnafold.dot"

        cmd = [
            self.config.rnafold,
            "--noPS",  # skip PostScript output
            "-p",      # compute partition function + base pair probabilities
            "-i", str(fasta_path),
        ]
        result = self._run_cmd(cmd)
        if result.returncode != 0:
            return ToolResult(
                success=False,
                error_message=f"RNAfold failed: {result.stderr[:300]}",
                runtime_seconds=result.runtime_seconds,
            )

        # Parse RNAfold output
        dot_bracket, mfe = self._parse_output(result.stdout)
        dot_file.write_text(f">{records[0].header}\n{seq}\n{dot_bracket} ({mfe:.2f})\n")

        # RNAfold writes dp.ps for base pair probabilities in cwd
        bpp_ps = self.work_dir / f"{records[0].id}_dp.ps"
        if not bpp_ps.exists():
            # sometimes it uses the sequence name differently
            candidates = list(self.work_dir.glob("*_dp.ps"))
            bpp_ps = candidates[0] if candidates else None

        return ToolResult(
            success=True,
            output_files={"dot": dot_file, "bpp_ps": bpp_ps},
            metrics={"dot_bracket": dot_bracket, "mfe": mfe, "length": len(seq)},
            runtime_seconds=result.runtime_seconds,
        )

    def _parse_output(self, stdout: str) -> tuple[str, float]:
        """Extract dot-bracket and MFE from RNAfold stdout."""
        lines = [l.strip() for l in stdout.strip().splitlines() if l.strip()]
        # RNAfold output format:
        # >header
        # SEQUENCE
        # (((....))) (-12.30)
        for line in lines:
            if "(" in line and ")" in line and " (" in line:
                # Last parenthesized value is MFE
                parts = line.rsplit(" (", 1)
                if len(parts) == 2:
                    dot_bracket = parts[0].strip()
                    mfe_str = parts[1].rstrip(")")
                    try:
                        mfe = float(mfe_str)
                        return dot_bracket, mfe
                    except ValueError:
                        continue
        logger.warning("Could not parse RNAfold output, returning empty structure")
        return "", 0.0
