"""Infernal (cmscan) wrapper for Rfam family identification and MSA building."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from rnapipey.config import ToolsConfig
from rnapipey.tools.base import BaseTool, ToolResult
from rnapipey.utils import which

logger = logging.getLogger("rnapipey")


class InfernalTool(BaseTool):
    config: ToolsConfig

    @property
    def name(self) -> str:
        return "infernal"

    def check(self) -> bool:
        if not which(self.config.cmscan):
            return False
        if not self.config.rfam_cm or not Path(self.config.rfam_cm).exists():
            logger.warning("Rfam.cm not found at: %s", self.config.rfam_cm)
            return False
        return True

    def run(self, **kwargs: Any) -> ToolResult:
        fasta_path: Path = kwargs["fasta_path"]
        tblout = self.work_dir / "cmscan_tblout.txt"
        output = self.work_dir / "cmscan_output.txt"

        cmd = [
            self.config.cmscan,
            "--cut_ga",
            "--rfam",
            "--nohmmonly",
            "--fmt", "2",
            "--tblout", str(tblout),
            "-o", str(output),
        ]
        if self.config.rfam_clanin and Path(self.config.rfam_clanin).exists():
            cmd.extend(["--clanin", str(self.config.rfam_clanin)])
        cmd.extend([self.config.rfam_cm, str(fasta_path)])

        result = self._run_cmd(cmd)
        if result.returncode != 0:
            return ToolResult(
                success=False,
                error_message=f"cmscan failed: {result.stderr[:300]}",
                runtime_seconds=result.runtime_seconds,
            )

        family, evalue = self._parse_tblout(tblout)
        output_files: dict[str, Path | None] = {"tblout": tblout, "output": output}
        metrics: dict[str, Any] = {"rfam_family": family, "evalue": evalue}

        # If a family was found, try to build an MSA
        msa_path = None
        if family:
            msa_path = self._build_msa(family, fasta_path)
            output_files["msa"] = msa_path

        return ToolResult(
            success=True,
            output_files=output_files,
            metrics=metrics,
            runtime_seconds=result.runtime_seconds,
        )

    def _parse_tblout(self, tblout: Path) -> tuple[str | None, float | None]:
        """Parse cmscan tblout for top hit."""
        if not tblout.exists():
            return None, None
        for line in tblout.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) >= 16:
                family = fields[1]  # target name (Rfam accession)
                evalue = float(fields[15])
                logger.info("Rfam hit: %s (E-value: %e)", family, evalue)
                return family, evalue
        return None, None

    def _build_msa(self, family: str, fasta_path: Path) -> Path | None:
        """Fetch the CM for the family and align the query to build an MSA."""
        cm_out = self.work_dir / f"{family}.cm"
        msa_out = self.work_dir / "alignment.sto"

        # cmfetch: extract single family CM
        fetch_result = self._run_cmd([
            self.config.cmfetch,
            "-o", str(cm_out),
            self.config.rfam_cm,
            family,
        ])
        if fetch_result.returncode != 0:
            logger.warning("cmfetch failed for %s", family)
            return None

        # cmalign: align query to the family CM
        align_result = self._run_cmd([
            self.config.cmalign,
            "--outformat", "Stockholm",
            "-o", str(msa_out),
            str(cm_out),
            str(fasta_path),
        ])
        if align_result.returncode != 0:
            logger.warning("cmalign failed for %s", family)
            return None

        logger.info("Built MSA: %s", msa_out)
        return msa_out
