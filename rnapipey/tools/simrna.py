"""SimRNA wrapper for physics-based RNA 3D structure prediction."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from rnapipey.config import SimRNAConfig
from rnapipey.tools.base import BaseTool, ToolResult
from rnapipey.utils import read_fasta, which

logger = logging.getLogger("rnapipey")


class SimRNATool(BaseTool):
    config: SimRNAConfig

    @property
    def name(self) -> str:
        return "simrna"

    def check(self) -> bool:
        binary = self.config.binary
        if not binary:
            return False
        return Path(binary).exists() or which(binary) is not None

    def run(self, **kwargs: Any) -> ToolResult:
        fasta_path: Path = kwargs["fasta_path"]
        secondary_structure: str = kwargs.get("secondary_structure", "")

        records = read_fasta(fasta_path)
        if not records:
            return ToolResult(success=False, error_message="Empty FASTA file")

        seq = records[0].sequence

        # Write SimRNA input file: sequence on first line, secondary structure on second
        input_file = self.work_dir / "input.seq"
        ss = secondary_structure if secondary_structure else "." * len(seq)
        input_file.write_text(f"{seq}\n{ss}\n")

        # Write restraints file from secondary structure
        restraints_file = None
        if secondary_structure:
            restraints_file = self._generate_restraints(secondary_structure)

        # Run SimRNA
        cmd = [
            self.config.binary,
            "-s", str(input_file),
            "-o", str(self.work_dir / "simrna_run"),
            "-n", str(self.config.steps),
            "-R", str(self.config.replicas),
        ]
        if self.config.data_dir:
            cmd.extend(["-E", str(self.config.data_dir)])
        if restraints_file:
            cmd.extend(["-r", str(restraints_file)])

        result = self._run_cmd(cmd, timeout=86400)  # can take hours
        if result.returncode != 0:
            return ToolResult(
                success=False,
                error_message=f"SimRNA failed: {result.stderr[:300]}",
                runtime_seconds=result.runtime_seconds,
            )

        # Process trajectory: extract lowest energy frames, cluster
        trafl_files = list(self.work_dir.glob("*.trafl"))
        pdb_files = self._cluster_and_extract(trafl_files)

        return ToolResult(
            success=len(pdb_files) > 0,
            output_files={
                "pdb": pdb_files[0] if pdb_files else None,
                "trafl": trafl_files[0] if trafl_files else None,
                "all_pdbs": pdb_files,
            },
            metrics={"n_clusters": len(pdb_files)},
            runtime_seconds=result.runtime_seconds,
        )

    def _generate_restraints(self, ss: str) -> Path:
        """Convert dot-bracket to SimRNA distance restraints."""
        restraints_file = self.work_dir / "restraints.txt"
        lines = []
        stack: list[int] = []
        for i, ch in enumerate(ss):
            if ch == "(":
                stack.append(i)
            elif ch == ")" and stack:
                j = stack.pop()
                # SimRNA restraint format: DIST chain res1 atom1 chain res2 atom2 dist_min dist_max weight
                # Use N1/N3 atoms for base-pairing distance restraints
                lines.append(
                    f"DIST A {j+1} N1 A {i+1} N3 5.0 10.0 1.0"
                )
        restraints_file.write_text("\n".join(lines) + "\n")
        return restraints_file

    def _cluster_and_extract(self, trafl_files: list[Path]) -> list[Path]:
        """Extract top N structures from SimRNA trajectory files.

        SimRNA provides clustering scripts (SimRNA_trafl2pdbs). If not available,
        we extract the lowest-energy frames directly from the .trafl file.
        """
        pdb_files: list[Path] = []
        if not trafl_files:
            return pdb_files

        trafl = trafl_files[0]
        # Check if SimRNA_trafl2pdbs is available
        trafl2pdbs = which("SimRNA_trafl2pdbs")
        if trafl2pdbs:
            # Use the official clustering tool
            result = self._run_cmd([
                trafl2pdbs,
                str(trafl),
                str(self.config.clustering_top_n),
            ])
            pdb_files = sorted(self.work_dir.glob("simrna_run*.pdb"))
        else:
            # Fallback: extract low-energy frames manually
            # trafl format: each frame is multi-line, energy on the header line
            logger.warning(
                "SimRNA_trafl2pdbs not found. "
                "Trajectory files saved but PDB extraction skipped. "
                "Run SimRNA_trafl2pdbs manually on: %s", trafl
            )

        return pdb_files[: self.config.clustering_top_n]
