"""Shared utilities: FASTA I/O, subprocess helpers, logging."""

from __future__ import annotations

import logging
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger("rnapipey")


def setup_logging(log_file: Path | None = None, verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    fmt = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
    handlers: list[logging.Handler] = [logging.StreamHandler()]
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    logging.basicConfig(level=level, format=fmt, handlers=handlers)


@dataclass
class FastaRecord:
    header: str
    sequence: str

    @property
    def id(self) -> str:
        return self.header.split()[0]


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header = ""
    seq_lines: list[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                records.append(FastaRecord(header, "".join(seq_lines)))
            header = line[1:].strip()
            seq_lines = []
        else:
            seq_lines.append(line.upper())
    if header:
        records.append(FastaRecord(header, "".join(seq_lines)))
    return records


def write_fasta(records: list[FastaRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for rec in records:
        lines.append(f">{rec.header}")
        # wrap at 80 chars
        seq = rec.sequence
        for i in range(0, len(seq), 80):
            lines.append(seq[i : i + 80])
    path.write_text("\n".join(lines) + "\n")


@dataclass
class CmdResult:
    returncode: int
    stdout: str
    stderr: str
    runtime_seconds: float


def run_cmd(
    cmd: list[str],
    work_dir: Path | None = None,
    timeout: int = 86400,
    stdout_file: Path | None = None,
    stderr_file: Path | None = None,
) -> CmdResult:
    """Run a subprocess with logging and timing."""
    logger.info("Running: %s", " ".join(cmd))
    start = time.monotonic()
    proc = subprocess.run(
        cmd,
        cwd=work_dir,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    elapsed = time.monotonic() - start
    logger.debug("Finished in %.1fs (exit %d)", elapsed, proc.returncode)
    if stdout_file:
        stdout_file.parent.mkdir(parents=True, exist_ok=True)
        stdout_file.write_text(proc.stdout)
    if stderr_file:
        stderr_file.parent.mkdir(parents=True, exist_ok=True)
        stderr_file.write_text(proc.stderr)
    if proc.returncode != 0:
        logger.warning("Command failed (exit %d): %s", proc.returncode, proc.stderr[:500])
    return CmdResult(
        returncode=proc.returncode,
        stdout=proc.stdout,
        stderr=proc.stderr,
        runtime_seconds=elapsed,
    )


def which(name: str) -> str | None:
    return shutil.which(name)


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path
