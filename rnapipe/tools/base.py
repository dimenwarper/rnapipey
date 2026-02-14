"""Base class for external tool wrappers."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from rnapipe.utils import CmdResult, run_cmd


@dataclass
class ToolResult:
    """Uniform result from any external tool."""

    success: bool
    output_files: dict[str, Path | None] = field(default_factory=dict)
    metrics: dict[str, Any] = field(default_factory=dict)
    error_message: str = ""
    runtime_seconds: float = 0.0


class BaseTool(ABC):
    """Abstract base for wrapping an external bioinformatics tool."""

    def __init__(self, config: Any, work_dir: Path, logs_dir: Path | None = None):
        self.config = config
        self.work_dir = work_dir
        self.logs_dir = logs_dir or work_dir
        self.work_dir.mkdir(parents=True, exist_ok=True)

    @property
    @abstractmethod
    def name(self) -> str: ...

    @abstractmethod
    def check(self) -> bool:
        """Return True if the tool is installed and accessible."""

    @abstractmethod
    def run(self, **kwargs: Any) -> ToolResult:
        """Execute the tool. Inputs vary by tool."""

    def _run_cmd(self, cmd: list[str], **kwargs: Any) -> CmdResult:
        """Run a subprocess with stdout/stderr captured to log files."""
        return run_cmd(
            cmd,
            work_dir=kwargs.get("work_dir", self.work_dir),
            stdout_file=self.logs_dir / f"{self.name}.stdout",
            stderr_file=self.logs_dir / f"{self.name}.stderr",
            **{k: v for k, v in kwargs.items() if k != "work_dir"},
        )
