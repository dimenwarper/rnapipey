"""Configuration loading and validation."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger("rnapipey")

DEFAULT_CONFIG = Path(__file__).parent.parent / "configs" / "default.yaml"


@dataclass
class RhoFoldConfig:
    script: str = ""
    model_dir: str = ""
    device: str = "cuda:0"


@dataclass
class SimRNAConfig:
    binary: str = ""
    data_dir: str = ""
    replicas: int = 10
    steps: int = 10_000_000
    clustering_top_n: int = 5


@dataclass
class ProtenixConfig:
    binary: str = "protenix"
    model: str = ""
    data_dir: str = "~/.protenix"


@dataclass
class RNAdvisorConfig:
    docker: bool = True
    metrics: list[str] = field(default_factory=lambda: ["rsRNASP", "DFIRE", "RASP", "MCQ"])


@dataclass
class ToolsConfig:
    cmscan: str = "cmscan"
    cmfetch: str = "cmfetch"
    cmalign: str = "cmalign"
    rfam_cm: str = ""
    rfam_clanin: str = ""
    rnafold: str = "RNAfold"
    spotrna: str = ""
    rhofold: RhoFoldConfig = field(default_factory=RhoFoldConfig)
    simrna: SimRNAConfig = field(default_factory=SimRNAConfig)
    protenix: ProtenixConfig = field(default_factory=ProtenixConfig)
    rnadvisor: RNAdvisorConfig = field(default_factory=RNAdvisorConfig)


@dataclass
class EnsembleConfig:
    nstruct: int = 1
    cluster: bool = True
    cluster_cutoff: float = 5.0  # RMSD in Angstroms


@dataclass
class PipelineConfig:
    tools: ToolsConfig = field(default_factory=ToolsConfig)
    ensemble: EnsembleConfig = field(default_factory=EnsembleConfig)


def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge override into base."""
    merged = base.copy()
    for key, val in override.items():
        if key in merged and isinstance(merged[key], dict) and isinstance(val, dict):
            merged[key] = _deep_merge(merged[key], val)
        else:
            merged[key] = val
    return merged


def _dict_to_dataclass(cls: type, data: dict[str, Any]) -> Any:
    """Recursively convert a dict to a nested dataclass."""
    field_types = {f.name: f.type for f in cls.__dataclass_fields__.values()}
    kwargs = {}
    for key, val in data.items():
        if key not in field_types:
            continue
        ft = field_types[key]
        # Resolve string type annotations to actual classes
        if isinstance(ft, str):
            ft = eval(ft, globals(), {
                "RhoFoldConfig": RhoFoldConfig,
                "SimRNAConfig": SimRNAConfig,
                "ProtenixConfig": ProtenixConfig,
                "RNAdvisorConfig": RNAdvisorConfig,
                "ToolsConfig": ToolsConfig,
                "EnsembleConfig": EnsembleConfig,
            })
        if isinstance(val, dict) and hasattr(ft, "__dataclass_fields__"):
            kwargs[key] = _dict_to_dataclass(ft, val)
        else:
            kwargs[key] = val
    return cls(**kwargs)


def load_config(path: Path | None = None) -> PipelineConfig:
    """Load config from YAML, merging with defaults."""
    base_data: dict[str, Any] = {}
    if DEFAULT_CONFIG.exists():
        base_data = yaml.safe_load(DEFAULT_CONFIG.read_text()) or {}

    if path:
        user_data = yaml.safe_load(path.read_text()) or {}
        base_data = _deep_merge(base_data, user_data)

    return _dict_to_dataclass(PipelineConfig, base_data)
