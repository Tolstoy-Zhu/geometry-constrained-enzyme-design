"""Configuration utilities for AquaEsterNet-StructMine.

This file intentionally uses English filenames and UTF-8 encoding so that the
project can be opened correctly in PyCharm on Windows, macOS, Linux, or WSL2.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict

import yaml

from .models import ProjectConfig


def project_root() -> Path:
    """Return project root based on this source file location."""
    return Path(__file__).resolve().parents[2]


def load_config(config_path: str | Path) -> Dict[str, Any]:
    """Load and validate a YAML configuration file."""
    path = Path(config_path)
    if not path.is_absolute():
        path = project_root() / path

    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with path.open("r", encoding="utf-8") as f:
        raw = yaml.safe_load(f)
    if not isinstance(raw, dict):
        raise ValueError(f"Config root must be a mapping: {path}")
    return ProjectConfig.model_validate(raw).model_dump(mode="python", exclude_none=False)


def resolve_path(path_like: str | Path) -> Path:
    """Resolve a path relative to project root if it is not absolute."""
    path = Path(path_like)
    if path.is_absolute():
        return path
    return project_root() / path
