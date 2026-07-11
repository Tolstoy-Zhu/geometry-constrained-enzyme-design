"""Search database configuration for AquaEsterNet-StructMine."""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from .config import resolve_path
from .foldseek_runner import FOLDSEEK_FIELDS

PROSTT5_FIELDS: List[str] = [
    "query",
    "target",
    "evalue",
    "bits",
    "prob",
    "fident",
    "alnlen",
    "qlen",
    "tlen",
    "qcov",
    "tcov",
]


@dataclass
class SearchDatabase:
    """One Foldseek-searchable data source."""

    name: str
    enabled: bool = True
    kind: str = "foldseek"
    target_db: Optional[Path] = None
    annotation_source: str = "auto"
    foldseek_database: str = ""
    archive_url: str = ""
    version_url: str = ""
    archive_sha256: str = ""
    archive_size: Optional[int] = None
    archive_prefix: str = ""
    input_fasta: Optional[Path] = None
    structure_dir: Optional[Path] = None
    result_tsv: Optional[Path] = None
    search_engine: str = "foldseek"
    auto_prepare: bool = False
    prostt5_model: Optional[Path] = None
    use_gpu: bool = False
    cluster_search: Optional[bool] = None
    search_extra_args: List[str] = field(default_factory=list)
    prepare_extra_args: List[str] = field(default_factory=list)
    format_fields: List[str] = field(default_factory=lambda: list(FOLDSEEK_FIELDS))

    @property
    def safe_name(self) -> str:
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", self.name.strip())
        return safe or "database"


def _as_bool(value: Any, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _as_list(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [x.strip() for x in value.split() if x.strip()]
    return [str(x) for x in value]


def _as_fields(value: Any, kind: str) -> List[str]:
    if value:
        if isinstance(value, str):
            return [x.strip() for x in value.split(",") if x.strip()]
        return [str(x).strip() for x in value if str(x).strip()]
    if kind in {"ncbi_fasta", "prostt5_fasta"}:
        return list(PROSTT5_FIELDS)
    return list(FOLDSEEK_FIELDS)


def _as_optional_bool(value: Any) -> Optional[bool]:
    if value is None:
        return None
    return _as_bool(value)


def _as_optional_int(value: Any) -> Optional[int]:
    if value in (None, ""):
        return None
    result = int(value)
    if result <= 0:
        raise ValueError("archive_size must be a positive integer")
    return result


def _resolve_optional_path(value: Any) -> Optional[Path]:
    if value in (None, ""):
        return None
    return resolve_path(value)


def _infer_annotation_source(entry: Dict[str, Any], kind: str) -> str:
    explicit = str(entry.get("annotation_source", "")).strip().lower()
    if explicit:
        return explicit

    text = " ".join(
        str(entry.get(k, ""))
        for k in ("name", "kind", "foldseek_database", "target_db", "input_fasta", "structure_dir")
    ).lower()
    if kind in {"alphafold", "afdb"} or "alphafold" in text or "afdb" in text:
        if "swiss" in text:
            return "swissprot"
        return "alphafold"
    if kind in {"swissprot", "uniprot"} or "swiss" in text:
        return "swissprot"
    if kind in {"ncbi", "ncbi_fasta"} or "ncbi" in text:
        return "ncbi"
    if "pdb" in text:
        return "pdb"
    return "auto"


def _iter_database_entries(cfg: Dict[str, Any]) -> Iterable[Dict[str, Any]]:
    entries = cfg.get("search_databases")
    if isinstance(entries, dict):
        for name, entry in entries.items():
            if isinstance(entry, dict):
                item = dict(entry)
                item.setdefault("name", name)
                yield item
        return

    if isinstance(entries, list):
        for entry in entries:
            if isinstance(entry, dict):
                yield dict(entry)
        return

    foldseek_cfg = cfg.get("foldseek", {})
    yield {
        "name": "pdb",
        "enabled": True,
        "kind": "foldseek",
        "target_db": foldseek_cfg.get("target_db", "databases/pdb"),
        "annotation_source": "pdb",
    }


def load_search_databases(cfg: Dict[str, Any]) -> List[SearchDatabase]:
    """Load enabled search database definitions from config."""
    databases: List[SearchDatabase] = []
    for raw in _iter_database_entries(cfg):
        name = str(raw.get("name", "")).strip()
        if not name:
            raise ValueError("Each search database entry must have a name")

        enabled = _as_bool(raw.get("enabled", True), default=True)
        if not enabled:
            continue

        kind = str(raw.get("kind", raw.get("type", "foldseek"))).strip().lower()
        annotation_source = _infer_annotation_source(raw, kind)

        databases.append(
            SearchDatabase(
                name=name,
                enabled=enabled,
                kind=kind,
                target_db=_resolve_optional_path(raw.get("target_db")),
                annotation_source=annotation_source,
                foldseek_database=str(raw.get("foldseek_database", "")).strip(),
                archive_url=str(raw.get("archive_url", "")).strip(),
                version_url=str(raw.get("version_url", "")).strip(),
                archive_sha256=str(raw.get("archive_sha256", "")).strip().lower(),
                archive_size=_as_optional_int(raw.get("archive_size")),
                archive_prefix=str(raw.get("archive_prefix", "")).strip(),
                input_fasta=_resolve_optional_path(raw.get("input_fasta")),
                structure_dir=_resolve_optional_path(raw.get("structure_dir")),
                result_tsv=_resolve_optional_path(raw.get("result_tsv")),
                search_engine=str(raw.get("search_engine", "foldseek")).strip().lower(),
                auto_prepare=_as_bool(raw.get("auto_prepare", False), default=False),
                prostt5_model=_resolve_optional_path(raw.get("prostt5_model")),
                use_gpu=_as_bool(raw.get("use_gpu", False), default=False),
                cluster_search=_as_optional_bool(raw.get("cluster_search")),
                search_extra_args=_as_list(raw.get("search_extra_args")),
                prepare_extra_args=_as_list(raw.get("prepare_extra_args")),
                format_fields=_as_fields(raw.get("format_fields"), kind),
            )
        )

    if not databases:
        raise ValueError("No enabled search databases were configured")
    return databases
