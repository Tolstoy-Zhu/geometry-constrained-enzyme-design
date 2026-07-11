"""Versioned run directories, manifests, and progress event streams."""
from __future__ import annotations

import hashlib
import platform
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import Bio
import numpy
import pandas

from .io_utils import ensure_dir
from .models import ProgressEvent, RunManifest, RunSummary

ProgressCallback = Callable[[ProgressEvent], None]


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def make_run_id(label: str = "run") -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "-", label).strip("-") or "run"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    return f"{stamp}_{safe}"


def mark_run_cancelled(run_dir: str | Path, reason: str = "Cancelled by user") -> bool:
    """Atomically mark an interrupted run and append a cancellation event."""
    run_dir = Path(run_dir).resolve()
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.exists():
        return False
    manifest = RunManifest.model_validate_json(manifest_path.read_text(encoding="utf-8"))
    if manifest.status in {"completed", "failed", "cancelled"}:
        return manifest.status == "cancelled"
    manifest.status = "cancelled"
    manifest.completed_at = datetime.now(timezone.utc)
    manifest.error = reason
    tmp = manifest_path.with_suffix(".json.tmp")
    tmp.write_text(manifest.model_dump_json(indent=2), encoding="utf-8")
    tmp.replace(manifest_path)
    event = ProgressEvent(stage="run", state="cancelled", message=reason)
    with (run_dir / "events.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(event.model_dump_json() + "\n")
    return True


class RunContext:
    """Persist enough state to inspect or resume a pipeline run."""

    def __init__(
        self,
        run_root: str | Path,
        config_path: str | Path,
        config_snapshot: dict[str, Any],
        run_id: str | None = None,
        label: str = "run",
        callback: ProgressCallback | None = None,
    ) -> None:
        self.run_id = run_id or make_run_id(label)
        self.run_dir = ensure_dir(Path(run_root) / self.run_id)
        self.events_path = self.run_dir / "events.jsonl"
        self.manifest_path = self.run_dir / "run_manifest.json"
        self.callback = callback
        config_path = Path(config_path).resolve()
        self.manifest = RunManifest(
            run_id=self.run_id,
            config_path=str(config_path),
            config_sha256=file_sha256(config_path),
            config_snapshot=config_snapshot,
            software_versions={
                "python": platform.python_version(),
                "platform": platform.platform(),
                "biopython": Bio.__version__,
                "numpy": numpy.__version__,
                "pandas": pandas.__version__,
            },
        )
        self._write_manifest()

    def _write_manifest(self) -> None:
        tmp = self.manifest_path.with_suffix(".json.tmp")
        tmp.write_text(self.manifest.model_dump_json(indent=2), encoding="utf-8")
        tmp.replace(self.manifest_path)

    def emit(
        self,
        stage: str,
        state: str,
        message: str,
        progress: float | None = None,
        artifact: str | Path | None = None,
    ) -> ProgressEvent:
        event = ProgressEvent(
            stage=stage,
            state=state,
            message=message,
            progress=progress,
            artifact=str(artifact) if artifact else None,
        )
        with self.events_path.open("a", encoding="utf-8") as handle:
            handle.write(event.model_dump_json() + "\n")
        if self.callback:
            self.callback(event)
        return event

    def start(self) -> None:
        self.manifest.status = "running"
        self.manifest.started_at = datetime.now(timezone.utc)
        self._write_manifest()
        self.emit("run", "started", "Pipeline run started", progress=0.0)

    def set_query(
        self,
        path: str | Path,
        pdb_id: str = "",
        chain: str = "",
        residue_range: str = "",
        source_url: str = "",
    ) -> None:
        path = Path(path).resolve()
        self.manifest.query_path = str(path)
        self.manifest.query_sha256 = file_sha256(path)
        self.manifest.query_pdb_id = pdb_id
        self.manifest.query_chain = chain
        self.manifest.query_residue_range = residue_range
        self.manifest.query_source_url = source_url
        self._write_manifest()

    def set_prediction_lineage(
        self,
        parent_prediction_id: str | None,
        catalytic_reference_sites: dict[str, int] | None,
    ) -> None:
        self.manifest.parent_prediction_id = parent_prediction_id or ""
        self.manifest.catalytic_reference_sites = dict(catalytic_reference_sites or {})
        self._write_manifest()

    def add_artifact(self, name: str, path: str | Path) -> None:
        path = Path(path).resolve()
        self.manifest.artifacts[name] = str(path)
        self._write_manifest()

    def set_software_version(self, name: str, version: str) -> None:
        self.manifest.software_versions[name] = version
        self._write_manifest()

    def set_database_version(self, name: str, version: str) -> None:
        self.manifest.database_versions[name] = version
        self._write_manifest()

    def add_warning(self, message: str) -> None:
        if message not in self.manifest.warnings:
            self.manifest.warnings.append(message)
        self._write_manifest()
        self.emit("run", "warning", message)

    def complete(self, raw: int, filtered: int, annotated: int) -> RunSummary:
        self.manifest.status = "completed"
        self.manifest.completed_at = datetime.now(timezone.utc)
        self.manifest.raw_hit_count = raw
        self.manifest.filtered_hit_count = filtered
        self.manifest.annotated_hit_count = annotated
        self._write_manifest()
        self.emit("run", "completed", "Pipeline run completed", progress=1.0)
        return self.summary()

    def fail(self, exc: BaseException) -> None:
        self.manifest.status = "failed"
        self.manifest.completed_at = datetime.now(timezone.utc)
        self.manifest.error = f"{type(exc).__name__}: {exc}"
        self._write_manifest()
        self.emit("run", "failed", self.manifest.error)

    def summary(self) -> RunSummary:
        return RunSummary(
            run_id=self.run_id,
            status=self.manifest.status,
            run_dir=self.run_dir,
            manifest_path=self.manifest_path,
            artifacts={key: Path(value) for key, value in self.manifest.artifacts.items()},
            raw_hit_count=self.manifest.raw_hit_count,
            filtered_hit_count=self.manifest.filtered_hit_count,
            annotated_hit_count=self.manifest.annotated_hit_count,
            warnings=list(self.manifest.warnings),
        )
