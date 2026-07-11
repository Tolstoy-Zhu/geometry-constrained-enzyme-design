"""Shared AlphaFold2 prediction and prediction-to-mining service layer."""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from .colabfold_runtime import (
    ColabFoldRuntime,
    PredictionCancelled,
    ProcessFileLock,
    build_msa_arguments,
    build_prediction_arguments,
    parse_prediction_models,
    preset_parameters,
)
from .config import resolve_path
from .models import (
    PredictionManifest,
    PredictionRequest,
    PredictionSummary,
    ProgressEvent,
    RunManifest,
    RunSummary,
)
from .pipeline import run_pipeline
from .prediction_sequence import FastaRecord, map_8hfw_catalytic_sites, parse_single_fasta, read_single_fasta, validate_sequence
from .run_context import make_run_id


class PrivacyConsentRequired(PermissionError):
    pass


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    tmp.replace(path)


class PredictionContext:
    def __init__(
        self,
        prediction_root: Path,
        prediction_id: str,
        record: FastaRecord,
        request: PredictionRequest,
        source: dict[str, str],
    ) -> None:
        self.prediction_id = prediction_id
        self.prediction_dir = prediction_root / prediction_id
        self.prediction_dir.mkdir(parents=True, exist_ok=False)
        self.manifest_path = self.prediction_dir / "prediction_manifest.json"
        self.events_path = self.prediction_dir / "events.jsonl"
        self.cancel_path = self.prediction_dir / "cancel.requested"
        self.log_path = self.prediction_dir / "colabfold.log"
        parameters = preset_parameters(request.preset)
        self.manifest = PredictionManifest(
            prediction_id=prediction_id,
            sequence_name=record.name,
            sequence_length=len(record.sequence),
            sequence_sha256=record.sha256,
            source_type=source.get("source_type", "fasta"),
            source_run_id=source.get("source_run_id", ""),
            source_accession=source.get("source_accession", ""),
            preset=request.preset,
            parameters=parameters.to_dict(),
            use_templates=request.use_templates,
            public_msa_confirmed=request.confirm_public_msa,
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
    ) -> None:
        event = ProgressEvent(
            stage=stage,
            state=state,
            message=message,
            progress=progress,
            artifact=str(artifact) if artifact else None,
        )
        with self.events_path.open("a", encoding="utf-8") as handle:
            handle.write(event.model_dump_json() + "\n")

    def add_artifact(self, name: str, path: str | Path) -> None:
        self.manifest.artifacts[name] = str(Path(path).resolve())
        self._write_manifest()

    def start(self) -> None:
        self.manifest.status = "running"
        self.manifest.started_at = datetime.now(timezone.utc)
        self._write_manifest()
        self.emit("prediction", "started", "AlphaFold2 prediction started", 0.0)

    def fail(self, exc: BaseException) -> None:
        self.manifest.status = "cancelled" if isinstance(exc, PredictionCancelled) else "failed"
        self.manifest.completed_at = datetime.now(timezone.utc)
        self.manifest.error = f"{type(exc).__name__}: {exc}"
        self._write_manifest()
        self.emit("prediction", self.manifest.status, self.manifest.error)

    def complete(self) -> PredictionSummary:
        self.manifest.status = "completed"
        self.manifest.completed_at = datetime.now(timezone.utc)
        self._write_manifest()
        self.emit("prediction", "completed", "AlphaFold2 prediction completed", 1.0, self.manifest.rank_001_path)
        return self.summary()

    def summary(self) -> PredictionSummary:
        rank_path = Path(self.manifest.rank_001_path) if self.manifest.rank_001_path else None
        return PredictionSummary(
            prediction_id=self.prediction_id,
            status=self.manifest.status,
            prediction_dir=self.prediction_dir,
            manifest_path=self.manifest_path,
            rank_001_path=rank_path,
            model_count=len(self.manifest.models),
            mean_plddt=self.manifest.mean_plddt,
            ptm=self.manifest.ptm,
            child_run_ids=list(self.manifest.child_run_ids),
            warnings=list(self.manifest.warnings),
        )


def _candidate_record(request: PredictionRequest) -> tuple[FastaRecord, dict[str, str]]:
    run_root = resolve_path(request.candidate_run_root)
    manifest_path = run_root / str(request.candidate_run_id) / "run_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"Candidate run manifest not found: {manifest_path}")
    manifest = RunManifest.model_validate_json(manifest_path.read_text(encoding="utf-8"))
    annotated_path = Path(manifest.artifacts.get("annotated_tsv", ""))
    if not annotated_path.is_file():
        raise FileNotFoundError(f"Candidate table not found: {annotated_path}")
    candidates = pd.read_csv(annotated_path, sep="\t")
    accession = str(request.accession).strip().upper()
    identifier_columns = [
        column
        for column in (
            "target",
            "target_uniprot_accession",
            "uniprot_accessions",
            "ncbi_protein_accession",
            "pdb_id",
        )
        if column in candidates.columns
    ]
    matches = pd.Series(False, index=candidates.index)
    for column in identifier_columns:
        values = candidates[column].fillna("").astype(str).str.upper()
        matches |= values.str.contains(accession, regex=False)
    selected = candidates[matches]
    if selected.empty:
        raise ValueError(f"Accession {accession} was not found in run {request.candidate_run_id}")
    sequence = str(selected.iloc[0].get("final_protein_sequence", ""))
    record = validate_sequence(accession, sequence)
    return record, {
        "source_type": "candidate",
        "source_run_id": str(request.candidate_run_id),
        "source_accession": accession,
    }


def resolve_prediction_input(request: PredictionRequest) -> tuple[FastaRecord, dict[str, str]]:
    if request.fasta_text:
        return parse_single_fasta(request.fasta_text), {"source_type": "fasta_text"}
    if request.fasta_path:
        path = resolve_path(request.fasta_path)
        return read_single_fasta(path), {"source_type": "fasta_file"}
    return _candidate_record(request)


def _copy_tree_contents(source: Path, destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for item in source.iterdir():
        target = destination / item.name
        if item.is_dir():
            shutil.copytree(item, target, dirs_exist_ok=True)
        else:
            shutil.copy2(item, target)


def _cache_key(record: FastaRecord, use_templates: bool) -> str:
    material = f"colabfold-1.6.1|{record.sha256}|templates={int(use_templates)}"
    return hashlib.sha256(material.encode("ascii")).hexdigest()


def _cache_complete(cache_dir: Path) -> bool:
    return (cache_dir / "cache_manifest.json").is_file() and any(cache_dir.rglob("*.a3m"))


def _msa_failure_detail(log_path: Path, start_offset: int) -> str:
    if not log_path.is_file():
        return "ColabFold MSA stage completed without an A3M file"
    with log_path.open("r", encoding="utf-8", errors="replace") as handle:
        handle.seek(start_offset)
        attempt_log = handle.read()
    details: list[str] = []
    for line in attempt_log.splitlines():
        if "Server didn't reply with json:" in line or "Could not get MSA/templates" in line:
            message = line.split(" ", 2)[-1].strip()
            if message and message not in details:
                details.append(message)
    if details:
        return "Public MSA service error: " + " | ".join(details[-2:])
    return "ColabFold MSA stage completed without an A3M file"


def _prepare_msa(
    context: PredictionContext,
    request: PredictionRequest,
    runtime: ColabFoldRuntime,
    input_fasta: Path,
    output_dir: Path,
    cache_root: Path,
) -> Path:
    cache_key = _cache_key(read_single_fasta(input_fasta), request.use_templates)
    cache_dir = cache_root / cache_key
    context.manifest.msa_cache_key = cache_key
    context.manifest.msa_cache_hit = _cache_complete(cache_dir)
    context._write_manifest()
    if context.manifest.msa_cache_hit:
        _copy_tree_contents(cache_dir, output_dir)
        msa_path = sorted(output_dir.rglob("*.a3m"))[0]
        context.emit("msa", "completed", "Reused sequence-hash MSA cache", 0.30, msa_path)
        return msa_path

    if not request.confirm_public_msa:
        raise PrivacyConsentRequired(
            "Public MSA consent is required because this sequence will be sent to the ColabFold MMseqs2 service"
        )

    msa_lock = ProcessFileLock(cache_root / ".public_msa.lock", "public MSA request")
    with msa_lock:
        if _cache_complete(cache_dir):
            _copy_tree_contents(cache_dir, output_dir)
            context.manifest.msa_cache_hit = True
            context._write_manifest()
            return sorted(output_dir.rglob("*.a3m"))[0]
        attempt_root = context.prediction_dir / "msa_download"
        attempt_root.mkdir(parents=True, exist_ok=True)
        last_error: BaseException | None = None
        retry_base_seconds = max(0.0, float(os.environ.get("AQUA_MSA_RETRY_BASE_SECONDS", "15")))
        for attempt in range(1, 4):
            context.emit("msa", "running", f"Public MSA request {attempt}/3", 0.08 + attempt * 0.05)
            try:
                log_offset = context.log_path.stat().st_size if context.log_path.is_file() else 0
                arguments = build_msa_arguments(
                    runtime.path(input_fasta),
                    runtime.path(attempt_root),
                    request.use_templates,
                )
                runtime.run(
                    arguments,
                    context.log_path,
                    context.cancel_path,
                    context.prediction_dir / "colabfold_msa.pid",
                    timeout=60 * 60,
                )
                if not any(attempt_root.rglob("*.a3m")):
                    raise RuntimeError(_msa_failure_detail(context.log_path, log_offset))
                last_error = None
                break
            except PredictionCancelled:
                raise
            except BaseException as exc:  # noqa: BLE001 - retry all external MSA failures.
                last_error = exc
                context.emit("msa", "warning", f"MSA attempt {attempt} failed: {type(exc).__name__}: {exc}")
                if attempt < 3:
                    time.sleep(retry_base_seconds * (2 ** (attempt - 1)))
        if last_error is not None:
            raise RuntimeError(f"Public MSA failed after 3 attempts: {last_error}") from last_error

        cache_root.mkdir(parents=True, exist_ok=True)
        temporary_cache = cache_root / f".{cache_key}.{context.prediction_id}.tmp"
        shutil.rmtree(temporary_cache, ignore_errors=True)
        shutil.copytree(attempt_root, temporary_cache)
        _write_json_atomic(
            temporary_cache / "cache_manifest.json",
            {
                "cache_key": cache_key,
                "sequence_sha256": context.manifest.sequence_sha256,
                "use_templates": request.use_templates,
                "created_at": datetime.now(timezone.utc).isoformat(),
                "service": context.manifest.public_msa_service,
            },
        )
        if cache_dir.exists():
            shutil.rmtree(temporary_cache)
        else:
            temporary_cache.replace(cache_dir)
        _copy_tree_contents(cache_dir, output_dir)

    msa_path = sorted(output_dir.rglob("*.a3m"))[0]
    context.emit("msa", "completed", "Public MSA retrieved and cached", 0.30, msa_path)
    return msa_path


def _write_quality_artifacts(context: PredictionContext, scores_path: Path) -> tuple[Path, Path]:
    scores = json.loads(scores_path.read_text(encoding="utf-8"))
    plddt_path = context.prediction_dir / "rank_001_plddt.json"
    pae_path = context.prediction_dir / "rank_001_pae.json"
    _write_json_atomic(plddt_path, {"plddt": scores.get("plddt", [])})
    _write_json_atomic(
        pae_path,
        {"pae": scores.get("pae", []), "max_pae": scores.get("max_pae")},
    )
    return plddt_path, pae_path


def run_prediction(request: PredictionRequest) -> PredictionSummary:
    record, source = resolve_prediction_input(request)
    prediction_root = resolve_path(request.output_root)
    prediction_root.mkdir(parents=True, exist_ok=True)
    prediction_id = request.prediction_id or make_run_id(f"AF2-{record.name}")
    context = PredictionContext(prediction_root, prediction_id, record, request, source)
    input_fasta = context.prediction_dir / "input.fasta"
    input_fasta.write_text(record.to_fasta(), encoding="ascii")
    context.manifest.input_fasta_path = str(input_fasta.resolve())
    context.add_artifact("input_fasta", input_fasta)
    context.add_artifact("events", context.events_path)
    context.add_artifact("log", context.log_path)

    runtime = ColabFoldRuntime(
        backend=request.runtime_backend,
        executable=request.colabfold_executable,
        wsl_distribution=request.wsl_distribution,
    )
    parameters = preset_parameters(request.preset)
    output_dir = context.prediction_dir / "colabfold"
    cache_root = resolve_path(request.msa_cache_dir)

    try:
        context.start()
        context.manifest.runtime_versions = runtime.versions()
        context.manifest.runtime_patches = runtime.patch_provenance()
        context._write_manifest()
        msa_path = _prepare_msa(context, request, runtime, input_fasta, output_dir, cache_root)
        context.manifest.msa_path = str(msa_path.resolve())
        context.add_artifact("msa", msa_path)

        gpu_lock = ProcessFileLock(prediction_root / ".gpu.lock", "GPU prediction")
        with gpu_lock:
            gpu = runtime.gpu_details()
            context.manifest.runtime_versions["jax_devices"] = "; ".join(gpu.get("devices", []))
            context._write_manifest()
            context.emit("inference", "running", "Running local GPU AlphaFold2-PTM inference", 0.38)
            arguments = build_prediction_arguments(
                runtime.path(input_fasta),
                runtime.path(output_dir),
                parameters,
                request.use_templates,
            )
            runtime.run(
                arguments,
                context.log_path,
                context.cancel_path,
                context.prediction_dir / "colabfold_inference.pid",
                timeout=24 * 60 * 60,
            )
        context.emit("inference", "completed", "Local GPU inference completed", 0.86)

        context.manifest.model_parameter_provenance = runtime.parameter_provenance()
        context.manifest.model_parameter_checksums = runtime.parameter_checksums()
        if not context.manifest.model_parameter_checksums:
            context.manifest.warnings.append("AlphaFold2 model parameter checksums were not available")
        context._write_manifest()

        models = parse_prediction_models(output_dir)
        expected_models = parameters.num_models * parameters.num_seeds
        if len(models) < expected_models:
            raise RuntimeError(
                f"Expected {expected_models} ranked models but parsed {len(models)} from {output_dir}"
            )
        rank_one = models[0]
        canonical_pdb = context.prediction_dir / "rank_001.pdb"
        shutil.copy2(rank_one.pdb_path, canonical_pdb)
        canonical_scores = context.prediction_dir / "rank_001_scores.json"
        shutil.copy2(rank_one.scores_path, canonical_scores)
        plddt_path, pae_path = _write_quality_artifacts(context, canonical_scores)

        mapping = map_8hfw_catalytic_sites(record.sequence, override=request.catalytic_override)
        context.manifest.models = models
        context.manifest.rank_001_path = str(canonical_pdb.resolve())
        context.manifest.mean_plddt = rank_one.mean_plddt
        context.manifest.ptm = rank_one.ptm
        context.manifest.catalytic_mapping = mapping
        if not mapping["complete"]:
            context.manifest.warnings.append(
                "One or more 8HFW catalytic residues could not be mapped with matching identity; catalytic evidence remains incomplete"
            )
        context.add_artifact("rank_001_pdb", canonical_pdb)
        context.add_artifact("rank_001_scores", canonical_scores)
        context.add_artifact("plddt", plddt_path)
        context.add_artifact("pae", pae_path)
        context.emit("quality", "completed", "Ranking, confidence, PAE, and catalytic mapping parsed", 0.96)
        return context.complete()
    except BaseException as exc:
        context.fail(exc)
        raise


def request_prediction_cancel(prediction_dir: str | Path) -> Path:
    prediction_dir = Path(prediction_dir).resolve()
    cancel_path = prediction_dir / "cancel.requested"
    cancel_path.write_text(datetime.now(timezone.utc).isoformat(), encoding="utf-8")
    return cancel_path


def _append_prediction_child_run(prediction_manifest_path: Path, run_id: str) -> None:
    manifest = PredictionManifest.model_validate_json(prediction_manifest_path.read_text(encoding="utf-8"))
    if run_id not in manifest.child_run_ids:
        manifest.child_run_ids.append(run_id)
    tmp = prediction_manifest_path.with_suffix(".json.tmp")
    tmp.write_text(manifest.model_dump_json(indent=2), encoding="utf-8")
    tmp.replace(prediction_manifest_path)


def run_prediction_and_mine(
    request: PredictionRequest,
    mining_config: str = "configs/alphafold_swissprot.yaml",
    mining_output_root: str | None = None,
) -> tuple[PredictionSummary, RunSummary]:
    prediction_summary = run_prediction(request)
    if prediction_summary.rank_001_path is None:
        raise RuntimeError("Prediction completed without rank_001.pdb")
    manifest = PredictionManifest.model_validate_json(
        prediction_summary.manifest_path.read_text(encoding="utf-8")
    )
    catalytic_sites = {
        label: int(position)
        for label, position in manifest.catalytic_mapping.get("valid_positions", {}).items()
    }
    mining = run_pipeline(
        config_path=mining_config,
        output_root=mining_output_root,
        query_structure=str(prediction_summary.rank_001_path),
        query_chain="A",
        parent_prediction_id=prediction_summary.prediction_id,
        catalytic_reference_sites=catalytic_sites,
    )
    _append_prediction_child_run(prediction_summary.manifest_path, mining.run_id)
    updated = PredictionManifest.model_validate_json(
        prediction_summary.manifest_path.read_text(encoding="utf-8")
    )
    return prediction_summary.model_copy(update={"child_run_ids": list(updated.child_run_ids)}), mining
