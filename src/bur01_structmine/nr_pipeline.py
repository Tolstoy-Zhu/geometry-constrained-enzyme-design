"""NCBI NR sequence expansion, catalytic mapping, and novelty audit services."""
from __future__ import annotations

import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import yaml

from .config import resolve_path
from .io_utils import ensure_dir, write_excel_safely, write_fasta
from .mmseqs_runner import (
    MMSEQS_FIELDS,
    get_mmseqs_version,
    run_easy_cluster,
    run_easy_search,
)
from .models import (
    NrExpansionRequest,
    NrNoveltyManifest,
    NrNoveltyRequest,
    NrRunManifest,
    NrRunSummary,
    ProgressEvent,
    RunManifest,
)
from .ncbi_blast import NcbiRemoteBlastClient, hits_to_rows
from .prediction_sequence import (
    FastaRecord,
    SequenceValidationError,
    parse_single_fasta,
    read_single_fasta,
    safe_sequence_name,
    validate_sequence,
)
from .results_parser import parse_ncbi_protein_accession
from .run_context import file_sha256, make_run_id
from .scoring import has_serine_hydrolase_motif
from .structure_search import load_ca_trace


DEFAULT_NR_CONFIG: dict[str, Any] = {
    "database": {
        "scope": "nr_remote_snapshot",
        "name": "NCBI nr remote snapshot",
        "path": "",
        "manifest": "databases/nr_demo/latest.json",
    },
    "runtime": {
        "backend": "wsl",
        "executable": "mmseqs",
        "wsl_distribution": "Ubuntu-Foldseek",
        "wsl_executable": "/opt/colabfold/env/bin/mmseqs",
    },
    "search": {
        "paper_preset": "hong2025",
        "iterations": 2,
        "max_seqs": 2500,
        "search_type": 1,
        "evalue": 1e-5,
        "threads": 16,
        "cluster_min_sequence_identity": 0.5,
    },
    "filter": {
        "max_evalue": 1e-5,
        "min_qcov": 0.60,
        "min_tcov": 0.40,
        "min_length": 180,
        "max_length": 500,
        "max_x_fraction": 0.05,
        "required_catalytic_identities": 3,
    },
    "output": {"run_root": "results/nr_runs"},
}


NUMERIC_MMSEQS_COLUMNS = (
    "pident",
    "fident",
    "nident",
    "alnlen",
    "qlen",
    "tlen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "qcov",
    "tcov",
    "bits",
    "evalue",
)


def _deep_merge(base: dict[str, Any], update: dict[str, Any]) -> dict[str, Any]:
    result = json.loads(json.dumps(base))
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def load_nr_config(path: str | Path) -> tuple[Path, dict[str, Any]]:
    config_path = resolve_path(path)
    if not config_path.is_file():
        raise FileNotFoundError(config_path)
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise ValueError("NR configuration must be a YAML mapping")
    config = _deep_merge(DEFAULT_NR_CONFIG, raw)
    scope = str(config["database"].get("scope", "")).strip()
    if scope not in {"nr_remote_snapshot", "full_nr", "custom_mmseqs"}:
        raise ValueError(f"Unsupported NR database scope: {scope}")
    search = config["search"]
    if int(search["iterations"]) < 1 or int(search["max_seqs"]) < 1:
        raise ValueError("NR search iterations and max_seqs must be positive")
    identity = float(search["cluster_min_sequence_identity"])
    if not 0.0 < identity <= 1.0:
        raise ValueError("cluster_min_sequence_identity must be in (0, 1]")
    if str(search.get("paper_preset", "")).strip() == "hong2025":
        expected = {
            "iterations": 2,
            "max_seqs": 2500,
            "search_type": 1,
        }
        mismatches = [
            f"{key}={search.get(key)!r} (expected {value})"
            for key, value in expected.items()
            if int(search.get(key, -1)) != value
        ]
        if abs(identity - 0.5) > 1e-9:
            mismatches.append(
                "cluster_min_sequence_identity="
                f"{identity!r} (expected 0.5)"
            )
        if mismatches:
            raise ValueError(
                "hong2025 preset parameters are fixed: " + "; ".join(mismatches)
            )
    return config_path, config


def _write_json_atomic(path: Path, payload: dict[str, Any]) -> None:
    ensure_dir(path.parent)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    tmp.replace(path)


def _events_append(path: Path, event: ProgressEvent) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(event.model_dump_json() + "\n")


class NrRunContext:
    def __init__(
        self,
        request: NrExpansionRequest,
        config_path: Path,
        config: dict[str, Any],
        scope: str,
    ) -> None:
        self.nr_run_id = request.nr_run_id or make_run_id("nr")
        self.run_dir = resolve_path(request.output_root) / self.nr_run_id
        self.run_dir.mkdir(parents=True, exist_ok=False)
        self.manifest_path = self.run_dir / "nr_run_manifest.json"
        self.events_path = self.run_dir / "events.jsonl"
        self.manifest = NrRunManifest(
            nr_run_id=self.nr_run_id,
            parent_structure_run_id=request.source_run_id,
            execution_mode=request.execution_mode,
            database_scope=scope,
            database_name=str(config["database"].get("name") or "NCBI nr"),
            config_path=str(config_path.resolve()),
            config_sha256=file_sha256(config_path),
            config_snapshot=config,
            paper_preset=str(config["search"].get("paper_preset", "hong2025")),
            resolved_parameters=dict(config["search"]),
        )
        self._write()

    def _write(self) -> None:
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
        _events_append(
            self.events_path,
            ProgressEvent(
                stage=stage,
                state=state,
                message=message,
                progress=progress,
                artifact=str(artifact) if artifact else None,
            ),
        )

    def start(self) -> None:
        self.manifest.status = "running"
        self.manifest.started_at = datetime.now(timezone.utc)
        self._write()
        self.emit("nr_run", "started", "NR sequence expansion started", 0.0)

    def artifact(self, name: str, path: str | Path) -> None:
        self.manifest.artifacts[name] = str(Path(path).resolve())
        self._write()

    def warning(self, message: str) -> None:
        if message not in self.manifest.warnings:
            self.manifest.warnings.append(message)
            self._write()
        self.emit("nr_run", "warning", message)

    def complete(self) -> NrRunSummary:
        self.manifest.status = "completed"
        self.manifest.completed_at = datetime.now(timezone.utc)
        self._write()
        self.emit("nr_run", "completed", "NR sequence expansion completed", 1.0)
        return self.summary()

    def fail(self, exc: BaseException) -> None:
        self.manifest.status = "failed"
        self.manifest.completed_at = datetime.now(timezone.utc)
        self.manifest.error = f"{type(exc).__name__}: {exc}"
        self._write()
        self.emit("nr_run", "failed", self.manifest.error)

    def summary(self) -> NrRunSummary:
        return NrRunSummary(
            nr_run_id=self.nr_run_id,
            status=self.manifest.status,
            run_dir=self.run_dir,
            manifest_path=self.manifest_path,
            artifacts={key: Path(value) for key, value in self.manifest.artifacts.items()},
            candidate_count=self.manifest.candidate_count,
            af2_ready_count=self.manifest.af2_ready_count,
            warnings=list(self.manifest.warnings),
        )


def _structure_manifest(request: NrExpansionRequest) -> tuple[Path, RunManifest]:
    path = (
        resolve_path(request.structure_run_root)
        / request.source_run_id
        / "run_manifest.json"
    )
    if not path.is_file():
        raise FileNotFoundError(f"Structure run manifest not found: {path}")
    return path, RunManifest.model_validate_json(path.read_text(encoding="utf-8"))


def _clean_sequence(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return re.sub(r"[^A-Za-zXx]", "", str(value or "")).upper()


def _first_sequence(*values: object) -> str:
    for value in values:
        sequence = _clean_sequence(value)
        if sequence:
            return sequence
    return ""


def _extract_structure_seeds(
    request: NrExpansionRequest,
    out_fasta: Path,
    out_tsv: Path,
) -> tuple[list[FastaRecord], list[str]]:
    _, manifest = _structure_manifest(request)
    warnings: list[str] = []
    sequence_by_target: dict[str, str] = {}
    annotated_path = Path(manifest.artifacts.get("annotated_tsv", ""))
    if annotated_path.is_file():
        annotated = pd.read_csv(annotated_path, sep="\t")
        for _, row in annotated.iterrows():
            target = str(row.get("target", "")).strip()
            sequence = _first_sequence(
                row.get("final_protein_sequence", ""),
                row.get("tseq", ""),
                row.get("ncbi_protein_sequence", ""),
            )
            if target and sequence:
                sequence_by_target[target] = sequence

    filtered_path = Path(manifest.artifacts.get("filtered_tsv", ""))
    filtered = pd.read_csv(filtered_path, sep="\t") if filtered_path.is_file() else pd.DataFrame()
    rows: list[dict[str, object]] = []
    records: list[FastaRecord] = []
    seen_hashes: set[str] = set()

    query_fasta_path = Path(manifest.artifacts.get("query_fasta", ""))
    if query_fasta_path.is_file():
        query = read_single_fasta(query_fasta_path)
        records.append(FastaRecord("Bur01_8HFW_A", query.sequence))
        seen_hashes.add(query.sha256)
        rows.append(
            {
                "seed_id": "Bur01_8HFW_A",
                "source_db": "query",
                "target": manifest.query_pdb_id or "8HFW",
                "sequence_sha256": query.sha256,
                "length": len(query.sequence),
            }
        )

    for index, row in filtered.iterrows():
        target = str(row.get("target", "")).strip()
        sequence = _first_sequence(
            row.get("tseq", ""),
            sequence_by_target.get(target, ""),
        )
        if not sequence:
            warnings.append(f"No sequence was available for structural hit {target or index}")
            continue
        try:
            record = validate_sequence(target or f"seed_{index + 1}", sequence)
        except SequenceValidationError as exc:
            warnings.append(f"Rejected structural seed {target or index}: {exc}")
            continue
        if record.sha256 in seen_hashes:
            continue
        seen_hashes.add(record.sha256)
        seed_id = f"seed_{len(records):05d}_{safe_sequence_name(target)}"
        records.append(FastaRecord(seed_id, record.sequence))
        rows.append(
            {
                "seed_id": seed_id,
                "source_db": row.get("source_db", ""),
                "target": target,
                "sequence_sha256": record.sha256,
                "length": len(record.sequence),
            }
        )

    if not records:
        raise ValueError(
            "No valid sequence seeds could be extracted; rerun Foldseek with target sequence output"
        )
    write_fasta([(record.name, record.sequence) for record in records], out_fasta)
    pd.DataFrame(rows).to_csv(out_tsv, sep="\t", index=False)
    return records, warnings


def load_mmseqs_tsv(path: str | Path, fields: Iterable[str] = MMSEQS_FIELDS) -> pd.DataFrame:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    fields = list(fields)
    if path.stat().st_size == 0:
        return pd.DataFrame(columns=fields)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        first = handle.readline().rstrip("\r\n")
    observed = len(first.split("\t"))
    columns = fields[:observed]
    if observed > len(columns):
        columns.extend(f"extra_{index}" for index in range(observed - len(columns)))
    frame = pd.read_csv(path, sep="\t", names=columns, comment="#")
    if not frame.empty and str(frame.iloc[0].get("query", "")).lower() == "query":
        frame = pd.read_csv(path, sep="\t", comment="#")
    for column in NUMERIC_MMSEQS_COLUMNS:
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    return frame


def _accession(value: object) -> str:
    text = str(value or "").strip()
    return parse_ncbi_protein_accession(text) or text.split(maxsplit=1)[0]


def deduplicate_hits(raw: pd.DataFrame) -> pd.DataFrame:
    if raw.empty:
        return raw.copy()
    frame = raw.copy()
    frame["target_accession"] = frame["target"].map(_accession)
    sequences = [
        _first_sequence(row.get("tseq", ""), row.get("target_sequence", ""))
        for _, row in frame.iterrows()
    ]
    frame["sequence_sha256"] = [
        hashlib.sha256(sequence.encode("ascii")).hexdigest() if sequence else ""
        for sequence in sequences
    ]
    frame = frame.sort_values(
        [column for column in ("evalue", "bits", "qcov", "tcov") if column in frame],
        ascending=[
            True if column == "evalue" else False
            for column in ("evalue", "bits", "qcov", "tcov")
            if column in frame
        ],
        na_position="last",
    )
    frame = frame.drop_duplicates("target_accession", keep="first").reset_index(drop=True)
    frame["_dedup_sequence_key"] = frame.apply(
        lambda row: (
            str(row["sequence_sha256"])
            if row["sequence_sha256"]
            else f"accession:{row['target_accession']}"
        ),
        axis=1,
    )
    frame = frame.drop_duplicates("_dedup_sequence_key", keep="first").drop(
        columns=["_dedup_sequence_key"]
    )
    frame = frame.reset_index(drop=True)
    frame.insert(0, "sequence_id", [f"nr_{index + 1:07d}" for index in range(len(frame))])
    return frame


def _read_fasta(path: str | Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = ""
    sequence: list[str] = []
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                records[name] = "".join(sequence).upper()
            name = line[1:].split(maxsplit=1)[0]
            sequence = []
        else:
            sequence.append(line)
    if name:
        records[name] = "".join(sequence).upper()
    return records


def _write_hit_fasta(frame: pd.DataFrame, path: Path) -> None:
    records: list[tuple[str, str]] = []
    for _, row in frame.iterrows():
        sequence = _first_sequence(
            row.get("tseq", ""),
            row.get("target_sequence", ""),
        )
        if sequence:
            records.append((str(row["sequence_id"]), sequence))
    if not records:
        raise ValueError("No complete target sequences were returned by MMseqs2")
    write_fasta(records, path)


def _reference_site_positions() -> dict[str, tuple[int, str]]:
    trace = load_ca_trace(resolve_path("data/query/Bur01m.pdb"), chain_id="A")
    index_by_number = {
        residue.residue_number: index + 1 for index, residue in enumerate(trace.residues)
    }
    sites: dict[str, tuple[int, str]] = {}
    for label, number, amino_acid in (
        ("S111", 111, "S"),
        ("D241", 241, "D"),
        ("H274", 274, "H"),
    ):
        query_position = index_by_number.get(number)
        if query_position is None:
            raise ValueError(f"8HFW reference residue {label} is missing from chain A")
        sites[label] = (query_position, amino_acid)
    return sites


def map_catalytic_alignment(row: pd.Series | None) -> dict[str, Any]:
    sites = _reference_site_positions()
    if row is None:
        return {
            "identity_count": 0,
            "complete": False,
            "status": "unmapped",
            "mapping": {
                label: {
                    "target_position": None,
                    "expected_amino_acid": expected,
                    "observed_amino_acid": None,
                    "status": "unmapped",
                }
                for label, (_, expected) in sites.items()
            },
        }
    qaln = str(row.get("qaln", "") or "")
    taln = str(row.get("taln", "") or "")
    qpos = int(row.get("qstart", 1) or 1) - 1
    tpos = int(row.get("tstart", 1) or 1) - 1
    mapped: dict[str, dict[str, Any]] = {}
    target_by_query: dict[int, tuple[int, str]] = {}
    for qchar, tchar in zip(qaln, taln):
        if qchar != "-":
            qpos += 1
        if tchar != "-":
            tpos += 1
        if qchar != "-" and tchar != "-":
            target_by_query[qpos] = (tpos, tchar.upper())
    identities = 0
    for label, (query_position, expected) in sites.items():
        target = target_by_query.get(query_position)
        if target is None:
            status = "unmapped"
            target_position = None
            observed = None
        else:
            target_position, observed = target
            status = "mapped" if observed == expected else "identity_mismatch"
            identities += int(status == "mapped")
        mapped[label] = {
            "query_sequence_position": query_position,
            "target_position": target_position,
            "expected_amino_acid": expected,
            "observed_amino_acid": observed,
            "status": status,
        }
    return {
        "identity_count": identities,
        "complete": identities == 3,
        "status": "complete" if identities == 3 else ("partial" if identities else "unmapped"),
        "mapping": mapped,
    }


def _screen_candidates(
    representatives: dict[str, str],
    deduplicated: pd.DataFrame,
    cluster_map: pd.DataFrame,
    direct: pd.DataFrame,
    filters: dict[str, Any],
) -> pd.DataFrame:
    by_sequence_id = deduplicated.set_index("sequence_id", drop=False)
    direct_by_target = {
        str(row["target"]): row for _, row in direct.iterrows()
    }
    cluster_sizes = (
        cluster_map.groupby("representative").size().to_dict()
        if not cluster_map.empty and {"representative", "member"} <= set(cluster_map.columns)
        else {}
    )
    rows: list[dict[str, Any]] = []
    for representative_id, sequence in representatives.items():
        if representative_id not in by_sequence_id.index:
            continue
        source = by_sequence_id.loc[representative_id]
        if isinstance(source, pd.DataFrame):
            source = source.iloc[0]
        alignment = direct_by_target.get(representative_id)
        mapping = map_catalytic_alignment(alignment)
        length = len(sequence)
        x_fraction = sequence.count("X") / length if length else 1.0
        evalue = float(source.get("evalue")) if pd.notna(source.get("evalue")) else None
        qcov = float(source.get("qcov")) if pd.notna(source.get("qcov")) else None
        tcov = float(source.get("tcov")) if pd.notna(source.get("tcov")) else None
        reasons: list[str] = []
        if evalue is None or evalue > float(filters["max_evalue"]):
            reasons.append("evalue")
        if qcov is None or qcov < float(filters["min_qcov"]):
            reasons.append("query_coverage")
        if tcov is None or tcov < float(filters["min_tcov"]):
            reasons.append("target_coverage")
        if not int(filters["min_length"]) <= length <= int(filters["max_length"]):
            reasons.append("length")
        if x_fraction > float(filters["max_x_fraction"]):
            reasons.append("ambiguous_sequence")
        if mapping["identity_count"] < int(filters["required_catalytic_identities"]):
            reasons.append("catalytic_triad")
        af2_ready = not reasons
        target = str(source.get("target", representative_id))
        title = str(source.get("theader", "") or "")
        organism_match = re.search(r"\[([^\[\]]+)\]\s*$", title)
        rows.append(
            {
                "target": target,
                "target_accession": str(source.get("target_accession", "") or _accession(target)),
                "sequence_id": representative_id,
                "source_seed": str(source.get("query", "")),
                "target_title": title,
                "organism": organism_match.group(1) if organism_match else "",
                "target_sequence": sequence,
                "target_length": length,
                "cluster_size": int(cluster_sizes.get(representative_id, 1)),
                "evalue": evalue,
                "bits": source.get("bits"),
                "pident": source.get("pident"),
                "fident": source.get("fident"),
                "qcov": qcov,
                "tcov": tcov,
                "direct_evalue": alignment.get("evalue") if alignment is not None else None,
                "direct_fident": alignment.get("fident") if alignment is not None else None,
                "direct_qcov": alignment.get("qcov") if alignment is not None else None,
                "direct_tcov": alignment.get("tcov") if alignment is not None else None,
                "motif_gxsxg_like": has_serine_hydrolase_motif(sequence),
                "mapped_catalytic_identity_count": mapping["identity_count"],
                "catalytic_mapping_status": mapping["status"],
                "catalytic_mapping": json.dumps(mapping["mapping"], ensure_ascii=False),
                "sequence_screen_status": "pass" if af2_ready else "hold:" + ",".join(reasons),
                "af2_ready": af2_ready,
                "wase_screen_score": None,
                "wase_evidence_coverage": None,
                "wase_tier": "Not evaluated - structure required",
                "wase_scope": "not_evaluated_sequence_only",
                "wase_missing_evidence": "predicted_or_experimental_structure;tm_score;rmsd;catalytic_geometry",
            }
        )
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame = frame.sort_values(
        [
            "af2_ready",
            "mapped_catalytic_identity_count",
            "motif_gxsxg_like",
            "qcov",
            "tcov",
            "bits",
        ],
        ascending=[False, False, False, False, False, False],
        na_position="last",
    ).reset_index(drop=True)
    frame.insert(0, "nr_rank", range(1, len(frame) + 1))
    return frame


def _database_from_config(
    request: NrExpansionRequest,
    config: dict[str, Any],
) -> tuple[Path, dict[str, Any]]:
    database = config["database"]
    explicit = request.database_path or str(database.get("path", "") or "")
    if explicit:
        path = resolve_path(explicit)
        snapshot = {
            "scope": request.database_scope or database["scope"],
            "database_path": str(path.resolve()),
        }
        return path, snapshot
    pointer_path = resolve_path(database.get("manifest", "databases/nr_demo/latest.json"))
    if not pointer_path.is_file():
        raise FileNotFoundError(
            f"NR demo database pointer not found: {pointer_path}; run nr-demo-build first"
        )
    pointer = json.loads(pointer_path.read_text(encoding="utf-8"))
    database_path = Path(pointer.get("fasta_path", ""))
    if not database_path.is_file():
        raise FileNotFoundError(f"NR demo FASTA from pointer is missing: {database_path}")
    snapshot_path = Path(pointer.get("manifest_path", ""))
    snapshot = (
        json.loads(snapshot_path.read_text(encoding="utf-8"))
        if snapshot_path.is_file()
        else pointer
    )
    return database_path, snapshot


def _link_parent_nr_run(request: NrExpansionRequest, nr_run_id: str) -> None:
    path, manifest = _structure_manifest(request)
    if nr_run_id not in manifest.child_nr_run_ids:
        manifest.child_nr_run_ids.append(nr_run_id)
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(manifest.model_dump_json(indent=2), encoding="utf-8")
        tmp.replace(path)


def run_nr_expansion(request: NrExpansionRequest) -> NrRunSummary:
    config_path, config = load_nr_config(request.config_path)
    scope = request.database_scope or str(config["database"]["scope"])
    context = NrRunContext(request, config_path, config, scope)
    try:
        context.start()
        snapshot_path = context.run_dir / "config.snapshot.yaml"
        snapshot_path.write_text(
            yaml.safe_dump(config, sort_keys=False, allow_unicode=True),
            encoding="utf-8",
        )
        context.artifact("config_snapshot", snapshot_path)
        context.artifact("events", context.events_path)

        seeds_fasta = context.run_dir / "seeds.fasta"
        seeds_tsv = context.run_dir / "seeds.tsv"
        seeds, seed_warnings = _extract_structure_seeds(request, seeds_fasta, seeds_tsv)
        context.manifest.seed_count = len(seeds)
        context._write()
        context.artifact("seeds_fasta", seeds_fasta)
        context.artifact("seeds_tsv", seeds_tsv)
        for warning in seed_warnings:
            context.warning(warning)
        context.emit("seeds", "completed", f"Prepared {len(seeds)} unique sequence seeds", 0.12)

        database_path, database_snapshot = _database_from_config(request, config)
        context.manifest.database_snapshot = database_snapshot
        runtime = config["runtime"]
        search = config["search"]
        context.manifest.software_versions["mmseqs"] = get_mmseqs_version(
            str(runtime.get("executable", "mmseqs")),
            backend=str(runtime.get("backend", "wsl")),
            wsl_distribution=str(runtime.get("wsl_distribution", "Ubuntu-Foldseek")),
            wsl_executable=str(
                runtime.get("wsl_executable", "/opt/colabfold/env/bin/mmseqs")
            ),
        )
        context._write()

        raw_tsv = context.run_dir / "nr_hits.raw.tsv"
        log_path = context.run_dir / "mmseqs.log"
        context.emit("search", "running", "Searching configured NR sequence database", 0.20)
        run_easy_search(
            seeds_fasta,
            database_path,
            raw_tsv,
            context.run_dir / "tmp" / "search",
            threads=int(search["threads"]),
            iterations=int(search["iterations"]),
            max_seqs=int(search["max_seqs"]),
            search_type=int(search["search_type"]),
            evalue=float(search["evalue"]),
            executable=str(runtime.get("executable", "mmseqs")),
            backend=str(runtime.get("backend", "wsl")),
            wsl_distribution=str(runtime.get("wsl_distribution", "Ubuntu-Foldseek")),
            wsl_executable=str(
                runtime.get("wsl_executable", "/opt/colabfold/env/bin/mmseqs")
            ),
            log_path=log_path,
        )
        raw = load_mmseqs_tsv(raw_tsv)
        context.manifest.raw_hit_count = len(raw)
        context.artifact("raw_hits_tsv", raw_tsv)
        context.artifact("mmseqs_log", log_path)
        if raw.empty:
            raise RuntimeError("MMseqs2 NR search returned no hits")

        deduplicated = deduplicate_hits(raw)
        dedup_tsv = context.run_dir / "nr_hits.deduplicated.tsv"
        deduplicated.to_csv(dedup_tsv, sep="\t", index=False)
        context.manifest.deduplicated_hit_count = len(deduplicated)
        context.artifact("deduplicated_hits_tsv", dedup_tsv)

        targets_fasta = context.run_dir / "nr_hits.fasta"
        _write_hit_fasta(deduplicated, targets_fasta)
        context.artifact("hit_sequences_fasta", targets_fasta)
        context.emit(
            "search",
            "completed",
            f"Retained {len(deduplicated)} unique target accessions",
            0.45,
        )

        context.emit("cluster", "running", "Clustering NR targets at 50% identity", 0.50)
        cluster_prefix = context.run_dir / "clusters" / "nr_s50"
        representatives_fasta, cluster_tsv = run_easy_cluster(
            targets_fasta,
            cluster_prefix,
            context.run_dir / "tmp" / "cluster",
            min_sequence_identity=float(search["cluster_min_sequence_identity"]),
            threads=int(search["threads"]),
            executable=str(runtime.get("executable", "mmseqs")),
            backend=str(runtime.get("backend", "wsl")),
            wsl_distribution=str(runtime.get("wsl_distribution", "Ubuntu-Foldseek")),
            wsl_executable=str(
                runtime.get("wsl_executable", "/opt/colabfold/env/bin/mmseqs")
            ),
            log_path=log_path,
        )
        representatives = _read_fasta(representatives_fasta)
        cluster_map = pd.read_csv(
            cluster_tsv,
            sep="\t",
            names=["representative", "member"],
        )
        context.manifest.representative_count = len(representatives)
        context.artifact("representatives_fasta", representatives_fasta)
        context.artifact("clusters_tsv", cluster_tsv)
        context.emit(
            "cluster",
            "completed",
            f"Produced {len(representatives)} cluster representatives",
            0.65,
        )

        direct_tsv = context.run_dir / "bur01_vs_representatives.tsv"
        canonical_query = context.run_dir / "bur01_query.fasta"
        query_record = read_single_fasta(resolve_path("data/query/Bur01m_seed.fasta"))
        write_fasta([("Bur01_8HFW_A", query_record.sequence)], canonical_query)
        direct_fields = MMSEQS_FIELDS
        run_easy_search(
            canonical_query,
            representatives_fasta,
            direct_tsv,
            context.run_dir / "tmp" / "direct",
            threads=int(search["threads"]),
            iterations=1,
            max_seqs=max(len(representatives), int(search["max_seqs"])),
            search_type=1,
            evalue=max(10.0, float(search["evalue"])),
            fields=direct_fields,
            executable=str(runtime.get("executable", "mmseqs")),
            backend=str(runtime.get("backend", "wsl")),
            wsl_distribution=str(runtime.get("wsl_distribution", "Ubuntu-Foldseek")),
            wsl_executable=str(
                runtime.get("wsl_executable", "/opt/colabfold/env/bin/mmseqs")
            ),
            log_path=log_path,
        )
        direct = load_mmseqs_tsv(direct_tsv, direct_fields)
        context.artifact("direct_bur01_alignments_tsv", direct_tsv)

        candidates = _screen_candidates(
            representatives,
            deduplicated,
            cluster_map,
            direct,
            config["filter"],
        )
        candidates_tsv = context.run_dir / "nr_candidates.sequence_screen.tsv"
        candidates_xlsx = context.run_dir / "nr_candidates.sequence_screen.xlsx"
        candidates_fasta = context.run_dir / "nr_candidates.fasta"
        af2_queue = context.run_dir / "af2_ready.tsv"
        candidates.to_csv(candidates_tsv, sep="\t", index=False)
        write_excel_safely(candidates, candidates_xlsx)
        write_fasta(
            [
                (str(row["target_accession"] or row["sequence_id"]), str(row["target_sequence"]))
                for _, row in candidates.iterrows()
            ],
            candidates_fasta,
        )
        if "af2_ready" in candidates:
            candidates[candidates["af2_ready"].fillna(False).astype(bool)].to_csv(
                af2_queue,
                sep="\t",
                index=False,
            )
        else:
            candidates.head(0).to_csv(af2_queue, sep="\t", index=False)
        context.manifest.candidate_count = len(candidates)
        context.manifest.af2_ready_count = int(candidates.get("af2_ready", pd.Series(dtype=bool)).sum())
        context.artifact("candidates_tsv", candidates_tsv)
        context.artifact("candidates_xlsx", candidates_xlsx)
        context.artifact("candidates_fasta", candidates_fasta)
        context.artifact("af2_queue_tsv", af2_queue)
        context._write()
        context.emit(
            "screen",
            "completed",
            f"{context.manifest.af2_ready_count} representatives are AF2-ready; WASE remains unevaluated",
            0.94,
            candidates_xlsx,
        )
        checksum_path = context.run_dir / "SHA256SUMS"
        context.artifact("run_checksums", checksum_path)
        _link_parent_nr_run(request, context.nr_run_id)
        summary = context.complete()
        write_nr_checksums(context.run_dir)
        return summary
    except BaseException as exc:
        context.fail(exc)
        try:
            write_nr_checksums(context.run_dir)
        except OSError:
            pass
        raise


def build_nr_demo_snapshot(
    query_fasta: str | Path,
    *,
    contact_email: str,
    confirm_ncbi: bool,
    max_hits: int = 2500,
    evalue: float = 1e-5,
    output_root: str | Path = "databases/nr_demo",
    cache_dir: str | Path = ".cache/nr_blast",
    api_key: str = "",
    client: NcbiRemoteBlastClient | None = None,
) -> Path:
    record = read_single_fasta(resolve_path(query_fasta))
    client = client or NcbiRemoteBlastClient(
        contact_email,
        resolve_path(cache_dir),
    )
    hits, provenance = client.search(
        record,
        max_hits=max_hits,
        evalue=evalue,
        confirm_remote=confirm_ncbi,
    )
    all_hits = client.fetch_target_sequences(hits, api_key=api_key)
    retrieved_hits = [hit for hit in all_hits if hit.target_sequence]
    if not retrieved_hits:
        raise RuntimeError("NCBI nr search returned no retrievable protein sequences")
    snapshot_id = (
        datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        + "_"
        + record.sha256[:12]
        + f"_n{len(retrieved_hits)}"
    )
    snapshot_dir = ensure_dir(resolve_path(output_root) / snapshot_id)
    fasta_path = snapshot_dir / "nr_remote_snapshot.fasta"
    hits_path = snapshot_dir / "nr_remote_hits.tsv"
    manifest_path = snapshot_dir / "nr_snapshot_manifest.json"
    write_fasta(
        [
            (
                " ".join(
                    value
                    for value in (
                        hit.target_accession,
                        hit.target_title,
                        f"[{hit.organism}]" if hit.organism else "",
                    )
                    if value
                ),
                hit.target_sequence,
            )
            for hit in retrieved_hits
        ],
        fasta_path,
    )
    pd.DataFrame(hits_to_rows(all_hits)).to_csv(hits_path, sep="\t", index=False)
    manifest = {
        "project_version": "0.4.0",
        "database_scope": "nr_remote_snapshot",
        "database_name": "NCBI nr",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "query_name": record.name,
        "query_sha256": record.sha256,
        "raw_hit_count": len(all_hits),
        "sequence_count": len(retrieved_hits),
        "sequence_fetch_error_count": sum(
            bool(hit.annotation_error) for hit in all_hits
        ),
        "max_hits": max_hits,
        "evalue": evalue,
        "fasta_path": str(fasta_path.resolve()),
        "fasta_sha256": file_sha256(fasta_path),
        "hits_path": str(hits_path.resolve()),
        "hits_sha256": file_sha256(hits_path),
        "remote_provenance": provenance,
        "limitations": [
            "This is a query-derived NCBI nr hit snapshot, not the complete nr database.",
            "It is suitable for local workflow validation and demonstrations only.",
        ],
    }
    _write_json_atomic(manifest_path, manifest)
    _write_json_atomic(
        resolve_path(output_root) / "latest.json",
        {
            "snapshot_id": snapshot_id,
            "fasta_path": str(fasta_path.resolve()),
            "manifest_path": str(manifest_path.resolve()),
        },
    )
    return manifest_path


def _novelty_record(request: NrNoveltyRequest) -> FastaRecord:
    if request.fasta_text:
        return parse_single_fasta(request.fasta_text)
    return read_single_fasta(resolve_path(str(request.fasta_path)))


def run_nr_novelty(request: NrNoveltyRequest) -> tuple[NrNoveltyManifest, Path]:
    record = _novelty_record(request)
    novelty_id = request.novelty_id or make_run_id(f"nr-novelty-{record.name}")
    run_dir = resolve_path(request.output_root) / novelty_id
    run_dir.mkdir(parents=True, exist_ok=False)
    manifest_path = run_dir / "nr_novelty_manifest.json"
    manifest = NrNoveltyManifest(
        novelty_id=novelty_id,
        status="running",
        started_at=datetime.now(timezone.utc),
        backend=request.backend,
        database_scope=request.database_scope,
        sequence_name=record.name,
        sequence_sha256=record.sha256,
        parameters={"max_hits": request.max_hits, "evalue": request.evalue},
    )
    _write_json_atomic(manifest_path, json.loads(manifest.model_dump_json()))
    try:
        query_path = run_dir / "query.fasta"
        query_path.write_text(record.to_fasta(), encoding="ascii")
        hits_path = run_dir / "nr_novelty_hits.tsv"
        if request.backend == "ncbi_remote":
            email = request.contact_email or os.environ.get("AQUA_NCBI_EMAIL", "")
            client = NcbiRemoteBlastClient(email, resolve_path(".cache/nr_blast"))
            hits, provenance = client.search(
                record,
                max_hits=request.max_hits,
                evalue=request.evalue,
                confirm_remote=request.confirm_ncbi,
            )
            frame = pd.DataFrame(hits_to_rows(hits))
            manifest.database_snapshot = provenance
            manifest.software_versions["search_engine"] = "NCBI BLAST URL API"
        else:
            raw_path = run_dir / "nr_novelty.raw.tsv"
            run_easy_search(
                query_path,
                resolve_path(str(request.database_path)),
                raw_path,
                run_dir / "tmp",
                iterations=1,
                max_seqs=request.max_hits,
                search_type=1,
                evalue=request.evalue,
                executable=request.mmseqs_executable,
                backend=request.runtime_backend,
                wsl_distribution=request.wsl_distribution,
                wsl_executable=request.wsl_executable,
            )
            frame = load_mmseqs_tsv(raw_path)
            frame["target_accession"] = frame["target"].map(_accession)
            manifest.software_versions["mmseqs"] = get_mmseqs_version(
                request.mmseqs_executable,
                backend=request.runtime_backend,
                wsl_distribution=request.wsl_distribution,
                wsl_executable=request.wsl_executable,
            )
            manifest.database_snapshot = {
                "scope": request.database_scope,
                "database_path": str(resolve_path(str(request.database_path))),
            }
        frame.to_csv(hits_path, sep="\t", index=False)
        manifest.artifacts = {
            "query_fasta": str(query_path.resolve()),
            "hits_tsv": str(hits_path.resolve()),
        }
        if not frame.empty:
            nearest = frame.sort_values(
                ["evalue", "bits"],
                ascending=[True, False],
                na_position="last",
            ).iloc[0]
            identity = nearest.get("fident")
            if pd.isna(identity):
                pident = nearest.get("pident")
                identity = float(pident) / 100.0 if pd.notna(pident) else None
            coverage = nearest.get("qcov")
            manifest.nearest_nr_accession = str(
                nearest.get("target_accession", "") or _accession(nearest.get("target", ""))
            )
            manifest.nearest_nr_identity = float(identity) if identity is not None and pd.notna(identity) else None
            manifest.nearest_nr_coverage = float(coverage) if coverage is not None and pd.notna(coverage) else None
            manifest.nearest_nr_evalue = float(nearest["evalue"]) if pd.notna(nearest.get("evalue")) else None
            manifest.nearest_nr_taxid = str(nearest.get("taxid", "") or "")
            manifest.nearest_nr_title = str(
                nearest.get("target_title", "") or nearest.get("theader", "") or ""
            )
            manifest.nearest_nr_organism = str(nearest.get("organism", "") or "")
            if not manifest.nearest_nr_organism:
                organism_match = re.search(
                    r"\[([^\[\]]+)\]\s*$",
                    manifest.nearest_nr_title,
                )
                if organism_match:
                    manifest.nearest_nr_organism = organism_match.group(1)
            if (
                manifest.nearest_nr_identity is not None
                and manifest.nearest_nr_identity >= 0.999
                and manifest.nearest_nr_coverage is not None
                and manifest.nearest_nr_coverage >= 0.99
            ):
                manifest.self_hit_warning = (
                    "The nearest NR match is effectively identical; this is expected for an "
                    "NR-derived natural sequence and is not novelty evidence."
                )
        manifest.status = "completed"
        manifest.completed_at = datetime.now(timezone.utc)
        checksum_path = run_dir / "SHA256SUMS"
        manifest.artifacts["run_checksums"] = str(checksum_path.resolve())
        _write_json_atomic(manifest_path, json.loads(manifest.model_dump_json()))
        write_nr_checksums(run_dir)
        return manifest, manifest_path
    except BaseException as exc:
        manifest.status = "failed"
        manifest.completed_at = datetime.now(timezone.utc)
        manifest.error = f"{type(exc).__name__}: {exc}"
        _write_json_atomic(manifest_path, json.loads(manifest.model_dump_json()))
        try:
            write_nr_checksums(run_dir)
        except OSError:
            pass
        raise


def write_nr_checksums(run_dir: str | Path) -> Path:
    run_dir = Path(run_dir)
    checksum_path = run_dir / "SHA256SUMS"
    lines: list[str] = []
    for path in sorted(item for item in run_dir.rglob("*") if item.is_file()):
        if path == checksum_path or "tmp" in path.parts:
            continue
        lines.append(f"{file_sha256(path)}  {path.relative_to(run_dir).as_posix()}")
    checksum_path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return checksum_path
