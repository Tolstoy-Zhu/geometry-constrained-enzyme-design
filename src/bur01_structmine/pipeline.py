"""Versioned AquaEsterNet-StructMine search, annotation, and export pipeline."""
from __future__ import annotations

import json
import hashlib
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from .annotator import annotate_hits, export_annotation_outputs
from .config import load_config, resolve_path
from .datasources import SearchDatabase, load_search_databases
from .foldseek_runner import (
    FoldseekUnavailable,
    create_foldseek_database,
    download_foldseek_database,
    get_foldseek_version,
    run_easy_search,
)
from .io_utils import ensure_dir, write_excel_safely, write_fasta
from .models import RunSummary
from .ncbi_commands import write_ncbi_download_script
from .pdb_utils import extract_sequence_from_pdb
from .results_parser import filter_hits, load_foldseek_tsv
from .run_context import ProgressCallback, RunContext
from .scoring import add_candidate_scores
from .structure_search import LOCAL_TM_FIELDS, run_local_tm_search

LOCAL_ENGINES = {"python_ca", "local_ca", "tmtools_tmalign", "tmalign"}


def _database_fingerprint(source: SearchDatabase) -> str:
    path = source.structure_dir if source.search_engine in LOCAL_ENGINES else source.target_db
    if path is None:
        return "unconfigured"
    candidates: list[Path]
    if path.is_dir():
        candidates = sorted(item for item in path.rglob("*") if item.is_file())
    else:
        candidates = sorted(item for item in path.parent.glob(path.name + "*") if item.is_file())
    digest = hashlib.sha256()
    total = 0
    for item in candidates:
        stat = item.stat()
        total += stat.st_size
        digest.update(f"{item.name}\0{stat.st_size}\0{stat.st_mtime_ns}\n".encode("utf-8"))
    label = source.foldseek_database or source.kind
    version_path = Path(f"{path}.version")
    database_version = version_path.read_text(encoding="utf-8").strip() if version_path.is_file() else "unknown"
    database_version = " ".join(database_version.split())
    return (
        f"{label};version={database_version};path={path};files={len(candidates)};"
        f"bytes={total};metadata_sha256={digest.hexdigest()}"
    )


def export_query_sequence(query_pdb: Path, out_fasta: Path, chain_id: str | None = None) -> None:
    """Extract the selected query chain and write a provenance FASTA."""
    sequences = extract_sequence_from_pdb(query_pdb, chain_id=chain_id)
    if chain_id and chain_id not in sequences:
        raise ValueError(f"Query chain {chain_id!r} was not found in {query_pdb}")
    records = [(f"{query_pdb.stem}_chain_{chain}", sequence) for chain, sequence in sequences.items()]
    if not records:
        raise ValueError(f"No protein sequence could be extracted from {query_pdb}")
    write_fasta(records, out_fasta)


def _source_raw_tsv(out_prefix: Path, source: SearchDatabase, multiple_sources: bool) -> Path:
    if not multiple_sources:
        return out_prefix.with_suffix(".search.raw.tsv")
    return Path(f"{out_prefix}.{source.safe_name}.search.raw.tsv")


def _as_bool(value: object, default: bool = False) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _runtime_kwargs(fs_cfg: dict[str, Any]) -> dict[str, str]:
    return {
        "backend": str(fs_cfg.get("backend", "auto")),
        "wsl_distribution": str(fs_cfg.get("wsl_distribution", "Ubuntu-Foldseek")),
        "wsl_executable": str(fs_cfg.get("wsl_executable", "/opt/foldseek/bin/foldseek")),
    }


def _prepare_search_database(source: SearchDatabase, fs_cfg: dict[str, Any]) -> None:
    if source.search_engine in LOCAL_ENGINES:
        if not source.structure_dir:
            raise ValueError(f"Search database {source.name} uses {source.search_engine} but has no structure_dir")
        if not source.structure_dir.exists():
            raise FileNotFoundError(f"Structure directory not found: {source.structure_dir}")
        return

    executable = str(fs_cfg.get("executable", "foldseek"))
    tmp_dir = resolve_path(fs_cfg.get("tmp_dir", "tmp"))
    runtime = _runtime_kwargs(fs_cfg)
    if not source.target_db:
        raise ValueError(f"Search database {source.name} has no target_db")

    if source.foldseek_database:
        download_foldseek_database(
            source.foldseek_database,
            source.target_db,
            tmp_dir,
            executable=executable,
            archive_url=source.archive_url,
            version_url=source.version_url,
            archive_sha256=source.archive_sha256,
            archive_size=source.archive_size,
            archive_prefix=source.archive_prefix,
            **runtime,
        )
        return
    if source.input_fasta:
        if not source.prostt5_model:
            raise ValueError(
                f"Search database {source.name} uses input_fasta but has no prostt5_model"
            )
        create_foldseek_database(
            source.input_fasta,
            source.target_db,
            executable=executable,
            prostt5_model=source.prostt5_model,
            use_gpu=source.use_gpu,
            extra_args=source.prepare_extra_args,
            **runtime,
        )
        return
    if source.structure_dir:
        create_foldseek_database(
            source.structure_dir,
            source.target_db,
            executable=executable,
            use_gpu=source.use_gpu,
            extra_args=source.prepare_extra_args,
            **runtime,
        )


def prepare_search_databases(
    config_path: str = "config.yaml",
    only: set[str] | None = None,
    dry_run: bool = False,
) -> list[SearchDatabase]:
    """Prepare enabled databases or inspect the exact planned operation."""
    cfg = load_config(config_path)
    fs_cfg = cfg.get("foldseek", {})
    sources = load_search_databases(cfg)
    selected = [source for source in sources if only is None or source.name in only]
    missing = sorted(only - {source.name for source in sources}) if only else []
    if missing:
        raise ValueError(f"Unknown enabled search database(s): {', '.join(missing)}")
    for source in selected:
        if dry_run:
            action = "use-tm-align" if source.search_engine in LOCAL_ENGINES else "inspect"
            detail = source.foldseek_database or source.input_fasta or source.structure_dir or source.target_db or ""
            if source.foldseek_database:
                action = "download-foldseek-database"
            elif source.input_fasta:
                action = "create-from-fasta"
            elif source.structure_dir and source.search_engine not in LOCAL_ENGINES:
                action = "create-from-structures"
            print(f"[DRY-RUN] {source.name}: {action} -> {detail}")
        else:
            _prepare_search_database(source, fs_cfg)
    return selected


def _run_or_load_source(
    source: SearchDatabase,
    query_pdb: Path,
    query_chain: str | None,
    out_prefix: Path,
    fs_cfg: dict[str, Any],
    skip_search: bool,
    multiple_sources: bool,
    catalytic_reference_sites: dict[str, int] | None,
) -> pd.DataFrame:
    if skip_search:
        if not source.result_tsv:
            raise ValueError(f"Search database {source.name} has no result_tsv for skipped search")
        raw_tsv = source.result_tsv
        fields = source.format_fields
        print(f"[INFO] Reuse search TSV for {source.name}: {raw_tsv}")
    elif source.search_engine in LOCAL_ENGINES:
        if not source.structure_dir:
            raise ValueError(f"Search database {source.name} has no structure_dir")
        raw_tsv = _source_raw_tsv(out_prefix, source, multiple_sources)
        site_map = (
            {"S111": 111, "D241": 241, "H274": 274}
            if catalytic_reference_sites is None
            else catalytic_reference_sites
        )
        catalytic_residues = tuple(
            (int(position), label[0].upper())
            for label, position in site_map.items()
            if label and label[0].upper() in {"S", "D", "H"}
        )
        run_local_tm_search(
            query_pdb,
            source.structure_dir,
            raw_tsv,
            query_chain=query_chain,
            catalytic_query_residues=catalytic_residues,
        )
        fields = LOCAL_TM_FIELDS
    else:
        if not source.target_db:
            raise ValueError(f"Search database {source.name} has no target_db")
        raw_tsv = _source_raw_tsv(out_prefix, source, multiple_sources)
        run_easy_search(
            query_pdb=query_pdb,
            target_db=source.target_db,
            out_tsv=raw_tsv,
            tmp_dir=resolve_path(fs_cfg.get("tmp_dir", "tmp")),
            executable=str(fs_cfg.get("executable", "foldseek")),
            max_seqs=int(fs_cfg.get("max_seqs", 1000)),
            threads=int(fs_cfg.get("threads", 8)),
            format_fields=source.format_fields,
            cluster_search=source.cluster_search,
            extra_args=source.search_extra_args,
            **_runtime_kwargs(fs_cfg),
        )
        fields = source.format_fields

    raw_df = load_foldseek_tsv(raw_tsv, fields=fields)
    raw_df.insert(0, "source_db", source.name)
    raw_df.insert(1, "source_kind", source.kind)
    raw_df.insert(2, "annotation_source", source.annotation_source)
    raw_df.insert(3, "target_database", str(source.target_db or source.structure_dir or ""))
    if "search_engine" not in raw_df.columns:
        raw_df.insert(4, "search_engine", source.search_engine or "foldseek")
    if "metric_source" not in raw_df.columns:
        raw_df.insert(5, "metric_source", "Foldseek" if source.search_engine not in LOCAL_ENGINES else "TM-align")
    if "query_catalytic_residues" not in raw_df.columns:
        labels = [
            f"{label[0].upper()}{int(position)}"
            for label, position in (catalytic_reference_sites or {}).items()
            if label
        ]
        raw_df["query_catalytic_residues"] = ";".join(labels)
    return raw_df


def _write_latest_pointer(run_root: Path, summary: RunSummary) -> None:
    pointer = run_root / "latest.json"
    tmp = pointer.with_suffix(".json.tmp")
    tmp.write_text(
        json.dumps(
            {
                "run_id": summary.run_id,
                "status": summary.status,
                "run_dir": str(summary.run_dir.resolve()),
                "manifest": str(summary.manifest_path.resolve()),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    tmp.replace(pointer)


def _link_parent_prediction(parent_prediction_id: str | None, run_id: str) -> None:
    if not parent_prediction_id:
        return
    from .models import PredictionManifest

    manifest_path = resolve_path("results/predictions") / parent_prediction_id / "prediction_manifest.json"
    if not manifest_path.is_file():
        return
    manifest = PredictionManifest.model_validate_json(manifest_path.read_text(encoding="utf-8"))
    if run_id in manifest.child_run_ids:
        return
    manifest.child_run_ids.append(run_id)
    tmp = manifest_path.with_suffix(".json.tmp")
    tmp.write_text(manifest.model_dump_json(indent=2), encoding="utf-8")
    tmp.replace(manifest_path)


def run_pipeline(
    config_path: str = "config.yaml",
    skip_foldseek: bool = False,
    input_tsv: str | None = None,
    run_id: str | None = None,
    output_root: str | None = None,
    legacy_output: bool = False,
    progress_callback: ProgressCallback | None = None,
    query_structure: str | None = None,
    query_chain: str | None = None,
    parent_prediction_id: str | None = None,
    catalytic_reference_sites: dict[str, int] | None = None,
) -> RunSummary:
    """Run search, filtering, annotation, WASE-screen, and versioned export."""
    cfg = load_config(config_path)
    config_file = resolve_path(config_path)
    query_cfg = cfg["query"]
    if query_structure:
        query_pdb = resolve_path(query_structure)
        query_cfg["pdb_file"] = str(query_pdb)
        if parent_prediction_id:
            query_cfg["pdb_id"] = None
            query_cfg["residue_range"] = None
            query_cfg["source_url"] = None
    else:
        query_pdb = resolve_path(query_cfg["pdb_file"])
    selected_query_chain = query_chain if query_chain is not None else query_cfg.get("chain")
    query_cfg["chain"] = selected_query_chain
    cfg["lineage"] = {
        "parent_prediction_id": parent_prediction_id,
        "catalytic_reference_sites": catalytic_reference_sites or {},
    }
    configured_prefix = resolve_path(cfg["output"]["out_prefix"])
    run_root = resolve_path(output_root or cfg["output"].get("run_root", "results/runs"))
    label = configured_prefix.name
    context = RunContext(
        run_root=run_root,
        config_path=config_file,
        config_snapshot=cfg,
        run_id=run_id,
        label=label,
        callback=progress_callback,
    )
    out_prefix = configured_prefix if legacy_output else context.run_dir / configured_prefix.name
    ensure_dir(out_prefix.parent)

    try:
        context.start()
        context.set_query(
            query_pdb,
            pdb_id=str(query_cfg.get("pdb_id") or ""),
            chain=str(selected_query_chain or ""),
            residue_range=str(query_cfg.get("residue_range") or ""),
            source_url=str(query_cfg.get("source_url") or ""),
        )
        context.set_prediction_lineage(parent_prediction_id, catalytic_reference_sites)
        snapshot_path = context.run_dir / "config.snapshot.yaml"
        snapshot_path.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8")
        context.add_artifact("config_snapshot", snapshot_path)
        context.add_artifact("events", context.events_path)

        query_fasta = out_prefix.with_suffix(".query_from_pdb.fasta")
        export_query_sequence(query_pdb, query_fasta, chain_id=selected_query_chain)
        context.add_artifact("query_fasta", query_fasta)
        context.emit("query", "completed", "Query provenance and sequence prepared", progress=0.08, artifact=query_fasta)

        fs_cfg = cfg.get("foldseek", {})
        sources = load_search_databases(cfg)
        skip_search = skip_foldseek
        if skip_foldseek and input_tsv:
            sources = [
                SearchDatabase(
                    name="provided_tsv",
                    kind="foldseek",
                    annotation_source="auto",
                    result_tsv=resolve_path(input_tsv),
                )
            ]
        if not skip_search:
            for source in sources:
                if source.auto_prepare or _as_bool(fs_cfg.get("prepare_databases", False)):
                    context.emit("database", "running", f"Preparing {source.name}", progress=0.12)
                    _prepare_search_database(source, fs_cfg)

        for source in sources:
            context.set_database_version(source.name, _database_fingerprint(source))

        if any(source.search_engine not in LOCAL_ENGINES for source in sources) and not skip_search:
            try:
                context.set_software_version(
                    "foldseek",
                    get_foldseek_version(
                        str(fs_cfg.get("executable", "foldseek")),
                        **_runtime_kwargs(fs_cfg),
                    ),
                )
            except FoldseekUnavailable:
                raise

        context.emit("search", "running", "Running configured structure searches", progress=0.18)
        raw_tables = [
            _run_or_load_source(
                source,
                query_pdb,
                selected_query_chain,
                out_prefix,
                fs_cfg,
                skip_search,
                multiple_sources=len(sources) > 1,
                catalytic_reference_sites=catalytic_reference_sites,
            )
            for source in sources
        ]
        raw_df = pd.concat(raw_tables, ignore_index=True) if raw_tables else pd.DataFrame()
        raw_tsv = out_prefix.with_suffix(".search.raw.tsv")
        raw_xlsx = out_prefix.with_suffix(".search.raw.xlsx")
        raw_df.to_csv(raw_tsv, sep="\t", index=False)
        write_excel_safely(raw_df, raw_xlsx)
        context.add_artifact("raw_tsv", raw_tsv)
        context.add_artifact("raw_xlsx", raw_xlsx)
        context.emit("search", "completed", f"Structure search produced {len(raw_df)} rows", progress=0.42, artifact=raw_xlsx)

        filter_cfg = cfg.get("filter", {})
        filtered_df = filter_hits(
            raw_df,
            min_alntmscore=float(filter_cfg.get("min_alntmscore", 0.45)),
            max_rmsd=float(filter_cfg.get("max_rmsd", 4.0)),
            min_qcov=float(filter_cfg.get("min_qcov", 0.60)),
            min_tcov=float(filter_cfg.get("min_tcov", 0.40)),
            min_alnlen=int(filter_cfg.get("min_alnlen", 180)),
            max_sequence_identity=float(filter_cfg.get("max_sequence_identity", 0.95)),
        )
        filtered_tsv = out_prefix.with_suffix(".search.filtered.tsv")
        filtered_xlsx = out_prefix.with_suffix(".search.filtered.xlsx")
        filtered_df.to_csv(filtered_tsv, sep="\t", index=False)
        write_excel_safely(filtered_df, filtered_xlsx)
        context.add_artifact("filtered_tsv", filtered_tsv)
        context.add_artifact("filtered_xlsx", filtered_xlsx)
        context.emit("filter", "completed", f"{len(filtered_df)} rows passed structural thresholds", progress=0.55, artifact=filtered_xlsx)

        annotation_cfg = cfg.get("annotation", {})
        annotation_limit = int(annotation_cfg.get("max_records", 200))
        annotation_input = filtered_df.head(annotation_limit).copy()
        if len(filtered_df) > annotation_limit:
            context.add_warning(
                f"Annotation was limited to the top {annotation_limit} of {len(filtered_df)} filtered hits"
            )
        context.emit("annotation", "running", "Retrieving cached cross-database annotations", progress=0.62)
        annotated_df = annotate_hits(
            annotation_input,
            enable_rcsb=_as_bool(annotation_cfg.get("enable_rcsb", True), True),
            enable_uniprot=_as_bool(annotation_cfg.get("enable_uniprot", True), True),
            enable_ncbi=_as_bool(annotation_cfg.get("enable_ncbi", True), True),
            enable_alphafold=_as_bool(annotation_cfg.get("enable_alphafold", True), True),
            request_sleep_seconds=float(annotation_cfg.get("request_sleep_seconds", 0.15)),
            cache_dir=resolve_path(annotation_cfg.get("cache_dir", ".cache/annotations")),
        )
        context.emit("annotation", "completed", f"Annotated {len(annotated_df)} candidate rows", progress=0.78)

        annotated_df = add_candidate_scores(annotated_df)
        annotated_xlsx = out_prefix.with_suffix(".annotated.xlsx")
        annotated_tsv = out_prefix.with_suffix(".annotated.tsv")
        annotated_df.to_csv(annotated_tsv, sep="\t", index=False)
        export_annotation_outputs(annotated_df, out_prefix)
        ncbi_script = out_prefix.with_suffix(".download_ncbi_cds.sh")
        write_ncbi_download_script(annotated_df, ncbi_script)
        candidate_fasta = out_prefix.with_suffix(".candidate_proteins.fasta")
        context.add_artifact("annotated_xlsx", annotated_xlsx)
        context.add_artifact("annotated_tsv", annotated_tsv)
        context.add_artifact("candidate_fasta", candidate_fasta)
        context.add_artifact("ncbi_script", ncbi_script)
        context.emit("score", "completed", "WASE-screen v0.1 scoring completed", progress=0.92, artifact=annotated_xlsx)

        summary = context.complete(len(raw_df), len(filtered_df), len(annotated_df))
        _write_latest_pointer(run_root, summary)
        _link_parent_prediction(parent_prediction_id, summary.run_id)
        print("\n[DONE] AquaEsterNet-StructMine v0.4 pipeline finished.")
        print(f"Run ID:           {summary.run_id}")
        print(f"Run directory:    {summary.run_dir}")
        print(f"Annotated table:  {annotated_xlsx}")
        return summary
    except BaseException as exc:
        context.fail(exc)
        raise


def reannotate_run(run_dir: str | Path) -> pd.DataFrame:
    """Refresh annotation and WASE-screen outputs for an existing run directory."""
    from .models import ProgressEvent, RunManifest

    run_dir = Path(run_dir).resolve()
    manifest_path = run_dir / "run_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Run manifest not found: {manifest_path}")
    manifest = RunManifest.model_validate_json(manifest_path.read_text(encoding="utf-8"))
    cfg = manifest.config_snapshot
    filtered_path = Path(manifest.artifacts.get("filtered_tsv", ""))
    if not filtered_path.exists():
        raise FileNotFoundError(f"Filtered result table not found: {filtered_path}")
    filtered_df = pd.read_csv(filtered_path, sep="\t")
    annotation_cfg = cfg.get("annotation", {})
    annotated_df = annotate_hits(
        filtered_df.head(int(annotation_cfg.get("max_records", 200))),
        enable_rcsb=_as_bool(annotation_cfg.get("enable_rcsb", True), True),
        enable_uniprot=_as_bool(annotation_cfg.get("enable_uniprot", True), True),
        enable_ncbi=_as_bool(annotation_cfg.get("enable_ncbi", True), True),
        enable_alphafold=_as_bool(annotation_cfg.get("enable_alphafold", True), True),
        request_sleep_seconds=float(annotation_cfg.get("request_sleep_seconds", 0.15)),
        cache_dir=resolve_path(annotation_cfg.get("cache_dir", ".cache/annotations")),
    )
    annotated_df = add_candidate_scores(annotated_df)
    annotated_tsv = Path(manifest.artifacts.get("annotated_tsv", run_dir / "candidates.annotated.tsv"))
    annotated_xlsx = Path(manifest.artifacts.get("annotated_xlsx", annotated_tsv.with_suffix(".xlsx")))
    prefix_text = str(annotated_tsv)
    suffix = ".annotated.tsv"
    out_prefix = Path(prefix_text[: -len(suffix)]) if prefix_text.endswith(suffix) else annotated_tsv.with_suffix("")
    annotated_df.to_csv(annotated_tsv, sep="\t", index=False)
    write_excel_safely(annotated_df, annotated_xlsx)
    export_annotation_outputs(annotated_df, out_prefix)
    manifest.annotated_hit_count = len(annotated_df)
    manifest.artifacts["annotated_tsv"] = str(annotated_tsv.resolve())
    manifest.artifacts["annotated_xlsx"] = str(annotated_xlsx.resolve())
    tmp_manifest = manifest_path.with_suffix(".json.tmp")
    tmp_manifest.write_text(manifest.model_dump_json(indent=2), encoding="utf-8")
    tmp_manifest.replace(manifest_path)
    event = ProgressEvent(stage="annotation", state="completed", message=f"Re-annotated {len(annotated_df)} candidates", progress=1.0)
    with (run_dir / "events.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(event.model_dump_json() + "\n")
    return annotated_df
