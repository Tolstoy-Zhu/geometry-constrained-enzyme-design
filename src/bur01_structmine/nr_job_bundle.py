"""Export credential-free full-NR jobs and import verified server results."""
from __future__ import annotations

import json
import os
import re
import shutil
import stat
import zipfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory
from typing import Any

import pandas as pd
import yaml

from .config import resolve_path
from .io_utils import ensure_dir, write_excel_safely, write_fasta
from .models import NrExpansionRequest, NrRunManifest, NrRunSummary
from .nr_pipeline import (
    NrRunContext,
    _accession,
    _extract_structure_seeds,
    _link_parent_nr_run,
    _read_fasta,
    _screen_candidates,
    load_mmseqs_tsv,
    load_nr_config,
    write_nr_checksums,
)
from .prediction_sequence import read_single_fasta
from .run_context import file_sha256, make_run_id


MIN_SERVER_RAM_GB = 128
MIN_SERVER_FREE_DISK_GB = 500
MAX_IMPORT_MEMBERS = 20_000
MAX_IMPORT_UNCOMPRESSED_BYTES = 20 * 1024**3
SAFE_RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")


def _write_text(path: Path, text: str, *, executable: bool = False) -> None:
    path.write_text(text.replace("\r\n", "\n"), encoding="utf-8", newline="\n")
    if executable and os.name != "nt":
        path.chmod(path.stat().st_mode | 0o111)


def _bundle_checksums(bundle_dir: Path) -> Path:
    path = bundle_dir / "BUNDLE_SHA256SUMS"
    lines = [
        f"{file_sha256(item)}  {item.relative_to(bundle_dir).as_posix()}"
        for item in sorted(bundle_dir.rglob("*"))
        if item.is_file() and item != path
    ]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return path


def _archive_directory(directory: Path, destination: Path) -> Path:
    with zipfile.ZipFile(
        destination,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=6,
    ) as archive:
        for path in sorted(item for item in directory.rglob("*") if item.is_file()):
            archive.write(path, arcname=path.relative_to(directory).as_posix())
    return destination


def _resource_script() -> str:
    return f"""#!/usr/bin/env bash
set -euo pipefail

MIN_RAM_GB="${{MIN_RAM_GB:-{MIN_SERVER_RAM_GB}}}"
MIN_FREE_DISK_GB="${{MIN_FREE_DISK_GB:-{MIN_SERVER_FREE_DISK_GB}}}"
TARGET_DIR="${{1:-$PWD}}"
RAM_GB="$(awk '/MemTotal/ {{printf "%d", $2/1024/1024}}' /proc/meminfo)"
FREE_GB="$(df -Pk "$TARGET_DIR" | awk 'NR==2 {{printf "%d", $4/1024/1024}}')"

printf 'RAM_GB=%s\\nFREE_DISK_GB=%s\\n' "$RAM_GB" "$FREE_GB"
if [[ "${{AQUA_SKIP_RESOURCE_CHECK:-0}}" != "1" ]]; then
  (( RAM_GB >= MIN_RAM_GB )) || {{
    echo "ERROR: full NR requires at least $MIN_RAM_GB GiB RAM for this preset." >&2
    exit 2
  }}
  (( FREE_GB >= MIN_FREE_DISK_GB )) || {{
    echo "ERROR: full NR requires at least $MIN_FREE_DISK_GB GiB free disk." >&2
    exit 3
  }}
fi
"""


def _prepare_database_script() -> str:
    return """#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/job.env"
MMSEQS="${MMSEQS:-mmseqs}"
NR_DB="${NR_DB:?Set NR_DB in job.env or the environment}"
TMP_DIR="${TMP_DIR:-$SCRIPT_DIR/work/mmseqs_tmp}"
THREADS="${THREADS:-64}"

mkdir -p "$(dirname "$NR_DB")" "$TMP_DIR"
"$SCRIPT_DIR/resource_precheck.sh" "$(dirname "$NR_DB")"
if [[ -f "${NR_DB}.dbtype" ]]; then
  echo "NR database already exists: $NR_DB"
  exit 0
fi
echo "Explicitly downloading and building the complete NCBI NR database."
"$MMSEQS" databases NR "$NR_DB" "$TMP_DIR" --threads "$THREADS"
test -f "${NR_DB}.dbtype"
"""


def _run_script() -> str:
    return """#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/job.env"
MMSEQS="${MMSEQS:-mmseqs}"
NR_DB="${NR_DB:?Set NR_DB in job.env or the environment}"
THREADS="${THREADS:-64}"
TMP_DIR="${TMP_DIR:-$SCRIPT_DIR/work/mmseqs_tmp}"
OUT_DIR="${OUT_DIR:-$SCRIPT_DIR/output}"
NR_SNAPSHOT_DATE="${NR_SNAPSHOT_DATE:?Set NR_SNAPSHOT_DATE to the source NR snapshot date}"
NR_RUN_ID="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["nr_run_id"])' "$SCRIPT_DIR/job_metadata.json")"

mkdir -p "$OUT_DIR"
"$SCRIPT_DIR/resource_precheck.sh" "$OUT_DIR"
test -f "${NR_DB}.dbtype" || {
  echo "ERROR: full NR database is missing. Run prepare_nr_database.sh explicitly." >&2
  exit 4
}
mkdir -p "$TMP_DIR/search" "$TMP_DIR/cluster" "$TMP_DIR/direct"
cp "$SCRIPT_DIR/job_metadata.json" "$OUT_DIR/job_metadata.json"

"$MMSEQS" easy-search \
  "$SCRIPT_DIR/seeds.fasta" "$NR_DB" "$OUT_DIR/nr_hits.raw.tsv" "$TMP_DIR/search" \
  --num-iterations 2 --max-seqs 2500 --search-type 1 -e 1e-5 \
  --threads "$THREADS" \
  --format-output query,target,pident,fident,nident,alnlen,qlen,tlen,qstart,qend,tstart,tend,qcov,tcov,bits,evalue,qaln,taln,tseq,theader

python3 "$SCRIPT_DIR/nr_server_postprocess.py" deduplicate \
  "$OUT_DIR/nr_hits.raw.tsv" \
  "$OUT_DIR/nr_hits.deduplicated.tsv" \
  "$OUT_DIR/nr_hits.fasta"

"$MMSEQS" easy-cluster \
  "$OUT_DIR/nr_hits.fasta" "$OUT_DIR/nr_s50" "$TMP_DIR/cluster" \
  --min-seq-id 0.5 --threads "$THREADS"

REP_COUNT="$(grep -c '^>' "$OUT_DIR/nr_s50_rep_seq.fasta" || true)"
(( REP_COUNT > 0 )) || {
  echo "ERROR: clustering produced no representative sequences." >&2
  exit 5
}
DIRECT_MAX_SEQS="$REP_COUNT"
(( DIRECT_MAX_SEQS >= 2500 )) || DIRECT_MAX_SEQS=2500

"$MMSEQS" easy-search \
  "$SCRIPT_DIR/bur01_query.fasta" "$OUT_DIR/nr_s50_rep_seq.fasta" \
  "$OUT_DIR/bur01_vs_representatives.tsv" "$TMP_DIR/direct" \
  --num-iterations 1 --max-seqs "$DIRECT_MAX_SEQS" --search-type 1 -e 10 \
  --threads "$THREADS" \
  --format-output query,target,pident,fident,nident,alnlen,qlen,tlen,qstart,qend,tstart,tend,qcov,tcov,bits,evalue,qaln,taln,tseq,theader

python3 "$SCRIPT_DIR/nr_server_postprocess.py" database-manifest \
  "$OUT_DIR" "$NR_DB" "$NR_SNAPSHOT_DATE" --mmseqs "$MMSEQS"
python3 "$SCRIPT_DIR/nr_server_postprocess.py" finalize \
  "$OUT_DIR" "$NR_RUN_ID" "NCBI NR full database" "$NR_SNAPSHOT_DATE" \
  --mmseqs "$MMSEQS"
python3 "$SCRIPT_DIR/nr_server_postprocess.py" archive \
  "$OUT_DIR" "$SCRIPT_DIR/${NR_RUN_ID}_server_result.zip"
echo "Result archive: $SCRIPT_DIR/${NR_RUN_ID}_server_result.zip"
"""


def _slurm_script() -> str:
    return """#!/usr/bin/env bash
#SBATCH --job-name=aqua_nr_v04
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=aqua_nr_%j.log
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bash "$SCRIPT_DIR/run_nr_pipeline.sh"
"""


def _pbs_script() -> str:
    return """#!/usr/bin/env bash
#PBS -N aqua_nr_v04
#PBS -l select=1:ncpus=64:mem=128gb
#PBS -l walltime=48:00:00
#PBS -j oe
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PBS_O_WORKDIR"
bash "$SCRIPT_DIR/run_nr_pipeline.sh"
"""


def _bundle_readme(nr_run_id: str) -> str:
    return f"""# AquaEsterNet-StructMine v0.4 full-NR job

Run ID: `{nr_run_id}`

## 中文

1. 编辑 `job.env`，设置服务器上的 `NR_DB`、`TMP_DIR`、线程数和输出目录。
2. 首次建库时显式运行 `bash prepare_nr_database.sh`。普通检索不会自动下载 NR。
3. 直接运行 `bash run_nr_pipeline.sh`，或使用 `sbatch submit.slurm` / `qsub submit.pbs`。
4. 将生成的 `{nr_run_id}_server_result.zip` 下载回本机，用 `nr-import-job` 导入。
5. 作业包不含 SSH 密码、令牌或服务器地址。`NR_DB` 必须是完整 NCBI NR，而不是本机演示子集。

固定复现参数：MMseqs2 两轮迭代、`max-seqs=2500`、`search-type=1`、50% 序列一致性聚类。

## English

1. Edit `job.env` and set the server-side `NR_DB`, `TMP_DIR`, threads, and output path.
2. For the initial database build, explicitly run `bash prepare_nr_database.sh`. A normal search never downloads NR.
3. Run `bash run_nr_pipeline.sh`, `sbatch submit.slurm`, or `qsub submit.pbs`.
4. Return `{nr_run_id}_server_result.zip` and import it with `nr-import-job`.
5. The bundle contains no SSH credentials, tokens, or server address. `NR_DB` must be the complete NCBI NR database.

The reproducibility preset is fixed to two MMseqs2 iterations, `max-seqs=2500`,
`search-type=1`, and 50% sequence-identity clustering.
"""


def export_nr_job_bundle(request: NrExpansionRequest) -> NrRunSummary:
    """Create a complete, credential-free job bundle without executing full NR."""
    config_path, config = load_nr_config(request.config_path)
    scope = request.database_scope or str(config["database"]["scope"])
    if scope != "full_nr":
        raise ValueError("nr-export-job requires a configuration with database.scope=full_nr")
    request = request.model_copy(
        update={
            "execution_mode": "job_bundle",
            "database_scope": "full_nr",
            "nr_run_id": request.nr_run_id or make_run_id("nr-full"),
        }
    )
    context = NrRunContext(request, config_path, config, "full_nr")
    try:
        context.start()
        bundle_dir = ensure_dir(context.run_dir / "job_bundle")
        seeds_fasta = bundle_dir / "seeds.fasta"
        seeds_tsv = bundle_dir / "seeds.tsv"
        seeds, warnings = _extract_structure_seeds(request, seeds_fasta, seeds_tsv)
        context.manifest.seed_count = len(seeds)
        for warning in warnings:
            context.warning(warning)

        query = read_single_fasta(resolve_path("data/query/Bur01m_seed.fasta"))
        write_fasta([("Bur01_8HFW_A", query.sequence)], bundle_dir / "bur01_query.fasta")
        shutil.copy2(
            resolve_path("scripts/nr_server_postprocess.py"),
            bundle_dir / "nr_server_postprocess.py",
        )
        metadata = {
            "project_version": "0.4.0",
            "nr_run_id": context.nr_run_id,
            "parent_structure_run_id": request.source_run_id,
            "database_scope": "full_nr",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "paper_preset": "hong2025",
            "parameters": config["search"],
            "seed_count": len(seeds),
            "seed_fasta_sha256": file_sha256(seeds_fasta),
        }
        (bundle_dir / "job_metadata.json").write_text(
            json.dumps(metadata, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        (bundle_dir / "config.snapshot.yaml").write_text(
            yaml.safe_dump(config, sort_keys=False, allow_unicode=True),
            encoding="utf-8",
        )
        _write_text(
            bundle_dir / "job.env",
            """# Server-local paths only. Do not place credentials in this file.
MMSEQS="${MMSEQS:-mmseqs}"
NR_DB="${NR_DB:-/data/databases/mmseqs/nr}"
TMP_DIR="${TMP_DIR:-$SCRIPT_DIR/work/mmseqs_tmp}"
OUT_DIR="${OUT_DIR:-$SCRIPT_DIR/output}"
THREADS="${THREADS:-64}"
# Required for auditable imports, for example 2026-07-24.
NR_SNAPSHOT_DATE="${NR_SNAPSHOT_DATE:-}"
""",
        )
        _write_text(bundle_dir / "resource_precheck.sh", _resource_script(), executable=True)
        _write_text(
            bundle_dir / "prepare_nr_database.sh",
            _prepare_database_script(),
            executable=True,
        )
        _write_text(bundle_dir / "run_nr_pipeline.sh", _run_script(), executable=True)
        _write_text(bundle_dir / "submit.slurm", _slurm_script(), executable=True)
        _write_text(bundle_dir / "submit.pbs", _pbs_script(), executable=True)
        _write_text(bundle_dir / "README_CN_EN.md", _bundle_readme(context.nr_run_id))
        checksums = _bundle_checksums(bundle_dir)
        archive_path = context.run_dir / f"{context.nr_run_id}_full_nr_job.zip"
        _archive_directory(bundle_dir, archive_path)

        context.manifest.database_snapshot = {
            "scope": "full_nr",
            "database_name": "NCBI NR",
            "database_path": "server-managed; not stored in this bundle",
            "validated": False,
        }
        context.manifest.artifacts.update(
            {
                "job_bundle_dir": str(bundle_dir.resolve()),
                "job_bundle_zip": str(archive_path.resolve()),
                "bundle_checksums": str(checksums.resolve()),
                "seeds_fasta": str(seeds_fasta.resolve()),
                "seeds_tsv": str(seeds_tsv.resolve()),
            }
        )
        context.warning(
            "Full NR has not been executed; import and validate a server result before release."
        )
        checksum_path = context.run_dir / "SHA256SUMS"
        context.artifact("run_checksums", checksum_path)
        context._write()
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


def _safe_relative_path(value: str) -> PurePosixPath:
    normalized = value.replace("\\", "/")
    relative = PurePosixPath(normalized)
    if (
        not normalized
        or "\x00" in normalized
        or relative.is_absolute()
        or any(part in {"", ".", ".."} or ":" in part for part in relative.parts)
    ):
        raise ValueError(f"Unsafe path in result archive: {value}")
    return relative


def _contained_target(root: Path, relative: PurePosixPath) -> Path:
    resolved_root = root.resolve()
    target = resolved_root.joinpath(*relative.parts).resolve()
    if not target.is_relative_to(resolved_root):
        raise ValueError(f"Unsafe path in result archive: {relative}")
    return target


def _validated_run_id(value: object, label: str = "nr_run_id") -> str:
    run_id = str(value or "").strip()
    if not SAFE_RUN_ID.fullmatch(run_id):
        raise ValueError(f"Invalid {label}: {run_id!r}")
    return run_id


def _safe_extract(archive_path: Path, destination: Path) -> list[Path]:
    extracted: list[Path] = []
    with zipfile.ZipFile(archive_path) as archive:
        members = archive.infolist()
        if len(members) > MAX_IMPORT_MEMBERS:
            raise ValueError("Server result archive contains too many members")
        total_size = sum(member.file_size for member in members)
        if total_size > MAX_IMPORT_UNCOMPRESSED_BYTES:
            raise ValueError("Server result archive is too large after extraction")
        seen: set[PurePosixPath] = set()
        for member in members:
            pure = _safe_relative_path(member.filename.rstrip("/"))
            if pure in seen:
                raise ValueError(f"Duplicate path in result archive: {member.filename}")
            seen.add(pure)
            mode = member.external_attr >> 16
            if stat.S_ISLNK(mode):
                raise ValueError(f"Symbolic links are not allowed: {member.filename}")
            target = _contained_target(destination, pure)
            if member.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            with archive.open(member) as source, target.open("wb") as output:
                shutil.copyfileobj(source, output)
            extracted.append(target)
    return extracted


def _find_single(root: Path, name: str) -> Path:
    matches = [path for path in root.rglob(name) if path.is_file()]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {name}, found {len(matches)}")
    return matches[0]


def _verify_result_checksums(output_dir: Path) -> None:
    checksum_path = output_dir / "RESULT_SHA256SUMS"
    if not checksum_path.is_file():
        raise ValueError("Server result is missing RESULT_SHA256SUMS")
    for line_number, line in enumerate(
        checksum_path.read_text(encoding="ascii").splitlines(),
        start=1,
    ):
        if not line.strip():
            continue
        parts = line.split(maxsplit=1)
        if len(parts) != 2 or len(parts[0]) != 64:
            raise ValueError(f"Malformed checksum line {line_number}")
        try:
            relative = _safe_relative_path(parts[1].strip())
            target = _contained_target(output_dir, relative)
        except ValueError as exc:
            raise ValueError(f"Unsafe checksum path: {parts[1].strip()}") from exc
        if not target.is_file():
            raise ValueError(f"Checksummed result is missing: {relative}")
        if file_sha256(target) != parts[0].lower():
            raise ValueError(f"Checksum mismatch: {relative}")


def _validate_database_manifest(
    path: Path,
    server_manifest: dict[str, Any],
) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("database_scope") != "full_nr":
        raise ValueError("Database manifest is not marked as full_nr")
    snapshot_date = str(payload.get("database_snapshot_date", ""))
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", snapshot_date):
        raise ValueError("Database manifest has no valid snapshot date")
    if snapshot_date != str(server_manifest.get("database_snapshot_date", "")):
        raise ValueError("Database and server manifests use different snapshot dates")
    if payload.get("source_command") != "mmseqs databases NR":
        raise ValueError("Database manifest does not identify the official NR build command")
    if int(payload.get("file_count", 0)) < 1 or int(payload.get("total_bytes", 0)) < 1:
        raise ValueError("Database manifest has invalid file-count or byte-size evidence")
    if not re.fullmatch(r"[0-9a-fA-F]{64}", str(payload.get("dbtype_sha256", ""))):
        raise ValueError("Database manifest has no valid dbtype SHA-256")
    server_version = str(server_manifest.get("mmseqs_version", "")).strip()
    database_version = str(payload.get("mmseqs_version", "")).strip()
    if not database_version or (server_version and database_version != server_version):
        raise ValueError("Database and server manifests use different MMseqs2 versions")
    return payload


def _write_imported_candidates(
    run_dir: Path,
    manifest: NrRunManifest,
    config: dict[str, Any],
    imported_dir: Path,
) -> None:
    deduplicated = pd.read_csv(
        imported_dir / "nr_hits.deduplicated.tsv",
        sep="\t",
    )
    if "target_accession" not in deduplicated:
        deduplicated["target_accession"] = deduplicated["target"].map(_accession)
    representatives_path = imported_dir / "nr_s50_rep_seq.fasta"
    representatives = _read_fasta(representatives_path)
    cluster_path = imported_dir / "nr_s50_cluster.tsv"
    cluster_map = pd.read_csv(
        cluster_path,
        sep="\t",
        names=["representative", "member"],
    )
    direct_path = imported_dir / "bur01_vs_representatives.tsv"
    direct = load_mmseqs_tsv(direct_path)
    candidates = _screen_candidates(
        representatives,
        deduplicated,
        cluster_map,
        direct,
        config["filter"],
    )
    if candidates.empty:
        raise ValueError("Imported full-NR result produced no sequence representatives")

    candidates_tsv = run_dir / "nr_candidates.sequence_screen.tsv"
    candidates_xlsx = run_dir / "nr_candidates.sequence_screen.xlsx"
    candidates_fasta = run_dir / "nr_candidates.fasta"
    af2_queue = run_dir / "af2_ready.tsv"
    candidates.to_csv(candidates_tsv, sep="\t", index=False)
    write_excel_safely(candidates, candidates_xlsx)
    write_fasta(
        [
            (str(row["target_accession"] or row["sequence_id"]), str(row["target_sequence"]))
            for _, row in candidates.iterrows()
        ],
        candidates_fasta,
    )
    candidates[candidates["af2_ready"].fillna(False).astype(bool)].to_csv(
        af2_queue,
        sep="\t",
        index=False,
    )
    manifest.candidate_count = len(candidates)
    manifest.af2_ready_count = int(candidates["af2_ready"].fillna(False).sum())
    manifest.artifacts.update(
        {
            "imported_server_dir": str(imported_dir.resolve()),
            "candidates_tsv": str(candidates_tsv.resolve()),
            "candidates_xlsx": str(candidates_xlsx.resolve()),
            "candidates_fasta": str(candidates_fasta.resolve()),
            "af2_queue_tsv": str(af2_queue.resolve()),
        }
    )


def import_nr_job_bundle(
    archive_path: str | Path,
    *,
    output_root: str | Path = "results/nr_runs",
    structure_run_root: str | Path = "results/runs",
) -> NrRunSummary:
    """Verify and import a full-NR server result archive."""
    archive_path = resolve_path(archive_path)
    if not archive_path.is_file():
        raise FileNotFoundError(archive_path)
    with TemporaryDirectory(prefix="aqua_nr_import_") as temp:
        extracted_root = Path(temp)
        _safe_extract(archive_path, extracted_root)
        server_manifest_path = _find_single(extracted_root, "server_result_manifest.json")
        imported_source = server_manifest_path.parent
        _verify_result_checksums(imported_source)
        server_manifest = json.loads(server_manifest_path.read_text(encoding="utf-8"))
        if server_manifest.get("project_version") != "0.4.0":
            raise ValueError("Server result was not produced by AquaEsterNet v0.4")
        if server_manifest.get("database_scope") != "full_nr":
            raise ValueError("Server result is not marked as a complete NR run")
        nr_run_id = _validated_run_id(server_manifest.get("nr_run_id"))
        metadata_path = imported_source / "job_metadata.json"
        if not metadata_path.is_file():
            raise ValueError("Server result is missing job_metadata.json")
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        if metadata.get("nr_run_id") != nr_run_id:
            raise ValueError("Server result and job metadata use different nr_run_id values")
        if metadata.get("database_scope") != "full_nr":
            raise ValueError("Imported job metadata is not marked full_nr")
        parent_run_id = _validated_run_id(
            metadata.get("parent_structure_run_id"),
            "parent_structure_run_id",
        )
        database_manifest_path = imported_source / "database_manifest.json"
        database_manifest = _validate_database_manifest(
            database_manifest_path,
            server_manifest,
        )

        run_dir = ensure_dir(resolve_path(output_root) / nr_run_id)
        local_manifest_path = run_dir / "nr_run_manifest.json"
        if not local_manifest_path.is_file():
            manifest = NrRunManifest(
                nr_run_id=nr_run_id,
                parent_structure_run_id=parent_run_id,
                execution_mode="import",
                database_scope="full_nr",
                database_name="NCBI NR",
                paper_preset=str(metadata.get("paper_preset", "hong2025")),
                resolved_parameters=dict(metadata.get("parameters") or {}),
            )
        else:
            manifest = NrRunManifest.model_validate_json(
                local_manifest_path.read_text(encoding="utf-8")
            )
            if manifest.database_scope != "full_nr":
                raise ValueError("The local run ID exists but is not a full-NR job")
            if (
                parent_run_id != manifest.parent_structure_run_id
            ):
                raise ValueError(
                    "Imported job metadata does not match the local parent structure run"
                )

        config_path = Path(manifest.config_path) if manifest.config_path else resolve_path(
            "configs/nr_full_server.yaml"
        )
        _, config = load_nr_config(config_path)
        archive_hash = file_sha256(archive_path)
        imported_dir = run_dir / f"server_import_{archive_hash[:12]}"
        if imported_dir.exists():
            raise FileExistsError(
                f"This server result has already been imported: {imported_dir}"
            )
        shutil.copytree(imported_source, imported_dir)
        _write_imported_candidates(run_dir, manifest, config, imported_dir)

        manifest.execution_mode = "import"
        manifest.status = "completed"
        manifest.started_at = manifest.started_at or datetime.now(timezone.utc)
        manifest.completed_at = datetime.now(timezone.utc)
        manifest.database_snapshot = {
            "scope": "full_nr",
            "database_name": server_manifest.get("database_label", "NCBI NR"),
            "snapshot_date": server_manifest.get("database_snapshot_date", ""),
            "server_completed_at": server_manifest.get("completed_at", ""),
            "server_host": server_manifest.get("host", ""),
            "validated": bool(server_manifest.get("database_snapshot_date")),
            "result_archive_sha256": archive_hash,
            "database_manifest_sha256": file_sha256(
                imported_dir / "database_manifest.json"
            ),
            "database_file_count": int(database_manifest["file_count"]),
            "database_total_bytes": int(database_manifest["total_bytes"]),
            "database_dbtype_sha256": str(database_manifest["dbtype_sha256"]),
        }
        manifest.software_versions["mmseqs"] = str(
            server_manifest.get("mmseqs_version", "")
        )
        manifest.raw_hit_count = int(server_manifest.get("raw_hit_count", 0))
        manifest.deduplicated_hit_count = int(
            server_manifest.get("deduplicated_hit_count", 0)
        )
        manifest.representative_count = int(
            server_manifest.get("representative_count", 0)
        )
        manifest.error = ""
        checksums = run_dir / "SHA256SUMS"
        manifest.artifacts["run_checksums"] = str(checksums.resolve())
        manifest.artifacts["database_manifest"] = str(
            (imported_dir / "database_manifest.json").resolve()
        )
        local_manifest_path.write_text(
            manifest.model_dump_json(indent=2),
            encoding="utf-8",
        )
        write_nr_checksums(run_dir)
        link_request = NrExpansionRequest(
            source_run_id=manifest.parent_structure_run_id,
            config_path=str(config_path),
            nr_run_id=nr_run_id,
            output_root=str(output_root),
            structure_run_root=str(structure_run_root),
            database_scope="full_nr",
            execution_mode="import",
        )
        _link_parent_nr_run(link_request, nr_run_id)
        return NrRunSummary(
            nr_run_id=nr_run_id,
            status=manifest.status,
            run_dir=run_dir,
            manifest_path=local_manifest_path,
            artifacts={key: Path(value) for key, value in manifest.artifacts.items()},
            candidate_count=manifest.candidate_count,
            af2_ready_count=manifest.af2_ready_count,
            warnings=manifest.warnings,
        )
