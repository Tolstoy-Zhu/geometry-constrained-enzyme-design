#!/usr/bin/env python3
"""Dependency-free post-processing bundled with full-NR server jobs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import re
import subprocess
import zipfile
from datetime import datetime, timezone
from pathlib import Path


FIELDS = [
    "query",
    "target",
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
    "qaln",
    "taln",
    "tseq",
    "theader",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_fasta(records: list[tuple[str, str]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def accession(value: str) -> str:
    for pattern in (
        r"\b(?:ref|gb|emb|dbj)\|([^|]+)\|",
        r"\b([A-Z]{2}_[0-9]+(?:\.[0-9]+)?)\b",
        r"\b([A-Z]{3,4}[0-9]{5,9}(?:\.[0-9]+)?)\b",
    ):
        match = re.search(pattern, value, flags=re.IGNORECASE)
        if match:
            return match.group(1).upper()
    return value.split(maxsplit=1)[0].upper()


def deduplicate(raw_path: Path, out_tsv: Path, out_fasta: Path) -> None:
    best: dict[str, dict[str, str]] = {}
    with raw_path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", fieldnames=FIELDS)
        for row in reader:
            target = row.get("target", "").strip()
            sequence = "".join(char for char in row.get("tseq", "").upper() if char.isalpha())
            if not target or not sequence:
                continue
            row["tseq"] = sequence
            key = accession(target)
            current = best.get(key)
            rank = (float(row.get("evalue") or "inf"), -float(row.get("bits") or 0.0))
            current_rank = (
                float(current.get("evalue") or "inf"),
                -float(current.get("bits") or 0.0),
            ) if current else None
            if current is None or rank < current_rank:
                best[key] = row

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    records: list[tuple[str, str]] = []
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["sequence_id", "sequence_sha256", *FIELDS],
            extrasaction="ignore",
        )
        writer.writeheader()
        seen_sequences: set[str] = set()
        retained_index = 0
        for row in sorted(
            best.values(),
            key=lambda item: (
                float(item.get("evalue") or "inf"),
                -float(item.get("bits") or 0.0),
            ),
        ):
            sequence_hash = hashlib.sha256(row["tseq"].upper().encode("ascii")).hexdigest()
            if sequence_hash in seen_sequences:
                continue
            seen_sequences.add(sequence_hash)
            retained_index += 1
            index = retained_index
            sequence_id = f"nr_{index:09d}"
            writer.writerow(
                {
                    "sequence_id": sequence_id,
                    "sequence_sha256": sequence_hash,
                    **row,
                }
            )
            records.append((sequence_id, row["tseq"]))
    write_fasta(records, out_fasta)


def count_fasta(path: Path) -> int:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return sum(line.startswith(">") for line in handle)


def count_rows(path: Path, has_header: bool = False) -> int:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        count = sum(bool(line.strip()) for line in handle)
    return max(0, count - int(has_header))


def write_database_manifest(
    output_dir: Path,
    database_prefix: Path,
    database_snapshot_date: str,
    executable: str,
) -> None:
    files = sorted(
        path
        for path in database_prefix.parent.glob(database_prefix.name + "*")
        if path.is_file()
    )
    dbtype = Path(f"{database_prefix}.dbtype")
    if not dbtype.is_file():
        raise SystemExit(f"Complete NR database is missing {dbtype.name}")
    payload = {
        "database_scope": "full_nr",
        "database_name": "NCBI NR",
        "database_snapshot_date": database_snapshot_date,
        "database_prefix_name": database_prefix.name,
        "mmseqs_version": mmseqs_version(executable),
        "file_count": len(files),
        "total_bytes": sum(path.stat().st_size for path in files),
        "dbtype_sha256": sha256(dbtype),
        "source_command": "mmseqs databases NR",
    }
    (output_dir / "database_manifest.json").write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )


def mmseqs_version(executable: str) -> str:
    command = [executable, "version"]
    if os.name == "nt" and Path(executable).suffix.lower() in {".cmd", ".bat"}:
        command = ["cmd.exe", "/d", "/c", executable, "version"]
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        errors="replace",
        check=False,
    )
    return (result.stdout or result.stderr).strip()


def finalize(
    output_dir: Path,
    nr_run_id: str,
    database_label: str,
    database_snapshot_date: str,
    executable: str,
) -> None:
    raw = output_dir / "nr_hits.raw.tsv"
    deduplicated = output_dir / "nr_hits.deduplicated.tsv"
    representatives = output_dir / "nr_s50_rep_seq.fasta"
    clusters = output_dir / "nr_s50_cluster.tsv"
    direct = output_dir / "bur01_vs_representatives.tsv"
    database_manifest = output_dir / "database_manifest.json"
    required = [
        raw,
        deduplicated,
        representatives,
        clusters,
        direct,
        database_manifest,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("Missing server result files: " + ", ".join(missing))

    manifest = {
        "nr_run_id": nr_run_id,
        "project_version": "0.4.0",
        "database_scope": "full_nr",
        "database_label": database_label,
        "database_snapshot_date": database_snapshot_date,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "host": platform.node(),
        "platform": platform.platform(),
        "mmseqs_version": mmseqs_version(executable),
        "raw_hit_count": count_rows(raw),
        "deduplicated_hit_count": count_rows(deduplicated, has_header=True),
        "representative_count": count_fasta(representatives),
        "environment": {
            "scheduler_job_id": os.environ.get(
                "SLURM_JOB_ID", os.environ.get("PBS_JOBID", "")
            )
        },
    }
    manifest_path = output_dir / "server_result_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    checksum_path = output_dir / "RESULT_SHA256SUMS"
    lines = [
        f"{sha256(path)}  {path.name}"
        for path in sorted(output_dir.iterdir())
        if path.is_file() and path != checksum_path
    ]
    checksum_path.write_text("\n".join(lines) + "\n", encoding="ascii")


def archive_results(output_dir: Path, archive_path: Path) -> None:
    required = output_dir / "RESULT_SHA256SUMS"
    if not required.is_file():
        raise SystemExit("Run finalize before creating the result archive")
    archive_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(
        archive_path,
        "w",
        compression=zipfile.ZIP_DEFLATED,
        compresslevel=6,
    ) as archive:
        for path in sorted(item for item in output_dir.iterdir() if item.is_file()):
            archive.write(path, arcname=path.name)


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    dedupe = subparsers.add_parser("deduplicate")
    dedupe.add_argument("raw", type=Path)
    dedupe.add_argument("out_tsv", type=Path)
    dedupe.add_argument("out_fasta", type=Path)
    finish = subparsers.add_parser("finalize")
    finish.add_argument("output_dir", type=Path)
    finish.add_argument("nr_run_id")
    finish.add_argument("database_label")
    finish.add_argument("database_snapshot_date")
    finish.add_argument("--mmseqs", default="mmseqs")
    archive = subparsers.add_parser("archive")
    archive.add_argument("output_dir", type=Path)
    archive.add_argument("archive_path", type=Path)
    database = subparsers.add_parser("database-manifest")
    database.add_argument("output_dir", type=Path)
    database.add_argument("database_prefix", type=Path)
    database.add_argument("database_snapshot_date")
    database.add_argument("--mmseqs", default="mmseqs")
    args = parser.parse_args()
    if args.command == "deduplicate":
        deduplicate(args.raw, args.out_tsv, args.out_fasta)
    elif args.command == "finalize":
        finalize(
            args.output_dir,
            args.nr_run_id,
            args.database_label,
            args.database_snapshot_date,
            args.mmseqs,
        )
    elif args.command == "database-manifest":
        write_database_manifest(
            args.output_dir,
            args.database_prefix,
            args.database_snapshot_date,
            args.mmseqs,
        )
    else:
        archive_results(args.output_dir, args.archive_path)


if __name__ == "__main__":
    main()
