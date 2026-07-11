"""Write a reproducible metadata manifest for a configured Foldseek database."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from bur01_structmine.config import load_config, project_root
from bur01_structmine.datasources import load_search_databases
from bur01_structmine.foldseek_runner import get_foldseek_version
from bur01_structmine.run_context import file_sha256


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/alphafold_swissprot.yaml")
    parser.add_argument("--source", required=True)
    parser.add_argument("--full-hash", action="store_true", help="Hash every database file, including large payloads")
    parser.add_argument("--out")
    args = parser.parse_args()

    cfg = load_config(args.config)
    sources = {source.name: source for source in load_search_databases(cfg)}
    if args.source not in sources:
        raise SystemExit(f"Configured source not found: {args.source}")
    source = sources[args.source]
    if source.target_db is None:
        raise SystemExit(f"Source has no Foldseek target_db: {source.name}")
    files = sorted(path for path in source.target_db.parent.glob(source.target_db.name + "*") if path.is_file())
    if not files:
        raise SystemExit(f"No database files found for prefix: {source.target_db}")

    metadata_digest = hashlib.sha256()
    records = []
    total_bytes = 0
    for path in files:
        stat = path.stat()
        total_bytes += stat.st_size
        metadata_line = f"{path.name}\0{stat.st_size}\0{stat.st_mtime_ns}\n"
        metadata_digest.update(metadata_line.encode("utf-8"))
        should_hash = args.full_hash or stat.st_size <= 16 * 1024 * 1024
        records.append(
            {
                "name": path.name,
                "size_bytes": stat.st_size,
                "modified_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
                "sha256": file_sha256(path) if should_hash else None,
            }
        )

    fs_cfg = cfg.get("foldseek", {})
    foldseek_version = get_foldseek_version(
        str(fs_cfg.get("executable", "foldseek")),
        backend=str(fs_cfg.get("backend", "auto")),
        wsl_distribution=str(fs_cfg.get("wsl_distribution", "Ubuntu-Foldseek")),
        wsl_executable=str(fs_cfg.get("wsl_executable", "/opt/foldseek/bin/foldseek")),
    )
    version_path = Path(f"{source.target_db}.version")
    database_version = version_path.read_text(encoding="utf-8").strip() if version_path.is_file() else "unknown"
    provenance_path = (
        project_root()
        / "databases"
        / "environment"
        / f"{source.target_db.name}_download_provenance.json"
    )
    manifest = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "source_name": source.name,
        "database_reference": source.foldseek_database,
        "target_prefix": str(source.target_db.resolve()),
        "foldseek_version": foldseek_version,
        "database_version": " ".join(database_version.split()),
        "download_command": f"foldseek databases {source.foldseek_database} <target> <tmp>",
        "archive_url": source.archive_url or None,
        "archive_size": source.archive_size,
        "archive_sha256": source.archive_sha256 or None,
        "download_provenance": str(provenance_path.resolve()) if provenance_path.is_file() else None,
        "file_count": len(records),
        "total_bytes": total_bytes,
        "metadata_sha256": metadata_digest.hexdigest(),
        "large_file_hash_policy": "full" if args.full_hash else "files <= 16 MiB",
        "files": records,
    }
    out = Path(args.out) if args.out else project_root() / "databases" / "environment" / f"{source.safe_name}_database_manifest.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    tmp.replace(out)
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
