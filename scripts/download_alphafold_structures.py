"""Download AlphaFold DB structures for a small local Foldseek database.

Examples:
  python scripts/download_alphafold_structures.py --accessions Q6RYA0 Q9FVW3 --format pdb --out-dir data/alphafold/structures
  python scripts/download_alphafold_structures.py --accession-file accessions.txt --format cif --out-dir data/alphafold/structures
  python scripts/download_alphafold_structures.py --uniprot-query "reviewed:true AND (protein_name:esterase OR protein_name:hydrolase)" --max-records 50 --format pdb
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Iterable, List

import requests

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from bur01_structmine.alphafold_client import AlphaFoldClient  # noqa: E402


def read_accessions(path: Path) -> List[str]:
    out: List[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            out.append(line)
    return out


def search_uniprot_accessions(query: str, max_records: int) -> List[str]:
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "fields": "accession",
        "format": "tsv",
        "size": str(max_records),
    }
    response = requests.get(url, params=params, timeout=60)
    response.raise_for_status()
    accessions: List[str] = []
    for line in response.text.splitlines()[1:]:
        acc = line.strip().split("\t")[0]
        if acc:
            accessions.append(acc)
    return accessions


def choose_url(record, file_format: str) -> str:
    if file_format == "pdb":
        return record.pdb_url
    if file_format == "cif":
        return record.cif_url
    if file_format == "bcif":
        return record.bcif_url
    raise ValueError(file_format)


def download_file(url: str, out_file: Path, timeout: int = 120) -> None:
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_bytes(response.content)


def write_manifest(rows: Iterable[dict], path: Path) -> None:
    columns = [
        "accession",
        "model_entity_id",
        "global_plddt",
        "is_reviewed",
        "organism",
        "protein_description",
        "structure_file",
        "pdb_url",
        "cif_url",
        "pae_url",
    ]
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(str(row.get(col, "")).replace("\t", " ") for col in columns))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download AlphaFold DB structures")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--accessions", nargs="+", help="UniProt accessions")
    source.add_argument("--accession-file", type=Path, help="Text file with one UniProt accession per line")
    source.add_argument("--uniprot-query", help="UniProt search query, e.g. reviewed:true AND protein_name:esterase")
    parser.add_argument("--max-records", type=int, default=50, help="Maximum UniProt records for --uniprot-query")
    parser.add_argument("--format", choices=["pdb", "cif", "bcif"], default="pdb", help="Structure file format")
    parser.add_argument("--out-dir", type=Path, default=Path("data/alphafold/structures"), help="Output structure directory")
    parser.add_argument("--manifest", type=Path, default=Path("data/alphafold/alphafold_manifest.tsv"), help="Output manifest TSV")
    parser.add_argument("--sleep-seconds", type=float, default=0.15, help="Delay between AlphaFold API calls")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.accessions:
        accessions = args.accessions
    elif args.accession_file:
        accessions = read_accessions(args.accession_file)
    else:
        accessions = search_uniprot_accessions(args.uniprot_query, args.max_records)

    if not accessions:
        raise SystemExit("No UniProt accessions were found.")

    client = AlphaFoldClient(sleep_seconds=args.sleep_seconds)
    rows = []
    downloaded = 0
    suffix = {"pdb": ".pdb", "cif": ".cif", "bcif": ".bcif"}[args.format]
    for accession in accessions:
        record = client.fetch(accession)
        if not record:
            print(f"[WARN] No AlphaFold DB record for {accession}", file=sys.stderr)
            continue
        url = choose_url(record, args.format)
        if not url:
            print(f"[WARN] No {args.format} URL for {accession}", file=sys.stderr)
            continue
        out_file = args.out_dir / f"{record.model_entity_id or accession}{suffix}"
        download_file(url, out_file)
        rows.append(
            {
                "accession": accession,
                "model_entity_id": record.model_entity_id,
                "global_plddt": record.global_plddt,
                "is_reviewed": record.is_reviewed,
                "organism": record.organism,
                "protein_description": record.protein_description,
                "structure_file": out_file.as_posix(),
                "pdb_url": record.pdb_url,
                "cif_url": record.cif_url,
                "pae_url": record.pae_url,
            }
        )
        downloaded += 1
        time.sleep(args.sleep_seconds)

    if downloaded == 0:
        raise SystemExit("No AlphaFold structures were downloaded.")
    write_manifest(rows, args.manifest)
    print(f"[OK] Downloaded {downloaded} AlphaFold structure file(s) to {args.out_dir}")
    print(f"[OK] Manifest: {args.manifest}")


if __name__ == "__main__":
    try:
        main()
    except requests.RequestException as exc:
        print(f"Download failed: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
