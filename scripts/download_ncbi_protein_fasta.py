"""Download NCBI Protein FASTA records for AquaEsterNet-StructMine.

Examples:
  python scripts/download_ncbi_protein_fasta.py --query "esterase[Protein Name] AND bacteria[Organism]" --max-records 200 --out data/ncbi/ncbi_candidates.fasta
  python scripts/download_ncbi_protein_fasta.py --accessions NP_414543.1 WP_012345678.1 --out data/ncbi/ncbi_candidates.fasta
  python scripts/download_ncbi_protein_fasta.py --accession-file accessions.txt --out data/ncbi/ncbi_candidates.fasta
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Iterable, List

import requests

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _request_text(endpoint: str, params: dict, timeout: int = 60) -> str:
    response = requests.get(f"{EUTILS_BASE}/{endpoint}", params=params, timeout=timeout)
    response.raise_for_status()
    return response.text


def search_ncbi_protein(query: str, max_records: int) -> List[str]:
    """Search NCBI Protein and return internal NCBI IDs."""
    text = _request_text(
        "esearch.fcgi",
        {
            "db": "protein",
            "term": query,
            "retmode": "json",
            "retmax": str(max_records),
            "sort": "relevance",
        },
    )
    data = json.loads(text)
    return data.get("esearchresult", {}).get("idlist", [])


def read_accessions(path: Path) -> List[str]:
    """Read one accession per line, ignoring blanks and comments."""
    out: List[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            out.append(line)
    return out


def batched(values: List[str], batch_size: int) -> Iterable[List[str]]:
    for i in range(0, len(values), batch_size):
        yield values[i:i + batch_size]


def fetch_fasta(ids_or_accessions: List[str], batch_size: int = 100, sleep_seconds: float = 0.34) -> str:
    """Fetch FASTA records from NCBI Protein."""
    chunks: List[str] = []
    for batch in batched(ids_or_accessions, batch_size):
        text = _request_text(
            "efetch.fcgi",
            {
                "db": "protein",
                "id": ",".join(batch),
                "rettype": "fasta",
                "retmode": "text",
            },
        )
        chunks.append(text.strip())
        time.sleep(sleep_seconds)
    return "\n".join(chunk for chunk in chunks if chunk)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Download NCBI Protein FASTA records")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--query", help="NCBI Protein search query")
    source.add_argument("--accessions", nargs="+", help="NCBI protein accessions")
    source.add_argument("--accession-file", type=Path, help="Text file with one accession per line")
    parser.add_argument("--max-records", type=int, default=200, help="Maximum records for --query")
    parser.add_argument("--out", type=Path, default=Path("data/ncbi/ncbi_candidates.fasta"), help="Output FASTA path")
    parser.add_argument("--batch-size", type=int, default=100, help="NCBI efetch batch size")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.query:
        ids = search_ncbi_protein(args.query, args.max_records)
    elif args.accession_file:
        ids = read_accessions(args.accession_file)
    else:
        ids = args.accessions

    if not ids:
        raise SystemExit("No NCBI Protein records found.")

    fasta = fetch_fasta(ids, batch_size=args.batch_size)
    if not fasta.startswith(">"):
        raise SystemExit("NCBI did not return FASTA content.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(fasta + "\n", encoding="utf-8")
    print(f"[OK] Wrote {len(ids)} NCBI Protein records to {args.out}")


if __name__ == "__main__":
    try:
        main()
    except requests.RequestException as exc:
        print(f"NCBI request failed: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
