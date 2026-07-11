"""Validate real AquaEsterNet-StructMine search outputs after a Foldseek run."""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd


BASE_COLUMNS = {
    "source_db",
    "source_kind",
    "annotation_source",
    "target",
    "target_uniprot_accession",
    "ncbi_protein_accession",
    "alphafold_model_entity_id",
    "alphafold_pdb_url",
    "alphafold_pae_url",
    "final_protein_sequence",
    "synthesis_recommendation",
}


def fail(message: str) -> None:
    raise SystemExit(message)


def nonempty_count(series: pd.Series) -> int:
    return int(series.fillna("").astype(str).str.strip().ne("").sum())


def validate_required_columns(df: pd.DataFrame, columns: Iterable[str]) -> None:
    missing = sorted(set(columns) - set(df.columns))
    if missing:
        fail(f"Missing required column(s): {', '.join(missing)}")


def validate_source(df: pd.DataFrame, source: str, allow_missing_alphafold_metadata: bool = False) -> None:
    subset = df[df["source_db"].fillna("").astype(str) == source]
    if subset.empty:
        observed = ", ".join(sorted(df["source_db"].fillna("").astype(str).unique()))
        fail(f"Required source_db '{source}' was not found. Observed: {observed}")

    annotation_sources = set(subset["annotation_source"].fillna("").astype(str))
    if "alphafold" in annotation_sources or "swissprot" in annotation_sources:
        if nonempty_count(subset["target_uniprot_accession"]) == 0:
            fail(f"Source '{source}' has no parsed UniProt accession.")
        if (
            not allow_missing_alphafold_metadata
            and nonempty_count(subset["alphafold_model_entity_id"]) == 0
            and nonempty_count(subset["alphafold_pdb_url"]) == 0
        ):
            fail(f"Source '{source}' has no AlphaFold model metadata.")

    if "ncbi" in annotation_sources:
        has_ncbi_accession = nonempty_count(subset["ncbi_protein_accession"]) > 0
        has_sequence = nonempty_count(subset["final_protein_sequence"]) > 0
        if not has_ncbi_accession and not has_sequence:
            fail(f"Source '{source}' has neither NCBI protein accessions nor final sequences.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate real Foldseek search outputs")
    parser.add_argument(
        "--annotated",
        type=Path,
        default=Path("results/Bur01_vs_database.annotated.xlsx"),
        help="Annotated Excel workbook from run_pipeline.py",
    )
    parser.add_argument("--require-source", action="append", default=[], help="Required source_db value")
    parser.add_argument("--min-rows", type=int, default=1, help="Minimum annotated row count")
    parser.add_argument(
        "--allow-missing-alphafold-metadata",
        action="store_true",
        help="Do not fail AlphaFold/Swiss-Prot sources without AlphaFold API metadata",
    )
    parser.add_argument(
        "--require-final-sequence",
        action="store_true",
        help="Require at least one non-empty final_protein_sequence value",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.annotated.exists():
        fail(f"Annotated workbook not found: {args.annotated}")

    df = pd.read_excel(args.annotated)
    validate_required_columns(df, BASE_COLUMNS)
    if len(df) < args.min_rows:
        fail(f"Annotated workbook has {len(df)} row(s), expected at least {args.min_rows}.")

    for source in args.require_source:
        validate_source(df, source, allow_missing_alphafold_metadata=args.allow_missing_alphafold_metadata)

    if args.require_final_sequence and nonempty_count(df["final_protein_sequence"]) == 0:
        fail("No final_protein_sequence values were found.")

    print(f"[OK] Validated {len(df)} annotated row(s) from {args.annotated}")
    if args.require_source:
        print(f"[OK] Required sources present: {', '.join(args.require_source)}")


if __name__ == "__main__":
    main()
