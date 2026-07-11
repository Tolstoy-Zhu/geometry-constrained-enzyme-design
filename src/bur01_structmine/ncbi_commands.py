"""Generate NCBI Datasets commands for downloading gene/CDS packages.

The exact natural CDS is often obtained by following RefSeq/GeneID accessions
from UniProt. This helper creates commands that the user can run manually after
installing NCBI Datasets CLI.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Set

import pandas as pd


def _split_multi(value: str) -> List[str]:
    if not value or str(value).lower() == "nan":
        return []
    out: List[str] = []
    for part in str(value).replace(";", ",").split(","):
        x = part.strip()
        if x:
            out.append(x)
    return out


def collect_refseq_accessions(df: pd.DataFrame) -> Set[str]:
    """Collect RefSeq accessions from an annotated hit table."""
    accessions: Set[str] = set()
    if "RefSeq_accessions" not in df.columns:
        return accessions
    for value in df["RefSeq_accessions"].dropna():
        for acc in _split_multi(value):
            accessions.add(acc)
    return accessions


def collect_gene_ids(df: pd.DataFrame) -> Set[str]:
    """Collect NCBI Gene IDs from an annotated hit table."""
    gene_ids: Set[str] = set()
    if "NCBI_GeneID" not in df.columns:
        return gene_ids
    for value in df["NCBI_GeneID"].dropna():
        for gid in _split_multi(value):
            gene_ids.add(gid)
    return gene_ids


def collect_ncbi_protein_accessions(df: pd.DataFrame) -> Set[str]:
    """Collect NCBI protein accessions from an annotated hit table."""
    accessions: Set[str] = set()
    if "ncbi_protein_accession" not in df.columns:
        return accessions
    for value in df["ncbi_protein_accession"].dropna():
        for acc in _split_multi(value):
            accessions.add(acc)
    return accessions


def write_ncbi_download_script(df: pd.DataFrame, out_sh: str | Path) -> None:
    """Write a Linux/WSL shell script for NCBI Datasets downloads."""
    out_sh = Path(out_sh)
    out_sh.parent.mkdir(parents=True, exist_ok=True)

    refseqs = sorted(collect_refseq_accessions(df))
    gene_ids = sorted(collect_gene_ids(df))
    protein_accessions = sorted(collect_ncbi_protein_accessions(df))

    with out_sh.open("w", encoding="utf-8") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write("set -e\n")
        f.write("mkdir -p results/ncbi_downloads\n\n")
        f.write("# Install NCBI Datasets CLI first: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/\n")
        f.write("# Commands below are generated automatically from UniProt cross-references.\n\n")

        if gene_ids:
            f.write("# Download by NCBI Gene ID. This is usually best for CDS retrieval.\n")
            for gid in gene_ids:
                safe = gid.replace(".", "_")
                f.write(
                    f"datasets download gene gene-id {gid} --include gene,cds,protein "
                    f"--filename results/ncbi_downloads/gene_{safe}.zip\n"
                )
            f.write("\n")

        if refseqs:
            f.write("# RefSeq protein accessions. If gene-id command fails, search these in NCBI Protein.\n")
            for acc in refseqs:
                f.write(f"# RefSeq protein accession: {acc}\n")
            f.write("\n")

        if protein_accessions:
            f.write("# NCBI protein accessions from NCBI-derived Foldseek hits.\n")
            for acc in protein_accessions:
                f.write(f"# NCBI protein accession: {acc}\n")

        if not gene_ids and not refseqs and not protein_accessions:
            f.write("echo 'No NCBI GeneID or RefSeq accession was found in the annotation table.'\n")

    out_sh.chmod(0o755)
