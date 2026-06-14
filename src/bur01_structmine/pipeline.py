"""Main Bur01-StructMine V1 pipeline."""
from __future__ import annotations

from pathlib import Path

from .annotator import annotate_hits, export_annotation_outputs
from .config import load_config, resolve_path
from .foldseek_runner import run_easy_search
from .io_utils import ensure_dir, write_excel_safely, write_fasta
from .ncbi_commands import write_ncbi_download_script
from .pdb_utils import extract_sequence_from_pdb
from .results_parser import filter_hits, load_foldseek_tsv
from .scoring import add_candidate_scores


def export_query_sequence(query_pdb: Path, out_fasta: Path) -> None:
    """Extract query PDB sequence and write FASTA."""
    seqs = extract_sequence_from_pdb(query_pdb)
    records = [(f"{query_pdb.stem}_chain_{chain}", seq) for chain, seq in seqs.items()]
    write_fasta(records, out_fasta)


def run_pipeline(config_path: str = "config.yaml", skip_foldseek: bool = False, input_tsv: str | None = None) -> None:
    """Run the full V1 pipeline.

    Parameters
    ----------
    config_path:
        YAML configuration file.
    skip_foldseek:
        If True, skip Foldseek execution and parse an existing TSV file.
    input_tsv:
        Existing Foldseek TSV file. Required when skip_foldseek=True.
    """
    cfg = load_config(config_path)

    query_pdb = resolve_path(cfg["query"]["pdb_file"])
    out_prefix = resolve_path(cfg["output"]["out_prefix"])
    ensure_dir(out_prefix.parent)

    # Step 0: export query sequence, useful for records and debugging.
    export_query_sequence(query_pdb, out_prefix.with_suffix(".query_from_pdb.fasta"))

    # Step 1: run Foldseek or use an existing TSV.
    raw_tsv = out_prefix.with_suffix(".foldseek.raw.tsv")
    if skip_foldseek:
        if not input_tsv:
            raise ValueError("input_tsv is required when skip_foldseek=True")
        raw_tsv = resolve_path(input_tsv)
        print(f"[INFO] Skip Foldseek. Use existing TSV: {raw_tsv}")
    else:
        fs_cfg = cfg["foldseek"]
        run_easy_search(
            query_pdb=query_pdb,
            target_db=resolve_path(fs_cfg["target_db"]),
            out_tsv=raw_tsv,
            tmp_dir=resolve_path(fs_cfg.get("tmp_dir", "tmp")),
            executable=fs_cfg.get("executable", "foldseek"),
            max_seqs=int(fs_cfg.get("max_seqs", 1000)),
            threads=int(fs_cfg.get("threads", 8)),
        )

    # Step 2: parse and filter hits.
    raw_df = load_foldseek_tsv(raw_tsv)
    write_excel_safely(raw_df, out_prefix.with_suffix(".foldseek.raw.xlsx"))

    flt_cfg = cfg.get("filter", {})
    filtered_df = filter_hits(
        raw_df,
        min_alntmscore=float(flt_cfg.get("min_alntmscore", 0.45)),
        max_rmsd=float(flt_cfg.get("max_rmsd", 4.0)),
        min_qcov=float(flt_cfg.get("min_qcov", 0.60)),
        min_tcov=float(flt_cfg.get("min_tcov", 0.40)),
        min_alnlen=int(flt_cfg.get("min_alnlen", 180)),
        max_sequence_identity=float(flt_cfg.get("max_sequence_identity", 0.95)),
    )
    filtered_tsv = out_prefix.with_suffix(".foldseek.filtered.tsv")
    filtered_df.to_csv(filtered_tsv, sep="\t", index=False)
    write_excel_safely(filtered_df, out_prefix.with_suffix(".foldseek.filtered.xlsx"))

    # Step 3: annotate hits with RCSB and UniProt.
    ann_cfg = cfg.get("annotation", {})
    annotated_df = annotate_hits(
        filtered_df,
        enable_rcsb=bool(ann_cfg.get("enable_rcsb", True)),
        enable_uniprot=bool(ann_cfg.get("enable_uniprot", True)),
        request_sleep_seconds=float(ann_cfg.get("request_sleep_seconds", 0.15)),
    )

    # Step 4: add simple candidate scores and export final files.
    annotated_df = add_candidate_scores(annotated_df)
    export_annotation_outputs(annotated_df, out_prefix)

    # Step 5: write NCBI download helper script.
    write_ncbi_download_script(annotated_df, out_prefix.with_suffix(".download_ncbi_cds.sh"))

    print("\n[DONE] Bur01-StructMine V1 pipeline finished.")
    print(f"Raw table:        {out_prefix}.foldseek.raw.xlsx")
    print(f"Filtered table:   {out_prefix}.foldseek.filtered.xlsx")
    print(f"Annotated table:  {out_prefix}.annotated.xlsx")
    print(f"Protein FASTA:    {out_prefix}.candidate_proteins.fasta")
    print(f"NCBI script:      {out_prefix}.download_ncbi_cds.sh")
