"""Annotate Foldseek hits by RCSB and UniProt."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List

import pandas as pd

from .io_utils import write_excel_safely, write_fasta
from .rcsb_client import RcsbClient
from .results_parser import parse_pdb_chain, parse_uniprot_from_afdb
from .uniprot_client import UniProtClient


def _split_accessions(acc_string: str) -> List[str]:
    if not acc_string or str(acc_string).lower() == "nan":
        return []
    return [x.strip() for x in str(acc_string).replace(";", ",").split(",") if x.strip()]


def annotate_hits(
    df: pd.DataFrame,
    enable_rcsb: bool = True,
    enable_uniprot: bool = True,
    request_sleep_seconds: float = 0.15,
) -> pd.DataFrame:
    """Annotate Foldseek hits.

    For PDB hits: target -> PDB ID/chain -> RCSB -> UniProt -> NCBI refs.
    For AlphaFold DB hits: target -> UniProt accession -> UniProt -> NCBI refs.
    """
    rcsb = RcsbClient(sleep_seconds=request_sleep_seconds) if enable_rcsb else None
    uniprot = UniProtClient(sleep_seconds=request_sleep_seconds) if enable_uniprot else None

    records: List[Dict] = []
    for _, row in df.iterrows():
        rec: Dict = row.to_dict()
        target = str(row.get("target", ""))

        pdb_id, chain = parse_pdb_chain(target)
        af_uniprot = parse_uniprot_from_afdb(target)

        rec["pdb_id"] = pdb_id or ""
        rec["chain"] = chain or ""
        rec["afdb_uniprot_accession"] = af_uniprot or ""

        # RCSB annotation for PDB hits.
        if enable_rcsb and rcsb and pdb_id:
            entities = rcsb.get_polymer_entities(pdb_id)
            ent = rcsb.choose_entity_for_chain(entities, chain)
            if ent:
                rec["rcsb_entity_id"] = ent.entity_id
                rec["rcsb_protein_name"] = ent.protein_name
                rec["rcsb_organism"] = ent.organism
                rec["rcsb_protein_sequence"] = ent.protein_sequence
                rec["uniprot_accessions"] = ent.uniprot_accessions
            else:
                rec.setdefault("uniprot_accessions", "")
        else:
            rec.setdefault("uniprot_accessions", "")

        # AlphaFold DB IDs usually already contain UniProt accession.
        if af_uniprot and not rec.get("uniprot_accessions"):
            rec["uniprot_accessions"] = af_uniprot

        # UniProt annotation.
        all_accs = _split_accessions(rec.get("uniprot_accessions", ""))
        up_records = []
        if enable_uniprot and uniprot:
            for acc in all_accs:
                up = uniprot.fetch(acc)
                if up:
                    up_records.append(up)

        if up_records:
            rec["uniprot_protein_name"] = "; ".join(sorted({x.protein_name for x in up_records if x.protein_name}))
            rec["uniprot_organism"] = "; ".join(sorted({x.organism for x in up_records if x.organism}))
            rec["gene_names"] = "; ".join(sorted({x.gene_names for x in up_records if x.gene_names}))
            rec["uniprot_protein_sequence"] = up_records[0].protein_sequence
            rec["EMBL_accessions"] = "; ".join(sorted({x.embl_accessions for x in up_records if x.embl_accessions}))
            rec["RefSeq_accessions"] = "; ".join(sorted({x.refseq_accessions for x in up_records if x.refseq_accessions}))
            rec["NCBI_GeneID"] = "; ".join(sorted({x.geneid for x in up_records if x.geneid}))
        else:
            rec.setdefault("uniprot_protein_name", "")
            rec.setdefault("uniprot_organism", "")
            rec.setdefault("gene_names", "")
            rec.setdefault("uniprot_protein_sequence", "")
            rec.setdefault("EMBL_accessions", "")
            rec.setdefault("RefSeq_accessions", "")
            rec.setdefault("NCBI_GeneID", "")

        # Choose final protein sequence with priority: UniProt > RCSB > empty.
        final_seq = rec.get("uniprot_protein_sequence") or rec.get("rcsb_protein_sequence") or ""
        rec["final_protein_sequence"] = final_seq
        rec["final_protein_length"] = len(final_seq) if final_seq else ""
        rec["has_unknown_X"] = "X" in final_seq if final_seq else ""

        records.append(rec)

    return pd.DataFrame(records)


def export_annotation_outputs(df: pd.DataFrame, out_prefix: str | Path) -> None:
    """Export annotated Excel and protein FASTA."""
    out_prefix = Path(out_prefix)
    write_excel_safely(df, f"{out_prefix}.annotated.xlsx")

    fasta_records = []
    for _, row in df.iterrows():
        seq = row.get("final_protein_sequence", "")
        if not seq or str(seq).lower() == "nan":
            continue
        rec_id = row.get("uniprot_accessions") or row.get("target")
        desc = row.get("uniprot_protein_name") or row.get("rcsb_protein_name") or "candidate"
        fasta_records.append((f"{rec_id} {desc}".strip(), seq))
    write_fasta(fasta_records, f"{out_prefix}.candidate_proteins.fasta")
