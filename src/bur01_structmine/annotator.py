"""Cross-database annotation with batching, cache, and explicit failure state."""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List

import pandas as pd

from .alphafold_client import AlphaFoldClient, AlphaFoldRecord
from .config import project_root
from .io_utils import write_excel_safely, write_fasta
from .ncbi_client import NcbiProteinClient, NcbiProteinRecord
from .rcsb_client import RcsbClient
from .results_parser import parse_ncbi_protein_accession, parse_pdb_chain, parse_uniprot_accession
from .uniprot_client import UniProtClient, UniProtRecord

ALPHAFOLD_COLUMNS = (
    "alphafold_model_entity_id",
    "alphafold_latest_version",
    "alphafold_model_created_date",
    "alphafold_sequence_version_date",
    "alphafold_global_plddt",
    "alphafold_is_reviewed",
    "alphafold_organism",
    "alphafold_protein_description",
    "alphafold_pdb_url",
    "alphafold_cif_url",
    "alphafold_bcif_url",
    "alphafold_pae_url",
    "alphafold_plddt_url",
)


def _split_accessions(value: Any) -> List[str]:
    if value is None or str(value).strip().lower() in {"", "nan", "none"}:
        return []
    return [part.strip() for part in str(value).replace(";", ",").split(",") if part.strip()]


def _join(values: Iterable[str]) -> str:
    return "; ".join(sorted({str(value) for value in values if str(value).strip()}))


def _set_alphafold(rec: Dict[str, Any], record: AlphaFoldRecord | None) -> None:
    values = {
        "alphafold_model_entity_id": record.model_entity_id if record else "",
        "alphafold_latest_version": record.latest_version if record else "",
        "alphafold_model_created_date": record.model_created_date if record else "",
        "alphafold_sequence_version_date": record.sequence_version_date if record else "",
        "alphafold_global_plddt": record.global_plddt if record else "",
        "alphafold_is_reviewed": record.is_reviewed if record else "",
        "alphafold_organism": record.organism if record else "",
        "alphafold_protein_description": record.protein_description if record else "",
        "alphafold_pdb_url": record.pdb_url if record else "",
        "alphafold_cif_url": record.cif_url if record else "",
        "alphafold_bcif_url": record.bcif_url if record else "",
        "alphafold_pae_url": record.pae_url if record else "",
        "alphafold_plddt_url": record.plddt_url if record else "",
    }
    rec.update(values)


def _set_ncbi(rec: Dict[str, Any], record: NcbiProteinRecord | None) -> None:
    rec["ncbi_protein_name"] = record.protein_name if record else ""
    rec["ncbi_organism"] = record.organism if record else ""
    rec["ncbi_protein_sequence"] = record.protein_sequence if record else ""
    rec["ncbi_fasta_header"] = record.fasta_header if record else ""


def _set_uniprot(rec: Dict[str, Any], records: list[UniProtRecord]) -> None:
    if not records:
        rec.update(
            {
                "uniprot_protein_name": "",
                "uniprot_organism": "",
                "gene_names": "",
                "uniprot_entry_type": "",
                "uniprot_is_reviewed": "",
                "uniprot_protein_sequence": "",
                "EMBL_accessions": "",
                "RefSeq_accessions": "",
                "NCBI_GeneID": "",
            }
        )
        return
    rec["uniprot_protein_name"] = _join(record.protein_name for record in records)
    rec["uniprot_organism"] = _join(record.organism for record in records)
    rec["gene_names"] = _join(record.gene_names for record in records)
    rec["uniprot_entry_type"] = _join(record.entry_type for record in records)
    rec["uniprot_is_reviewed"] = all(record.is_reviewed.lower() == "true" for record in records)
    rec["uniprot_protein_sequence"] = records[0].protein_sequence
    rec["EMBL_accessions"] = _join(record.embl_accessions for record in records)
    rec["RefSeq_accessions"] = _join(record.refseq_accessions for record in records)
    rec["NCBI_GeneID"] = _join(record.geneid for record in records)


def annotate_hits(
    df: pd.DataFrame,
    enable_rcsb: bool = True,
    enable_uniprot: bool = True,
    enable_ncbi: bool = True,
    enable_alphafold: bool = True,
    request_sleep_seconds: float = 0.15,
    cache_dir: str | Path | None = None,
) -> pd.DataFrame:
    """Annotate search hits and preserve per-row error provenance."""
    root_cache = Path(cache_dir) if cache_dir else project_root() / ".cache" / "annotations"
    rcsb = RcsbClient(request_sleep_seconds, cache_dir=root_cache / "rcsb") if enable_rcsb else None
    uniprot = UniProtClient(request_sleep_seconds, cache_dir=root_cache / "uniprot") if enable_uniprot else None
    ncbi = NcbiProteinClient(max(request_sleep_seconds, 0.34), cache_dir=root_cache / "ncbi") if enable_ncbi else None
    alphafold = AlphaFoldClient(request_sleep_seconds, cache_dir=root_cache / "alphafold") if enable_alphafold else None

    records: list[Dict[str, Any]] = []
    for _, row in df.iterrows():
        rec = row.to_dict()
        target = str(row.get("target", ""))
        annotation_source = str(row.get("annotation_source", "auto") or "auto").lower()
        target_uniprot = parse_uniprot_accession(target)
        target_ncbi = parse_ncbi_protein_accession(target)
        pdb_id, chain = (None, None)
        if annotation_source in {"pdb", "auto"} and not target_uniprot:
            pdb_id, chain = parse_pdb_chain(target)
        rec["pdb_id"] = pdb_id or ""
        rec["chain"] = chain or rec.get("target_chain", "") or ""
        rec["target_uniprot_accession"] = target_uniprot or ""
        rec["ncbi_protein_accession"] = target_ncbi or ""
        rec["uniprot_accessions"] = rec.get("uniprot_accessions", "") or ""
        rec["_annotation_errors"] = []

        if rcsb and pdb_id:
            entities = rcsb.get_polymer_entities(pdb_id)
            entity = rcsb.choose_entity_for_chain(entities, chain)
            if entity:
                rec["rcsb_entity_id"] = entity.entity_id
                rec["rcsb_protein_name"] = entity.protein_name
                rec["rcsb_organism"] = entity.organism
                rec["rcsb_protein_sequence"] = entity.protein_sequence
                rec["uniprot_accessions"] = entity.uniprot_accessions
            else:
                rec["_annotation_errors"].append(
                    f"RCSB {pdb_id}: {next(iter(rcsb.errors.values()), 'polymer entity not returned')}"
                )
        rec.setdefault("rcsb_entity_id", "")
        rec.setdefault("rcsb_protein_name", "")
        rec.setdefault("rcsb_organism", "")
        rec.setdefault("rcsb_protein_sequence", "")
        if annotation_source in {"alphafold", "swissprot", "uniprot", "auto"} and target_uniprot:
            rec["uniprot_accessions"] = rec["uniprot_accessions"] or target_uniprot
        records.append(rec)

    all_uniprot = [accession for rec in records for accession in _split_accessions(rec["uniprot_accessions"])]
    all_ncbi = [str(rec["ncbi_protein_accession"]) for rec in records if rec["ncbi_protein_accession"]]
    all_alphafold = [
        (_split_accessions(rec["uniprot_accessions"]) or [rec["target_uniprot_accession"]])[0]
        for rec in records
        if str(rec.get("annotation_source", "auto")).lower() in {"alphafold", "swissprot", "uniprot", "auto"}
        and (_split_accessions(rec["uniprot_accessions"]) or [rec["target_uniprot_accession"]])[0]
    ]
    uniprot_records = uniprot.fetch_many(all_uniprot) if uniprot else {}
    ncbi_records = ncbi.fetch_many(all_ncbi) if ncbi else {}
    alphafold_records = alphafold.fetch_many(all_alphafold) if alphafold else {}

    for rec in records:
        accessions = _split_accessions(rec["uniprot_accessions"])
        selected_uniprot = [uniprot_records[accession] for accession in accessions if uniprot_records.get(accession)]
        _set_uniprot(rec, selected_uniprot)
        if uniprot:
            rec["_annotation_errors"].extend(
                f"UniProt {accession}: {uniprot.errors[accession]}"
                for accession in accessions
                if accession in uniprot.errors
            )

        af_accession = (accessions or [rec.get("target_uniprot_accession", "")])[0]
        af_record = alphafold_records.get(af_accession) if alphafold else None
        _set_alphafold(rec, af_record)
        if alphafold and af_accession in alphafold.errors:
            rec["_annotation_errors"].append(f"AlphaFold {af_accession}: {alphafold.errors[af_accession]}")

        ncbi_accession = str(rec.get("ncbi_protein_accession", ""))
        _set_ncbi(rec, ncbi_records.get(ncbi_accession) if ncbi else None)
        if ncbi and ncbi_accession in ncbi.errors:
            rec["_annotation_errors"].append(f"NCBI {ncbi_accession}: {ncbi.errors[ncbi_accession]}")

        final_sources = (
            ("UniProt", rec.get("uniprot_protein_sequence", "")),
            ("RCSB", rec.get("rcsb_protein_sequence", "")),
            ("NCBI", rec.get("ncbi_protein_sequence", "")),
        )
        source, final_sequence = next(
            ((name, sequence) for name, sequence in final_sources if sequence),
            ("", ""),
        )
        rec["final_protein_sequence"] = final_sequence
        rec["final_sequence_source"] = source
        rec["final_protein_length"] = len(final_sequence) if final_sequence else ""
        rec["has_unknown_X"] = "X" in final_sequence if final_sequence else ""
        errors = list(dict.fromkeys(rec.pop("_annotation_errors")))
        rec["annotation_errors"] = "; ".join(errors)
        if not any((enable_rcsb, enable_uniprot, enable_ncbi, enable_alphafold)):
            rec["annotation_status"] = "disabled"
        elif errors and final_sequence:
            rec["annotation_status"] = "partial"
        elif errors:
            rec["annotation_status"] = "failed"
        else:
            rec["annotation_status"] = "complete"

    return pd.DataFrame(records)


def export_annotation_outputs(df: pd.DataFrame, out_prefix: str | Path) -> None:
    """Export annotated Excel and deduplicated candidate protein FASTA."""
    out_prefix = Path(out_prefix)
    write_excel_safely(df, f"{out_prefix}.annotated.xlsx")
    fasta_records = []
    seen: set[tuple[str, str]] = set()
    for _, row in df.iterrows():
        sequence = row.get("final_protein_sequence", "")
        if not sequence or str(sequence).lower() == "nan":
            continue
        record_id = row.get("uniprot_accessions") or row.get("ncbi_protein_accession") or row.get("target")
        description = (
            row.get("uniprot_protein_name")
            or row.get("rcsb_protein_name")
            or row.get("ncbi_protein_name")
            or "candidate"
        )
        key = (str(record_id), str(sequence))
        if key in seen:
            continue
        seen.add(key)
        fasta_records.append((f"{record_id} {description}".strip(), str(sequence)))
    write_fasta(fasta_records, f"{out_prefix}.candidate_proteins.fasta")
