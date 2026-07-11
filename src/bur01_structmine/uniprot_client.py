"""Cached UniProtKB REST API client with batch retrieval."""
from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional

from .config import project_root
from .http_utils import CachedSession


@dataclass
class UniProtRecord:
    accession: str = ""
    entry_type: str = ""
    is_reviewed: str = ""
    protein_name: str = ""
    organism: str = ""
    gene_names: str = ""
    protein_sequence: str = ""
    embl_accessions: str = ""
    refseq_accessions: str = ""
    geneid: str = ""


class UniProtClient:
    """Retrieve reviewed status, annotation, sequence, and cross-references."""

    def __init__(
        self,
        sleep_seconds: float = 0.15,
        timeout: int = 30,
        cache_dir: str | Path | None = None,
    ) -> None:
        self.sleep_seconds = sleep_seconds
        self._cache: Dict[str, Optional[UniProtRecord]] = {}
        self.errors: Dict[str, str] = {}
        self.http = CachedSession(
            cache_dir or project_root() / ".cache" / "annotations" / "uniprot",
            timeout=timeout,
        )

    @staticmethod
    def _extract_recommended_name(data: dict) -> str:
        desc = data.get("proteinDescription", {}) or {}
        rec = desc.get("recommendedName", {}) or {}
        full = rec.get("fullName", {}) or {}
        value = full.get("value", "")
        if value:
            return value
        sub = desc.get("submissionNames", []) or []
        return sub[0].get("fullName", {}).get("value", "") if sub else ""

    @staticmethod
    def _extract_gene_names(data: dict) -> str:
        names = []
        for gene in data.get("genes", []) or []:
            if gene.get("geneName", {}).get("value"):
                names.append(gene["geneName"]["value"])
            names.extend(
                synonym["value"]
                for synonym in gene.get("synonyms", []) or []
                if synonym.get("value")
            )
        return ",".join(sorted(set(names)))

    @staticmethod
    def _extract_crossrefs(data: dict) -> Dict[str, str]:
        out = {"EMBL": [], "RefSeq": [], "GeneID": []}
        for ref in data.get("uniProtKBCrossReferences", []) or []:
            database = ref.get("database", "")
            identifier = ref.get("id", "")
            if database in out and identifier:
                out[database].append(identifier)
        return {key: ",".join(sorted(set(values))) for key, values in out.items()}

    def _to_record(self, data: dict) -> UniProtRecord:
        accession = str(data.get("primaryAccession", ""))
        xrefs = self._extract_crossrefs(data)
        entry_type = str(data.get("entryType", ""))
        return UniProtRecord(
            accession=accession,
            entry_type=entry_type,
            is_reviewed=str("reviewed" in entry_type.lower()),
            protein_name=self._extract_recommended_name(data),
            organism=(data.get("organism", {}) or {}).get("scientificName", ""),
            gene_names=self._extract_gene_names(data),
            protein_sequence=(data.get("sequence", {}) or {}).get("value", ""),
            embl_accessions=xrefs.get("EMBL", ""),
            refseq_accessions=xrefs.get("RefSeq", ""),
            geneid=xrefs.get("GeneID", ""),
        )

    def fetch_many(self, accessions: Iterable[str], batch_size: int = 100) -> Dict[str, Optional[UniProtRecord]]:
        normalized = list(dict.fromkeys(str(accession).strip().upper() for accession in accessions if str(accession).strip()))
        missing = [accession for accession in normalized if accession not in self._cache]
        for start in range(0, len(missing), batch_size):
            batch = missing[start:start + batch_size]
            query = " OR ".join(f"(accession:{accession})" for accession in batch)
            data = self.http.get_json(
                "https://rest.uniprot.org/uniprotkb/search",
                params={"query": query, "format": "json", "size": len(batch)},
            )
            observed: set[str] = set()
            if isinstance(data, dict):
                for item in data.get("results", []) or []:
                    record = self._to_record(item)
                    if record.accession:
                        self._cache[record.accession] = record
                        observed.add(record.accession)
            for accession in batch:
                if accession not in observed:
                    self._cache[accession] = None
                    self.errors[accession] = self.http.last_error or "UniProt record was not returned"
            if batch and self.sleep_seconds:
                time.sleep(self.sleep_seconds)
        return {accession: self._cache.get(accession) for accession in normalized}

    def fetch(self, accession: str) -> Optional[UniProtRecord]:
        accession = accession.strip().upper()
        if not accession:
            return None
        return self.fetch_many([accession]).get(accession)
