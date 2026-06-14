"""UniProt REST API client."""
from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Dict, Optional

import requests


@dataclass
class UniProtRecord:
    accession: str = ""
    protein_name: str = ""
    organism: str = ""
    gene_names: str = ""
    protein_sequence: str = ""
    embl_accessions: str = ""
    refseq_accessions: str = ""
    geneid: str = ""


class UniProtClient:
    """Small wrapper for UniProtKB REST API."""

    def __init__(self, sleep_seconds: float = 0.15, timeout: int = 30):
        self.sleep_seconds = sleep_seconds
        self.timeout = timeout
        self._cache: Dict[str, Optional[UniProtRecord]] = {}

    def _get_json(self, accession: str) -> Optional[dict]:
        url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        try:
            r = requests.get(url, timeout=self.timeout)
            if r.status_code != 200:
                return None
            return r.json()
        except requests.RequestException:
            return None

    @staticmethod
    def _extract_recommended_name(data: dict) -> str:
        desc = data.get("proteinDescription", {}) or {}
        rec = desc.get("recommendedName", {}) or {}
        full = rec.get("fullName", {}) or {}
        value = full.get("value", "")
        if value:
            return value
        sub = desc.get("submissionNames", []) or []
        if sub:
            return sub[0].get("fullName", {}).get("value", "")
        return ""

    @staticmethod
    def _extract_gene_names(data: dict) -> str:
        genes = data.get("genes", []) or []
        names = []
        for g in genes:
            if "geneName" in g and g["geneName"].get("value"):
                names.append(g["geneName"]["value"])
            for syn in g.get("synonyms", []) or []:
                if syn.get("value"):
                    names.append(syn["value"])
        return ",".join(sorted(set(names)))

    @staticmethod
    def _extract_crossrefs(data: dict) -> Dict[str, str]:
        refs = data.get("uniProtKBCrossReferences", []) or []
        out = {"EMBL": [], "RefSeq": [], "GeneID": []}
        for ref in refs:
            db = ref.get("database", "")
            rid = ref.get("id", "")
            if db in out and rid:
                out[db].append(rid)
        return {k: ",".join(sorted(set(v))) for k, v in out.items()}

    def fetch(self, accession: str) -> Optional[UniProtRecord]:
        """Fetch one UniProt record by accession."""
        accession = accession.strip()
        if not accession:
            return None
        if accession in self._cache:
            return self._cache[accession]

        data = self._get_json(accession)
        if not data:
            self._cache[accession] = None
            return None

        xrefs = self._extract_crossrefs(data)
        record = UniProtRecord(
            accession=accession,
            protein_name=self._extract_recommended_name(data),
            organism=(data.get("organism", {}) or {}).get("scientificName", ""),
            gene_names=self._extract_gene_names(data),
            protein_sequence=(data.get("sequence", {}) or {}).get("value", ""),
            embl_accessions=xrefs.get("EMBL", ""),
            refseq_accessions=xrefs.get("RefSeq", ""),
            geneid=xrefs.get("GeneID", ""),
        )
        self._cache[accession] = record
        time.sleep(self.sleep_seconds)
        return record
