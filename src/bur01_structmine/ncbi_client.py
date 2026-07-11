"""Cached NCBI Protein E-utilities client."""
from __future__ import annotations

import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional

from .config import project_root
from .http_utils import CachedSession


@dataclass
class NcbiProteinRecord:
    accession: str = ""
    protein_name: str = ""
    organism: str = ""
    protein_sequence: str = ""
    fasta_header: str = ""


class NcbiProteinClient:
    """Retrieve NCBI Protein FASTA records with rate limiting and cache."""

    def __init__(
        self,
        sleep_seconds: float = 0.34,
        timeout: int = 30,
        cache_dir: str | Path | None = None,
    ) -> None:
        self.sleep_seconds = sleep_seconds
        self._cache: Dict[str, Optional[NcbiProteinRecord]] = {}
        self.errors: Dict[str, str] = {}
        self.http = CachedSession(
            cache_dir or project_root() / ".cache" / "annotations" / "ncbi",
            timeout=timeout,
        )

    @staticmethod
    def _parse_fasta(accession: str, fasta_text: str) -> NcbiProteinRecord:
        lines = [line.strip() for line in fasta_text.splitlines() if line.strip()]
        header = lines[0][1:] if lines else ""
        sequence = "".join(lines[1:])
        parts = header.split(maxsplit=1)
        name = parts[1] if len(parts) == 2 else ""
        if name.endswith("]"):
            name = re.sub(r"\s+\[[^\[\]]+\]\s*$", "", name)
        match = re.search(r"\[([^\[\]]+)\]\s*$", header)
        return NcbiProteinRecord(
            accession=accession,
            protein_name=name,
            organism=match.group(1) if match else "",
            protein_sequence=sequence,
            fasta_header=header,
        )

    def fetch_many(self, accessions: Iterable[str]) -> Dict[str, Optional[NcbiProteinRecord]]:
        normalized = list(dict.fromkeys(str(accession).strip() for accession in accessions if str(accession).strip()))
        endpoint = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        for accession in normalized:
            if accession in self._cache:
                continue
            text = self.http.get_text(
                endpoint,
                params={"db": "protein", "id": accession, "rettype": "fasta", "retmode": "text"},
            )
            if text and text.startswith(">"):
                self._cache[accession] = self._parse_fasta(accession, text)
            else:
                self._cache[accession] = None
                self.errors[accession] = self.http.last_error or "NCBI Protein FASTA was not returned"
            if self.sleep_seconds:
                time.sleep(self.sleep_seconds)
        return {accession: self._cache.get(accession) for accession in normalized}

    def fetch(self, accession: str) -> Optional[NcbiProteinRecord]:
        accession = accession.strip()
        if not accession:
            return None
        return self.fetch_many([accession]).get(accession)
