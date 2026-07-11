"""Cached AlphaFold DB prediction API client."""
from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Optional

from .config import project_root
from .http_utils import CachedSession


@dataclass
class AlphaFoldRecord:
    accession: str = ""
    model_entity_id: str = ""
    latest_version: str = ""
    model_created_date: str = ""
    sequence_version_date: str = ""
    global_plddt: str = ""
    is_reviewed: str = ""
    organism: str = ""
    protein_description: str = ""
    pdb_url: str = ""
    cif_url: str = ""
    bcif_url: str = ""
    pae_url: str = ""
    plddt_url: str = ""


class AlphaFoldClient:
    """Retrieve AlphaFold model provenance and confidence metadata."""

    def __init__(
        self,
        sleep_seconds: float = 0.15,
        timeout: int = 30,
        cache_dir: str | Path | None = None,
    ) -> None:
        self.sleep_seconds = sleep_seconds
        self._cache: Dict[str, Optional[AlphaFoldRecord]] = {}
        self.errors: Dict[str, str] = {}
        self.http = CachedSession(
            cache_dir or project_root() / ".cache" / "annotations" / "alphafold",
            timeout=timeout,
        )

    @staticmethod
    def _choose_best(data: list) -> dict:
        latest = [item for item in data if item.get("latestVersion")]
        return latest[0] if latest else data[0]

    @staticmethod
    def _to_record(accession: str, item: dict) -> AlphaFoldRecord:
        return AlphaFoldRecord(
            accession=accession,
            model_entity_id=str(item.get("modelEntityId", "")),
            latest_version=str(item.get("latestVersion", "")),
            model_created_date=str(item.get("modelCreatedDate", "")),
            sequence_version_date=str(item.get("sequenceVersionDate", "")),
            global_plddt=str(item.get("globalMetricValue", "")),
            is_reviewed=str(item.get("isReviewed", item.get("isUniProtReviewed", ""))),
            organism=str(item.get("organismScientificName", "")),
            protein_description=str(item.get("uniprotDescription", "")),
            pdb_url=str(item.get("pdbUrl", "")),
            cif_url=str(item.get("cifUrl", "")),
            bcif_url=str(item.get("bcifUrl", "")),
            pae_url=str(item.get("paeDocUrl", "")),
            plddt_url=str(item.get("plddtDocUrl", "")),
        )

    def fetch_many(self, accessions: Iterable[str]) -> Dict[str, Optional[AlphaFoldRecord]]:
        normalized = list(dict.fromkeys(str(accession).strip().upper() for accession in accessions if str(accession).strip()))
        for accession in normalized:
            if accession in self._cache:
                continue
            data = self.http.get_json(f"https://alphafold.ebi.ac.uk/api/prediction/{accession}")
            if isinstance(data, list) and data:
                self._cache[accession] = self._to_record(accession, self._choose_best(data))
            else:
                self._cache[accession] = None
                self.errors[accession] = self.http.last_error or "AlphaFold prediction was not returned"
            if self.sleep_seconds:
                time.sleep(self.sleep_seconds)
        return {accession: self._cache.get(accession) for accession in normalized}

    def fetch(self, accession: str) -> Optional[AlphaFoldRecord]:
        accession = accession.strip().upper()
        if not accession:
            return None
        return self.fetch_many([accession]).get(accession)
