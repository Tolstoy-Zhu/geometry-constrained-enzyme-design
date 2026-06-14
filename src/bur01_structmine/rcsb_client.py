"""RCSB PDB REST API client.

This module converts PDB hits into protein annotation, source organism,
protein sequence, and UniProt accession when available.
"""
from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Dict, List, Optional

import requests


@dataclass
class RcsbEntity:
    pdb_id: str = ""
    entity_id: str = ""
    auth_asym_ids: str = ""
    asym_ids: str = ""
    protein_name: str = ""
    organism: str = ""
    protein_sequence: str = ""
    uniprot_accessions: str = ""


class RcsbClient:
    """Small wrapper for RCSB REST endpoints."""

    def __init__(self, sleep_seconds: float = 0.15, timeout: int = 30):
        self.sleep_seconds = sleep_seconds
        self.timeout = timeout
        self._entry_cache: Dict[str, List[RcsbEntity]] = {}

    def _get_json(self, url: str) -> Optional[dict]:
        try:
            r = requests.get(url, timeout=self.timeout)
            if r.status_code != 200:
                return None
            return r.json()
        except requests.RequestException:
            return None

    def get_polymer_entities(self, pdb_id: str) -> List[RcsbEntity]:
        """Return polymer entities for a PDB entry."""
        pdb_id = pdb_id.upper()
        if pdb_id in self._entry_cache:
            return self._entry_cache[pdb_id]

        entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        entry = self._get_json(entry_url)
        if not entry:
            self._entry_cache[pdb_id] = []
            return []

        entity_ids = (
            entry.get("rcsb_entry_container_identifiers", {})
            .get("polymer_entity_ids", [])
        )

        entities: List[RcsbEntity] = []
        for entity_id in entity_ids:
            url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
            ent = self._get_json(url)
            if not ent:
                continue

            ids = ent.get("rcsb_polymer_entity_container_identifiers", {})
            refs = ids.get("reference_sequence_identifiers", []) or []
            uniprot_accs = []
            for ref in refs:
                db_name = str(ref.get("database_name", "")).lower()
                if db_name == "uniprot":
                    acc = ref.get("database_accession")
                    if acc:
                        uniprot_accs.append(acc)

            source_orgs = ent.get("rcsb_entity_source_organism", []) or []
            organism = "; ".join(
                sorted({x.get("scientific_name", "") for x in source_orgs if x.get("scientific_name")})
            )

            entities.append(
                RcsbEntity(
                    pdb_id=pdb_id,
                    entity_id=str(entity_id),
                    auth_asym_ids=",".join(ids.get("auth_asym_ids", []) or []),
                    asym_ids=",".join(ids.get("asym_ids", []) or []),
                    protein_name=ent.get("rcsb_polymer_entity", {}).get("pdbx_description", ""),
                    organism=organism,
                    protein_sequence=ent.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", ""),
                    uniprot_accessions=",".join(sorted(set(uniprot_accs))),
                )
            )
            time.sleep(self.sleep_seconds)

        self._entry_cache[pdb_id] = entities
        return entities

    @staticmethod
    def choose_entity_for_chain(entities: List[RcsbEntity], chain: Optional[str]) -> Optional[RcsbEntity]:
        """Choose an entity that contains the requested chain."""
        if not entities:
            return None
        if chain:
            for ent in entities:
                chains = set((ent.auth_asym_ids + "," + ent.asym_ids).split(","))
                if chain in chains:
                    return ent
        return entities[0]
