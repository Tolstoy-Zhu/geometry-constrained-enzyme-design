"""Coordinate-only ligand contacts and insertion-code-safe site selections."""
from __future__ import annotations

import hashlib
import io
import json
from dataclasses import dataclass

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFIO, MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa

from .models import SiteSelection

WATER = {"HOH", "WAT", "DOD"}
IONS = {"NA", "K", "CL", "CA", "MG", "ZN", "FE", "MN", "CU", "CO", "NI", "CD", "IOD", "BR", "F", "LI", "CS", "SR"}


@dataclass
class StructureDocument:
    structure: object
    sha256: str

    @property
    def models(self) -> list[str]:
        return [str(model.id) for model in self.structure]

    def model(self, model_id: str):
        return next(model for model in self.structure if str(model.id) == str(model_id))

    def chains(self, model_id: str) -> list[str]:
        return [chain.id for chain in self.model(model_id)]

    def model_cif(self, model_id: str) -> str:
        writer = MMCIFIO()
        writer.set_structure(self.model(model_id))
        output = io.StringIO()
        writer.save(output)
        return output.getvalue()


def parse_structure_bytes(content: bytes, filename: str) -> StructureDocument:
    if len(content) > 30 * 1024 * 1024:
        raise ValueError("Structure exceeds 30 MiB")
    suffix = filename.rsplit(".", 1)[-1].lower()
    if suffix not in {"pdb", "cif", "mmcif"}:
        raise ValueError("Expected PDB/CIF/mmCIF")
    parser = PDBParser(QUIET=True) if suffix == "pdb" else MMCIFParser(QUIET=True)
    structure = parser.get_structure("uploaded", io.StringIO(content.decode("utf-8-sig")))
    if not list(structure.get_atoms()):
        raise ValueError("No coordinate atoms found")
    if not all(np.isfinite(atom.coord).all() for atom in structure.get_atoms()):
        raise ValueError("Non-finite atom coordinates")
    return StructureDocument(structure, hashlib.sha256(content).hexdigest())


def residue_key(residue) -> str:
    return json.dumps([residue.parent.id, residue.id[0], int(residue.id[1]), residue.id[2].strip()], separators=(",", ":"))


def residue_info(residue) -> dict:
    return {"residue_id": residue_key(residue), "chain_id": residue.parent.id,
            "residue_number": int(residue.id[1]), "insertion_code": residue.id[2].strip(), "residue_name": residue.resname,
            "label": f"{residue.parent.id}:{residue.resname}{residue.id[1]}{residue.id[2].strip()}"}


def nonprotein_components(doc: StructureDocument, model_id: str) -> pd.DataFrame:
    rows = []
    for residue in doc.model(model_id).get_residues():
        if is_aa(residue, standard=False):
            continue
        kind = "water" if residue.resname in WATER else "ion" if residue.resname in IONS else "ligand_or_cofactor"
        rows.append({**residue_info(residue), "kind": kind, "default_selected": kind == "ligand_or_cofactor",
                     "heavy_atoms": len(_heavy_atoms(residue))})
    return pd.DataFrame(rows, columns=["residue_id", "chain_id", "residue_number", "insertion_code", "residue_name", "label", "kind", "default_selected", "heavy_atoms"])


def _heavy_atoms(residue) -> list:
    # Biopython chooses the highest-occupancy altloc; never mix conformers.
    return [atom for atom in residue.get_atoms() if str(atom.element).upper() not in {"H", "D"}]


def ligand_contacts(doc: StructureDocument, model_id: str, chains: list[str], ligand_ids: list[str], cutoff: float = 4.0) -> pd.DataFrame:
    if not 0 < cutoff <= 20:
        raise ValueError("Contact cutoff must be in (0, 20] Angstrom")
    model = doc.model(model_id)
    residues = {residue_key(residue): residue for residue in model.get_residues()}
    if not set(ligand_ids).issubset(residues):
        raise ValueError("Selected ligand is absent from the selected model")
    rows = []
    for ligand_id in ligand_ids:
        ligand = residues[ligand_id]
        ligand_atoms = _heavy_atoms(ligand)
        if not ligand_atoms:
            continue
        coords = np.array([atom.coord for atom in ligand_atoms])
        for residue in model.get_residues():
            if residue.parent.id not in chains or not is_aa(residue, standard=False):
                continue
            atoms = _heavy_atoms(residue)
            if not atoms:
                continue
            distances = np.linalg.norm(np.array([atom.coord for atom in atoms])[:, None, :] - coords[None, :, :], axis=2)
            i, j = np.unravel_index(np.argmin(distances), distances.shape)
            if distances[i, j] <= cutoff:
                rows.append({**residue_info(residue), "model_id": str(model_id), "ligand_id": ligand_id,
                             "protein_atom": atoms[i].name, "ligand_atom": ligand_atoms[j].name,
                             "minimum_distance_angstrom": float(distances[i, j]), "heavy_atom_pairs": int((distances <= cutoff).sum()),
                             "cutoff_angstrom": cutoff, "contact_type": "distance_proximity_not_confirmed_hbond"})
    columns = ["residue_id", "chain_id", "residue_number", "insertion_code", "residue_name", "label", "model_id", "ligand_id",
               "protein_atom", "ligand_atom", "minimum_distance_angstrom", "heavy_atom_pairs", "cutoff_angstrom", "contact_type"]
    return pd.DataFrame(rows, columns=columns).sort_values("minimum_distance_angstrom", ignore_index=True)


def protein_sites(doc: StructureDocument, model_id: str, chains: list[str]) -> pd.DataFrame:
    return pd.DataFrame([residue_info(residue) for residue in doc.model(model_id).get_residues()
                         if residue.parent.id in chains and is_aa(residue, standard=False)])


def site_selections(doc: StructureDocument, model_id: str, marks: dict[str, str], ligand_ids: list[str], cutoff: float,
                    substrate_name: str = "", smiles: str = "", conditions: str = "") -> list[SiteSelection]:
    residues = {residue_key(residue): residue for residue in doc.model(model_id).get_residues() if is_aa(residue, standard=False)}
    contacts = ligand_contacts(doc, model_id, doc.chains(model_id), ligand_ids, cutoff)
    distances = contacts.groupby("residue_id").minimum_distance_angstrom.min().to_dict() if not contacts.empty else {}
    output = []
    for key, mark in marks.items():
        if key not in residues:
            raise ValueError(f"Unknown protein site: {key}")
        info = residue_info(residues[key])
        output.append(SiteSelection(structure_sha256=doc.sha256, model_id=str(model_id),
                                    **{col: info[col] for col in ("chain_id", "residue_number", "insertion_code", "residue_name")},
                                    mark=mark, ligand_ids=ligand_ids, cutoff_angstrom=cutoff, minimum_distance_angstrom=distances.get(key),
                                    substrate_name=substrate_name, smiles=smiles, conditions=conditions))
    return output


def viewer_selection(residue_id: str) -> dict:
    chain, hetero, number, insertion = json.loads(residue_id)
    return {"chain": chain, "resi": number, "icode": insertion, "hetflag": bool(hetero.strip())}
