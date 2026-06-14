"""PDB parsing utilities."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List

from Bio.PDB import PDBParser, PPBuilder

AA3_TO_AA1: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEC": "U", "PYL": "O",
}


def extract_sequence_from_pdb(pdb_file: str | Path, chain_id: str | None = None) -> Dict[str, str]:
    """Extract amino-acid sequences from a PDB file.

    Returns a dictionary: {chain_id: sequence}.
    """
    pdb_file = Path(pdb_file)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("query", str(pdb_file))
    model = next(structure.get_models())

    pp_builder = PPBuilder()
    results: Dict[str, str] = {}

    for chain in model:
        if chain_id and chain.id != chain_id:
            continue

        peptides = pp_builder.build_peptides(chain)
        seq = "".join(str(pp.get_sequence()) for pp in peptides)

        # Fallback: build sequence from residues when PPBuilder fails.
        if not seq:
            chars: List[str] = []
            for residue in chain.get_residues():
                resname = residue.get_resname().upper()
                if resname in AA3_TO_AA1:
                    chars.append(AA3_TO_AA1[resname])
            seq = "".join(chars)

        if seq:
            results[chain.id] = seq

    return results
