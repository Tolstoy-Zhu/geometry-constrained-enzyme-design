"""Structure-overlay helpers shared by the Streamlit workbench and tests."""
from __future__ import annotations

import io
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select
from tmtools import tm_align

from .structure_search import (
    DEFAULT_CATALYTIC_QUERY_RESIDUES,
    CaTrace,
    compare_ca_traces,
    load_ca_trace,
    load_ca_traces,
)


@dataclass(frozen=True)
class StructureOverlay:
    query_pdb: str
    target_pdb: str
    query_chain: str
    target_chain: str
    qtmscore: float
    ttmscore: float
    rmsd: float
    aligned_query: str
    alignment_match_line: str
    aligned_target: str
    mapped_catalytic_residues: str
    mapped_target_residue_numbers: tuple[int, ...]
    query_residue_numbers: tuple[int, ...]


class _ChainSelect(Select):
    def __init__(self, chain_id: str) -> None:
        self.chain_id = chain_id

    def accept_chain(self, chain) -> bool:  # noqa: ANN001 - Biopython callback protocol.
        return chain.id == self.chain_id


def _parser(path: Path):
    return MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)


def _chain_as_pdb(
    path: Path,
    chain_id: str,
    rotation: np.ndarray | None = None,
    translation: np.ndarray | None = None,
) -> str:
    structure = _parser(path).get_structure(path.stem, str(path))
    selected = False
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            selected = True
            if rotation is not None and translation is not None:
                for atom in chain.get_atoms():
                    atom.set_coord(rotation @ np.asarray(atom.coord, dtype=float) + translation)
    if not selected:
        raise ValueError(f"Chain {chain_id!r} was not found in {path}")
    output = io.StringIO()
    writer = PDBIO()
    writer.set_structure(structure)
    writer.save(output, select=_ChainSelect(chain_id))
    return output.getvalue()


def _choose_target(path: Path, chain_id: str | None) -> CaTrace:
    if chain_id:
        return load_ca_trace(path, chain_id=chain_id)
    traces = load_ca_traces(path)
    return max(traces, key=lambda trace: len(trace.residues))


def align_structure_overlay(
    query_path: str | Path,
    target_path: str | Path,
    query_chain: str = "A",
    target_chain: str | None = None,
    catalytic_query_residues: tuple[tuple[int, str], ...] = DEFAULT_CATALYTIC_QUERY_RESIDUES,
) -> StructureOverlay:
    """Superpose the selected query chain onto a selected target chain with TM-align."""
    query_path = Path(query_path).resolve()
    target_path = Path(target_path).resolve()
    query = load_ca_trace(query_path, chain_id=query_chain)
    target = _choose_target(target_path, target_chain)
    alignment = tm_align(query.coords, target.coords, query.sequence, target.sequence)
    metrics = compare_ca_traces(query, target, catalytic_query_residues)
    query_pdb = _chain_as_pdb(
        query_path,
        query.chain_id,
        rotation=np.asarray(alignment.u, dtype=float),
        translation=np.asarray(alignment.t, dtype=float),
    )
    target_pdb = _chain_as_pdb(target_path, target.chain_id)
    mapped_numbers = tuple(
        int(match.group(1))
        for label in str(metrics["mapped_catalytic_residues"]).split(";")
        if (match := re.search(r"(-?\d+)", label))
    )
    return StructureOverlay(
        query_pdb=query_pdb,
        target_pdb=target_pdb,
        query_chain=query.chain_id,
        target_chain=target.chain_id,
        qtmscore=float(alignment.tm_norm_chain1),
        ttmscore=float(alignment.tm_norm_chain2),
        rmsd=float(alignment.rmsd),
        aligned_query=alignment.seqxA,
        alignment_match_line=alignment.seqM,
        aligned_target=alignment.seqyA,
        mapped_catalytic_residues=str(metrics["mapped_catalytic_residues"]),
        mapped_target_residue_numbers=mapped_numbers,
        query_residue_numbers=tuple(residue.residue_number for residue in query.residues),
    )
