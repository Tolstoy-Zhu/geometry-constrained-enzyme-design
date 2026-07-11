"""Chain-aware small-set structure search backed by TM-align."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser

from .foldseek_runner import FOLDSEEK_FIELDS
from .io_utils import ensure_dir
from .pdb_utils import AA3_TO_AA1

STRUCTURE_SUFFIXES = {".pdb", ".ent", ".cif", ".mmcif"}
DEFAULT_CATALYTIC_QUERY_RESIDUES = ((111, "S"), (241, "D"), (274, "H"))
LOCAL_TM_EXTRA_FIELDS = [
    "search_engine",
    "metric_source",
    "alignment_version",
    "query_chain",
    "target_chain",
    "aligned_query",
    "alignment_match_line",
    "aligned_target",
    "query_structure_path",
    "target_structure_path",
    "query_catalytic_residues",
    "mapped_catalytic_residues",
    "mapped_catalytic_identity_count",
    "catalytic_ser_his_distance",
    "catalytic_asp_his_distance",
    "catalytic_geometry_contacts",
]
LOCAL_TM_FIELDS = list(FOLDSEEK_FIELDS) + LOCAL_TM_EXTRA_FIELDS


@dataclass
class TraceResidue:
    chain_id: str
    residue_number: int
    insertion_code: str
    amino_acid: str
    ca: np.ndarray
    atoms: dict[str, np.ndarray]

    @property
    def label(self) -> str:
        suffix = self.insertion_code.strip()
        return f"{self.amino_acid}{self.residue_number}{suffix}"


@dataclass
class CaTrace:
    """One structure chain represented by its ordered C-alpha trace."""

    name: str
    path: Path
    chain_id: str
    sequence: str
    coords: np.ndarray
    residues: list[TraceResidue]


def iter_structure_files(structure_dir: str | Path) -> Iterable[Path]:
    """Yield supported structure files from a directory."""
    structure_dir = Path(structure_dir)
    if not structure_dir.exists():
        raise FileNotFoundError(f"Structure directory not found: {structure_dir}")
    for path in sorted(structure_dir.rglob("*")):
        if path.is_file() and path.suffix.lower() in STRUCTURE_SUFFIXES:
            yield path


def _parse_structure(path: Path):
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    return parser.get_structure(path.stem, str(path))


def load_ca_traces(
    path: str | Path,
    chain_id: str | None = None,
    name: str | None = None,
) -> list[CaTrace]:
    """Load each protein chain independently from a PDB or mmCIF file."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)

    structure = _parse_structure(path)
    model = next(structure.get_models())
    traces: list[CaTrace] = []
    chains = [chain for chain in model if chain_id is None or chain.id == chain_id]
    if chain_id is not None and not chains:
        raise ValueError(f"Chain {chain_id!r} was not found in {path}")

    for chain in chains:
        residues: list[TraceResidue] = []
        for residue in chain.get_residues():
            resname = residue.get_resname().upper()
            if resname not in AA3_TO_AA1 or "CA" not in residue:
                continue
            hetero, number, insertion = residue.id
            if hetero.strip():
                continue
            atoms = {
                atom.name: np.asarray(atom.coord, dtype=float)
                for atom in residue.get_unpacked_list()
            }
            residues.append(
                TraceResidue(
                    chain_id=chain.id,
                    residue_number=int(number),
                    insertion_code=str(insertion or "").strip(),
                    amino_acid=AA3_TO_AA1[resname],
                    ca=np.asarray(residue["CA"].coord, dtype=float),
                    atoms=atoms,
                )
            )
        if not residues:
            continue
        base_name = name or path.stem
        trace_name = base_name if len(chains) == 1 else f"{base_name}_chain_{chain.id}"
        traces.append(
            CaTrace(
                name=trace_name,
                path=path,
                chain_id=chain.id,
                sequence="".join(residue.amino_acid for residue in residues),
                coords=np.vstack([residue.ca for residue in residues]),
                residues=residues,
            )
        )

    if not traces:
        requested = f" chain {chain_id}" if chain_id else ""
        raise ValueError(f"No C-alpha protein trace found in{requested} {path}")
    return traces


def load_ca_trace(path: str | Path, name: str | None = None, chain_id: str | None = None) -> CaTrace:
    """Load one chain, rejecting implicit concatenation of multi-chain files."""
    traces = load_ca_traces(path, chain_id=chain_id, name=name)
    if chain_id is None and len(traces) > 1:
        available = ", ".join(trace.chain_id for trace in traces)
        raise ValueError(f"Multiple protein chains found in {path}: {available}. Select a chain explicitly.")
    return traces[0]


def _alignment_map(aligned_query: str, aligned_target: str) -> tuple[dict[int, int], int, int]:
    query_index = -1
    target_index = -1
    mapping: dict[int, int] = {}
    aligned = 0
    matches = 0
    for query_aa, target_aa in zip(aligned_query, aligned_target):
        if query_aa != "-":
            query_index += 1
        if target_aa != "-":
            target_index += 1
        if query_aa != "-" and target_aa != "-":
            mapping[query_index] = target_index
            aligned += 1
            matches += int(query_aa == target_aa)
    return mapping, aligned, matches


def _minimum_atom_distance(
    first: TraceResidue | None,
    first_atoms: tuple[str, ...],
    second: TraceResidue | None,
    second_atoms: tuple[str, ...],
) -> float | None:
    if first is None or second is None:
        return None
    distances = [
        float(np.linalg.norm(first.atoms[a] - second.atoms[b]))
        for a in first_atoms
        for b in second_atoms
        if a in first.atoms and b in second.atoms
    ]
    return min(distances) if distances else None


def _map_catalytic_residues(
    query: CaTrace,
    target: CaTrace,
    mapping: dict[int, int],
    catalytic_query_residues: tuple[tuple[int, str], ...],
) -> dict[str, object]:
    query_index_by_number = {
        residue.residue_number: index for index, residue in enumerate(query.residues)
    }
    mapped: list[TraceResidue | None] = []
    identity_count = 0
    labels: list[str] = []
    for residue_number, expected_aa in catalytic_query_residues:
        query_index = query_index_by_number.get(residue_number)
        target_index = mapping.get(query_index) if query_index is not None else None
        residue = target.residues[target_index] if target_index is not None else None
        mapped.append(residue)
        labels.append(residue.label if residue else "missing")
        identity_count += int(residue is not None and residue.amino_acid == expected_aa)

    mapped_by_identity = {
        expected_aa: residue
        for (_, expected_aa), residue in zip(catalytic_query_residues, mapped)
    }
    serine = mapped_by_identity.get("S")
    aspartate = mapped_by_identity.get("D")
    histidine = mapped_by_identity.get("H")
    ser_his = _minimum_atom_distance(serine, ("OG",), histidine, ("NE2", "ND1"))
    asp_his = _minimum_atom_distance(aspartate, ("OD1", "OD2"), histidine, ("ND1", "NE2"))
    contacts = sum(
        int(distance is not None and 2.0 <= distance <= 5.0)
        for distance in (ser_his, asp_his)
    )
    return {
        "query_catalytic_residues": ";".join(
            f"{amino_acid}{residue_number}"
            for residue_number, amino_acid in catalytic_query_residues
        ),
        "mapped_catalytic_residues": ";".join(labels),
        "mapped_catalytic_identity_count": identity_count,
        "catalytic_ser_his_distance": ser_his,
        "catalytic_asp_his_distance": asp_his,
        "catalytic_geometry_contacts": contacts,
    }


def compare_ca_traces(
    query: CaTrace,
    target: CaTrace,
    catalytic_query_residues: tuple[tuple[int, str], ...] = DEFAULT_CATALYTIC_QUERY_RESIDUES,
) -> dict[str, object]:
    """Align two chains with TM-align and return explicit, non-fabricated metrics."""
    try:
        from tmtools import tm_align
    except ImportError as exc:  # pragma: no cover - dependency check is covered by doctor.
        raise RuntimeError("tmtools 0.3.0 is required for the local TM-align preview") from exc

    result = tm_align(query.coords, target.coords, query.sequence, target.sequence)
    mapping, alnlen, matches = _alignment_map(result.seqxA, result.seqyA)
    qlen = len(query.residues)
    tlen = len(target.residues)
    qtmscore = float(result.tm_norm_chain1)
    ttmscore = float(result.tm_norm_chain2)
    row: dict[str, object] = {
        "query": query.name,
        "target": target.name,
        "evalue": np.nan,
        "bits": np.nan,
        "prob": np.nan,
        "alntmscore": min(qtmscore, ttmscore),
        "qtmscore": qtmscore,
        "ttmscore": ttmscore,
        "rmsd": float(result.rmsd),
        "fident": matches / alnlen if alnlen else 0.0,
        "alnlen": alnlen,
        "qlen": qlen,
        "tlen": tlen,
        "qcov": alnlen / qlen if qlen else 0.0,
        "tcov": alnlen / tlen if tlen else 0.0,
        "search_engine": "tmtools_tmalign",
        "metric_source": "TM-align",
        "alignment_version": "tmtools-0.3.0/TM-align-20210224",
        "query_chain": query.chain_id,
        "target_chain": target.chain_id,
        "aligned_query": result.seqxA,
        "alignment_match_line": result.seqM,
        "aligned_target": result.seqyA,
        "query_structure_path": str(query.path.resolve()),
        "target_structure_path": str(target.path.resolve()),
    }
    row.update(_map_catalytic_residues(query, target, mapping, catalytic_query_residues))
    return row


def run_local_tm_search(
    query_pdb: str | Path,
    structure_dir: str | Path,
    out_tsv: str | Path,
    query_chain: str | None = None,
    max_targets: int | None = None,
    catalytic_query_residues: tuple[tuple[int, str], ...] = DEFAULT_CATALYTIC_QUERY_RESIDUES,
) -> Path:
    """Compare one query chain against each chain in a small structure directory."""
    out_tsv = Path(out_tsv)
    ensure_dir(out_tsv.parent)
    query = load_ca_trace(query_pdb, name=Path(query_pdb).stem, chain_id=query_chain)

    rows: list[dict[str, object]] = []
    for path in iter_structure_files(structure_dir):
        try:
            for target in load_ca_traces(path):
                rows.append(compare_ca_traces(query, target, catalytic_query_residues))
                if max_targets is not None and len(rows) >= max_targets:
                    break
        except Exception as exc:  # noqa: BLE001 - preserve batch progress and surface the skipped file.
            print(f"[WARN] Skipping {path}: {exc}")
        if max_targets is not None and len(rows) >= max_targets:
            break

    if not rows:
        raise ValueError(f"No comparable structure chains were found in {structure_dir}")

    df = pd.DataFrame(rows).reindex(columns=LOCAL_TM_FIELDS)
    df = df.sort_values(["alntmscore", "qcov", "tcov"], ascending=False)
    df.to_csv(out_tsv, sep="\t", index=False, header=False)
    return out_tsv


def run_local_ca_search(
    query_pdb: str | Path,
    structure_dir: str | Path,
    out_tsv: str | Path,
    max_targets: int | None = None,
    query_chain: str | None = None,
    catalytic_query_residues: tuple[tuple[int, str], ...] = DEFAULT_CATALYTIC_QUERY_RESIDUES,
) -> Path:
    """Compatibility alias for the v0.2 TM-align implementation."""
    return run_local_tm_search(
        query_pdb=query_pdb,
        structure_dir=structure_dir,
        out_tsv=out_tsv,
        query_chain=query_chain,
        max_targets=max_targets,
        catalytic_query_residues=catalytic_query_residues,
    )
