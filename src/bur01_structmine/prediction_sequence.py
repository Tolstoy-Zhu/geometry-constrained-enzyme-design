"""FASTA validation and 8HFW catalytic-site mapping for AF2 predictions."""
from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from Bio.Align import PairwiseAligner

from .config import project_root
from .structure_search import load_ca_trace

CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
ALLOWED_AA = CANONICAL_AA | {"X"}
REFERENCE_SITES = {"S111": (111, "S"), "D241": (241, "D"), "H274": (274, "H")}


class SequenceValidationError(ValueError):
    """One or more FASTA validation errors collected before execution."""

    def __init__(self, errors: list[str]) -> None:
        self.errors = errors
        super().__init__("; ".join(errors))


@dataclass(frozen=True)
class FastaRecord:
    name: str
    sequence: str

    @property
    def sha256(self) -> str:
        return hashlib.sha256(self.sequence.encode("ascii")).hexdigest()

    def to_fasta(self) -> str:
        lines = [f">{self.name}"]
        lines.extend(self.sequence[index : index + 80] for index in range(0, len(self.sequence), 80))
        return "\n".join(lines) + "\n"


def safe_sequence_name(value: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip()).strip("_.-")
    return clean[:80] or "query"


def validate_sequence(name: str, sequence: str) -> FastaRecord:
    clean = re.sub(r"\s+", "", sequence).upper()
    errors: list[str] = []
    if not clean:
        errors.append("sequence is empty")
    if len(clean) < 20 or len(clean) > 1000:
        errors.append(f"sequence length must be 20-1000 aa (received {len(clean)})")
    invalid = sorted(set(clean) - ALLOWED_AA)
    if invalid:
        errors.append(f"unsupported residue symbols: {', '.join(invalid)}")
    x_fraction = clean.count("X") / len(clean) if clean else 0.0
    if x_fraction > 0.05:
        errors.append(f"X content must not exceed 5% (received {x_fraction:.1%})")
    if errors:
        raise SequenceValidationError(errors)
    return FastaRecord(safe_sequence_name(name), clean)


def parse_single_fasta(text: str) -> FastaRecord:
    errors: list[str] = []
    lines = [line.strip() for line in text.replace("\ufeff", "").splitlines() if line.strip()]
    headers = [index for index, line in enumerate(lines) if line.startswith(">")]
    if not lines:
        raise SequenceValidationError(["FASTA input is empty"])
    if not headers:
        raise SequenceValidationError(["FASTA header beginning with '>' is required"])
    if len(headers) != 1 or headers[0] != 0:
        errors.append(f"exactly one FASTA record is required (received {len(headers)})")
    header = lines[0][1:].strip() if lines[0].startswith(">") else ""
    if not header:
        errors.append("FASTA header is empty")
    sequence_lines = [line for line in lines[1:] if not line.startswith(">")]
    if any(line.startswith(">") for line in lines[1:]):
        errors.append("multiple FASTA records are not supported in v0.3")
    if errors:
        raise SequenceValidationError(list(dict.fromkeys(errors)))
    return validate_sequence(header.split()[0], "".join(sequence_lines))


def read_single_fasta(path: str | Path) -> FastaRecord:
    path = Path(path)
    if not path.is_file():
        raise SequenceValidationError([f"FASTA file not found: {path}"])
    return parse_single_fasta(path.read_text(encoding="utf-8"))


def _alignment_index_map(reference: str, target: str) -> tuple[dict[int, int], float, float]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(reference, target)[0]
    mapping: dict[int, int] = {}
    matches = 0
    aligned_count = 0
    for (ref_start, ref_end), (target_start, target_end) in zip(*alignment.aligned):
        block_length = min(ref_end - ref_start, target_end - target_start)
        for offset in range(block_length):
            ref_index = int(ref_start + offset)
            target_index = int(target_start + offset)
            mapping[ref_index] = target_index
            matches += int(reference[ref_index] == target[target_index])
            aligned_count += 1
    identity = matches / aligned_count if aligned_count else 0.0
    return mapping, identity, float(alignment.score)


def map_8hfw_catalytic_sites(
    target_sequence: str,
    override: dict[str, int] | None = None,
    reference_pdb: str | Path | None = None,
) -> dict[str, Any]:
    """Map S111/D241/H274 from experimental 8HFW chain A to a new sequence."""
    target = validate_sequence("target", target_sequence).sequence
    reference_path = Path(reference_pdb) if reference_pdb else project_root() / "data" / "query" / "Bur01m.pdb"
    trace = load_ca_trace(reference_path, chain_id="A")
    reference_index = {
        residue.residue_number: index for index, residue in enumerate(trace.residues)
    }
    alignment_map, identity, alignment_score = _alignment_index_map(trace.sequence, target)

    mappings: dict[str, dict[str, Any]] = {}
    valid_positions: dict[str, int] = {}
    missing: list[str] = []
    for label, (reference_number, expected_aa) in REFERENCE_SITES.items():
        source = "manual" if override and label in override else "global_alignment"
        if source == "manual":
            target_position = int(override[label])
        else:
            ref_index = reference_index.get(reference_number)
            target_index = alignment_map.get(ref_index) if ref_index is not None else None
            target_position = target_index + 1 if target_index is not None else 0
        observed = target[target_position - 1] if 1 <= target_position <= len(target) else ""
        if not target_position:
            status = "unmapped"
            missing.append(f"{label}:unmapped")
        elif observed != expected_aa:
            status = "identity_mismatch"
            missing.append(f"{label}:expected_{expected_aa}_observed_{observed or 'missing'}")
        else:
            status = "mapped"
            valid_positions[label] = target_position
        mappings[label] = {
            "reference_residue": label,
            "predicted_position": target_position or None,
            "expected_amino_acid": expected_aa,
            "observed_amino_acid": observed or None,
            "status": status,
            "source": source,
        }

    return {
        "reference_structure": "PDB 8HFW chain A",
        "reference_path": str(reference_path.resolve()),
        "alignment_method": "Biopython PairwiseAligner global",
        "alignment_score": round(alignment_score, 3),
        "aligned_identity": round(identity, 6),
        "sites": mappings,
        "valid_positions": valid_positions,
        "complete": len(valid_positions) == len(REFERENCE_SITES),
        "missing_evidence": missing,
    }
