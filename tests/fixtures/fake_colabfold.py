"""Small deterministic colabfold_batch stand-in used by non-GPU CI tests."""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path


AA1_TO_AA3 = {
    "A": "ALA", "C": "CYS", "D": "ASP", "E": "GLU", "F": "PHE",
    "G": "GLY", "H": "HIS", "I": "ILE", "K": "LYS", "L": "LEU",
    "M": "MET", "N": "ASN", "P": "PRO", "Q": "GLN", "R": "ARG",
    "S": "SER", "T": "THR", "V": "VAL", "W": "TRP", "Y": "TYR",
    "X": "GLY",
}


def option(name: str, default: str) -> str:
    try:
        return sys.argv[sys.argv.index(name) + 1]
    except (ValueError, IndexError):
        return default


def read_sequence(path: Path) -> tuple[str, str]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    return lines[0].lstrip(">") or "query", "".join(lines[1:]).upper()


def write_pdb(path: Path, sequence: str, confidence: float) -> None:
    lines = []
    for index, amino_acid in enumerate(sequence, start=1):
        x = float(index) * 1.5
        lines.append(
            f"ATOM  {index:5d}  CA  {AA1_TO_AA3[amino_acid]:>3s} A{index:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{confidence:6.2f}           C"
        )
    lines.extend(["TER", "END"])
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def main() -> int:
    input_path = Path(sys.argv[1])
    output_dir = Path(sys.argv[2])
    output_dir.mkdir(parents=True, exist_ok=True)
    name, sequence = read_sequence(input_path)
    safe_name = "".join(character if character.isalnum() else "_" for character in name)
    sleep_seconds = float(os.environ.get("FAKE_COLABFOLD_SLEEP", "0"))
    if sleep_seconds:
        time.sleep(sleep_seconds)
    if "--msa-only" in sys.argv:
        fail_attempts = int(os.environ.get("FAKE_COLABFOLD_FAIL_MSA_ATTEMPTS", "0"))
        counter_path = os.environ.get("FAKE_COLABFOLD_COUNTER")
        if fail_attempts and counter_path:
            counter = Path(counter_path)
            count = int(counter.read_text(encoding="ascii")) if counter.exists() else 0
            counter.write_text(str(count + 1), encoding="ascii")
            if count < fail_attempts:
                print("simulated public MSA failure", file=sys.stderr)
                return 7
        (output_dir / f"{safe_name}.a3m").write_text(
            f">101\n{sequence}\n",
            encoding="ascii",
        )
        return 0

    num_models = int(option("--num-models", "5"))
    num_seeds = int(option("--num-seeds", "1"))
    num_relax = int(option("--num-relax", "0"))
    total = num_models * num_seeds
    pae = [[round(abs(row - column) / 4.0, 2) for column in range(len(sequence))] for row in range(len(sequence))]
    for rank in range(1, total + 1):
        model = (rank - 1) % num_models + 1
        seed = (rank - 1) // num_models
        tag = f"rank_{rank:03d}_alphafold2_ptm_model_{model}_seed_{seed:03d}"
        confidence = 88.0 - rank
        state = "relaxed" if rank <= num_relax else "unrelaxed"
        write_pdb(output_dir / f"{safe_name}_{state}_{tag}.pdb", sequence, confidence)
        (output_dir / f"{safe_name}_scores_{tag}.json").write_text(
            json.dumps(
                {
                    "plddt": [confidence] * len(sequence),
                    "ptm": round(0.86 - rank / 100.0, 3),
                    "max_pae": max(max(row) for row in pae),
                    "pae": pae,
                }
            ),
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
