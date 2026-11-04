"""Pinned local inverse-folding models. Network access is restricted to explicit setup."""
from __future__ import annotations

import hashlib
import os
import tempfile
from pathlib import Path

import numpy as np
import requests

from .config import project_root

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
REVISION = "26ec57ac976ade5379920dbd43c7f97a91cf82de"
SOURCE = f"https://raw.githubusercontent.com/dauparas/LigandMPNN/{REVISION}/"
WEIGHTS = "https://files.ipd.uw.edu/pub/ligandmpnn/"
CHECKPOINTS = {"protein_mpnn": "proteinmpnn_v_48_020.pt", "ligand_mpnn": "ligandmpnn_v_32_010_25.pt"}
FILES = {
    "model_utils.py": (SOURCE, "5ad8f6992f83cdd1c97103e19c80833afff2bff201662a23bcb06eae669f0f95"),
    "LICENSE": (SOURCE, "8af7ba31868e08e8949e20ce2ef3c0b1d0ef78c9358c91563a60a9fe3b36c9ea"),
    "proteinmpnn_v_48_020.pt": (WEIGHTS, "c9cb4a671d79604111231f8dbfc7c590e06f1197453b7a6854ac6661a642f5bd"),
    "ligandmpnn_v_32_010_25.pt": (WEIGHTS, "161cd264061fda9680cbb940255522ae42f2966c552d045d87913d9452a80970"),
}


def model_dir() -> Path:
    return Path(os.environ.get("AQUA_DESIGN_MODELS") or project_root() / ".tools" / "aqua-design")


def verified_bytes(name: str) -> bytes:
    path = model_dir() / name
    if not path.is_file():
        raise ValueError(f"Missing {name}. Run aquaesternet-structmine design-setup first.")
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != FILES[name][1]:
        raise ValueError(f"Checksum mismatch for {name}. Run design-setup to repair it.")
    return data


def install_models() -> Path:
    directory = model_dir()
    directory.mkdir(parents=True, exist_ok=True)
    for name, (base, expected) in FILES.items():
        try:
            verified_bytes(name)
            continue
        except ValueError:
            pass
        with requests.get(base + name, stream=True, timeout=(12, 60)) as response:
            response.raise_for_status()
            digest = hashlib.sha256()
            with tempfile.NamedTemporaryFile(dir=directory, delete=False) as target:
                tmp = Path(target.name)
            try:
                with tmp.open("wb") as target:
                    for block in response.iter_content(256 * 1024):
                        digest.update(block)
                        target.write(block)
                if digest.hexdigest() != expected:
                    raise ValueError(f"Unexpected SHA-256 for {name} from the upstream source.")
                os.replace(tmp, directory / name)
            finally:
                tmp.unlink(missing_ok=True)
    return directory


def conditional_scores(trace, positions, ligand_xyz, ligand_types, *, seed, device):
    """Score each requested position last, conditioned on every other parent residue."""
    import io
    try:
        import torch
    except ImportError as exc:
        raise ValueError("Install PyTorch in this Python environment before running a design.") from exc
    if device not in {"cpu", "cuda"} or (device == "cuda" and not torch.cuda.is_available()):
        raise ValueError(f"Requested PyTorch device is unavailable: {device}")
    backend = "ligand_mpnn" if len(ligand_types) else "protein_mpnn"
    weight = CHECKPOINTS[backend]
    # Execute exactly the pinned, verified bytes; never import editable third-party code by path.
    from types import ModuleType
    upstream = ModuleType("aqua_verified_ligandmpnn")
    exec(compile(verified_bytes("model_utils.py"), "LigandMPNN/model_utils.py", "exec"), upstream.__dict__)
    checkpoint = torch.load(io.BytesIO(verified_bytes(weight)), map_location="cpu", weights_only=True)
    n_context = checkpoint.get("atom_context_num", 1)
    model = upstream.ProteinMPNN(k_neighbors=checkpoint["num_edges"], atom_context_num=n_context,
                                 model_type=backend, device=torch.device(device)).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()

    # Both checkpoints use N, CA, C, O in this order, and the same 20-letter alphabet plus X.
    xyz = np.array([[r.atoms[name] for name in ("N", "CA", "C", "O")] for r in trace.residues], dtype=np.float32)
    length = len(trace.sequence)
    numbers = np.array([r.residue_number for r in trace.residues])
    indices = numbers + np.cumsum(np.r_[0, numbers[1:] == numbers[:-1]])
    def tensor(array, dtype):
        return torch.as_tensor(array, dtype=dtype, device=device).unsqueeze(0)
    features = {
        "X": tensor(xyz, torch.float32), "S": tensor([ALPHABET.index(a) for a in trace.sequence], torch.long),
        "R_idx": tensor(indices, torch.long), "chain_labels": tensor(np.zeros(length), torch.long),
        "mask": tensor(np.ones(length), torch.float32), "chain_mask": tensor(np.ones(length), torch.float32),
        "symmetry_residues": [[]],
    }
    if backend == "ligand_mpnn":
        # Upstream defines each residue's local context by nearest heavy atoms to virtual C-beta.
        n, ca, c = xyz[:, 0], xyz[:, 1], xyz[:, 2]
        b, d = ca - n, c - ca
        cb = ca - 0.58273431 * np.cross(b, d) + 0.56802827 * b - 0.54067466 * d
        distances = np.linalg.norm(cb[:, None] - ligand_xyz[None], axis=-1)
        near = np.argsort(distances, axis=1, kind="stable")[:, :n_context]
        count = near.shape[1]
        coords = np.zeros((length, n_context, 3), dtype=np.float32)
        elements = np.zeros((length, n_context), dtype=int)
        valid = np.zeros((length, n_context), dtype=np.float32)
        coords[:, :count], elements[:, :count], valid[:, :count] = ligand_xyz[near], ligand_types[near], 1
        features.update(Y=tensor(coords, torch.float32), Y_t=tensor(elements, torch.long), Y_m=tensor(valid, torch.float32))

    # A fixed random decoding order is shared by all sites, except that the evaluated site is last.
    order = np.random.default_rng(seed).permutation(length) + 1
    output, masked_output = [], []
    with torch.inference_mode():
        for start in range(0, len(positions), 4):
            sites = np.asarray(positions[start:start + 4])
            priority = np.tile(order, (len(sites), 1))
            priority[np.arange(len(sites)), sites] = length + 1
            features["batch_size"] = len(sites)
            features["randn"] = torch.as_tensor(priority, dtype=torch.float32, device=device)
            result = model.score(features, use_sequence=True)
            rows = torch.arange(len(sites), device=device)
            target = torch.as_tensor(sites, device=device)
            if not torch.equal(result["decoding_order"][:, -1], target):
                raise RuntimeError("Target site was not decoded last.")
            output.extend(result["log_probs"][rows, target, :20].cpu().numpy().astype(float))
            if backend == "ligand_mpnn":
                masked = dict(features, Y_m=torch.zeros_like(features["Y_m"]))
                control = model.score(masked, use_sequence=True)
                masked_output.extend(control["log_probs"][rows, target, :20].cpu().numpy().astype(float))
    scores = np.asarray(output)
    control_scores = np.asarray(masked_output) if masked_output else None
    if not np.isfinite(scores).all() or (control_scores is not None and not np.isfinite(control_scores).all()):
        raise ValueError("Non-finite output from the inverse-folding model.")
    return scores, control_scores, {"model": backend, "checkpoint": weight, "checkpoint_sha256": FILES[weight][1],
                    "model_code_revision": REVISION, "torch_version": str(torch.__version__), "device": device,
                    "seed": seed, "ligand_heavy_atoms": len(ligand_types)}
