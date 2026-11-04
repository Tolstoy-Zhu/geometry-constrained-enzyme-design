"""Isolated inference worker. Accepts local paths only; downloads belong to properties-setup."""

from __future__ import annotations

import argparse
import ast
import importlib.metadata
import json
import math
import os
from pathlib import Path

PROTEIN_CODES = "AILVFWYNCQMSTDERHKGPOUXBZ"  # Maranas group EnzRank, index 0 is left padding
CANONICAL = frozenset("ACDEFGHIKLMNPQRSTVWY")


def swi_parameters(path):
    # Parse only literal fitted coefficients from the hash-checked author distribution.
    params = {}
    for line in ast.parse(Path(path).read_text(encoding="utf-8")).body:
        if isinstance(line, ast.Assign) and len(line.targets) == 1 and isinstance(line.targets[0], ast.Name):
            name = line.targets[0].id
            if name in {"weights", "A", "B"}:
                params[name] = ast.literal_eval(line.value)
    if set(params) != {"weights", "A", "B"} or set(params["weights"]) != CANONICAL:
        raise ValueError("Invalid SWI model parameters.")
    return params


def predict_solubility(sequence, params):
    swi = math.fsum(params["weights"][aa] for aa in sequence) / len(sequence)
    probability = 1.0 / (1.0 + math.exp(-(params["A"] * swi + params["B"])))
    return {"swi_index": swi, "swi_soluble_score": probability}


def thermal(sequences, root, device):
    import torch
    from peft import PeftModel, get_peft_model_state_dict
    from safetensors.torch import load_file
    from transformers import AutoTokenizer, EsmForSequenceClassification

    if device == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA was requested but is unavailable in the property runtime.")
    torch.set_num_threads(min(6, os.cpu_count() or 1))
    tokenizer = AutoTokenizer.from_pretrained(str(root / "esm"), local_files_only=True)
    base = EsmForSequenceClassification.from_pretrained(
        str(root / "esm"), num_labels=1, local_files_only=True
    )
    model = PeftModel.from_pretrained(base, str(root / "adapter"), local_files_only=True)
    # The ESM-2 backbone alone has a randomly initialized regression head. Verify the full fitted adapter and head.
    published = load_file(str(root / "adapter/adapter_model.safetensors"))
    loaded = get_peft_model_state_dict(model)
    if any(k not in loaded or not torch.equal(v, loaded[k].cpu()) for k, v in published.items()):
        raise ValueError("TmProt regression head or LoRA weights were not fully restored.")
    model = model.to(device).eval()
    out = {k: None for k in sequences}
    eligible = [
        k for k, s in sequences.items() if len(s) <= 510
    ]  # Two of the authors' 512 tokens are special tokens.
    with torch.inference_mode():
        for idx in range(0, len(eligible), 2):
            keys = eligible[idx : idx + 2]
            batch = tokenizer(
                [sequences[k] for k in keys], return_tensors="pt", padding=True, truncation=False
            )
            if batch.input_ids.shape[1] > 512:
                raise ValueError("Unexpected ESM token count.")
            scores = model(**batch.to(device)).logits.flatten().cpu().tolist()
            if len(scores) != len(keys) or not all(math.isfinite(v) for v in scores):
                raise ValueError("Invalid TmProt regression output.")
            out.update(zip(keys, scores))
    return out


def predict_pairs(sequences, substrates, root):
    import numpy as np
    import tensorflow as tf
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator

    tf.config.set_visible_devices([], "GPU")
    fingerprint = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    chemical, fps, seen = [], [], set()
    for name, smiles in substrates:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or mol.GetNumHeavyAtoms() == 0 or len(Chem.GetMolFrags(mol)) != 1:
            raise ValueError(f"Provide one valid, connected substrate SMILES for {name}.")
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        if canonical in seen:
            raise ValueError(f"Duplicate substrate structure: {name}.")
        seen.add(canonical)
        chemical.append({"name": name, "input_smiles": smiles, "canonical_smiles": canonical})
        fps.append(fingerprint.GetFingerprintAsNumPy(mol))

    exported = tf.saved_model.load(str(root / "enzyme"))
    predict = exported.signatures["serving_default"]
    inputs = predict.structured_input_signature[1]
    by_width = {int(spec.shape[-1]): key for key, spec in inputs.items()}
    if len(inputs) != 2 or set(by_width) != {2048, 2500}:
        raise ValueError("Unexpected EnzRank input signature.")
    code = {aa: i for i, aa in enumerate(PROTEIN_CODES, 1)}
    padded = {}
    for key, seq in sequences.items():
        if len(seq) > 2500:
            raise ValueError("EnzRank supports up to 2500 residues without truncation.")
        array = np.zeros(2500, dtype=np.int32)
        array[-len(seq) :] = [code[a] for a in seq]
        padded[key] = array
    combinations = [(key, j) for key in sequences for j in range(len(chemical))]
    output = []
    for start in range(0, len(combinations), 16):
        batch = combinations[start : start + 16]
        args = {}
        for width, values in ((2048, [fps[j] for _, j in batch]), (2500, [padded[key] for key, _ in batch])):
            key = by_width[width]
            args[key] = tf.convert_to_tensor(np.asarray(values), dtype=inputs[key].dtype)
        result = predict(**args)
        if len(result) != 1:
            raise ValueError("Unexpected EnzRank output signature.")
        values = next(iter(result.values())).numpy()
        if values.shape != (len(batch), 1) or not np.all((values >= 0.0) & (values <= 1.0)):
            raise ValueError("Invalid EnzRank output shape or score.")
        for (key, j), value in zip(batch, values[:, 0]):
            output.append(
                {"sequence_sha256": key, "substrate": chemical[j]["name"], "enzrank_score": float(value)}
            )
    return {"substrates": chemical, "pairs": output}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=["thermal", "pairs"])
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("models", type=Path)
    parser.add_argument("device", choices=["cpu", "cuda"])
    args = parser.parse_args()
    request = json.loads(args.input.read_text(encoding="utf-8"))
    if args.mode == "thermal":
        result = {"tm": thermal(request["sequences"], args.models, args.device)}
        packages = ["torch", "transformers", "peft", "safetensors"]
    else:
        result = predict_pairs(request["sequences"], request["substrates"], args.models)
        packages = ["tensorflow", "rdkit", "numpy"]
    result["versions"] = {k: importlib.metadata.version(k) for k in packages}
    args.output.write_text(json.dumps(result, ensure_ascii=False, allow_nan=False), encoding="utf-8")


if __name__ == "__main__":
    main()
