"""Local, parent-matched property screening of one or more AquaDesign runs."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import subprocess
import tempfile
from pathlib import Path

from .config import project_root
from .property_assets import (
    ASSETS,
    ENZ_REV,
    ESM_REV,
    SWI_REV,
    TM_REV,
    file_hash,
    runtime_python,
    verify_assets,
)
from .property_inference import CANONICAL, predict_solubility, swi_parameters
from .run_context import make_run_id
from .structure_search import load_ca_trace

DEFAULT_SUBSTRATES = (("octanoic_acid", "CCCCCCCC(=O)O"), ("ethanol", "CCO"))


def profiles_root():
    return Path(os.environ.get("AQUA_PROPERTY_ROOT") or project_root() / "results/properties")


def sequence_hash(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def read_table(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_table(path, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_substrates(specs):
    pairs = []
    for spec in specs:
        name, sep, smiles = spec.partition("=")
        name, smiles = name.strip(), smiles.strip()
        if not sep or not name or not smiles or any(ord(c) < 32 for c in name + smiles):
            raise ValueError("Use a separate name=SMILES entry for each substrate.")
        if name in {key for key, _ in pairs}:
            raise ValueError(f"Duplicate substrate name: {name}")
        pairs.append((name, smiles))
    if not 1 <= len(pairs) <= 32:
        raise ValueError("Supply 1–32 individual substrates.")
    return pairs


def collect_runs(directories):
    sequences, samples, members, sources = {}, {}, [], {}
    ref = load_ca_trace(project_root() / "data/query/Bur01m.pdb", chain_id="A")
    for directory in map(Path, directories):
        label = directory.name
        if not label or label in sources:
            raise ValueError(f"Each input design run must have a distinct name: {label}")
        manifest_path, table_path = directory / "manifest.json", directory / "candidates.tsv"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        parent_path = (directory / manifest["parent_file"]).resolve()
        if parent_path.parent != directory.resolve() or file_hash(parent_path) != manifest["parent_sha256"]:
            raise ValueError(f"Invalid parent coordinates in {label}.")
        trace = load_ca_trace(parent_path, chain_id=manifest["chain"])
        parent = trace.sequence
        parent_hash = sequence_hash(parent)
        if set(parent) - CANONICAL or parent_hash != manifest["parent_sequence_sha256"]:
            raise ValueError(f"Parent sequence mismatch in {label}.")
        sequences[parent_hash] = parent
        samples.setdefault(parent_hash, {"id": f"parent_{parent_hash[:16]}", "parent": False, "runs": []})[
            "parent"
        ] = True
        samples[parent_hash]["runs"].append(label)
        seen = set()
        rows = read_table(table_path)
        if not rows or manifest["candidates"] != len(rows):
            raise ValueError(f"Candidate count differs from the manifest in {label}.")
        for row in rows:
            sequence, digest = row["sequence"], row["sequence_sha256"]
            if (
                not sequence
                or set(sequence) - CANONICAL
                or len(sequence) != len(parent)
                or sequence_hash(sequence) != digest
                or row["candidate_id"] != f"AD_{digest[:16]}"
                or digest in seen
                or digest == parent_hash
            ):
                raise ValueError(
                    f"Inconsistent candidate sequence or identifier in {label}: {row['candidate_id']}"
                )
            seen.add(digest)
            changed = [i for i, (a, b) in enumerate(zip(sequence, parent)) if a != b]
            edits = len(changed)
            if edits != int(row["substitutions"]):
                raise ValueError(f"Incorrect mutation count in {label}: {row['candidate_id']}")
            sequence_label = ",".join(f"{parent[i]}{i + 1}{sequence[i]}" for i in changed)
            source_label = ",".join(trace.residues[i].label + sequence[i] for i in changed)
            if sequence_label != row["mutations_sequence"] or source_label != row["mutations_pdb"]:
                raise ValueError(f"Inconsistent mutation annotation in {label}: {row['candidate_id']}")
            ref_label = (
                ",".join(ref.residues[i].label + sequence[i] for i in changed)
                if parent == ref.sequence
                else ""
            )
            sequences[digest] = sequence
            sample = samples.setdefault(digest, {"id": row["candidate_id"], "parent": False, "runs": []})
            sample["runs"].append(label)
            members.append(
                dict(
                    run=label,
                    candidate_id=row["candidate_id"],
                    mutations_pdb=row["mutations_pdb"],
                    mutations_sequence=sequence_label,
                    mutations_8hfw=ref_label,
                    substitutions=edits,
                    sequence_sha256=digest,
                    parent_sha256=parent_hash,
                )
            )
        sources[label] = {
            "manifest_sha256": file_hash(manifest_path),
            "table_sha256": file_hash(table_path),
            "parent_sha256": parent_hash,
            "design_model": manifest.get("model", {}).get("model"),
            "candidate_count": len(rows),
            "numbering_reference": "8HFW:A" if parent == ref.sequence else "input_sequence",
        }
    if not members:
        raise ValueError("Select at least one non-empty AquaDesign run.")
    return sequences, samples, members, sources


def run_worker(mode, request, folder, assets, device):
    python = runtime_python()
    if not python.is_file():
        raise ValueError("Property runtime not installed. Run properties-setup --install-runtime.")
    path = folder / (mode + ".json")
    env = dict(
        os.environ,
        HF_HUB_OFFLINE="1",
        TRANSFORMERS_OFFLINE="1",
        HF_HUB_DISABLE_TELEMETRY="1",
        USE_TF="0",
        TF_CPP_MIN_LOG_LEVEL="3",
        TF_NUM_INTRAOP_THREADS="4",
        TF_NUM_INTEROP_THREADS="1",
        TOKENIZERS_PARALLELISM="false",
        PYTHONNOUSERSITE="1",
    )
    # Keep editable packages and unrelated site directories out of the model interpreter.
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    result = subprocess.run(
        [
            str(python),
            str(Path(__file__).with_name("property_inference.py")),
            mode,
            str(request),
            str(path),
            str(assets),
            device,
        ],
        env=env,
        capture_output=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    if result.returncode:
        raise ValueError(
            f"{mode} inference failed: " + (result.stderr or result.stdout or "No runtime message.")[-1600:]
        )
    return json.loads(path.read_text(encoding="utf-8"))


def assemble(sequences, samples, members, substrates, tms, pair_result, params):
    if set(tms) != set(sequences):
        raise ValueError("TmProt returned the wrong set of sequences.")
    chemistry = pair_result["substrates"]
    if [(r["name"], r["input_smiles"]) for r in chemistry] != list(substrates):
        raise ValueError("EnzRank substrate identities do not match the request.")
    pairs = pair_result["pairs"]
    by_pair = {(r["sequence_sha256"], r["substrate"]): r["enzrank_score"] for r in pairs}
    expected = {(h, name) for h in sequences for name, _ in substrates}
    if (
        set(by_pair) != expected
        or len(by_pair) != len(pairs)
        or not all(0 <= v <= 1 for v in by_pair.values())
    ):
        raise ValueError("Incomplete, duplicate or invalid enzyme–substrate predictions.")
    proteins = {}
    for digest, sequence in sequences.items():
        tm = tms[digest]
        if (tm is None) != (len(sequence) > 510) or (tm is not None and not math.isfinite(tm)):
            raise ValueError(f"Invalid TmProt result for {digest}.")
        sample = samples[digest]
        proteins[digest] = dict(
            sample_id=sample["id"],
            is_parent=sample["parent"],
            source_runs="|".join(dict.fromkeys(sample["runs"])),
            sequence_sha256=digest,
            length=len(sequence),
            **predict_solubility(sequence, params),
            tmprot_tm_c=tm,
            tmprot_status="ok" if tm is not None else "exceeds_510_residues",
            sequence=sequence,
        )
    pair_table = []
    for h, entry in proteins.items():
        for chem in chemistry:
            name = chem["name"]
            pair_table.append(
                dict(
                    sample_id=entry["sample_id"],
                    sequence_sha256=h,
                    substrate=name,
                    input_smiles=chem["input_smiles"],
                    canonical_smiles=chem["canonical_smiles"],
                    enzrank_score=by_pair[h, name],
                )
            )
    comparisons = []
    for member in members:
        h, p = member["sequence_sha256"], member["parent_sha256"]
        a, b = proteins[h], proteins[p]
        for name, _ in substrates:
            delta_tm = (
                a["tmprot_tm_c"] - b["tmprot_tm_c"]
                if a["tmprot_tm_c"] is not None and b["tmprot_tm_c"] is not None
                else None
            )
            comparisons.append(
                dict(
                    **member,
                    substrate=name,
                    swi_soluble_score=a["swi_soluble_score"],
                    parent_swi_soluble_score=b["swi_soluble_score"],
                    delta_swi=a["swi_soluble_score"] - b["swi_soluble_score"],
                    tmprot_tm_c=a["tmprot_tm_c"],
                    parent_tmprot_tm_c=b["tmprot_tm_c"],
                    delta_tm_c=delta_tm,
                    enzrank_score=by_pair[h, name],
                    parent_enzrank_score=by_pair[p, name],
                    delta_enzrank=by_pair[h, name] - by_pair[p, name],
                )
            )
    return list(proteins.values()), pair_table, comparisons


def profile_runs(directories, substrates=DEFAULT_SUBSTRATES, output=None, *, device="cpu", progress=None):
    if device not in {"cpu", "cuda"}:
        raise ValueError("Choose cpu or cuda for TmProt.")
    substrates = parse_substrates([name + "=" + smiles for name, smiles in substrates])
    sequences, samples, members, sources = collect_runs(directories)
    target = Path(output) if output is not None else profiles_root() / make_run_id("properties")
    if target.exists():
        raise FileExistsError(f"Output already exists: {target}")

    def status(message):
        if progress is not None:
            progress(message)

    status("正在核对模型权重 / Checking model weights")
    assets = verify_assets()
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=target.parent) as temp:
        folder = Path(temp)
        report = folder / "complete"
        report.mkdir()
        request = folder / "input.json"
        request.write_text(
            json.dumps({"sequences": sequences, "substrates": substrates}, ensure_ascii=False),
            encoding="utf-8",
        )
        status("EnzRank · 正在比较酶与底物 / Scoring enzyme–substrate pairs")
        pairs = run_worker("pairs", request, folder, assets, device)
        status("TmProt · 正在预测熔解温度 / Predicting melting temperatures")
        thermal = run_worker("thermal", request, folder, assets, device)
        status("SWI · 正在计算可溶表达分数及亲本差值 / Comparing to each parent")
        proteins, pair_table, comparisons = assemble(
            sequences,
            samples,
            members,
            substrates,
            thermal["tm"],
            pairs,
            swi_parameters(assets / "swi/swi.py"),
        )
        for filename, rows in (
            ("protein_predictions.tsv", proteins),
            ("substrate_predictions.tsv", pair_table),
            ("candidate_comparison.tsv", comparisons),
        ):
            write_table(report / filename, rows)
        manifest = dict(
            method="aqua_properties_v1",
            run_id=target.name,
            source_runs=sources,
            unique_sequences=len(sequences),
            design_records=len(members),
            pair_count=len(pair_table),
            substrates=pairs["substrates"],
            thermal_device=device,
            versions=thermal["versions"] | pairs["versions"],
            revisions={"SWI": SWI_REV, "TmProt": TM_REV, "ESM2": ESM_REV, "EnzRank": ENZ_REV},
            asset_sha256=ASSETS,
            output_sha256={f.name: file_hash(f) for f in report.glob("*.tsv")},
        )
        (report / "manifest.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        from .property_report import write_viewer

        write_viewer(report)
        report.rename(target)
    return target
