"""Generate a portable result viewer from verified property tables."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


def report_data(directory):
    base = Path(directory)
    manifest = json.loads((base / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("method") != "aqua_properties_v1":
        raise ValueError("Unsupported property report format.")

    def table(name):
        path = base / name
        if hashlib.sha256(path.read_bytes()).hexdigest() != manifest["output_sha256"].get(name):
            raise ValueError(f"Property results have changed: {name}")
        with path.open(encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    proteins = {r["sequence_sha256"]: r for r in table("protein_predictions.tsv")}
    pairs = {
        (r["sequence_sha256"], r["substrate"]): float(r["enzrank_score"])
        for r in table("substrate_predictions.tsv")
    }
    substrate_names = [r["name"] for r in manifest["substrates"]]

    def profile(h):
        p = proteins[h]
        return {
            "swi": float(p["swi_soluble_score"]),
            "tm": float(p["tmprot_tm_c"]) if p["tmprot_tm_c"] else None,
            "enz": {s: pairs[h, s] for s in substrate_names},
            "length": int(p["length"]),
        }

    rows, parents = {}, {}
    mode = manifest.get("comparison_mode", "design_parent")
    if mode == "absolute":
        for source, info in manifest["source_runs"].items():
            for record in info["records"]:
                cid = record["id"]
                rows[(source, cid)] = dict(id=cid, run=source, **profile(record["sequence_sha256"]))
    for entry in table("candidate_comparison.tsv") if mode != "absolute" else []:
        run, cid, h, parent = (entry[k] for k in ("run", "candidate_id", "sequence_sha256", "parent_sha256"))
        key = (run, cid)
        if key not in rows:
            rows[key] = dict(
                id=cid,
                run=run,
                mutations=entry.get("mutations_8hfw") or entry["mutations_sequence"],
                reference="8HFW:A" if entry.get("mutations_8hfw") else "sequence",
                original_mutations=entry["mutations_pdb"],
                sequence_mutations=entry["mutations_sequence"],
                edits=int(entry["substitutions"]),
                parent=parent,
                **profile(h),
            )
        if parent not in parents:
            parents[parent] = profile(parent)
    return {
        "mode": mode,
        "runId": manifest["run_id"],
        "sources": manifest["source_runs"],
        "rows": list(rows.values()),
        "parents": parents,
        "substrates": manifest["substrates"],
        "pairs": manifest["pair_count"],
    }


def write_viewer(directory):
    payload = json.dumps(report_data(directory), ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    payload = payload.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    template = Path(__file__).with_name("property_viewer.html").read_text(encoding="utf-8")
    if template.count("__AQUA_PROPERTY_DATA__") != 1:
        raise ValueError("Invalid property viewer template.")
    path = Path(directory) / "index.html"
    path.write_text(template.replace("__AQUA_PROPERTY_DATA__", payload), encoding="utf-8")
    return path
