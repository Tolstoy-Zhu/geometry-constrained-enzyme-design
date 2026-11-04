import json

import pytest

from bur01_structmine import property_screen as ps
from bur01_structmine.config import project_root
from bur01_structmine.property_assets import file_hash
from bur01_structmine.property_report import report_data
from bur01_structmine.structure_search import load_ca_trace


def make_run(root, name, shift=0):
    source = project_root() / "data/query/Bur01m.pdb"
    run = root / name
    run.mkdir()
    path = run / "parent.pdb"
    lines = []
    for line in source.read_text().splitlines():
        if line.startswith("ATOM"):
            line = line[:22] + f"{int(line[22:26]) + shift:4d}" + line[26:]
        lines.append(line)
    path.write_text("\n".join(lines) + "\n")
    tr = load_ca_trace(path, chain_id="A")
    before = tr.sequence[0]
    variant_aa = "C" if before != "C" else "A"
    variant = variant_aa + tr.sequence[1:]
    h = ps.sequence_hash(variant)
    row = dict(
        candidate_id="AD_" + h[:16],
        sequence=variant,
        sequence_sha256=h,
        substitutions=1,
        mutations_sequence=before + "1" + variant_aa,
        mutations_pdb=tr.residues[0].label + variant_aa,
    )
    ps.write_table(run / "candidates.tsv", [row])
    (run / "manifest.json").write_text(
        json.dumps(
            dict(
                parent_file=path.name,
                parent_sha256=file_hash(path),
                parent_sequence_sha256=ps.sequence_hash(tr.sequence),
                chain="A",
                candidates=1,
                model={"model": "protein_mpnn"},
            )
        )
    )
    return run


def test_exact_parent_and_mutation_numbering(tmp_path):
    a, b = make_run(tmp_path, "first"), make_run(tmp_path, "second", shift=-11)
    sequences, _, members, _ = ps.collect_runs([a, b])
    assert len(sequences) == 2
    assert members[0]["mutations_pdb"] != members[1]["mutations_pdb"]
    assert members[0]["mutations_8hfw"] == members[1]["mutations_8hfw"] == members[0]["mutations_pdb"]
    f = b / "candidates.tsv"
    row = ps.read_table(f)[0]
    row["mutations_sequence"] = "X999Z"
    ps.write_table(f, [row])
    with pytest.raises(ValueError, match="mutation annotation"):
        ps.collect_runs([b])


def test_scores_bind_to_ids_and_to_each_parent_not_row_order(tmp_path, monkeypatch):
    a, b = make_run(tmp_path, "first"), make_run(tmp_path, "second", shift=-11)
    assets = tmp_path / "weights"
    (assets / "swi").mkdir(parents=True)
    from bur01_structmine.property_inference import CANONICAL

    (assets / "swi/swi.py").write_text(f"weights={dict.fromkeys(CANONICAL, 0.75)!r}\nA=1.0\nB=0.0\n")
    monkeypatch.setattr(ps, "verify_assets", lambda: assets)

    def fixture_worker(mode, request, *_):
        spec = json.loads(request.read_text(encoding="utf-8"))
        keys = list(spec["sequences"])
        if mode == "thermal":
            return {"tm": {keys[1]: 42.0, keys[0]: 40.0}, "versions": {}}
        subs = spec["substrates"]
        return {
            "substrates": [dict(name=n, input_smiles=s, canonical_smiles=s) for n, s in subs],
            "pairs": list(
                reversed(
                    [
                        dict(sequence_sha256=h, substrate=n, enzrank_score=0.1 + 0.2 * i + 0.1 * j)
                        for i, h in enumerate(keys)
                        for j, (n, _) in enumerate(subs)
                    ]
                )
            ),
            "versions": {},
        }

    monkeypatch.setattr(ps, "run_worker", fixture_worker)
    output = ps.profile_runs([a, b], output=tmp_path / "report")
    frame = ps.read_table(output / "candidate_comparison.tsv")
    assert len(frame) == 4
    assert all(float(r["delta_tm_c"]) == 2 for r in frame)
    assert all(float(r["delta_enzrank"]) == pytest.approx(0.2) for r in frame)
    assert all(float(r["delta_swi"]) == 0 for r in frame)
    data = report_data(output)
    assert len(data["rows"]) == 2 and len(data["parents"]) == 1
    assert "__AQUA_PROPERTY_DATA__" not in (output / "index.html").read_text(encoding="utf-8")
    with pytest.raises(FileExistsError):
        ps.profile_runs([a, b], output=output)
    with (output / "substrate_predictions.tsv").open("a") as fh:
        fh.write("corrupt")
    with pytest.raises(ValueError, match="changed"):
        report_data(output)


def test_incomplete_inference_does_not_publish_output(tmp_path, monkeypatch):
    run = make_run(tmp_path, "input")
    monkeypatch.setattr(ps, "verify_assets", lambda: tmp_path)

    def fail(*args):
        raise ValueError("model unavailable")

    monkeypatch.setattr(ps, "run_worker", fail)
    with pytest.raises(ValueError, match="model unavailable"):
        ps.profile_runs([run], output=tmp_path / "report")
    assert not (tmp_path / "report").exists()


def test_invalid_substrates_and_overlength_missing_values():
    assert ps.parse_substrates(["acid=CC(=O)O"]) == [("acid", "CC(=O)O")]
    with pytest.raises(ValueError, match="Duplicate"):
        ps.parse_substrates(["x=CC", "x=CCO"])
    with pytest.raises(ValueError):
        ps.parse_substrates(["water"])
    from bur01_structmine.property_inference import CANONICAL

    seq = "A" * 511
    props = dict(weights=dict.fromkeys(CANONICAL, 0.75), A=1.0, B=0.0)
    ids = {"hash": {"id": "parent", "parent": True, "runs": ["r"]}}
    q = {
        "substrates": [{"name": "a", "input_smiles": "CC", "canonical_smiles": "CC"}],
        "pairs": [{"sequence_sha256": "hash", "substrate": "a", "enzrank_score": 0.0}],
    }
    rows, pairs, _ = ps.assemble({"hash": seq}, ids, [], [("a", "CC")], {"hash": None}, q, props)
    assert rows[0]["tmprot_status"] == "exceeds_510_residues" and pairs[0]["enzrank_score"] == 0
    with pytest.raises(ValueError, match="TmProt"):
        ps.assemble({"hash": seq}, ids, [], [("a", "CC")], {"hash": 0}, q, props)
    with pytest.raises(ValueError, match="Incomplete"):
        ps.assemble({"hash": "ACDE"}, ids, [], [("a", "CC")], {"hash": 40.0}, dict(q, pairs=[]), props)
