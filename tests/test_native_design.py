from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from bur01_structmine import design
from bur01_structmine.config import project_root
from bur01_structmine.design_models import ALPHABET, verified_bytes

PARENT = project_root() / "data/query/Bur01m.pdb"


def test_policy_preserves_numbering_contact_motif_and_manual_mask():
    parent = design.protein_trace(PARENT, "A")
    rows, sites, triad = design.site_policy(parent, np.empty((0, 3)), None, ["A240"], 16, 10.0)
    assert [rows[i]["position"] for i in triad] == [100, 230, 263]
    assert [rows[i]["residue"] for i in triad] == ["S111", "D241", "H274"]
    assert len(sites) == 16
    assert all(rows[i]["core_distance_A"] > 4 and rows[i]["status"] == "designable" for i in sites)
    assert "user_fixed" in next(r["status"] for r in rows if r["residue"] == "A240")
    assert sum("serine_motif" in r["status"] for r in rows) == 5
    with pytest.raises(ValueError, match="Unknown fixed"):
        design.site_policy(parent, np.empty((0, 3)), None, ["S999"], 16, 10.0)
    with pytest.raises(ValueError, match="outside"):
        design.site_policy(parent, np.full((1, 3), 1000), None, [], 16, 10.0)


def test_msa_deletes_insertions_ignores_gaps_and_collapses_duplicate_homologs(tmp_path):
    msa = tmp_path / "alignment.a3m"
    msa.write_text(">parent\nACD-E\n>h1\nAaeCD-F\n>h2\nACD-F\n>h3\nAC--G\n>h4\nACD-H\n")
    n, valid, fractions = design.alignment_support(msa, "ACDE")
    assert n == 3
    assert list(valid) == [3, 3, 2, 3]
    assert list(fractions) == [1, 1, 1, 0]
    with pytest.raises(ValueError, match="parent"):
        design.alignment_support(msa, "ACDF")


def test_same_checkpoint_mask_contrast_changes_prioritisation():
    trace = design.protein_trace(PARENT, "A")
    i = next(i for i, r in enumerate(trace.residues) if r.label == "A240")
    values = np.full((1, 20), -8.0)
    for aa, score in {"A": -4.0, "C": -1.0, "F": -2.0, "H": -3.0}.items():
        values[0, ALPHABET.index(aa)] = score
    masked = values.copy()
    masked[0, ALPHABET.index("F")] = -6
    options = {}
    for priority in ("compatibility", "ligand_response"):
        variants, landscape = design.propose(trace, [i], values, masked, 1, 8, priority)
        options[priority] = variants[0]
        assert len(landscape) == 19
    assert options["compatibility"]["mutations_pdb"] == "A240C"
    assert options["ligand_response"]["mutations_pdb"] == "A240F"
    assert options["ligand_response"]["ligand_context_delta"] == 4.0


def test_design_writes_only_complete_runs_and_rejects_parent_coordinates(tmp_path, monkeypatch):
    def local_scores(trace, positions, *_, **__):
        scores = np.tile(np.arange(20), (len(positions), 1))
        return scores, None, {"model": "fixture"}
    monkeypatch.setattr(design, "conditional_scores", local_scores)
    run = tmp_path / "run"
    design.run_design(PARENT, run, count=7, site_limit=8)
    table = pd.read_csv(run / "candidates.tsv", sep="\t")
    assert len(table) == table.sequence_sha256.nunique() == 7
    assert set(table.substitutions) == {1, 2}
    assert table.ligand_context_delta.isna().all()
    assert all(design.sha(seq.encode()) == h for seq, h in zip(table.sequence, table.sequence_sha256))
    with pytest.raises(FileExistsError):
        design.run_design(PARENT, run)
    with pytest.raises(ValueError, match="do not match"):
        design.review_design(run, table.iloc[0].candidate_id, PARENT)
    with pytest.raises(ValueError, match="selected candidate"):
        design.run_design(PARENT, tmp_path / "bad", expected_sequence="M" * 283)
    assert not (tmp_path / "bad").exists()
    assert json.loads((run / "manifest.json").read_text())["model"]["model"] == "fixture"

    # Synthetic coordinate attribution fixture: modify residue labels only, never use it as a fold prediction.
    from Bio.Data.PDBData import protein_letters_1to3
    first = table.iloc[0]
    positions = {int(m[1:-1]): protein_letters_1to3[m[-1]] for m in first.mutations_pdb.split(",")}
    directory = tmp_path / "fixture_prediction"
    directory.mkdir()
    model = directory / "rank_001.pdb"
    lines = []
    for line in PARENT.read_text().splitlines():
        if line.startswith("ATOM  ") and int(line[22:26]) in positions:
            line = line[:17] + positions[int(line[22:26])] + line[20:]
        lines.append(line)
    model.write_text("\n".join(lines) + "\n")
    pm = directory / "prediction_manifest.json"
    pred = dict(status="completed", sequence_sha256=first.sequence_sha256, prediction_id="unit-test", ptm=0.71)
    pm.write_text(json.dumps(pred))
    (directory / "rank_001_plddt.json").write_text(json.dumps({"plddt": [87] * 283}))
    report = design.review_design(run, first.candidate_id, model, origin="predicted", prediction_manifest=pm)
    assert report["af2_prediction_id"] == "unit-test"
    assert report["triad_plddt"] == [87.0] * 3
    assert report["qTM"] == pytest.approx(1.0)
    pm.write_text(json.dumps(dict(pred, sequence_sha256="wrong")))
    with pytest.raises(ValueError, match="AF2 manifest"):
        design.review_design(run, first.candidate_id, model, origin="predicted", prediction_manifest=pm)
    assert len(list((run / "reviews").glob("*.json"))) == 1


def test_network_free_inference_fails_on_missing_asset_and_rejects_tampering(tmp_path, monkeypatch):
    monkeypatch.setenv("AQUA_DESIGN_MODELS", str(tmp_path))
    with pytest.raises(ValueError, match="Missing"):
        verified_bytes("model_utils.py")
    (tmp_path / "model_utils.py").write_text("invalid")
    with pytest.raises(ValueError, match="Checksum"):
        verified_bytes("model_utils.py")


def test_ligand_selection_requires_coordinates_and_uses_exact_residue(tmp_path):
    ligand_pdb = tmp_path / "ligand.pdb"
    ligand_pdb.write_text("HETATM    1  C1  ETH B  17      11.000  12.000  13.000  1.00 20.00           C  \n"
                          "HETATM    2  O   HOH C   1      11.000  12.000  14.000  1.00 20.00           O  \nEND\n")
    coords, elements = design.ligand_context(ligand_pdb, ["B:ETH:17"])
    assert list(elements) == [6] and coords.shape == (1, 3)
    with pytest.raises(ValueError, match="Unknown ligand"):
        design.ligand_context(ligand_pdb, ["C:HOH:1"])
    with pytest.raises(ValueError, match="Ligand-response"):
        design.run_design(PARENT, tmp_path / "output", priority="ligand_response")


def test_failed_wsl_probe_preserves_localised_message(monkeypatch):
    import subprocess
    from bur01_structmine.colabfold_runtime import ColabFoldUnavailable, _probe_json

    failure = subprocess.CompletedProcess(["wsl"], 1, "虚拟机无法启动".encode("utf-16"), b"")
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: failure)
    with pytest.raises(ColabFoldUnavailable, match="虚拟机无法启动"):
        _probe_json(["wsl"])
