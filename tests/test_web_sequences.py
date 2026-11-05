import json
from threading import Event

import pytest
from streamlit.testing.v1 import AppTest

from bur01_structmine import property_screen as ps
from bur01_structmine.property_fasta import collect_fasta, parse_property_fasta
from bur01_structmine.property_inference import CANONICAL
from bur01_structmine.property_jobs import PropertyJobs
from bur01_structmine.property_report import report_data


def test_web_opens_new_submission_with_no_saved_runs(tmp_path, monkeypatch):
    from bur01_structmine.config import project_root

    monkeypatch.setenv("AQUA_PROPERTY_ROOT", str(tmp_path / "properties"))
    monkeypatch.setenv("AQUA_DESIGN_ROOT", str(tmp_path / "designs"))
    test = AppTest.from_file(str(project_root() / "ui/properties_app.py")).run(timeout=30)
    assert not test.exception
    assert any(button.label == "提交本地预测" for button in test.button)
    test.text_area(key="web_fasta").input(">parent\nACDE\n>candidate\nACGE")
    test.checkbox(key="web_parent").check().run(timeout=30)
    assert not test.exception
    assert any("已识别 2 条记录" in line.value for line in test.caption)
    test.text_area(key="web_fasta").input(">parent\nACDE\n>short\nACD").run(timeout=30)
    assert any("长度不同" in item.value for item in test.error)


def test_invalid_uploaded_structure_gives_an_input_error(tmp_path):
    from bur01_structmine.design import ligand_groups

    path = tmp_path / "no-coordinates.pdb"
    path.write_text("END\n")
    with pytest.raises(ValueError, match="原子记录"):
        ligand_groups(path)


def test_parent_is_explicit_and_numbering_uses_input_sequence():
    fasta = "\ufeff>WT with description\r\nac de\r\n>variant\r\nac ge\r\n>alias\r\nac de"
    seqs, _, members, sources = collect_fasta(fasta, True)
    assert len(seqs) == 2
    assert [(m["candidate_id"], m["mutations_sequence"], m["substitutions"]) for m in members] == [
        ("variant", "D3G", 1),
        ("alias", "", 0),
    ]
    assert all(not m["mutations_pdb"] and not m["mutations_8hfw"] for m in members)
    assert [r["id"] for r in sources["FASTA"]["records"]] == ["WT", "variant", "alias"]
    assert not collect_fasta(fasta, False)[2]
    assert parse_property_fasta("ac\nd e")[0].sequence == "ACDE"
    with pytest.raises(ValueError, match="长度不同"):
        collect_fasta(">parent\nACDE\n>short\nACD", True)
    with pytest.raises(ValueError, match="至少一条候选"):
        collect_fasta("ACDE", True)


@pytest.mark.parametrize(
    "fasta",
    [
        "",
        ">no-sequence",
        ">x\nACX",
        ">x\nA\n>x\nC",
        "A\n>x\nCD",
        ">p\n" + "A" * 2501,
        "\n".join(f">s{i}\nAC" for i in range(101)),
    ],
)
def test_fasta_input_failures(fasta):
    with pytest.raises(ValueError):
        parse_property_fasta(fasta)


@pytest.mark.parametrize("parent", [False, True])
def test_fasta_results_retain_names_absolute_scores_and_downloads(tmp_path, monkeypatch, parent):
    (tmp_path / "swi").mkdir()
    (tmp_path / "swi/swi.py").write_text(f"weights={dict.fromkeys(CANONICAL, 0.75)!r}\nA=1\nB=0\n")
    monkeypatch.setattr(ps, "verify_assets", lambda: tmp_path)

    def inference(mode, request, *_):
        spec = json.loads(request.read_text())
        if mode == "thermal":
            return {"tm": {k: 42.0 for k in spec["sequences"]}, "versions": {}}
        return {
            "substrates": [dict(name=n, input_smiles=s, canonical_smiles=s) for n, s in spec["substrates"]],
            "pairs": [
                dict(sequence_sha256=k, substrate=n, enzrank_score=0.25)
                for k in spec["sequences"]
                for n, _ in spec["substrates"]
            ],
            "versions": {},
        }

    monkeypatch.setattr(ps, "run_worker", inference)
    fasta = ">WT\nACDE\n>alias\nACDE\n>variant\nACGE\n"
    out = ps.profile_fasta(fasta, [("ethanol", "CCO")], tmp_path / "report", first_as_parent=parent)
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["unique_sequences"] == 2 and manifest["input_records"] == 3
    assert manifest["pair_count"] == 2 and manifest["design_records"] == 0
    assert "input.fasta" in manifest["output_sha256"] and (out / "input.fasta").read_text() == fasta
    assert (out / "candidate_comparison.tsv").exists() == parent
    data = report_data(out)
    assert data["mode"] == ("first_as_parent" if parent else "absolute")
    assert [r["id"] for r in data["rows"]] == (["alias", "variant"] if parent else ["WT", "alias", "variant"])
    assert all(r["enz"]["ethanol"] == 0.25 for r in data["rows"])
    assert len(data["parents"]) == int(parent)


def test_queue_survives_client_navigation_and_isolates_job_failures(tmp_path, monkeypatch):
    from bur01_structmine import property_jobs

    monkeypatch.setattr(property_jobs, "profiles_root", lambda: tmp_path)
    queue = PropertyJobs()
    started, release, ran = Event(), Event(), Event()

    def first(output, progress):
        progress("TmProt")
        started.set()
        assert release.wait(5)
        raise ValueError("bad input")

    def second(output, progress):
        ran.set()

    try:
        a = queue.submit("first", first)
        assert started.wait(5)
        b = queue.submit("second", second)
        assert [(j.key, j.state) for j in queue.snapshot()] == [(b, "queued"), (a, "running")]
    finally:
        release.set()
        queue._pool.shutdown(wait=True)
    assert ran.is_set()
    jobs = queue.snapshot()
    assert [j.state for j in jobs] == ["complete", "failed"] and jobs[1].message == "bad input"
