from __future__ import annotations

import pandas as pd

from bur01_structmine.config import project_root
from bur01_structmine.models import RunManifest
from bur01_structmine.pipeline import reannotate_run


def test_reannotation_creates_new_run_and_preserves_history(monkeypatch, tmp_path):
    source = tmp_path / "runs" / "old_run"
    source.mkdir(parents=True)
    filtered = source / "filtered.tsv"
    pd.DataFrame([{"target": "8HFW_A", "pdb_id": "8HFW", "qtmscore": 1, "ttmscore": 1,
                   "qcov": 1, "tcov": 1, "rmsd": 0, "final_protein_sequence": "GASAG"}]).to_csv(filtered, sep="\t", index=False)
    manifest = RunManifest(run_id="old_run", config_path="config.yaml", status="completed", query_pdb_id="8HFW", query_chain="A",
                           query_path=str(project_root() / "data/query/Bur01m.pdb"), artifacts={"filtered_tsv": str(filtered)})
    (source / "run_manifest.json").write_text(manifest.model_dump_json(), encoding="utf-8")
    before = {path.name: path.read_bytes() for path in source.iterdir()}
    monkeypatch.setattr("bur01_structmine.pipeline.annotate_hits", lambda frame, **kwargs: frame)
    result = reannotate_run(source)
    assert before == {path.name: path.read_bytes() for path in source.iterdir()}
    assert result.iloc[0].source_run_id == "old_run"
    child = next(path for path in source.parent.iterdir() if path != source)
    assert (child / "run_manifest.json").is_file()
    assert (child / "reannotation_source.json").is_file()
    assert (child / "candidates.annotated.tsv").is_file()
