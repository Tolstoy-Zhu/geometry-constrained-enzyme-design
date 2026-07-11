from __future__ import annotations

from pathlib import Path

from streamlit.testing.v1 import AppTest

from bur01_structmine.config import project_root


def test_streamlit_workbench_loads_six_views(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setenv("AQUA_RUN_ROOT", str(tmp_path / "runs"))
    monkeypatch.setenv("AQUA_PREDICTION_ROOT", str(tmp_path / "predictions"))
    monkeypatch.setenv("AQUA_REMOTE_ACCESS_PATH", str(tmp_path / "remote_access.json"))
    app = project_root() / "ui" / "app.py"
    test = AppTest.from_file(str(app), default_timeout=60).run(timeout=60)
    assert not test.exception
    assert [tab.label for tab in test.tabs] == ["运行", "AF预测", "候选", "结构", "历史", "诊断"]
    assert any(title.value == "AquaEsterNet-StructMine" for title in test.title)
