from __future__ import annotations

import json
from pathlib import Path

import requests
from streamlit.testing.v1 import AppTest

from bur01_structmine import demo_service
from bur01_structmine.config import project_root
from bur01_structmine.run_context import file_sha256


def test_demo_five_steps_no_compute_or_http(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Read-only demo attempted a process or external request")
    monkeypatch.setattr("subprocess.run", forbidden)
    monkeypatch.setattr("subprocess.Popen", forbidden)
    monkeypatch.setattr(requests.Session, "request", forbidden)
    root = project_root()
    manifest_path = root / "data/demo/demo_manifest.json"
    before = file_sha256(manifest_path)
    demo = AppTest.from_file(str(root / "ui/demo.py"), default_timeout=60).run()
    for step in range(5):
        demo.radio[0].set_value(step).run()
        assert not demo.exception
    assert before == file_sha256(manifest_path)


def test_demo_stop_rejects_unrelated_pid(monkeypatch, tmp_path):
    path = tmp_path / "state.json"
    path.write_text(json.dumps({"pid": 42, "token": "demo-test"}))
    monkeypatch.setattr(demo_service, "state_path", lambda: path)
    monkeypatch.setattr(demo_service, "_process_command", lambda pid: "unrelated program")
    assert demo_service.stop_demo()["status"] == "not_running_or_identity_mismatch"
    assert path.exists()


def test_demo_asset_hashes():
    root = project_root()
    manifest = json.loads((root / "data/demo/demo_manifest.json").read_text())
    assert all(item["passed"] for item in manifest["benchmark_checks"])
    for relative, digest in manifest["files"].items():
        assert file_sha256(root / relative) == digest
    assert not any("results/" in name for name in manifest["files"])
