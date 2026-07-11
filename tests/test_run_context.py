from __future__ import annotations

import json
from pathlib import Path

from bur01_structmine.run_context import RunContext, mark_run_cancelled


def test_mark_run_cancelled_updates_manifest_and_events(tmp_path: Path) -> None:
    config = tmp_path / "config.yaml"
    config.write_text("query: {}\n", encoding="utf-8")
    context = RunContext(tmp_path / "runs", config, {"query": {}}, run_id="cancel-test")
    context.start()
    assert mark_run_cancelled(context.run_dir) is True
    manifest = json.loads(context.manifest_path.read_text(encoding="utf-8"))
    assert manifest["status"] == "cancelled"
    assert manifest["completed_at"]
    last_event = json.loads(context.events_path.read_text(encoding="utf-8").splitlines()[-1])
    assert last_event["state"] == "cancelled"
