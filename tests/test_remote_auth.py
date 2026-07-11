import base64
import hashlib
import json
from pathlib import Path

from streamlit.testing.v1 import AppTest

from bur01_structmine.config import project_root
from bur01_structmine.remote_auth import load_remote_access, verify_remote_password


def _record(password: str, *, enabled: bool = True, iterations: int = 200_000) -> dict[str, object]:
    salt = b"0123456789abcdef"
    digest = hashlib.pbkdf2_hmac("sha256", password.encode(), salt, iterations, dklen=32)
    return {
        "enabled": enabled,
        "iterations": iterations,
        "salt_b64": base64.b64encode(salt).decode(),
        "hash_b64": base64.b64encode(digest).decode(),
    }


def test_remote_password_verification() -> None:
    record = _record("correct-horse-battery-staple")

    assert verify_remote_password("correct-horse-battery-staple", record)
    assert not verify_remote_password("wrong", record)
    assert not verify_remote_password("", record)


def test_remote_password_rejects_unsafe_or_malformed_records() -> None:
    assert not verify_remote_password("secret", _record("secret", iterations=10))
    assert not verify_remote_password("secret", {"iterations": 200_000})
    assert not verify_remote_password("x" * 257, _record("secret"))


def test_remote_access_loader_requires_enabled_valid_json(tmp_path) -> None:
    path = tmp_path / "access.json"
    path.write_text(json.dumps(_record("secret", enabled=False)), encoding="utf-8")
    assert load_remote_access(path) is None

    expected = _record("secret")
    path.write_text(json.dumps(expected), encoding="utf-8")
    assert load_remote_access(path) == expected

    path.write_text("not-json", encoding="utf-8")
    assert load_remote_access(path) is None


def test_streamlit_remote_gate_unlocks_the_workbench(monkeypatch, tmp_path: Path) -> None:
    password = "expert-demo-password"
    access_path = tmp_path / "access.json"
    access_path.write_text(json.dumps(_record(password)), encoding="utf-8")
    monkeypatch.setenv("AQUA_REMOTE_ACCESS_PATH", str(access_path))
    monkeypatch.setenv("AQUA_RUN_ROOT", str(tmp_path / "runs"))
    monkeypatch.setenv("AQUA_PREDICTION_ROOT", str(tmp_path / "predictions"))

    test = AppTest.from_file(str(project_root() / "ui" / "app.py"), default_timeout=60).run(timeout=60)
    assert not test.exception
    assert not test.tabs
    assert any(item.value == "在线演示访问 / Secure remote access" for item in test.subheader)

    test.text_input[0].input(password)
    test.button[0].click().run(timeout=60)

    assert not test.exception
    assert [tab.label for tab in test.tabs] == ["运行", "AF预测", "候选", "结构", "历史", "诊断"]
