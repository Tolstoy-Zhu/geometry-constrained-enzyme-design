from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import pytest


@pytest.fixture
def fake_mmseqs_executable(tmp_path: Path) -> Path:
    source = Path(__file__).parent / "fixtures" / "fake_mmseqs.py"
    if os.name == "nt":
        wrapper = tmp_path / "fake_mmseqs.cmd"
        wrapper.write_text(
            f'@echo off\r\n"{sys.executable}" "{source}" %*\r\n',
            encoding="ascii",
        )
        return wrapper
    destination = tmp_path / "fake_mmseqs"
    shutil.copy2(source, destination)
    destination.chmod(0o755)
    return destination
