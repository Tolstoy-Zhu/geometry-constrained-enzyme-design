from __future__ import annotations

import hashlib
import io
import tarfile

import pytest

from bur01_structmine.foldseek_runner import (
    _extract_pinned_database_archive,
    _validate_pinned_archive,
)


def _write_fixture_archive(path) -> dict[str, bytes]:
    files = {
        "fixture_db": b"coordinates",
        "fixture_db.dbtype": b"\x00\x00\x00\x00",
        "fixture_db.index": b"index",
    }
    with tarfile.open(path, "w:gz") as bundle:
        for name, payload in files.items():
            member = tarfile.TarInfo(name)
            member.size = len(payload)
            bundle.addfile(member, io.BytesIO(payload))
    return files


def test_validate_and_extract_pinned_database_archive(tmp_path):
    archive = tmp_path / "fixture_db.tar.gz"
    files = _write_fixture_archive(archive)
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()

    assert _validate_pinned_archive(archive, archive.stat().st_size, digest) == digest

    database_dir = tmp_path / "databases"
    database_dir.mkdir()
    target = database_dir / "renamed_db"
    _extract_pinned_database_archive(archive, target, "fixture_db")

    for source_name, payload in files.items():
        suffix = source_name.removeprefix("fixture_db")
        assert (database_dir / f"renamed_db{suffix}").read_bytes() == payload


def test_pinned_database_archive_rejects_bad_size(tmp_path):
    archive = tmp_path / "fixture_db.tar.gz"
    _write_fixture_archive(archive)

    with pytest.raises(RuntimeError, match="size mismatch"):
        _validate_pinned_archive(archive, archive.stat().st_size + 1, "0" * 64)
