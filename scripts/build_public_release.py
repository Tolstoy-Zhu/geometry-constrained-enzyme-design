"""Build the sanitized v0.3 release archive after local acceptance passes."""
from __future__ import annotations

import json
import re
import zipfile
from datetime import datetime, timezone
from pathlib import Path

from bur01_structmine.config import project_root
from bur01_structmine.run_context import file_sha256

EXCLUDED_PARTS = {
    ".cache",
    ".git",
    ".obsidian",
    ".pytest_cache",
    "__pycache__",
    "databases",
    "Github-Bur01",
    "results",
    "tmp",
}
EXCLUDED_SUFFIXES = {".a3m", ".appx", ".npz", ".zip"}
TEXT_SUFFIXES = {".bat", ".json", ".md", ".ps1", ".py", ".sh", ".txt", ".yaml", ".yml"}
MAX_PUBLIC_FILE_SIZE = 25 * 1024 * 1024


def public_files(root: Path):
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(root)
        if any(part in EXCLUDED_PARTS or part.endswith(".egg-info") for part in relative.parts):
            continue
        if path.suffix.lower() in EXCLUDED_SUFFIXES:
            continue
        if path.stat().st_size > MAX_PUBLIC_FILE_SIZE:
            raise RuntimeError(f"Public file exceeds 25 MB: {relative}")
        yield path, relative


def audit_public_text(root: Path, files: list[tuple[Path, Path]]) -> None:
    forbidden = [
        re.compile(r"[A-Za-z]:\\Users\\", re.IGNORECASE),
        re.compile(re.escape(str(root)), re.IGNORECASE),
        re.compile(r"\\\\wsl(?:\.localhost)?\\", re.IGNORECASE),
    ]
    findings: list[str] = []
    for path, relative in files:
        if path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for pattern in forbidden:
            if pattern.search(text):
                findings.append(f"{relative}: {pattern.pattern}")
    if findings:
        raise RuntimeError("Machine-specific paths remain in public files:\n" + "\n".join(findings))


def main() -> int:
    root = project_root()
    validation_path = root / "docs" / "validation" / "v0.3_validation_summary.json"
    if not validation_path.is_file():
        raise RuntimeError(f"Missing sanitized validation summary: {validation_path}")
    validation = json.loads(validation_path.read_text(encoding="utf-8"))
    if not validation.get("passed"):
        raise RuntimeError("Sanitized validation summary is not passing")
    screenshots = [
        root / "docs" / "screenshots" / "af_prediction_desktop.png",
        root / "docs" / "screenshots" / "af_prediction_mobile.png",
    ]
    missing_screenshots = [str(path) for path in screenshots if not path.is_file()]
    if missing_screenshots:
        raise RuntimeError("Missing public UI screenshot(s): " + ", ".join(missing_screenshots))

    files = list(public_files(root))
    audit_public_text(root, files)
    stamp = datetime.now().strftime("%Y%m%d")
    archive = root.parent / f"AquaEsterNet-StructMine_{stamp}_v0.3_verified.zip"
    if archive.exists():
        archive = root.parent / f"AquaEsterNet-StructMine_{stamp}_v0.3_verified_{datetime.now().strftime('%H%M%S')}.zip"

    index = {
        "project_version": "0.3.0",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "sanitized": True,
        "excludes": sorted(EXCLUDED_PARTS),
        "validation": validation,
        "files": [
            {"path": str(relative).replace("\\", "/"), "sha256": file_sha256(path), "size_bytes": path.stat().st_size}
            for path, relative in files
        ],
    }
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as bundle:
        for path, relative in files:
            bundle.write(path, Path("AquaEsterNet-StructMine") / relative)
        bundle.writestr(
            "AquaEsterNet-StructMine/VERIFICATION_INDEX.json",
            json.dumps(index, indent=2, ensure_ascii=False),
        )
    digest = file_sha256(archive)
    checksum_path = archive.with_suffix(".sha256")
    checksum_path.write_text(f"{digest}  {archive.name}\n", encoding="ascii")
    print(f"Archive: {archive}")
    print(f"SHA-256: {digest}")
    print(f"Checksum: {checksum_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
