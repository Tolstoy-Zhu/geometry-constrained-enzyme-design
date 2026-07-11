"""Apply the recorded ColabFold 1.6.1 result-download timeout fix."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
from datetime import datetime, timezone
from pathlib import Path


ORIGINAL = "requests.get(f'{host_url}/result/download/{ID}', timeout=6.02, headers=headers)"
PATCHED = "requests.get(f'{host_url}/result/download/{ID}', timeout=(6.02, 180), headers=headers)"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    version = importlib.metadata.version("colabfold")
    if version != "1.6.1":
        raise RuntimeError(f"Timeout patch is validated only for ColabFold 1.6.1, found {version}")
    distribution = importlib.metadata.distribution("colabfold")
    package_root = Path(distribution.locate_file("colabfold"))
    target = package_root / "colabfold.py"
    text = target.read_text(encoding="utf-8")
    before = digest(target)
    provenance_path = Path("/opt/colabfold/COLABFOLD_PATCH_PROVENANCE.json")
    if PATCHED in text:
        status = "already_applied"
        try:
            previous = json.loads(provenance_path.read_text(encoding="utf-8"))
            original_sha256 = previous.get("original_sha256")
        except (OSError, json.JSONDecodeError):
            original_sha256 = None
    elif ORIGINAL in text:
        original_sha256 = before
        updated = text.replace(ORIGINAL, PATCHED, 1)
        temporary = target.with_suffix(".py.tmp")
        temporary.write_text(updated, encoding="utf-8")
        temporary.replace(target)
        status = "applied"
    else:
        raise RuntimeError("Expected ColabFold result-download code was not found; refusing an unverified patch")
    after = digest(target)
    provenance = {
        "component": "colabfold",
        "version": version,
        "patch": "result download read timeout 6.02s -> 180s",
        "algorithm_changes": False,
        "status": status,
        "target": str(target),
        "original_sha256": original_sha256,
        "patched_sha256": after,
        "applied_at": datetime.now(timezone.utc).isoformat(),
    }
    temporary = provenance_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    temporary.replace(provenance_path)
    print(json.dumps(provenance, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
