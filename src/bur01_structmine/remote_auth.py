"""Password verification for the optional remote Streamlit workbench."""
from __future__ import annotations

import base64
import hashlib
import hmac
import json
from pathlib import Path
from typing import Any, Mapping


MIN_ITERATIONS = 100_000
MAX_ITERATIONS = 2_000_000
MAX_PASSWORD_LENGTH = 256


def load_remote_access(path: str | Path) -> dict[str, Any] | None:
    """Load an enabled remote-access record, rejecting malformed local state."""

    record_path = Path(path)
    if not record_path.is_file():
        return None
    try:
        record = json.loads(record_path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(record, dict) or record.get("enabled") is not True:
        return None
    return record


def verify_remote_password(password: str, record: Mapping[str, Any]) -> bool:
    """Verify a PBKDF2-SHA256 password record in constant time."""

    if not password or len(password) > MAX_PASSWORD_LENGTH:
        return False
    try:
        iterations = int(record["iterations"])
        if not MIN_ITERATIONS <= iterations <= MAX_ITERATIONS:
            return False
        salt = base64.b64decode(str(record["salt_b64"]), validate=True)
        expected = base64.b64decode(str(record["hash_b64"]), validate=True)
    except (KeyError, TypeError, ValueError):
        return False
    if len(salt) < 16 or len(expected) != 32:
        return False
    actual = hashlib.pbkdf2_hmac("sha256", password.encode("utf-8"), salt, iterations, dklen=32)
    return hmac.compare_digest(actual, expected)
