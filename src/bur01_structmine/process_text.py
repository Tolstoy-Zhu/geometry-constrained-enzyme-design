"""Text exchanged with local tools and displayed in the workbench."""

from __future__ import annotations

import codecs
import locale
import os
import subprocess
from pathlib import Path


def utf8_process_env() -> dict[str, str]:
    # A redirected file handle does not set the child process's own encoding.
    return {**os.environ, "PYTHONIOENCODING": "utf-8", "WSL_UTF8": "1"}


def decode_process_output(data: bytes) -> str:
    """WSL host diagnostics can be UTF-16 even when the Linux tool uses UTF-8."""
    if data.startswith((codecs.BOM_UTF16_LE, codecs.BOM_UTF16_BE)):
        return data.decode("utf-16", errors="backslashreplace")
    if b"\x00" in data:
        return data.decode("utf-16-le", errors="backslashreplace")
    try:
        return data.decode("utf-8-sig")
    except UnicodeDecodeError:
        # Older native tools use the Windows ANSI code page for piped output.
        encoding = "mbcs" if os.name == "nt" else locale.getpreferredencoding(False)
        return data.decode(encoding, errors="backslashreplace")


def capture_process_text(
    command: list[str] | str, *, timeout: float | None = None, check: bool = False
) -> subprocess.CompletedProcess[str]:
    raw = subprocess.run(command, capture_output=True, env=utf8_process_env(), timeout=timeout, check=False)
    result = subprocess.CompletedProcess(
        raw.args, raw.returncode, decode_process_output(raw.stdout), decode_process_output(raw.stderr)
    )
    if check:
        result.check_returncode()
    return result


def read_process_log(path: Path) -> str:
    data = path.read_bytes()
    if b"\x00" not in data or data.startswith(codecs.BOM_UTF16_BE):
        return decode_process_output(data)

    # Old Windows logs can contain UTF-8 command lines around UTF-16LE WSL output.
    # Only an aligned UTF-16 newline ends a wide line (e.g., 上 also contains 0x0a).
    lines = []
    start = 0
    previous_wide = False
    while start < len(data):
        lf = data.find(b"\n", start)
        end = len(data) if lf < 0 else lf + 1
        line = data[start:end]
        wide = b"\x00" in line or (lf >= 0 and data[end : end + 1] == b"\x00")
        if previous_wide and not wide:
            try:
                line.decode("utf-8")
            except UnicodeDecodeError:
                wide = True
        if wide:
            lf = data.find(b"\n\x00", start)
            while lf != -1 and (lf - start) % 2:
                lf = data.find(b"\n\x00", lf + 1)
            end = len(data) if lf < 0 else lf + 2
            lines.append(
                data[start:end].decode("utf-16-le", errors="backslashreplace").removeprefix("\ufeff")
            )
        else:
            lines.append(decode_process_output(line))
        previous_wide = wide
        start = end
    return "".join(lines)
