"""On-demand loopback-only demo process, without WSL, GPU or autostart."""
from __future__ import annotations

import json
import os
import socket
import subprocess
import sys
import time
import urllib.request
import uuid
import webbrowser
from pathlib import Path

from .config import project_root


def state_path() -> Path:
    return project_root() / ".cache" / "demo_service" / "state.json"


def _process_command(pid: int) -> str:
    if os.name == "nt":
        command = f"(Get-CimInstance Win32_Process -Filter 'ProcessId = {int(pid)}').CommandLine"
        result = subprocess.run(["powershell.exe", "-NoProfile", "-Command", command], capture_output=True, text=True,
                                errors="replace", timeout=15, creationflags=subprocess.CREATE_NO_WINDOW)
        return result.stdout.strip()
    try:
        return Path(f"/proc/{pid}/cmdline").read_text().replace("\0", " ")
    except OSError:
        return ""


def _owned(state: dict) -> bool:
    command = _process_command(int(state.get("pid", 0)))
    token = state.get("token", "")
    script = str((project_root() / "ui" / "demo.py").resolve())
    return bool(token and token in command and script in command and "streamlit" in command)


def start_demo(port: int | None = None, no_browser: bool = False) -> dict:
    path = state_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        state = json.loads(path.read_text(encoding="utf-8"))
        if _owned(state):
            if not no_browser:
                webbrowser.open(state["url"])
            return state
    with socket.socket() as listener:
        listener.bind(("127.0.0.1", port or 0))
        selected_port = listener.getsockname()[1]
    token = "demo-" + uuid.uuid4().hex
    command = [sys.executable, "-m", "streamlit", "run", str(project_root() / "ui" / "demo.py"),
               "--server.address", "127.0.0.1", "--server.port", str(selected_port), "--server.headless", "true",
               "--browser.gatherUsageStats", "false", "--server.fileWatcherType", "none", "--server.baseUrlPath", token]
    environment = dict(os.environ, PYTHONPATH=str(project_root() / "src"), CUDA_VISIBLE_DEVICES="", AQUA_DEMO_MODE="1")
    log = path.parent / f"{token}.log"
    with log.open("w", encoding="utf-8") as handle:
        process = subprocess.Popen(command, cwd=project_root(), env=environment, stdout=handle, stderr=subprocess.STDOUT,
                                   creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0)
    state = {"pid": process.pid, "token": token, "port": selected_port, "url": f"http://127.0.0.1:{selected_port}/{token}/", "log": str(log)}
    path.write_text(json.dumps(state, indent=2), encoding="utf-8")
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    for _ in range(60):
        if process.poll() is not None:
            raise RuntimeError(f"Demo failed to start. See {log}")
        try:
            with opener.open(state["url"] + "_stcore/health", timeout=1) as response:
                if response.status == 200:
                    if not no_browser:
                        webbrowser.open(state["url"])
                    return state
        except OSError:
            time.sleep(0.25)
    stop_demo()
    raise RuntimeError(f"Demo startup timed out. See {log}")


def stop_demo() -> dict:
    path = state_path()
    if not path.exists():
        return {"status": "not_running"}
    state = json.loads(path.read_text(encoding="utf-8"))
    if not _owned(state):
        return {"status": "not_running_or_identity_mismatch", "pid": state.get("pid")}
    pid = int(state["pid"])
    if os.name == "nt":
        subprocess.run(["taskkill.exe", "/PID", str(pid), "/T", "/F"], check=True, capture_output=True, creationflags=subprocess.CREATE_NO_WINDOW)
    else:
        import signal
        os.kill(pid, signal.SIGTERM)
    for _ in range(30):
        if not _owned(state):
            path.unlink(missing_ok=True)
            return {"status": "stopped", "pid": pid}
        time.sleep(0.1)
    raise RuntimeError("Demo process did not stop")
