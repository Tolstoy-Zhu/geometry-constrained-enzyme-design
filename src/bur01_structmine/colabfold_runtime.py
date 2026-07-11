"""Pinned LocalColabFold command execution, locking, and result parsing."""
from __future__ import annotations

import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from statistics import fmean
from typing import Any, Callable, Sequence

from .models import PredictionModelRecord, PredictionPreset

ProgressWriter = Callable[[str, str, str, float | None, str | Path | None], None]


class ColabFoldUnavailable(RuntimeError):
    pass


class PredictionCancelled(RuntimeError):
    pass


class GpuOutOfMemoryError(RuntimeError):
    pass


class ResourceBusyError(RuntimeError):
    pass


@dataclass(frozen=True)
class PredictionParameters:
    num_models: int
    num_recycle: int
    num_seeds: int
    num_relax: int
    random_seed: int = 0
    model_type: str = "alphafold2_ptm"

    def to_dict(self) -> dict[str, Any]:
        return {
            "model_type": self.model_type,
            "num_models": self.num_models,
            "num_recycle": self.num_recycle,
            "num_seeds": self.num_seeds,
            "num_relax": self.num_relax,
            "random_seed": self.random_seed,
            "rank_by": "plddt",
            "unified_memory": False,
            "relax_device": "cpu",
        }


PRESETS: dict[PredictionPreset, PredictionParameters] = {
    "quick": PredictionParameters(num_models=1, num_recycle=3, num_seeds=1, num_relax=0),
    "standard": PredictionParameters(num_models=5, num_recycle=3, num_seeds=1, num_relax=1),
    "high": PredictionParameters(num_models=5, num_recycle=12, num_seeds=3, num_relax=1),
}


def preset_parameters(preset: PredictionPreset) -> PredictionParameters:
    try:
        return PRESETS[preset]
    except KeyError as exc:  # pragma: no cover - Pydantic prevents this through public APIs.
        raise ValueError(f"Unknown prediction preset: {preset}") from exc


def to_wsl_path(path: str | Path) -> str:
    """Convert a local Windows path to the default WSL mount representation."""
    raw_value = str(path)
    match = re.match(r"^([A-Za-z]):[\\/](.*)$", raw_value)
    if match:
        drive, remainder = match.groups()
        normalized = remainder.replace("\\", "/")
        return f"/mnt/{drive.lower()}/{normalized}"
    value = str(Path(path).resolve())
    return value.replace("\\", "/")


def process_exists(pid: int) -> bool:
    if pid <= 0:
        return False
    if os.name == "nt":
        import ctypes

        process_query_limited_information = 0x1000
        handle = ctypes.windll.kernel32.OpenProcess(  # type: ignore[attr-defined]
            process_query_limited_information,
            False,
            pid,
        )
        if not handle:
            return False
        ctypes.windll.kernel32.CloseHandle(handle)  # type: ignore[attr-defined]
        return True
    try:
        os.kill(pid, 0)
    except OSError:
        return False
    return True


class ProcessFileLock:
    """Atomic process lock with stale-owner recovery."""

    def __init__(self, path: str | Path, purpose: str) -> None:
        self.path = Path(path)
        self.purpose = purpose
        self.acquired = False

    def acquire(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "pid": os.getpid(),
            "purpose": self.purpose,
            "created_at": datetime.now(timezone.utc).isoformat(),
        }
        for _ in range(2):
            try:
                descriptor = os.open(self.path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            except FileExistsError:
                try:
                    owner = json.loads(self.path.read_text(encoding="utf-8"))
                    owner_pid = int(owner.get("pid", -1))
                except (OSError, ValueError, json.JSONDecodeError):
                    owner_pid = -1
                    owner = {}
                if process_exists(owner_pid):
                    raise ResourceBusyError(
                        f"{self.purpose} is already in use by process {owner_pid}"
                    )
                try:
                    self.path.unlink()
                except FileNotFoundError:
                    pass
                continue
            with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, indent=2)
            self.acquired = True
            return
        raise ResourceBusyError(f"Could not acquire {self.purpose} lock: {self.path}")

    def release(self) -> None:
        if not self.acquired:
            return
        try:
            owner = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            owner = {}
        if int(owner.get("pid", -1)) == os.getpid():
            self.path.unlink(missing_ok=True)
        self.acquired = False

    def __enter__(self) -> "ProcessFileLock":
        self.acquire()
        return self

    def __exit__(self, *_: object) -> None:
        self.release()


def build_msa_arguments(input_path: str, output_dir: str, use_templates: bool) -> list[str]:
    arguments = [input_path, output_dir, "--msa-only", "--model-type", "alphafold2_ptm"]
    if use_templates:
        arguments.append("--templates")
    return arguments


def build_prediction_arguments(
    input_path: str,
    output_dir: str,
    parameters: PredictionParameters,
    use_templates: bool,
) -> list[str]:
    arguments = [
        input_path,
        output_dir,
        "--model-type",
        parameters.model_type,
        "--num-models",
        str(parameters.num_models),
        "--num-recycle",
        str(parameters.num_recycle),
        "--num-seeds",
        str(parameters.num_seeds),
        "--random-seed",
        str(parameters.random_seed),
        "--num-relax",
        str(parameters.num_relax),
        "--rank",
        "plddt",
        "--disable-unified-memory",
    ]
    if use_templates:
        arguments.append("--templates")
    return arguments


@dataclass(frozen=True)
class ColabFoldRuntime:
    backend: str = "wsl"
    executable: str = "/opt/colabfold/bin/colabfold_batch"
    wsl_distribution: str = "Ubuntu-Foldseek"

    def _command(
        self,
        arguments: Sequence[str],
        linux_pid_path: str | None = None,
    ) -> list[str]:
        if self.backend == "native":
            if Path(self.executable).suffix.lower() == ".py":
                return [sys.executable, self.executable, *arguments]
            return [self.executable, *arguments]
        if self.backend != "wsl":
            raise ValueError(f"Unsupported ColabFold backend: {self.backend}")
        if not shutil.which("wsl.exe"):
            raise ColabFoldUnavailable("wsl.exe was not found")
        command = shlex.join([self.executable, *arguments])
        if linux_pid_path:
            command = f"echo $$ > {shlex.quote(linux_pid_path)}; exec {command}"
        return [
            "wsl.exe",
            "-d",
            self.wsl_distribution,
            "--",
            "bash",
            "-lc",
            command,
        ]

    def path(self, path: str | Path) -> str:
        return to_wsl_path(path) if self.backend == "wsl" else str(Path(path).resolve())

    def _terminate(self, process: subprocess.Popen[str], pid_path: Path | None) -> None:
        if self.backend == "wsl" and pid_path and pid_path.exists():
            try:
                linux_pid = int(pid_path.read_text(encoding="utf-8").strip())
                subprocess.run(
                    ["wsl.exe", "-d", self.wsl_distribution, "--", "kill", "-TERM", str(linux_pid)],
                    capture_output=True,
                    timeout=10,
                    check=False,
                )
            except (OSError, ValueError, subprocess.SubprocessError):
                pass
        process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()

    def run(
        self,
        arguments: Sequence[str],
        log_path: str | Path,
        cancel_path: str | Path,
        pid_path: str | Path,
        timeout: float | None = None,
    ) -> None:
        log_path = Path(log_path)
        cancel_path = Path(cancel_path)
        pid_path = Path(pid_path)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        pid_path.unlink(missing_ok=True)
        linux_pid_path = self.path(pid_path) if self.backend == "wsl" else None
        command = self._command(arguments, linux_pid_path=linux_pid_path)
        started = time.monotonic()
        creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
        with log_path.open("a", encoding="utf-8") as log_handle:
            log_handle.write(f"\n$ {shlex.join(command)}\n")
            log_handle.flush()
            process = subprocess.Popen(
                command,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                errors="replace",
                creationflags=creationflags,
            )
            while process.poll() is None:
                if cancel_path.exists():
                    self._terminate(process, pid_path)
                    raise PredictionCancelled("Prediction cancelled by user")
                if timeout is not None and time.monotonic() - started > timeout:
                    self._terminate(process, pid_path)
                    raise TimeoutError(f"ColabFold command exceeded {timeout:.0f} seconds")
                time.sleep(0.5)
            return_code = int(process.returncode or 0)
        pid_path.unlink(missing_ok=True)
        if return_code == 0:
            return
        tail = log_path.read_text(encoding="utf-8", errors="replace")[-20000:]
        lowered = tail.lower()
        if any(marker in lowered for marker in ("resource_exhausted", "out of memory", "cuda_error_out_of_memory", "failed to allocate")):
            raise GpuOutOfMemoryError(
                "ColabFold exhausted GPU memory. Reduce the sequence length or use the quick preset."
            )
        raise RuntimeError(f"ColabFold exited with code {return_code}; see {log_path}")

    def gpu_details(self) -> dict[str, Any]:
        if os.environ.get("AQUA_SKIP_GPU_CHECK") == "1":
            return {"platform": "test", "devices": ["fake-gpu"]}
        python = "/opt/colabfold/env/bin/python" if self.backend == "wsl" else shutil.which("python") or "python"
        script = (
            "import json,jax; d=jax.devices(); "
            "print(json.dumps({'jax':jax.__version__,'devices':[str(x) for x in d],"
            "'platforms':[x.platform for x in d]}))"
        )
        if self.backend == "wsl":
            command = ["wsl.exe", "-d", self.wsl_distribution, "--", python, "-c", script]
        else:
            command = [python, "-c", script]
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=60, check=False)
        except (OSError, subprocess.SubprocessError) as exc:
            raise ColabFoldUnavailable(str(exc)) from exc
        if result.returncode != 0:
            raise ColabFoldUnavailable((result.stderr or result.stdout).strip())
        details = json.loads(result.stdout.strip().splitlines()[-1])
        if "gpu" not in details.get("platforms", []):
            raise ColabFoldUnavailable("JAX did not detect a GPU; CPU fallback is disabled")
        return details

    def versions(self) -> dict[str, str]:
        if self.backend == "native" and os.environ.get("AQUA_FAKE_COLABFOLD_VERSION"):
            return {"colabfold": os.environ["AQUA_FAKE_COLABFOLD_VERSION"]}
        python = "/opt/colabfold/env/bin/python" if self.backend == "wsl" else shutil.which("python") or "python"
        script = (
            "import importlib.metadata as m,json; "
            "print(json.dumps({k:m.version(k) for k in ['colabfold','jax','openmm','alphafold-colabfold']}))"
        )
        command = (
            ["wsl.exe", "-d", self.wsl_distribution, "--", python, "-c", script]
            if self.backend == "wsl"
            else [python, "-c", script]
        )
        result = subprocess.run(command, capture_output=True, text=True, timeout=60, check=False)
        if result.returncode != 0:
            raise ColabFoldUnavailable((result.stderr or result.stdout).strip())
        return json.loads(result.stdout.strip().splitlines()[-1])

    def patch_provenance(self) -> list[dict[str, Any]]:
        if self.backend == "native":
            return []
        result = subprocess.run(
            [
                "wsl.exe",
                "-d",
                self.wsl_distribution,
                "--",
                "cat",
                "/opt/colabfold/COLABFOLD_PATCH_PROVENANCE.json",
            ],
            capture_output=True,
            text=True,
            errors="replace",
            timeout=30,
            check=False,
        )
        if result.returncode != 0:
            return []
        return [json.loads(result.stdout)]

    def parameter_checksums(self) -> dict[str, str]:
        if self.backend == "native" and os.environ.get("AQUA_FAKE_COLABFOLD_VERSION"):
            return {"fake_parameters": "test-only"}
        if self.backend != "wsl":
            return {}
        provenance = self.parameter_provenance()
        recorded = provenance.get("parameter_files", {})
        if isinstance(recorded, dict) and recorded:
            return {str(name): str(checksum) for name, checksum in recorded.items()}

        script = (
            "import hashlib,json,pathlib; "
            "files=sorted(pathlib.Path('/opt/colabfold/data').rglob('params_model_*.npz')); "
            "print(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}))"
        )
        result = subprocess.run(
            [
                "wsl.exe",
                "-d",
                self.wsl_distribution,
                "--",
                "/opt/colabfold/env/bin/python",
                "-c",
                script,
            ],
            capture_output=True,
            text=True,
            errors="replace",
            timeout=30 * 60,
            check=False,
        )
        if result.returncode != 0:
            raise ColabFoldUnavailable((result.stderr or result.stdout).strip())
        return {str(name): str(checksum) for name, checksum in json.loads(result.stdout).items()}

    def parameter_provenance(self) -> dict[str, Any]:
        if self.backend == "native" and os.environ.get("AQUA_FAKE_COLABFOLD_VERSION"):
            return {"source_url": "test-only", "archive_sha256": "test-only"}
        if self.backend != "wsl":
            return {}
        result = subprocess.run(
            [
                "wsl.exe",
                "-d",
                self.wsl_distribution,
                "--",
                "cat",
                "/opt/colabfold/ALPHAFOLD_PARAMETER_PROVENANCE.json",
            ],
            capture_output=True,
            text=True,
            errors="replace",
            timeout=30,
            check=False,
        )
        if result.returncode != 0:
            return {}
        return json.loads(result.stdout)


def parse_prediction_models(output_dir: str | Path) -> list[PredictionModelRecord]:
    output_dir = Path(output_dir)
    score_files = sorted(output_dir.rglob("*scores_rank_*.json"))
    records: list[PredictionModelRecord] = []
    for scores_path in score_files:
        rank_match = re.search(r"rank_(\d{3})", scores_path.name)
        if not rank_match:
            continue
        rank = int(rank_match.group(1))
        rank_token = f"rank_{rank:03d}"
        pdb_candidates = [path for path in output_dir.rglob(f"*{rank_token}*.pdb") if "_relaxed_" in path.name]
        relaxed = bool(pdb_candidates)
        if not pdb_candidates:
            pdb_candidates = [path for path in output_dir.rglob(f"*{rank_token}*.pdb") if "_unrelaxed_" in path.name]
        if not pdb_candidates:
            continue
        scores = json.loads(scores_path.read_text(encoding="utf-8"))
        plddt = [float(value) for value in scores.get("plddt", [])]
        mean_plddt = fmean(plddt) if plddt else None
        ptm = scores.get("ptm")
        max_pae = scores.get("max_pae")
        suffix = scores_path.stem.split(rank_token, 1)[-1].lstrip("_")
        records.append(
            PredictionModelRecord(
                rank=rank,
                model_name=suffix or scores_path.stem,
                pdb_path=str(pdb_candidates[0].resolve()),
                scores_path=str(scores_path.resolve()),
                mean_plddt=mean_plddt,
                ptm=float(ptm) if ptm is not None else None,
                max_pae=float(max_pae) if max_pae is not None else None,
                relaxed=relaxed,
            )
        )
    return sorted(records, key=lambda record: record.rank)
