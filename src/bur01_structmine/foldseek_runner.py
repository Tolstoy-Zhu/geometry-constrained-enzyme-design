"""Foldseek command wrapper with native and WSL2 execution backends."""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import tarfile
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, List, Sequence
from urllib.parse import urlparse

from .io_utils import ensure_dir
from .process_text import capture_process_text, utf8_process_env

FOLDSEEK_FIELDS: List[str] = [
    "query",
    "target",
    "evalue",
    "bits",
    "prob",
    "alntmscore",
    "qtmscore",
    "ttmscore",
    "rmsd",
    "fident",
    "alnlen",
    "qlen",
    "tlen",
    "qcov",
    "tcov",
    "tseq",
    "theader",
]


@dataclass(frozen=True)
class FoldseekRuntime:
    backend: str
    executable: str
    wsl_distribution: str = "Ubuntu-Foldseek"


class FoldseekUnavailable(RuntimeError):
    pass


def _wsl_foldseek_works(distribution: str, executable: str) -> bool:
    if not shutil.which("wsl.exe"):
        return False
    try:
        result = capture_process_text(
            ["wsl.exe", "-d", distribution, "--", executable, "version"],
            timeout=20,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return result.returncode == 0


def resolve_foldseek_runtime(
    executable: str = "foldseek",
    backend: str = "auto",
    wsl_distribution: str = "Ubuntu-Foldseek",
    wsl_executable: str = "/opt/foldseek/bin/foldseek",
) -> FoldseekRuntime:
    """Resolve a usable Foldseek runtime without shell-string execution."""
    backend = str(backend or "auto").strip().lower()
    if backend not in {"auto", "native", "wsl"}:
        raise ValueError(f"Unsupported Foldseek backend: {backend}")

    native_path = shutil.which(executable)
    if backend in {"auto", "native"} and native_path:
        return FoldseekRuntime("native", native_path, wsl_distribution)
    if backend == "native":
        raise FoldseekUnavailable(f"Native Foldseek executable not found: {executable}")

    if backend in {"auto", "wsl"} and _wsl_foldseek_works(wsl_distribution, wsl_executable):
        return FoldseekRuntime("wsl", wsl_executable, wsl_distribution)
    raise FoldseekUnavailable(
        "Foldseek is unavailable. Run `python start_here.py --action doctor` and configure "
        f"the native executable or WSL distribution {wsl_distribution!r}."
    )


def check_foldseek_executable(
    executable: str = "foldseek",
    backend: str = "auto",
    wsl_distribution: str = "Ubuntu-Foldseek",
    wsl_executable: str = "/opt/foldseek/bin/foldseek",
) -> bool:
    try:
        resolve_foldseek_runtime(executable, backend, wsl_distribution, wsl_executable)
    except FoldseekUnavailable:
        return False
    return True


def _wsl_path(path: str | Path, runtime: FoldseekRuntime) -> str:
    resolved = str(Path(path).resolve())
    result = capture_process_text(
        ["wsl.exe", "-d", runtime.wsl_distribution, "--", "wslpath", "-a", resolved],
        timeout=20,
        check=True,
    )
    converted = result.stdout.strip()
    if not converted:
        raise RuntimeError(f"WSL path conversion returned no value for {resolved}")
    return converted


def _runtime_path(path: str | Path, runtime: FoldseekRuntime) -> str:
    return _wsl_path(path, runtime) if runtime.backend == "wsl" else str(Path(path).resolve())


def run_command(cmd: Iterable[str], runtime: FoldseekRuntime) -> None:
    command = [str(part) for part in cmd]
    if runtime.backend == "wsl":
        command = ["wsl.exe", "-d", runtime.wsl_distribution, "--", runtime.executable, *command]
    else:
        command = [runtime.executable, *command]
    print("[CMD]", " ".join(command))
    subprocess.run(command, env=utf8_process_env(), check=True)


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_pinned_archive(path: Path, expected_size: int, expected_sha256: str) -> str:
    observed_size = path.stat().st_size
    if observed_size != expected_size:
        raise RuntimeError(
            f"Database archive size mismatch for {path}: expected {expected_size}, observed {observed_size}"
        )
    observed_sha256 = _file_sha256(path)
    if observed_sha256.lower() != expected_sha256.lower():
        raise RuntimeError(
            f"Database archive SHA-256 mismatch for {path}: "
            f"expected {expected_sha256}, observed {observed_sha256}"
        )
    return observed_sha256


def _download_with_curl(url: str, destination: Path) -> None:
    curl = shutil.which("curl.exe") or shutil.which("curl")
    if not curl:
        raise RuntimeError("curl is required for resumable pinned database downloads")
    partial = destination.with_name(destination.name + ".part")
    command = [
        curl,
        "--fail",
        "--location",
        "--retry",
        "10",
        "--retry-all-errors",
        "--continue-at",
        "-",
        "--output",
        str(partial),
        url,
    ]
    print("[CMD]", " ".join(command))
    subprocess.run(command, env=utf8_process_env(), check=True)
    partial.replace(destination)


def _extract_pinned_database_archive(archive: Path, target_db: Path, archive_prefix: str) -> None:
    if not archive_prefix or Path(archive_prefix).name != archive_prefix:
        raise ValueError(f"Invalid database archive prefix: {archive_prefix!r}")
    with tempfile.TemporaryDirectory(prefix=f".{target_db.name}.extract-", dir=target_db.parent) as temp_name:
        staging = Path(temp_name)
        with tarfile.open(archive, mode="r:gz") as bundle:
            members = bundle.getmembers()
            if not members:
                raise RuntimeError(f"Database archive is empty: {archive}")
            for member in members:
                if not member.isfile() or Path(member.name).name != member.name:
                    raise RuntimeError(f"Unsafe or unsupported archive member: {member.name}")
                if not member.name.startswith(archive_prefix):
                    raise RuntimeError(
                        f"Archive member {member.name!r} does not match prefix {archive_prefix!r}"
                    )
                source = bundle.extractfile(member)
                if source is None:
                    raise RuntimeError(f"Could not read archive member: {member.name}")
                with source, (staging / member.name).open("wb") as destination:
                    shutil.copyfileobj(source, destination, length=1024 * 1024)

        staged = sorted(staging.iterdir())
        destinations = [
            target_db.parent / f"{target_db.name}{path.name[len(archive_prefix):]}" for path in staged
        ]
        existing = [path for path in destinations if path.exists()]
        if existing:
            raise FileExistsError(
                "Refusing to replace existing database files: " + ", ".join(str(path) for path in existing)
            )
        for source, destination in zip(staged, destinations):
            source.replace(destination)

    if not Path(f"{target_db}.dbtype").is_file():
        raise RuntimeError(f"Installed Foldseek database is missing its dbtype file: {target_db}.dbtype")


def _write_download_provenance(
    target_db: Path,
    archive: Path,
    archive_url: str,
    archive_sha256: str,
    version_url: str,
    version_text: str,
) -> Path:
    out = target_db.parent / "environment" / f"{target_db.name}_download_provenance.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    record = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "target_prefix": str(target_db.resolve()),
        "archive_url": archive_url,
        "archive_path": str(archive.resolve()),
        "archive_size": archive.stat().st_size,
        "archive_sha256": archive_sha256,
        "version_url": version_url,
        "database_version": version_text,
    }
    tmp = out.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(record, indent=2, ensure_ascii=False), encoding="utf-8")
    tmp.replace(out)
    return out


def get_foldseek_version(
    executable: str = "foldseek",
    backend: str = "auto",
    wsl_distribution: str = "Ubuntu-Foldseek",
    wsl_executable: str = "/opt/foldseek/bin/foldseek",
) -> str:
    runtime = resolve_foldseek_runtime(executable, backend, wsl_distribution, wsl_executable)
    command = [runtime.executable, "version"]
    if runtime.backend == "wsl":
        command = ["wsl.exe", "-d", runtime.wsl_distribution, "--", *command]
    result = capture_process_text(command, timeout=20, check=True)
    return (result.stdout or result.stderr).strip()


def run_easy_search(
    query_pdb: str | Path,
    target_db: str | Path,
    out_tsv: str | Path,
    tmp_dir: str | Path,
    executable: str = "foldseek",
    max_seqs: int = 1000,
    threads: int = 8,
    format_fields: Sequence[str] | None = None,
    cluster_search: bool | None = None,
    extra_args: Sequence[str] | None = None,
    backend: str = "auto",
    wsl_distribution: str = "Ubuntu-Foldseek",
    wsl_executable: str = "/opt/foldseek/bin/foldseek",
) -> Path:
    """Run Foldseek easy-search through the configured backend."""
    query_pdb = Path(query_pdb)
    target_db = Path(target_db)
    out_tsv = Path(out_tsv)
    tmp_dir = Path(tmp_dir)
    ensure_dir(out_tsv.parent)
    ensure_dir(tmp_dir)
    if not query_pdb.exists():
        raise FileNotFoundError(f"Query PDB not found: {query_pdb}")

    runtime = resolve_foldseek_runtime(executable, backend, wsl_distribution, wsl_executable)
    cmd = [
        "easy-search",
        _runtime_path(query_pdb, runtime),
        _runtime_path(target_db, runtime),
        _runtime_path(out_tsv, runtime),
        _runtime_path(tmp_dir, runtime),
        "--max-seqs",
        str(max_seqs),
        "--threads",
        str(threads),
        "--format-output",
        ",".join(format_fields or FOLDSEEK_FIELDS),
    ]
    if cluster_search is not None:
        cmd.extend(["--cluster-search", "1" if cluster_search else "0"])
    if extra_args:
        cmd.extend(str(value) for value in extra_args)
    run_command(cmd, runtime)
    return out_tsv


def download_foldseek_database(
    database_name: str,
    target_db: str | Path,
    tmp_dir: str | Path,
    executable: str = "foldseek",
    backend: str = "auto",
    wsl_distribution: str = "Ubuntu-Foldseek",
    wsl_executable: str = "/opt/foldseek/bin/foldseek",
    archive_url: str = "",
    version_url: str = "",
    archive_sha256: str = "",
    archive_size: int | None = None,
    archive_prefix: str = "",
) -> Path:
    target_db = Path(target_db)
    tmp_dir = Path(tmp_dir)
    ensure_dir(target_db.parent)
    ensure_dir(tmp_dir)
    dbtype = Path(f"{target_db}.dbtype")
    if archive_url:
        if not archive_sha256 or archive_size is None or not archive_prefix:
            raise ValueError(
                "Pinned database downloads require archive_sha256, archive_size, and archive_prefix"
            )
        cache_dir = ensure_dir(tmp_dir / "official_downloads")
        archive_name = Path(urlparse(archive_url).path).name
        if not archive_name:
            raise ValueError(f"Could not determine archive filename from URL: {archive_url}")
        archive = cache_dir / archive_name
        if not archive.exists():
            _download_with_curl(archive_url, archive)
        observed_sha256 = _validate_pinned_archive(archive, archive_size, archive_sha256)

        version_file = Path(f"{target_db}.version")
        if version_url and not version_file.exists():
            cached_version = cache_dir / Path(urlparse(version_url).path).name
            if not cached_version.exists():
                _download_with_curl(version_url, cached_version)
            shutil.copy2(cached_version, version_file)
        version_text = version_file.read_text(encoding="utf-8").strip() if version_file.exists() else ""

        if not dbtype.exists():
            _extract_pinned_database_archive(archive, target_db, archive_prefix)
        _write_download_provenance(
            target_db,
            archive,
            archive_url,
            observed_sha256,
            version_url,
            version_text,
        )
        return target_db

    runtime = resolve_foldseek_runtime(executable, backend, wsl_distribution, wsl_executable)
    run_command(
        [
            "databases",
            database_name,
            _runtime_path(target_db, runtime),
            _runtime_path(tmp_dir, runtime),
        ],
        runtime,
    )
    return target_db


def create_foldseek_database(
    input_path: str | Path,
    target_db: str | Path,
    executable: str = "foldseek",
    prostt5_model: str | Path | None = None,
    use_gpu: bool = False,
    extra_args: Sequence[str] | None = None,
    backend: str = "auto",
    wsl_distribution: str = "Ubuntu-Foldseek",
    wsl_executable: str = "/opt/foldseek/bin/foldseek",
) -> Path:
    input_path = Path(input_path)
    target_db = Path(target_db)
    ensure_dir(target_db.parent)
    if not input_path.exists():
        raise FileNotFoundError(f"Foldseek database input not found: {input_path}")
    runtime = resolve_foldseek_runtime(executable, backend, wsl_distribution, wsl_executable)
    cmd = ["createdb", _runtime_path(input_path, runtime), _runtime_path(target_db, runtime)]
    if prostt5_model:
        cmd.extend(["--prostt5-model", _runtime_path(prostt5_model, runtime)])
    if use_gpu:
        cmd.extend(["--gpu", "1"])
    if extra_args:
        cmd.extend(str(value) for value in extra_args)
    run_command(cmd, runtime)
    return target_db
