"""MMseqs2 command construction and native/WSL execution."""
from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from .io_utils import ensure_dir
from .process_text import capture_process_text, utf8_process_env


MMSEQS_FIELDS: tuple[str, ...] = (
    "query",
    "target",
    "pident",
    "fident",
    "nident",
    "alnlen",
    "qlen",
    "tlen",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "qcov",
    "tcov",
    "bits",
    "evalue",
    "qaln",
    "taln",
    "tseq",
    "theader",
)


@dataclass(frozen=True)
class MmseqsRuntime:
    backend: str
    executable: str
    wsl_distribution: str = "Ubuntu-Foldseek"


class MmseqsUnavailable(RuntimeError):
    pass


def _wsl_works(distribution: str, executable: str) -> bool:
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


def resolve_mmseqs_runtime(
    executable: str = "mmseqs",
    backend: str = "auto",
    wsl_distribution: str = "Ubuntu-Foldseek",
    wsl_executable: str = "/opt/colabfold/env/bin/mmseqs",
) -> MmseqsRuntime:
    backend = str(backend or "auto").strip().lower()
    if backend not in {"auto", "native", "wsl"}:
        raise ValueError(f"Unsupported MMseqs2 backend: {backend}")

    native = shutil.which(executable)
    if backend in {"auto", "native"} and native:
        return MmseqsRuntime("native", native, wsl_distribution)
    if backend == "native":
        raise MmseqsUnavailable(f"Native MMseqs2 executable not found: {executable}")
    if _wsl_works(wsl_distribution, wsl_executable):
        return MmseqsRuntime("wsl", wsl_executable, wsl_distribution)
    raise MmseqsUnavailable(
        f"MMseqs2 is unavailable in WSL distribution {wsl_distribution!r}; "
        f"expected {wsl_executable}"
    )


def _wsl_path(path: str | Path, runtime: MmseqsRuntime) -> str:
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


def runtime_path(path: str | Path, runtime: MmseqsRuntime) -> str:
    return _wsl_path(path, runtime) if runtime.backend == "wsl" else str(Path(path).resolve())


def build_easy_search_args(
    query: str,
    target: str,
    output: str,
    tmp: str,
    *,
    threads: int,
    iterations: int,
    max_seqs: int,
    search_type: int,
    evalue: float,
    fields: Sequence[str] = MMSEQS_FIELDS,
    extra_args: Sequence[str] | None = None,
) -> list[str]:
    if threads < 1 or iterations < 1 or max_seqs < 1:
        raise ValueError("threads, iterations, and max_seqs must be positive")
    args = [
        "easy-search",
        query,
        target,
        output,
        tmp,
        "--threads",
        str(threads),
        "--num-iterations",
        str(iterations),
        "--max-seqs",
        str(max_seqs),
        "--search-type",
        str(search_type),
        "-e",
        str(evalue),
        "--format-output",
        ",".join(fields),
    ]
    if extra_args:
        args.extend(str(value) for value in extra_args)
    return args


def _command(runtime: MmseqsRuntime, args: Iterable[str]) -> list[str]:
    parts = [str(value) for value in args]
    if runtime.backend == "wsl":
        return ["wsl.exe", "-d", runtime.wsl_distribution, "--", runtime.executable, *parts]
    if os.name == "nt" and Path(runtime.executable).suffix.lower() in {".cmd", ".bat"}:
        return ["cmd.exe", "/d", "/c", runtime.executable, *parts]
    return [runtime.executable, *parts]


def run_mmseqs_command(
    runtime: MmseqsRuntime,
    args: Iterable[str],
    *,
    log_path: str | Path | None = None,
) -> None:
    command = _command(runtime, args)
    print("[CMD]", " ".join(command))
    if log_path is None:
        result = capture_process_text(command)
        if result.returncode:
            detail = (result.stderr or result.stdout).strip()
            raise subprocess.CalledProcessError(result.returncode, command, result.stdout, detail)
        return

    log_path = Path(log_path)
    ensure_dir(log_path.parent)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write("[CMD] " + " ".join(command) + "\n")
        handle.flush()
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, env=utf8_process_env(), check=True)


def get_mmseqs_version(
    executable: str = "mmseqs",
    backend: str = "auto",
    wsl_distribution: str = "Ubuntu-Foldseek",
    wsl_executable: str = "/opt/colabfold/env/bin/mmseqs",
) -> str:
    runtime = resolve_mmseqs_runtime(executable, backend, wsl_distribution, wsl_executable)
    result = capture_process_text(
        _command(runtime, ["version"]),
        timeout=20,
        check=True,
    )
    return (result.stdout or result.stderr).strip()


def run_easy_search(
    query_fasta: str | Path,
    target_database: str | Path,
    out_tsv: str | Path,
    tmp_dir: str | Path,
    *,
    threads: int = 16,
    iterations: int = 2,
    max_seqs: int = 2500,
    search_type: int = 1,
    evalue: float = 1e-5,
    fields: Sequence[str] = MMSEQS_FIELDS,
    extra_args: Sequence[str] | None = None,
    executable: str = "mmseqs",
    backend: str = "auto",
    wsl_distribution: str = "Ubuntu-Foldseek",
    wsl_executable: str = "/opt/colabfold/env/bin/mmseqs",
    log_path: str | Path | None = None,
) -> Path:
    query_fasta = Path(query_fasta)
    out_tsv = Path(out_tsv)
    tmp_dir = Path(tmp_dir)
    if not query_fasta.is_file():
        raise FileNotFoundError(query_fasta)
    ensure_dir(out_tsv.parent)
    ensure_dir(tmp_dir)
    runtime = resolve_mmseqs_runtime(executable, backend, wsl_distribution, wsl_executable)
    args = build_easy_search_args(
        runtime_path(query_fasta, runtime),
        runtime_path(target_database, runtime),
        runtime_path(out_tsv, runtime),
        runtime_path(tmp_dir, runtime),
        threads=threads,
        iterations=iterations,
        max_seqs=max_seqs,
        search_type=search_type,
        evalue=evalue,
        fields=fields,
        extra_args=extra_args,
    )
    run_mmseqs_command(runtime, args, log_path=log_path)
    if not out_tsv.is_file():
        raise RuntimeError(f"MMseqs2 did not create its expected output: {out_tsv}")
    return out_tsv


def run_easy_cluster(
    input_fasta: str | Path,
    output_prefix: str | Path,
    tmp_dir: str | Path,
    *,
    min_sequence_identity: float = 0.5,
    threads: int = 16,
    extra_args: Sequence[str] | None = None,
    executable: str = "mmseqs",
    backend: str = "auto",
    wsl_distribution: str = "Ubuntu-Foldseek",
    wsl_executable: str = "/opt/colabfold/env/bin/mmseqs",
    log_path: str | Path | None = None,
) -> tuple[Path, Path]:
    input_fasta = Path(input_fasta)
    output_prefix = Path(output_prefix)
    tmp_dir = Path(tmp_dir)
    if not input_fasta.is_file():
        raise FileNotFoundError(input_fasta)
    ensure_dir(output_prefix.parent)
    ensure_dir(tmp_dir)
    runtime = resolve_mmseqs_runtime(executable, backend, wsl_distribution, wsl_executable)
    args = [
        "easy-cluster",
        runtime_path(input_fasta, runtime),
        runtime_path(output_prefix, runtime),
        runtime_path(tmp_dir, runtime),
        "--min-seq-id",
        str(min_sequence_identity),
        "--threads",
        str(threads),
    ]
    if extra_args:
        args.extend(str(value) for value in extra_args)
    run_mmseqs_command(runtime, args, log_path=log_path)
    representatives = Path(f"{output_prefix}_rep_seq.fasta")
    clusters = Path(f"{output_prefix}_cluster.tsv")
    if not representatives.is_file() or not clusters.is_file():
        raise RuntimeError(
            "MMseqs2 clustering did not create representative FASTA and cluster TSV"
        )
    return representatives, clusters


def download_nr_database(
    target_database: str | Path,
    tmp_dir: str | Path,
    *,
    threads: int = 64,
    executable: str = "mmseqs",
    backend: str = "native",
    wsl_distribution: str = "Ubuntu-Foldseek",
    wsl_executable: str = "/opt/colabfold/env/bin/mmseqs",
    log_path: str | Path | None = None,
) -> Path:
    """Explicitly download and build NCBI NR; never called by a normal search."""
    target_database = Path(target_database)
    tmp_dir = Path(tmp_dir)
    ensure_dir(target_database.parent)
    ensure_dir(tmp_dir)
    runtime = resolve_mmseqs_runtime(executable, backend, wsl_distribution, wsl_executable)
    args = [
        "databases",
        "NR",
        runtime_path(target_database, runtime),
        runtime_path(tmp_dir, runtime),
        "--threads",
        str(threads),
    ]
    run_mmseqs_command(runtime, args, log_path=log_path)
    if not Path(f"{target_database}.dbtype").is_file():
        raise RuntimeError(f"NCBI NR MMseqs2 database was not created: {target_database}")
    return target_database
