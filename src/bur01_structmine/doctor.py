"""Environment diagnostics for reproducible AquaEsterNet-StructMine runs."""
from __future__ import annotations

import importlib.metadata
import json
import os
import platform
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import requests

from .config import load_config, resolve_path
from .colabfold_runtime import ColabFoldRuntime, ColabFoldUnavailable, process_exists
from .datasources import load_search_databases
from .foldseek_runner import FoldseekUnavailable, get_foldseek_version
from .mmseqs_runner import MmseqsUnavailable, get_mmseqs_version
from .run_context import file_sha256
from .structure_search import load_ca_trace

LOCAL_ENGINES = {"tmtools_tmalign", "tmalign", "python_ca", "local_ca"}


@dataclass(frozen=True)
class DoctorCheck:
    name: str
    category: str
    status: str
    message: str
    details: dict[str, Any]


@dataclass(frozen=True)
class DoctorReport:
    project_version: str
    platform: str
    checks: list[DoctorCheck]

    @property
    def passed(self) -> bool:
        return not any(check.status == "fail" for check in self.checks)

    def to_dict(self) -> dict[str, Any]:
        return {
            "project_version": self.project_version,
            "platform": self.platform,
            "passed": self.passed,
            "summary": {
                status: sum(check.status == status for check in self.checks)
                for status in ("pass", "warn", "fail")
            },
            "checks": [asdict(check) for check in self.checks],
        }

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2, ensure_ascii=False)


def _check(name: str, category: str, status: str, message: str, **details: Any) -> DoctorCheck:
    return DoctorCheck(name, category, status, message, details)


def _feature_state(feature_name: str) -> tuple[str, str]:
    if os.name != "nt":
        return "not-applicable", "Windows optional features are not used on this platform"
    try:
        result = subprocess.run(
            ["dism.exe", "/Online", "/English", "/Get-FeatureInfo", f"/FeatureName:{feature_name}"],
            capture_output=True,
            text=True,
            errors="replace",
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return "unknown", f"{type(exc).__name__}: {exc}"
    output = (result.stdout or "") + (result.stderr or "")
    for line in output.splitlines():
        if line.strip().lower().startswith("state"):
            return line.split(":", 1)[-1].strip().lower(), output.strip()
    return "unknown", output.strip()


def _cpu_has_avx2() -> bool | None:
    try:
        import numpy as np

        features = getattr(np._core._multiarray_umath, "__cpu_features__", {})
        if "AVX2" in features:
            return bool(features["AVX2"])
    except (AttributeError, ImportError):
        pass
    return None


def _drive_checks(paths: list[Path]) -> list[DoctorCheck]:
    checks: list[DoctorCheck] = []
    roots = {Path(path.anchor) for path in paths if path.anchor}
    for root in sorted(roots, key=str):
        root_label = str(root).rstrip("\\/")
        try:
            usage = shutil.disk_usage(root)
            free_gb = usage.free / (1024**3)
            status = "pass" if free_gb >= 50 else "warn"
            checks.append(
                _check(
                    f"disk_{root_label}",
                    "storage",
                    status,
                    f"{free_gb:.1f} GB free on {root}",
                    root=str(root),
                    free_bytes=usage.free,
                    total_bytes=usage.total,
                )
            )
        except OSError as exc:
            checks.append(_check(f"disk_{root}", "storage", "fail", str(exc), root=str(root)))
    return checks


def _database_exists(path: Path) -> bool:
    return path.exists() or Path(f"{path}.dbtype").exists() or Path(f"{path}.index").exists()


def _nr_snapshot_check() -> DoctorCheck:
    pointer_path = resolve_path("databases/nr_demo/latest.json")
    if not pointer_path.is_file():
        return _check(
            "nr_demo_snapshot",
            "nr",
            "warn",
            "No local NR hit snapshot is installed; run nr-demo-build with explicit NCBI consent",
            path=str(pointer_path),
            database_scope="not_installed",
        )
    try:
        pointer = json.loads(pointer_path.read_text(encoding="utf-8"))
        fasta_path = Path(pointer["fasta_path"])
        manifest_path = Path(pointer["manifest_path"])
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        checksum_ok = (
            fasta_path.is_file()
            and file_sha256(fasta_path) == manifest.get("fasta_sha256")
        )
        scope_ok = manifest.get("database_scope") == "nr_remote_snapshot"
        return _check(
            "nr_demo_snapshot",
            "nr",
            "pass" if checksum_ok and scope_ok else "fail",
            (
                f"Local NCBI NR hit snapshot: {manifest.get('sequence_count', 0)} sequences"
                if checksum_ok and scope_ok
                else "Local NR snapshot provenance or checksum is invalid"
            ),
            path=str(fasta_path),
            database_scope=manifest.get("database_scope", ""),
            created_at=manifest.get("created_at", ""),
            sequence_count=manifest.get("sequence_count", 0),
            checksum_ok=checksum_ok,
        )
    except (OSError, KeyError, json.JSONDecodeError) as exc:
        return _check(
            "nr_demo_snapshot",
            "nr",
            "fail",
            f"{type(exc).__name__}: {exc}",
            path=str(pointer_path),
        )


def _nr_job_bundle_check() -> DoctorCheck:
    run_root = resolve_path("results/nr_runs")
    bundles = list(run_root.glob("*/job_bundle/BUNDLE_SHA256SUMS"))
    if not bundles:
        return _check(
            "full_nr_job_bundles",
            "nr",
            "warn",
            "No complete-NR server job bundle has been exported",
            run_root=str(run_root),
            minimum_server_ram_gb=128,
            recommended_free_disk_gb=500,
        )
    invalid: list[str] = []
    for checksum_path in bundles:
        bundle_dir = checksum_path.parent
        try:
            for line in checksum_path.read_text(encoding="ascii").splitlines():
                if not line.strip():
                    continue
                checksum, relative = line.split(maxsplit=1)
                target = bundle_dir / relative.strip()
                if not target.is_file() or file_sha256(target) != checksum:
                    invalid.append(str(target))
        except (OSError, ValueError):
            invalid.append(str(checksum_path))
    return _check(
        "full_nr_job_bundles",
        "nr",
        "pass" if not invalid else "fail",
        (
            f"{len(bundles)} complete-NR server job bundle(s) passed local integrity checks"
            if not invalid
            else f"{len(invalid)} full-NR bundle file(s) failed integrity checks"
        ),
        bundle_count=len(bundles),
        invalid=invalid[:20],
        minimum_server_ram_gb=128,
        recommended_free_disk_gb=500,
    )


def run_doctor(config_path: str | Path = "config.yaml", check_network: bool = True) -> DoctorReport:
    checks: list[DoctorCheck] = []
    try:
        cfg = load_config(config_path)
        checks.append(_check("configuration", "project", "pass", "Configuration is valid"))
    except Exception as exc:  # noqa: BLE001 - doctor must aggregate all configuration failures.
        checks.append(_check("configuration", "project", "fail", f"{type(exc).__name__}: {exc}"))
        return DoctorReport("0.4.0", platform.platform(), checks)
    configured_sources = load_search_databases(cfg)
    needs_foldseek = any(source.search_engine not in LOCAL_ENGINES for source in configured_sources)

    python_ok = tuple(int(part) for part in platform.python_version_tuple()[:2]) >= (3, 10)
    checks.append(
        _check(
            "python",
            "runtime",
            "pass" if python_ok else "fail",
            f"Python {platform.python_version()}",
            executable=os.sys.executable,
        )
    )

    required_packages = {
        "pydantic": "2.13.4",
        "tmtools": "0.3.0",
        "streamlit": "1.59.1",
        "plotly": "6.9.0",
        "py3Dmol": "2.5.5",
        "pandas": "2.0.0",
        "biopython": "1.83",
    }
    for package, minimum in required_packages.items():
        try:
            version = importlib.metadata.version(package)
            status = "pass"
            message = f"{package} {version}"
        except importlib.metadata.PackageNotFoundError:
            version = ""
            status = "fail"
            message = f"{package} is not installed"
        checks.append(_check(f"package_{package}", "runtime", status, message, minimum=minimum, version=version))

    avx2 = _cpu_has_avx2()
    checks.append(
        _check(
            "cpu_avx2",
            "hardware",
            "pass" if avx2 else ("warn" if avx2 is None else "fail"),
            "AVX2 is available" if avx2 else ("AVX2 detection was inconclusive" if avx2 is None else "AVX2 is unavailable"),
        )
    )

    for feature in ("Microsoft-Windows-Subsystem-Linux", "VirtualMachinePlatform"):
        state, output = _feature_state(feature)
        enabled = state == "enabled"
        checks.append(
            _check(
                f"windows_feature_{feature}",
                "wsl",
                "pass" if enabled else ("fail" if needs_foldseek else "warn"),
                f"{feature}: {state}",
                state=state,
                diagnostic=output[-2000:],
            )
        )

    wsl_path = shutil.which("wsl.exe")
    checks.append(
        _check(
            "wsl_executable",
            "wsl",
            "pass" if wsl_path else "fail",
            f"WSL executable: {wsl_path}" if wsl_path else "wsl.exe was not found",
        )
    )
    configured_distro = str(cfg.get("foldseek", {}).get("wsl_distribution", "Ubuntu-Foldseek"))
    distributions: list[str] = []
    if wsl_path:
        try:
            result = subprocess.run(
                ["wsl.exe", "--list", "--quiet"],
                capture_output=True,
                text=True,
                errors="replace",
                timeout=30,
                check=False,
            )
            distributions = [line.replace("\x00", "").strip() for line in result.stdout.splitlines() if line.replace("\x00", "").strip()]
        except (OSError, subprocess.SubprocessError):
            distributions = []
    checks.append(
        _check(
            "wsl_distribution",
            "wsl",
            "pass" if configured_distro in distributions else ("fail" if needs_foldseek else "warn"),
            f"{configured_distro}: {'installed' if configured_distro in distributions else 'not installed'}",
            installed=distributions,
        )
    )

    colabfold_runtime = ColabFoldRuntime(
        backend="wsl",
        executable="/opt/colabfold/bin/colabfold_batch",
        wsl_distribution=configured_distro,
    )
    if configured_distro in distributions:
        try:
            mmseqs_version = get_mmseqs_version(
                backend="wsl",
                wsl_distribution=configured_distro,
                wsl_executable="/opt/colabfold/env/bin/mmseqs",
            )
            checks.append(
                _check(
                    "mmseqs2",
                    "nr",
                    "pass" if mmseqs_version == "18.8cc5c" else "warn",
                    f"MMseqs2 {mmseqs_version}",
                    expected_version="18.8cc5c",
                    executable="/opt/colabfold/env/bin/mmseqs",
                )
            )
        except (MmseqsUnavailable, OSError, subprocess.SubprocessError) as exc:
            checks.append(
                _check("mmseqs2", "nr", "fail", f"{type(exc).__name__}: {exc}")
            )

        try:
            blast_result = subprocess.run(
                [
                    "wsl.exe",
                    "-d",
                    configured_distro,
                    "--",
                    "/opt/ncbi-blast/bin/blastp",
                    "-version",
                ],
                capture_output=True,
                text=True,
                errors="replace",
                timeout=30,
                check=False,
            )
            blast_text = (blast_result.stdout or blast_result.stderr).strip()
            checks.append(
                _check(
                    "ncbi_blast_plus",
                    "nr",
                    "pass"
                    if blast_result.returncode == 0 and "2.17.0" in blast_text
                    else "warn",
                    blast_text or "NCBI BLAST+ is not installed at the locked path",
                    expected_version="2.17.0",
                    executable="/opt/ncbi-blast/bin/blastp",
                )
            )
        except (OSError, subprocess.SubprocessError) as exc:
            checks.append(
                _check(
                    "ncbi_blast_plus",
                    "nr",
                    "warn",
                    f"{type(exc).__name__}: {exc}",
                )
            )

        try:
            gpu_result = subprocess.run(
                [
                    "wsl.exe",
                    "-d",
                    configured_distro,
                    "--",
                    "nvidia-smi",
                    "--query-gpu=name,memory.total,driver_version,compute_cap",
                    "--format=csv,noheader",
                ],
                capture_output=True,
                text=True,
                errors="replace",
                timeout=30,
                check=False,
            )
            gpu_text = gpu_result.stdout.strip()
            checks.append(
                _check(
                    "wsl_gpu",
                    "prediction",
                    "pass" if gpu_result.returncode == 0 and gpu_text else "fail",
                    gpu_text or (gpu_result.stderr.strip() or "WSL GPU was not detected"),
                )
            )
        except (OSError, subprocess.SubprocessError) as exc:
            checks.append(_check("wsl_gpu", "prediction", "fail", str(exc)))

        try:
            versions = colabfold_runtime.versions()
            version_ok = versions.get("colabfold") == "1.6.1"
            checks.append(
                _check(
                    "colabfold",
                    "prediction",
                    "pass" if version_ok else "fail",
                    ", ".join(f"{name} {version}" for name, version in versions.items()),
                    versions=versions,
                )
            )
            patches = colabfold_runtime.patch_provenance()
            checks.append(
                _check(
                    "colabfold_network_patch",
                    "prediction",
                    "pass" if patches and patches[0].get("patched_sha256") else "warn",
                    "Recorded result-download timeout patch" if patches else "ColabFold network patch provenance is missing",
                    patches=patches,
                )
            )
            gpu_details = colabfold_runtime.gpu_details()
            checks.append(
                _check(
                    "jax_cuda",
                    "prediction",
                    "pass",
                    "; ".join(gpu_details.get("devices", [])),
                    **gpu_details,
                )
            )
        except (ColabFoldUnavailable, OSError, subprocess.SubprocessError, json.JSONDecodeError) as exc:
            checks.append(_check("colabfold", "prediction", "warn", f"{type(exc).__name__}: {exc}"))

        parameter_result = subprocess.run(
            [
                "wsl.exe",
                "-d",
                configured_distro,
                "--",
                "bash",
                "-lc",
                "find /opt/colabfold/data -maxdepth 3 -type f -name 'params_model_*.npz' -print -quit",
            ],
            capture_output=True,
            text=True,
            errors="replace",
            timeout=30,
            check=False,
        )
        parameter_path = parameter_result.stdout.strip() if parameter_result.returncode == 0 else ""
        checks.append(
            _check(
                "alphafold_parameters",
                "prediction",
                "pass" if parameter_path else "warn",
                parameter_path or "AlphaFold2 parameters are not cached yet; ColabFold downloads them on the first prediction",
            )
        )
        provenance_result = subprocess.run(
            [
                "wsl.exe",
                "-d",
                configured_distro,
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
        try:
            parameter_provenance = json.loads(provenance_result.stdout)
        except json.JSONDecodeError:
            parameter_provenance = {}
        provenance_ok = bool(parameter_provenance.get("archive_sha256")) and len(
            parameter_provenance.get("parameter_files", {})
        ) == 5
        checks.append(
            _check(
                "alphafold_parameter_provenance",
                "prediction",
                "pass" if provenance_ok else "warn",
                "AlphaFold2-PTM archive and five parameter checksums recorded"
                if provenance_ok
                else "AlphaFold2 parameter provenance is not complete yet",
                provenance=parameter_provenance,
            )
        )

    fs_cfg = cfg.get("foldseek", {})
    try:
        version = get_foldseek_version(
            str(fs_cfg.get("executable", "foldseek")),
            backend=str(fs_cfg.get("backend", "auto")),
            wsl_distribution=str(fs_cfg.get("wsl_distribution", "Ubuntu-Foldseek")),
            wsl_executable=str(fs_cfg.get("wsl_executable", "/opt/foldseek/bin/foldseek")),
        )
        checks.append(_check("foldseek", "runtime", "pass", f"Foldseek {version}", version=version))
    except (FoldseekUnavailable, OSError, subprocess.SubprocessError) as exc:
        checks.append(_check("foldseek", "runtime", "fail" if needs_foldseek else "warn", str(exc)))

    query_cfg = cfg["query"]
    query_path = resolve_path(query_cfg["pdb_file"])
    try:
        trace = load_ca_trace(query_path, chain_id=query_cfg.get("chain"))
        checks.append(
            _check(
                "query_structure",
                "project",
                "pass",
                f"{query_cfg.get('pdb_id') or query_path.stem} chain {trace.chain_id}: {len(trace.residues)} residues",
                path=str(query_path),
                sha256=file_sha256(query_path),
                pdb_id=query_cfg.get("pdb_id", ""),
                chain=trace.chain_id,
                residue_range=query_cfg.get("residue_range", ""),
            )
        )
    except Exception as exc:  # noqa: BLE001
        checks.append(_check("query_structure", "project", "fail", f"{type(exc).__name__}: {exc}", path=str(query_path)))

    prediction_root = resolve_path("results/predictions")
    nr_run_root = resolve_path("results/nr_runs")
    data_paths = [
        resolve_path("databases"),
        resolve_path(cfg["output"].get("run_root", "results/runs")),
        prediction_root,
        nr_run_root,
    ]
    if os.name == "nt" and Path("D:/").exists():
        data_paths.append(Path("D:/Apps"))
    for source in configured_sources:
        path = source.structure_dir if source.search_engine in LOCAL_ENGINES else source.target_db
        if path is None:
            checks.append(_check(f"database_{source.safe_name}", "database", "fail", f"{source.name} has no usable path"))
            continue
        data_paths.append(path)
        exists = path.is_dir() if source.search_engine in LOCAL_ENGINES else _database_exists(path)
        checks.append(
            _check(
                f"database_{source.safe_name}",
                "database",
                "pass" if exists else "fail",
                f"{source.name}: {'ready' if exists else 'missing'} at {path}",
                path=str(path),
                search_engine=source.search_engine,
            )
        )

    output_root = resolve_path(cfg["output"].get("run_root", "results/runs"))
    try:
        output_root.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(prefix="doctor-", dir=output_root, delete=True):
            pass
        checks.append(_check("output_permissions", "storage", "pass", f"Output directory is writable: {output_root}"))
    except OSError as exc:
        checks.append(_check("output_permissions", "storage", "fail", str(exc), path=str(output_root)))

    checks.extend(_drive_checks(data_paths))
    checks.append(_nr_snapshot_check())
    checks.append(_nr_job_bundle_check())
    ncbi_email = os.environ.get("AQUA_NCBI_EMAIL", "").strip()
    checks.append(
        _check(
            "ncbi_contact_email",
            "nr",
            "pass" if "@" in ncbi_email else "warn",
            (
                "AQUA_NCBI_EMAIL is configured"
                if "@" in ncbi_email
                else "AQUA_NCBI_EMAIL is required before a real remote NR request"
            ),
            configured=bool(ncbi_email),
        )
    )

    gpu_lock = prediction_root / ".gpu.lock"
    if gpu_lock.exists():
        try:
            owner = json.loads(gpu_lock.read_text(encoding="utf-8"))
            owner_pid = int(owner.get("pid", -1))
            if not process_exists(owner_pid):
                raise OSError("stale lock owner")
            checks.append(
                _check(
                    "gpu_lock",
                    "prediction",
                    "warn",
                    f"GPU is reserved by process {owner_pid}",
                    owner=owner,
                )
            )
        except (OSError, ValueError, json.JSONDecodeError):
            checks.append(
                _check(
                    "gpu_lock",
                    "prediction",
                    "warn",
                    f"Stale GPU lock can be recovered automatically: {gpu_lock}",
                )
            )
    else:
        checks.append(_check("gpu_lock", "prediction", "pass", "GPU prediction lock is available"))

    if check_network:
        endpoints = {
            "uniprot": "https://rest.uniprot.org/uniprotkb/P00520.json",
            "rcsb": "https://data.rcsb.org/rest/v1/core/entry/8HFW",
            "alphafold": "https://alphafold.ebi.ac.uk/api/prediction/Q6RYA0",
            "colabfold_msa": "https://api.colabfold.com",
            "ncbi_blast_policy": "https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html",
        }
        session = requests.Session()
        session.headers.update({"User-Agent": "AquaEsterNet-StructMine/0.4.0 doctor"})
        for name, url in endpoints.items():
            try:
                response = session.get(url, timeout=15)
                if response.status_code >= 500:
                    response.raise_for_status()
                checks.append(_check(f"network_{name}", "network", "pass", f"{name} API reachable", status_code=response.status_code))
            except requests.RequestException as exc:
                checks.append(_check(f"network_{name}", "network", "warn", f"{type(exc).__name__}: {exc}", url=url))

    return DoctorReport("0.4.0", platform.platform(), checks)


def print_doctor(report: DoctorReport) -> None:
    labels = {"pass": "OK", "warn": "WARN", "fail": "FAIL"}
    for check in report.checks:
        print(f"[{labels[check.status]}] {check.category}/{check.name}: {check.message}")
    summary = report.to_dict()["summary"]
    print(f"\nDoctor summary: {summary['pass']} passed, {summary['warn']} warnings, {summary['fail']} failed")
