"""Capture and validate the NR workbench with an isolated, sanitized fixture."""
from __future__ import annotations

import argparse
import json
import os
import shutil
import socket
import subprocess
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
from PIL import Image, ImageStat
from playwright.sync_api import Browser, Page, sync_playwright

from bur01_structmine.config import project_root
from bur01_structmine.models import NrRunManifest, ProgressEvent, RunManifest
from bur01_structmine.run_context import file_sha256

EDGE = Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe")
FIXTURE_RUN_ID = "ui-validation-8hfw-structure"
FIXTURE_NR_RUN_ID = "ui-validation-nr-sequence-screen"
DESKTOP_VIEWPORT = {"width": 1440, "height": 900}
MOBILE_VIEWPORT = {"width": 390, "height": 844}


def _free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as handle:
        handle.bind(("127.0.0.1", 0))
        return int(handle.getsockname()[1])


def _wait_for_server(url: str, process: subprocess.Popen[str]) -> None:
    deadline = time.monotonic() + 90
    health_url = f"{url}/_stcore/health"
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise RuntimeError(f"Streamlit exited early with code {process.returncode}")
        try:
            with urllib.request.urlopen(health_url, timeout=2) as response:
                if response.status == 200:
                    return
        except OSError:
            time.sleep(0.4)
    raise TimeoutError(f"Streamlit did not become ready: {health_url}")


def _safe_remove_fixture(path: Path, root: Path) -> None:
    resolved = path.resolve()
    allowed_parent = (root / "results" / "validation").resolve()
    if resolved.parent != allowed_parent or not resolved.name.startswith("ui_nr_fixture_"):
        raise ValueError(f"Refusing to remove unexpected validation path: {resolved}")
    if resolved.exists():
        shutil.rmtree(resolved)


def _mutate(sequence: str, replacements: dict[int, str]) -> str:
    residues = list(sequence)
    for one_based_position, residue in replacements.items():
        residues[one_based_position - 1] = residue
    return "".join(residues)


def _prepare_fixture(root: Path, fixture_root: Path) -> dict[str, Path]:
    run_root = fixture_root / "runs"
    nr_root = fixture_root / "nr_runs"
    prediction_root = fixture_root / "predictions"
    run_dir = run_root / FIXTURE_RUN_ID
    nr_dir = nr_root / FIXTURE_NR_RUN_ID
    run_dir.mkdir(parents=True, exist_ok=True)
    nr_dir.mkdir(parents=True, exist_ok=True)
    prediction_root.mkdir(parents=True, exist_ok=True)

    query_fasta = root / "data" / "query" / "Bur01m_seed.fasta"
    query_pdb = root / "data" / "query" / "Bur01m.pdb"
    fasta_lines = query_fasta.read_text(encoding="ascii").splitlines()
    query_sequence = "".join(line.strip() for line in fasta_lines if not line.startswith(">"))
    fixture_query = run_dir / "query.fasta"
    shutil.copy2(query_fasta, fixture_query)

    run_manifest = RunManifest(
        run_id=FIXTURE_RUN_ID,
        status="completed",
        config_path=str((root / "configs" / "benchmark_foldseek.yaml").resolve()),
        query_pdb_id="8HFW",
        query_chain="A",
        query_structure_path=str(query_pdb.resolve()),
        artifacts={"query_fasta": str(fixture_query.resolve())},
        child_nr_run_ids=[FIXTURE_NR_RUN_ID],
    )
    (run_dir / "run_manifest.json").write_text(
        run_manifest.model_dump_json(indent=2),
        encoding="utf-8",
    )

    rows = [
        {
            "nr_rank": 1,
            "target_accession": "8HFW_A",
            "organism": "Burkholderia sp. (self control)",
            "evalue": 0.0,
            "bits": 612.0,
            "fident": 1.0,
            "qcov": 1.0,
            "tcov": 1.0,
            "cluster_size": 7,
            "motif_gxsxg_like": True,
            "mapped_catalytic_identity_count": 3,
            "sequence_screen_status": "pass_sequence_prescreen",
            "af2_ready": True,
            "wase_tier": "Not evaluated - structure required",
            "wase_scope": "not_evaluated_sequence_only",
            "target_sequence": query_sequence,
        },
        {
            "nr_rank": 2,
            "target_accession": "Q9FVW3",
            "organism": "Arabidopsis thaliana",
            "evalue": 2.1e-74,
            "bits": 286.0,
            "fident": 0.43,
            "qcov": 0.93,
            "tcov": 0.88,
            "cluster_size": 12,
            "motif_gxsxg_like": True,
            "mapped_catalytic_identity_count": 3,
            "sequence_screen_status": "pass_sequence_prescreen",
            "af2_ready": True,
            "wase_tier": "Not evaluated - structure required",
            "wase_scope": "not_evaluated_sequence_only",
            "target_sequence": _mutate(query_sequence, {16: "A", 58: "V", 196: "L"}),
        },
        {
            "nr_rank": 3,
            "target_accession": "Q6RYA0",
            "organism": "Oryza sativa",
            "evalue": 7.6e-61,
            "bits": 247.0,
            "fident": 0.39,
            "qcov": 0.90,
            "tcov": 0.86,
            "cluster_size": 9,
            "motif_gxsxg_like": True,
            "mapped_catalytic_identity_count": 3,
            "sequence_screen_status": "pass_sequence_prescreen",
            "af2_ready": True,
            "wase_tier": "Not evaluated - structure required",
            "wase_scope": "not_evaluated_sequence_only",
            "target_sequence": _mutate(query_sequence, {22: "G", 83: "A", 220: "V"}),
        },
        {
            "nr_rank": 4,
            "target_accession": "DEMO_NR_004",
            "organism": "Uncultured bacterium",
            "evalue": 1.4e-38,
            "bits": 171.0,
            "fident": 0.31,
            "qcov": 0.84,
            "tcov": 0.78,
            "cluster_size": 4,
            "motif_gxsxg_like": True,
            "mapped_catalytic_identity_count": 2,
            "sequence_screen_status": "hold_catalytic_triad_identity",
            "af2_ready": False,
            "wase_tier": "Not evaluated - structure required",
            "wase_scope": "not_evaluated_sequence_only",
            "target_sequence": _mutate(query_sequence, {111: "A", 142: "V"}),
        },
        {
            "nr_rank": 5,
            "target_accession": "P00520",
            "organism": "Rous sarcoma virus (negative control)",
            "evalue": 0.018,
            "bits": 34.0,
            "fident": 0.16,
            "qcov": 0.31,
            "tcov": 0.28,
            "cluster_size": 1,
            "motif_gxsxg_like": False,
            "mapped_catalytic_identity_count": 0,
            "sequence_screen_status": "reject_default_sequence_threshold",
            "af2_ready": False,
            "wase_tier": "Not evaluated - structure required",
            "wase_scope": "not_evaluated_sequence_only",
            "target_sequence": "A" * len(query_sequence),
        },
    ]
    candidates = pd.DataFrame(rows)
    candidates_tsv = nr_dir / "nr_candidates.tsv"
    candidates_xlsx = nr_dir / "nr_candidates.xlsx"
    candidates_fasta = nr_dir / "af2_ready.fasta"
    candidates.to_csv(candidates_tsv, sep="\t", index=False)
    candidates.to_excel(candidates_xlsx, index=False)
    candidates_fasta.write_text(
        "".join(
            f">{row['target_accession']}\n{row['target_sequence']}\n"
            for row in rows
            if row["af2_ready"]
        ),
        encoding="ascii",
    )

    timestamp = datetime.now(timezone.utc)
    event_specs = [
        ("seed_export", "completed", "Exported 8HFW chain A and structure-derived seeds", 0.12),
        ("mmseqs_search", "completed", "Completed MMseqs2 iterative search (2 rounds)", 0.36),
        ("deduplicate", "completed", "Merged accessions and identical sequence hashes", 0.52),
        ("cluster", "completed", "Clustered hits at 50% sequence identity", 0.68),
        ("catalytic_mapping", "completed", "Mapped Ser111/Asp241/His274 to representatives", 0.82),
        ("sequence_prescreen", "completed", "Applied transparent sequence-only evidence filters", 0.92),
        ("af2_queue", "completed", "Prepared AF2-ready queue; WASE remains unevaluated", 1.0),
    ]
    events_path = nr_dir / "events.jsonl"
    events_path.write_text(
        "\n".join(
            ProgressEvent(
                timestamp=timestamp,
                stage=stage,
                state=state,
                message=message,
                progress=progress,
            ).model_dump_json()
            for stage, state, message, progress in event_specs
        )
        + "\n",
        encoding="utf-8",
    )

    nr_manifest = NrRunManifest(
        nr_run_id=FIXTURE_NR_RUN_ID,
        status="completed",
        started_at=timestamp,
        completed_at=timestamp,
        parent_structure_run_id=FIXTURE_RUN_ID,
        execution_mode="local",
        database_scope="custom_mmseqs",
        database_name="Sanitized deterministic UI validation fixture",
        database_snapshot={
            "scope": "custom_mmseqs",
            "snapshot_date": timestamp.date().isoformat(),
            "is_real_ncbi_result": False,
        },
        config_path=str((root / "configs" / "nr_local_demo.yaml").resolve()),
        paper_preset="hong2025",
        resolved_parameters={
            "iterations": 2,
            "max_seqs": 2500,
            "search_type": 1,
            "cluster_min_sequence_identity": 0.5,
        },
        software_versions={"mmseqs": "18.8cc5c", "blastp": "2.17.0+"},
        seed_count=1,
        raw_hit_count=48,
        deduplicated_hit_count=31,
        representative_count=len(rows),
        candidate_count=len(rows),
        af2_ready_count=sum(bool(row["af2_ready"]) for row in rows),
        artifacts={
            "candidates_tsv": str(candidates_tsv.resolve()),
            "candidates_xlsx": str(candidates_xlsx.resolve()),
            "candidates_fasta": str(candidates_fasta.resolve()),
            "events": str(events_path.resolve()),
        },
        warnings=[
            "界面验证数据，不是真实 NCBI NR 命中；真实查询必须逐次确认序列传输。"
        ],
    )
    (nr_dir / "nr_run_manifest.json").write_text(
        nr_manifest.model_dump_json(indent=2),
        encoding="utf-8",
    )
    return {
        "run_root": run_root,
        "nr_root": nr_root,
        "prediction_root": prediction_root,
        "manifest": nr_dir / "nr_run_manifest.json",
    }


def _page_variance(path: Path) -> float:
    image = Image.open(path).convert("RGB")
    return float(sum(ImageStat.Stat(image).var) / 3.0)


def _assert_layout(page: Page) -> dict[str, Any]:
    dimensions = page.evaluate(
        """() => ({
          innerWidth: window.innerWidth,
          innerHeight: window.innerHeight,
          scrollWidth: document.documentElement.scrollWidth,
          bodyScrollWidth: document.body.scrollWidth
        })"""
    )
    if dimensions["scrollWidth"] > dimensions["innerWidth"] + 2:
        raise AssertionError(f"Horizontal overflow detected: {dimensions}")
    page.keyboard.press("Tab")
    focus = page.evaluate(
        """() => {
          const el = document.activeElement;
          const r = el?.getBoundingClientRect();
          return {tag: el?.tagName || '', role: el?.getAttribute('role') || '',
                  visible: !!r && r.width > 0 && r.height > 0};
        }"""
    )
    if focus["tag"] in {"", "BODY"} or not focus["visible"]:
        raise AssertionError(f"Keyboard focus is not visible: {focus}")
    return {"dimensions": dimensions, "keyboard_focus": focus}


def _open_nr_page(browser: Browser, url: str, viewport: dict[str, int]) -> tuple[Any, Page]:
    context = browser.new_context(viewport=viewport, device_scale_factor=2)
    page = context.new_page()
    page.goto(url, wait_until="domcontentloaded", timeout=60_000)
    page.locator('[data-testid="stAppViewContainer"]').wait_for(timeout=60_000)
    page.get_by_role("tab", name="NR挖掘").click()
    page.get_by_text("NCBI NR 序列扩增", exact=True).wait_for()
    page.get_by_text("NR 处理轨迹", exact=True).wait_for()
    page.get_by_text("序列预筛证据", exact=True).wait_for()
    page.wait_for_timeout(1_500)
    return context, page


def _scroll_to(page: Page, locator: Any, block: str = "start") -> None:
    locator.wait_for()
    locator.evaluate(
        f"element => element.scrollIntoView({{block: '{block}', inline: 'nearest'}})"
    )
    page.wait_for_timeout(700)


def _capture(page: Page, path: Path, expected_size: tuple[int, int]) -> dict[str, Any]:
    page.screenshot(path=str(path), full_page=False)
    with Image.open(path) as image:
        size = image.size
    if size != expected_size:
        raise AssertionError(f"Unexpected screenshot dimensions for {path.name}: {size}")
    variance = _page_variance(path)
    if variance < 8:
        raise AssertionError(f"Screenshot appears blank: {path} ({variance:.3f})")
    return {
        "path": str(path.relative_to(project_root())).replace("\\", "/"),
        "pixels": {"width": size[0], "height": size[1]},
        "pixel_variance": round(variance, 3),
        "sha256": file_sha256(path),
    }


def _capture_desktop(
    browser: Browser,
    url: str,
    screenshot_dir: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    context, page = _open_nr_page(browser, url, DESKTOP_VIEWPORT)
    page_errors: list[str] = []
    page.on("pageerror", lambda error: page_errors.append(str(error)))
    expected = (DESKTOP_VIEWPORT["width"] * 2, DESKTOP_VIEWPORT["height"] * 2)
    captures: list[dict[str, Any]] = []

    heading = page.get_by_text("NCBI NR 序列扩增", exact=True)
    _scroll_to(page, heading)
    captures.append(
        _capture(
            page,
            screenshot_dir / "nr_workflow_01_workbench.png",
            expected,
        )
    )

    trace = page.get_by_text("NR 处理轨迹", exact=True)
    _scroll_to(page, trace, block="center")
    captures.append(
        _capture(
            page,
            screenshot_dir / "nr_workflow_02_processing_trace.png",
            expected,
        )
    )

    evidence = page.get_by_text("序列预筛证据", exact=True)
    _scroll_to(page, evidence, block="center")
    table = page.locator('[data-testid="stDataFrame"]:visible').last
    table.wait_for()
    if table.bounding_box() is None:
        raise AssertionError("NR candidate table is not visible")
    captures.append(
        _capture(
            page,
            screenshot_dir / "nr_workflow_03_candidates.png",
            expected,
        )
    )

    _scroll_to(page, heading)
    page.get_by_text("完整 NR 服务器作业", exact=True).click()
    page.get_by_text(
        "完整 NR 仅在服务器/HPC 上运行；导回且通过校验前不会发布 v0.4.0。",
        exact=True,
    ).wait_for()
    captures.append(
        _capture(
            page,
            screenshot_dir / "nr_workflow_04_server_job.png",
            expected,
        )
    )
    layout = _assert_layout(page)
    if page_errors:
        raise AssertionError("Browser page errors: " + " | ".join(page_errors))
    context.close()
    return captures, layout


def _capture_mobile(
    browser: Browser,
    url: str,
    screenshot_dir: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    context, page = _open_nr_page(browser, url, MOBILE_VIEWPORT)
    page_errors: list[str] = []
    page.on("pageerror", lambda error: page_errors.append(str(error)))
    heading = page.get_by_text("NCBI NR 序列扩增", exact=True)
    _scroll_to(page, heading)
    page.wait_for_timeout(1_000)
    capture = _capture(
        page,
        screenshot_dir / "nr_workflow_05_mobile.png",
        (MOBILE_VIEWPORT["width"] * 2, MOBILE_VIEWPORT["height"] * 2),
    )
    layout = _assert_layout(page)
    if page_errors:
        raise AssertionError("Browser page errors: " + " | ".join(page_errors))
    context.close()
    return capture, layout


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--edge", default=str(EDGE))
    parser.add_argument(
        "--output",
        default="docs/validation/v0.4_nr_ui_validation.json",
    )
    args = parser.parse_args()

    root = project_root()
    edge = Path(args.edge)
    if not edge.is_file():
        raise FileNotFoundError(f"Microsoft Edge was not found: {edge}")
    fixture_root = (
        root
        / "results"
        / "validation"
        / f"ui_nr_fixture_{os.getpid()}_{int(time.time())}"
    )
    screenshot_dir = root / "docs" / "screenshots"
    screenshot_dir.mkdir(parents=True, exist_ok=True)
    fixture = _prepare_fixture(root, fixture_root)

    port = _free_port()
    url = f"http://127.0.0.1:{port}"
    env = os.environ.copy()
    env.update(
        {
            "PYTHONPATH": str(root / "src"),
            "AQUA_TEST_MODE": "1",
            "AQUA_RUN_ROOT": str(fixture["run_root"]),
            "AQUA_NR_RUN_ROOT": str(fixture["nr_root"]),
            "AQUA_PREDICTION_ROOT": str(fixture["prediction_root"]),
            "AQUA_REMOTE_ACCESS_PATH": str(fixture_root / "no_remote_access.json"),
            "STREAMLIT_BROWSER_GATHER_USAGE_STATS": "false",
        }
    )
    log_path = fixture_root / "streamlit.log"
    log_handle = log_path.open("w", encoding="utf-8")
    creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    process = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "streamlit",
            "run",
            str(root / "ui" / "app.py"),
            "--server.headless=true",
            f"--server.port={port}",
            "--server.address=127.0.0.1",
            "--browser.gatherUsageStats=false",
        ],
        cwd=root,
        env=env,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        text=True,
        creationflags=creationflags,
    )
    service_stopped = False
    report: dict[str, Any]
    try:
        _wait_for_server(url, process)
        with sync_playwright() as playwright:
            browser = playwright.chromium.launch(
                executable_path=str(edge),
                headless=True,
                args=[
                    "--enable-webgl",
                    "--ignore-gpu-blocklist",
                    "--use-angle=swiftshader",
                ],
            )
            try:
                desktop, desktop_layout = _capture_desktop(
                    browser,
                    url,
                    screenshot_dir,
                )
                mobile, mobile_layout = _capture_mobile(
                    browser,
                    url,
                    screenshot_dir,
                )
            finally:
                browser.close()
        report = {
            "project_version": "0.4.0",
            "passed": True,
            "captured_at": datetime.now(timezone.utc).isoformat(),
            "data_scope": "sanitized deterministic UI validation fixture",
            "real_ncbi_result": False,
            "scientific_guardrail": (
                "Sequence-only NR evidence; WASE-screen remains unevaluated "
                "until AF2/Foldseek/TM-align validation."
            ),
            "viewports": {
                "desktop_css": DESKTOP_VIEWPORT,
                "desktop_output_pixels": {"width": 2880, "height": 1800},
                "mobile_css": MOBILE_VIEWPORT,
                "mobile_output_pixels": {"width": 780, "height": 1688},
            },
            "screenshots": [*desktop, mobile],
            "layout_checks": {
                "desktop": desktop_layout,
                "mobile": mobile_layout,
            },
            "temporary_service": {
                "binding": "localhost ephemeral port",
                "hidden_window": bool(creationflags),
                "devtunnel_started": False,
            },
        }
    finally:
        if process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=15)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=5)
        log_handle.close()
        service_stopped = process.poll() is not None
        _safe_remove_fixture(fixture_root, root)

    report["temporary_service"]["stopped_after_capture"] = service_stopped
    output = Path(args.output)
    if not output.is_absolute():
        output = root / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
