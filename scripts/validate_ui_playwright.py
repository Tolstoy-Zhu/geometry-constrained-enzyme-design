"""Exercise the AF prediction workbench at desktop and mobile viewports."""
from __future__ import annotations

import argparse
import base64
import io
import json
import math
import re
from pathlib import Path
from typing import Any

from PIL import Image, ImageStat
from playwright.sync_api import Page, sync_playwright

from bur01_structmine.config import project_root

EDGE = Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe")
PREDICTION_ID = "20260711_AF2_8HFW_standard_final_v3"
VIEWPORTS = {
    "desktop": {"width": 1440, "height": 900},
    "mobile": {"width": 390, "height": 844},
}


def _pixel_variance(png: bytes) -> float:
    image = Image.open(io.BytesIO(png)).convert("RGB")
    return float(sum(ImageStat.Stat(image).var) / 3.0)


def _capture_viewport_png(page: Page) -> bytes:
    session = page.context.new_cdp_session(page)
    try:
        result = session.send(
            "Page.captureScreenshot",
            {"format": "png", "fromSurface": True, "captureBeyondViewport": False},
        )
    finally:
        session.detach()
    return base64.b64decode(result["data"])


def _box_variance(png: bytes, box: dict[str, float], minimum: int = 2) -> float:
    image = Image.open(io.BytesIO(png)).convert("RGB")
    left = max(0, math.floor(box["x"]))
    top = max(0, math.floor(box["y"]))
    right = min(image.width, math.ceil(box["x"] + box["width"]))
    bottom = min(image.height, math.ceil(box["y"] + box["height"]))
    if right - left < minimum or bottom - top < minimum:
        raise AssertionError("Rendered element is not sufficiently visible in the viewport")
    crop = image.crop((left, top, right, bottom))
    return float(sum(ImageStat.Stat(crop).var) / 3.0)


def _scroll(page: Page, locator: Any) -> None:
    locator.wait_for()
    locator.evaluate("element => element.scrollIntoView({block: 'center', inline: 'nearest'})")
    page.wait_for_timeout(250)


def _locator_variance(page: Page, locator: Any) -> float:
    _scroll(page, locator)
    box = locator.bounding_box()
    if not box or box["width"] < 2 or box["height"] < 2:
        raise AssertionError("Rendered element has no measurable pixels")
    return _box_variance(_capture_viewport_png(page), box)


def _click_af_tab(page: Page) -> None:
    tab = page.get_by_role("tab", name=re.compile(r"^AF(?:预测| Prediction)$"))
    tab.click()
    page.get_by_text(re.compile(r"AlphaFold2 (?:单体预测|monomer prediction)"), exact=True).wait_for()


def _choose_prediction(page: Page, prediction_id: str) -> None:
    select = page.get_by_label(re.compile(r"^(?:选择预测|Select prediction)$"))
    _scroll(page, select)
    select.click()
    page.get_by_role("option", name=re.compile(re.escape(prediction_id))).click()
    page.get_by_text("5", exact=True).first.wait_for()


def _fill_public_query(page: Page, fasta: str) -> float:
    textarea = page.get_by_label(re.compile(r"^(?:单条蛋白 FASTA|Single-protein FASTA)$"))
    textarea.fill(fasta)
    textarea.press("Tab")
    page.wait_for_timeout(1200)
    _click_af_tab(page)
    page.get_by_text(re.compile(r"8HFW (?:催化位点映射|catalytic-site mapping)"), exact=True).wait_for()
    mapping_table = page.locator('[data-testid="stDataFrame"]:visible').first
    mapping_table.wait_for()
    variance = _locator_variance(page, mapping_table)
    if variance < 3.0:
        raise AssertionError(f"Catalytic mapping table is visually blank (variance={variance:.3f})")
    return variance


def _nonblank_structure_canvas(page: Page) -> tuple[int, float]:
    heading = page.get_by_text(
        re.compile(r"^(?:按 pLDDT 着色的结构|Structure colored by pLDDT)$"),
        exact=True,
    )
    heading.wait_for()
    heading.evaluate("element => element.scrollIntoView({block: 'start', inline: 'nearest'})")
    page.wait_for_timeout(1800)
    canvas_boxes: list[dict[str, float]] = []
    for frame in page.frames:
        locator = frame.locator("canvas")
        if locator.count() and frame != page.main_frame:
            frame_element = frame.frame_element()
            frame_element.evaluate(
                "element => element.scrollIntoView({block: 'center', inline: 'nearest'})"
            )
            page.wait_for_timeout(300)
        for index in range(locator.count()):
            canvas = locator.nth(index)
            box = canvas.bounding_box()
            if box and box["width"] >= 200 and box["height"] >= 200:
                canvas_boxes.append(box)
    if not canvas_boxes:
        raise AssertionError("No rendered py3Dmol canvas was found")
    variances: list[float] = []
    viewport_png = _capture_viewport_png(page)
    for box in canvas_boxes:
        try:
            variances.append(_box_variance(viewport_png, box, minimum=200))
        except AssertionError:
            continue
    if not variances:
        raise AssertionError("The py3Dmol canvas was not sufficiently visible in the viewport")
    variance = max(variances)
    if variance < 3.0:
        raise AssertionError(f"py3Dmol canvas is visually blank (variance={variance:.3f})")
    return len(canvas_boxes), variance


def _validate_charts(page: Page) -> tuple[int, float]:
    _scroll(page, page.get_by_text(re.compile(r"^(?:pLDDT 曲线|pLDDT profile)$"), exact=True))
    page.wait_for_timeout(800)
    charts = page.locator('[data-testid="stPlotlyChart"]:visible')
    count = charts.count()
    if count < 2:
        raise AssertionError(f"Expected pLDDT and PAE Plotly charts, found {count}")
    variance = min(_locator_variance(page, charts.nth(index)) for index in range(2))
    if variance < 3.0:
        raise AssertionError(f"A prediction chart is visually blank (variance={variance:.3f})")
    return count, variance


def _layout_checks(page: Page) -> dict[str, Any]:
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


def _exercise_viewport(
    browser: Any,
    url: str,
    name: str,
    viewport: dict[str, int],
    fasta: str,
    screenshot_dir: Path,
    prediction_id: str,
) -> dict[str, Any]:
    context = browser.new_context(viewport=viewport, device_scale_factor=1)
    page = context.new_page()
    page_errors: list[str] = []
    page.on("pageerror", lambda error: page_errors.append(str(error)))
    page.goto(url, wait_until="domcontentloaded", timeout=60_000)
    page.locator('[data-testid="stAppViewContainer"]').wait_for(timeout=60_000)
    _click_af_tab(page)
    mapping_variance = _fill_public_query(page, fasta)
    page.get_by_label(
        re.compile(r"^(?:单条蛋白 FASTA|Single-protein FASTA)$")
    ).evaluate("element => { element.scrollTop = 0; element.setSelectionRange(0, 0); }")

    workbench = page.get_by_text(
        re.compile(r"AlphaFold2 (?:单体预测|monomer prediction)"), exact=True
    )
    _scroll(page, workbench)
    page.wait_for_timeout(750)
    screenshot = screenshot_dir / f"af_prediction_{name}.png"
    page.screenshot(path=str(screenshot), full_page=False)
    with Image.open(screenshot) as image:
        if image.size != (viewport["width"], viewport["height"]):
            raise AssertionError(f"Unexpected screenshot dimensions: {image.size}")

    _choose_prediction(page, prediction_id)
    page.get_by_text(re.compile(r"^(?:模型排名|Model ranking)$"), exact=True).wait_for()
    chart_count, chart_variance = _validate_charts(page)
    canvas_count, canvas_variance = _nonblank_structure_canvas(page)
    layout = _layout_checks(page)
    if page_errors:
        raise AssertionError("Browser page errors: " + " | ".join(page_errors))
    context.close()
    return {
        "viewport": viewport,
        "screenshot": str(screenshot.relative_to(project_root())).replace("\\", "/"),
        "plotly_chart_count": chart_count,
        "catalytic_mapping_table_pixel_variance": round(mapping_variance, 3),
        "minimum_chart_pixel_variance": round(chart_variance, 3),
        "structure_canvas_count": canvas_count,
        "structure_canvas_pixel_variance": round(canvas_variance, 3),
        **layout,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", default="http://127.0.0.1:8501")
    parser.add_argument("--edge", default=str(EDGE))
    parser.add_argument("--prediction-id", default=PREDICTION_ID)
    parser.add_argument("--output", default="results/validation/ui_acceptance_report.json")
    args = parser.parse_args()

    root = project_root()
    screenshot_dir = root / "docs" / "screenshots"
    screenshot_dir.mkdir(parents=True, exist_ok=True)
    fasta = (root / "data" / "query" / "Bur01m_seed.fasta").read_text(encoding="ascii")
    edge = Path(args.edge)
    if not edge.is_file():
        raise FileNotFoundError(f"Microsoft Edge was not found: {edge}")

    results: dict[str, Any] = {}
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(
            executable_path=str(edge),
            headless=True,
            args=["--enable-webgl", "--ignore-gpu-blocklist", "--use-angle=swiftshader"],
        )
        try:
            for name, viewport in VIEWPORTS.items():
                results[name] = _exercise_viewport(
                    browser,
                    args.url,
                    name,
                    viewport,
                    fasta,
                    screenshot_dir,
                    args.prediction_id,
                )
        finally:
            browser.close()

    report = {
        "project_version": "0.4.0",
        "passed": True,
        "prediction_id": args.prediction_id,
        "viewports": results,
    }
    output = Path(args.output)
    if not output.is_absolute():
        output = root / output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    print(f"Report: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
