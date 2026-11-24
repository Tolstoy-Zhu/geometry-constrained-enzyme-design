"""Local-only browser acceptance for the v0.5 evidence workbench."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import time
from pathlib import Path
from urllib.parse import urlparse

from PIL import Image, ImageStat
from playwright.sync_api import sync_playwright

from bur01_structmine.config import project_root

EDGE = r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"


def _canvas_check(page):
    for frame in page.frames[1:]:
        canvas = frame.locator("canvas").first
        if canvas.count() and canvas.is_visible():
            canvas.scroll_into_view_if_needed()
            page.wait_for_timeout(400)
            data = canvas.screenshot()
            variance = sum(ImageStat.Stat(Image.open(io.BytesIO(data)).convert("RGB")).var) / 3
            if variance > 3:
                box = canvas.bounding_box()
                page.mouse.move(box["x"] + box["width"] * .5, box["y"] + box["height"] * .5)
                page.mouse.down()
                page.mouse.move(box["x"] + box["width"] * .7, box["y"] + box["height"] * .6, steps=8)
                page.mouse.up()
                page.wait_for_timeout(500)
                moved = canvas.screenshot()
                assert hashlib.sha256(data).digest() != hashlib.sha256(moved).digest(), "Structure did not rotate"
                return variance
    raise AssertionError("No nonblank molecular canvas")


def nonblank(page):
    from playwright.sync_api import Error
    for attempt in range(10):
        try:
            return _canvas_check(page)
        except Error:
            if attempt == 9:
                raise
            page.wait_for_timeout(800)


def choose_step(page, step, mobile):
    expand = page.locator('[data-testid="stExpandSidebarButton"]')
    if mobile and expand.is_visible():
        expand.click()
        page.wait_for_timeout(600)
    radio = page.get_by_role("radio", name=re.compile(f"^{step} "))
    if not radio.is_visible():
        page.locator('[data-testid="stExpandSidebarButton"]').click()
    radio.locator("xpath=ancestor::label[1]").click()
    page.wait_for_timeout(1200)
    if mobile:
        close_sidebar(page)


def close_sidebar(page):
    close = page.locator('[data-testid="stSidebarCollapseButton"] button')
    box = close.bounding_box() if close.count() else None
    if box and 0 <= box["x"] < page.viewport_size["width"]:
        close.click()
        page.wait_for_timeout(600)


def check_errors(page):
    assert page.locator('[data-testid="stException"]').count() == 0, page.locator('[data-testid="stException"]').all_text_contents()


def wait_ready(page):
    page.wait_for_timeout(300)
    page.wait_for_function("document.querySelector('[data-testid=stApp]')?.getAttribute('data-test-script-state') === 'notRunning'", timeout=60000)
    page.wait_for_timeout(600)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--workbench-url")
    parser.add_argument("--only-workbench", action="store_true")
    args = parser.parse_args()
    out = project_root() / "docs/validation/v0.5"
    out.mkdir(parents=True, exist_ok=True)
    report = {"demo": [], "external_requests": [], "workbench": []}
    if args.only_workbench and (out / "ui_report.json").is_file():
        previous = json.loads((out / "ui_report.json").read_text(encoding="utf-8"))
        report["demo"] = previous.get("demo", [])
        report["external_requests"] = previous.get("external_requests", [])
    with sync_playwright() as p:
        browser = p.chromium.launch(executable_path=EDGE, headless=True, args=["--enable-webgl", "--use-gl=angle", "--use-angle=swiftshader-webgl", "--disable-gpu-sandbox"])
        for device, width, height in ([] if args.only_workbench else [("desktop", 1440, 900), ("mobile", 390, 844)]):
            context = browser.new_context(viewport={"width": width, "height": height}, device_scale_factor=1)
            def route(request):
                if urlparse(request.request.url).hostname not in {"127.0.0.1", "localhost"}:
                    report["external_requests"].append(request.request.url)
                    request.abort()
                else:
                    request.continue_()
            context.route("**/*", route)
            page = context.new_page()
            page.goto(args.url, wait_until="networkidle")
            page.get_by_role("heading", name="AquaEsterNet-StructMine", exact=True).wait_for(timeout=30000)
            for language in ["zh", "en"]:
                if language == "en":
                    expand = page.locator('[data-testid="stExpandSidebarButton"]')
                    if device == "mobile" and expand.is_visible():
                        expand.click()
                        page.wait_for_timeout(600)
                    page.get_by_text("English", exact=True).click()
                    page.wait_for_timeout(700)
                for step in range(1, 6):
                    choose_step(page, step, device == "mobile")
                    check_errors(page)
                    variance = nonblank(page) if step in {1, 3} else None
                    if step == 5:
                        label = "计算 Top-k 命中率" if language == "zh" else "Calculate Top-k hit rate"
                        page.get_by_text(label, exact=True).click()
                        page.wait_for_timeout(800)
                        check_errors(page)
                    page.keyboard.press("Tab")
                    focus = page.evaluate("document.activeElement.tagName")
                    horizontal = page.evaluate("document.documentElement.scrollWidth > innerWidth + 2")
                    assert not horizontal, "Page overflows viewport"
                    page.screenshot(path=str(out / f"demo_{device}_{language}_{step}.png"), full_page=False)
                    report["demo"].append({"device": device, "language": language, "step": step, "canvas_variance": variance, "focus_tag": focus, "page_overflow": horizontal})
            context.close()
        if not args.only_workbench:
            (out / "ui_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        if args.workbench_url:
            views = [(device, width, height, language) for device, width, height in [("desktop", 1440, 900), ("mobile", 390, 844)] for language in ["zh", "en"]]
            for device, width, height, language in views:
                print(f"Checking workbench: {device}/{language}", flush=True)
                context = browser.new_context(viewport={"width": width, "height": height}, accept_downloads=True)
                page = context.new_page()
                page.goto(args.workbench_url, wait_until="networkidle")
                page.get_by_role("tab", name="实验评估", exact=True).wait_for(timeout=60000)
                if language == "en":
                    expand = page.locator('[data-testid="stExpandSidebarButton"]')
                    if expand.is_visible():
                        expand.click()
                        page.wait_for_timeout(600)
                    page.get_by_text("English", exact=True).click()
                    page.wait_for_timeout(1200)
                if device == "mobile":
                    close_sidebar(page)
                page.get_by_role("tab", name="实验评估" if language == "zh" else "Experimental assessment", exact=True).click()
                page.locator('[data-testid="stFileUploader"]:visible input[type="file"]').set_input_files(str(project_root() / "data/examples/simulated/experimental_measurements.csv"))
                export_label = "导出完整评估 Excel" if language == "zh" else "Export full evaluation Excel"
                page.get_by_text(export_label, exact=True).wait_for()
                check_errors(page)
                charts = page.locator('[data-testid="stPlotlyChart"]:visible')
                assert charts.count() == 2
                chart_variances = []
                for chart_index in range(charts.count()):
                    chart = charts.nth(chart_index)
                    chart.scroll_into_view_if_needed()
                    pixels = chart.screenshot()
                    chart_variances.append(sum(ImageStat.Stat(Image.open(io.BytesIO(pixels)).convert("RGB")).var) / 3)
                assert min(chart_variances) > 3
                top_k = page.get_by_role("checkbox", name="计算 Top-k 命中率" if language == "zh" else "Calculate Top-k hit rate", exact=True)
                top_k.focus()
                top_k.press("Space")
                wait_ready(page)
                assert top_k.is_checked()
                with page.expect_download() as download:
                    page.get_by_text(export_label, exact=True).click()
                download.value.save_as(out / f"evaluation_{device}_{language}.xlsx")
                page.screenshot(path=str(out / f"workbench_{device}_{language}_experiments.png"))
                page.get_by_role("tab", name="结构" if language == "zh" else "Structure", exact=True).click()
                page.get_by_text("底物与位点" if language == "zh" else "Substrates and sites", exact=True).first.click()
                page.locator('[data-testid="stFileUploader"]:visible input[type="file"]').set_input_files(str(project_root() / "data/benchmark/structures/1XKL.pdb"))
                wait_ready(page)
                check_errors(page)
                variance = nonblank(page)
                focus_label = "关注" if language == "zh" else "Focus"
                page.get_by_label(focus_label, exact=True).click()
                page.get_by_role("option", name="A:GLU3", exact=True).click()
                wait_ready(page)
                with page.expect_download() as download:
                    page.get_by_text("导出位点与底物记录" if language == "zh" else "Export sites and substrate record", exact=True).click()
                export_path = out / f"sites_{device}_{language}.json"
                download.value.save_as(export_path)
                sites_export = json.loads(export_path.read_text(encoding="utf-8"))
                assert len(sites_export["sites"]) == 1
                assert sites_export["coordinate_status"] == "coordinates_present"
                page.screenshot(path=str(out / f"workbench_{device}_{language}_sites.png"))
                report["workbench"].append({"device": device, "language": language, "import_export": True, "keyboard_toggle": True,
                                             "chart_variances": chart_variances, "site_marks_exported": len(sites_export["sites"]), "canvas_variance": variance})
                context.close()
        browser.close()
    assert not report["external_requests"], report["external_requests"]
    (out / "ui_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
