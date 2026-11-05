"""Browser entry points for submitting sequences, designing candidates and following inference."""

from __future__ import annotations

import json
import os
import tempfile

import pandas as pd
import streamlit as st

from .config import project_root
from .design_view import design_root, render_inputs
from .property_assets import ASSETS, asset_root, runtime_python
from .property_fasta import collect_fasta
from .property_jobs import PropertyJobs
from .property_screen import DEFAULT_SUBSTRATES, parse_substrates, profile_fasta, profile_runs, profiles_root
from .property_view import show_report

BUR01_RUNS = ("20261105_Bur01_backbone", "20261105_Bur01_octanoate_ethanol")

STYLE = """
<style>
 [data-testid="stAppViewContainer"] {background:#f7fafb;color:#203746}
 [data-testid="stMainBlockContainer"] {max-width:1300px;padding-top:2.8rem;padding-bottom:3rem}
 [data-testid="stHeader"] {background:#f7fafbf5}
 [data-testid="stAppDeployButton"], #MainMenu {display:none}
 [data-testid="stVerticalBlockBorderWrapper"] {background:white;border-color:#cfdee3}
 .aqua-heading {display:flex;align-items:baseline;gap:22px;padding:12px 0 13px;border-bottom:2px solid #167c81}
 .aqua-heading .wordmark {font:600 clamp(25px,3.2vw,42px)/1.2 Georgia,serif;color:#224674;letter-spacing:-1.6px}
 .wordmark i {color:#167c81}
 .aqua-heading h1 {font:650 clamp(21px,2.7vw,34px)/1.3 "Microsoft YaHei","Segoe UI",sans-serif;color:#203746;margin:0;padding:0}
 .aqua-sub {color:#526b79;font-size:1rem;line-height:1.55;margin:9px 0 2px}
 .aqua-models {display:grid;grid-template-columns:5.5rem 1fr;gap:6px 12px;font-size:.92rem;color:#385461;margin:.1rem 0 .65rem}
 .aqua-models dt {font:600 .93rem/1.55 Consolas,monospace;color:#224674}.aqua-models dd {margin:0}
 [data-testid="stButton"] button[kind="primary"],[data-testid="stFormSubmitButton"] button[kind="primary"] {background:#146f76;color:white;border-color:#146f76}
 [data-testid="stTextArea"] textarea {font-family:Consolas,"Microsoft YaHei",monospace;line-height:1.5}
 [data-testid="stMarkdownContainer"] h3 {color:#224674}
 .aqua-triad {display:flex;align-items:center;gap:.55rem;color:#167c81;font:600 1rem Consolas,monospace;margin:.5rem 0 .75rem}
 .aqua-triad i {height:2px;background:#a9c9cd;width:2.6rem}
 @media(max-width:760px){[data-testid="stMainBlockContainer"]{padding:3.5rem .9rem 2rem}.aqua-heading{flex-wrap:wrap;gap:.3rem 1rem}}
</style>
"""


def heading(title, subtitle):
    st.html(STYLE)
    st.html(
        f'<div class="aqua-heading"><span class="wordmark"><i>Aqua</i>EsterNet</span>'
        f'<h1>{title}</h1></div><p class="aqua-sub">{subtitle}</p>'
    )


@st.cache_resource
def jobs():
    return PropertyJobs()


def go(section):
    st.switch_page(st.session_state["_aqua_pages"][section])


def load_example():
    st.session_state["web_fasta"] = (project_root() / "data/query/Bur01m_seed.fasta").read_text(encoding="utf-8")
    st.session_state["_web_upload_revision"] = st.session_state.get("_web_upload_revision", 0) + 1


def prepare_batch(names):
    st.session_state["web_source"] = "设计批次"
    st.session_state["web_batches"] = list(names)


def predict_page():
    heading("酶性质预测", "提交蛋白序列和底物，在本机计算可溶表达倾向、熔解温度与酶–底物匹配分数。")
    available = sorted((p.parent.name for p in design_root().glob("*/manifest.json")), reverse=True)
    pending = st.session_state.pop("_web_use_design", None)
    if pending:
        prepare_batch([pending])
    st.session_state.setdefault("web_source", "序列输入")
    source = st.segmented_control("输入来源", ["序列输入", "设计批次"], required=True, key="web_source")
    left, right = st.columns([1.42, 1], gap="large")
    text, selected, parent, input_error = "", [], False, None
    with left.container(border=True):
        if source == "序列输入":
            st.subheader("蛋白序列")
            st.button("载入 Bur01 亲本", on_click=load_example, key="web_example")
            text = st.text_area(
                "粘贴序列或多条 FASTA",
                key="web_fasta",
                height=145,
                placeholder=">蛋白名称\n氨基酸序列",
                help="支持 20 种标准氨基酸，单次最多 100 条。",
            )
            with st.expander("或上传 FASTA 文件"):
                uploaded = st.file_uploader(
                    "选择序列文件",
                    type=["fa", "fas", "faa", "fasta", "txt"],
                    max_upload_size=2,
                    key=f"web_fasta_file_{st.session_state.get('_web_upload_revision', 0)}",
                )
            if uploaded is not None:
                try:
                    text = uploaded.getvalue().decode("utf-8-sig")
                    st.caption(f"本次使用上传文件 {uploaded.name}")
                except UnicodeDecodeError:
                    input_error = "请将 FASTA 文件保存为 UTF-8 编码后重新上传。"
            parent = st.checkbox("以首条序列为亲本，比较等长突变体", key="web_parent")
            if text and not input_error:
                try:
                    sequences, _, _, sources = collect_fasta(text, parent)
                    n = len(sources["FASTA"]["records"])
                    st.caption(f"已识别 {n} 条记录 · {len(sequences)} 条独立序列")
                    if any(len(seq) > 510 for seq in sequences.values()):
                        st.info("超过 510 个残基的序列仅计算可溶表达与底物匹配，Tm 留空。")
                except ValueError as exc:
                    input_error = str(exc)
        else:
            st.subheader("设计批次")
            if set(BUR01_RUNS).issubset(available):
                st.button(
                    "载入 Bur01 的 24 条候选", on_click=prepare_batch, args=(BUR01_RUNS,), key="web_example_batch"
                )
            selected = st.multiselect(
                "参与预测的批次", available, key="web_batches", placeholder="选择一个或多个设计批次"
            )
            st.caption("程序从每个批次的坐标读取亲本，核对候选的序列与替换编号。")
            if not available and st.button("先生成候选"):
                go("design")
        if input_error:
            st.error(input_error)
    with right.container(border=True):
        st.subheader("底物与模型")
        st.session_state.setdefault("web_chem", "\n".join(f"{n}={s}" for n, s in DEFAULT_SUBSTRATES))
        raw = st.text_area(
            "底物名称与 SMILES",
            height=93,
            key="web_chem",
            help="每行一种底物，按 名称=SMILES 填写。不同分子分别计算。",
        )
        st.html(
            '<dl class="aqua-models"><dt>SWI</dt><dd>大肠杆菌可溶表达倾向</dd>'
            "<dt>TmProt</dt><dd>蛋白熔解温度（°C）</dd>"
            "<dt>EnzRank</dt><dd>每种底物的独立配对分数</dd></dl>"
        )
        ready = runtime_python().is_file() and all((asset_root() / p).is_file() for p in ASSETS)
        if not ready:
            st.info(
                "本机尚未装全性质模型。请先运行 python -m bur01_structmine.cli properties-setup --install-runtime。"
            )
        if st.button("提交本地预测", type="primary", width="stretch", disabled=not ready or bool(input_error)):
            try:
                substrates = parse_substrates([s.strip() for s in raw.splitlines() if s.strip()])
                if source == "序列输入":
                    _, _, _, src = collect_fasta(text, parent)
                    records = src["FASTA"]["records"]
                    title = f"{records[0]['id']} · {len(records)} 条序列"
                    job = jobs().submit(
                        title, profile_fasta, text=text, first_as_parent=parent, substrates=substrates
                    )
                else:
                    if not selected:
                        raise ValueError("请至少选择一个设计批次。")
                    job = jobs().submit(
                        "、".join(selected),
                        profile_runs,
                        directories=[design_root() / s for s in selected],
                        substrates=substrates,
                    )
            except ValueError as exc:
                st.error(str(exc))
            else:
                st.session_state["_web_select_job"] = job
                go("results")


def design_page():
    heading("结构约束设计", "选择 Bur01 亲本或上传同折叠结构，指定可参与设计的配体和固定残基，再生成候选。")
    left, right = st.columns([1.45, 1], gap="large")
    new_run = None
    with left.container(border=True):
        st.subheader("亲本与设计参数")
        try:
            with tempfile.TemporaryDirectory() as tmp:
                new_run = render_inputs(pd.DataFrame(), tmp, compact=True)
        except (ValueError, OSError, RuntimeError, KeyError) as exc:
            st.error(str(exc))
    with right.container(border=True):
        st.subheader("设计约束")
        st.caption("Bur01 参考结构 · 8HFW A 链")
        st.html('<div class="aqua-triad">Ser111<i></i>Asp241<i></i>His274</div>')
        st.write("固定催化三联体、Ser 中心模体和紧邻催化残基的位点。其他亲本按结构比对映射编号。")
        st.caption(
            "ProteinMPNN 评价给定骨架上的替换。选中有坐标的配体后，LigandMPNN 比较保留与遮蔽配体时的序列偏好。"
        )
        st.divider()
        st.subheader("候选与后续预测")
        runs = sorted((p.parent.name for p in design_root().glob("*/manifest.json")), reverse=True)
        if not runs:
            st.caption("生成的候选会出现在这里。")
            return
        if new_run:
            st.session_state["web_design_result"] = new_run.name
        selected = st.selectbox("已保存的设计", runs, key="web_design_result")
        run = design_root() / selected
        rows = pd.read_csv(run / "candidates.tsv", sep="\t")
        st.dataframe(
            rows[["mutations_pdb", "substitutions"]].rename(
                columns={"mutations_pdb": "结构编号替换", "substitutions": "替换数"}
            ),
            width="stretch",
            hide_index=True,
            height=180,
        )
        if st.button("用本批候选预测性质", type="primary", width="stretch"):
            st.session_state["_web_use_design"] = selected
            go("predict")
        st.download_button(
            "下载候选 FASTA", (run / "candidates.fasta").read_bytes(), file_name=selected + ".fasta"
        )


@st.fragment(run_every=2)
def watch_job(key):
    job = next((j for j in jobs().snapshot() if j.key == key), None)
    if not job or job.state not in {"running", "queued"}:
        st.rerun()
    indicator = st.status(job.message, state="running", expanded=True)
    indicator.write(job.title)
    indicator.caption("本机正在处理任务。可以切换页面，稍后从任务记录打开结果。")


def results_page():
    heading("任务与结果", "查看本机预测进度，按序列和底物检查模型输出，下载原始表格。")
    current = {j.key: j for j in jobs().snapshot()}
    saved = sorted(profiles_root().glob("*/manifest.json"), key=lambda p: p.stat().st_mtime_ns, reverse=True)
    reports = {p.parent.name: p.parent for p in saved if (p.parent / "index.html").is_file()}
    keys = list(current) + [key for key in reports if key not in current]
    if not keys:
        st.info("还没有预测记录。")
        if st.button("提交第一条序列", type="primary"):
            go("predict")
        return
    pending = st.session_state.pop("_web_select_job", None)
    if pending in keys:
        st.query_params["run"] = pending
    requested = st.query_params.get("run")
    if requested in keys and st.session_state.get("web_report_choice") != requested:
        st.session_state["web_report_choice"] = requested
    states = {"queued": "排队中", "running": "计算中", "failed": "失败", "complete": "已完成"}

    def label(key):
        if key in current:
            job = current[key]
            return f"{key}　·　{states[job.state]}　·　{job.title}"
        return f"{key}　·　已保存"

    def select_result():
        st.query_params["run"] = st.session_state["web_report_choice"]

    selected = st.selectbox("选择任务", keys, format_func=label, key="web_report_choice", on_change=select_result)
    if selected != requested:
        select_result()
    job = current.get(selected)
    if job and job.state in {"running", "queued"}:
        watch_job(selected)
    elif job and job.state == "failed":
        st.error("本次预测失败。 " + job.message)
        if st.button("返回修改输入"):
            go("predict")
    elif selected in reports:
        manifest = json.loads((reports[selected] / "manifest.json").read_text(encoding="utf-8"))
        st.caption(f"已完成 · {manifest['unique_sequences']} 条独立序列 · {manifest['pair_count']} 个酶–底物配对")
        show_report(reports[selected])


def mining_page():
    """The original workbench continues in its entrypoint after navigation."""


def navigate(*, standalone):
    # Keep editable inputs when a user checks another page and returns to the draft.
    for key in ("web_fasta", "web_chem", "web_parent", "web_source", "web_batches", "web_report_choice"):
        if key in st.session_state:
            st.session_state[key] = st.session_state[key]
    default_models = standalone or os.environ.get("AQUA_MODEL_HOME") == "1"
    pages = {
        "predict": st.Page(predict_page, title="性质预测", url_path="predict", default=default_models),
        "design": st.Page(design_page, title="结构约束设计", url_path="design"),
        "results": st.Page(results_page, title="任务与结果", url_path="results"),
        "mine": st.Page(
            "app.py" if standalone else mining_page,
            title="结构挖掘与 AF2",
            url_path="mining",
            default=not default_models,
        ),
    }
    st.session_state["_aqua_pages"] = pages
    page = st.navigation(list(pages.values()), position="top")
    page.run()
    return page == pages["mine"]
