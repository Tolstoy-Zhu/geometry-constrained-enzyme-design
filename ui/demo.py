"""Read-only public Bur01 demonstration. No production service imports."""
from __future__ import annotations

import json

import pandas as pd
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

from bur01_structmine.config import project_root
from bur01_structmine.local_viewer import offline_viewer_html
from bur01_structmine.run_context import file_sha256
from bur01_structmine.visualization import align_structure_overlay
from bur01_structmine.workbench_views import experimental_view, score_details, structure_sites_view, tr

ROOT = project_root()
DATA = ROOT / "data" / "demo"
st.set_page_config(page_title="AquaEsterNet | Bur01", layout="wide", initial_sidebar_state="expanded")
st.markdown("""<style>
.stApp{background:#F7F9F8;color:#17211E;font-family:'IBM Plex Sans','Microsoft YaHei',sans-serif}
.block-container{padding:2rem 2.5rem 3rem;max-width:1400px}
h1{font-size:2rem!important;letter-spacing:0!important}h2,h3{letter-spacing:0!important}
button:focus-visible,input:focus-visible{outline:3px solid #C77D00!important;outline-offset:2px}
[data-testid="stMetricValue"]{font-size:1.4rem;white-space:normal}
.stButton button,.stDownloadButton button{border-radius:6px}iframe{max-width:100%}
.stAppDeployButton{display:none}
@media(max-width:640px){.block-container{padding:4rem .8rem 2rem}h1{font-size:1.4rem!important}
[data-testid="stHorizontalBlock"]{flex-wrap:wrap!important}[data-testid="stColumn"]{min-width:140px!important}}
</style>""", unsafe_allow_html=True)
language = st.sidebar.segmented_control("Language", ["中文", "English"], default="中文", key="demo_language")
st.session_state.language = "en" if language == "English" else "zh"
st.sidebar.title("Bur01 / 8HFW")
steps = [tr("1 输入结构", "1 Input structure"), tr("2 候选筛选", "2 Candidate filtering"),
         tr("3 结构比较", "3 Structure comparison"), tr("4 评分解释", "4 Score explanation"),
         tr("5 实验评估", "5 Experimental assessment")]
step = st.sidebar.radio(tr("演示步骤", "Demonstration steps"), list(range(5)), format_func=lambda i: steps[i])
st.sidebar.caption(tr("只读演示 · 不提交计算任务", "Read-only demo · no task submission"))
st.title("AquaEsterNet-StructMine")
st.caption(tr("Bur01 结构挖掘", "Bur01 structure mining") + "  |  v0.5.0.dev0  |  PDB 8HFW:A")
if not (DATA / "demo_manifest.json").exists():
    st.error(tr("演示资料尚未准备，需先完成固定基准验证。", "Demo data is not prepared; verified benchmark required."))
    st.stop()
manifest = json.loads((DATA / "demo_manifest.json").read_text(encoding="utf-8"))
for relative, digest in manifest["files"].items():
    path = (ROOT / relative).resolve()
    if not path.is_relative_to(ROOT) or not path.is_file() or file_sha256(path) != digest:
        st.error(tr("演示资料校验失败：", "Demo checksum failed: ") + relative)
        st.stop()
table = pd.read_csv(DATA / "candidates.tsv", sep="\t")
query = ROOT / "data" / "query" / "Bur01m.pdb"
st.subheader(steps[step])
if step == 0:
    columns = st.columns(3)
    columns[0].metric(tr("结构", "Structure"), "8HFW:A")
    columns[1].metric(tr("残基范围", "Residue range"), "12-294")
    columns[2].metric(tr("催化位点", "Catalytic sites"), "S111 / D241 / H274")
    viewer = py3Dmol.view(width="100%", height=450)
    viewer.addModel(query.read_text(), "pdb")
    viewer.setStyle({}, {"cartoon": {"color": "#0F6B5B"}})
    viewer.addStyle({"resi": [111, 241, 274]}, {"stick": {"color": "#B94A48", "radius": 0.25}})
    active_site = st.segmented_control(tr("催化残基", "Catalytic residue"), ["Ser111", "Asp241", "His274"], default="Ser111")
    if active_site:
        position = {"Ser111": 111, "Asp241": 241, "His274": 274}[active_site]
        viewer.addLabel(active_site, {"fontSize": 13, "backgroundColor": "#17211E"}, {"chain": "A", "resi": position})
    viewer.setBackgroundColor("#F7F9F8")
    viewer.zoomTo()
    components.html(offline_viewer_html(viewer), height=470)
    st.caption(tr("公开实验结构；本页不进行新的结构预测。", "Public experimental structure; no new prediction on this page."))
elif step == 1:
    st.caption(tr("真实 Foldseek 固定基准计算记录", "Real Foldseek fixed-benchmark computation") + " · " + manifest["source_run_id"])
    minimum = st.slider("min(qTM,tTM)", 0.0, 1.0, 0.0, 0.05)
    only_pass = st.checkbox(tr("仅显示通过默认结构阈值的候选", "Only candidates passing default structural thresholds"))
    shown = table[table[["qtmscore", "ttmscore"]].min(axis=1) >= minimum]
    if only_pass:
        shown = shown[shown.default_filter_pass]
    st.dataframe(shown[["benchmark_label", "target", "qtmscore", "ttmscore", "rmsd", "qcov", "tcov", "default_filter_pass", "filter_exclusion_reasons"]], hide_index=True, width="stretch")
    st.success(tr("自比对、阳性排序、P00520 排除检查：通过", "Self-alignment, positive ranking and P00520 rejection: passed"))
    st.json(manifest["benchmark_checks"], expanded=False)
elif step in {2, 3}:
    selected = st.selectbox(tr("候选结构", "Candidate structure"), list(range(len(table))), format_func=lambda i: str(table.iloc[i].target))
    row = table.iloc[selected]
    if step == 3:
        score_details(row)
    else:
        target = ROOT / row.demo_structure_path
        overlay = align_structure_overlay(query, target, query_chain="A", target_chain=str(row.target_chain))
        columns = st.columns(3)
        columns[0].metric("qTM", f"{overlay.qtmscore:.3f}")
        columns[1].metric("tTM", f"{overlay.ttmscore:.3f}")
        columns[2].metric("RMSD", f"{overlay.rmsd:.3f} Å")
        st.caption(tr("当前叠合：本地 TM-align；候选筛选页保留 Foldseek 原始指标。", "Current overlay: local TM-align; candidate filtering retains original Foldseek metrics."))
        viewer = py3Dmol.view(width="100%", height=440)
        viewer.addModel(overlay.query_pdb, "pdb")
        viewer.addModel(overlay.target_pdb, "pdb")
        viewer.setStyle({"model": 0}, {"cartoon": {"color": "#0F6B5B"}})
        viewer.setStyle({"model": 1}, {"cartoon": {"color": "#C77D00", "opacity": 0.7}})
        viewer.addStyle({"model": 0, "resi": [111, 241, 274]}, {"stick": {"color": "#B94A48"}})
        viewer.setBackgroundColor("#F7F9F8")
        viewer.zoomTo()
        components.html(offline_viewer_html(viewer), height=460)
        st.write("Ser111 / Asp241 / His274 → " + overlay.mapped_catalytic_residues)
        st.code(overlay.aligned_query + "\n" + overlay.alignment_match_line + "\n" + overlay.aligned_target, language=None, wrap_lines=False)
        with st.expander(tr("底物与位点", "Substrates and sites")):
            structure_sites_view(query, key="demo_sites", read_only=True)
else:
    experimental_view(demo=True, key="demo_experiments")
with st.expander(tr("计算来源", "Computation provenance")):
    st.json(manifest, expanded=False)
