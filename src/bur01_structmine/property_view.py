"""Streamlit controls shared by the standalone page and the design workbench."""

import streamlit as st

from .design_view import design_root
from .property_assets import asset_root, runtime_python
from .property_screen import parse_substrates, profile_runs, profiles_root


def show_properties(current_run=None):
    root = design_root()
    runs = sorted(p.parent.name for p in root.glob("*/manifest.json"))
    reports = sorted(
        (p.parent for p in profiles_root().glob("*/manifest.json") if (p.parent / "index.html").is_file()),
        reverse=True,
    )
    with st.expander("本地计算设置 · 选择设计批次和底物", expanded=not reports):
        preferred = [current_run.name] if current_run and current_run.name in runs else runs[-2:]
        selected = st.multiselect("设计批次", runs, default=preferred, key="property_inputs")
        text = st.text_area(
            "底物名称与 SMILES（每行一种）", "octanoic_acid=CCCCCCCC(=O)O\nethanol=CCO", key="property_chem"
        )
        ready = runtime_python().is_file() and (asset_root() / "esm/model.safetensors").is_file()
        if not ready:
            st.info("请在项目环境运行 python -m bur01_structmine.cli properties-setup --install-runtime。")
        if st.button(
            "用本机 SWI、TmProt 和 EnzRank 重新计算", disabled=not (ready and selected), key="property_run"
        ):
            msg = st.empty()
            try:
                substrates = parse_substrates([line.strip() for line in text.splitlines() if line.strip()])
                report = profile_runs([root / name for name in selected], substrates, progress=msg.info)
                st.session_state["property_report"] = report.name
                st.rerun()
            except (ValueError, OSError, KeyError) as exc:
                msg.error(str(exc))
    if not reports:
        st.info("完成一次本地预测后，此处将显示候选与亲本的对照。")
        return
    names = [p.name for p in reports]
    chosen = st.selectbox("已保存的性质预测", names, key="property_report")
    show_report(profiles_root() / chosen)


def show_report(directory):
    html = (directory / "index.html").read_text(encoding="utf-8")
    st.iframe(html, height="content")
    with st.expander("下载结果与模型记录"):
        for name, label in (
            ("index.html", "交互结果网页"),
            ("protein_predictions.tsv", "序列性质原始表"),
            ("substrate_predictions.tsv", "酶–底物配对表"),
            ("candidate_comparison.tsv", "候选与亲本差值表"),
            ("input.fasta", "本次提交的序列"),
            ("manifest.json", "模型版本与文件记录"),
        ):
            if (directory / name).is_file():
                st.download_button(
                    label, (directory / name).read_bytes(), file_name=name, key="property_download_" + name
                )
