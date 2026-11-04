"""Streamlit controls for running native designs and handing sequences to existing AF2 tools."""
from __future__ import annotations

import json
import os
import re
import tempfile
from pathlib import Path

import pandas as pd
import streamlit as st

from .config import project_root
from .design import ligand_groups, review_design, run_design
from .design_models import CHECKPOINTS, model_dir
from .prediction_sequence import validate_sequence
from .run_context import make_run_id


def tr(zh, en):
    return zh if st.session_state.get("language", "zh") == "zh" else en


def design_root():
    return Path(os.environ.get("AQUA_DESIGN_ROOT") or project_root() / "results/designs")


def choose_run(key):
    root = design_root()
    runs = sorted((p.parent.name for p in root.glob("*/manifest.json")), reverse=True)
    if not runs:
        st.caption(tr("尚无本地设计记录。", "No local design runs yet."))
        return None
    return root / st.selectbox(tr("设计批次", "Design run"), runs, key=key)


def select_design_for_prediction():
    run = choose_run("af_design_run")
    if run:
        df = pd.read_csv(run / "candidates.tsv", sep="\t")
        cid = st.selectbox(tr("选择设计序列", "Select a designed sequence"), df.candidate_id, key="af_design_id")
        row = df[df.candidate_id == cid].iloc[0]
        return validate_sequence(cid, row.sequence)
    return None


def upload_file(upload, directory, basename):
    if upload.size > 20 * 1024 * 1024:
        raise ValueError("Maximum file size is 20 MiB.")
    dest = Path(directory) / (basename + Path(upload.name).suffix.lower())
    dest.write_bytes(upload.getvalue())
    return dest


def render_local_design(candidates, predictions):
    with st.expander(tr("本地酶序列设计 · AquaDesign", "Local enzyme sequence design · AquaDesign")):
        if not st.toggle(tr("打开序列设计工作台", "Open the design workbench"), key="open_aquadesign"):
            return
        st.caption(tr("先按 Bur01 催化核心保护位点，再用本地模型比较允许的替换。含配体时可比较保留和遮蔽配体原子的两次推理。",
                      "Protect the Bur01 catalytic core and score allowed substitutions locally. With a ligand, compare the same model with and without its atom context."))
        try:
            with tempfile.TemporaryDirectory() as temp:
                render_inputs(candidates, temp)
        except (ValueError, OSError, RuntimeError) as exc:
            st.error(str(exc))
        st.markdown(tr("**候选与结构复核**", "**Candidates and structure review**"))
        run = choose_run("design_results_run")
        if run:
            render_results(run, predictions)


def render_inputs(candidates, temp):
    parents = {"8HFW · Bur01": (project_root() / "data/query/Bur01m.pdb", "A", None)}
    for _, row in candidates.iterrows():
        path = row.get("target_structure_path")
        seq = row.get("final_protein_sequence")
        if isinstance(path, str) and path and Path(path).is_file() and isinstance(seq, str) and seq:
            saved_chain = row.get("target_chain")
            saved_chain = saved_chain if isinstance(saved_chain, str) and saved_chain.strip() else "A"
            label = str(row.get("target", "")) + " · " + saved_chain
            parents[label] = (Path(path), saved_chain, seq)
    parent_id = st.selectbox(tr("亲本（含已保存坐标的候选）", "Parent (including candidates with saved coordinates)"),
                             list(parents), key="design_parent")
    parent, chain, sequence = parents[parent_id]
    upload = st.file_uploader(tr("或上传亲本坐标", "Or upload parent coordinates"), type=["pdb", "cif", "mmcif"], key="design_parent_file")
    if upload:
        parent = upload_file(upload, temp, "parent")
        parent_id, chain, sequence = Path(upload.name).name, "A", None
    chain = st.text_input(tr("亲本蛋白链", "Parent protein chain"), value=chain, key="design_chain_" + parent_id)
    groups = ligand_groups(parent)
    ligands = st.multiselect(tr("参与设计的配体（留空使用 ProteinMPNN）", "Ligands for LigandMPNN (empty for ProteinMPNN)"),
                             list(groups), key="design_ligands")
    priority = "compatibility"
    if ligands:
        priority = st.selectbox(tr("候选排序", "Candidate ranking"), ["compatibility", "ligand_response"],
                                format_func=lambda x: tr("骨架相容性", "Structural compatibility") if x == "compatibility" else tr("配体遮蔽对照", "Ligand masking contrast"),
                                key="design_priority")
    msa = st.file_uploader(tr("可选同源比对（首条为亲本）", "Optional homolog alignment (parent first)"),
                           type=["a3m", "fa", "fas", "faa", "fasta"], key="design_msa")
    fixed = st.text_input(tr("额外保护残基，如 A240、L192（空格或逗号分隔）", "Additional fixed residues, e.g. A240 L192"), key="design_fixed")
    c1, c2, c3 = st.columns(3)
    sites = c1.number_input(tr("最多设计位点", "Maximum sites"), 1, 48, 16, key="design_sites")
    edits = c2.number_input(tr("每条最多替换", "Maximum edits per design"), 1, 3, 2, key="design_edits")
    count = c3.number_input(tr("候选数量上限", "Candidate limit"), 1, 100, 12, key="design_count")
    required = CHECKPOINTS["ligand_mpnn" if ligands else "protein_mpnn"]
    ready = (model_dir() / required).is_file() and (model_dir() / "model_utils.py").is_file()
    if not ready:
        st.info(tr("请先在项目环境安装 PyTorch，并运行 aquaesternet-structmine design-setup 下载权重。",
                   "Install PyTorch and run aquaesternet-structmine design-setup in the project environment."))
    if st.button(tr("生成设计序列", "Generate designs"), disabled=not ready, key="design_generate"):
        with st.spinner(tr("正在进行本地推理…", "Running local inference…")):
            run = run_design(parent, design_root() / make_run_id("design"), chain=chain, ligands=ligands,
                             fixed=re.findall(r"[^,，、\s]+", fixed), msa=upload_file(msa, temp, "msa") if msa else None,
                             site_limit=sites, max_edits=edits, count=count, priority=priority, parent_id=parent_id,
                             expected_sequence=sequence)
        st.session_state["design_results_run"] = run.name
        st.success(tr("已保存本次设计。", "Design run saved."))


def render_results(run, predictions):
    df = pd.read_csv(run / "candidates.tsv", sep="\t")
    manifest = json.loads((run / "manifest.json").read_text(encoding="utf-8"))
    st.caption(f"{manifest['model']['model']} · {manifest['priority']} · {', '.join(manifest['triad_labels'])}")
    st.dataframe(df[["candidate_id", "mutations_pdb", "substitutions", "additive_log_odds", "ligand_context_delta",
                     "gravy", "estimated_charge_pH7"]], hide_index=True, width="stretch")
    st.caption(tr("相容性与配体遮蔽差值描述模型对序列的偏好。辛酸乙酯生成速率和合成/水解比需单独测定。",
                  "Scores describe sequence preference. Ethyl octanoate formation and synthesis/hydrolysis ratios require measurements."))
    for ui_col, name in zip(st.columns(3), ("candidates.fasta", "sites.tsv", "substitutions.tsv")):
        ui_col.download_button(name, (run / name).read_bytes(), file_name=name, key=f"design_dl_{name}")
    cid = st.selectbox(tr("需复核的设计", "Design for structural review"), df.candidate_id, key="design_review_id")
    digest = df.loc[df.candidate_id == cid, "sequence_sha256"].iloc[0]
    matching = {p["prediction_id"]: Path(p["_manifest_path"]) for p in predictions
                if p.get("status") == "completed" and p.get("sequence_sha256") == digest}
    if matching:
        pid = st.selectbox(tr("同序列的已完成 AF2 任务", "Completed AF2 runs for this sequence"), list(matching), key="design_af2_match")
        if st.button(tr("复核本机 AF2 模型", "Review saved AF2 model"), key="design_link_af2"):
            try:
                pm = matching[pid]
                review_design(run, cid, pm.parent / "rank_001.pdb", "A", "predicted", pm)
            except (OSError, ValueError) as exc:
                st.error(str(exc))
    st.caption(tr("可在“AF预测”页选择“本地设计序列”进行预测，或上传已有坐标。", "Choose Local design in AF Prediction, or upload existing coordinates."))
    upload = st.file_uploader(tr("设计序列对应的坐标", "Coordinates for this designed sequence"), type=["pdb", "cif", "mmcif"], key="design_review_file")
    if upload:
        chain = st.text_input(tr("待复核结构的链", "Chain to review"), value="A", key="design_review_chain")
        origin = st.selectbox(tr("坐标来源", "Coordinate origin"), ["unspecified", "predicted", "experimental", "template_model"], key="design_origin")
        if st.button(tr("核对序列并比对结构", "Verify sequence and compare structures"), key="design_check"):
            try:
                with tempfile.TemporaryDirectory() as temp:
                    review_design(run, cid, upload_file(upload, temp, "target"), chain, origin)
            except (ValueError, OSError) as exc:
                st.error(str(exc))
    reports = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((run / "reviews").glob("*.json"))]
    if reports:
        st.dataframe(pd.DataFrame(reports), hide_index=True, width="stretch")
