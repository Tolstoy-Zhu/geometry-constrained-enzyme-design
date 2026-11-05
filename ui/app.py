"""Bilingual local workbench for AquaEsterNet-StructMine."""
from __future__ import annotations

import html
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components
import yaml

from bur01_structmine.config import load_config, project_root
from bur01_structmine.design_view import render_local_design, tr as model_text
from bur01_structmine.structmine_model_tabs import (
    candidate_property_button, open_model_tab, show_property_tab, workbench_tabs,
)
from bur01_structmine.doctor import run_doctor
from bur01_structmine.http_utils import CachedSession
from bur01_structmine.models import PredictionManifest, ProjectConfig
from bur01_structmine.pipeline import reannotate_run
from bur01_structmine.prediction import request_prediction_cancel
from bur01_structmine.prediction_sequence import (
    SequenceValidationError,
    map_8hfw_catalytic_sites,
    parse_single_fasta,
    validate_sequence,
)
from bur01_structmine.results_parser import parse_pdb_chain
from bur01_structmine.remote_auth import load_remote_access, verify_remote_password
from bur01_structmine.run_context import make_run_id, mark_run_cancelled
from bur01_structmine.visualization import StructureOverlay, align_structure_overlay

ROOT = project_root()
RUN_ROOT = Path(os.environ.get("AQUA_RUN_ROOT", ROOT / "results" / "runs"))
UI_CONFIG_ROOT = ROOT / "results" / "ui_configs"
UI_STRUCTURE_CACHE = ROOT / ".cache" / "ui_structures"
PREDICTION_ROOT = Path(os.environ.get("AQUA_PREDICTION_ROOT", ROOT / "results" / "predictions"))
NR_RUN_ROOT = Path(os.environ.get("AQUA_NR_RUN_ROOT", ROOT / "results" / "nr_runs"))
UI_PREDICTION_INPUTS = ROOT / "results" / "ui_prediction_inputs"
UI_NR_INPUTS = ROOT / "results" / "ui_nr_inputs"
REMOTE_ACCESS_PATH = Path(os.environ.get("AQUA_REMOTE_ACCESS_PATH", ROOT / ".cache" / "remote_ui" / "access.json"))
FRAGMENT_RUN_EVERY = None if os.environ.get("AQUA_TEST_MODE") == "1" else 2

COLORS = {
    "background": "#F7F9F8",
    "ink": "#17211E",
    "green": "#0F6B5B",
    "amber": "#C77D00",
    "red": "#B94A48",
    "line": "#D6DEDA",
    "muted": "#61706A",
    "panel": "#EEF3F1",
}

CONFIG_OPTIONS = {
    "benchmark_foldseek": "configs/benchmark_foldseek.yaml",
    "alphafold_swissprot": "configs/alphafold_swissprot.yaml",
    "pdb": "config.yaml",
    "benchmark_tmtools": "configs/benchmark_tmtools.yaml",
}

TEXT = {
    "zh": {
        "settings": "运行设置",
        "language": "语言",
        "query": "查询结构",
        "chain": "查询链",
        "database": "搜索数据库",
        "benchmark_foldseek": "固定基准 · Foldseek",
        "alphafold_swissprot": "AlphaFold / Swiss-Prot",
        "pdb": "PDB 结构库",
        "benchmark_tmtools": "固定基准 · TM-align 预览",
        "min_tm": "最低 TM-score",
        "max_rmsd": "最高 RMSD (Å)",
        "min_qcov": "最低查询覆盖度",
        "min_tcov": "最低目标覆盖度",
        "min_alnlen": "最低对齐长度",
        "max_identity": "最高序列一致性",
        "current_run": "当前运行",
        "run_tab": "运行",
        "candidates_tab": "候选",
        "structure_tab": "结构",
        "history_tab": "历史",
        "doctor_tab": "诊断",
        "workbench": "运行工作台",
        "preflight": "预检",
        "ready": "可运行",
        "blocked": "存在堵点",
        "warnings": "警告",
        "run": "运行",
        "cancel": "取消",
        "reannotate": "重新注释",
        "open_folder": "打开结果目录",
        "query_metric": "查询",
        "database_metric": "数据库",
        "environment_metric": "环境",
        "candidates_metric": "候选数",
        "not_run": "尚未运行",
        "running": "运行中",
        "completed": "已完成",
        "failed": "失败",
        "cancelled": "已取消",
        "recent_events": "实时进度",
        "process_log": "运行日志",
        "latest_results": "最新结果",
        "no_results": "尚无可显示的候选结果。",
        "source_filter": "来源",
        "tier_filter": "WASE-screen 层级",
        "reviewed_filter": "Reviewed 状态",
        "evidence_filter": "最低证据完整度",
        "any": "全部",
        "reviewed": "Reviewed",
        "unreviewed": "Unreviewed",
        "unknown": "未知",
        "candidate_table": "候选表",
        "evidence_map": "结构证据图",
        "download_excel": "导出 Excel",
        "download_tsv": "导出 TSV",
        "download_fasta": "导出 FASTA",
        "select_candidate": "候选结构",
        "overlay": "查询与候选叠合",
        "alignment": "催化核心残基对齐",
        "query_legend": "8HFW 查询",
        "target_legend": "候选结构",
        "catalytic_legend": "催化残基",
        "target_chain": "目标链",
        "history": "运行历史",
        "run_details": "运行详情",
        "refresh": "刷新",
        "network": "检查网络",
        "rerun_doctor": "重新诊断",
        "diagnostic_results": "诊断结果",
        "download_json": "导出诊断 JSON",
        "no_structure": "未找到可加载的候选结构文件。",
        "start_error": "无法启动运行",
        "annotation_done": "重新注释已完成。",
        "status": "状态",
        "search_engine": "搜索引擎",
        "run_count": "运行数",
    },
    "en": {
        "settings": "Run settings",
        "language": "Language",
        "query": "Query structure",
        "chain": "Query chain",
        "database": "Search database",
        "benchmark_foldseek": "Fixed benchmark · Foldseek",
        "alphafold_swissprot": "AlphaFold / Swiss-Prot",
        "pdb": "PDB structure database",
        "benchmark_tmtools": "Fixed benchmark · TM-align preview",
        "min_tm": "Minimum TM-score",
        "max_rmsd": "Maximum RMSD (Å)",
        "min_qcov": "Minimum query coverage",
        "min_tcov": "Minimum target coverage",
        "min_alnlen": "Minimum alignment length",
        "max_identity": "Maximum sequence identity",
        "current_run": "Current run",
        "run_tab": "Run",
        "candidates_tab": "Candidates",
        "structure_tab": "Structure",
        "history_tab": "History",
        "doctor_tab": "Doctor",
        "workbench": "Run workbench",
        "preflight": "Preflight",
        "ready": "Ready",
        "blocked": "Blocked",
        "warnings": "Warnings",
        "run": "Run",
        "cancel": "Cancel",
        "reannotate": "Re-annotate",
        "open_folder": "Open results folder",
        "query_metric": "Query",
        "database_metric": "Database",
        "environment_metric": "Environment",
        "candidates_metric": "Candidates",
        "not_run": "Not run",
        "running": "Running",
        "completed": "Completed",
        "failed": "Failed",
        "cancelled": "Cancelled",
        "recent_events": "Live progress",
        "process_log": "Run log",
        "latest_results": "Latest results",
        "no_results": "No candidate results are available yet.",
        "source_filter": "Source",
        "tier_filter": "WASE-screen tier",
        "reviewed_filter": "Reviewed status",
        "evidence_filter": "Minimum evidence coverage",
        "any": "Any",
        "reviewed": "Reviewed",
        "unreviewed": "Unreviewed",
        "unknown": "Unknown",
        "candidate_table": "Candidate table",
        "evidence_map": "Structural evidence map",
        "download_excel": "Export Excel",
        "download_tsv": "Export TSV",
        "download_fasta": "Export FASTA",
        "select_candidate": "Candidate structure",
        "overlay": "Query-candidate overlay",
        "alignment": "Catalytic-core alignment",
        "query_legend": "8HFW query",
        "target_legend": "Candidate structure",
        "catalytic_legend": "Catalytic residues",
        "target_chain": "Target chain",
        "history": "Run history",
        "run_details": "Run details",
        "refresh": "Refresh",
        "network": "Check network",
        "rerun_doctor": "Run doctor",
        "diagnostic_results": "Diagnostic results",
        "download_json": "Export doctor JSON",
        "no_structure": "No loadable candidate structure was found.",
        "start_error": "Run could not start",
        "annotation_done": "Re-annotation completed.",
        "status": "Status",
        "search_engine": "Search engine",
        "run_count": "Runs",
    },
}

TEXT["zh"].update(
    {
        "settings": "运行设置",
        "language": "语言",
        "query": "查询结构",
        "chain": "查询链",
        "database": "搜索数据库",
        "benchmark_foldseek": "固定基准 · Foldseek",
        "alphafold_swissprot": "AlphaFold / Swiss-Prot",
        "pdb": "PDB 结构库",
        "benchmark_tmtools": "固定基准 · TM-align 预览",
        "min_tm": "最低 TM-score",
        "max_rmsd": "最高 RMSD (Å)",
        "min_qcov": "最低查询覆盖度",
        "min_tcov": "最低目标覆盖度",
        "min_alnlen": "最低对齐长度",
        "max_identity": "最高序列一致性",
        "current_run": "当前运行",
        "run_tab": "运行",
        "prediction_tab": "AF预测",
        "candidates_tab": "候选",
        "structure_tab": "结构",
        "history_tab": "历史",
        "doctor_tab": "诊断",
        "workbench": "运行工作台",
        "preflight": "预检",
        "ready": "可运行",
        "blocked": "存在堵点",
        "warnings": "警告",
        "run": "运行",
        "cancel": "取消",
        "reannotate": "重新注释",
        "open_folder": "打开结果目录",
        "query_metric": "查询",
        "database_metric": "数据库",
        "environment_metric": "环境",
        "candidates_metric": "候选数",
        "not_run": "尚未运行",
        "running": "运行中",
        "completed": "已完成",
        "failed": "失败",
        "cancelled": "已取消",
        "recent_events": "实时进度",
        "process_log": "运行日志",
        "latest_results": "最新结果",
        "no_results": "尚无可显示的候选结果。",
        "source_filter": "来源",
        "tier_filter": "WASE-screen 层级",
        "reviewed_filter": "Reviewed 状态",
        "evidence_filter": "最低证据完整度",
        "any": "全部",
        "reviewed": "Reviewed",
        "unreviewed": "Unreviewed",
        "unknown": "未知",
        "candidate_table": "候选表",
        "evidence_map": "结构证据图",
        "download_excel": "导出 Excel",
        "download_tsv": "导出 TSV",
        "download_fasta": "导出 FASTA",
        "select_candidate": "候选结构",
        "overlay": "查询与候选叠合",
        "alignment": "催化核心残基对齐",
        "query_legend": "查询结构",
        "target_legend": "候选结构",
        "catalytic_legend": "催化残基",
        "target_chain": "目标链",
        "history": "运行历史",
        "run_details": "运行详情",
        "refresh": "刷新",
        "network": "检查网络",
        "rerun_doctor": "重新诊断",
        "diagnostic_results": "诊断结果",
        "download_json": "导出诊断 JSON",
        "no_structure": "未找到可加载的候选结构文件。",
        "start_error": "无法启动运行",
        "annotation_done": "重新注释已完成。",
        "status": "状态",
        "search_engine": "搜索引擎",
        "run_count": "运行数",
        "af_workbench": "AlphaFold2 单体预测",
        "input_mode": "输入方式",
        "paste_fasta": "粘贴 FASTA",
        "upload_fasta": "上传 FASTA",
        "candidate_input": "已有候选",
        "fasta_text": "单条蛋白 FASTA",
        "preset": "质量档",
        "quick": "快速",
        "standard": "标准",
        "high": "高采样",
        "use_templates": "启用模板搜索",
        "privacy_confirm": "我确认该序列会在本次任务中发送到公共 ColabFold MSA 服务。",
        "privacy_required": "未命中本地 MSA 缓存时，必须逐次确认第三方序列传输。",
        "map_preview": "8HFW 催化位点映射",
        "manual_sites": "手动修正位点",
        "predict_only": "仅预测",
        "predict_mine": "预测并挖掘",
        "predictions": "预测记录",
        "select_prediction": "选择预测",
        "model_ranking": "模型排名",
        "mean_plddt": "平均 pLDDT",
        "ptm": "pTM",
        "pae": "PAE 热图",
        "plddt_curve": "pLDDT 曲线",
        "predicted_structure": "按 pLDDT 着色的结构",
        "no_prediction": "尚无预测结果。",
        "model_count": "模型数",
        "retrigger": "重新预测",
        "accession": "候选 accession",
        "service_mode": "后台常驻 · 自动恢复",
        "online_mode": "HTTPS 在线演示 · 密码保护",
    }
)
TEXT["en"].update(
    {
        "benchmark_foldseek": "Fixed benchmark · Foldseek",
        "benchmark_tmtools": "Fixed benchmark · TM-align preview",
        "max_rmsd": "Maximum RMSD (Å)",
        "prediction_tab": "AF Prediction",
        "af_workbench": "AlphaFold2 monomer prediction",
        "input_mode": "Input source",
        "paste_fasta": "Paste FASTA",
        "upload_fasta": "Upload FASTA",
        "candidate_input": "Existing candidate",
        "fasta_text": "Single-protein FASTA",
        "preset": "Quality preset",
        "quick": "Quick",
        "standard": "Standard",
        "high": "High sampling",
        "use_templates": "Enable template search",
        "privacy_confirm": "I confirm that this sequence will be sent to the public ColabFold MSA service for this job.",
        "privacy_required": "Third-party sequence transfer requires confirmation for every MSA cache miss.",
        "map_preview": "8HFW catalytic-site mapping",
        "manual_sites": "Manual site correction",
        "predict_only": "Predict",
        "predict_mine": "Predict and mine",
        "predictions": "Prediction history",
        "select_prediction": "Select prediction",
        "model_ranking": "Model ranking",
        "mean_plddt": "Mean pLDDT",
        "ptm": "pTM",
        "pae": "PAE heatmap",
        "plddt_curve": "pLDDT profile",
        "predicted_structure": "Structure colored by pLDDT",
        "no_prediction": "No prediction results are available yet.",
        "model_count": "Models",
        "retrigger": "Predict again",
        "accession": "Candidate accession",
        "service_mode": "Persistent local service · auto-recovery",
        "online_mode": "HTTPS expert demo · password protected",
    }
)
TEXT["zh"].update(
    {
        "nr_tab": "NR挖掘",
        "nr_workbench": "NCBI NR 序列扩增",
        "nr_scope": "数据库范围",
        "nr_local_scope": "本机真实命中子集",
        "nr_full_scope": "完整 NR 服务器作业",
        "nr_source_run": "结构种子运行",
        "nr_run_history": "NR 运行记录",
        "nr_build_snapshot": "构建本机 NR 子集",
        "nr_expand_local": "运行本机子集",
        "nr_export_job": "导出服务器作业",
        "nr_import_job": "导入服务器结果",
        "nr_contact_email": "NCBI 联系邮箱",
        "nr_consent": "我确认本次会将 Bur01 序列发送至 NCBI BLAST 公共服务",
        "nr_snapshot_missing": "尚未安装本机 NR 命中子集。先填写邮箱并逐次确认。",
        "nr_server_note": "完整 NR 仅在服务器/HPC 上运行；导回且通过校验前不会发布 v0.4.0。",
        "nr_sequence_evidence": "序列预筛证据",
        "nr_wase_notice": "未经 AF2/Foldseek/TM-align 复核，WASE-screen 保持未评估。",
        "nr_af2_ready": "AF2-ready",
        "nr_representatives": "代表序列",
        "nr_send_af2": "送入 AF2 复核",
        "nr_candidate_input": "NR 候选",
        "nr_result_zip": "服务器结果 ZIP",
        "nr_min_qcov": "最低查询覆盖度",
        "nr_triad": "催化三联体",
        "nr_accession": "NR accession",
        "nr_process": "NR 后台任务",
        "nr_no_runs": "尚无 NR 运行记录。",
        "nr_snapshot": "本机快照",
        "nr_integrity": "可追溯性",
        "nr_not_structural": "仅序列证据",
        "nr_pipeline_trace": "NR 处理轨迹",
    }
)
TEXT["en"].update(
    {
        "nr_tab": "NR Mining",
        "nr_workbench": "NCBI NR sequence expansion",
        "nr_scope": "Database scope",
        "nr_local_scope": "Real local hit subset",
        "nr_full_scope": "Complete-NR server job",
        "nr_source_run": "Structure seed run",
        "nr_run_history": "NR run history",
        "nr_build_snapshot": "Build local NR subset",
        "nr_expand_local": "Run local subset",
        "nr_export_job": "Export server job",
        "nr_import_job": "Import server result",
        "nr_contact_email": "NCBI contact email",
        "nr_consent": "I confirm that the Bur01 sequence may be sent to the public NCBI BLAST service for this request",
        "nr_snapshot_missing": "No local NR hit subset is installed. Enter an email and confirm each transfer first.",
        "nr_server_note": "Complete NR runs only on a server/HPC; v0.4.0 is not released before verified import.",
        "nr_sequence_evidence": "Sequence-screen evidence",
        "nr_wase_notice": "WASE-screen remains unevaluated until AF2/Foldseek/TM-align validation.",
        "nr_af2_ready": "AF2-ready",
        "nr_representatives": "Representatives",
        "nr_send_af2": "Send to AF2 validation",
        "nr_candidate_input": "NR candidate",
        "nr_result_zip": "Server result ZIP",
        "nr_min_qcov": "Minimum query coverage",
        "nr_triad": "Catalytic triad",
        "nr_accession": "NR accession",
        "nr_process": "NR background task",
        "nr_no_runs": "No NR run is available yet.",
        "nr_snapshot": "Local snapshot",
        "nr_integrity": "Traceability",
        "nr_not_structural": "Sequence evidence only",
        "nr_pipeline_trace": "NR processing trace",
    }
)


st.set_page_config(page_title="AquaEsterNet-StructMine", layout="wide")
st.markdown(
    f"""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;500;600&display=swap');
    :root {{ --bg:{COLORS['background']}; --ink:{COLORS['ink']}; --green:{COLORS['green']}; --amber:{COLORS['amber']}; --red:{COLORS['red']}; --line:{COLORS['line']}; --muted:{COLORS['muted']}; }}
    html, body, .stApp, input, textarea, button {{ font-family: 'IBM Plex Sans', 'Segoe UI', sans-serif; color: var(--ink); letter-spacing: 0; }}
    [data-testid="stIconMaterial"], [data-testid="stAlertDynamicIcon"] {{ font-family: 'Material Symbols Rounded' !important; }}
    [data-testid="stAppDeployButton"], #MainMenu {{ display:none !important; }}
    .stApp {{ background: var(--bg); }}
    [data-testid="stSidebar"] {{ background: #EEF3F1; border-right: 1px solid var(--line); }}
    [data-testid="stSidebar"] h1 {{ font-size: 1rem; line-height: 1.2; }}
    h1 {{ font-size: 1.75rem !important; line-height: 1.2 !important; font-weight: 600 !important; }}
    h2 {{ font-size: 1.2rem !important; line-height: 1.3 !important; font-weight: 600 !important; }}
    h3 {{ font-size: 1rem !important; line-height: 1.3 !important; font-weight: 600 !important; }}
    button, [data-testid="stBaseButton-primary"], [data-testid="stBaseButton-secondary"] {{ border-radius: 4px !important; }}
    [data-testid="stBaseButton-primary"] {{ background: var(--green) !important; border-color: var(--green) !important; }}
    [data-testid="stBaseButton-primary"]:disabled {{ background:#7FA99F !important; color:#FFFFFF !important; border-color:#7FA99F !important; opacity:.72; }}
    [data-testid="stMetric"] {{ padding: .55rem 0; border-top: 1px solid var(--line); border-bottom: 1px solid var(--line); }}
    [data-testid="stMetricValue"] {{ font-family: 'IBM Plex Mono', Consolas, monospace; font-size: 1.15rem; }}
    [data-baseweb="tab-list"] {{ gap: 1.25rem; border-bottom: 1px solid var(--line); }}
    [data-baseweb="tab"] {{ padding-left: 0; padding-right: 0; }}
    [data-baseweb="tab-highlight"] {{ background-color: var(--green); }}
    [data-testid="stDataFrame"] {{ border: 1px solid var(--line); border-radius: 4px; overflow: hidden; }}
    .status-line {{ display:flex; gap:.7rem; align-items:center; font-size:.92rem; color:var(--muted); margin:.2rem 0 .65rem; }}
    .status-dot {{ width:.55rem; height:.55rem; border-radius:50%; display:inline-block; flex:0 0 auto; }}
    .legend-row {{ display:flex; flex-wrap:wrap; gap:1rem; margin:.25rem 0 .7rem; color:var(--muted); font-size:.86rem; }}
    .swatch {{ display:inline-block; width:.8rem; height:.8rem; margin-right:.35rem; vertical-align:-.08rem; border:1px solid rgba(0,0,0,.12); }}
    .alignment-scroll {{ overflow-x:auto; border-top:1px solid var(--line); border-bottom:1px solid var(--line); padding:.8rem 0; background:#FFFFFF; }}
    .alignment-row {{ font-family:'IBM Plex Mono', Consolas, monospace; white-space:pre; line-height:1.65; font-size:.78rem; min-width:max-content; }}
    .alignment-label {{ display:inline-block; width:5.5rem; color:var(--muted); font-family:'IBM Plex Sans', sans-serif; }}
    .cat-residue {{ color:#FFFFFF; background:var(--red); font-weight:600; padding:1px 0; }}
    .muted {{ color:var(--muted); }}
    @media (max-width: 640px) {{
      .block-container {{ padding: 1rem .75rem 3rem !important; }}
      h1 {{ font-size:1.35rem !important; }}
      [data-baseweb="tab-list"] {{ overflow-x:auto; gap:1rem; }}
      [data-testid="stHorizontalBlock"] {{ gap:.45rem; }}
      .alignment-label {{ width:4.2rem; }}
      .alignment-row {{ font-size:.7rem; }}
    }}
    </style>
    """,
    unsafe_allow_html=True,
)

remote_access = load_remote_access(REMOTE_ACCESS_PATH)
if remote_access is not None and not st.session_state.get("remote_authenticated", False):
    st.title("AquaEsterNet-StructMine")
    st.subheader("在线演示访问 / Secure remote access")
    st.caption("请输入共享访问密码。模型计算仍在所有者电脑上运行，序列和结果不会托管到本页面服务器。")
    with st.form("remote_access_login"):
        remote_password = st.text_input("访问密码 / Access password", type="password", max_chars=256)
        login_submitted = st.form_submit_button("进入工作台 / Open workbench", type="primary")
    if login_submitted:
        if verify_remote_password(remote_password, remote_access):
            st.session_state.remote_authenticated = True
            st.session_state.remote_auth_failures = 0
            st.rerun()
        else:
            st.session_state.remote_auth_failures = st.session_state.get("remote_auth_failures", 0) + 1
            time.sleep(min(3.0, float(st.session_state.remote_auth_failures)))
            st.error("访问密码不正确 / Incorrect access password")
    st.stop()


def tx(key: str) -> str:
    return TEXT[st.session_state.get("language", "zh")][key]


def doctor_message(check: dict[str, Any]) -> str:
    if st.session_state.get("language", "zh") != "zh":
        return check["message"]
    name = check.get("name", "")
    translations = {
        "windows_feature_Microsoft-Windows-Subsystem-Linux": "Windows 的 WSL 功能尚未启用",
        "windows_feature_VirtualMachinePlatform": "Windows 的虚拟机平台功能尚未启用",
        "wsl_distribution": "尚未安装 Ubuntu-Foldseek 发行版",
        "foldseek": "Foldseek 当前不可用",
        "database_benchmark_v01": "基准数据库尚未创建；Foldseek 就绪后运行时会自动创建",
        "colabfold": "LocalColabFold 尚未就绪",
        "jax_cuda": "JAX 未检测到 CUDA GPU",
    }
    if name in translations:
        return translations[name]
    if name.startswith("database_") and check.get("status") == "fail":
        path = check.get("details", {}).get("path", "")
        return f"目标数据库尚未就绪：{path}"
    return check["message"]


def _clean_value(value: Any) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"nan", "none"} else text


def scan_runs() -> list[dict[str, Any]]:
    runs: list[dict[str, Any]] = []
    if not RUN_ROOT.exists():
        return runs
    for manifest_path in RUN_ROOT.glob("*/run_manifest.json"):
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        manifest["_manifest_path"] = str(manifest_path)
        manifest["_run_dir"] = str(manifest_path.parent)
        runs.append(manifest)
    return sorted(runs, key=lambda item: item.get("created_at", ""), reverse=True)


def scan_predictions() -> list[dict[str, Any]]:
    predictions: list[dict[str, Any]] = []
    if not PREDICTION_ROOT.exists():
        return predictions
    for manifest_path in PREDICTION_ROOT.glob("*/prediction_manifest.json"):
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        manifest["_manifest_path"] = str(manifest_path)
        manifest["_prediction_dir"] = str(manifest_path.parent)
        predictions.append(manifest)
    return sorted(predictions, key=lambda item: item.get("created_at", ""), reverse=True)


def scan_nr_runs() -> list[dict[str, Any]]:
    nr_runs: list[dict[str, Any]] = []
    if not NR_RUN_ROOT.exists():
        return nr_runs
    for manifest_path in NR_RUN_ROOT.glob("*/nr_run_manifest.json"):
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        manifest["_manifest_path"] = str(manifest_path)
        manifest["_nr_run_dir"] = str(manifest_path.parent)
        nr_runs.append(manifest)
    return sorted(nr_runs, key=lambda item: item.get("created_at", ""), reverse=True)


def read_nr_candidates(manifest: dict[str, Any] | None) -> pd.DataFrame:
    if not manifest:
        return pd.DataFrame()
    path_text = _clean_value(manifest.get("artifacts", {}).get("candidates_tsv", ""))
    path = Path(path_text) if path_text else Path("__missing_nr_candidates__")
    if not path.is_file():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        return pd.DataFrame()


def read_nr_events(manifest: dict[str, Any] | None) -> list[dict[str, Any]]:
    if not manifest:
        return []
    run_dir_text = _clean_value(manifest.get("_nr_run_dir", ""))
    run_dir = Path(run_dir_text) if run_dir_text else None
    path = run_dir / "events.jsonl" if run_dir else None
    if path is None or not path.is_file():
        return []
    events: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        try:
            events.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return events


def read_prediction_events(prediction_id: str | None) -> list[dict[str, Any]]:
    if not prediction_id:
        return []
    path = PREDICTION_ROOT / prediction_id / "events.jsonl"
    if not path.is_file():
        return []
    events: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        try:
            events.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return events


def read_candidates(manifest: dict[str, Any] | None) -> pd.DataFrame:
    if not manifest:
        return pd.DataFrame()
    artifacts = manifest.get("artifacts", {})
    path = Path(artifacts.get("annotated_tsv", ""))
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t")
    except (OSError, pd.errors.ParserError, pd.errors.EmptyDataError):
        return pd.DataFrame()


def read_events(run_id: str | None) -> list[dict[str, Any]]:
    if not run_id:
        return []
    path = RUN_ROOT / run_id / "events.jsonl"
    if not path.exists():
        return []
    events: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        try:
            events.append(json.loads(line))
        except json.JSONDecodeError:
            continue
    return events


@st.cache_data(ttl=30, show_spinner=False)
def cached_doctor(config_path: str, modified: float, network: bool) -> dict[str, Any]:
    del modified
    return run_doctor(config_path, check_network=network).to_dict()


def doctor_for(config_path: str, network: bool = False) -> dict[str, Any]:
    if os.environ.get("AQUA_TEST_MODE") == "1":
        checks = [
            {
                "name": name,
                "category": category,
                "status": "pass",
                "message": "isolated UI test",
                "details": {},
            }
            for name, category in (
                ("configuration", "project"),
                ("wsl_gpu", "prediction"),
                ("colabfold", "prediction"),
                ("jax_cuda", "prediction"),
                ("mmseqs2", "nr"),
                ("nr_demo_snapshot", "nr"),
            )
        ]
        return {
            "project_version": "0.4.0",
            "platform": "isolated-test",
            "passed": True,
            "summary": {"pass": len(checks), "warn": 0, "fail": 0},
            "checks": checks,
        }
    path = ROOT / config_path
    return cached_doctor(config_path, path.stat().st_mtime, network)


def write_runtime_config(
    base_config: str,
    query_path: str,
    chain: str,
    thresholds: dict[str, Any],
    run_id: str,
) -> Path:
    cfg = load_config(base_config)
    cfg["query"]["pdb_file"] = query_path
    cfg["query"]["chain"] = chain
    cfg["filter"].update(thresholds)
    cfg["output"]["run_root"] = str(RUN_ROOT)
    if base_config == "configs/benchmark_foldseek.yaml":
        cfg["search_databases"][0]["auto_prepare"] = True
    ProjectConfig.model_validate(cfg)
    UI_CONFIG_ROOT.mkdir(parents=True, exist_ok=True)
    path = UI_CONFIG_ROOT / f"{run_id}.yaml"
    tmp = path.with_suffix(".yaml.tmp")
    tmp.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True), encoding="utf-8")
    tmp.replace(path)
    return path


def start_run(config_path: Path, run_id: str) -> None:
    run_dir = RUN_ROOT / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    log_path = run_dir / "ui_process.log"
    log_handle = log_path.open("w", encoding="utf-8")
    command = [
        sys.executable,
        "-m",
        "bur01_structmine.cli",
        "run",
        "--config",
        str(config_path),
        "--run-id",
        run_id,
    ]
    creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    process = subprocess.Popen(
        command,
        cwd=ROOT,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        creationflags=creationflags,
    )
    st.session_state.ui_process = process
    st.session_state.ui_log_handle = log_handle
    st.session_state.ui_run_id = run_id
    st.session_state.ui_log_path = str(log_path)
    st.session_state.ui_fragment_state = "running"


def process_state() -> tuple[str, int | None]:
    process = st.session_state.get("ui_process")
    if process is None:
        return "idle", None
    code = process.poll()
    if code is None:
        return "running", None
    handle = st.session_state.get("ui_log_handle")
    if handle and not handle.closed:
        handle.close()
    run_id = st.session_state.get("ui_run_id")
    manifest_path = RUN_ROOT / run_id / "run_manifest.json" if run_id else None
    if manifest_path and manifest_path.exists():
        try:
            manifest_status = json.loads(manifest_path.read_text(encoding="utf-8")).get("status")
            if manifest_status == "cancelled":
                return "cancelled", code
        except (OSError, json.JSONDecodeError):
            pass
    return ("completed" if code == 0 else "failed"), code


def cancel_run() -> None:
    process = st.session_state.get("ui_process")
    if process is None or process.poll() is not None:
        return
    process.terminate()
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        process.kill()
    handle = st.session_state.get("ui_log_handle")
    if handle and not handle.closed:
        handle.close()
    run_id = st.session_state.get("ui_run_id")
    if run_id:
        mark_run_cancelled(RUN_ROOT / run_id)


def start_prediction_process(
    prediction_id: str,
    preset: str,
    use_templates: bool,
    consent: bool,
    catalytic_sites: dict[str, int],
    fasta_path: Path | None = None,
    candidate_run_id: str | None = None,
    candidate_nr_run_id: str | None = None,
    accession: str | None = None,
    mine_after: bool = False,
) -> None:
    log_root = ROOT / "results" / "ui_process_logs"
    log_root.mkdir(parents=True, exist_ok=True)
    log_path = log_root / f"{prediction_id}.log"
    log_handle = log_path.open("w", encoding="utf-8")
    command = [
        sys.executable,
        "-m",
        "bur01_structmine.cli",
        "predict-and-mine" if mine_after else "predict",
    ]
    if fasta_path is not None:
        command.extend(["--fasta", str(fasta_path)])
    elif candidate_nr_run_id is not None:
        command.extend(["--nr-run", str(candidate_nr_run_id), "--accession", str(accession)])
    else:
        command.extend(["--candidate-run", str(candidate_run_id), "--accession", str(accession)])
    command.extend(
        [
            "--preset",
            preset,
            "--prediction-id",
            prediction_id,
            "--prediction-output-root",
            str(PREDICTION_ROOT),
        ]
    )
    if mine_after:
        command.extend(["--mining-config", "configs/alphafold_swissprot.yaml"])
    if use_templates:
        command.append("--templates")
    if consent:
        command.append("--confirm-public-msa")
    for label, position in catalytic_sites.items():
        command.extend(["--catalytic-site", f"{label}={position}"])
    creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    process = subprocess.Popen(
        command,
        cwd=ROOT,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        creationflags=creationflags,
    )
    st.session_state.prediction_process = process
    st.session_state.prediction_log_handle = log_handle
    st.session_state.prediction_process_log = str(log_path)
    st.session_state.active_prediction_id = prediction_id
    st.session_state.prediction_fragment_state = "running"


def prediction_process_state() -> tuple[str, int | None]:
    process = st.session_state.get("prediction_process")
    if process is None:
        return "idle", None
    code = process.poll()
    if code is None:
        return "running", None
    handle = st.session_state.get("prediction_log_handle")
    if handle and not handle.closed:
        handle.close()
    prediction_id = st.session_state.get("active_prediction_id")
    manifest_path = PREDICTION_ROOT / prediction_id / "prediction_manifest.json" if prediction_id else None
    if manifest_path and manifest_path.is_file():
        try:
            status = json.loads(manifest_path.read_text(encoding="utf-8")).get("status", "")
            if status in {"completed", "failed", "cancelled"}:
                return status, code
        except (OSError, json.JSONDecodeError):
            pass
    return ("completed" if code == 0 else "failed"), code


def cancel_prediction_process() -> None:
    process = st.session_state.get("prediction_process")
    prediction_id = st.session_state.get("active_prediction_id")
    if process is None or process.poll() is not None:
        return
    prediction_dir = PREDICTION_ROOT / prediction_id if prediction_id else None
    if prediction_dir and prediction_dir.exists():
        request_prediction_cancel(prediction_dir)
        try:
            process.wait(timeout=8)
        except subprocess.TimeoutExpired:
            process.terminate()
    else:
        process.terminate()
    handle = st.session_state.get("prediction_log_handle")
    if handle and not handle.closed:
        handle.close()


def start_nr_process(task_id: str, command: list[str]) -> None:
    log_root = ROOT / "results" / "ui_process_logs"
    log_root.mkdir(parents=True, exist_ok=True)
    log_path = log_root / f"{task_id}.log"
    log_handle = log_path.open("w", encoding="utf-8")
    creationflags = subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
    process = subprocess.Popen(
        [sys.executable, "-m", "bur01_structmine.cli", *command],
        cwd=ROOT,
        stdout=log_handle,
        stderr=subprocess.STDOUT,
        creationflags=creationflags,
    )
    st.session_state.nr_process = process
    st.session_state.nr_log_handle = log_handle
    st.session_state.nr_process_log = str(log_path)
    st.session_state.nr_task_id = task_id


def nr_process_state() -> tuple[str, int | None]:
    process = st.session_state.get("nr_process")
    if process is None:
        return "idle", None
    code = process.poll()
    if code is None:
        return "running", None
    handle = st.session_state.get("nr_log_handle")
    if handle and not handle.closed:
        handle.close()
    return ("completed" if code == 0 else "failed"), code


def cancel_nr_process() -> None:
    process = st.session_state.get("nr_process")
    if process is None or process.poll() is not None:
        return
    if os.name == "nt":
        subprocess.run(
            ["taskkill.exe", "/PID", str(process.pid), "/T", "/F"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            creationflags=subprocess.CREATE_NO_WINDOW,
            check=False,
        )
    else:
        process.terminate()
    try:
        process.wait(timeout=8)
    except subprocess.TimeoutExpired:
        process.kill()
    handle = st.session_state.get("nr_log_handle")
    if handle and not handle.closed:
        handle.close()


@st.fragment(run_every=FRAGMENT_RUN_EVERY)
def render_nr_progress_fragment() -> None:
    nr_state, nr_exit_code = nr_process_state()
    if nr_state == "idle":
        return
    color = COLORS["amber"] if nr_state == "running" else (
        COLORS["green"] if nr_state == "completed" else COLORS["red"]
    )
    st.markdown(
        (
            '<div class="status-line">'
            f'<span class="status-dot" style="background:{color}"></span>'
            f'<strong>{tx("nr_process")}</strong><span>{tx(nr_state)}</span>'
            f'<span>exit={nr_exit_code if nr_exit_code is not None else "..."}</span>'
            "</div>"
        ),
        unsafe_allow_html=True,
    )
    log_path_text = st.session_state.get("nr_process_log", "")
    log_path = Path(log_path_text) if log_path_text else None
    if log_path and log_path.is_file():
        lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
        st.code("\n".join(lines[-18:]) or tx("running"), language="text")


@st.fragment(run_every=FRAGMENT_RUN_EVERY)
def render_progress_fragment() -> None:
    fragment_state, fragment_exit_code = process_state()
    active_run_id = st.session_state.get("ui_run_id")
    events = read_events(active_run_id)
    if fragment_state == "idle" and not events:
        return
    latest_event = events[-1] if events else {}
    progress = float(latest_event.get("progress") or 0.0)
    state_key = fragment_state if fragment_state in TEXT[st.session_state.language] else "running"
    st.progress(progress, text=latest_event.get("message") or tx(state_key))
    st.subheader(tx("recent_events"))
    event_rows = [
        {
            "time": event.get("timestamp", "")[11:19],
            "stage": event.get("stage", ""),
            "state": event.get("state", ""),
            "message": event.get("message", ""),
        }
        for event in events[-8:]
    ]
    if event_rows:
        st.dataframe(
            pd.DataFrame(event_rows),
            hide_index=True,
            width="stretch",
            height=min(310, 38 * (len(event_rows) + 1)),
        )
    log_text_path = _clean_value(st.session_state.get("ui_log_path", ""))
    log_path = Path(log_text_path) if log_text_path else None
    if log_path and log_path.is_file():
        log_text = log_path.read_text(encoding="utf-8", errors="replace")
        with st.expander(tx("process_log"), expanded=fragment_state == "failed"):
            st.code("\n".join(log_text.splitlines()[-80:]), language="text")
    if fragment_state == "failed":
        st.error(f"{tx('failed')} · exit code {fragment_exit_code}")
    previous_state = st.session_state.get("ui_fragment_state")
    st.session_state.ui_fragment_state = fragment_state
    if previous_state == "running" and fragment_state != "running":
        st.rerun(scope="app")


@st.fragment(run_every=FRAGMENT_RUN_EVERY)
def render_prediction_progress_fragment() -> None:
    fragment_state, exit_code = prediction_process_state()
    prediction_id = st.session_state.get("active_prediction_id")
    events = read_prediction_events(prediction_id)
    if fragment_state == "idle" and not events:
        return
    latest = events[-1] if events else {}
    progress = float(latest.get("progress") or 0.0)
    st.progress(progress, text=latest.get("message") or tx("running"))
    if events:
        rows = [
            {
                "time": event.get("timestamp", "")[11:19],
                "stage": event.get("stage", ""),
                "state": event.get("state", ""),
                "message": event.get("message", ""),
            }
            for event in events[-8:]
        ]
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch", height=min(310, 38 * (len(rows) + 1)))
    log_text = _clean_value(st.session_state.get("prediction_process_log", ""))
    log_path = Path(log_text) if log_text else None
    if log_path and log_path.is_file():
        with st.expander(tx("process_log"), expanded=fragment_state == "failed"):
            st.code("\n".join(log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-80:]), language="text")
    if fragment_state == "failed":
        st.error(f"{tx('failed')} · exit code {exit_code}")
    previous = st.session_state.get("prediction_fragment_state")
    st.session_state.prediction_fragment_state = fragment_state
    if previous == "running" and fragment_state != "running":
        st.rerun(scope="app")


def open_folder(path: str | Path) -> None:
    resolved = Path(path).resolve()
    if os.name == "nt":
        os.startfile(resolved)  # type: ignore[attr-defined]
    elif sys.platform == "darwin":
        subprocess.Popen(["open", str(resolved)])
    else:
        subprocess.Popen(["xdg-open", str(resolved)])


def candidate_structure_path(row: pd.Series) -> Path | None:
    local = Path(_clean_value(row.get("target_structure_path"))) if _clean_value(row.get("target_structure_path")) else None
    if local and local.exists():
        return local
    target = _clean_value(row.get("target"))
    for path in (ROOT / "data" / "benchmark" / "structures").glob("*.pdb"):
        if any(label in target.upper() and label in path.name.upper() for label in ("8HFW", "1XKL", "6TJ2", "Q6RYA0", "Q9FVW3", "P60174", "P00520")):
            return path
    url = _clean_value(row.get("alphafold_pdb_url"))
    if not url:
        pdb_id = _clean_value(row.get("pdb_id"))
        if not pdb_id:
            pdb_id, _ = parse_pdb_chain(target)
        if pdb_id:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    if not url:
        return None
    UI_STRUCTURE_CACHE.mkdir(parents=True, exist_ok=True)
    safe_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", target) or "candidate"
    path = UI_STRUCTURE_CACHE / f"{safe_name}.pdb"
    if not path.exists():
        session = CachedSession(ROOT / ".cache" / "ui_structure_http", timeout=60)
        content = session.get_text(url)
        if not content or "ATOM" not in content:
            return None
        tmp = path.with_suffix(".pdb.tmp")
        tmp.write_text(content, encoding="utf-8")
        tmp.replace(path)
    return path


@st.cache_data(show_spinner=False)
def cached_overlay(
    query: str,
    query_mtime: int,
    target: str,
    target_mtime: int,
    query_chain: str,
    target_chain: str,
    catalytic_query_residues: tuple[tuple[int, str], ...],
) -> StructureOverlay:
    del query_mtime, target_mtime
    return align_structure_overlay(
        query,
        target,
        query_chain=query_chain,
        target_chain=target_chain or None,
        catalytic_query_residues=catalytic_query_residues,
    )


def render_overlay(overlay: StructureOverlay, catalytic_numbers: tuple[int, ...]) -> None:
    viewer = py3Dmol.view(width="100%", height=520)
    viewer.addModel(overlay.query_pdb, "pdb")
    viewer.addModel(overlay.target_pdb, "pdb")
    viewer.setStyle({"model": 0}, {"cartoon": {"color": COLORS["green"], "opacity": 0.88}})
    viewer.setStyle({"model": 1}, {"cartoon": {"color": COLORS["amber"], "opacity": 0.74}})
    viewer.addStyle({"model": 0, "resi": list(catalytic_numbers)}, {"stick": {"color": COLORS["red"], "radius": 0.22}, "sphere": {"color": COLORS["red"], "scale": 0.28}})
    if overlay.mapped_target_residue_numbers:
        viewer.addStyle(
            {"model": 1, "resi": list(overlay.mapped_target_residue_numbers)},
            {"stick": {"color": COLORS["red"], "radius": 0.18}},
        )
    viewer.setBackgroundColor(COLORS["background"])
    viewer.zoomTo()
    viewer.spin("y", 0.15)
    components.html(viewer._make_html(), height=540, scrolling=False)


def alignment_html(overlay: StructureOverlay, catalytic_numbers: tuple[int, ...]) -> str:
    catalytic_number_set = set(catalytic_numbers)
    query_index = -1
    highlighted_columns: set[int] = set()
    query_chars: list[str] = []
    for column, amino_acid in enumerate(overlay.aligned_query):
        residue_number = None
        if amino_acid != "-":
            query_index += 1
            if query_index < len(overlay.query_residue_numbers):
                residue_number = overlay.query_residue_numbers[query_index]
        escaped = html.escape(amino_acid)
        if residue_number in catalytic_number_set:
            highlighted_columns.add(column)
            query_chars.append(f'<span class="cat-residue" title="{amino_acid}{residue_number}">{escaped}</span>')
        else:
            query_chars.append(escaped)
    match_chars = [f'<span class="cat-residue">{html.escape(char)}</span>' if index in highlighted_columns else html.escape(char) for index, char in enumerate(overlay.alignment_match_line)]
    target_chars = [f'<span class="cat-residue">{html.escape(char)}</span>' if index in highlighted_columns else html.escape(char) for index, char in enumerate(overlay.aligned_target)]
    return (
        '<div class="alignment-scroll">'
        f'<div class="alignment-row"><span class="alignment-label">8HFW:A</span>{"".join(query_chars)}</div>'
        f'<div class="alignment-row"><span class="alignment-label"></span>{"".join(match_chars)}</div>'
        f'<div class="alignment-row"><span class="alignment-label">Target:{html.escape(overlay.target_chain)}</span>{"".join(target_chars)}</div>'
        '</div>'
    )


def candidate_accession(row: pd.Series) -> str:
    for column in ("target_uniprot_accession", "uniprot_accessions", "ncbi_protein_accession", "pdb_id", "target"):
        value = _clean_value(row.get(column))
        if value:
            return value.split(";")[0].strip()
    return ""


def render_prediction_structure(path: Path) -> None:
    viewer = py3Dmol.view(width="100%", height=520)
    viewer.addModel(path.read_text(encoding="utf-8", errors="replace"), "pdb")
    viewer.setStyle(
        {"model": 0},
        {
            "cartoon": {
                "colorscheme": {"prop": "b", "gradient": "roygb", "min": 50, "max": 90}
            }
        },
    )
    viewer.setBackgroundColor(COLORS["background"])
    viewer.zoomTo()
    viewer.spin("y", 0.12)
    components.html(viewer._make_html(), height=540, scrolling=False)


if "language" not in st.session_state:
    st.session_state.language = "zh"

st.sidebar.title("AquaEsterNet-StructMine")
language_label = st.sidebar.segmented_control(
    "Language",
    options=("中文", "English"),
    default="中文" if st.session_state.language == "zh" else "English",
    label_visibility="collapsed",
    width="stretch",
)
st.session_state.language = "zh" if language_label == "中文" else "en"
if remote_access is not None and st.sidebar.button("退出在线会话 / Sign out", width="stretch"):
    st.session_state.remote_authenticated = False
    st.rerun()

st.sidebar.header(tx("settings"))
query_path = st.sidebar.text_input(tx("query"), value="data/query/Bur01m.pdb")
query_chain = st.sidebar.text_input(tx("chain"), value="A", max_chars=4)
config_key = st.sidebar.selectbox(
    tx("database"),
    options=list(CONFIG_OPTIONS),
    format_func=lambda key: tx(key),
)
base_config = CONFIG_OPTIONS[config_key]
base_cfg = load_config(base_config)
base_filter = base_cfg["filter"]
min_tm = st.sidebar.number_input(tx("min_tm"), min_value=0.0, max_value=1.0, value=float(base_filter["min_alntmscore"]), step=0.05)
max_rmsd = st.sidebar.number_input(tx("max_rmsd"), min_value=0.1, max_value=20.0, value=float(base_filter["max_rmsd"]), step=0.5)
min_qcov = st.sidebar.number_input(tx("min_qcov"), min_value=0.0, max_value=1.0, value=float(base_filter["min_qcov"]), step=0.05)
min_tcov = st.sidebar.number_input(tx("min_tcov"), min_value=0.0, max_value=1.0, value=float(base_filter["min_tcov"]), step=0.05)
min_alnlen = st.sidebar.number_input(tx("min_alnlen"), min_value=1, max_value=5000, value=int(base_filter["min_alnlen"]), step=10)
max_identity = st.sidebar.number_input(tx("max_identity"), min_value=0.0, max_value=1.0, value=float(base_filter["max_sequence_identity"]), step=0.05)

runs = scan_runs()
run_by_id = {run["run_id"]: run for run in runs}
run_options = list(run_by_id)
selected_run_id = st.sidebar.selectbox(
    tx("current_run"),
    options=run_options,
    index=0 if run_options else None,
    format_func=lambda run_id: f"{run_by_id[run_id].get('status', '?')} · {run_id}",
    disabled=not run_options,
    placeholder=tx("not_run"),
)
selected_manifest = run_by_id.get(selected_run_id) if selected_run_id else None
candidates = read_candidates(selected_manifest)
predictions = scan_predictions()
prediction_by_id = {prediction["prediction_id"]: prediction for prediction in predictions}
nr_runs = scan_nr_runs()
nr_run_by_id = {run["nr_run_id"]: run for run in nr_runs}

preflight = doctor_for(base_config, network=False)
preflight_fails = [check for check in preflight["checks"] if check["status"] == "fail"]
preflight_warnings = [check for check in preflight["checks"] if check["status"] == "warn"]
blocking_fails = [
    check
    for check in preflight_fails
    if not (config_key == "benchmark_foldseek" and check["name"] == "database_benchmark_v01")
]
state, exit_code = process_state()
prediction_state, prediction_exit_code = prediction_process_state()
nr_state, nr_exit_code = nr_process_state()

st.title("AquaEsterNet-StructMine")
service_label = f"  |  {tx('service_mode')}" if os.environ.get("AQUA_UI_SERVICE") == "1" else ""
if remote_access is not None:
    service_label += f"  |  {tx('online_mode')}"
st.caption(f"v0.4.0  |  Structure + NR + Local AlphaFold2-PTM  |  PDB 8HFW chain A  |  WASE-screen v0.1.0{service_label}")

tabs, design_tab, property_tab = workbench_tabs(
    [
        tx("run_tab"),
        tx("nr_tab"),
        tx("prediction_tab"),
        tx("candidates_tab"),
        tx("structure_tab"),
        tx("history_tab"),
        tx("doctor_tab"),
    ],
    after=tx("structure_tab"),
)

with tabs[0]:
    st.subheader(tx("workbench"))
    metric_columns = st.columns(4)
    metric_columns[0].metric(tx("query_metric"), f"8HFW:{query_chain or 'A'}")
    metric_columns[1].metric(tx("database_metric"), tx(config_key))
    metric_columns[2].metric(tx("environment_metric"), tx("ready") if not blocking_fails else tx("blocked"))
    metric_columns[3].metric(tx("candidates_metric"), len(candidates))

    dot_color = COLORS["green"] if not blocking_fails else COLORS["red"]
    status_text = tx("ready") if not blocking_fails else f"{tx('blocked')} · {len(blocking_fails)}"
    st.markdown(
        f'<div class="status-line"><span class="status-dot" style="background:{dot_color}"></span><strong>{tx("preflight")}</strong><span>{status_text}</span><span>{tx("warnings")}: {len(preflight_warnings)}</span></div>',
        unsafe_allow_html=True,
    )
    if blocking_fails:
        for check in blocking_fails:
            st.error(doctor_message(check), icon=":material/error:")
    for check in preflight_fails:
        if check not in blocking_fails:
            st.warning(doctor_message(check), icon=":material/database:")

    controls = st.columns([1, 1, 1, 1.35])
    run_clicked = controls[0].button(
        tx("run"),
        type="primary",
        icon=":material/play_arrow:",
        width="stretch",
        disabled=state == "running" or bool(blocking_fails),
    )
    if controls[1].button(tx("cancel"), icon=":material/stop:", width="stretch", disabled=state != "running"):
        cancel_run()
        st.rerun()
    reannotate_clicked = controls[2].button(
        tx("reannotate"),
        icon=":material/sync:",
        width="stretch",
        disabled=selected_manifest is None or state == "running",
    )
    if controls[3].button(
        tx("open_folder"),
        icon=":material/folder_open:",
        width="stretch",
        disabled=selected_manifest is None or remote_access is not None,
    ) and selected_manifest:
        open_folder(selected_manifest["_run_dir"])

    if run_clicked:
        try:
            run_id = make_run_id("ui")
            runtime_config = write_runtime_config(
                base_config,
                query_path,
                query_chain or "A",
                {
                    "min_alntmscore": min_tm,
                    "max_rmsd": max_rmsd,
                    "min_qcov": min_qcov,
                    "min_tcov": min_tcov,
                    "min_alnlen": int(min_alnlen),
                    "max_sequence_identity": max_identity,
                },
                run_id,
            )
            start_run(runtime_config, run_id)
            st.rerun()
        except Exception as exc:  # noqa: BLE001 - surface configuration and process errors in the workbench.
            st.error(f"{tx('start_error')}: {type(exc).__name__}: {exc}")

    if reannotate_clicked and selected_manifest:
        try:
            with st.spinner(tx("reannotate")):
                reannotate_run(selected_manifest["_run_dir"])
            st.success(tx("annotation_done"))
            st.rerun()
        except Exception as exc:  # noqa: BLE001
            st.error(f"{type(exc).__name__}: {exc}")

    render_progress_fragment()

    st.subheader(tx("latest_results"))
    if candidates.empty:
        st.info(tx("no_results"))
    else:
        preview_columns = [column for column in ("target", "source_db", "alntmscore", "qcov", "tcov", "wase_screen_score", "wase_tier", "annotation_status") if column in candidates.columns]
        st.dataframe(candidates[preview_columns].head(12), hide_index=True, width="stretch", height=420)

with tabs[1]:
    st.subheader(tx("nr_workbench"))
    st.caption(tx("nr_wase_notice"))
    nr_checks = {
        check["name"]: check
        for check in preflight["checks"]
        if check["name"] in {"mmseqs2", "ncbi_blast_plus", "nr_demo_snapshot"}
    }
    snapshot_ready = nr_checks.get("nr_demo_snapshot", {}).get("status") == "pass"
    mmseqs_ready = nr_checks.get("mmseqs2", {}).get("status") == "pass"

    scope_label = st.segmented_control(
        tx("nr_scope"),
        options=[tx("nr_local_scope"), tx("nr_full_scope")],
        default=tx("nr_local_scope"),
        width="stretch",
    )
    source_run_id = st.selectbox(
        tx("nr_source_run"),
        options=run_options,
        index=0 if run_options else None,
        format_func=lambda value: f"{run_by_id[value].get('status', '?')} | {value}",
        disabled=not run_options or nr_state == "running",
        placeholder=tx("not_run"),
        key="nr_source_structure_run",
    )

    if scope_label == tx("nr_local_scope"):
        snapshot_columns = st.columns([1.25, 1.5, 1])
        nr_email = snapshot_columns[0].text_input(
            tx("nr_contact_email"),
            value=os.environ.get("AQUA_NCBI_EMAIL", ""),
            placeholder="name@example.org",
        )
        nr_consent = snapshot_columns[1].checkbox(tx("nr_consent"), value=False)
        build_snapshot = snapshot_columns[2].button(
            tx("nr_build_snapshot"),
            icon=":material/cloud_download:",
            width="stretch",
            disabled=(
                nr_state == "running"
                or snapshot_ready
                or "@" not in nr_email
                or not nr_consent
            ),
        )
        if snapshot_ready:
            st.success(tx("nr_snapshot"), icon=":material/verified:")
        else:
            st.info(tx("nr_snapshot_missing"), icon=":material/info:")

        local_controls = st.columns([1, 1, 1.2])
        run_local = local_controls[0].button(
            tx("nr_expand_local"),
            type="primary",
            icon=":material/hub:",
            width="stretch",
            disabled=(
                nr_state == "running"
                or not source_run_id
                or not snapshot_ready
                or not mmseqs_ready
            ),
        )
        if local_controls[1].button(
            tx("cancel"),
            icon=":material/stop:",
            width="stretch",
            disabled=nr_state != "running",
            key="cancel_nr_local",
        ):
            cancel_nr_process()
            st.rerun()
        local_controls[2].metric(
            tx("nr_integrity"),
            "nr_remote_snapshot" if snapshot_ready else tx("blocked"),
        )

        if build_snapshot:
            task_id = make_run_id("nr-demo-build")
            start_nr_process(
                task_id,
                [
                    "nr-demo-build",
                    "--query-fasta",
                    "data/query/Bur01m_seed.fasta",
                    "--email",
                    nr_email,
                    "--max-hits",
                    "2500",
                    "--confirm-ncbi",
                ],
            )
            st.rerun()
        if run_local and source_run_id:
            nr_run_id = make_run_id("nr-local")
            start_nr_process(
                nr_run_id,
                [
                    "nr-expand",
                    "--source-run",
                    source_run_id,
                    "--config",
                    "configs/nr_local_demo.yaml",
                    "--nr-run-id",
                    nr_run_id,
                ],
            )
            st.rerun()
    else:
        st.info(tx("nr_server_note"), icon=":material/dns:")
        server_controls = st.columns([1, 1, 1.35])
        export_job = server_controls[0].button(
            tx("nr_export_job"),
            type="primary",
            icon=":material/package_2:",
            width="stretch",
            disabled=nr_state == "running" or not source_run_id,
        )
        if server_controls[1].button(
            tx("cancel"),
            icon=":material/stop:",
            width="stretch",
            disabled=nr_state != "running",
            key="cancel_nr_server",
        ):
            cancel_nr_process()
            st.rerun()
        server_controls[2].metric(tx("nr_integrity"), "full_nr | SHA-256")
        if export_job and source_run_id:
            nr_run_id = make_run_id("nr-full")
            start_nr_process(
                nr_run_id,
                [
                    "nr-export-job",
                    "--source-run",
                    source_run_id,
                    "--config",
                    "configs/nr_full_server.yaml",
                    "--nr-run-id",
                    nr_run_id,
                ],
            )
            st.rerun()

        import_columns = st.columns([2.5, 1])
        server_archive = import_columns[0].file_uploader(
            tx("nr_result_zip"),
            type=["zip"],
            key="nr_server_result_upload",
        )
        import_clicked = import_columns[1].button(
            tx("nr_import_job"),
            icon=":material/upload_file:",
            width="stretch",
            disabled=nr_state == "running" or server_archive is None,
        )
        if import_clicked and server_archive is not None:
            UI_NR_INPUTS.mkdir(parents=True, exist_ok=True)
            uploaded_path = UI_NR_INPUTS / (
                f"{make_run_id('nr-import')}_{Path(server_archive.name).name}"
            )
            uploaded_path.write_bytes(server_archive.getvalue())
            start_nr_process(
                make_run_id("nr-import"),
                ["nr-import-job", "--bundle", str(uploaded_path)],
            )
            st.rerun()

    render_nr_progress_fragment()

    st.subheader(tx("nr_run_history"))
    if not nr_runs:
        st.info(tx("nr_no_runs"))
    else:
        nr_run_ids = list(nr_run_by_id)
        preferred_nr_run = st.session_state.get("prediction_nr_run_id", "")
        nr_default = nr_run_ids.index(preferred_nr_run) if preferred_nr_run in nr_run_ids else 0
        selected_nr_run_id = st.selectbox(
            tx("nr_run_history"),
            nr_run_ids,
            index=nr_default,
            format_func=lambda value: (
                f"{nr_run_by_id[value].get('database_scope', '?')} | "
                f"{nr_run_by_id[value].get('status', '?')} | {value}"
            ),
            key="selected_nr_run",
        )
        selected_nr_manifest = nr_run_by_id[selected_nr_run_id]
        selected_nr_candidates = read_nr_candidates(selected_nr_manifest)
        nr_metrics = st.columns(4)
        nr_metrics[0].metric(tx("status"), selected_nr_manifest.get("status", ""))
        nr_metrics[1].metric(
            tx("nr_representatives"),
            selected_nr_manifest.get("representative_count", 0),
        )
        nr_metrics[2].metric(
            tx("nr_af2_ready"),
            selected_nr_manifest.get("af2_ready_count", 0),
        )
        nr_metrics[3].metric(
            tx("nr_scope"),
            selected_nr_manifest.get("database_scope", ""),
        )
        for warning in selected_nr_manifest.get("warnings", []):
            st.warning(warning)

        selected_nr_events = read_nr_events(selected_nr_manifest)
        if selected_nr_events:
            st.markdown(f"**{tx('nr_pipeline_trace')}**")
            nr_event_rows = [
                {
                    "time": _clean_value(event.get("timestamp", ""))[11:19],
                    "stage": _clean_value(event.get("stage", "")),
                    "state": _clean_value(event.get("state", "")),
                    "progress": (
                        f"{float(event.get('progress')) * 100:.0f}%"
                        if event.get("progress") is not None
                        else ""
                    ),
                    "message": _clean_value(event.get("message", "")),
                }
                for event in selected_nr_events[-10:]
            ]
            st.dataframe(
                pd.DataFrame(nr_event_rows),
                hide_index=True,
                width="stretch",
                height=min(340, 38 * (len(nr_event_rows) + 1)),
            )

        if selected_nr_candidates.empty:
            st.info(tx("no_results"))
        else:
            filter_columns = st.columns(3)
            minimum_nr_qcov = filter_columns[0].slider(
                tx("nr_min_qcov"),
                min_value=0.0,
                max_value=1.0,
                value=0.60,
                step=0.05,
            )
            triad_only = filter_columns[1].toggle(tx("nr_triad"), value=False)
            af2_only = filter_columns[2].toggle(tx("nr_af2_ready"), value=False)
            filtered_nr = selected_nr_candidates.copy()
            filtered_nr = filtered_nr[
                pd.to_numeric(filtered_nr.get("qcov", 0), errors="coerce").fillna(0)
                >= minimum_nr_qcov
            ]
            if triad_only and "mapped_catalytic_identity_count" in filtered_nr:
                filtered_nr = filtered_nr[
                    pd.to_numeric(
                        filtered_nr["mapped_catalytic_identity_count"],
                        errors="coerce",
                    ).fillna(0)
                    == 3
                ]
            if af2_only and "af2_ready" in filtered_nr:
                filtered_nr = filtered_nr[
                    filtered_nr["af2_ready"].fillna(False).astype(bool)
                ]
            nr_display_columns = [
                column
                for column in (
                    "nr_rank",
                    "target_accession",
                    "organism",
                    "evalue",
                    "fident",
                    "qcov",
                    "mapped_catalytic_identity_count",
                    "af2_ready",
                    "wase_tier",
                    "sequence_screen_status",
                )
                if column in filtered_nr.columns
            ]
            nr_column_config: dict[str, Any] = {
                "nr_rank": st.column_config.NumberColumn("排名", format="%d", width="small"),
                "target_accession": st.column_config.TextColumn("accession", width="medium"),
                "organism": st.column_config.TextColumn(
                    "物种" if st.session_state.language == "zh" else "Organism",
                    width="medium",
                ),
                "evalue": st.column_config.NumberColumn(
                    "E-value",
                    format="%.2e",
                    width="small",
                ),
                "fident": st.column_config.NumberColumn(
                    "一致性" if st.session_state.language == "zh" else "Identity",
                    format="%.2f",
                    width="small",
                ),
                "qcov": st.column_config.NumberColumn(
                    "查询覆盖" if st.session_state.language == "zh" else "Query coverage",
                    format="%.2f",
                    width="small",
                ),
                "mapped_catalytic_identity_count": st.column_config.NumberColumn(
                    "三联体" if st.session_state.language == "zh" else "Triad",
                    format="%d / 3",
                    width="small",
                ),
                "af2_ready": st.column_config.CheckboxColumn("AF2-ready", width="small"),
                "wase_tier": st.column_config.TextColumn(
                    "WASE-screen",
                    width="medium",
                ),
                "sequence_screen_status": st.column_config.TextColumn(
                    "预筛状态" if st.session_state.language == "zh" else "Prescreen status",
                    width="medium",
                ),
            }
            st.markdown(f"**{tx('nr_sequence_evidence')}**")
            nr_display = filtered_nr[nr_display_columns].copy()
            if st.session_state.language == "zh" and "wase_tier" in nr_display:
                nr_display["wase_tier"] = nr_display["wase_tier"].replace(
                    {"Not evaluated - structure required": "未评估：待结构复核"}
                )
            st.dataframe(
                nr_display,
                hide_index=True,
                width="stretch",
                height=440,
                column_config=nr_column_config,
            )

            action_columns = st.columns([2.4, 1, 1, 1])
            nr_candidate_labels = [
                f"{_clean_value(row.get('target_accession')) or _clean_value(row.get('sequence_id'))} | "
                f"{_clean_value(row.get('organism')) or _clean_value(row.get('target'))}"
                for _, row in filtered_nr.iterrows()
            ]
            selected_nr_candidate = action_columns[0].selectbox(
                tx("nr_accession"),
                range(len(filtered_nr)),
                format_func=nr_candidate_labels.__getitem__,
                disabled=not nr_candidate_labels,
            )
            nr_row = filtered_nr.iloc[selected_nr_candidate] if selected_nr_candidate is not None else {}
            candidate_property_button(
                action_columns[2], _clean_value(nr_row.get("target_accession")) or _clean_value(nr_row.get("sequence_id")),
                _clean_value(nr_row.get("target_sequence")), "nr_to_properties",
            )
            if action_columns[1].button(
                tx("nr_send_af2"),
                icon=":material/model_training:",
                width="stretch",
                disabled=not nr_candidate_labels,
            ):
                st.session_state.prediction_nr_accession = (
                    _clean_value(nr_row.get("target_accession"))
                    or _clean_value(nr_row.get("sequence_id"))
                )
                st.session_state.prediction_nr_run_id = selected_nr_run_id
                st.info(tx("prediction_tab"))
            if action_columns[3].button(
                tx("open_folder"),
                icon=":material/folder_open:",
                width="stretch",
                disabled=remote_access is not None,
                key="open_nr_folder",
            ):
                open_folder(selected_nr_manifest["_nr_run_dir"])

            artifacts = selected_nr_manifest.get("artifacts", {})
            download_columns = st.columns(3)
            for column, key, label, mime in (
                (
                    download_columns[0],
                    "candidates_xlsx",
                    "download_excel",
                    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                ),
                (download_columns[1], "candidates_tsv", "download_tsv", "text/tab-separated-values"),
                (download_columns[2], "candidates_fasta", "download_fasta", "text/plain"),
            ):
                path_text = _clean_value(artifacts.get(key, ""))
                path = Path(path_text) if path_text else Path("__missing_nr_artifact__")
                column.download_button(
                    tx(label),
                    data=path.read_bytes() if path.is_file() else b"",
                    file_name=path.name if path.is_file() else key,
                    mime=mime,
                    icon=":material/download:",
                    width="stretch",
                    disabled=not path.is_file(),
                )

with tabs[2]:
    st.subheader(tx("af_workbench"))
    prediction_checks = {
        check["name"]: check for check in preflight["checks"] if check["name"] in {"wsl_gpu", "colabfold", "jax_cuda"}
    }
    prediction_blockers = [check for check in prediction_checks.values() if check["status"] != "pass"]
    if prediction_blockers:
        for check in prediction_blockers:
            st.warning(doctor_message(check), icon=":material/warning:")

    native_design_label = "本地设计序列" if st.session_state.language == "zh" else "Local design"
    input_labels = [
        tx("paste_fasta"),
        tx("upload_fasta"),
        tx("candidate_input"),
        tx("nr_candidate_input"),
        native_design_label,
    ]
    if st.session_state.get("prediction_nr_accession"):
        preferred_mode = tx("nr_candidate_input")
    elif st.session_state.get("prediction_candidate_accession"):
        preferred_mode = tx("candidate_input")
    else:
        preferred_mode = tx("paste_fasta")
    input_mode = st.segmented_control(
        tx("input_mode"),
        options=input_labels,
        default=preferred_mode,
        width="stretch",
    )
    prediction_record = None
    candidate_run_id = None
    candidate_nr_run_id = None
    candidate_accession_value = None
    input_error = ""
    if input_mode == tx("paste_fasta"):
        fasta_text = st.text_area(
            tx("fasta_text"),
            value=st.session_state.get("prediction_fasta_text", ""),
            height=180,
            placeholder=">protein_name\nMST...",
        )
        if fasta_text.strip():
            try:
                prediction_record = parse_single_fasta(fasta_text)
            except SequenceValidationError as exc:
                input_error = str(exc)
    elif input_mode == native_design_label:
        from bur01_structmine.design_view import select_design_for_prediction
        prediction_record = select_design_for_prediction()
    elif input_mode == tx("upload_fasta"):
        upload = st.file_uploader(tx("fasta_text"), type=["fa", "faa", "fasta", "txt"])
        if upload is not None:
            try:
                prediction_record = parse_single_fasta(upload.getvalue().decode("utf-8-sig"))
            except (UnicodeDecodeError, SequenceValidationError) as exc:
                input_error = str(exc)
    elif input_mode == tx("candidate_input"):
        if candidates.empty or not selected_manifest:
            st.info(tx("no_results"))
        else:
            candidate_labels = [
                f"{candidate_accession(row)} · {_clean_value(row.get('target'))}"
                for _, row in candidates.iterrows()
            ]
            preselected_accession = st.session_state.get("prediction_candidate_accession", "")
            default_index = next(
                (index for index, label in enumerate(candidate_labels) if preselected_accession and label.startswith(preselected_accession)),
                0,
            )
            selected_candidate_label = st.selectbox(tx("accession"), candidate_labels, index=default_index)
            selected_candidate_row = candidates.iloc[candidate_labels.index(selected_candidate_label)]
            candidate_accession_value = candidate_accession(selected_candidate_row)
            candidate_run_id = selected_manifest["run_id"]
            try:
                prediction_record = validate_sequence(
                    candidate_accession_value,
                    _clean_value(selected_candidate_row.get("final_protein_sequence")),
                )
            except SequenceValidationError as exc:
                input_error = str(exc)
    else:
        if not nr_runs:
            st.info(tx("nr_no_runs"))
        else:
            nr_run_ids = list(nr_run_by_id)
            preferred_nr_run = st.session_state.get("prediction_nr_run_id", "")
            default_nr_index = (
                nr_run_ids.index(preferred_nr_run)
                if preferred_nr_run in nr_run_ids
                else 0
            )
            candidate_nr_run_id = st.selectbox(
                tx("nr_run_history"),
                nr_run_ids,
                index=default_nr_index,
                format_func=lambda value: (
                    f"{nr_run_by_id[value].get('database_scope', '?')} | {value}"
                ),
                key="prediction_nr_run_select",
            )
            source_nr_candidates = read_nr_candidates(
                nr_run_by_id[candidate_nr_run_id]
            )
            if source_nr_candidates.empty:
                st.info(tx("no_results"))
            else:
                nr_labels = [
                    f"{_clean_value(row.get('target_accession')) or _clean_value(row.get('sequence_id'))} | "
                    f"{_clean_value(row.get('organism')) or _clean_value(row.get('target'))}"
                    for _, row in source_nr_candidates.iterrows()
                ]
                preferred_accession = st.session_state.get(
                    "prediction_nr_accession",
                    "",
                )
                default_candidate_index = next(
                    (
                        index
                        for index, label in enumerate(nr_labels)
                        if preferred_accession and label.startswith(preferred_accession)
                    ),
                    0,
                )
                selected_nr_label = st.selectbox(
                    tx("nr_accession"),
                    nr_labels,
                    index=default_candidate_index,
                    key="prediction_nr_accession_select",
                )
                selected_nr_row = source_nr_candidates.iloc[
                    nr_labels.index(selected_nr_label)
                ]
                candidate_accession_value = (
                    _clean_value(selected_nr_row.get("target_accession"))
                    or _clean_value(selected_nr_row.get("sequence_id"))
                )
                try:
                    prediction_record = validate_sequence(
                        candidate_accession_value,
                        _clean_value(selected_nr_row.get("target_sequence")),
                    )
                except SequenceValidationError as exc:
                    input_error = str(exc)

    if input_error:
        st.error(input_error, icon=":material/error:")

    mapping = None
    manual_overrides: dict[str, int] = {}
    if prediction_record is not None:
        try:
            mapping = map_8hfw_catalytic_sites(prediction_record.sequence)
            mapping_rows = []
            for label, site in mapping["sites"].items():
                mapping_rows.append(
                    {
                        "site": label,
                        "position": site.get("predicted_position"),
                        "observed": site.get("observed_amino_acid"),
                        "status": site.get("status"),
                    }
                )
            st.markdown(f"**{tx('map_preview')}**")
            st.dataframe(pd.DataFrame(mapping_rows), hide_index=True, width="stretch", height=176)
            manual_enabled = st.toggle(tx("manual_sites"), value=False)
            if manual_enabled:
                site_columns = st.columns(3)
                for column, label in zip(site_columns, ("S111", "D241", "H274")):
                    automatic = mapping["sites"][label].get("predicted_position") or 1
                    manual_overrides[label] = int(
                        column.number_input(
                            label,
                            min_value=1,
                            max_value=len(prediction_record.sequence),
                            value=int(automatic),
                            step=1,
                        )
                    )
        except Exception as exc:  # noqa: BLE001 - mapping failure must not hide the prediction controls.
            st.warning(f"{tx('map_preview')}: {type(exc).__name__}: {exc}")

    settings_columns = st.columns(2)
    preset_label = settings_columns[0].segmented_control(
        tx("preset"),
        options=[tx("quick"), tx("standard"), tx("high")],
        default=tx("standard"),
        width="stretch",
    )
    preset_lookup = {tx("quick"): "quick", tx("standard"): "standard", tx("high"): "high"}
    use_templates = settings_columns[1].toggle(tx("use_templates"), value=False)
    consent = st.checkbox(tx("privacy_confirm"), value=False)
    st.caption(tx("privacy_required"))

    prediction_controls = st.columns([1, 1.25, 1, 1.2])
    prediction_disabled = (
        prediction_record is None
        or not consent
        or prediction_state == "running"
        or bool(prediction_blockers)
    )
    predict_clicked = prediction_controls[0].button(
        tx("predict_only"),
        type="primary",
        icon=":material/model_training:",
        width="stretch",
        disabled=prediction_disabled,
    )
    predict_mine_clicked = prediction_controls[1].button(
        tx("predict_mine"),
        icon=":material/account_tree:",
        width="stretch",
        disabled=prediction_disabled,
    )
    if prediction_controls[2].button(
        tx("cancel"),
        icon=":material/stop:",
        width="stretch",
        disabled=prediction_state != "running",
        key="cancel_prediction",
    ):
        cancel_prediction_process()
        st.rerun()
    selected_prediction_for_open = predictions[0] if predictions else None
    if prediction_controls[3].button(
        tx("open_folder"),
        icon=":material/folder_open:",
        width="stretch",
        disabled=selected_prediction_for_open is None or remote_access is not None,
        key="open_prediction_folder",
    ) and selected_prediction_for_open:
        open_folder(selected_prediction_for_open["_prediction_dir"])

    if predict_clicked or predict_mine_clicked:
        try:
            prediction_id = make_run_id(f"AF2-{prediction_record.name}")
            fasta_path = None
            if input_mode in {tx("paste_fasta"), tx("upload_fasta"), native_design_label}:
                UI_PREDICTION_INPUTS.mkdir(parents=True, exist_ok=True)
                fasta_path = UI_PREDICTION_INPUTS / f"{prediction_id}.fasta"
                fasta_path.write_text(prediction_record.to_fasta(), encoding="ascii")
            start_prediction_process(
                prediction_id=prediction_id,
                preset=preset_lookup[preset_label],
                use_templates=use_templates,
                consent=consent,
                catalytic_sites=manual_overrides,
                fasta_path=fasta_path,
                candidate_run_id=candidate_run_id,
                candidate_nr_run_id=candidate_nr_run_id,
                accession=candidate_accession_value,
                mine_after=predict_mine_clicked,
            )
            st.session_state.pop("prediction_candidate_accession", None)
            st.session_state.pop("prediction_nr_accession", None)
            st.session_state.pop("prediction_nr_run_id", None)
            st.rerun()
        except Exception as exc:  # noqa: BLE001
            st.error(f"{tx('start_error')}: {type(exc).__name__}: {exc}")

    render_prediction_progress_fragment()

    st.subheader(tx("predictions"))
    if not predictions:
        st.info(tx("no_prediction"))
    else:
        prediction_ids = list(prediction_by_id)
        selected_prediction_id = st.selectbox(
            tx("select_prediction"),
            prediction_ids,
            format_func=lambda value: f"{prediction_by_id[value].get('status', '?')} · {value}",
        )
        selected_prediction = prediction_by_id[selected_prediction_id]
        quality_metrics = st.columns(4)
        quality_metrics[0].metric(tx("status"), selected_prediction.get("status", ""))
        quality_metrics[1].metric(tx("model_count"), len(selected_prediction.get("models", [])))
        quality_metrics[2].metric(
            tx("mean_plddt"),
            f"{float(selected_prediction.get('mean_plddt')):.1f}" if selected_prediction.get("mean_plddt") is not None else "—",
        )
        quality_metrics[3].metric(
            tx("ptm"),
            f"{float(selected_prediction.get('ptm')):.3f}" if selected_prediction.get("ptm") is not None else "—",
        )
        if selected_prediction.get("error"):
            st.error(selected_prediction["error"])
        for warning in selected_prediction.get("warnings", []):
            st.warning(warning)
        model_rows = [
            {
                "rank": model.get("rank"),
                "model": model.get("model_name"),
                "mean_pLDDT": model.get("mean_plddt"),
                "pTM": model.get("ptm"),
                "max_PAE": model.get("max_pae"),
                "relaxed": model.get("relaxed"),
            }
            for model in selected_prediction.get("models", [])
        ]
        if model_rows:
            st.markdown(f"**{tx('model_ranking')}**")
            st.dataframe(pd.DataFrame(model_rows), hide_index=True, width="stretch", height=min(300, 38 * (len(model_rows) + 1)))
        artifacts = selected_prediction.get("artifacts", {})
        scores_path = Path(_clean_value(artifacts.get("rank_001_scores", ""))) if _clean_value(artifacts.get("rank_001_scores", "")) else None
        if scores_path and scores_path.is_file():
            scores = json.loads(scores_path.read_text(encoding="utf-8"))
            chart_columns = st.columns(2)
            plddt = [float(value) for value in scores.get("plddt", [])]
            if plddt:
                plddt_figure = go.Figure(go.Scatter(x=list(range(1, len(plddt) + 1)), y=plddt, mode="lines", line={"color": COLORS["green"], "width": 2}))
                plddt_figure.update_layout(title=tx("plddt_curve"), xaxis_title="Residue", yaxis_title="pLDDT", yaxis_range=[0, 100], margin={"l": 20, "r": 20, "t": 45, "b": 20}, paper_bgcolor=COLORS["background"], plot_bgcolor="#FFFFFF")
                chart_columns[0].plotly_chart(plddt_figure, width="stretch", config={"displaylogo": False})
            pae = scores.get("pae", [])
            if pae:
                step = max(1, len(pae) // 500)
                reduced_pae = [row[::step] for row in pae[::step]]
                pae_figure = go.Figure(go.Heatmap(z=reduced_pae, colorscale="Viridis_r", zmin=0, zmax=30, colorbar={"title": "Å"}))
                pae_figure.update_layout(title=tx("pae"), xaxis_title="Residue", yaxis_title="Residue", margin={"l": 20, "r": 20, "t": 45, "b": 20}, paper_bgcolor=COLORS["background"])
                chart_columns[1].plotly_chart(pae_figure, width="stretch", config={"displaylogo": False})
        rank_path_text = _clean_value(selected_prediction.get("rank_001_path", ""))
        rank_path = Path(rank_path_text) if rank_path_text else None
        if rank_path and rank_path.is_file():
            st.markdown(f"**{tx('predicted_structure')}**")
            render_prediction_structure(rank_path)

with tabs[3]:
    st.subheader(tx("candidate_table"))
    if candidates.empty:
        st.info(tx("no_results"))
    else:
        filtered = candidates.copy()
        filter_columns = st.columns(4)
        source_values = sorted(_clean_value(value) for value in filtered.get("source_db", pd.Series(dtype=str)).dropna().unique())
        selected_sources = filter_columns[0].multiselect(tx("source_filter"), source_values, default=source_values)
        tier_values = sorted(_clean_value(value) for value in filtered.get("wase_tier", pd.Series(dtype=str)).dropna().unique())
        selected_tiers = filter_columns[1].multiselect(tx("tier_filter"), tier_values, default=tier_values)
        reviewed_choice = filter_columns[2].selectbox(tx("reviewed_filter"), [tx("any"), tx("reviewed"), tx("unreviewed"), tx("unknown")])
        min_evidence = filter_columns[3].slider(tx("evidence_filter"), min_value=0, max_value=100, value=0, step=5)
        tm_filter = st.slider(tx("min_tm"), min_value=0.0, max_value=1.0, value=0.0, step=0.05, key="candidate_tm")
        if selected_sources and "source_db" in filtered:
            filtered = filtered[filtered["source_db"].fillna("").astype(str).isin(selected_sources)]
        if selected_tiers and "wase_tier" in filtered:
            filtered = filtered[filtered["wase_tier"].fillna("").astype(str).isin(selected_tiers)]
        tm_series = pd.to_numeric(filtered.get("alntmscore", 0), errors="coerce").fillna(0)
        filtered = filtered[tm_series >= tm_filter]
        evidence_series = pd.to_numeric(filtered.get("wase_evidence_coverage", 0), errors="coerce").fillna(0)
        filtered = filtered[evidence_series >= min_evidence]
        if reviewed_choice != tx("any") and "uniprot_is_reviewed" in filtered:
            reviewed = filtered["uniprot_is_reviewed"].fillna("").astype(str).str.lower()
            if reviewed_choice == tx("reviewed"):
                filtered = filtered[reviewed.isin({"true", "reviewed", "1"})]
            elif reviewed_choice == tx("unreviewed"):
                filtered = filtered[reviewed.isin({"false", "unreviewed", "0"})]
            else:
                filtered = filtered[~reviewed.isin({"true", "reviewed", "1", "false", "unreviewed", "0"})]

        display_columns = [
            column
            for column in (
                "rank_after_filter",
                "target",
                "source_db",
                "search_engine",
                "alntmscore",
                "qtmscore",
                "ttmscore",
                "qcov",
                "tcov",
                "rmsd",
                "wase_screen_score",
                "wase_evidence_coverage",
                "wase_tier",
                "uniprot_is_reviewed",
                "annotation_status",
            )
            if column in filtered.columns
        ]
        st.dataframe(filtered[display_columns], hide_index=True, width="stretch", height=480)

        if not filtered.empty and "wase_screen_score" in filtered.columns:
            chart = filtered.copy()
            chart["min_tm"] = pd.concat(
                [
                    pd.to_numeric(chart.get("qtmscore", chart.get("alntmscore", 0)), errors="coerce"),
                    pd.to_numeric(chart.get("ttmscore", chart.get("alntmscore", 0)), errors="coerce"),
                ],
                axis=1,
            ).min(axis=1)
            fig = px.scatter(
                chart,
                x="min_tm",
                y="wase_screen_score",
                color="wase_tier" if "wase_tier" in chart else None,
                hover_name="target",
                hover_data=[column for column in ("source_db", "qcov", "tcov", "rmsd", "wase_evidence_coverage") if column in chart],
                color_discrete_map={
                    "High priority": COLORS["green"],
                    "Medium priority": COLORS["amber"],
                    "Low priority": COLORS["red"],
                    "Incomplete evidence": COLORS["muted"],
                },
                labels={"min_tm": "min(qTM, tTM)", "wase_screen_score": "WASE-screen"},
            )
            fig.update_layout(
                title=tx("evidence_map"),
                paper_bgcolor=COLORS["background"],
                plot_bgcolor="#FFFFFF",
                font={"family": "IBM Plex Sans", "color": COLORS["ink"]},
                margin={"l": 20, "r": 20, "t": 50, "b": 20},
                legend_title_text="",
            )
            fig.update_xaxes(range=[0, 1], gridcolor=COLORS["line"])
            fig.update_yaxes(range=[0, 100], gridcolor=COLORS["line"])
            st.plotly_chart(fig, width="stretch", config={"displaylogo": False})

        if selected_manifest:
            artifacts = selected_manifest.get("artifacts", {})
            downloads = st.columns(3)
            for column, artifact_key, label_key, mime in (
                (downloads[0], "annotated_xlsx", "download_excel", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"),
                (downloads[1], "annotated_tsv", "download_tsv", "text/tab-separated-values"),
                (downloads[2], "candidate_fasta", "download_fasta", "text/plain"),
            ):
                path_text = _clean_value(artifacts.get(artifact_key, ""))
                path = Path(path_text) if path_text else Path("__missing_artifact__")
                column.download_button(
                    tx(label_key),
                    data=path.read_bytes() if path.is_file() else b"",
                    file_name=path.name if path.is_file() else artifact_key,
                    mime=mime,
                    icon=":material/download:",
                    width="stretch",
                    disabled=not path.is_file(),
                )

        repredict_columns = st.columns([3, 1, 1])
        repredict_labels = [
            f"{candidate_accession(row)} · {_clean_value(row.get('target'))}"
            for _, row in filtered.iterrows()
        ]
        repredict_choice = repredict_columns[0].selectbox(
            tx("accession"),
            range(len(filtered)),
            format_func=repredict_labels.__getitem__,
            disabled=filtered.empty,
            key="candidate_repredict_choice",
        )
        chosen_row = filtered.iloc[repredict_choice] if repredict_choice is not None else {}
        candidate_property_button(repredict_columns[2], candidate_accession(chosen_row),
                                  _clean_value(chosen_row.get("final_protein_sequence")), "candidate_to_properties")
        if repredict_columns[1].button(
            tx("retrigger"),
            icon=":material/model_training:",
            width="stretch",
            disabled=filtered.empty,
        ):
            st.session_state.prediction_candidate_accession = candidate_accession(chosen_row)
            st.info(tx("prediction_tab"))

if design_tab.open:
    with design_tab:
        render_local_design(candidates, predictions, open_model_tab)
if property_tab.open:
    with property_tab:
        show_property_tab()

with tabs[4]:
    if st.button(model_text("打开序列设计", "Open sequence design"), key="structure_to_design"):
        open_model_tab("design")
    st.subheader(tx("overlay"))
    if candidates.empty or not selected_manifest:
        st.info(tx("no_results"))
    else:
        labels = [f"{index + 1}. {_clean_value(row.get('target'))}" for index, (_, row) in enumerate(candidates.iterrows())]
        selected_label = st.selectbox(tx("select_candidate"), labels)
        selected_index = labels.index(selected_label)
        row = candidates.iloc[selected_index]
        target_path = candidate_structure_path(row)
        query_structure = Path(selected_manifest.get("query_path", ROOT / "data" / "query" / "Bur01m.pdb"))
        target_chain = _clean_value(row.get("target_chain")) or _clean_value(row.get("chain"))
        catalytic_site_map = selected_manifest.get("catalytic_reference_sites") or {"S111": 111, "D241": 241, "H274": 274}
        catalytic_query_residues = tuple(
            (int(position), label[0].upper())
            for label, position in catalytic_site_map.items()
            if label and label[0].upper() in {"S", "D", "H"}
        )
        catalytic_numbers = tuple(position for position, _ in catalytic_query_residues)
        if target_path is None or not target_path.exists() or not query_structure.exists():
            st.warning(tx("no_structure"))
        else:
            try:
                overlay = cached_overlay(
                    str(query_structure),
                    query_structure.stat().st_mtime_ns,
                    str(target_path),
                    target_path.stat().st_mtime_ns,
                    selected_manifest.get("query_chain", "A") or "A",
                    target_chain,
                    catalytic_query_residues,
                )
                metrics = st.columns(4)
                metrics[0].metric("qTM", f"{overlay.qtmscore:.3f}")
                metrics[1].metric("tTM", f"{overlay.ttmscore:.3f}")
                metrics[2].metric("RMSD", f"{overlay.rmsd:.2f} Å")
                metrics[3].metric(tx("target_chain"), overlay.target_chain)
                st.markdown(
                    f'<div class="legend-row"><span><i class="swatch" style="background:{COLORS["green"]}"></i>{tx("query_legend")}</span><span><i class="swatch" style="background:{COLORS["amber"]}"></i>{tx("target_legend")}</span><span><i class="swatch" style="background:{COLORS["red"]}"></i>{tx("catalytic_legend")}</span></div>',
                    unsafe_allow_html=True,
                )
                render_overlay(overlay, catalytic_numbers)
                st.subheader(tx("alignment"))
                site_labels = " · ".join(f"{amino_acid}{position}" for position, amino_acid in catalytic_query_residues)
                st.markdown(site_labels + "  →  " + html.escape(overlay.mapped_catalytic_residues))
                st.markdown(alignment_html(overlay, catalytic_numbers), unsafe_allow_html=True)
            except Exception as exc:  # noqa: BLE001
                st.error(f"{type(exc).__name__}: {exc}")

with tabs[5]:
    st.subheader(tx("history"))
    if not runs:
        st.info(tx("not_run"))
    else:
        history_rows = []
        for run in runs:
            config_snapshot = run.get("config_snapshot", {})
            sources = config_snapshot.get("search_databases") or []
            if isinstance(sources, dict):
                source_list = list(sources.values())
            else:
                source_list = sources
            engines = ", ".join(sorted({_clean_value(source.get("search_engine")) or "foldseek" for source in source_list if isinstance(source, dict)}))
            history_rows.append(
                {
                    "run_id": run.get("run_id", ""),
                    tx("status"): run.get("status", ""),
                    tx("search_engine"): engines,
                    tx("candidates_metric"): run.get("annotated_hit_count", 0),
                    "created_at": run.get("created_at", ""),
                }
            )
        st.dataframe(pd.DataFrame(history_rows), hide_index=True, width="stretch", height=min(520, 38 * (len(history_rows) + 1)))
        detail_run_id = st.selectbox(tx("run_details"), run_options, format_func=lambda run_id: f"{run_by_id[run_id].get('status', '?')} · {run_id}", key="history_run")
        detail = run_by_id[detail_run_id]
        detail_columns = st.columns(4)
        detail_columns[0].metric(tx("status"), detail.get("status", ""))
        detail_columns[1].metric("Raw", detail.get("raw_hit_count", 0))
        detail_columns[2].metric("Filtered", detail.get("filtered_hit_count", 0))
        detail_columns[3].metric("Annotated", detail.get("annotated_hit_count", 0))
        if detail.get("error"):
            st.error(detail["error"])
        for warning in detail.get("warnings", []):
            st.warning(warning)
        history_actions = st.columns(2)
        if history_actions[0].button(
            tx("open_folder"),
            icon=":material/folder_open:",
            width="stretch",
            key="history_open",
            disabled=remote_access is not None,
        ):
            open_folder(detail["_run_dir"])
        if history_actions[1].button(tx("reannotate"), icon=":material/sync:", width="stretch", key="history_annotate"):
            try:
                with st.spinner(tx("reannotate")):
                    reannotate_run(detail["_run_dir"])
                st.success(tx("annotation_done"))
                st.rerun()
            except Exception as exc:  # noqa: BLE001
                st.error(f"{type(exc).__name__}: {exc}")

with tabs[6]:
    st.subheader(tx("diagnostic_results"))
    diagnostic_controls = st.columns([1, 1, 2])
    network_check = diagnostic_controls[0].toggle(tx("network"), value=False)
    if diagnostic_controls[1].button(tx("rerun_doctor"), icon=":material/medical_services:", width="stretch"):
        cached_doctor.clear()
        st.session_state.doctor_report = doctor_for(base_config, network=network_check)
    diagnostic = st.session_state.get("doctor_report") or preflight
    summary = diagnostic["summary"]
    diagnostic_metrics = st.columns(3)
    diagnostic_metrics[0].metric("OK", summary["pass"])
    diagnostic_metrics[1].metric("WARN", summary["warn"])
    diagnostic_metrics[2].metric("FAIL", summary["fail"])
    doctor_rows = [
        {
            tx("status"): check["status"].upper(),
            "category": check["category"],
            "check": check["name"],
            "message": doctor_message(check),
        }
        for check in diagnostic["checks"]
    ]
    st.dataframe(pd.DataFrame(doctor_rows), hide_index=True, width="stretch", height=560)
    st.download_button(
        tx("download_json"),
        data=json.dumps(diagnostic, indent=2, ensure_ascii=False).encode("utf-8"),
        file_name=f"doctor_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
        mime="application/json",
        icon=":material/download:",
    )
