"""Reusable bilingual evidence, experiment and coordinate views."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import plotly.express as px
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

from .experiments import (COMPARISON, evaluate_experiments, evaluation_excel, read_experiments,
                          safe_table, simulated_experiments, template_bytes)
from .local_viewer import offline_viewer_html
from .scoring import read_score_display
from .sites import (ligand_contacts, nonprotein_components, parse_structure_bytes, protein_sites,
                    site_selections, viewer_selection)


def tr(zh: str, en: str) -> str:
    return zh if st.session_state.get("language", "zh") == "zh" else en


def score_details(row: pd.Series, key: str = "score") -> None:
    display = read_score_display(row)
    columns = st.columns(3)
    columns[0].metric("WASE-screen", f"{float(row.get('wase_screen_score', 0)):.1f} / 100")
    missing = display["information_missingness"]
    columns[1].metric(tr("信息缺失程度", "Information missingness"), f"{missing:.1f}%" if missing is not None else "N/A")
    columns[2].metric(tr("预测误差", "Prediction error"), tr("尚未估计", "Not estimated"))
    st.caption(tr("预筛分不代表活性、S/H 或实验成功概率。", "Screening score does not predict activity, S/H or success probability.") + f"  v{display['version']}")
    st.caption(tr("信息缺失程度仅针对当前评分字段；尚未评估的机制证据仍列在缺失项中。", "Missingness covers rubric fields only; unassessed mechanistic evidence remains listed below."))
    names = [("structure_fit_score", tr("结构匹配", "Structure fit"), 45),
             ("catalytic_core_score", tr("催化核心", "Catalytic core"), 30),
             ("annotation_confidence_score", tr("注释可靠性", "Annotation confidence"), 15),
             ("experimental_readiness_score", tr("实验准备度", "Experimental readiness"), 10)]
    detail = []
    try:
        explanation = json.loads(row.get("wase_score_explanation", "{}"))
    except (ValueError, TypeError):
        explanation = {}
    for field, name, maximum in names:
        item = explanation.get(field.removesuffix("_score"), {})
        chinese_basis = {"structure_fit_score": "25×min(qTM,tTM)+10×sqrt(qcov×tcov)+10×clamp(1−RMSD/5)",
                         "catalytic_core_score": "映射Ser中心模体5；S/D/H身份各5；两项几何距离各5",
                         "annotation_confidence_score": "相关功能5；审阅或实验结构5；结构可信度5",
                         "experimental_readiness_score": "序列4；长度2；结构覆盖2；可追溯ID 2"}
        basis = tr(chinese_basis[field], item["basis"]) if item.get("basis") else tr("历史记录未保存明细", "Details absent in historical record")
        detail.append({tr("项目", "Component"): name, tr("得分", "Score"): row.get(field), tr("满分", "Maximum"): maximum,
                       tr("计算依据", "Basis"): basis,
                       tr("缺失项", "Missing"): "; ".join(item.get("missing", []))})
    st.dataframe(pd.DataFrame(detail), hide_index=True, width="stretch")
    exclusions = row.get("wase_high_priority_exclusion_reasons", "")
    if isinstance(exclusions, str) and exclusions:
        labels = {"score_below_70": "总分低于70", "coverage_below_75": "可评估权重低于75", "min_tm_below_0.5": "双向TM最小值低于0.5", "catalytic_below_20": "催化核心低于20"}
        st.warning(tr("未列为高优先级的原因：", "High-priority exclusions: ") + tr("；".join(labels.get(code, code) for code in exclusions.split(";")), exclusions))
    geometry = row.get("catalytic_geometry_status", "not_assessable")
    status_labels = {"compatible": tr("已计算，距离符合条件", "Assessed, distances in range"),
                     "out_of_range": tr("已计算，距离不符合条件", "Assessed, distances out of range"),
                     "not_assessable": tr("无法完整评估", "Not fully assessable")}
    st.write(tr("催化几何：", "Catalytic geometry: ") + status_labels.get(str(geometry), str(geometry)))
    fields = [field for field in ("query_catalytic_residues", "mapped_catalytic_residues", "mapped_catalytic_ser_motif",
                                   "has_GXSXG_like_motif", "catalytic_ser_his_distance", "catalytic_asp_his_distance",
                                   "catalytic_geometry_missing", "catalytic_evidence_status", "source_run_id", "wase_score_delta") if field in row]
    with st.expander(tr("位点依据与缺失证据", "Site evidence and missing evidence")):
        st.dataframe(pd.DataFrame({"field": fields, "value": [str(row[field]) for field in fields]}), hide_index=True, width="stretch")
        st.write(str(row.get("wase_missing_evidence", "")))


def experimental_view(*, demo: bool = False, key: str = "experiments", candidates: pd.DataFrame | None = None, candidate_run_id: str = "") -> None:
    st.subheader(tr("实验评估", "Experimental assessment"))
    if not demo:
        st.download_button(tr("下载 Excel 模板", "Download Excel template"), template_bytes(), "experiment_template.xlsx",
                           icon=":material/download:", key=f"{key}_template")
        upload = st.file_uploader(tr("实验数据 Excel / CSV", "Experimental data Excel / CSV"), type=["xlsx", "csv"], key=f"{key}_upload")
        simulated = st.toggle(tr("模拟数据示例", "Simulated example"), value=False, key=f"{key}_simulated")
        if upload:
            try:
                raw = read_experiments(upload.getvalue(), upload.name)
            except Exception as exc:
                st.error(str(exc))
                return
        elif simulated:
            raw = simulated_experiments()
        else:
            st.info(tr("尚未导入实验记录", "No experimental records imported"))
            return
    else:
        raw = simulated_experiments()
    if "data_kind" in raw and raw.data_kind.eq("simulated").any():
        st.warning(tr("模拟数据，仅用于展示评估方法；不是 Bur01 实验结果。", "SIMULATED DATA: method demonstration only, not Bur01 experimental results."))
    if not demo and candidates is not None and not candidates.empty:
        if st.checkbox(tr("关联当前候选的 WASE-screen 分数", "Link WASE-screen from current candidates"), key=f"{key}_link"):
            lookup = {}
            versions = {}
            ambiguous = set()
            for _, candidate in candidates.iterrows():
                for field in ("target", "target_uniprot_accession"):
                    identifier = str(candidate.get(field, ""))
                    if identifier in {"", "nan"}:
                        continue
                    value = candidate.get("wase_screen_score")
                    if identifier in lookup and lookup[identifier] != value:
                        ambiguous.add(identifier)
                    lookup[identifier] = value
                    versions[identifier] = candidate.get("wase_version", "unknown")
            linked = raw.get("candidate_id", pd.Series(index=raw.index, dtype=str)).map({k: v for k, v in lookup.items() if k not in ambiguous})
            raw = raw.copy()
            existing = pd.to_numeric(raw.get("wase_screen_score", pd.Series(index=raw.index, dtype=float)), errors="coerce")
            raw["wase_screen_score"] = existing.fillna(linked)
            for field in ("wase_source_run_id", "wase_version"):
                if field not in raw:
                    raw[field] = ""
            newly_linked = existing.isna() & linked.notna()
            raw.loc[newly_linked, "wase_source_run_id"] = candidate_run_id
            raw.loc[newly_linked, "wase_version"] = raw.loc[newly_linked, "candidate_id"].map(versions)
            st.caption(tr("成功匹配行数：", "Matched rows: ") + str(int(linked.notna().sum())))
    st.dataframe(raw, hide_index=True, width="stretch", height=250)
    result = evaluate_experiments(raw)
    if not result.errors.empty:
        st.error(tr("以下行未纳入评估，原始记录仍保留。", "Rows below are excluded; original records are retained."))
        st.dataframe(result.errors, hide_index=True, width="stretch")
    st.caption(tr("有效记录 / 原始记录：", "Valid / original rows: ") + f"{result.summary.valid_rows} / {result.summary.raw_rows}")
    if result.aggregates.empty:
        return
    st.markdown("**" + tr("重复测量汇总", "Replicate summary") + "**")
    st.dataframe(result.aggregates, hide_index=True, width="stretch")
    groups = list(result.valid.groupby(COMPARISON, sort=False, dropna=False))
    group_index = st.selectbox(tr("指标、单位与条件", "Endpoint, unit and conditions"), list(range(len(groups))),
                               format_func=lambda i: " | ".join(str(value) for value in groups[i][0]), key=f"{key}_group")
    group_key, group = groups[group_index]
    threshold_enabled = st.checkbox(tr("计算 Top-k 命中率", "Calculate Top-k hit rate"), key=f"{key}_top_enabled")
    k, threshold, higher = 5, None, True
    if threshold_enabled:
        columns = st.columns(3)
        k = columns[0].number_input("Top-k", 1, 1000, 3, key=f"{key}_k")
        threshold = columns[1].number_input(tr("成功阈值", "Success threshold") + f" ({group_key[1]})", value=30.0, key=f"{key}_threshold")
        higher = columns[2].selectbox(tr("成功条件", "Success condition"), [">=", "<="], key=f"{key}_direction") == ">="
    selected = evaluate_experiments(group.drop(columns=["source_row"], errors="ignore"), top_k=k, success_threshold=threshold, higher_is_success=higher)
    st.dataframe(selected.comparisons, hide_index=True, width="stretch")
    paired = selected.aggregates.dropna(subset=["predicted_value"])
    if not paired.empty:
        chart = px.scatter(paired, x="predicted_value", y="mean", error_y="std", hover_name="sample_id", color="prediction_source",
                           labels={"predicted_value": tr("预测值", "Predicted") + f" ({group_key[1]})", "mean": tr("实测均值", "Measured mean") + f" ({group_key[1]})"})
        chart.update_layout(margin=dict(l=15, r=15, t=15, b=15), height=350)
        st.plotly_chart(chart, width="stretch", key=f"{key}_scatter")
    rankings = selected.aggregates.dropna(subset=["wase_screen_score"]).copy()
    if not rankings.empty and group_key[0] in {"activity", "relative_activity", "residual_activity", "S/H"}:
        rankings["measured_rank"] = rankings["mean"].rank(ascending=False)
        rankings["wase_rank"] = rankings.wase_screen_score.rank(ascending=False)
        chart = px.scatter(rankings, x="wase_rank", y="measured_rank", hover_name="sample_id", text="sample_id",
                           labels={"wase_rank": "WASE-screen " + tr("排序", "rank"), "measured_rank": tr("实测排序", "Measured rank")})
        chart.update_traces(marker=dict(color="#0F6B5B", size=10), textposition="top center")
        chart.update_layout(height=320, margin=dict(l=15, r=15, t=20, b=15))
        st.plotly_chart(chart, width="stretch", key=f"{key}_ranks")
    if not result.unmatched.empty:
        with st.expander(tr("未匹配记录与原因", "Unmatched records and reasons")):
            st.dataframe(result.unmatched, hide_index=True, width="stretch")
    if not demo:
        st.download_button(tr("导出完整评估 Excel", "Export full evaluation Excel"), evaluation_excel(result), "experimental_evaluation.xlsx",
                           icon=":material/download:", key=f"{key}_export")
        st.download_button(tr("导出当前分组评估（含阈值）", "Export selected group with threshold"), evaluation_excel(selected), "selected_group.xlsx",
                           icon=":material/download:", key=f"{key}_export_group")


def structure_sites_view(default_path: Path, *, key: str = "sites", read_only: bool = False) -> None:
    st.subheader(tr("底物与位点", "Substrates and sites"))
    uploaded = None if read_only else st.file_uploader(tr("上传 PDB / CIF / mmCIF", "Upload PDB / CIF / mmCIF"), type=["pdb", "cif", "mmcif"], key=f"{key}_upload")
    try:
        content = uploaded.getvalue() if uploaded else default_path.read_bytes()
        doc = parse_structure_bytes(content, uploaded.name if uploaded else default_path.name)
    except Exception as exc:
        st.error(str(exc))
        return
    columns = st.columns(2)
    model_id = columns[0].selectbox(tr("模型", "Model"), doc.models, key=f"{key}_model")
    chains = columns[1].multiselect(tr("蛋白链", "Protein chains"), doc.chains(model_id), default=doc.chains(model_id)[:1], key=f"{key}_{doc.sha256}_{model_id}_chains")
    components_table = nonprotein_components(doc, model_id)
    component_labels = dict(zip(components_table.residue_id, components_table.label))
    ligand_ids = st.multiselect(tr("已有配体 / 辅因子", "Existing ligands / cofactors"), list(component_labels),
                               default=components_table.loc[components_table.default_selected, "residue_id"].tolist(),
                               format_func=component_labels.get, key=f"{key}_{doc.sha256}_{model_id}_ligands")
    cutoff = st.slider(tr("重原子接触距离 (Å)", "Heavy-atom contact cutoff (Å)"), 1.0, 10.0, 4.0, 0.5, key=f"{key}_cutoff")
    if not ligand_ids:
        st.info(tr("没有选定配体坐标；底物信息仅作记录，不生成结合姿势。", "No ligand coordinates selected. Substrate metadata only; no binding pose is generated."))
    contacts = ligand_contacts(doc, model_id, chains, ligand_ids, cutoff)
    sites = protein_sites(doc, model_id, chains)
    if sites.empty:
        st.warning(tr("当前链没有蛋白残基", "No protein residues in selected chains"))
        return
    labels = dict(zip(sites.residue_id, sites.label))
    triad_default = sites.loc[(sites.chain_id == "A") & sites.residue_number.isin([111, 241, 274]) & sites.insertion_code.eq(""), "residue_id"].tolist() if not uploaded else []
    triad = st.multiselect(tr("催化位点（结构编号）", "Catalytic sites (structure numbering)"), list(labels), default=triad_default,
                          format_func=labels.get, key=f"{key}_{doc.sha256}_{model_id}_triad")
    st.caption(tr("距离邻近不等于已确认氢键。", "Distance proximity is not a confirmed hydrogen bond."))
    selection = st.dataframe(contacts, hide_index=True, width="stretch", on_select="rerun", selection_mode="multi-row", key=f"{key}_contacts")
    selected_ids = [contacts.iloc[index].residue_id for index in selection.selection.rows if index < len(contacts)]
    highlight = st.multiselect(tr("高亮残基", "Highlight residues"), list(labels), format_func=labels.get, key=f"{key}_{doc.sha256}_{model_id}_highlight")
    viewer = py3Dmol.view(width="100%", height=440)
    viewer.addModel(doc.model_cif(model_id), "cif")
    viewer.setStyle({}, {})
    viewer.setStyle({"chain": chains, "hetflag": False}, {"cartoon": {"color": "#0F6B5B"}})
    for residue_id in ligand_ids:
        viewer.addStyle(viewer_selection(residue_id), {"stick": {"colorscheme": "orangeCarbon", "radius": 0.22}})
    for residue_id in set(contacts.residue_id) | set(highlight) | set(selected_ids):
        viewer.addStyle(viewer_selection(residue_id), {"stick": {"color": "#C77D00", "radius": 0.15}})
    for residue_id in triad:
        viewer.addStyle(viewer_selection(residue_id), {"stick": {"color": "#B94A48", "radius": 0.24}})
    for residue_id in set(highlight) | set(selected_ids):
        viewer.addStyle(viewer_selection(residue_id), {"sphere": {"color": "#C77D00", "scale": 0.35}})
        viewer.addLabel(labels[residue_id], {"fontSize": 12, "backgroundOpacity": 0.75}, viewer_selection(residue_id))
    viewer.setBackgroundColor("#F7F9F8")
    viewer.zoomTo({"chain": chains})
    components.html(offline_viewer_html(viewer), height=460)
    if read_only:
        return
    st.markdown("**" + tr("位点标记", "Site marks") + "**")
    marks = {}
    for mark, zh, en in [("keep", "保留", "Keep"), ("focus", "关注", "Focus"), ("exclude", "排除", "Exclude")]:
        values = st.multiselect(tr(zh, en), list(labels), format_func=labels.get, key=f"{key}_{doc.sha256}_{model_id}_{mark}")
        for value in values:
            if value in marks:
                st.error(tr("同一位点有冲突标记", "Conflicting marks for one site") + ": " + labels[value])
                return
            marks[value] = mark
    substrate = st.text_input(tr("底物名称", "Substrate name"), key=f"{key}_name")
    smiles = st.text_input("SMILES", key=f"{key}_smiles")
    conditions = st.text_input(tr("实验条件", "Experimental conditions"), key=f"{key}_conditions")
    records = site_selections(doc, model_id, marks, ligand_ids, cutoff, substrate, smiles, conditions)
    metadata = {"structure_sha256": doc.sha256, "model_id": model_id, "chains": chains, "cutoff_angstrom": cutoff,
                "ligand_ids": ligand_ids, "substrate_name": substrate, "smiles": smiles, "conditions": conditions,
                "coordinate_status": "coordinates_present" if ligand_ids else "metadata_only",
                "sites": [record.model_dump() for record in records]}
    st.download_button(tr("导出位点与底物记录", "Export sites and substrate record"), json.dumps(metadata, ensure_ascii=False, indent=2),
                       "site_selection.json", icon=":material/download:", key=f"{key}_export")
    st.download_button(tr("导出接触距离 CSV", "Export contact distances CSV"), safe_table(contacts).to_csv(index=False).encode("utf-8-sig"),
                       "contacts.csv", icon=":material/download:", key=f"{key}_contacts_export")
