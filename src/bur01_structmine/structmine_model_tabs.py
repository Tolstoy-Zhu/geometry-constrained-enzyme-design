"""Model actions in the original StructureMine tabs, using the current mining session."""

import streamlit as st

from .design_view import tr
from .model_web import keep_property_draft, prediction_form, prediction_history
from .prediction_sequence import safe_sequence_name
from .property_fasta import parse_property_fasta


def open_model_tab(section):
    # Apply the request before the native tab widgets are created on the next run.
    st.session_state["_structmine_destination"] = section
    st.rerun()


def workbench_tabs(labels, after):
    design = tr("序列设计", "Sequence design")
    properties = tr("性质预测", "Property prediction")
    split = labels.index(after) + 1
    names = [*labels[:split], design, properties, *labels[split:]]
    keep_property_draft()
    destination = st.session_state.pop("_structmine_destination", None)
    if "structmine_tab" not in st.session_state and st.query_params.get("run"):
        destination = destination or "results"
    if destination:
        st.session_state["structmine_tab"] = design if destination == "design" else properties
        if destination in {"predict", "results"}:
            st.session_state["native_property_mode"] = destination
    elif st.session_state.get("structmine_tab", names[0]) not in names:
        st.session_state["structmine_tab"] = names[0]
    tabs = st.tabs(names, key="structmine_tab", on_change="rerun")
    design_tab, property_tab = tabs[split : split + 2]
    if not property_tab.open and "run" in st.query_params:
        del st.query_params["run"]
    return [*tabs[:split], *tabs[split + 2 :]], design_tab, property_tab


def show_property_tab():
    modes = {"predict": tr("提交预测", "Submit prediction"), "results": tr("任务与结果", "Jobs and results")}
    st.session_state.setdefault("native_property_mode", "predict")
    mode = st.segmented_control(
        tr("性质预测", "Property prediction"),
        list(modes),
        format_func=modes.__getitem__,
        required=True,
        key="native_property_mode",
    )
    if mode == "results":
        prediction_history(open_model_tab)
    else:
        if "run" in st.query_params:
            del st.query_params["run"]
        prediction_form(open_model_tab)


def candidate_property_button(column, name, sequence, key):
    fasta = ""
    hint = tr(
        "这条候选尚无完整序列，可在性质预测页上传 FASTA。", "Supply a full sequence in Property prediction."
    )
    if sequence:
        try:
            records = parse_property_fasta(f">{safe_sequence_name(name)}\n{sequence}")
            if len(records) != 1:
                raise ValueError(
                    tr("这条候选包含多条序列，请核对原记录。", "This candidate contains multiple sequences.")
                )
            record = records[0]
            fasta = record.to_fasta()
            hint = f"{name} · {len(record.sequence)} aa"
        except ValueError as exc:
            hint = str(exc)
    if column.button(
        tr("送入性质预测", "Send to property prediction"),
        key=key,
        disabled=not fasta,
        help=hint,
        width="stretch",
    ):
        st.session_state.pop("_web_use_design", None)
        st.session_state["_web_use_sequence"] = fasta
        open_model_tab("predict")
