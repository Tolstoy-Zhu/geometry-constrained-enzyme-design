from argparse import ArgumentParser

import pytest
from streamlit.testing.v1 import AppTest

from bur01_structmine.property_cli import add_property_commands


def candidate_workbench():
    import streamlit as st

    from bur01_structmine.structmine_model_tabs import (
        candidate_property_button,
        show_property_tab,
        workbench_tabs,
    )

    original, _, properties = workbench_tabs(["候选", "结构", "历史"], after="结构")
    with original[0]:
        candidate_property_button(st, "hit-1", st.session_state["example_sequence"], "handoff")
    if properties.open:
        with properties:
            show_property_tab()


def test_mined_sequence_replaces_old_parent_draft(tmp_path, monkeypatch):
    monkeypatch.setenv("AQUA_DESIGN_ROOT", str(tmp_path / "designs"))
    test = AppTest.from_function(candidate_workbench, default_timeout=30)
    test.session_state["example_sequence"] = "ACDE" * 76
    test.session_state["web_fasta"] = ">previous\nAGDE"
    test.session_state["web_parent"] = True
    test.session_state["web_source"] = "设计批次"
    test.run()
    assert [tab.label for tab in test.tabs] == ["候选", "结构", "序列设计", "性质预测", "历史"]
    assert not test.text_area
    test.button(key="handoff").click().run()
    assert not test.exception
    fasta = test.text_area(key="web_fasta").value.splitlines()
    assert fasta[0] == ">hit-1"
    assert "".join(fasta[1:]) == "ACDE" * 76
    assert test.session_state["web_source"] == "序列输入"
    assert not test.checkbox(key="web_parent").value
    assert test.session_state["_web_upload_revision"] == 1
    assert any("304 aa" in line.value for line in test.caption)


@pytest.mark.parametrize("sequence", ["", "ACXDE", "ACDE\n>other\nAAAA"])
def test_only_one_supported_mining_sequence_can_be_transferred(sequence):
    test = AppTest.from_function(candidate_workbench, default_timeout=30)
    test.session_state["example_sequence"] = sequence
    test.run()
    assert not test.exception
    assert test.button(key="handoff").disabled


def test_web_launcher_opens_the_original_app_and_keeps_the_existing_property_route():
    parser = ArgumentParser()
    add_property_commands(parser.add_subparsers(dest="command"))
    native = parser.parse_args(["web"])
    previous = parser.parse_args(["properties-web"])
    assert (native.web_app, native.port) == ("app.py", 8501)
    assert (previous.web_app, previous.port) == ("properties_app.py", 8765)
