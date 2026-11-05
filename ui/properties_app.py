"""Local AquaEsterNet browser application. Launch with 'aquaesternet-structmine web'."""

import streamlit as st

from bur01_structmine.model_web import navigate

st.set_page_config(page_title="AquaEsterNet · 模型工作台", layout="wide", initial_sidebar_state="collapsed")
st.session_state["_aqua_local_shell"] = True
navigate(standalone=True)
