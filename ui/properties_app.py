"""Standalone property pages retained for existing property-workbench links."""

import streamlit as st

from bur01_structmine.model_web import navigate

st.set_page_config(page_title="AquaEsterNet · 模型工作台", layout="wide", initial_sidebar_state="collapsed")
navigate()
