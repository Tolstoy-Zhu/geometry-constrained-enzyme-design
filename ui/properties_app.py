"""Local interactive property workbench. Launch with 'aquaesternet-structmine properties-web'."""

import streamlit as st

from bur01_structmine.property_view import show_properties

st.set_page_config(page_title="AquaEsterNet · 本地性质预测", layout="wide", initial_sidebar_state="collapsed")
st.title("AquaEsterNet 本地性质预测")
show_properties()
