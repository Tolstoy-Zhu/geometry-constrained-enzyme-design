@echo off
REM Run this file from Windows Command Prompt inside the project folder.
python -m venv .venv
call .venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
python run_demo_no_foldseek.py
pause
