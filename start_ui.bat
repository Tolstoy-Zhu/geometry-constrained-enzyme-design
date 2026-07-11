@echo off
setlocal
cd /d "%~dp0"
if not defined AQUA_PYTHON if exist "%~dp0.venv\Scripts\python.exe" set "AQUA_PYTHON=%~dp0.venv\Scripts\python.exe"
if not defined AQUA_PYTHON set "AQUA_PYTHON=python"
"%AQUA_PYTHON%" start_here.py --action ui
if errorlevel 1 pause
endlocal
