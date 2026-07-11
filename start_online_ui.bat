@echo off
setlocal
cd /d "%~dp0"
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0scripts\windows\remote_ui_service.ps1" -Action Start -OpenBrowser
if errorlevel 1 (
  pause
  exit /b 1
)
start "" "%~dp0.cache\remote_ui\SHARE_WITH_EXPERTS.txt"
endlocal
