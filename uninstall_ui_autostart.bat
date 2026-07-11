@echo off
setlocal
cd /d "%~dp0"
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0scripts\windows\ui_service.ps1" -Action Uninstall
if errorlevel 1 pause
endlocal
