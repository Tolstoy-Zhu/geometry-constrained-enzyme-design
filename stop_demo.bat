@echo off
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0scripts\windows\demo_service.ps1" -Action stop
if errorlevel 1 pause
