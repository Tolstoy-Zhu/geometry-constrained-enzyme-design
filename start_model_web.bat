@echo off
setlocal
cd /d "%~dp0"
call "%~dp0start_properties_web.bat"
if errorlevel 1 pause
endlocal
