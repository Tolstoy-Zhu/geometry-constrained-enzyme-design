@echo off
setlocal
cd /d "%~dp0"
powershell.exe -NoProfile -Command "$p=$env:AQUA_PYTHON; if (-not $p) { $p=[Environment]::GetEnvironmentVariable('AQUA_PYTHON','User') }; if (-not $p -and (Test-Path -LiteralPath '.venv\Scripts\python.exe')) { $p=(Resolve-Path -LiteralPath '.venv\Scripts\python.exe').Path }; if (-not $p) { $p='python' }; & $p -m bur01_structmine.cli properties-web; exit $LASTEXITCODE"
exit /b %ERRORLEVEL%
