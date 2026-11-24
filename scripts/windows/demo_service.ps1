param([ValidateSet('start', 'stop')][string]$Action = 'start')
$ErrorActionPreference = 'Stop'
$Root = [IO.Path]::GetFullPath((Join-Path $PSScriptRoot '..\..'))
$Candidates = @($env:AQUA_PYTHON, [Environment]::GetEnvironmentVariable('AQUA_PYTHON', 'User'), (Join-Path $Root '.venv\Scripts\python.exe'))
$Python = $Candidates | Where-Object { $_ -and (Test-Path -LiteralPath $_) } | Select-Object -First 1
if (-not $Python) { throw 'Python runtime not found. Set AQUA_PYTHON to the project interpreter.' }
$env:PYTHONPATH = Join-Path $Root 'src'
Push-Location -LiteralPath $Root
try {
    $Command = if ($Action -eq 'start') { 'demo' } else { 'stop-demo' }
    & $Python -m bur01_structmine.cli $Command
    exit $LASTEXITCODE
} finally { Pop-Location }
