[CmdletBinding()]
param(
    [string]$Python = $(if ($env:AQUA_PYTHON) { $env:AQUA_PYTHON } else { "python" }),
    [switch]$SkipWslSetup,
    [switch]$SkipSwissProtDownload
)

$ErrorActionPreference = "Stop"
$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
Set-Location $ProjectRoot

function Invoke-Checked([string]$Executable, [string[]]$Arguments) {
    Write-Host "[RUN] $Executable $($Arguments -join ' ')"
    & $Executable @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed with exit code ${LASTEXITCODE}: $Executable $($Arguments -join ' ')"
    }
}

if (-not (Test-Path -LiteralPath $Python -PathType Leaf)) {
    throw "Python environment not found: $Python"
}

if (-not $SkipWslSetup) {
    & (Join-Path $PSScriptRoot "setup_wsl_foldseek.ps1")
    if ($LASTEXITCODE -ne 0) {
        throw "WSL/Foldseek setup failed with exit code $LASTEXITCODE"
    }
}

Invoke-Checked $Python @("-m", "pytest", "-q")
Invoke-Checked $Python @("scripts\prepare_benchmark.py")
Invoke-Checked $Python @("-m", "bur01_structmine.cli", "benchmark", "--engine", "foldseek", "--no-prepare")

if (-not $SkipSwissProtDownload) {
    Invoke-Checked $Python @(
        "prepare_databases.py",
        "--config", "configs\alphafold_swissprot.yaml",
        "--only", "alphafold_swissprot"
    )
}
Invoke-Checked $Python @(
    "scripts\write_database_manifest.py",
    "--config", "configs\alphafold_swissprot.yaml",
    "--source", "alphafold_swissprot"
)

Invoke-Checked $Python @(
    "-m", "bur01_structmine.cli",
    "production-run",
    "--config", "configs\alphafold_swissprot.yaml"
)

$latest = Get-Content -Raw -LiteralPath "results\runs\latest.json" | ConvertFrom-Json
$manifest = Get-Content -Raw -LiteralPath $latest.manifest | ConvertFrom-Json
$annotated = $manifest.artifacts.annotated_xlsx
if (-not $annotated -or -not (Test-Path -LiteralPath $annotated -PathType Leaf)) {
    throw "Production annotated workbook is missing from the latest run manifest."
}
Invoke-Checked $Python @(
    "scripts\validate_real_search_outputs.py",
    "--annotated", $annotated,
    "--require-source", "alphafold_swissprot",
    "--require-final-sequence"
)

Write-Host "[OK] Full Foldseek acceptance completed."
Write-Host "[OK] Run manifest: $($latest.manifest)"
Write-Host "[OK] Annotated workbook: $annotated"
