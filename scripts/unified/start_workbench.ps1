[CmdletBinding()]
param([int]$Port = 8501, [switch]$NoBrowser)
Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"
$WorkspaceRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot "../..")).Path
$RuntimeConfig = Join-Path $WorkspaceRoot ".cache/unified/runtime.json"
if (-not (Test-Path -LiteralPath $RuntimeConfig)) {
    throw "Create .cache/unified/runtime.json using configs/unified/runtime.example.json."
}
$Config = Get-Content -LiteralPath $RuntimeConfig -Raw -Encoding UTF8 | ConvertFrom-Json
$PythonPath = $Config.python
if (-not (Test-Path -LiteralPath $PythonPath -PathType Leaf)) { throw "Configured Python does not exist." }
$Allowed = @("AQUA_PROPERTY_MODELS","AQUA_PROPERTY_PYTHON","AQUA_DESIGN_MODELS","AQUA_DESIGN_PYTHON",
    "AQUA_LEARNING_PYTHON","AQUA_MECHANISM_PYTHON","AQUA_FEATURE_PYTHON","AQUA_PREPARATION_ADAPTER")
foreach ($Property in $Config.environment.PSObject.Properties) {
    if ($Allowed -notcontains $Property.Name) { throw "Unsupported runtime setting: $($Property.Name)" }
    [Environment]::SetEnvironmentVariable($Property.Name, [string]$Property.Value, "Process")
}
$env:PYTHONPATH = Join-Path $WorkspaceRoot "src"
$env:PYTHONIOENCODING = "utf-8"
$AppPath = Join-Path $WorkspaceRoot "ui/app.py"
$Url = "http://127.0.0.1:$Port"
$Listener = Get-NetTCPConnection -LocalAddress "127.0.0.1" -LocalPort $Port -State Listen -ErrorAction SilentlyContinue | Select-Object -First 1
if ($null -ne $Listener) {
    $Running = Get-CimInstance Win32_Process -Filter "ProcessId = $($Listener.OwningProcess)"
    if ($Running.CommandLine -notmatch [Regex]::Escape($AppPath.Replace("\","/")) -and
        $Running.CommandLine -notmatch [Regex]::Escape($AppPath.Replace("/","\"))) {
        throw "The selected port belongs to another application or checkout. Use -Port to choose a different port."
    }
} else {
    $LogRoot = Join-Path $WorkspaceRoot ".cache/unified/web"
    New-Item -ItemType Directory -Force -Path $LogRoot | Out-Null
    $Stamp = Get-Date -Format "yyyyMMdd-HHmmss"
    $Arguments = @("-m","streamlit","run",('"' + $AppPath + '"'),
        "--server.address=127.0.0.1","--server.port=$Port","--server.headless=true",
        "--server.fileWatcherType=none","--server.enableStaticServing=true","--browser.gatherUsageStats=false")
    $Process = Start-Process -FilePath $PythonPath -ArgumentList $Arguments -WorkingDirectory $WorkspaceRoot -WindowStyle Hidden -PassThru `
        -RedirectStandardOutput (Join-Path $LogRoot "$Stamp.out.log") -RedirectStandardError (Join-Path $LogRoot "$Stamp.err.log")
    @{pid=$Process.Id;url=$Url;root=$WorkspaceRoot} | ConvertTo-Json | Set-Content -LiteralPath (Join-Path $LogRoot "process.json") -Encoding UTF8
}
$Ready = $false
for ($Attempt = 0; $Attempt -lt 30; $Attempt++) {
    try {
        $Health = Invoke-WebRequest -UseBasicParsing -Uri "$Url/_stcore/health" -TimeoutSec 1
        if ($Health.StatusCode -eq 200 -and $Health.Content.Trim() -eq "ok") { $Ready = $true; break }
    } catch {}
    Start-Sleep -Milliseconds 500
}
if (-not $Ready) { throw "Workbench startup failed. Check .cache/unified/web logs." }
if (-not $NoBrowser) { Start-Process $Url }
Write-Output "AquaEsterNet workbench ready: $Url"
