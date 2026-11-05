[CmdletBinding()]
param(
    [ValidateSet("Install", "Uninstall", "Start", "Stop", "Restart", "Run", "Status")]
    [string]$Action = "Status",
    [ValidateRange(1024, 65535)]
    [int]$Port = 8501,
    [switch]$OpenBrowser,
    [switch]$AsJson
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$TaskName = "AquaEsterNet-StructMine UI"
$ScriptPath = $MyInvocation.MyCommand.Path
$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$RuntimeRoot = Join-Path $ProjectRoot ".cache\ui_service"
$StatePath = Join-Path $RuntimeRoot "state.json"
$StopPath = Join-Path $RuntimeRoot "stop.requested"
$ServiceLog = Join-Path $RuntimeRoot "service.log"
$Url = "http://127.0.0.1:$Port"
$HealthUrl = "$Url/_stcore/health"

function Initialize-RuntimeRoot {
    New-Item -ItemType Directory -Force -Path $RuntimeRoot | Out-Null
}

function Write-ServiceLog {
    param([string]$Message)
    Initialize-RuntimeRoot
    $line = "{0} {1}" -f ([DateTime]::UtcNow.ToString("o")), $Message
    Add-Content -LiteralPath $ServiceLog -Value $line -Encoding UTF8
}

function Resolve-AquaPython {
    $candidates = @(
        $env:AQUA_PYTHON,
        [Environment]::GetEnvironmentVariable("AQUA_PYTHON", "User"),
        (Join-Path $ProjectRoot ".venv\Scripts\python.exe")
    )

    $pythonCommand = Get-Command python.exe -ErrorAction SilentlyContinue
    if ($null -ne $pythonCommand) {
        $candidates += $pythonCommand.Source
    }

    foreach ($candidate in $candidates) {
        if ([string]::IsNullOrWhiteSpace($candidate)) {
            continue
        }
        if (Test-Path -LiteralPath $candidate -PathType Leaf) {
            return (Resolve-Path -LiteralPath $candidate).Path
        }
    }

    throw "No usable Python interpreter was found. Set the user AQUA_PYTHON variable or create .venv."
}

function Test-UiHealth {
    try {
        $response = Invoke-WebRequest -UseBasicParsing -Uri $HealthUrl -TimeoutSec 3
        return $response.StatusCode -eq 200 -and $response.Content.Trim() -eq "ok"
    }
    catch {
        return $false
    }
}

function Get-UiListenerPid {
    $listener = Get-NetTCPConnection -LocalAddress "127.0.0.1" -LocalPort $Port -State Listen -ErrorAction SilentlyContinue |
        Select-Object -First 1
    if ($null -eq $listener) {
        return $null
    }
    return [int]$listener.OwningProcess
}

function Get-ProcessCommandLine {
    param([int]$ProcessId)
    $processInfo = Get-CimInstance Win32_Process -Filter "ProcessId = $ProcessId" -ErrorAction SilentlyContinue
    if ($null -eq $processInfo) {
        return $null
    }
    return [string]$processInfo.CommandLine
}

function Test-ProjectUiProcess {
    param([int]$ProcessId)
    $commandLine = Get-ProcessCommandLine -ProcessId $ProcessId
    if ([string]::IsNullOrWhiteSpace($commandLine)) {
        return $false
    }
    return $commandLine -match "streamlit" -and $commandLine -match "ui[\\/]app\.py"
}

function Test-ServiceProcess {
    param([int]$ProcessId)
    $commandLine = Get-ProcessCommandLine -ProcessId $ProcessId
    if ([string]::IsNullOrWhiteSpace($commandLine)) {
        return $false
    }
    return $commandLine -match [Regex]::Escape($ScriptPath) -and $commandLine -match "-Action\s+Run"
}

function Read-UiState {
    if (-not (Test-Path -LiteralPath $StatePath -PathType Leaf)) {
        return $null
    }
    try {
        return Get-Content -LiteralPath $StatePath -Raw -Encoding UTF8 | ConvertFrom-Json
    }
    catch {
        return $null
    }
}

function Write-UiState {
    param(
        [string]$Status,
        [Nullable[int]]$ServicePid,
        [Nullable[int]]$StreamlitPid,
        [int]$RestartCount,
        [string]$Message = ""
    )

    Initialize-RuntimeRoot
    $state = [ordered]@{
        status = $Status
        healthy = Test-UiHealth
        url = $Url
        port = $Port
        task_name = $TaskName
        project_root = $ProjectRoot
        service_pid = $ServicePid
        streamlit_pid = $StreamlitPid
        restart_count = $RestartCount
        message = $Message
        checked_at = [DateTime]::UtcNow.ToString("o")
    }
    $temporaryPath = "$StatePath.$PID.tmp"
    $state | ConvertTo-Json -Depth 4 | Set-Content -LiteralPath $temporaryPath -Encoding UTF8
    Move-Item -LiteralPath $temporaryPath -Destination $StatePath -Force
}

function Remove-OldUiLogs {
    Get-ChildItem -LiteralPath $RuntimeRoot -Filter "streamlit-*.log" -File -ErrorAction SilentlyContinue |
        Sort-Object LastWriteTime -Descending |
        Select-Object -Skip 20 |
        Remove-Item -Force -ErrorAction SilentlyContinue
}

function Start-StreamlitProcess {
    param([string]$Python)

    Initialize-RuntimeRoot
    Remove-OldUiLogs
    $stamp = Get-Date -Format "yyyyMMdd-HHmmss"
    $stdoutPath = Join-Path $RuntimeRoot "streamlit-$stamp.out.log"
    $stderrPath = Join-Path $RuntimeRoot "streamlit-$stamp.err.log"
    $env:AQUA_PYTHON = $Python
    $env:AQUA_UI_SERVICE = "1"
    $env:AQUA_UI_PORT = [string]$Port

    $arguments = @(
        "-m", "streamlit", "run", "ui/app.py",
        "--server.address=127.0.0.1",
        "--server.port=$Port",
        "--server.headless=true",
        "--server.fileWatcherType=none",
        "--server.enableStaticServing=true",
        "--browser.gatherUsageStats=false"
    )
    return Start-Process -FilePath $Python -ArgumentList $arguments -WorkingDirectory $ProjectRoot `
        -WindowStyle Hidden -RedirectStandardOutput $stdoutPath -RedirectStandardError $stderrPath -PassThru
}

function Start-Watchdog {
    Initialize-RuntimeRoot
    Remove-Item -LiteralPath $StopPath -Force -ErrorAction SilentlyContinue

    $task = Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
    if ($null -ne $task -and $task.State -ne "Disabled") {
        if ($task.State -ne "Running") {
            Start-ScheduledTask -TaskName $TaskName
        }
        return
    }

    $powershellPath = (Get-Process -Id $PID).Path
    $arguments = "-NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File `"$ScriptPath`" -Action Run -Port $Port"
    Start-Process -FilePath $powershellPath -ArgumentList $arguments -WorkingDirectory $ProjectRoot -WindowStyle Hidden | Out-Null
}

function Wait-ForUi {
    param([int]$TimeoutSeconds = 60)
    $deadline = (Get-Date).AddSeconds($TimeoutSeconds)
    while ((Get-Date) -lt $deadline) {
        if (Test-UiHealth) {
            return $true
        }
        Start-Sleep -Milliseconds 500
    }
    return $false
}

function Stop-UiProcesses {
    Initialize-RuntimeRoot
    New-Item -ItemType File -Force -Path $StopPath | Out-Null

    $listenerPid = Get-UiListenerPid
    if ($null -ne $listenerPid -and (Test-ProjectUiProcess -ProcessId $listenerPid)) {
        Stop-Process -Id $listenerPid -Force -ErrorAction SilentlyContinue
    }

    $state = Read-UiState
    if ($null -ne $state -and $null -ne $state.service_pid) {
        $servicePid = [int]$state.service_pid
        if ($servicePid -ne $PID -and (Test-ServiceProcess -ProcessId $servicePid)) {
            Stop-Process -Id $servicePid -Force -ErrorAction SilentlyContinue
        }
    }

    $task = Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
    if ($null -ne $task -and $task.State -eq "Running") {
        Stop-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
    }

    Write-UiState -Status "stopped" -ServicePid $null -StreamlitPid $null -RestartCount 0 -Message "Stopped by user"
    Write-ServiceLog "UI service stopped"
}

function Install-UiTask {
    Initialize-RuntimeRoot
    $powershellPath = (Get-Process -Id $PID).Path
    $taskArguments = "-NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File `"$ScriptPath`" -Action Run -Port $Port"
    $taskAction = New-ScheduledTaskAction -Execute $powershellPath -Argument $taskArguments -WorkingDirectory $ProjectRoot
    $identity = [Security.Principal.WindowsIdentity]::GetCurrent().Name
    $logonTrigger = New-ScheduledTaskTrigger -AtLogOn -User $identity
    $recoveryTrigger = New-ScheduledTaskTrigger -Once -At ((Get-Date).AddMinutes(1)) `
        -RepetitionInterval (New-TimeSpan -Minutes 1) -RepetitionDuration (New-TimeSpan -Days 3650)
    $principal = New-ScheduledTaskPrincipal -UserId $identity -LogonType Interactive -RunLevel Limited
    $settings = New-ScheduledTaskSettingsSet -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries `
        -StartWhenAvailable -RestartCount 999 -RestartInterval (New-TimeSpan -Minutes 1) `
        -ExecutionTimeLimit ([TimeSpan]::Zero) -MultipleInstances IgnoreNew
    $task = New-ScheduledTask -Action $taskAction -Trigger @($logonTrigger, $recoveryTrigger) `
        -Principal $principal -Settings $settings `
        -Description "Keeps the local AquaEsterNet-StructMine Streamlit workbench available on port $Port."
    Register-ScheduledTask -TaskName $TaskName -InputObject $task -Force | Out-Null
    Write-ServiceLog "Scheduled task installed for $identity on port $Port"
}

function Get-UiStatus {
    $task = Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
    $taskInfo = if ($null -ne $task) { Get-ScheduledTaskInfo -TaskName $TaskName -ErrorAction SilentlyContinue } else { $null }
    $state = Read-UiState
    $listenerPid = Get-UiListenerPid
    return [ordered]@{
        healthy = Test-UiHealth
        url = $Url
        port = $Port
        auto_start_installed = $null -ne $task
        task_state = if ($null -ne $task) { [string]$task.State } else { "NotInstalled" }
        last_task_result = if ($null -ne $taskInfo) { $taskInfo.LastTaskResult } else { $null }
        service_pid = if ($null -ne $state) { $state.service_pid } else { $null }
        streamlit_pid = $listenerPid
        restart_count = if ($null -ne $state) { $state.restart_count } else { 0 }
        state_message = if ($null -ne $state) { $state.message } else { "" }
        state_file = $StatePath
        service_log = $ServiceLog
    }
}

function Run-Watchdog {
    Initialize-RuntimeRoot
    if (Test-Path -LiteralPath $StopPath) {
        Write-ServiceLog "Watchdog start skipped because a stop request is present"
        return 0
    }

    $createdNew = $false
    $mutex = New-Object System.Threading.Mutex($true, "Local\AquaEsterNetStructMineUiService", [ref]$createdNew)
    if (-not $createdNew) {
        $mutex.Dispose()
        return 0
    }

    $ownedStreamlitPid = $null
    $previousState = Read-UiState
    $restartCount = if ($null -ne $previousState -and $null -ne $previousState.restart_count) {
        [int]$previousState.restart_count
    }
    else {
        0
    }
    $lastReportedStatus = ""
    try {
        $python = Resolve-AquaPython
        Write-ServiceLog "Watchdog started with Python $python"

        while (-not (Test-Path -LiteralPath $StopPath)) {
            if (Test-UiHealth) {
                $listenerPid = Get-UiListenerPid
                Write-UiState -Status "running" -ServicePid $PID -StreamlitPid $listenerPid `
                    -RestartCount $restartCount -Message "Healthy"
                if ($lastReportedStatus -ne "running") {
                    Write-ServiceLog "UI healthy at $Url (pid=$listenerPid)"
                    $lastReportedStatus = "running"
                }
                Start-Sleep -Seconds 5
                continue
            }

            $listenerPid = Get-UiListenerPid
            if ($null -ne $listenerPid) {
                Write-UiState -Status "blocked" -ServicePid $PID -StreamlitPid $listenerPid `
                    -RestartCount $restartCount -Message "Port $Port is occupied by an unhealthy process"
                if ($lastReportedStatus -ne "blocked") {
                    Write-ServiceLog "Port $Port is occupied but the health endpoint is unavailable"
                    $lastReportedStatus = "blocked"
                }
                Start-Sleep -Seconds 10
                continue
            }

            $restartCount += 1
            Write-UiState -Status "starting" -ServicePid $PID -StreamlitPid $null `
                -RestartCount $restartCount -Message "Starting Streamlit"
            Write-ServiceLog "Starting Streamlit attempt $restartCount"
            $streamlitProcess = Start-StreamlitProcess -Python $python
            $ownedStreamlitPid = $streamlitProcess.Id

            $ready = $false
            for ($attempt = 0; $attempt -lt 60; $attempt++) {
                if (Test-Path -LiteralPath $StopPath) {
                    break
                }
                if (Test-UiHealth) {
                    $ready = $true
                    break
                }
                if ($streamlitProcess.HasExited) {
                    break
                }
                Start-Sleep -Seconds 1
            }

            if ($ready) {
                $lastReportedStatus = ""
                continue
            }

            if (-not $streamlitProcess.HasExited) {
                Stop-Process -Id $streamlitProcess.Id -Force -ErrorAction SilentlyContinue
            }
            Write-UiState -Status "failed" -ServicePid $PID -StreamlitPid $null `
                -RestartCount $restartCount -Message "Streamlit did not become healthy"
            Write-ServiceLog "Streamlit attempt $restartCount failed; retrying"
            Start-Sleep -Seconds ([Math]::Min(30, 5 * $restartCount))
        }

        return 0
    }
    catch {
        Write-ServiceLog "Watchdog failed: $($_.Exception.Message)"
        Write-UiState -Status "failed" -ServicePid $PID -StreamlitPid $ownedStreamlitPid `
            -RestartCount $restartCount -Message $_.Exception.Message
        return 1
    }
    finally {
        if ($null -ne $ownedStreamlitPid -and (Test-ProjectUiProcess -ProcessId $ownedStreamlitPid)) {
            Stop-Process -Id $ownedStreamlitPid -Force -ErrorAction SilentlyContinue
        }
        $mutex.ReleaseMutex()
        $mutex.Dispose()
    }
}

switch ($Action) {
    "Install" {
        Install-UiTask
        Remove-Item -LiteralPath $StopPath -Force -ErrorAction SilentlyContinue
        Start-Watchdog
        if (-not (Wait-ForUi)) {
            throw "The UI service was installed but did not become healthy within 60 seconds."
        }
        if ($OpenBrowser) {
            Start-Process $Url
        }
        Write-Output "Persistent UI installed and running at $Url"
        exit 0
    }
    "Uninstall" {
        Stop-UiProcesses
        Unregister-ScheduledTask -TaskName $TaskName -Confirm:$false -ErrorAction SilentlyContinue
        Write-ServiceLog "Scheduled task removed"
        Write-Output "Persistent UI auto-start removed. Project files were not changed."
        exit 0
    }
    "Start" {
        Start-Watchdog
        if (-not (Wait-ForUi)) {
            throw "The UI did not become healthy within 60 seconds. Check $ServiceLog"
        }
        if ($OpenBrowser) {
            Start-Process $Url
        }
        Write-Output "UI ready at $Url"
        exit 0
    }
    "Stop" {
        Stop-UiProcesses
        Write-Output "UI stopped. Double-click start_ui.bat to start it again."
        exit 0
    }
    "Restart" {
        Stop-UiProcesses
        Start-Sleep -Seconds 2
        Start-Watchdog
        if (-not (Wait-ForUi)) {
            throw "The UI did not become healthy after restart. Check $ServiceLog"
        }
        if ($OpenBrowser) {
            Start-Process $Url
        }
        Write-Output "UI restarted at $Url"
        exit 0
    }
    "Run" {
        exit (Run-Watchdog)
    }
    "Status" {
        $status = Get-UiStatus
        if ($AsJson) {
            $status | ConvertTo-Json -Depth 4
        }
        else {
            Write-Output "AquaEsterNet-StructMine persistent UI"
            Write-Output ("Healthy: {0}" -f $status.healthy)
            Write-Output ("URL: {0}" -f $status.url)
            Write-Output ("Auto-start installed: {0}" -f $status.auto_start_installed)
            Write-Output ("Task state: {0}" -f $status.task_state)
            Write-Output ("Service PID: {0}" -f $status.service_pid)
            Write-Output ("Streamlit PID: {0}" -f $status.streamlit_pid)
            Write-Output ("Automatic restarts: {0}" -f $status.restart_count)
            Write-Output ("Log: {0}" -f $status.service_log)
        }
        exit 0
    }
}
