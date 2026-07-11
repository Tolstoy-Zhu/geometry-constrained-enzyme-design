[CmdletBinding()]
param(
    [ValidateSet("Install", "Uninstall", "Start", "Stop", "Restart", "Run", "Status")]
    [string]$Action = "Status",
    [ValidateSet("Quick", "Token")]
    [string]$Mode = "Quick",
    [string]$PublicUrl = "",
    [switch]$OpenBrowser,
    [switch]$ResetPassword,
    [switch]$AsJson
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$TaskName = "AquaEsterNet-StructMine Remote UI"
$ScriptPath = $MyInvocation.MyCommand.Path
$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$UiServiceScript = Join-Path $PSScriptRoot "ui_service.ps1"
$RuntimeRoot = Join-Path $ProjectRoot ".cache\remote_ui"
$ConfigPath = Join-Path $RuntimeRoot "config.json"
$AccessPath = Join-Path $RuntimeRoot "access.json"
$ShareJsonPath = Join-Path $RuntimeRoot "share.json"
$ShareTextPath = Join-Path $RuntimeRoot "SHARE_WITH_EXPERTS.txt"
$StatePath = Join-Path $RuntimeRoot "state.json"
$StopPath = Join-Path $RuntimeRoot "stop.requested"
$ServiceLog = Join-Path $RuntimeRoot "remote_service.log"
$OriginUrl = "http://127.0.0.1:8501"

function Initialize-RuntimeRoot {
    New-Item -ItemType Directory -Force -Path $RuntimeRoot | Out-Null
}

function Write-ServiceLog {
    param([string]$Message)
    Initialize-RuntimeRoot
    Add-Content -LiteralPath $ServiceLog -Encoding UTF8 `
        -Value ("{0} {1}" -f [DateTime]::UtcNow.ToString("o"), $Message)
}

function Resolve-Cloudflared {
    $candidates = @(
        $env:AQUA_CLOUDFLARED,
        [Environment]::GetEnvironmentVariable("AQUA_CLOUDFLARED", "User"),
        (Join-Path $ProjectRoot ".tools\cloudflared.exe")
    )
    $command = Get-Command cloudflared.exe -ErrorAction SilentlyContinue
    if ($null -ne $command) {
        $candidates += $command.Source
    }
    foreach ($candidate in $candidates) {
        if (-not [string]::IsNullOrWhiteSpace($candidate) -and (Test-Path -LiteralPath $candidate -PathType Leaf)) {
            return (Resolve-Path -LiteralPath $candidate).Path
        }
    }
    throw "cloudflared was not found. Set AQUA_CLOUDFLARED, add it to PATH, or place it under .tools."
}

function Read-JsonFile {
    param([string]$Path)
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
        return $null
    }
    try {
        return Get-Content -LiteralPath $Path -Raw -Encoding UTF8 | ConvertFrom-Json
    }
    catch {
        return $null
    }
}

function Write-JsonFile {
    param([string]$Path, [object]$Value)
    Initialize-RuntimeRoot
    $temporaryPath = "$Path.$PID.tmp"
    $Value | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $temporaryPath -Encoding UTF8
    Move-Item -LiteralPath $temporaryPath -Destination $Path -Force
}

function Protect-SecretFile {
    param([string]$Path)
    try {
        $identity = [Security.Principal.WindowsIdentity]::GetCurrent().Name
        $acl = Get-Acl -LiteralPath $Path
        $acl.SetAccessRuleProtection($true, $false)
        $rule = New-Object Security.AccessControl.FileSystemAccessRule($identity, "FullControl", "Allow")
        $acl.SetAccessRule($rule)
        Set-Acl -LiteralPath $Path -AclObject $acl
    }
    catch {
        Write-ServiceLog "Could not restrict ACL for ${Path}: $($_.Exception.Message)"
    }
}

function New-RemotePassword {
    $bytes = New-Object byte[] 24
    $rng = [Security.Cryptography.RandomNumberGenerator]::Create()
    try {
        $rng.GetBytes($bytes)
    }
    finally {
        $rng.Dispose()
    }
    $token = [Convert]::ToBase64String($bytes).TrimEnd("=").Replace("+", "-").Replace("/", "_")
    return "Aqua-$token"
}

function New-AccessRecord {
    param([string]$Password, [string]$TunnelMode, [string]$RemoteUrl)
    $iterations = 600000
    $salt = New-Object byte[] 16
    $rng = [Security.Cryptography.RandomNumberGenerator]::Create()
    try {
        $rng.GetBytes($salt)
    }
    finally {
        $rng.Dispose()
    }
    $derive = New-Object Security.Cryptography.Rfc2898DeriveBytes(
        $Password,
        $salt,
        $iterations,
        [Security.Cryptography.HashAlgorithmName]::SHA256
    )
    try {
        $digest = $derive.GetBytes(32)
    }
    finally {
        $derive.Dispose()
    }
    return [ordered]@{
        enabled = $true
        scheme = "pbkdf2-sha256"
        iterations = $iterations
        salt_b64 = [Convert]::ToBase64String($salt)
        hash_b64 = [Convert]::ToBase64String($digest)
        mode = $TunnelMode
        public_url = $RemoteUrl
        created_at = [DateTime]::UtcNow.ToString("o")
    }
}

function Write-ShareFiles {
    param([string]$Password, [string]$TunnelMode, [string]$RemoteUrl)
    $share = [ordered]@{
        public_url = $RemoteUrl
        access_password = $Password
        mode = $TunnelMode
        updated_at = [DateTime]::UtcNow.ToString("o")
    }
    Write-JsonFile -Path $ShareJsonPath -Value $share
    $lines = @(
        "AquaEsterNet-StructMine expert demonstration access",
        "====================================================",
        "",
        "Online URL: $RemoteUrl",
        "Access password: $Password",
        "",
        "Share these two lines only with invited experts.",
        "The owner computer must remain powered on, signed in, and connected to the Internet.",
        "Quick Tunnel URLs may change after the tunnel restarts. This file is updated automatically.",
        "Do not commit or upload this file."
    )
    $lines | Set-Content -LiteralPath $ShareTextPath -Encoding UTF8
    Protect-SecretFile -Path $ShareJsonPath
    Protect-SecretFile -Path $ShareTextPath
}

function Enable-RemoteAccess {
    param([string]$TunnelMode, [string]$RemoteUrl, [switch]$ForceNewPassword)

    $share = Read-JsonFile -Path $ShareJsonPath
    $access = Read-JsonFile -Path $AccessPath
    if ($ForceNewPassword -or $null -eq $share -or [string]::IsNullOrWhiteSpace([string]$share.access_password) -or $null -eq $access) {
        $password = New-RemotePassword
        $access = New-AccessRecord -Password $password -TunnelMode $TunnelMode -RemoteUrl $RemoteUrl
    }
    else {
        $password = [string]$share.access_password
        $access.enabled = $true
        $access.mode = $TunnelMode
        $access.public_url = $RemoteUrl
    }
    Write-JsonFile -Path $AccessPath -Value $access
    Write-ShareFiles -Password $password -TunnelMode $TunnelMode -RemoteUrl $RemoteUrl
    return $password
}

function Disable-RemoteAccess {
    $access = Read-JsonFile -Path $AccessPath
    if ($null -ne $access) {
        $access.enabled = $false
        $access | Add-Member -NotePropertyName "disabled_at" `
            -NotePropertyValue ([DateTime]::UtcNow.ToString("o")) -Force
        Write-JsonFile -Path $AccessPath -Value $access
    }
}

function Write-RemoteState {
    param(
        [string]$Status,
        [Nullable[int]]$ServicePid,
        [Nullable[int]]$TunnelPid,
        [int]$RestartCount,
        [string]$TunnelMode,
        [string]$RemoteUrl,
        [string]$Message = "",
        [string]$StdoutPath = "",
        [string]$StderrPath = ""
    )
    $state = [ordered]@{
        status = $Status
        local_healthy = Test-LocalUiHealth
        remote_healthy = if ([string]::IsNullOrWhiteSpace($RemoteUrl)) { $false } else { Test-RemoteHealth -RemoteUrl $RemoteUrl }
        mode = $TunnelMode
        public_url = $RemoteUrl
        origin_url = $OriginUrl
        service_pid = $ServicePid
        tunnel_pid = $TunnelPid
        restart_count = $RestartCount
        message = $Message
        stdout_path = $StdoutPath
        stderr_path = $StderrPath
        share_file = $ShareTextPath
        checked_at = [DateTime]::UtcNow.ToString("o")
    }
    Write-JsonFile -Path $StatePath -Value $state
}

function Test-LocalUiHealth {
    try {
        $response = Invoke-WebRequest -UseBasicParsing -Uri "$OriginUrl/_stcore/health" -TimeoutSec 3
        return $response.StatusCode -eq 200 -and $response.Content.Trim() -eq "ok"
    }
    catch {
        return $false
    }
}

function Test-RemoteHealth {
    param([string]$RemoteUrl)
    if ([string]::IsNullOrWhiteSpace($RemoteUrl)) {
        return $false
    }
    $handler = $null
    $client = $null
    try {
        Add-Type -AssemblyName System.Net.Http
        $handler = New-Object System.Net.Http.HttpClientHandler
        $handler.UseProxy = $false
        $client = New-Object System.Net.Http.HttpClient($handler)
        $client.Timeout = [TimeSpan]::FromSeconds(10)
        for ($attempt = 0; $attempt -lt 3; $attempt++) {
            $response = $null
            try {
                $response = $client.GetAsync("$($RemoteUrl.TrimEnd('/'))/_stcore/health").GetAwaiter().GetResult()
                $content = $response.Content.ReadAsStringAsync().GetAwaiter().GetResult()
                if ($response.IsSuccessStatusCode -and $content.Trim() -eq "ok") {
                    return $true
                }
            }
            catch {
                if ($attempt -eq 2) {
                    return $false
                }
            }
            finally {
                if ($null -ne $response) {
                    $response.Dispose()
                }
            }
            Start-Sleep -Milliseconds 500
        }
        return $false
    }
    catch {
        return $false
    }
    finally {
        if ($null -ne $client) {
            $client.Dispose()
        }
        if ($null -ne $handler) {
            $handler.Dispose()
        }
    }
}

function Get-ProcessCommandLine {
    param([int]$ProcessId)
    $info = Get-CimInstance Win32_Process -Filter "ProcessId = $ProcessId" -ErrorAction SilentlyContinue
    if ($null -eq $info) {
        return $null
    }
    return [string]$info.CommandLine
}

function Test-TunnelProcess {
    param([int]$ProcessId)
    $commandLine = Get-ProcessCommandLine -ProcessId $ProcessId
    return -not [string]::IsNullOrWhiteSpace($commandLine) -and $commandLine -match "cloudflared" -and $commandLine -match "tunnel"
}

function Test-RemoteServiceProcess {
    param([int]$ProcessId)
    $commandLine = Get-ProcessCommandLine -ProcessId $ProcessId
    return -not [string]::IsNullOrWhiteSpace($commandLine) -and `
        $commandLine -match [Regex]::Escape($ScriptPath) -and $commandLine -match "-Action\s+Run"
}

function Get-TunnelToken {
    param([object]$Config)
    if ($null -eq $Config.token_ciphertext -or [string]::IsNullOrWhiteSpace([string]$Config.token_ciphertext)) {
        throw "The named tunnel token is not configured."
    }
    $secure = ConvertTo-SecureString ([string]$Config.token_ciphertext)
    $credential = New-Object Management.Automation.PSCredential("token", $secure)
    return $credential.GetNetworkCredential().Password
}

function Write-RemoteConfig {
    param([string]$TunnelMode, [string]$RemoteUrl)
    $cloudflared = Resolve-Cloudflared
    $tokenCiphertext = ""
    if ($TunnelMode -eq "Token") {
        $token = if (-not [string]::IsNullOrWhiteSpace($env:AQUA_CLOUDFLARED_TOKEN)) {
            $env:AQUA_CLOUDFLARED_TOKEN
        }
        else {
            [Environment]::GetEnvironmentVariable("AQUA_CLOUDFLARED_TOKEN", "User")
        }
        if ([string]::IsNullOrWhiteSpace($token)) {
            throw "Set the user AQUA_CLOUDFLARED_TOKEN variable before installing a named tunnel."
        }
        if ($RemoteUrl -notmatch '^https://[A-Za-z0-9.-]+/?$') {
            throw "Token mode requires -PublicUrl with the HTTPS hostname configured in Cloudflare."
        }
        $tokenCiphertext = ConvertFrom-SecureString (ConvertTo-SecureString $token -AsPlainText -Force)
    }
    $config = [ordered]@{
        mode = $TunnelMode
        public_url = $RemoteUrl.TrimEnd("/")
        origin_url = $OriginUrl
        cloudflared_path = $cloudflared
        token_ciphertext = $tokenCiphertext
        cloudflared_version = (& $cloudflared --version | Out-String).Trim()
        installed_at = [DateTime]::UtcNow.ToString("o")
    }
    Write-JsonFile -Path $ConfigPath -Value $config
    Protect-SecretFile -Path $ConfigPath
}

function Start-CloudflaredProcess {
    param([object]$Config, [int]$Attempt)
    $stamp = Get-Date -Format "yyyyMMdd-HHmmss"
    $stdoutPath = Join-Path $RuntimeRoot "cloudflared-$stamp-$Attempt.out.log"
    $stderrPath = Join-Path $RuntimeRoot "cloudflared-$stamp-$Attempt.err.log"
    if ([string]$Config.mode -eq "Quick") {
        $arguments = @("tunnel", "--url", $OriginUrl, "--no-autoupdate", "--loglevel", "info")
    }
    else {
        $token = Get-TunnelToken -Config $Config
        $arguments = @("tunnel", "--no-autoupdate", "run", "--token", $token)
    }
    $process = Start-Process -FilePath ([string]$Config.cloudflared_path) -ArgumentList $arguments `
        -WorkingDirectory $ProjectRoot -WindowStyle Hidden -RedirectStandardOutput $stdoutPath `
        -RedirectStandardError $stderrPath -PassThru
    return [pscustomobject]@{ Process = $process; StdoutPath = $stdoutPath; StderrPath = $stderrPath }
}

function Find-QuickTunnelUrl {
    param([string[]]$Paths)
    foreach ($path in $Paths) {
        if ([string]::IsNullOrWhiteSpace($path) -or -not (Test-Path -LiteralPath $path -PathType Leaf)) {
            continue
        }
        $content = Get-Content -LiteralPath $path -Raw -Encoding UTF8 -ErrorAction SilentlyContinue
        if ([string]::IsNullOrWhiteSpace([string]$content)) {
            continue
        }
        $match = [Regex]::Match([string]$content, 'https://[a-z0-9-]+\.trycloudflare\.com')
        if ($match.Success) {
            return $match.Value
        }
    }
    return ""
}

function Wait-ForTunnelUrl {
    param([object]$Config, [object]$Tunnel, [int]$TimeoutSeconds = 75)
    if ([string]$Config.mode -eq "Token") {
        return [string]$Config.public_url
    }
    $deadline = (Get-Date).AddSeconds($TimeoutSeconds)
    while ((Get-Date) -lt $deadline) {
        $url = Find-QuickTunnelUrl -Paths @($Tunnel.StdoutPath, $Tunnel.StderrPath)
        if (-not [string]::IsNullOrWhiteSpace($url)) {
            return $url
        }
        if ($Tunnel.Process.HasExited) {
            break
        }
        Start-Sleep -Seconds 1
    }
    return ""
}

function Install-RemoteTask {
    $powershellPath = (Get-Process -Id $PID).Path
    $arguments = "-NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File `"$ScriptPath`" -Action Run"
    $taskAction = New-ScheduledTaskAction -Execute $powershellPath -Argument $arguments -WorkingDirectory $ProjectRoot
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
        -Description "Keeps the authenticated AquaEsterNet-StructMine HTTPS demonstration tunnel available."
    Register-ScheduledTask -TaskName $TaskName -InputObject $task -Force | Out-Null
    Write-ServiceLog "Remote scheduled task installed for $identity"
}

function Start-RemoteWatchdog {
    Initialize-RuntimeRoot
    Remove-Item -LiteralPath $StopPath -Force -ErrorAction SilentlyContinue
    & powershell.exe -NoProfile -ExecutionPolicy Bypass -File $UiServiceScript -Action Start | Out-Null
    if ($LASTEXITCODE -ne 0) {
        throw "The local UI service could not be started."
    }
    $task = Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
    if ($null -ne $task) {
        if ($task.State -ne "Running") {
            Start-ScheduledTask -TaskName $TaskName
        }
        return
    }
    $powershellPath = (Get-Process -Id $PID).Path
    $arguments = "-NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File `"$ScriptPath`" -Action Run"
    Start-Process -FilePath $powershellPath -ArgumentList $arguments -WorkingDirectory $ProjectRoot -WindowStyle Hidden | Out-Null
}

function Wait-ForRemote {
    param([int]$TimeoutSeconds = 120)
    $deadline = (Get-Date).AddSeconds($TimeoutSeconds)
    while ((Get-Date) -lt $deadline) {
        $state = Read-JsonFile -Path $StatePath
        if ($null -ne $state -and -not [string]::IsNullOrWhiteSpace([string]$state.public_url) -and `
            (Test-RemoteHealth -RemoteUrl ([string]$state.public_url))) {
            return $state
        }
        Start-Sleep -Seconds 1
    }
    return $null
}

function Stop-RemoteProcesses {
    Initialize-RuntimeRoot
    New-Item -ItemType File -Force -Path $StopPath | Out-Null
    Disable-RemoteAccess
    $state = Read-JsonFile -Path $StatePath
    Get-CimInstance Win32_Process -Filter "Name = 'cloudflared.exe'" -ErrorAction SilentlyContinue |
        Where-Object { $_.CommandLine -match "tunnel" -and $_.CommandLine -match [Regex]::Escape($OriginUrl) } |
        ForEach-Object { Stop-Process -Id $_.ProcessId -Force -ErrorAction SilentlyContinue }
    if ($null -ne $state -and $null -ne $state.service_pid) {
        $servicePid = [int]$state.service_pid
        if ($servicePid -ne $PID -and (Test-RemoteServiceProcess -ProcessId $servicePid)) {
            Stop-Process -Id $servicePid -Force -ErrorAction SilentlyContinue
        }
    }
    $task = Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
    if ($null -ne $task -and $task.State -eq "Running") {
        Stop-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
    }
    $modeValue = if ($null -ne $state) { [string]$state.mode } else { "" }
    $urlValue = if ($null -ne $state) { [string]$state.public_url } else { "" }
    Write-RemoteState -Status "stopped" -ServicePid $null -TunnelPid $null -RestartCount 0 `
        -TunnelMode $modeValue -RemoteUrl $urlValue -Message "Stopped by user"
    Write-ServiceLog "Remote UI stopped"
}

function Run-RemoteWatchdog {
    Initialize-RuntimeRoot
    if (Test-Path -LiteralPath $StopPath) {
        return 0
    }
    $createdNew = $false
    $mutex = New-Object Threading.Mutex($true, "Local\AquaEsterNetStructMineRemoteUi", [ref]$createdNew)
    if (-not $createdNew) {
        $mutex.Dispose()
        return 0
    }
    $tunnelProcess = $null
    $stdoutPath = ""
    $stderrPath = ""
    $remoteUrl = ""
    $restartCount = 0
    $config = $null
    $healthGraceUntil = Get-Date
    try {
        $config = Read-JsonFile -Path $ConfigPath
        if ($null -eq $config) {
            throw "Remote configuration is missing. Run install_online_ui.bat."
        }
        $previousState = Read-JsonFile -Path $StatePath
        if ($null -ne $previousState) {
            $restartCount = [int]$previousState.restart_count
            $remoteUrl = [string]$previousState.public_url
            $stdoutPath = [string]$previousState.stdout_path
            $stderrPath = [string]$previousState.stderr_path
            if ($null -ne $previousState.tunnel_pid -and (Test-TunnelProcess -ProcessId ([int]$previousState.tunnel_pid))) {
                $tunnelProcess = Get-Process -Id ([int]$previousState.tunnel_pid)
                $healthGraceUntil = $tunnelProcess.StartTime.AddSeconds(90)
            }
        }
        if ($null -ne $tunnelProcess -and [string]$config.mode -eq "Quick" -and [string]::IsNullOrWhiteSpace($remoteUrl)) {
            $remoteUrl = Find-QuickTunnelUrl -Paths @($stdoutPath, $stderrPath)
            if (-not [string]::IsNullOrWhiteSpace($remoteUrl)) {
                Enable-RemoteAccess -TunnelMode ([string]$config.mode) -RemoteUrl $remoteUrl | Out-Null
            }
        }
        Write-ServiceLog "Remote watchdog started in $($config.mode) mode"
        $remoteFailures = 0

        while (-not (Test-Path -LiteralPath $StopPath)) {
            if (-not (Test-LocalUiHealth)) {
                & powershell.exe -NoProfile -ExecutionPolicy Bypass -File $UiServiceScript -Action Start | Out-Null
                Start-Sleep -Seconds 3
                continue
            }

            if ($null -eq $tunnelProcess -or $tunnelProcess.HasExited) {
                $restartCount += 1
                Write-ServiceLog "Starting cloudflared attempt $restartCount"
                $tunnel = Start-CloudflaredProcess -Config $config -Attempt $restartCount
                $tunnelProcess = $tunnel.Process
                $healthGraceUntil = (Get-Date).AddSeconds(90)
                $stdoutPath = $tunnel.StdoutPath
                $stderrPath = $tunnel.StderrPath
                $remoteUrl = Wait-ForTunnelUrl -Config $config -Tunnel $tunnel
                if ([string]::IsNullOrWhiteSpace($remoteUrl)) {
                    if (-not $tunnelProcess.HasExited) {
                        Stop-Process -Id $tunnelProcess.Id -Force -ErrorAction SilentlyContinue
                    }
                    Write-RemoteState -Status "failed" -ServicePid $PID -TunnelPid $null `
                        -RestartCount $restartCount -TunnelMode ([string]$config.mode) -RemoteUrl "" `
                        -Message "cloudflared did not provide a public URL" -StdoutPath $stdoutPath -StderrPath $stderrPath
                    Start-Sleep -Seconds ([Math]::Min(30, 5 * $restartCount))
                    continue
                }
                Enable-RemoteAccess -TunnelMode ([string]$config.mode) -RemoteUrl $remoteUrl | Out-Null
                Write-ServiceLog "Tunnel URL ready: $remoteUrl"
                $remoteFailures = 0
            }

            if ([string]::IsNullOrWhiteSpace($remoteUrl)) {
                $remoteUrl = Find-QuickTunnelUrl -Paths @($stdoutPath, $stderrPath)
                if (-not [string]::IsNullOrWhiteSpace($remoteUrl)) {
                    Enable-RemoteAccess -TunnelMode ([string]$config.mode) -RemoteUrl $remoteUrl | Out-Null
                    Write-ServiceLog "Recovered tunnel URL from cloudflared logs: $remoteUrl"
                }
            }

            $remoteHealthy = Test-RemoteHealth -RemoteUrl $remoteUrl
            if ($remoteHealthy) {
                $remoteFailures = 0
                Write-RemoteState -Status "running" -ServicePid $PID -TunnelPid $tunnelProcess.Id `
                    -RestartCount $restartCount -TunnelMode ([string]$config.mode) -RemoteUrl $remoteUrl `
                    -Message "Authenticated HTTPS tunnel healthy" -StdoutPath $stdoutPath -StderrPath $stderrPath
            }
            else {
                if ((Get-Date) -lt $healthGraceUntil) {
                    Write-RemoteState -Status "starting" -ServicePid $PID -TunnelPid $tunnelProcess.Id `
                        -RestartCount $restartCount -TunnelMode ([string]$config.mode) -RemoteUrl $remoteUrl `
                        -Message "Waiting for public URL propagation" -StdoutPath $stdoutPath -StderrPath $stderrPath
                    Start-Sleep -Seconds 10
                    continue
                }
                $remoteFailures += 1
                Write-RemoteState -Status "degraded" -ServicePid $PID -TunnelPid $tunnelProcess.Id `
                    -RestartCount $restartCount -TunnelMode ([string]$config.mode) -RemoteUrl $remoteUrl `
                    -Message "Remote health check failed ($remoteFailures/3)" -StdoutPath $stdoutPath -StderrPath $stderrPath
                if ($remoteFailures -ge 3) {
                    Stop-Process -Id $tunnelProcess.Id -Force -ErrorAction SilentlyContinue
                    $tunnelProcess = $null
                    $remoteUrl = ""
                    $remoteFailures = 0
                }
            }
            Start-Sleep -Seconds 10
        }
        return 0
    }
    catch {
        Write-ServiceLog "Remote watchdog failed: $($_.Exception.Message)"
        Write-RemoteState -Status "failed" -ServicePid $PID `
            -TunnelPid $(if ($null -ne $tunnelProcess) { $tunnelProcess.Id } else { $null }) `
            -RestartCount $restartCount -TunnelMode $(if ($null -ne $config) { [string]$config.mode } else { "" }) `
            -RemoteUrl $remoteUrl -Message $_.Exception.Message -StdoutPath $stdoutPath -StderrPath $stderrPath
        return 1
    }
    finally {
        $mutex.ReleaseMutex()
        $mutex.Dispose()
    }
}

function Get-RemoteStatus {
    $task = Get-ScheduledTask -TaskName $TaskName -ErrorAction SilentlyContinue
    $state = Read-JsonFile -Path $StatePath
    $remoteUrl = if ($null -ne $state) { [string]$state.public_url } else { "" }
    return [ordered]@{
        local_healthy = Test-LocalUiHealth
        remote_healthy = if ([string]::IsNullOrWhiteSpace($remoteUrl)) { $false } else { Test-RemoteHealth -RemoteUrl $remoteUrl }
        public_url = $remoteUrl
        mode = if ($null -ne $state) { [string]$state.mode } else { "" }
        auto_start_installed = $null -ne $task
        task_state = if ($null -ne $task) { [string]$task.State } else { "NotInstalled" }
        service_pid = if ($null -ne $state) { $state.service_pid } else { $null }
        tunnel_pid = if ($null -ne $state) { $state.tunnel_pid } else { $null }
        restart_count = if ($null -ne $state) { $state.restart_count } else { 0 }
        message = if ($null -ne $state) { [string]$state.message } else { "" }
        share_file = $ShareTextPath
        service_log = $ServiceLog
    }
}

switch ($Action) {
    "Install" {
        Initialize-RuntimeRoot
        Write-RemoteConfig -TunnelMode $Mode -RemoteUrl $PublicUrl
        Enable-RemoteAccess -TunnelMode $Mode -RemoteUrl $PublicUrl -ForceNewPassword:$ResetPassword | Out-Null
        Install-RemoteTask
        Start-RemoteWatchdog
        $state = Wait-ForRemote
        if ($null -eq $state) {
            throw "The HTTPS tunnel did not become healthy. Check $ServiceLog"
        }
        if ($OpenBrowser) {
            Start-Process ([string]$state.public_url)
        }
        Write-Output "Online UI installed: $($state.public_url)"
        Write-Output "Share details: $ShareTextPath"
        exit 0
    }
    "Uninstall" {
        Stop-RemoteProcesses
        Unregister-ScheduledTask -TaskName $TaskName -Confirm:$false -ErrorAction SilentlyContinue
        Remove-Item -LiteralPath $ConfigPath, $AccessPath, $ShareJsonPath, $ShareTextPath -Force -ErrorAction SilentlyContinue
        Write-ServiceLog "Remote scheduled task and secrets removed"
        Write-Output "Online access removed. Local UI, project files, and results remain available."
        exit 0
    }
    "Start" {
        $config = Read-JsonFile -Path $ConfigPath
        if ($null -eq $config) {
            throw "Online access is not configured. Run install_online_ui.bat first."
        }
        Enable-RemoteAccess -TunnelMode ([string]$config.mode) -RemoteUrl ([string]$config.public_url) `
            -ForceNewPassword:$ResetPassword | Out-Null
        Start-RemoteWatchdog
        $state = Wait-ForRemote
        if ($null -eq $state) {
            throw "The HTTPS tunnel did not become healthy. Check $ServiceLog"
        }
        if ($OpenBrowser) {
            Start-Process ([string]$state.public_url)
        }
        Write-Output "Online UI ready: $($state.public_url)"
        Write-Output "Share details: $ShareTextPath"
        exit 0
    }
    "Stop" {
        Stop-RemoteProcesses
        Write-Output "Online access stopped. The local UI remains installed."
        exit 0
    }
    "Restart" {
        Stop-RemoteProcesses
        Start-Sleep -Seconds 2
        $config = Read-JsonFile -Path $ConfigPath
        Enable-RemoteAccess -TunnelMode ([string]$config.mode) -RemoteUrl ([string]$config.public_url) `
            -ForceNewPassword:$ResetPassword | Out-Null
        Start-RemoteWatchdog
        $state = Wait-ForRemote
        if ($null -eq $state) {
            throw "The HTTPS tunnel did not recover. Check $ServiceLog"
        }
        if ($OpenBrowser) {
            Start-Process ([string]$state.public_url)
        }
        Write-Output "Online UI restarted: $($state.public_url)"
        exit 0
    }
    "Run" {
        exit (Run-RemoteWatchdog)
    }
    "Status" {
        $status = Get-RemoteStatus
        if ($AsJson) {
            $status | ConvertTo-Json -Depth 4
        }
        else {
            Write-Output "AquaEsterNet-StructMine online UI"
            Write-Output ("Local UI healthy: {0}" -f $status.local_healthy)
            Write-Output ("Remote tunnel healthy: {0}" -f $status.remote_healthy)
            Write-Output ("Online URL: {0}" -f $status.public_url)
            Write-Output ("Mode: {0}" -f $status.mode)
            Write-Output ("Auto-start installed: {0}" -f $status.auto_start_installed)
            Write-Output ("Task state: {0}" -f $status.task_state)
            Write-Output ("Tunnel PID: {0}" -f $status.tunnel_pid)
            Write-Output ("Share file: {0}" -f $status.share_file)
            Write-Output ("Log: {0}" -f $status.service_log)
        }
        exit 0
    }
}
