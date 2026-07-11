[CmdletBinding()]
param(
    [string]$DistroName = "Ubuntu-Foldseek",
    [string]$InstallLocation = "D:\Apps\WSL\Ubuntu-Foldseek",
    [string]$UbuntuRootfs = "D:\Apps\WSL\downloads\ubuntu-jammy-wsl-amd64-ubuntu22.04lts.rootfs.tar.gz",
    [string]$FoldseekArchive = "D:\Apps\Foldseek\downloads\foldseek-linux-avx2.tar.gz"
)

$ErrorActionPreference = "Stop"
$ProgressPreference = "SilentlyContinue"
$ExpectedUbuntuSha256 = "1483cc5c1dce13064f774834cbffdff226559fd522a67a381a8ea77d63fb4109"
$ExpectedFoldseekSha256 = "3bafbc46cdfe6ac0d196edde5da40da99413b8410f12f0e2baee72d2edecc124"
$UbuntuSourceUrl = "https://cloud-images.ubuntu.com/wsl/jammy/current/ubuntu-jammy-wsl-amd64-ubuntu22.04lts.rootfs.tar.gz"
$FoldseekSourceUrl = "https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz"
$ProjectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path

function Assert-Administrator {
    $identity = [Security.Principal.WindowsIdentity]::GetCurrent()
    $principal = [Security.Principal.WindowsPrincipal]::new($identity)
    if (-not $principal.IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator)) {
        throw "Run this script from an elevated PowerShell window."
    }
}

function Assert-Sha256([string]$Path, [string]$Expected) {
    if (-not (Test-Path -LiteralPath $Path -PathType Leaf)) {
        throw "Required download is missing: $Path"
    }
    $actual = (Get-FileHash -LiteralPath $Path -Algorithm SHA256).Hash.ToLowerInvariant()
    if ($actual -ne $Expected.ToLowerInvariant()) {
        throw "SHA-256 mismatch for $Path. Expected $Expected, observed $actual"
    }
}

function Invoke-Wsl([string[]]$Arguments) {
    & wsl.exe @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "wsl.exe failed with exit code ${LASTEXITCODE}: $($Arguments -join ' ')"
    }
}

function Convert-ToWslMountPath([string]$Path) {
    $resolved = (Resolve-Path -LiteralPath $Path).Path
    if ($resolved -notmatch '^(?<drive>[A-Za-z]):\\(?<rest>.*)$') {
        throw "Expected an absolute Windows drive path, observed: $resolved"
    }
    $drive = $Matches.drive.ToLowerInvariant()
    $rest = $Matches.rest -replace '\\', '/'
    return "/mnt/$drive/$rest"
}

Assert-Administrator

$requiredFeatures = @(
    "Microsoft-Windows-Subsystem-Linux",
    "VirtualMachinePlatform"
)
foreach ($featureName in $requiredFeatures) {
    $feature = Get-WindowsOptionalFeature -Online -FeatureName $featureName
    if ($feature.State -ne "Enabled") {
        throw "Windows feature is not enabled: $featureName"
    }
}
if (Test-Path "HKLM:\SOFTWARE\Microsoft\Windows\CurrentVersion\Component Based Servicing\RebootPending") {
    throw "Windows still has a pending component restart. Restart Windows, then rerun this script."
}

Assert-Sha256 $UbuntuRootfs $ExpectedUbuntuSha256
Assert-Sha256 $FoldseekArchive $ExpectedFoldseekSha256

Invoke-Wsl @("--update", "--web-download")
Invoke-Wsl @("--set-default-version", "2")

$installed = @(& wsl.exe --list --quiet) | ForEach-Object { ($_ -replace "`0", "").Trim() } | Where-Object { $_ }
if ($DistroName -notin $installed) {
    if (Test-Path -LiteralPath $InstallLocation) {
        $existing = @(Get-ChildItem -LiteralPath $InstallLocation -Force -ErrorAction SilentlyContinue)
        if ($existing.Count -gt 0) {
            throw "Install location is not empty, but the distro is not registered: $InstallLocation"
        }
    }
    New-Item -ItemType Directory -Path $InstallLocation -Force | Out-Null
    Invoke-Wsl @("--import", $DistroName, $InstallLocation, $UbuntuRootfs, "--version", "2")
}

$linuxArchive = Convert-ToWslMountPath $FoldseekArchive
if (-not $linuxArchive) {
    throw "Could not convert the Foldseek archive path for WSL."
}

$linuxSetup = @'
set -euo pipefail
export DEBIAN_FRONTEND=noninteractive
apt-get update
apt-get install -y --no-install-recommends ca-certificates curl gzip tar sudo
printf '%s  %s\n' '__EXPECTED_SHA256__' '__ARCHIVE_PATH__' | sha256sum --check --status -
mkdir -p /opt/foldseek
tar -xzf '__ARCHIVE_PATH__' -C /opt/foldseek --strip-components=1
chmod 0755 /opt/foldseek/bin/foldseek
ln -sfn /opt/foldseek/bin/foldseek /usr/local/bin/foldseek
if ! id -u aqua >/dev/null 2>&1; then
  useradd --create-home --shell /bin/bash aqua
fi
usermod -aG sudo aqua
printf 'aqua ALL=(ALL) NOPASSWD:ALL\n' >/etc/sudoers.d/90-aqua
chmod 0440 /etc/sudoers.d/90-aqua
cat >/etc/wsl.conf <<'EOF'
[boot]
systemd=true

[user]
default=aqua

[automount]
enabled=true
options=metadata,umask=022,fmask=011
EOF
grep -qw avx2 /proc/cpuinfo
/opt/foldseek/bin/foldseek version
'@
$linuxSetupCommand = $linuxSetup.Replace("__ARCHIVE_PATH__", $linuxArchive).Replace("__EXPECTED_SHA256__", $ExpectedFoldseekSha256)

Invoke-Wsl @(
    "-d", $DistroName,
    "-u", "root",
    "--", "bash", "-lc", $linuxSetupCommand
)
Invoke-Wsl @("--terminate", $DistroName)

$foldseekVersion = (((& wsl.exe -d $DistroName -- /opt/foldseek/bin/foldseek version) -replace "`0", "") -join "`n").Trim()
$distroInfo = ((& wsl.exe --list --verbose) -replace "`0", "") -join "`n"
$wslVersionOutput = try { (((& wsl.exe --version) -replace "`0", "") -join "`n").Trim() } catch { "" }
$wslVersionValues = @([regex]::Matches($wslVersionOutput, '\d+(?:\.\d+)+(?:-\d+)?') | ForEach-Object { $_.Value })
$componentNames = @("wsl", "kernel", "wslg", "msrdc", "direct3d", "dxcore", "windows")
$wslComponents = [ordered]@{}
for ($index = 0; $index -lt [Math]::Min($componentNames.Count, $wslVersionValues.Count); $index++) {
    $wslComponents[$componentNames[$index]] = $wslVersionValues[$index]
}
$ubuntuOsRelease = (((& wsl.exe -d $DistroName -- cat /etc/os-release) -replace "`0", "") -join "`n").Trim()

$manifest = [ordered]@{
    created_at = (Get-Date).ToUniversalTime().ToString("o")
    distro_name = $DistroName
    distro_install_location = $InstallLocation
    distro_info = $distroInfo
    wsl_version = if ($wslComponents.wsl) { $wslComponents.wsl } else { "unknown" }
    wsl_kernel_version = if ($wslComponents.kernel) { $wslComponents.kernel } else { "unknown" }
    wsl_component_versions = $wslComponents
    ubuntu_os_release = $ubuntuOsRelease
    ubuntu = [ordered]@{
        source_url = $UbuntuSourceUrl
        local_archive = $UbuntuRootfs
        sha256 = $ExpectedUbuntuSha256
    }
    foldseek = [ordered]@{
        source_url = $FoldseekSourceUrl
        local_archive = $FoldseekArchive
        sha256 = $ExpectedFoldseekSha256
        version = $foldseekVersion
        executable = "/opt/foldseek/bin/foldseek"
        build = "linux-avx2"
    }
    database_root = (Join-Path $ProjectRoot "databases")
}
$manifestDir = Join-Path $ProjectRoot "databases\environment"
New-Item -ItemType Directory -Path $manifestDir -Force | Out-Null
$manifestPath = Join-Path $manifestDir "foldseek_install_manifest.json"
$manifest | ConvertTo-Json -Depth 6 | Set-Content -LiteralPath $manifestPath -Encoding UTF8

Write-Host "[OK] WSL distribution: $DistroName"
Write-Host "[OK] Install location: $InstallLocation"
Write-Host "[OK] Foldseek version: $foldseekVersion"
Write-Host "[OK] Manifest: $manifestPath"
