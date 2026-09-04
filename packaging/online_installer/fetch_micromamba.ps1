param(
    [Parameter(Mandatory=$true)][string]$Dest,
    [Parameter(Mandatory=$true)][string]$Expected
)

# DockStudio / Eric Studio - self-healing micromamba download with sha256 verify.
# Called when the bundled micromamba.exe is missing (e.g. removed by antivirus).
$ErrorActionPreference = 'Stop'
$exp = $Expected.Trim().ToLowerInvariant()
$ver = '2.9.0-0'
$base = "https://github.com/mamba-org/micromamba-releases/releases/download/$ver/micromamba-win-64.exe"
$urls = @(
    $base,
    "https://gh-proxy.com/$base",
    "https://ghfast.top/$base",
    "https://gh-proxy.com/https://github.com/mamba-org/micromamba-releases/releases/download/$ver/micromamba-win-64.exe"
)
$dir = Split-Path -Parent $Dest
if (-not (Test-Path $dir)) { New-Item -ItemType Directory -Path $dir -Force | Out-Null }

$lastErr = ''
foreach ($u in $urls) {
    try {
        Write-Host "Downloading micromamba: $u"
        [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
        Invoke-WebRequest -Uri $u -OutFile $Dest -UseBasicParsing -TimeoutSec 240
        if (-not (Test-Path $Dest)) { continue }
        $h = (Get-FileHash -Path $Dest -Algorithm SHA256).Hash.ToLowerInvariant()
        if ($h -eq $exp) {
            Write-Host "OK verified sha256=$h"
            exit 0
        } else {
            Write-Host "HASH MISMATCH expected=$exp got=$h"
            Remove-Item $Dest -Force -ErrorAction SilentlyContinue
        }
    } catch {
        $lastErr = $_.Exception.Message
        Write-Host "Download failed: $lastErr"
    }
}
Write-Host "ALL_DOWNLOADS_FAILED last=$lastErr"
exit 5
