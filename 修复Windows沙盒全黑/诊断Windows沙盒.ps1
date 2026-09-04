<#
    Diagnose-WindowsSandbox.ps1  (Eric Studio / DockStudio tools)
    Collects the minimum info needed to debug a black-screen / won't-start
    Windows Sandbox. Run:  powershell -ExecutionPolicy Bypass -File this.ps1
    Output: a .txt next to this script.
#>
$ErrorActionPreference = 'Continue'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$out  = Join-Path $here ('WindowsSandbox_Diagnose_' + (Get-Date -Format 'yyyyMMdd_HHmmss') + '.txt')
$L = New-Object System.Collections.Generic.List[string]

function Log($m){ $L.Add($m) }

Log '===== Windows Sandbox black-screen diagnostic ====='
Log ('Time: ' + (Get-Date))
Log ('User: ' + [Environment]::UserName)
Log ('Elevated: ' + ([Security.Principal.WindowsPrincipal][Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole([Security.Principal.WindowsBuiltInRole]::Administrator))

Log ''
Log '--- OS ---'
try { Get-ComputerInfo -Property OsName,OsVersion,OsBuildNumber,WindowsVersion | Format-List | Out-String -Width 200 | ForEach-Object { Log $_ } } catch { Log ('Get-ComputerInfo failed: ' + $_.Exception.Message) }

Log ''
Log '--- CPU virtualization (VM monitor mode extensions) ---'
try {
    $ci = Get-CimInstance Win32_Processor
    Log ('CPU: ' + $ci.Name)
    Log ('VirtualizationFirmwareEnabled: ' + $ci.VirtualizationFirmwareEnabled)
    Log ('SecondLevelAddressTranslationExtensions: ' + $ci.SecondLevelAddressTranslationExtensions)
} catch { Log ('CPU query failed: ' + $_.Exception.Message) }

Log ''
Log '--- Memory (free GB) ---'
try {
    $os = Get-CimInstance Win32_OperatingSystem
    Log ('TotalVisible: ' + [math]::Round($os.TotalVisibleMemorySize/1MB,1) + ' GB; Free: ' + [math]::Round($os.FreePhysicalMemory/1MB,1) + ' GB')
} catch {}

Log ''
Log '--- Optional features state ---'
try {
    Get-WindowsOptionalFeature -Online | Where-Object { $_.FeatureName -match 'Containers-DisposableClientVM|VirtualMachinePlatform|Microsoft-Hyper-V' } |
        Select-Object FeatureName,State | Format-Table -AutoSize | Out-String -Width 200 | ForEach-Object { Log $_ }
} catch { Log ('Feature query failed (need admin): ' + $_.Exception.Message) }

Log ''
Log '--- Sandbox related services ---'
try {
    Get-Service -ErrorAction SilentlyContinue | Where-Object { $_.Name -match 'vmcompute|vmms|HvHost|rdpinput' } |
        Select-Object Name,Status,StartType | Format-Table -AutoSize | Out-String -Width 200 | ForEach-Object { Log $_ }
} catch {}

Log ''
Log '--- GPU driver (host) ---'
try {
    Get-CimInstance Win32_VideoController | Select-Object Name,DriverVersion,DriverDate,Status | Format-Table -AutoSize | Out-String -Width 200 | ForEach-Object { Log $_ }
} catch { Log ('GPU query failed: ' + $_.Exception.Message) }

Log ''
Log '--- Recent errors (last 30 min, Hyper-V / Sandbox) ---'
try {
    $since = (Get-Date).AddMinutes(-30)
    Get-WinEvent -FilterHashtable @{ LogName='Microsoft-Windows-Hyper-V-Worker/Admin'; StartTime=$since } -MaxEvents 20 -ErrorAction SilentlyContinue |
        Select-Object TimeCreated,Id,LevelDisplayName,Message | Format-List | Out-String -Width 200 | ForEach-Object { Log $_ }
} catch { Log ('No Hyper-V-Worker events / need admin: ' + $_.Exception.Message) }
try {
    Get-WinEvent -FilterHashtable @{ LogName='Application'; StartTime=$since } -MaxEvents 40 -ErrorAction SilentlyContinue |
        Where-Object { $_.Message -match 'Sandbox|Hyper-V|vmcompute|VID' } |
        Select-Object TimeCreated,Id,ProviderName,Message | Format-List | Out-String -Width 200 | ForEach-Object { Log $_ }
} catch {}

$L | Set-Content -Path $out -Encoding UTF8
Write-Host ''
Write-Host ('Diagnostic saved to: ' + $out) -ForegroundColor Green
Write-Host 'Tip: if you cannot fix it yourself, attach this file when asking for support.'
