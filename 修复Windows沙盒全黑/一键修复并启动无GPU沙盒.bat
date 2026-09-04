@echo off
rem ============================================================
rem  Windows Sandbox black-screen repair helper (run as admin)
rem  Fixes feature state / corrupted image and offers the
rem  no-GPU sandbox config as a workaround.
rem  Usage: right-click -> Run as administrator
rem ============================================================
setlocal
net session >nul 2>&1 || (
    echo.
    echo   Please RIGHT-CLICK this file and choose
    echo   "Run as administrator".
    echo.
    pause
    exit /b 1
)

echo [1/4] Repair Windows image with DISM ...
DISM /Online /Cleanup-Image /RestoreHealth
echo.

echo [2/4] System File Checker (sfc) ...
sfc /scannow
echo.

echo [3/4] Make sure Sandbox / VirtualMachinePlatform / Hyper-V features exist ...
powershell -NoProfile -ExecutionPolicy Bypass -Command ^
  "Enable-WindowsOptionalFeature -Online -FeatureName Containers-DisposableClientVM -All -NoRestart; " ^
  "Enable-WindowsOptionalFeature -Online -FeatureName VirtualMachinePlatform -All -NoRestart; " ^
  "Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Hyper-V-All -All; " ^
  "Get-WindowsOptionalFeature -Online | Where-Object {$_.FeatureName -match 'Sandbox|VirtualMachinePlatform|Hyper-V'} | Select-Object FeatureName,State"

echo.
echo [4/4] Launch Sandbox with GPU disabled (black-screen workaround) ...
set "WSB=%~dp0Sandbox_NoGPU.wsb"
if exist "%WSB%" (
    start "" WindowsSandbox.exe "%WSB%"
) else (
    echo WSB config not found: %WSB%
)

echo.
echo Done. If the sandbox still fails, reboot and try again.
echo If you prefer hardware acceleration later, update the GPU driver and
echo remove the ^<vGPU^>disable^</vGPU^> line from Sandbox_NoGPU.wsb.
pause
endlocal
