@echo off
REM ============================================================
REM DockStudio (Eric Studio) - build Windows online installer
REM Works on Windows where NSIS is installed (makensis on PATH).
REM Usage: packaging\online_installer\build_online_installer.bat
REM Output: dist\DockStudio_Setup_1.0.0.exe
REM ============================================================
setlocal
cd /d "%~dp0..\.."

echo [1/2] staging payload...
python packaging\online_installer\prepare_payload.py
if errorlevel 1 exit /b 1

echo [2/2] compiling installer with makensis...
where makensis >nul 2>nul
if errorlevel 1 (
    echo ERROR: makensis not found. Install NSIS (choco install nsis) or add to PATH.
    exit /b 1
)
makensis packaging\online_installer\DockStudio_online.nsi
if errorlevel 1 exit /b 1
echo DONE: dist\DockStudio_Setup_1.0.0.exe
endlocal
