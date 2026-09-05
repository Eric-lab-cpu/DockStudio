@echo off
REM ============================================================
REM DockStudio (Eric Studio) - Build Windows Setup Installer
REM ============================================================
REM Requirements on the BUILD machine:
REM   1) Miniconda + env "dockstudio" (see packaging/environment.yml)
REM   2) Inno Setup 6  https://jrsoftware.org/isinfo.php  (ISCC on PATH)
REM
REM Produces: dist\installer\DockStudio_Setup_1.1.0.exe
REM ============================================================
setlocal enabledelayedexpansion
cd /d "%~dp0..\.."

echo [1/4] locate conda env / portable folder...
if not exist "dist\DockStudioPortable\pythonw.exe" (
    echo [info] portable folder missing - building it with conda-pack...
    call packaging\condapack\build_win_condapack.bat
    if errorlevel 1 exit /b 1
)

echo [2/4] locate Inno Setup compiler...
set ISCC=
where iscc >nul 2>nul && set ISCC=iscc
if "%ISCC%"=="" (
    if exist "%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe" set ISCC="%ProgramFiles(x86)%\Inno Setup 6\ISCC.exe"
)
if "%ISCC%"=="" if exist "%ProgramFiles%\Inno Setup 6\ISCC.exe" set ISCC="%ProgramFiles%\Inno Setup 6\ISCC.exe"
if "%ISCC%"=="" (
    echo ERROR: Inno Setup 6 not found. Install from https://jrsoftware.org/isinfo.php
    exit /b 1
)

echo [3/4] version...
set VER=1.1.0

echo [4/4] compile installer...
%ISCC% /DAPP_VERSION=%VER% /DPORT="dist\DockStudioPortable" packaging\installer\DockStudio_setup.iss
if errorlevel 1 exit /b 1

echo.
echo DONE: dist\installer\DockStudio_Setup_%VER%.exe
endlocal
