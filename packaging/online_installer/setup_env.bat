@echo off
rem ============================================================
rem  DockStudio environment bootstrap (Eric Studio, online)
rem
rem  Strategy (robust):
rem   1. micromamba creates ONLY a tiny python=3.10 + tk env
rem      (minimal conda surface -> fewer chances of extraction issues)
rem   2. Everything else is pip-installed from Windows wheels
rem      (plain unzip into site-packages, far more reliable than
rem      conda .conda extraction + linking on locked-down systems)
rem   3. PyMOL is OPTIONAL (conda-only on Windows); if it fails,
rem      DockStudio still runs - only 3D ray-traced PNGs are skipped.
rem ============================================================
setlocal EnableExtensions EnableDelayedExpansion

rem --- canonical paths ----------------------------------------
set "SCRIPTS=%~dp0"
for %%I in ("%SCRIPTS%..") do set "APP_ROOT=%%~fI"
set "ENV=%APP_ROOT%\env"
set "LOG=%APP_ROOT%\DockStudio_setup.log"
rem official micromamba 2.9.0-0 win-64 sha256
set "MM_SHA256=a6d804394b2418991c4e29562853eaace2f2ce9d9da661a98e74e02e8dbb44b0"

rem China-friendly mirrors (Tsinghua first)
set "CF_URL=https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge"
set "PIP_URL=https://pypi.tuna.tsinghua.edu.cn/simple"

rem Detect Windows Sandbox (WDAGUtilityAccount is its default user)
set "IS_SANDBOX=0"
if /i "%USERNAME%"=="WDAGUtilityAccount" set "IS_SANDBOX=1"

echo === DockStudio setup start: %date% %time% >> "%LOG%"

if "%IS_SANDBOX%"=="1" (
    echo.
    echo   *** Detected: Windows Sandbox ***
    echo   The Windows Sandbox uses a locked-down filesystem and antivirus
    echo   that frequently break large package extractions. The install may
    echo   fail here even though it works on a normal Windows system.
    echo   If the steps below fail, please install DockStudio on your real
    echo   Windows desktop instead of inside the Sandbox.
    echo.
)

rem --- locate bundled micromamba.exe --------------------------
set "MM="
if exist "%SCRIPTS%micromamba.exe"            set "MM=%SCRIPTS%micromamba.exe"
if not defined MM if exist "%APP_ROOT%\scripts\micromamba.exe" set "MM=%APP_ROOT%\scripts\micromamba.exe"
if not defined MM if exist "%APP_ROOT%\micromamba.exe"         set "MM=%APP_ROOT%\micromamba.exe"
if defined MM goto :HAVE_MM

echo.
echo   micromamba.exe was not found. Trying to re-download and verify...
powershell -NoProfile -ExecutionPolicy Bypass -File "%SCRIPTS%fetch_micromamba.ps1" ^
    -Dest "%APP_ROOT%\scripts\micromamba.exe" -Expected "%MM_SHA256%"
if exist "%APP_ROOT%\scripts\micromamba.exe" set "MM=%APP_ROOT%\scripts\micromamba.exe"
if defined MM goto :HAVE_MM

echo.
echo   ERROR: could not obtain micromamba.exe.
echo   micromamba-unavailable >> "%LOG%"
pause
exit /b 2

:HAVE_MM
echo Using micromamba: %MM% >> "%LOG%"

rem keep micromamba private to this installation
set "MAMBA_ROOT_PREFIX=%APP_ROOT%\.mamba"
set "CONDA_PKGS_DIRS=%APP_ROOT%\.pkgs"
set "PYTHONPATH=%APP_ROOT%\app"
set "PIP_INDEX_URL=%PIP_URL%"

if exist "%ENV%\python.exe" goto :READY

rem ============================================================
rem  STEP 1 - tiny base env: python + tk only (conda)
rem ============================================================
set /a STEP=0

:STEP1
set /a STEP+=1
echo.
echo [1/3] Creating base Python 3.10 env - attempt !STEP! of 2 ...
"%MM%" create -y -p "%ENV%" -c "%CF_URL%" -c conda-forge python=3.10 tk=8.6
if errorlevel 1 goto :STEP1_RETRY
goto :STEP1_OK

:STEP1_RETRY
if %STEP% GEQ 2 goto :BASE_FAILED
echo   Cleaning and retrying once...
if exist "%APP_ROOT%\.pkgs" rmdir /s /q "%APP_ROOT%\.pkgs"
if exist "%ENV%" rmdir /s /q "%ENV%"
if exist "%ENV%" (echo   Warning: could not fully remove old env & rmdir /s /q "%ENV%" 2>nul)
goto :STEP1

:STEP1_OK
echo Base environment created.
echo.

rem ============================================================
rem  STEP 2 - pip install scientific stack from Windows wheels
rem ============================================================
echo [2/3] Installing scientific stack via pip (Windows wheels) ...
rem NOTE: prody has no Windows pip wheel and is NOT required by the DockStudio
rem engine - it is intentionally omitted here.
"%ENV%\python.exe" -m pip install --no-input --disable-pip-version-check ^
    numpy pandas scipy pillow rdkit gemmi meeko plip ttkbootstrap
if errorlevel 1 goto :PIP_FAILED
echo Scientific stack installed.
echo.

rem ============================================================
rem  STEP 3 - verify runtime
rem ============================================================
echo [3/4] Verifying runtime - python + tkinter + DockStudio import ...
"%ENV%\python.exe" -c "import sys, tkinter, dockstudio, rdkit; print('RUNTIME_OK', sys.version.split()[0])"
if errorlevel 1 goto :VERIFY_FAILED

rem ============================================================
rem  OPTIONAL - PyMOL (conda-only on Windows; non-fatal)
rem  Uses the OFFICIAL conda-forge channel and a SHORT package-cache
rem  path to avoid Windows MAX_PATH problems with the deep qt-main
rem  header paths; retries a few times.
rem ============================================================
set "OLD_PKGS=%CONDA_PKGS_DIRS%"
set "CONDA_PKGS_DIRS=%USERPROFILE%\.ds_py_pkgs"
set /a PMTRY=0

:PYMOL_TRY
set /a PMTRY+=1
echo.
echo [optional] Installing PyMOL for 3D rendering - attempt !PMTRY! of 3 ...
echo            If this fails, DockStudio still works - 3D ray-traced
echo            PNGs are skipped; 2D images and .pml scripts still work.
"%MM%" install -y -p "%ENV%" -c conda-forge pymol-open-source
if errorlevel 1 goto :PYMOL_RETRY
echo PyMOL installed.
echo pymol-ok >> "%LOG%"
set "CONDA_PKGS_DIRS=%OLD_PKGS%"
goto :PYMOL_DONE

:PYMOL_RETRY
if %PMTRY% GEQ 3 goto :PYMOL_FAIL
echo   Retrying PyMOL with a clean cache ...
if exist "%CONDA_PKGS_DIRS%" rmdir /s /q "%CONDA_PKGS_DIRS%"
goto :PYMOL_TRY

:PYMOL_FAIL
echo.
echo Warning: PyMOL install failed after retries - continuing without 3D rendering.
echo pymol-optional-failed >> "%LOG%"
set "CONDA_PKGS_DIRS=%OLD_PKGS%"

:PYMOL_DONE
echo.
echo DockStudio environment ready.
echo done >> "%LOG%"
exit /b 0

:READY
echo Environment already present at %ENV%
exit /b 0

:BASE_FAILED
echo.
echo   ERROR: could not create even the minimal Python environment.
echo   base-env-failed >> "%LOG%"
echo.
echo   This almost always means the Windows Sandbox / antivirus is blocking
echo   package extraction. Recommended fixes:
echo     1. Install DockStudio on your REAL Windows system, not in the Sandbox.
echo     2. If on a real system: add an antivirus exclusion for this folder:
echo        %APP_ROOT%
echo     3. Re-run DockStudio.bat - it retries automatically.
echo.
echo   Free space:
powershell -NoProfile -Command "$d=Get-PSDrive -Name C; Write-Host ('Free on C: ' + [math]::Round($d.Free/1GB,1) + ' GB')" 2>nul
pause
exit /b 1

:PIP_FAILED
echo.
echo   ERROR: pip install failed.
echo   pip-failed >> "%LOG%"
echo   Tip: if you are in mainland China the Tsinghua PyPI mirror is used;
echo   if it is unreachable, temporarily disable PIP_INDEX_URL in
echo   setup_env.bat.
pause
exit /b 1

:VERIFY_FAILED
echo.
echo   ERROR: runtime verification failed - tkinter or DockStudio import.
echo   verify-failed >> "%LOG%"
pause
exit /b 1
