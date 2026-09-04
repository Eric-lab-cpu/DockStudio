@echo off
rem ============================================================
rem  DockStudio - optional PyMOL installer (Eric Studio)
rem  Use this if the automatic PyMOL step failed. It uses the
rem  official conda-forge channel + a short package-cache path to
rem  avoid Windows long-path issues with qt-main, and retries.
rem ============================================================
setlocal EnableExtensions EnableDelayedExpansion
set "ROOT=%~dp0"
set "SCRIPTS=%ROOT%scripts"
set "ENV=%ROOT%env"
set "LOG=%ROOT%DockStudio_setup.log"
set "MM=%SCRIPTS%micromamba.exe"

if not exist "%MM%" (
    echo ERROR: micromamba.exe not found in %SCRIPTS%
    pause
    exit /b 1
)
if not exist "%ENV%\python.exe" (
    echo ERROR: DockStudio environment not created yet. Run DockStudio.bat first.
    pause
    exit /b 1
)

set "MAMBA_ROOT_PREFIX=%ROOT%\.mamba"
set "CONDA_PKGS_DIRS=%USERPROFILE%\.ds_py_pkgs"
set /a PMTRY=0

:PYMOL_TRY
set /a PMTRY+=1
echo.
echo Installing PyMOL - attempt !PMTRY! of 3 ...
"%MM%" install -y -p "%ENV%" -c conda-forge pymol-open-source
if errorlevel 1 goto :PYMOL_RETRY
echo.
echo PyMOL installed successfully.
echo pymol-ok >> "%LOG%"
goto :END

:PYMOL_RETRY
if %PMTRY% GEQ 3 goto :PYMOL_FAIL
echo   Retrying with a clean cache ...
if exist "%CONDA_PKGS_DIRS%" rmdir /s /q "%CONDA_PKGS_DIRS%"
goto :PYMOL_TRY

:PYMOL_FAIL
echo.
echo PyMOL install failed after 3 attempts.
echo pymol-optional-failed >> "%LOG%"
echo DockStudio still works - only 3D ray-traced PNGs are skipped.

:END
pause
endlocal
