@echo off
rem ============================================================
rem  DockStudio launcher (Eric Studio)
rem  First run: creates the engine environment automatically.
rem  Runs the GUI in a console so startup errors are visible and
rem  the window does not silently flash away.
rem ============================================================
setlocal EnableExtensions
set "ROOT=%~dp0"
set "ENV=%ROOT%env"

if not exist "%ENV%\python.exe" (
    echo.
    echo   First run: preparing DockStudio engine environment.
    echo   This downloads PyMOL / RDKit / Meeko / PLIP etc. and may take
    echo   several minutes depending on your network. Please wait...
    echo.
    call "%ROOT%scripts\setup_env.bat"
)
if not exist "%ENV%\python.exe" (
    echo.
    echo   Environment setup FAILED. See messages above.
    pause
    exit /b 1
)

set "DOCKSTUDIO_VINA=%ROOT%vina\vina.exe"
set "PYTHONPATH=%ROOT%app"
set "PYTHONNOUSERSITE=1"
set "PATH=%ENV%\Scripts;%ENV%\Library\bin;%PATH%"
cd /d "%ROOT%app"

echo Starting DockStudio ...
"%ENV%\python.exe" -m dockstudio
if errorlevel 1 (
    echo.
    echo   DockStudio exited with an error. The message above explains why.
    pause
)
endlocal
