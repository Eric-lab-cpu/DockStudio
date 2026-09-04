@echo off
rem ============================================================
rem  DockStudio DEBUG launcher (Eric Studio)
rem  Same as DockStudio.bat but keeps a console window and shows
rem  the Python traceback if the GUI fails to start.
rem ============================================================
setlocal EnableExtensions
set "ROOT=%~dp0"
set "ENV=%ROOT%env"

if not exist "%ENV%\python.exe" (
    call "%ROOT%scripts\setup_env.bat"
    if errorlevel 1 (
        pause
        exit /b 1
    )
)

set "DOCKSTUDIO_VINA=%ROOT%vina\vina.exe"
set "PYTHONPATH=%ROOT%app"
set "PYTHONNOUSERSITE=1"
set "PATH=%ENV%\Scripts;%ENV%\Library\bin;%PATH%"
cd /d "%ROOT%app"

echo Starting DockStudio in debug mode...
"%ENV%\python.exe" -X faulthandler -m dockstudio
echo.
echo DockStudio exited with code %errorlevel%
pause
endlocal
