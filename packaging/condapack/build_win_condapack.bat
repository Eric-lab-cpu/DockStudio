@echo off
REM DockStudio conda-pack build (Windows) - produce a green portable folder
REM Usage:
REM   1) conda env create -f packaging/environment.yml
REM   2) conda activate dockstudio
REM   3) packaging\condapack\build_win_condapack.bat
REM Output: dist\DockStudioPortable  (copy anywhere, run DockStudio.bat)

setlocal enabledelayedexpansion
cd /d "%~dp0..\.."

set ENV_NAME=dockstudio
set OUT=dist\DockStudioPortable

echo [1/5] resolve conda env prefix...
for /f "delims=" %%i in ('conda env list') do (
    echo %%i | findstr /r /c:"^%ENV_NAME% " >nul
    if not errorlevel 1 for /f "tokens=2" %%p in ("%%i") do set ENV_PREFIX=%%p
)
if "%ENV_PREFIX%"=="" (
    echo ERROR: cannot find conda env "%ENV_NAME%". Run: conda env create -f packaging/environment.yml
    exit /b 1
)
echo env prefix: %ENV_PREFIX%

echo [2/5] conda-pack the environment...
if exist "%OUT%" rmdir /s /q "%OUT%"
mkdir "%OUT%"
conda pack -n %ENV_NAME% -o "%OUT%\dockstudio_env.tar.gz"
if errorlevel 1 exit /b 1

echo [3/5] extract...
tar -xzf "%OUT%\dockstudio_env.tar.gz" -C "%OUT%"
del "%OUT%\dockstudio_env.tar.gz"

echo [4/5] copy DockStudio python package...
xcopy dockstudio "%OUT%\dockstudio\" /E /I /Y >nul
xcopy examples\demo "%OUT%\examples\demo\" /E /I /Y >nul 2>nul

echo [5/5] write launcher...
(
echo @echo off
echo cd /d "%%~dp0"
echo set PYTHONNOUSERSITE=1
echo set PYTHONPATH=%%~dp0
echo start "" "%%~dp0pythonw.exe" -m dockstudio
) > "%OUT%\DockStudio.bat"

echo.
echo DONE: "%OUT%"
echo Copy this whole folder to any Windows PC and run DockStudio.bat
echo (no installation needed).
endlocal
