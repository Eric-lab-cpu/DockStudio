@echo off
REM DockStudio PyInstaller build (Windows, advanced)
REM Preconditions: conda env dockstudio active, pip install pyinstaller
cd /d "%~dp0..\.."
python -m PyInstaller --noconfirm packaging\pyinstaller\dockstudio.spec
echo Build finished: dist\DockStudio\
