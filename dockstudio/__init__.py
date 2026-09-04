"""DockStudio: reproducible, batch, GUI-first molecular docking platform.

Engines are importable without tkinter so the pipeline can be tested and
scripted on headless machines; the graphical application lives in
``dockstudio.gui``.
"""

from ._version import (APP_NAME, APP_NAME_ASCII, APP_TITLE, BRAND,
                       COPYRIGHT, COPYRIGHT_CN, __version__)

__all__ = ["APP_NAME", "APP_NAME_ASCII", "APP_TITLE", "BRAND",
           "COPYRIGHT", "COPYRIGHT_CN", "__version__"]
