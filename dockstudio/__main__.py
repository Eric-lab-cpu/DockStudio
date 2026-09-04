"""``python -m dockstudio`` launches the graphical application."""

import sys


def main() -> int:
    from dockstudio.gui.app import main as gui_main

    return gui_main(sys.argv)


if __name__ == "__main__":
    raise SystemExit(main())
