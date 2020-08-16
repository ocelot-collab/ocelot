"""Ocelot GUI tool"""

import sys
from PyQt5.QtWidgets import QApplication

from app.gui import GUIWindow


class OcelotWindow():
    
    def __init__(self):
        self.gui = GUIWindow()


def main():

    app = QApplication(sys.argv)
    window = OcelotWindow()
    window.gui.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
