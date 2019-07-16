"""GUI description"""

import os
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QDesktopWidget, QMainWindow, QAction
from app.ui_forms.main import *

from app.menu_file import *
from app.menu_edit import *
from app.menu_simulation import *
from app.lattice import *


class GUIWindow(QMainWindow, Ui_MainWindow):

    def __init__(self):

        super(GUIWindow, self).__init__()

        self.setupUi(self)
        self.initUI()
        #self.loadStyleSheet()

        self.lattice = GUILattice()


    def __del__(self):
        pass


    def initUI(self):

        # Init main window
        self.statusBar.showMessage('')
        self.menuBar.setNativeMenuBar(False)

        # Init menu actions
        self.menu_file = GUIMenuFile(self)
        self.menu_edit = GUIMenuEdit(self)
        self.menu_sim = GUIMenuSimulation(self)

        self.init_menu_actions()

        # Centering main window
        self.centering_window()
        self.init_central_widget()


    def centering_window(self):

        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())


    def init_menu_actions(self):

        # File menu
        self.action_new_lattice.triggered.connect(self.menu_file.new_lattice)
        self.action_open_lattice.triggered.connect(self.menu_file.dialog_open_lattice)
        self.action_save_lattice.triggered.connect(self.menu_file.dialog_save_lattice)
        self.action_exit.triggered.connect(qApp.quit)

        # Edit menu
        self.action_edit_lattice.triggered.connect(self.menu_edit.edit_lattice)

        # Simulation menu
        self.action_calc_twiss.triggered.connect(self.menu_sim.calc_twiss)


    def init_central_widget(self):
        """Central widget - grid layout"""

        self.layout = QtWidgets.QGridLayout()
        self.central_widget.setLayout(self.layout)


    def clean_central_widget(self):

        while self.layout.count():
            child = self.layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()


    def loadStyleSheet(self, filename="colinDark.css"):
        """
        Sets the dark GUI theme from a css file.

        :return:
        """
        root_dir = os.path.dirname(os.path.abspath(__file__))
        style_dir = os.path.join(root_dir, "ui_forms/ui_style/")
        try:
            self.cssfile = os.path.join(style_dir, filename)
            with open(self.cssfile, "r") as f:
                self.setStyleSheet(f.read())
        except IOError:
            print ('No style sheet found!')
        
    
    def error_window(self, title, msg, txt='Error'):
        """Error dialog window"""
    
        window = QMessageBox()
        window.setIcon(QMessageBox.Critical)
        window.setText(txt)
        window.setInformativeText(msg)
        window.setWindowTitle(title)
        window.exec()
        