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

        super().__init__()
        self.setupUi(self)
        
        self.initUI()
        #self.loadStyleSheet()
        
        self.tunable_elements = {'Bend':['angle', 'k1'], 'SBend':['angle', 'k1'], 'RBend':['angle', 'k1'], 'Quadrupole':['k1'], 'Drift':['l']}
        self.matchable_elements = {'Bend':'k1','SBend':'k1', 'RBend':'k1', 'Quadrupole':'k1', 'Drift':'l'}
        
        self.lattice = GUILattice()


    def initUI(self):

        # Init main window
        self.statusBar().showMessage('')
        self.setWindowTitle('Ocelot GUI')
        self.centering_window()

        # Init menu actions
        self.menu_file = GUIMenuFile(self)
        self.menu_edit = GUIMenuEdit(self)
        self.menu_sim = GUIMenuSimulation(self)
        
        # Init menu
        self.create_menu()

        
    def centering_window(self):

        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())


    def create_menu(self):
        
        menubar = self.menuBar()

        # File menu
        action_new_lattice = QAction('New', self)
        action_new_lattice.triggered.connect(self.menu_file.new_lattice)

        action_open_lattice = QAction('Open lattice', self)
        action_open_lattice.triggered.connect(self.menu_file.dialog_open_lattice)

        action_save_lattice = QAction('Save lattice', self)
        action_save_lattice.triggered.connect(self.menu_file.dialog_save_lattice)
        
        action_exit = QAction('Exit', self)
        action_exit.triggered.connect(qApp.quit)
        
        menu_file = menubar.addMenu('File')
        menu_file.addAction(action_new_lattice)
        menu_file.addAction(action_open_lattice)
        menu_file.addAction(action_save_lattice)
        menu_file.addAction(action_exit)

        # Edit menu
        action_edit_lattice = QAction('Edit lattice and parameters', self)
        action_edit_lattice.triggered.connect(self.menu_edit.edit_lattice)

        menu_edit = menubar.addMenu('Edit')
        menu_edit.addAction(action_edit_lattice)

        # Simulation menu
        action_calc_twiss = QAction('Twiss functions', self)
        action_calc_twiss.triggered.connect(self.menu_sim.calc_twiss)
        
        action_calc_params = QAction('Main parameters', self)
        action_calc_params.triggered.connect(self.menu_sim.calc_params)

        action_calc_matching = QAction('Matching', self)
        action_calc_matching.triggered.connect(self.menu_sim.matching)

        menu_sim = menubar.addMenu('Simulation')
        menu_sim.addAction(action_calc_twiss)
        menu_sim.addAction(action_calc_params)
        menu_sim.addAction(action_calc_matching)


    def init_central_widget(self):
        """Central widget - grid layout"""

        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)
        self.layout = QtWidgets.QGridLayout(central_widget)
        
    
    def error_window(self, title, msg, txt='Error'):
        """Error dialog window"""
    
        window = QMessageBox()
        window.setIcon(QMessageBox.Critical)
        window.setText(txt)
        window.setInformativeText(msg)
        window.setWindowTitle(title)
        window.exec()
        
    
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
