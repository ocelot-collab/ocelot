"""GUI description"""

from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QDesktopWidget, QMainWindow, QAction

from app.gui_menu_file import *
from app.gui_menu_edit import *
from app.gui_menu_simulation import *
from app.lattice import *


class GUIWindow(QMainWindow):

    def __init__(self):

        super().__init__()
        self.initUI()
        self.tunable_elements = {'Bend':['angle', 'k1'],'SBend':['angle', 'k1'], 'RBend':['angle', 'k1'], 'Quadrupole':['k1'], 'Drift':['l']}
        self.matchable_elements = {'Bend':['k1'],'SBend':['k1'], 'RBend':['k1'], 'Quadrupole':['k1'], 'Drift':['l']}
        self.lattice = GUILattice()


    def initUI(self):

        self.statusBar().showMessage('')
        self.setWindowTitle('Ocelot GUI')
        self.setGeometry(0, 0, 1300, 900)
        self.centering_window()

        self.menubar = self.menuBar()
        self.menu_file = GUIMenuFile(self)
        self.menu_edit = GUIMenuEdit(self)
        self.menu_sim = GUIMenuSimulation(self)
        self.create_menu()
        
        self.init_central_widget()  
        #self.show()

        
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
        action_calc_twiss = QAction('Main parameters and Twiss functions', self)
        action_calc_twiss.triggered.connect(self.menu_sim.calc_parameters)

        action_calc_matching = QAction('Matching', self)
        action_calc_matching.triggered.connect(self.menu_sim.matching)

        menu_sim = menubar.addMenu('Simulation')
        menu_sim.addAction(action_calc_twiss)
        menu_sim.addAction(action_calc_matching)


    def init_central_widget(self):
        """Central widget - grid layout"""

        central_widget = QtWidgets.QWidget()
        self.setCentralWidget(central_widget)
        self.layout = QtWidgets.QGridLayout(central_widget)
