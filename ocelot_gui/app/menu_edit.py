"""Menu Edit description"""

from copy import deepcopy

from app.action_edit_lattice import *


class GUIMenuEdit():

    def __init__(self, MainWindow):
        self.mw = MainWindow


    def edit_lattice(self):
        
        self.mw.init_central_widget()
        
        edit_lattice_form = EditLattice(self.mw)
        self.mw.layout.addWidget(edit_lattice_form, 0, 0)
