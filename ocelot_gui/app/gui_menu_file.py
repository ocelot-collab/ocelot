"""Menu File description"""

from PyQt5.QtWidgets import qApp, QFileDialog

from app.parser import *
from app.lattice import *


class GUIMenuFile():

    def __init__(self, MainWindow):
        self.mw = MainWindow


    def new_lattice(self):
        """Reset and init new lattice"""
        
        self.mw.init_central_widget()
        self.mw.lattice = GUILattice()
        self.mw.menu_edit.edit_lattice()


    def dialog_open_lattice(self):
        """Read lattice from file"""

        filename = QFileDialog.getOpenFileName(self.mw, 'Open Lattice', '', "Lattice Files (*.inp);;All Files (*)", options=QFileDialog.DontUseNativeDialog)[0]
        if filename == '':
            return 0

        with open(filename) as fp:
            lines = fp.read()

        lp = LatticeRead()
        self.mw.lattice = GUILattice()
        self.mw.lattice.tws, self.mw.lattice.beam, self.mw.lattice.elements, self.mw.lattice.cells, self.mw.lattice.cells_order = lp.parsing(lines)
        self.mw.menu_edit.edit_lattice()


    def dialog_save_lattice(self):
        """Save lattice to file"""

        filename = QFileDialog.getSaveFileName(self.mw, 'Save Lattice', '', "Lattice Files (*.inp);;All Files (*)", options=QFileDialog.DontUseNativeDialog)
        if filename[0] == '':
            return 0

        ls = LatticeSave()
        lines = ls.parsing(self.mw.lattice)

        if filename[1] == 'Lattice Files (*.inp)' and filename[0][-4:] != '.inp':
            filename = filename[0] + '.inp'
        else:
            filename = filename[0]

        with open(filename, 'w') as fp:
            fp.write(lines)
