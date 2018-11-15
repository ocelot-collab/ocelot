"""Menu Edit Lattice description"""

from copy import deepcopy

from app.ui_forms.edit_lattice import *
from app.parser import *


class EditLattice(QtWidgets.QWidget):

    def __init__(self, MainWindow):
        super().__init__()
        
        # init user interface
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        
        self.mw = MainWindow

        # init button actions
        self.ui.btn1.clicked.connect(self.action_update_parameters)
        self.ui.btn2.clicked.connect(self.action_cancel)
        
        self.load_parameters()
    
    
    def load_parameters(self):
        """Load all parameters to window"""

        parcer = LatticeSave()

        # elements
        lines = parcer.parsing_elements(self.mw.lattice.elements)
        lines = lines[0:-2]
        self.ui.edit_elements.setText(lines)

        # sequences
        lines = parcer.parsing_cell(self.mw.lattice.cell)
        lines = lines[0:-2]
        self.ui.edit_cells.setText(lines)

        # number of superperiods
        self.ui.edit_nsuperperiods.setValue(self.mw.lattice.nsuperperiods)

        # beam
        lines = parcer.parsing_beam(self.mw.lattice.beam)
        lines = lines[0:-2]
        self.ui.edit_beam.setText(lines)

        # tws
        lines = parcer.parsing_twiss(self.mw.lattice.tws0)
        lines = lines[0:-2]
        self.ui.edit_twiss.setText(lines)
        
        # periodic solution
        if self.mw.lattice.periodic_solution:
            self.ui.edit_periodic_solution.toggle()
            
            
    def action_update_parameters(self):
        
        parcer = LatticeRead()
        
        # elements and sequences
        lines = self.ui.edit_elements.toPlainText()
        lines += '\n'
        lines += self.ui.edit_cells.toPlainText()
        
        loc_dict = {}
        try:
            exec(lines, globals(), loc_dict)
        except Exception as err:
            self.mw.error_window('Edit Lattice Error', str(err))
            return
            
        if 'cell' in loc_dict:
            self.mw.lattice.cell = deepcopy(loc_dict['cell'])
            self.mw.lattice.elements = parcer.parsing_elements(self.mw.lattice.cell)
        else:
            self.mw.lattice.cell = ()
            self.mw.lattice.elements = {}
            self.mw.error_window('Edit Lattice Error', 'NO secuence named "cell"')
        
        # beam
        lines = self.ui.edit_beam.toPlainText()
        
        loc_dict = {}
        try:
            exec(lines, globals(), loc_dict)
        except Exception as err:
            self.mw.error_window('Edit Lattice Error', str(err))
            
        if 'beam' in loc_dict:
            self.mw.lattice.beam = deepcopy(loc_dict['beam'])
        else:
            self.mw.lattice.beam = Beam()
        
        # tws
        lines = self.ui.edit_twiss.toPlainText()
        
        loc_dict = {}
        try:
            exec(lines, globals(), loc_dict)
        except Exception as err:
            self.mw.error_window('Edit Lattice Error', str(err))
            
        if 'tws0' in loc_dict:
            self.mw.lattice.tws0 = deepcopy(loc_dict['tws0'])
        else:
            self.mw.lattice.tws0 = Twiss()
        
        # number of superperiods
        try:
            val = int(self.ui.edit_nsuperperiods.value())
            if val > 0:
                self.mw.lattice.nsuperperiods = val
        except ValueError:
            self.mw.lattice.nsuperperiods = 1

        # periodic solution
        if self.ui.edit_periodic_solution.checkState():
            self.mw.lattice.periodic_solution = True
        else:
            self.mw.lattice.periodic_solution = False
            
            
    def action_cancel(self):
        self.mw.init_central_widget()
    