"""Lattice edit functions"""

from ocelot import *


class GUILattice():

    def __init__(self):

        self.tws = Twiss()
        #self.tws.beta_x = 1.0
        #self.tws.beta_y = 1.0
        
        self.beam = Beam()
        #self.beam.E = 1.0
        
        self.elements = {}
        self.cells = {}
        self.cells_order = []

        self.nsuperperiods = 1
        
        self.method = MethodTM()
        self.method.global_method = TransferMap
        self.method.params[Sextupole] = KickTM
        #self.method.global_method = SecondTM

        self.lattice = None
        self.lattice_ename_sequence = []
        self.periodic_solution = False

        self.init_lattice()


    def init_lattice(self):
        """Create magnetic lattice from the last cell sequences"""

        sequence = []

        if self.cells_order != [] and self.cells_order[-1] in self.cells:
            for cell_element in self.cells[self.cells_order[-1]]:
                if cell_element in self.elements.keys():
                    sequence.append(self.elements[cell_element])
                    self.lattice_ename_sequence.append(cell_element)
                        
                if cell_element in self.cells.keys():
                    sequence.extend(self._convert_cell(self.cells[cell_element]))

        self.lattice = MagneticLattice(tuple(sequence),  method=self.method)


    def _convert_cell(self, cell):

        sequence = []
        for elem in cell:
            if elem in self.cells.keys():
                sequence.extend(self._convert_cell(self.cells[elem]))
            else:
                sequence.append(self.elements[elem])
                self.lattice_ename_sequence.append(elem)
        
        return sequence
