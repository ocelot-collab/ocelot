"""Lattice edit functions"""

from ocelot import *


class GUILattice():

    def __init__(self):

        self.tws0 = Twiss()
        self.beam = Beam()
        
        self.elements = {}
        self.cell = ()

        self.nsuperperiods = 1
        
        self.method = MethodTM()
        self.method.global_method = TransferMap
        self.method.params[Sextupole] = KickTM
        #self.method.global_method = SecondTM

        self.lattice = None
        self.periodic_solution = False


    def init_lattice(self):
        """Create magnetic lattice from the last cell sequences"""
        
        self.lattice = MagneticLattice(self.cell,  method=self.method)
