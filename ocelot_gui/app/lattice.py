"""Lattice edit functions"""

from ocelot import *


class GUILattice():

    def __init__(self):

        self.tws0 = Twiss()
        self.beam = Beam()
        
        self.elements = {}
        self.cell = ()

        self.nsuperperiods = 1

        # use for first order matrix tracking
        self.method = {'global': TransferMap, 'Sextupole': KickTM}

        # use for second order matrix tracking
        # self.method = {'global': SecondTM, 'Sextupole': KickTM}

        self.lattice = None
        self.periodic_solution = False

        self.tunable_elements = {'Bend':['angle', 'k1'],'SBend':['angle', 'k1'], 'RBend':['angle', 'k1'], 'Quadrupole':['k1'], 'Drift':['l']}
        self.matchable_elements = {'Bend':'k1','SBend':'k1', 'RBend':'k1', 'Quadrupole':'k1', 'Drift':'l'}


    def init_lattice(self):
        """Create magnetic lattice from the last cell sequences"""
        
        self.init_tuneablity()
        self.init_matchablity()
        self.lattice = MagneticLattice(self.cell,  method=self.method)


    def init_tuneablity(self):
        """Set possibility to tune element"""

        for elem_id, elem in self.elements.items():
            if elem.__class__.__name__ in self.tunable_elements:
                elem.is_tuneable = True
                elem.tune_params = self.tunable_elements[elem.__class__.__name__]
            else:
                elem.is_tuneable = False


    def init_matchablity(self):
        """Set possibility to match element"""

        for elem_id, elem in self.elements.items():
            if elem.__class__.__name__ in self.matchable_elements:
                elem.is_matchable = True
                elem.match_params = self.matchable_elements[elem.__class__.__name__]
            else:
                elem.is_matchable = False
