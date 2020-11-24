"""Calculation twiss parameters description"""

from ocelot.cpbd.transformations.optics import *


class TwissParameters():

    def __init__(self, lattice):
        
        self.lattice = lattice
        self.tws = None
        self.tws_step = 0.0


    def __del__(self):
        pass


    def calc_twiss(self):

        if self.lattice.periodic_solution:
            tws0 = periodic_twiss(Twiss(), lattice_transfer_map(self.lattice.lattice, self.lattice.tws0.E))
        else:
            tws0 = self.lattice.tws0

        if tws0 is not None:

            # fix only for periodic solution
            tws0.E = self.lattice.tws0.E

            if self.tws_step != 0.0:
                n = int(self.lattice.lattice.totalLen / self.tws_step) + 1
            else:
                n = None

            self.tws = twiss(self.lattice.lattice, tws0=tws0, nPoints=n)
        else:
            self.tws = None
