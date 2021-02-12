"""Test parameters description file"""

import pytest
import numpy as np

from ocelot import *
from ocelot.cpbd.beam import generate_parray
#from ocelot.utils.acc_utils import chicane_RTU
#from ocelot.cpbd.sc import LSC


class LogProc(PhysProc):
    def __init__(self):
        PhysProc.__init__(self)
        self.dz_list = []
        self.totalLen = 0.
        self.L_list = []
        self.s_list = []
        self.s_stop_list = []
        self.s_start_list = []

    def prepare(self, lat):
        self.totalLen = lat.totalLen

    def apply(self, p_array, dz):
        self.dz_list.append(dz)
        L = self.s_stop - self.s_start
        self.L_list.append(L)
        self.s_list.append(p_array.s)
        self.s_stop_list.append(self.s_stop)
        self.s_start_list.append(self.s_start)



"""Lattice elements definition"""

start = Marker()
m_extra1 = Marker()
D0 = Drift(l=0.3, eid='D0')

m_kick = Marker()

B = SBend(l=0.25, angle=-0.099484, e1=0.0, e2 = -0.099484, tilt=0.0, fint=0.0, eid='BL.48I.I1')

D1 = Drift(l=0.4, eid='D1')
ap = Aperture(xmax=0.01, ymax=0.01)
m1 = Marker()
D2 = Drift(l=0.5)
m2 = Marker()
m_extra2 = Marker()
stop = Marker()


"""pytest fixtures definition"""

@pytest.fixture(scope='module')
def cell():
    return (start, m_extra1, D0, m_kick, B, ap, D1, m1, D2, m2, m_extra2, stop)


@pytest.fixture(scope='module')
def method():
    
    m = MethodTM()
    m.global_method = SecondTM
    
    return m
    
    
@pytest.fixture(scope='module')
def lattice(cell, method):
    return MagneticLattice(cell, method=method)


@pytest.fixture(scope='function')
def p_array():

    np.random.seed(11)

    # generate beam file
    sigma_x = 0.000121407185261
    sigma_px = 1.80989470506e-05
    sigma_y = 0.000165584800564
    sigma_py = 4.00994225888e-05

    p_array = generate_parray(sigma_x=sigma_x, sigma_px=sigma_px, sigma_y=sigma_y, sigma_py=sigma_py,
                              sigma_tau=1.30190131e-04, sigma_p=3.09815718e-04, chirp=0.002, charge=0.5e-9,
                              nparticles=10, energy=0.13)

    p_array.rparticles *= 0.

    return p_array