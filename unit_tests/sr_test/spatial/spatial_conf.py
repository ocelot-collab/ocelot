"""Test parameters description file"""

from ocelot.cpbd.beam.core import Beam
from ocelot.cpbd.elements.undulator import Undulator
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.rad.screen import Screen


import pytest



"""Lattice elements definition"""

und = Undulator(Kx=0.43, nperiods=500, lperiod=0.007, eid="und")


"""pytest fixtures definition"""
    
@pytest.fixture(scope='module')
def lattice():
    return MagneticLattice((und))


@pytest.fixture(scope='function')
def beam():

    b = Beam()
    b.E = 2.5
    b.I = 0.1
    b.beta_x = 12.84
    b.beta_y = 6.11
    b.Dx = 0.526

    return b


@pytest.fixture(scope='function')
def screen():

    s = Screen()
    s.z = 100.0

    return s
