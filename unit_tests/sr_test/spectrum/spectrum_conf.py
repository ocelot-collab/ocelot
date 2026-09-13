"""Test parameters description file"""

from ocelot.cpbd.beam.core import Beam
from ocelot.cpbd.elements.undulator import Undulator
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.rad.screen import Screen


import pytest



"""Lattice elements defenition"""

und = Undulator(Kx=0.43, nperiods=280, lperiod=0.007, eid="und")


"""pytest fixtures defenition"""
    
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
    s.size_x = 0.0
    s.size_y = 0.0
    s.nx = 1
    s.ny = 1
    s.start_energy = 7300
    s.end_energy = 7900
    s.num_energy = 101

    return s
