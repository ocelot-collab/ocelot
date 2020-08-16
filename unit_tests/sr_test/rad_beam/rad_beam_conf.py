"""Test parameters description file"""

import pytest

from ocelot import *
from ocelot.rad import *
import numpy as np


"""Lattice elements defenition"""

und = Undulator(lperiod=0.4, nperiods=9, Kx=44.81)

"""pytest fixtures defenition"""
    
@pytest.fixture(scope='module')
def lattice():
    return MagneticLattice((und))


@pytest.fixture(scope='function')
def beam():

    chirp_coeff = 0.01 / 2.36
    sigma_tau = 100e-6 / 2.36

    tau = np.array([-0.7, -0.3, 0, 0.3, 0.7]) * sigma_tau

    p_array = ParticleArray(n=5)
    p_array.E = 0.6
    p_array.rparticles[4, :] = tau[:]
    p_array.rparticles[5, :] = chirp_coeff * tau / sigma_tau
    p_array.q_array[:] = 1e-10

    return p_array


@pytest.fixture(scope='function')
def screen():
    s = Screen()
    s.z = 1000.0
    s.size_x = 15
    s.size_y = 15
    s.nx = 101
    s.ny = 1
    s.start_energy = 0.00850446  # eV
    s.end_energy = 15e-3  # eV
    s.num_energy = 1

    return s
