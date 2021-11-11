"""Test parameters description"""

import pytest
import numpy as np

from ocelot import *
from ocelot import XYQuadrupole

"""lattice elements description"""

D = Drift(l=1, eid="D")
Qf = Quadrupole(l=0.3, k1=1, eid="Qf")
Qd = Quadrupole(l=0.3, k1=-1, eid="Qd")
qxy = XYQuadrupole(l=0.5, x_offs=0.01, y_offs=0.025, k1=-0.5, eid='QK')
octf = Octupole(l=0.3,  k3=300, eid="of")
octd = Octupole(l=0.3,  k3=-300, eid="of")


"""pytest fixtures description"""

@pytest.fixture(scope='module')
def cell():
    return (D, Qf, D, Qd, D, octf, D, qxy, D, octd, D, Qf, D, Qd, D)


@pytest.fixture(scope='module')
def lattice(cell):
    return MagneticLattice(cell, method={"global": SecondTM, Octupole: KickTM, "nkick": 5})

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
                              nparticles=1000, energy=0.13)

    return p_array