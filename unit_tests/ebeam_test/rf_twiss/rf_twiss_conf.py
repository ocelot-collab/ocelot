"""Test parameters description"""

import pytest
import numpy as np

from ocelot import *

"""lattice elements description"""

phi_bc2 = 0.033646252962410

l_bc2 = 0.5 * phi_bc2 / np.sin(phi_bc2)

ac_v = 0.02265625 # in GV

bb_393_b2 = Bend(l=l_bc2, angle=phi_bc2, e1=0.0, e2=phi_bc2, tilt=1.570796330, eid= 'bb_393_b2')
bb_402_b2 = Bend(l=l_bc2, angle=-phi_bc2, e1=-phi_bc2, e2=0.0, tilt=1.570796330, eid= 'bb_402_b2')
bb_404_b2 = Bend(l=l_bc2, angle=-phi_bc2, e1=0.0, e2=-phi_bc2, tilt=1.570796330,  eid= 'bb_404_b2')
bb_413_b2 = Bend(l=l_bc2, angle=phi_bc2, e1=phi_bc2, e2=0.0, tilt=1.570796330,  eid= 'bb_413_b2')

d10cm =   Drift(l=0.1, eid= 'd10cm')
cd850cm = Drift(l=8.5 / np.cos(phi_bc2), eid= 'cd850cm')
cd150cm = Drift(l=1.5, eid= 'cd150cm')
cd100cm = Drift(l=1.0, eid= 'cd100cm')
d34cm59 = Drift(l=0.3459, eid= 'd34cm59')
d13cm =   Drift(l=0.13, eid= 'd13cm')
d130cm =  Drift(l=1.3, eid= 'd130cm')

qd_415_b2 = Quadrupole(l=0.2, k1=0.3, tilt=0.0, eid= 'qd_415_b2')
qd_417_b2 = Quadrupole(l=0.2, k1=-0.2, tilt=0.0, eid= 'qd_417_b2')
qd_418_b2 = Quadrupole(l=0.2, k1=-0.5, tilt=0.0, eid= 'qd_418_b2')
q_249_l2 =  Quadrupole(l=0.3, k1=0.25, tilt=0.0, eid= 'q_249_l2')
q_261_l2 =  Quadrupole(l=0.3, k1=-0.29711100, tilt=0.0, eid= 'q_261_l2')

c_a3 = Cavity(l=1.0377000, phi=0.0, v = ac_v, freq=1.300e+009, eid= 'c_a3')


"""pytest fixtures descripteion"""

@pytest.fixture(scope='module')
def cell():

    bc2 = (d10cm, bb_393_b2, cd850cm, bb_402_b2, cd150cm, bb_404_b2, cd850cm, bb_413_b2, cd100cm)

    l3 = (d13cm, qd_415_b2, d130cm, qd_417_b2, d130cm, qd_418_b2, d130cm,
          c_a3, d34cm59, c_a3, d34cm59, c_a3, d34cm59, c_a3, d34cm59,
          c_a3, d34cm59, c_a3, d34cm59, c_a3, d34cm59, c_a3, d13cm,
          q_249_l2, d34cm59, c_a3, d34cm59, c_a3, d34cm59, c_a3, d34cm59, c_a3, d34cm59,
          c_a3, d34cm59, c_a3, d34cm59, c_a3, d34cm59, c_a3, d13cm, q_261_l2, d130cm)

    return (bc2, l3)


@pytest.fixture(scope='module')
def lattice(cell):
    return MagneticLattice(cell, method=MethodTM({'global': SecondTM}))


@pytest.fixture(scope='function')
def p_array():

    np.random.seed(11)
    n = 5000
    # generate beam file
    sigma_x = 0.000121407185261
    sigma_px = 1.80989470506e-05
    sigma_y = 0.000165584800564
    sigma_py = 4.00994225888e-05

    x = np.random.randn(n) * sigma_x
    px = np.random.randn(n) * sigma_px
    y = np.random.randn(n) * sigma_y
    py = np.random.randn(n) * sigma_py

    # covariance matrix for [tau, p] for beam compression in BC
    cov_t_p = [[1.30190131e-06, 2.00819771e-05], [2.00819771e-05, 3.09815718e-04]]
    long_dist = np.random.multivariate_normal((0,0), cov_t_p, n)
    tau = long_dist[:, 0]
    dp = long_dist[:, 1]

    p_array = ParticleArray(n=n)
    p_array.E = 0.130 # GeV
    p_array.rparticles[0] = x
    p_array.rparticles[1] = px
    p_array.rparticles[2] = y
    p_array.rparticles[3] = py
    p_array.rparticles[4] = tau
    p_array.rparticles[5] = dp

    Q = 5e-9

    p_array.q_array = np.ones(n) * Q / n

    return p_array