"""Test parameters description file"""

import pytest
import numpy as np

from ocelot import *
from ocelot.cpbd.r_matrix import rot_mtx
from ocelot.utils.acc_utils import chicane_RTU

"""Lattice elements definition"""

b1 = Bend(l=0.501471, angle=0.132729704703, e1=0.0, e2=0.132729704703,   tilt=0.0, eid="b")
b2 = Bend(l=0.501471, angle=-0.132729704703, e1=-0.132729704703, e2=0.0,  tilt=0.0, eid="b")
b3 = Bend(l=0.501471, angle=-0.132729704703, e1=0.0, e2=-0.132729704703,  tilt=0.0, eid="b")
b4 = Bend(l=0.501471, angle=0.132729704703, e1=0.132729704703, e2=0.0,   tilt=0.0, eid="b")
d1 = Drift(l=1.0)
d2 = Drift(l=1.5)


"""pytest fixtures definition"""

@pytest.fixture(scope='module')
def cell():
    return (Marker(eid="m1"), Drift(l=0.1), b1, d1, b2, d2, b3, d1, b4, d1, Marker(eid="m2"))


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