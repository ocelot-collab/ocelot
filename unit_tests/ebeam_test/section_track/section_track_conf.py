"""Test parameters description file"""

import pytest
import numpy as np

from ocelot import *
from accelerator.s2e_sections.sections import *
from ocelot.utils.section_track import *
import time

data_dir = "./unit_tests/ebeam_test/section_track/data"



@pytest.fixture(scope='module')
def all_sections():
    return [A1, AH1, LH, DL, BC0, L1, BC1]


@pytest.fixture(scope='module')
def tws0():
    tws = Twiss()
    tws.E = 0.005
    tws.beta_x = 0.286527307369
    tws.beta_y = 0.286527307369
    tws.alpha_x = -0.838833736086
    tws.alpha_y = -0.838833736086
    return tws


@pytest.fixture(scope='module')
def section_lat(all_sections, tws0):
    lats = SectionLattice(sequence=all_sections, tws0=tws0, data_dir=data_dir)
    return lats


@pytest.fixture(scope='function')
def p_array():

    np.random.seed(11)
    n = 20000
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