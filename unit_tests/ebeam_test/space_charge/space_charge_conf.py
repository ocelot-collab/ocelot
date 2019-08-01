"""Test parameters description file"""

import pytest
import numpy as np

from ocelot import *

"""Lattice elements defenition"""

D_14 = Drift(l=0.2216+0.0996, eid='D_14')
D_15 = Drift(l=0.3459, eid='D_15')
D_22 = Drift(l=0.2043, eid='D_22')
D_23 = Drift(l=0.085+0.4579+0.2211+0.085, eid='D_23')

phi1 = 18.7268
V1 = 18.50662e-3 / np.cos(phi1 * np.pi / 180.0)

C_A1_1_1_I1 = Cavity(l=1.0377, v=V1, freq=1.3e9, phi=phi1, eid='C.A1.1.1.I1')
C_A1_1_2_I1 = Cavity(l=1.0377, v=V1, freq=1.3e9, phi=phi1, eid='C.A1.1.2.I1')
C_A1_1_3_I1 = Cavity(l=1.0377, v=V1, freq=1.3e9, phi=phi1, eid='C.A1.1.3.I1')
C_A1_1_4_I1 = Cavity(l=1.0377, v=V1, freq=1.3e9, phi=phi1, eid='C.A1.1.4.I1')
C_A1_1_5_I1 = Cavity(l=1.0377, v=V1, freq=1.3e9, phi=phi1, eid='C.A1.1.5.I1')
C_A1_1_6_I1 = Cavity(l=1.0377, v=V1, freq=1.3e9, phi=phi1, eid='C.A1.1.6.I1')
C_A1_1_7_I1 = Cavity(l=1.0377, v=V1, freq=1.3e9, phi=phi1, eid='C.A1.1.7.I1')
C_A1_1_8_I1 = Cavity(l=1.0377, v=V1, freq=1.3e9, phi=phi1, eid='C.A1.1.8.I1')

phi13 = 180.0
V13 = -20.2E-3 / 8.0 / np.cos(phi13 * np.pi / 180.0)
C3_AH1_1_1_I1 = Cavity(l=0.346, v=V13, freq=3.9e9, phi=phi13, eid='C3.AH1.1.1.I1')
C3_AH1_1_2_I1 = Cavity(l=0.346, v=V13, freq=3.9e9, phi=phi13, eid='C3.AH1.1.2.I1')
C3_AH1_1_3_I1 = Cavity(l=0.346, v=V13, freq=3.9e9, phi=phi13, eid='C3.AH1.1.3.I1')
C3_AH1_1_4_I1 = Cavity(l=0.346, v=V13, freq=3.9e9, phi=phi13, eid='C3.AH1.1.4.I1')
C3_AH1_1_5_I1 = Cavity(l=0.346, v=V13, freq=3.9e9, phi=phi13, eid='C3.AH1.1.5.I1')
C3_AH1_1_6_I1 = Cavity(l=0.346, v=V13, freq=3.9e9, phi=phi13, eid='C3.AH1.1.6.I1')
C3_AH1_1_7_I1 = Cavity(l=0.346, v=V13, freq=3.9e9, phi=phi13, eid='C3.AH1.1.7.I1')
C3_AH1_1_8_I1 = Cavity(l=0.346, v=V13, freq=3.9e9, phi=phi13, eid='C3.AH1.1.8.I1')

Q_37_I1 = Quadrupole(l=0.3, k1=-1.537886, tilt=0.0, eid='Q.37.I1')
Q_38_I1 = Quadrupole(l=0.3, k1=1.435078, tilt=0.0, eid='Q.38.I1')


"""pytest fixtures defenition"""

@pytest.fixture(scope='module')
def cell():
    return (Marker(), D_14, C_A1_1_1_I1, D_15, C_A1_1_2_I1, D_15, C_A1_1_3_I1, D_15, C_A1_1_4_I1, D_15, C_A1_1_5_I1, D_15, C_A1_1_6_I1, D_15, C_A1_1_7_I1, D_15, C_A1_1_8_I1, D_22, Q_37_I1, D_23, Q_38_I1)


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

    sigma_x = 0.000231507245956
    sigma_px =0.000204206874319
    sigma_y = 0.000231583942392
    sigma_py =0.000204272734636
    n = 20000
    x = np.random.randn(n) * sigma_x
    px = np.random.randn(n) * sigma_px
    y = np.random.randn(n) * sigma_y
    py = np.random.randn(n) * sigma_py

    cov_t_p =  [[  6.89508231e-07,  -2.98688604e-07], [ -2.98688604e-07,   1.87434257e-07]]

    long_dist = np.random.multivariate_normal((0, 0), cov_t_p, n)
    tau = long_dist[:, 0]
    dp = long_dist[:, 1]

    p_array = ParticleArray(n=n)
    p_array.E = 0.0065 # GeV
    p_array.rparticles[0] = x
    p_array.rparticles[1] = px
    p_array.rparticles[2] = y
    p_array.rparticles[3] = py
    p_array.rparticles[4] = tau
    p_array.rparticles[5] = dp

    # creating charge array
    Q = 5e-10
    p_array.q_array = np.ones(n) * Q / n

    # beam transformation
    tws = Twiss()
    tws.beta_x  = 1.59966676201
    tws.beta_y  = 1.60018325757
    tws.alpha_x = -0.995487979563
    tws.alpha_y = -0.996116091572
    tws.mux = 0.0
    tws.muy = 0.0

    bt = BeamTransform(tws=tws)
    bt.apply(p_array, dz=0.0)

    return p_array