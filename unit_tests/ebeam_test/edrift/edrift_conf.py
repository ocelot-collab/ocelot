__author__ = "Nikita Kuklev"

import pytest
from ocelot import *

"""
Have 2 basic test configs: just drift or drift + symmetric lens

Latter is defined by any 2 of these parameters:
    mu = phase advance
    l = length
    f = focal length
    beta_star or beta_edge
    alpha_edge
"""

l = 2.0
mu = 0.25

# Solutions are analytic, saving reference set is not necessary
f = l / 4 / np.sin(mu * np.pi) ** 2
beta_star = l / 2 / np.tan(np.pi*mu)
beta_edge = beta_star + (l/2)**2 / beta_star
alpha_edge = - l / 2 / beta_star

d_matrix = np.eye(6)
d_matrix[0, 1] = d_matrix[2, 3] = l / 2

m_matrix = np.eye(6)
m_matrix[1, 0] = m_matrix[3, 2] = -1.0 / f / 2.0  # R21 and R43 for symmetric focusing


@pytest.fixture(scope='function')
def method():
    return MethodTM()


@pytest.fixture(scope='function')
def lattice(method):
    d = Drift(l=l/2, eid='D')
    m = Matrix(eid='M')
    m.r = m_matrix.copy()
    return MagneticLattice((m, d, d, m), method=method)


@pytest.fixture(scope='function')
def lattice_exact(method, edrift_method=1):
    d = EDrift(l=l/2, eid='D', method=edrift_method)
    m = Matrix(eid='M')
    m.r = m_matrix.copy()
    return MagneticLattice((m, d, d, m), method=method)


@pytest.fixture(scope='function')
def drift(method):
    d = Drift(l=l/2, eid='D')
    return MagneticLattice((d, d), method=method)


@pytest.fixture(scope='function')
def drift_exact(method, edrift_method=1):
    de = EDrift(l=l/2, eid='D', method=edrift_method)
    return MagneticLattice((de, de), method=method)
