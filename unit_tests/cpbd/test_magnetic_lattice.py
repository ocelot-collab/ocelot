import pytest
import numpy as np

from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.elements import Drift, SBend, Quadrupole, Cavity
from ocelot.cpbd.beam import Twiss


def test_magnetic_lattice_find_indices():
    d1 = Drift(l=0.1)
    sb = SBend(l=0.2, angle=3)
    d2 = Drift(l=0.2)
    sequence = [d1, sb, d2, sb]
    mlat = MagneticLattice(sequence)
    indices = mlat.find_indices(SBend)
    assert indices == [1, 3]


def test_magnetic_lattice_find_indices_by_predicate():
    d1 = Drift(l=0.1)
    sb = SBend(l=0.2, angle=3)
    d2 = Drift(l=0.2)

    sequence = [d1, sb, d2, sb]
    mlat = MagneticLattice(sequence)
    def pred(element):
        return element.l == 0.2 and isinstance(element, SBend)

    indices = mlat.find_indices_by_predicate(pred)
    assert indices == [1, 3]


def test_magnetic_lattice_periodic_twiss_with_cav():
    d2 = Drift(l=0.5)
    qf = Quadrupole(l=0.2, k1=0.3)
    qdh = Quadrupole(l=0.2 / 2, k1=-0.3)
    c = Cavity(l=0.5, v=0.01, freq=1.3e9, phi=0)
    cell = (qdh, d2, c, d2, qf, d2, c, d2, qdh)
    lat = MagneticLattice(cell)

    tws0 = Twiss()
    tws0.E = 0.005
    tws0.beta_x = 10
    tws0.beta_y = 10
    tw_p = lat.periodic_twiss(tws=Twiss(E=0.005))
    beta_x = 1.60031832321874
    beta_y = 1.8790818475215212
    alpha_x = -0.8882452561429366
    alpha_y = -0.9590585725276709
    gamma_x = 1.1178898654751601
    gamma_y = 1.0216656331766445
    E = 0.005
    s = 0.0
    ref = np.array([beta_x, beta_y, alpha_x, alpha_y, gamma_x, gamma_y, E, s])
    res = np.array([tw_p.beta_x, tw_p.beta_y, tw_p.alpha_x, tw_p.alpha_y, tw_p.gamma_x, tw_p.gamma_y, tw_p.E, tw_p.s])
    assert (ref == res).all()

def test_magnetic_lattice_periodic_twiss():
    d2 = Drift(l=0.5)
    qf = Quadrupole(l=0.2, k1=0.3)
    qdh = Quadrupole(l=0.2 / 2, k1=-0.3)
    c = Drift(l=0.5)
    cell = (qdh, d2, c, d2, qf, d2, c, d2, qdh)
    lat = MagneticLattice(cell)

    tws0 = Twiss()
    tws0.E = 0.005
    tws0.beta_x = 10
    tws0.beta_y = 10
    tw_p = lat.periodic_twiss(tws=Twiss(E=0.005))
    beta_x = 33.09417032395688
    beta_y = 36.43155958891432
    alpha_x = 0.0
    alpha_y = 5.675938272265563e-16
    gamma_x = 0.030216802240728775
    gamma_y = 0.027448728829722893
    E = 0.005
    s = 0.0
    ref = np.array([beta_x, beta_y, alpha_x, alpha_y, gamma_x, gamma_y, E, s])
    res = np.array([tw_p.beta_x, tw_p.beta_y, tw_p.alpha_x, tw_p.alpha_y, tw_p.gamma_x, tw_p.gamma_y, tw_p.E, tw_p.s])
    assert (ref == res).all()