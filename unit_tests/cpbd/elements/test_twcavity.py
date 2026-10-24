"""
Physics tests for the travelling-wave cavity.

The travelling-wave structure differs from the standing-wave one only in its
transverse block: the body has no ponderomotive focusing and all the radial
focusing sits in the edge lenses. The longitudinal map is shared through
``ocelot.cpbd.rf_maps``, so most of these tests pin the travelling-wave result
either to a closed-form reference or to the standing-wave cavity.
"""

import numpy as np
import pytest

from ocelot.common.globals import m_e_GeV, speed_of_light
from ocelot.cpbd.beam import ParticleArray
from ocelot.cpbd.elements import Cavity, Drift, TWCavity
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.navi import Navigator
from ocelot.cpbd.physics_proc import PhysProc
from ocelot.cpbd.track import track

E0 = 0.150          # GeV
L = 1.0             # m
V = 0.02            # GV
FREQ = 1.3e9        # Hz


def total_R(element, energy):
    """Product of the ENTRANCE -> MAIN -> EXIT matrices of one element."""
    total = np.eye(6)
    for R in element.R(energy):
        total = R @ total
    return total


def reference_R(z, v, energy, freq, phi_deg):
    """
    Closed-form travelling-wave matrix, written out independently of the
    implementation (Rosenzweig & Serafini, PRE 49, 1599 (1994), plus the
    standard RF longitudinal block).
    """
    phi = np.radians(phi_deg)
    de = v * np.cos(phi)
    gradient = de / z
    e1 = energy + de
    g0, g1 = energy / m_e_GeV, e1 / m_e_GeV
    b0, b1 = np.sqrt(1. - 1. / g0 ** 2), np.sqrt(1. - 1. / g1 ** 2)
    k = 2. * np.pi * freq / speed_of_light

    p0, p1 = np.sqrt(g0 ** 2 - 1.), np.sqrt(g1 ** 2 - 1.)
    body = np.eye(6)
    # at the zero crossing de is 1e-18 rather than 0, so g1 - g0 rounds to
    # exactly zero; the drift limit of R12 is the length
    body[0, 1] = body[2, 3] = (z * p0 * (np.arcsinh(p1) - np.arcsinh(p0)) / (g1 - g0)
                               if g1 != g0 else z)
    body[1, 1] = body[3, 3] = p0 / p1
    body[4, 5] = -z / (g1 * g1 * g0 * b1) * (g1 + g0) / (b1 + b0)
    # tau is c*delta_t, so the RF phase is phi-k*tau.
    if abs(g1 - g0) < 1e-8 * abs(g0) and abs(np.cos(phi)) < 1e-3:
        body[4, 4] = 1. - k * z * v * np.sin(phi) / (2. * m_e_GeV * g0 ** 3 * b0 ** 3)
    else:
        body[4, 4] = 1. + (k * z * v / m_e_GeV * np.sin(phi)
                           * (g0 * g1 * (b0 * b1 - 1.) + 1.) / (b1 * g1 * (g0 - g1) ** 2))
    body[5, 4] = k * np.sin(phi) * v / (g1 * b1 * m_e_GeV)
    body[5, 5] = g0 / g1 * b0 / b1

    entrance = np.eye(6)
    entrance[1, 0] = entrance[3, 2] = -gradient / (2. * energy)
    exit_ = np.eye(6)
    exit_[1, 0] = exit_[3, 2] = +gradient / (2. * e1)
    return exit_ @ body @ entrance


def jacobian(element, energy, eps=1e-9):
    """Numerical Jacobian of the active tracking map of a single element."""
    lattice = MagneticLattice((element,))
    num = np.zeros((6, 6))
    for j in range(6):
        parray = ParticleArray(n=1)
        parray.E = energy
        parray.rparticles[:] = 0.
        parray.rparticles[j, 0] = eps
        _, parray = track(lattice, parray, navi=Navigator(lattice), print_progress=False)
        num[:, j] = parray.rparticles[:, 0] / eps
    return num


def test_twcavity_constructor_is_silent(capsys):
    """The atom used to announce 'NOT FINISHED. IT WILL NOT WORK' on stdout."""
    TWCavity(l=L, v=V, phi=0., freq=FREQ, eid="TW")
    assert capsys.readouterr().out == ""


@pytest.mark.parametrize("phi", [0., 20., 90., 160., 180., 270.])
def test_twcavity_matches_closed_form(phi):
    cav = TWCavity(l=L, v=V, phi=phi, freq=FREQ, eid="TW")
    assert total_R(cav, E0) == pytest.approx(reference_R(L, V, E0, FREQ, phi), abs=1e-12)


@pytest.mark.parametrize("phi", [0., 35., 180.])
def test_twcavity_exit_edge_uses_the_exit_energy(phi):
    """
    The exit edge lens is +gamma'/(2 gamma1). The transfer-map chain hands the
    edge the energy *behind* the body, so the energy gain must not be added a
    second time inside the edge map.
    """
    cav = TWCavity(l=L, v=V, phi=phi, freq=FREQ, eid="TW")
    _, _, r_exit = cav.R(E0)
    gradient = V * np.cos(np.radians(phi)) / L
    assert r_exit[1, 0] == pytest.approx(gradient / (2. * (E0 + V * np.cos(np.radians(phi)))), rel=1e-12)
    assert r_exit[3, 2] == pytest.approx(r_exit[1, 0], rel=1e-12)


def test_twcavity_r12_is_the_length_at_the_zero_crossing():
    """
    cos(90 deg) is 6e-17 rather than 0, so a naive log(1 + de/E) underflows to
    zero and collapses R12. The body must still transport like a drift.
    """
    cav = TWCavity(l=L, v=V, phi=90., freq=FREQ, eid="TW")
    _, r_main, _ = cav.R(E0)
    assert r_main[0, 1] == pytest.approx(L, rel=1e-12)
    assert r_main[2, 3] == pytest.approx(L, rel=1e-12)


def test_twcavity_converges_to_a_drift_as_the_voltage_vanishes():
    """Deviation from a drift must be first order in the voltage."""
    drift = total_R(Drift(l=L), E0)
    block = np.ix_([0, 1, 2, 3], [0, 1, 2, 3])
    errors = []
    for v in (1e-5, 1e-6, 1e-7):
        M = total_R(TWCavity(l=L, v=v, phi=10., freq=FREQ), E0)
        errors.append(np.max(np.abs(M[block] - drift[block])))
    assert errors[0] / errors[1] == pytest.approx(10., rel=1e-3)
    assert errors[1] / errors[2] == pytest.approx(10., rel=1e-3)
    assert errors[-1] < 1e-6


@pytest.mark.parametrize("phi", [0., 25., 90., 200., 270.])
def test_twcavity_longitudinal_block_matches_standing_wave(phi):
    """The two families share ocelot.cpbd.rf_maps, so this must be exact."""
    kwargs = dict(l=L, v=V, phi=phi, freq=FREQ)
    tw = total_R(TWCavity(**kwargs), E0)
    sw = total_R(Cavity(**kwargs), E0)
    assert tw[4:6, 4:6] == pytest.approx(sw[4:6, 4:6], abs=0.)


@pytest.mark.parametrize("phi", [0., 25., 90., 180., 270.])
def test_twcavity_tracking_linearises_to_the_r_matrix(phi):
    kwargs = dict(l=L, v=V, phi=phi, freq=FREQ)
    R = total_R(TWCavity(**kwargs), E0)
    num = jacobian(TWCavity(**kwargs), E0)

    assert num == pytest.approx(R, abs=1e-6)


@pytest.mark.parametrize("phi,expected", [(0., V), (60., 0.5 * V), (90., 0.), (180., -V)])
def test_twcavity_energy_gain(phi, expected):
    cav = TWCavity(l=L, v=V, phi=phi, freq=FREQ)
    lattice = MagneticLattice((Drift(l=0.1), cav, Drift(l=0.1)))
    parray = ParticleArray(n=8)
    parray.E = E0
    parray.rparticles[:] = np.random.default_rng(3).normal(0., 1., (6, 8)) * 1e-4
    _, parray = track(lattice, parray, navi=Navigator(lattice), print_progress=False)
    assert parray.E - E0 == pytest.approx(expected, abs=1e-12)
    assert np.all(np.isfinite(parray.rparticles))


@pytest.mark.parametrize("phi", [0., 25., 90., 180., 270.])
def test_twcavity_longitudinal_tracking_matches_standing_wave(phi):
    """
    The second-order longitudinal map is shared, and neither family's edge maps
    touch tau or delta, so the longitudinal coordinates must come out equal.
    """
    kwargs = dict(l=L, v=V, phi=phi, freq=FREQ)
    start = np.random.default_rng(11).normal(0., 1., (6, 32)) * np.array(
        [[1e-4], [1e-5], [1e-4], [1e-5], [1e-3], [1e-3]])
    out = {}
    for cls in (TWCavity, Cavity):
        lattice = MagneticLattice((cls(**kwargs),))
        parray = ParticleArray(n=32)
        parray.E = E0
        parray.rparticles[:] = start
        _, parray = track(lattice, parray, navi=Navigator(lattice), print_progress=False)
        out[cls] = parray.rparticles.copy()
    assert out[TWCavity][4:6] == pytest.approx(out[Cavity][4:6], abs=0.)


def test_twcavity_slicing_preserves_the_full_map():
    """At relativistic energy, real slices agree within the RF approximation."""
    class CountSteps(PhysProc):
        def __init__(self):
            super().__init__()
            self.steps = []

        def apply(self, parray, dz):
            self.steps.append(dz)

    kwargs = dict(l=L, v=V, phi=25., freq=FREQ)
    start = np.random.default_rng(5).normal(0., 1., (6, 16)) * 1e-4
    results = []
    for unit_step, expected_steps in ((L, 1), (0.1, 10)):
        cavity = TWCavity(**kwargs)
        end = Drift(l=0.1)
        lattice = MagneticLattice((Drift(l=0.1), cavity, end))
        parray = ParticleArray(n=16)
        parray.E = E0
        parray.rparticles[:] = start
        navi = Navigator(lattice)
        navi.unit_step = unit_step
        process = CountSteps()
        navi.add_physics_proc(process, cavity, end)
        _, parray = track(lattice, parray, navi=navi, print_progress=False)
        assert len(process.steps) == expected_steps
        assert sum(process.steps) == pytest.approx(L)
        results.append((parray.E, parray.rparticles.copy()))
    assert results[0][0] == pytest.approx(results[1][0], rel=1e-12)
    assert results[0][1] == pytest.approx(results[1][1], abs=5e-8)


@pytest.mark.parametrize("energy", [0.0007, 0.001, 0.002, 0.005, 0.010, 0.050, 0.150, 1.000])
def test_twcavity_body_matches_elegant_default_rf_matrix(energy):
    """
    elegant's default RFCA matrix is the travelling-wave one (STANDING_WAVE
    defaults to 0, BODY_FOCUS_MODEL to None): the "old TRANSPORT matrix" block
    of compute_matrices.c. Transcribed here from its own variables, where
    sf1 = gamma0, f2 = gamma1, sf5 = p1, f4 = p0 + gamma0, f6 = gamma1 + p1
    and cscPhi0 / dgMax = 1 / dgamma.

    The body carries no focusing there either (R[0][0] = 1, R[1][0] unset), so
    the two codes should agree exactly, including at low energy where the
    beta -> 1 form of R12 and R22 is off by tens of percent.
    """
    phi = 0.
    g0 = energy / m_e_GeV
    p00 = np.sqrt(g0 ** 2 - 1.)
    dgamma = V * np.cos(np.radians(phi)) / m_e_GeV
    sf1 = np.sqrt(1. + p00 ** 2)
    f2 = sf1 + dgamma
    sf5 = np.sqrt(f2 ** 2 - 1.)
    f4 = p00 + sf1
    f6 = f2 + sf5
    r12_elegant = L * p00 * (np.log(f6) - np.log(f4)) / dgamma
    r22_elegant = p00 / sf5

    _, r_main, _ = TWCavity(l=L, v=V, phi=phi, freq=FREQ).R(energy)
    assert r_main[0, 0] == 1.
    assert r_main[1, 0] == 0.
    assert r_main[0, 1] == pytest.approx(r12_elegant, rel=1e-12)
    assert r_main[1, 1] == pytest.approx(r22_elegant, rel=1e-12)
    assert r_main[2, 3] == pytest.approx(r_main[0, 1], rel=1e-12)
    assert r_main[3, 3] == pytest.approx(r_main[1, 1], rel=1e-12)
