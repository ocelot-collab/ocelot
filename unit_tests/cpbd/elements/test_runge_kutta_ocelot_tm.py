import copy

import numpy as np
import pytest

from ocelot.cpbd.beam.particle import ParticleArray
from ocelot.cpbd.elements.bend import Bend
from ocelot.cpbd.elements.drift import Drift
from ocelot.cpbd.elements.quadrupole import Quadrupole
from ocelot.cpbd.elements.sextupole import Sextupole
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaGlobalTM, RungeKuttaOcelotTM, RungeKuttaTM
from ocelot.cpbd.transformations.second_order import SecondTM


def _probe_parray(energy=1.0):
    coords = np.array([
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0],
        [-1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 2.0e-5, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0e-4, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 2.0e-5, 0.0, 0.0],
    ])
    parray = ParticleArray(n=len(coords))
    parray.E = energy
    parray.rparticles[:] = coords.T
    return parray


def _track(element, tm_class, parray):
    tracked = copy.deepcopy(parray)
    element.set_tm(tm_class)
    for tm in element.tms:
        tm.apply(tracked)
    return tracked.rparticles.copy()


def test_runge_kutta_tm_is_legacy_global_tm():
    assert issubclass(RungeKuttaTM, RungeKuttaGlobalTM)


def test_global_runge_kutta_keeps_bend_reference_in_fixed_frame():
    bend = Bend(l=0.5, angle=0.01, e1=0.0, e2=0.0)
    tracked = _track(bend, RungeKuttaGlobalTM, _probe_parray())

    assert tracked[0, 0] != pytest.approx(0.0)
    assert tracked[1, 0] != pytest.approx(0.0)


def test_ocelot_runge_kutta_returns_bend_reference_to_zero():
    bend = Bend(l=0.5, angle=0.01, e1=0.0, e2=0.0)
    tracked = _track(bend, RungeKuttaOcelotTM, _probe_parray())

    assert tracked[:, 0] == pytest.approx(np.zeros(6), abs=1.0e-14)


@pytest.mark.parametrize(
    "element",
    [
        Bend(l=0.5, angle=0.01, e1=0.0, e2=0.0),
        Quadrupole(l=0.5, k1=0.7),
        Sextupole(l=0.2, k2=10.0),
    ],
)
def test_ocelot_runge_kutta_matches_second_order_for_small_amplitudes(element):
    parray = _probe_parray()
    second = _track(copy.deepcopy(element), SecondTM, parray)
    rk = _track(copy.deepcopy(element), RungeKuttaOcelotTM, parray)

    np.testing.assert_allclose(rk[:4], second[:4], rtol=0.0, atol=1.0e-7)


def test_ocelot_runge_kutta_tracks_relative_to_custom_field_reference():
    drift = Drift(l=0.3)
    drift.mag_field = lambda x, y, z: (0.0 * x, 0.02 + 0.0 * x, 0.0 * x)
    parray = _probe_parray()

    global_rk = _track(copy.deepcopy(drift), RungeKuttaGlobalTM, parray)
    ocelot_rk = _track(copy.deepcopy(drift), RungeKuttaOcelotTM, parray)

    assert global_rk[0, 0] != pytest.approx(0.0)
    assert global_rk[1, 0] != pytest.approx(0.0)
    assert ocelot_rk[:, 0] == pytest.approx(np.zeros(6), abs=1.0e-14)
