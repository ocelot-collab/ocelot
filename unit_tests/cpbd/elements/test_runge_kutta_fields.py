import numpy as np
import pytest

from ocelot.common.globals import m_e_GeV, speed_of_light
from ocelot.cpbd.elements.bend import Bend
from ocelot.cpbd.elements.drift import Drift
from ocelot.cpbd.elements.octupole import Octupole
from ocelot.cpbd.elements.quadrupole import Quadrupole
from ocelot.cpbd.elements.sextupole import Sextupole
from ocelot.cpbd.elements.solenoid import Solenoid
from ocelot.cpbd.elements.undulator import Undulator
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaTM


def _brho(energy):
    gamma = energy / m_e_GeV
    beta = np.sqrt(1. - 1. / (gamma * gamma))
    return beta * energy * 1e9 / speed_of_light


def _field_tuple(field, x=1.e-3, y=2.e-3, z=3.e-3):
    return field(np.array([x]), np.array([y]), np.array([z]))


def test_undulator_runge_kutta_uses_user_mag_field():
    def custom_field(x, y, z):
        return x + 1., y + 2., z + 3.

    und = Undulator(lperiod=0.03, nperiods=2, Kx=1.0, tm=RungeKuttaTM)
    und.mag_field = custom_field

    bx, by, bz = _field_tuple(und.tms[0].get_params(1.0).mag_field)

    assert bx == pytest.approx(np.array([1.001]))
    assert by == pytest.approx(np.array([2.002]))
    assert bz == pytest.approx(np.array([3.003]))


def test_drift_and_bend_can_use_user_mag_field():
    def custom_field(x, y, z):
        return 0.1 + x, 0.2 + y, 0.3 + z

    drift = Drift(l=0.2, tm=RungeKuttaTM)
    drift.mag_field = custom_field
    bend = Bend(l=0.5, angle=0.01, tm=RungeKuttaTM)
    bend.mag_field = custom_field

    for field in (
        drift.tms[0].get_params(1.0).mag_field,
        bend.element.create_runge_kutta_main_params(1.0).mag_field,
    ):
        bx, by, bz = _field_tuple(field)
        assert bx == pytest.approx(np.array([0.101]))
        assert by == pytest.approx(np.array([0.202]))
        assert bz == pytest.approx(np.array([0.303]))


def test_bend_default_runge_kutta_field_is_uniform_hard_edge_dipole():
    energy = 1.25
    bend = Bend(l=2.0, angle=0.1)
    field = bend.element.create_runge_kutta_main_params(energy).mag_field

    bx, by, bz = _field_tuple(field)

    assert bx == pytest.approx(np.array([0.]))
    assert by == pytest.approx(np.array([_brho(energy) * bend.angle / bend.l]))
    assert bz == pytest.approx(np.array([0.]))


def test_tilted_bend_default_runge_kutta_field_matches_csr_geometry():
    energy = 1.25
    bend = Bend(l=2.0, angle=0.1, tilt=np.pi / 2.)
    field = bend.element.create_runge_kutta_main_params(energy).mag_field

    bx, by, bz = _field_tuple(field)

    assert bx == pytest.approx(np.array([-_brho(energy) * bend.angle / bend.l]))
    assert by == pytest.approx(np.array([0.]), abs=1.e-16)
    assert bz == pytest.approx(np.array([0.]))


def test_quadrupole_default_runge_kutta_field_is_hard_edge_quadrupole():
    energy = 1.25
    quad = Quadrupole(l=0.4, k1=2.0, tm=RungeKuttaTM)
    field = quad.tms[0].get_params(energy).mag_field

    bx, by, bz = _field_tuple(field)

    assert bx == pytest.approx(np.array([_brho(energy) * quad.k1 * 2.e-3]))
    assert by == pytest.approx(np.array([_brho(energy) * quad.k1 * 1.e-3]))
    assert bz == pytest.approx(np.array([0.]))


def test_sextupole_default_runge_kutta_field_is_hard_edge_sextupole():
    energy = 1.25
    sext = Sextupole(l=0.3, k2=3.0, tm=RungeKuttaTM)
    field = sext.tms[0].get_params(energy).mag_field

    bx, by, bz = _field_tuple(field)

    assert bx == pytest.approx(np.array([_brho(energy) * sext.k2 * 1.e-3 * 2.e-3]))
    assert by == pytest.approx(np.array([_brho(energy) * 0.5 * sext.k2 * (1.e-6 - 4.e-6)]))
    assert bz == pytest.approx(np.array([0.]))


def test_octupole_default_runge_kutta_field_is_hard_edge_octupole():
    energy = 1.25
    octupole = Octupole(l=0.3, k3=4.0, tm=RungeKuttaTM)
    field = octupole.tms[0].get_params(energy).mag_field

    bx, by, bz = _field_tuple(field)

    assert bx == pytest.approx(np.array([_brho(energy) * octupole.k3 * (3. * 1.e-6 * 2.e-3 - 8.e-9) / 6.]))
    assert by == pytest.approx(np.array([_brho(energy) * octupole.k3 * (1.e-9 - 3. * 1.e-3 * 4.e-6) / 6.]))
    assert bz == pytest.approx(np.array([0.]))


def test_solenoid_default_runge_kutta_field_is_hard_edge_longitudinal_field():
    energy = 1.25
    solenoid = Solenoid(l=0.3, k=0.7, tm=RungeKuttaTM)
    field = solenoid.tms[0].get_params(energy).mag_field

    bx, by, bz = _field_tuple(field)

    assert bx == pytest.approx(np.array([0.]))
    assert by == pytest.approx(np.array([0.]))
    assert bz == pytest.approx(np.array([2. * solenoid.k * _brho(energy)]))
