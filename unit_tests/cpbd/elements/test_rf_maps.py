"""Independent time-of-flight and energy checks for both RF cavity models."""

import numpy as np
import pytest
from scipy.integrate import quad

from ocelot.common.globals import m_e_GeV, speed_of_light
from ocelot.cpbd.beam import ParticleArray
from ocelot.cpbd.elements import Cavity, TWCavity
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.rf_maps import rf_long_r_terms, rf_long_t_terms
from ocelot.cpbd.track import track


@pytest.mark.parametrize("energy", [0.0007, 0.002, 0.150, 2.4])
@pytest.mark.parametrize("voltage_fraction", [0., 1e-14, 0.1, -0.1])
@pytest.mark.parametrize("phase", [0., 25., 89.999999, 90., 180., 270., 270.000001])
def test_longitudinal_coefficients_from_flight_time_integrals(energy, voltage_fraction, phase):
    """Integrate derivatives of 1/beta, independently of the endpoint formulas."""
    length, freq = 0.73, 1.3e9
    voltage = energy * voltage_fraction
    phi = np.deg2rad(phase)
    k = 2. * np.pi * freq / speed_of_light
    g0, v = energy / m_e_GeV, voltage / m_e_GeV
    p0 = np.sqrt(g0 ** 2 - 1.)
    p1 = np.sqrt((g0 + v * np.cos(phi)) ** 2 - 1.)
    a = v * k * np.sin(phi)

    def integral(order, power):
        def integrand(u):
            g = g0 + u * v * np.cos(phi)
            derivative = -(g * g - 1.) ** -1.5 if order == 1 else 3. * g * (g * g - 1.) ** -2.5
            return u ** power * derivative
        return quad(integrand, 0., 1., epsabs=1e-24, epsrel=1e-11)[0]

    expected_r = (length * a * integral(1, 1), length * p0 * integral(1, 0), a / p1, p0 / p1)
    expected_t = (
        0.5 * length * p0 ** 2 * integral(2, 0),
        length * p0 * a * integral(2, 1),
        0.5 * length * (a ** 2 * integral(2, 2) - v * k ** 2 * np.cos(phi) * integral(1, 1)),
    )
    np.testing.assert_allclose(rf_long_r_terms(length, voltage, energy, freq, phi), expected_r, rtol=2e-12, atol=1e-25)
    np.testing.assert_allclose(rf_long_t_terms(length, voltage, energy, freq, phi), expected_t, rtol=2e-12, atol=1e-25)


@pytest.mark.parametrize("element_type", [Cavity, TWCavity])
@pytest.mark.parametrize("energy", [0.0007, 0.002, 0.150])
@pytest.mark.parametrize("phase", [0., 25., 90., 180., 270.])
def test_cavity_tracking_jacobian_and_rf_energy(element_type, energy, phase):
    """Check every matrix entry and energy gain in the time-like coordinate."""
    length, voltage, freq = 0.73, 0.1 * energy, 1.3e9
    cavity = element_type(l=length, v=voltage, phi=phase, freq=freq)
    lattice = MagneticLattice((cavity,))
    matrix = np.eye(6)
    for r in cavity.R(energy):
        matrix = r @ matrix
    step = 1e-7
    initial = np.column_stack((step * np.eye(6), -step * np.eye(6)))
    particles = ParticleArray(n=12)
    particles.E = energy
    particles.rparticles[:] = initial
    _, particles = track(lattice, particles, print_progress=False, calc_tws=False)
    jacobian = (particles.rparticles[:, :6] - particles.rparticles[:, 6:]) / (2. * step)
    np.testing.assert_allclose(jacobian, matrix, rtol=1e-9, atol=1e-9)

    # An RF oscillator advances by omega*delta_t, independent of beam beta.
    tau = np.array([-1e-3, 2e-3, 3e-3])
    delta = np.array([2e-3, -1e-3, 0.])
    particles = ParticleArray(n=3)
    particles.E = energy
    particles.rparticles[4] = tau
    particles.rparticles[5] = delta
    pc0 = np.sqrt(energy ** 2 - m_e_GeV ** 2)
    expected_energy = energy + pc0 * delta + voltage * np.cos(
        np.deg2rad(phase) - 2. * np.pi * freq * tau / speed_of_light)
    _, particles = track(lattice, particles, print_progress=False, calc_tws=False)
    pc1 = np.sqrt(particles.E ** 2 - m_e_GeV ** 2)
    np.testing.assert_allclose(particles.E + pc1 * particles.rparticles[5], expected_energy, rtol=2e-15)


@pytest.mark.parametrize("energy", [0.0007, 0.002, 0.150])
def test_both_zero_crossings_have_the_correct_parity(energy):
    args = dict(z=0.73, V=0.1 * energy, E=energy, freq=1.3e9)
    r_plus = np.array(rf_long_r_terms(**args, phi=np.pi / 2.))
    r_minus = np.array(rf_long_r_terms(**args, phi=-np.pi / 2.))
    t_plus = np.array(rf_long_t_terms(**args, phi=np.pi / 2.))
    t_minus = np.array(rf_long_t_terms(**args, phi=-np.pi / 2.))
    np.testing.assert_allclose(r_minus, r_plus * [-1., 1., -1., 1.], rtol=1e-14)
    np.testing.assert_allclose(t_minus, t_plus * [1., -1., 1.], rtol=1e-14)
