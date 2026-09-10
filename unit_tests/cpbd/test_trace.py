from copy import deepcopy

import numpy as np
import pytest

from ocelot.cpbd.beam import Particle, Twiss
from ocelot.cpbd.elements import Cavity, Drift, Quadrupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.optics import trace_obj, trace_z


@pytest.mark.parametrize("lengths", [(1.0,), (0.25, 0.5, 0.25)])
def test_trace_z_particle_matches_drift_at_each_requested_position(lengths):
    lattice = MagneticLattice([Drift(l=length) for length in lengths])
    particle = Particle(x=1e-3, px=2e-3, y=-3e-3, py=-1e-3,
                        tau=4e-3, s=3.0, E=1.0, q=1e-12)
    initial = deepcopy(particle)
    z = np.array([0.0, 0.125, 0.25, 0.25, 0.5, 0.75, 0.875, 1.0])

    samples = trace_z(lattice, particle, z)

    np.testing.assert_allclose([p.s for p in samples], initial.s + z, atol=1e-14)
    np.testing.assert_allclose([p.x for p in samples], initial.x + initial.px * z)
    np.testing.assert_allclose([p.y for p in samples], initial.y + initial.py * z)
    for name in ("px", "py", "tau", "p", "E", "q"):
        np.testing.assert_allclose([getattr(p, name) for p in samples], getattr(initial, name))
    assert vars(particle) == vars(initial)
    assert len({id(p) for p in samples}) == len(z)
    assert all(p is not particle for p in samples)


@pytest.mark.parametrize("z", [0.25, 0.75])
def test_trace_z_preserves_input_when_sampling_inside_or_after_first_element(z):
    lattice = MagneticLattice([Drift(l=0.5), Drift(l=0.5)])
    particle = Particle(x=1e-3, px=2e-3, E=1.0)
    initial = deepcopy(particle)

    trace_z(lattice, particle, [z])

    assert vars(particle) == vars(initial)


@pytest.mark.parametrize("n_points", [2, 5, 11])
def test_trace_obj_particle_endpoint_is_independent_of_sampling_density(n_points):
    lattice = MagneticLattice([Drift(l=1.0)])
    particle = Particle(px=1e-3, E=1.0)

    samples = trace_obj(lattice, particle, nPoints=n_points)

    z = np.linspace(0.0, 1.0, n_points)
    np.testing.assert_allclose([p.s for p in samples], z, atol=1e-14)
    np.testing.assert_allclose([p.x for p in samples], 1e-3 * z, atol=1e-14)


def test_trace_z_particle_endpoint_matches_full_map_after_acceleration():
    lattice = MagneticLattice([
        Drift(l=0.25),
        Cavity(l=0.25, v=0.1, phi=0.0),
        Quadrupole(l=0.25, k1=0.7),
        Drift(l=0.25),
    ])
    particle = Particle(x=1e-3, px=2e-3, y=-3e-3, py=-1e-3, E=0.1)
    initial = deepcopy(particle)
    coordinates = ("x", "px", "y", "py", "tau", "p")
    b, r, _ = lattice.transfer_maps(energy=initial.E)
    expected = r @ np.array([getattr(initial, name) for name in coordinates]) + b.ravel()

    samples = trace_z(lattice, particle, np.linspace(0.0, 1.0, 17))

    np.testing.assert_allclose([getattr(samples[-1], name) for name in coordinates],
                               expected, rtol=1e-12, atol=1e-14)
    assert samples[-1].s == pytest.approx(1.0)
    assert samples[-1].E == pytest.approx(0.2)
    assert vars(particle) == vars(initial)


def test_trace_z_twiss_matches_drift_optics_without_modifying_input():
    drift = Drift(l=0.5)
    lattice = MagneticLattice([drift, drift])
    seed = Twiss(beta_x=2.0, beta_y=3.0, alpha_x=0.2, alpha_y=-0.3, E=1.0, s=3.0)
    initial = deepcopy(seed)
    z = np.linspace(0.0, 1.0, 9)

    samples = trace_z(lattice, seed, z)

    np.testing.assert_allclose([t.s for t in samples], initial.s + z)
    for plane in ("x", "y"):
        beta = getattr(initial, "beta_" + plane)
        alpha = getattr(initial, "alpha_" + plane)
        gamma = (1 + alpha**2) / beta
        np.testing.assert_allclose([getattr(t, "beta_" + plane) for t in samples],
                                   beta - 2 * alpha * z + gamma * z**2)
        np.testing.assert_allclose([getattr(t, "alpha_" + plane) for t in samples],
                                   alpha - gamma * z)
    assert vars(seed) == vars(initial)
