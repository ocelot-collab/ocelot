"""XYQuadrupole against the revised-v4 paper's Hamilton equations.

The numerical oracle evolves both the reference and the fundamental matrix;
it does not reuse the element's analytical transport or arrival functions.
"""
from copy import deepcopy

import numpy as np
import pytest
from scipy.integrate import quad, solve_ivp

from ocelot.common.globals import m_e_GeV
from ocelot.cpbd.beam import Particle, ParticleArray, Twiss
from ocelot.cpbd.elements import Drift, Quadrupole, XYQuadrupole
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.navi import Navigator
from ocelot.cpbd.optics import lattice_transfer_map, trace_z, twiss
from ocelot.cpbd.r_matrix import rot_mtx
from ocelot.cpbd.track import tracking_step
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


def _integrated_map(k1, a, b, length, energy, start=0.0):
    """Integrate Eq. (29), with tau=-sigma, along the paraxial reference."""
    igamma2 = (m_e_GeV / energy)**2 if energy else 0.0
    beta2 = 1.0 - igamma2
    beta = np.sqrt(beta2)

    def rhs(s, state):
        x, px, y, py = state[:4]
        generator = np.zeros((6, 6))
        generator[0, 1] = generator[2, 3] = 1.0
        generator[1, 0], generator[3, 2] = -k1, k1
        generator[1, 5] = generator[4, 0] = k1 * x / beta
        generator[3, 5] = generator[4, 2] = -k1 * y / beta
        generator[4, 5] = -igamma2 / beta2
        return np.r_[px, -k1 * x, py, k1 * y,
                     (generator @ state[4:].reshape(6, 6)).ravel()]

    initial = np.r_[a, 0.0, b, 0.0, np.eye(6).ravel()]
    if start:
        before = solve_ivp(rhs, (0, start), initial, method="DOP853", rtol=2e-13, atol=2e-15)
        assert before.success
        initial[:4] = before.y[:4, -1]
    solution = solve_ivp(rhs, (start, start + length), initial,
                         method="DOP853", rtol=2e-13, atol=2e-15)
    assert solution.success
    return solution.y[4:, -1].reshape(6, 6)


@pytest.mark.parametrize("k1", [-2.0, -0.1, 0.0, 0.1, 2.0])
@pytest.mark.parametrize("energy", [0.0, 0.001, 0.13])
def test_map_matches_independent_hamilton_equations(k1, energy):
    element = XYQuadrupole(l=1.2, k1=k1, x_offs=0.02, y_offs=-0.03)
    expected = _integrated_map(k1, element.x_offs, element.y_offs, element.l, energy)
    np.testing.assert_allclose(element.R(energy)[0], expected, rtol=3e-12, atol=3e-14)


@pytest.mark.parametrize("k1", [-1.3, 0.0, 1.3])
@pytest.mark.parametrize("tilt", [0.0, 0.37])
def test_map_is_symplectic_in_ocelot_coordinates(k1, tilt):
    element = XYQuadrupole(l=0.9, k1=k1, x_offs=0.02, y_offs=0.03, tilt=tilt)
    r = element.R(0.002)[0]
    j = np.kron(np.eye(3), [[0.0, 1.0], [-1.0, 0.0]])
    j[4:, 4:] *= -1  # Ocelot uses tau=-sigma and the same energy deviation.
    np.testing.assert_allclose(r.T @ j @ r, j, atol=2e-15)


@pytest.mark.parametrize("k1", [-2.0, 0.0, 2.0])
@pytest.mark.parametrize("energy", [0.0, 0.001, 1.0])
def test_zero_offsets_reduce_to_ordinary_quadrupole(k1, energy):
    xy = XYQuadrupole(l=0.7, k1=k1, tilt=0.23)
    ordinary = Quadrupole(l=xy.l, k1=k1, tilt=xy.tilt)
    np.testing.assert_allclose(xy.R(energy), ordinary.R(energy), atol=2e-15)
    np.testing.assert_array_equal(xy.B(energy), np.zeros((1, 6, 1)))


def test_zero_gradient_is_drift_and_zero_length_is_identity():
    xy = XYQuadrupole(l=0.7, x_offs=0.02, y_offs=0.03)
    np.testing.assert_allclose(xy.R(0.001), Drift(l=xy.l).R(0.001), atol=1e-15)
    xy.l, xy.k1 = 0.0, 3.0
    np.testing.assert_array_equal(xy.R(0.001)[0], np.eye(6))
    np.testing.assert_array_equal(xy.get_section_tms(0)[0].get_params(0.001).R, np.eye(6))


@pytest.mark.parametrize("a,b,absent", [(0.02, 0.0, [2, 3]), (0.0, -0.03, [0, 1])])
def test_one_plane_offsets_generate_dispersion_only_in_that_plane(a, b, absent):
    xy = XYQuadrupole(l=0.8, k1=1.2, x_offs=a, y_offs=b)
    r = xy.R(0)[0]
    np.testing.assert_array_equal(r[absent, 5], 0.0)
    np.testing.assert_array_equal(r[4, absent], 0.0)
    assert np.linalg.norm(r[:4, 5]) > 0
    assert r[4, 5] > 0  # Geometric timing term has Ocelot's tau sign.
    reversed_offsets = XYQuadrupole(l=xy.l, k1=xy.k1, x_offs=-a, y_offs=-b).R(0)[0]
    np.testing.assert_array_equal(reversed_offsets[:4, 5], -r[:4, 5])
    assert reversed_offsets[4, 5] == r[4, 5]


@pytest.mark.parametrize("k1", [-1e-14, 1e-14, -0.099999, 0.099999, -0.100001, 0.100001, -1.0, 1.0])
def test_arrival_term_matches_quadrature_including_weak_gradient(k1):
    xy = XYQuadrupole(l=1.0, k1=k1, x_offs=0.02, y_offs=-0.03)

    def integrand(s, k):
        w = np.sqrt(abs(k))
        if k > 0:
            return s * np.cos(w*s) * np.sin(w*s) / w
        return s * np.cosh(w*s) * np.sinh(w*s) / w

    ix = quad(integrand, 0, xy.l, args=(k1,), epsabs=1e-14)[0]
    iy = quad(integrand, 0, xy.l, args=(-k1,), epsabs=1e-14)[0]
    expected = k1**2 / 2 * (xy.x_offs**2 * ix + xy.y_offs**2 * iy)
    assert xy.R(0)[0][4, 5] == pytest.approx(expected, rel=2e-13, abs=0)


@pytest.mark.parametrize("k1", [-1.2, 1.2])
@pytest.mark.parametrize("first_order_only", [False, True])
def test_interior_slice_matches_hamilton_equations(k1, first_order_only):
    xy = XYQuadrupole(l=1.3, k1=k1, x_offs=0.02, y_offs=-0.03)
    tm = xy.get_section_tms(0.4, start_l=0.6, first_order_only=first_order_only)[0]
    expected = _integrated_map(k1, xy.x_offs, xy.y_offs, 0.4, 0.002, start=0.6)
    np.testing.assert_allclose(tm.get_params(0.002).get_rotated_R(), expected, rtol=3e-12, atol=3e-14)
    assert tm.delta_length == 0.4
    assert tm.get_delta_e() == 0.0


@pytest.mark.parametrize("tilt", [0.0, 0.41])
@pytest.mark.parametrize("aligned", [False, True])
def test_slice_composition_preserves_full_affine_map(tilt, aligned):
    xy = XYQuadrupole(l=1.1, k1=1.3, x_offs=0.02, y_offs=-0.03, tilt=tilt,
                      dx=1e-3 if aligned else 0.0, dy=-2e-3 if aligned else 0.0)
    r, b = np.eye(6), np.zeros((6, 1))
    positions = [0.0, 0.07, 0.31, 0.72, xy.l]
    for start, end in zip(positions[:-1], positions[1:]):
        params = xy.get_section_tms(end-start, start_l=start)[0].get_params(0.001)
        section = params.get_rotated_R()
        r, b = section @ r, section @ b + params.B
    np.testing.assert_allclose(r, xy.R(0.001)[0], atol=2e-15)
    np.testing.assert_allclose(b, xy.B(0.001)[0], atol=2e-17)


def test_tilt_and_alignment_translate_the_rolled_map():
    xy = XYQuadrupole(l=0.8, k1=-1.2, x_offs=0.02, y_offs=-0.03, tilt=0.31, dx=1e-3, dy=-2e-3)
    unrolled = XYQuadrupole(l=xy.l, k1=xy.k1, x_offs=xy.x_offs, y_offs=xy.y_offs)
    expected_r = rot_mtx(-xy.tilt) @ unrolled.R(0.13)[0] @ rot_mtx(xy.tilt)
    offset = np.array([xy.dx, 0, xy.dy, 0, 0, 0])[:, None]
    np.testing.assert_allclose(xy.R(0.13)[0], expected_r, atol=1e-15)
    np.testing.assert_allclose(xy.B(0.13)[0], (np.eye(6) - expected_r) @ offset, atol=1e-17)


@pytest.mark.parametrize("steps", [1, 7, 50])
def test_navigator_particle_array_tracking_is_independent_of_step_size(steps):
    xy = XYQuadrupole(l=1.1, k1=-1.2, x_offs=0.02, y_offs=-0.03, tilt=0.27)
    lattice = MagneticLattice([xy], method={"global": SecondTM})
    particles = ParticleArray(n=8)
    particles.E = 0.13
    particles.rparticles[:] = np.random.default_rng(12).normal(size=(6, 8)) * 1e-4
    particles.rparticles[:, 0] = 0  # The reference remains the zero particle.
    expected = xy.R(particles.E)[0] @ particles.rparticles
    navigator = Navigator(lattice)
    for _ in range(steps):
        tracking_step(lattice, particles, xy.l / steps, navi=navigator)
    np.testing.assert_allclose(particles.rparticles, expected, rtol=2e-12, atol=2e-17)
    assert particles.s == pytest.approx(xy.l)
    assert particles.E == 0.13
    assert type(xy.tms[0]) is TransferMap


@pytest.mark.parametrize("steps", [30, 100])
def test_repeated_occurrences_reset_reference_and_tolerate_boundary_roundoff(steps):
    xy = XYQuadrupole(l=0.3, k1=1.2, x_offs=0.02, y_offs=0.03)
    lattice = MagneticLattice([Drift(l=0.1), xy, Drift(l=0.2), xy, Drift(l=0.1)])
    initial = np.array([1e-4, -2e-4, 3e-4, 1e-4, 2e-4, 1e-3])
    particles = ParticleArray(n=1)
    particles.E = 0.13
    particles.rparticles[:, 0] = initial
    expected = lattice_transfer_map(lattice, particles.E) @ initial
    navigator = Navigator(lattice)
    for _ in range(steps):
        tracking_step(lattice, particles, lattice.totalLen / steps, navi=navigator)
    np.testing.assert_allclose(particles.rparticles[:, 0], expected, rtol=2e-12, atol=2e-17)
    assert particles.s == pytest.approx(lattice.totalLen)


def test_particle_trace_and_twiss_use_the_same_reference_at_each_position():
    xy = XYQuadrupole(l=1.1, k1=1.2, x_offs=0.02, y_offs=-0.03)
    lattice = MagneticLattice([xy])
    p = Particle(x=1e-4, px=-2e-4, y=3e-4, py=1e-4, tau=2e-4, p=1e-3, E=0.13)
    initial = deepcopy(p)
    positions = np.linspace(0.0, xy.l, 11)
    names = ("x", "px", "y", "py", "tau", "p")
    seed = np.array([getattr(p, name) for name in names])
    for s, sample in zip(positions, trace_z(lattice, p, positions)):
        expected = np.eye(6) if s == 0 else _integrated_map(xy.k1, xy.x_offs, xy.y_offs, s, p.E)
        np.testing.assert_allclose([getattr(sample, name) for name in names], expected @ seed, atol=2e-16)
    assert vars(p) == vars(initial)
    tws = Twiss(beta_x=3.0, beta_y=4.0, E=p.E)
    end = twiss(lattice, tws, nPoints=11)[-1]
    r = lattice_transfer_map(lattice, p.E)
    np.testing.assert_allclose([end.Dx, end.Dxp, end.Dy, end.Dyp], r[:4, 5], atol=2e-15)


@pytest.mark.parametrize("name,value", [("k1", -0.7), ("x_offs", 0.04), ("y_offs", 0.05), ("l", 0.6), ("tilt", 0.2)])
def test_parameter_changes_invalidate_maps(name, value):
    xy = XYQuadrupole(l=0.9, k1=1.0, x_offs=0.02, y_offs=-0.03)
    before = xy.R(0.13)[0].copy()
    active, optics = xy.tms, xy.first_order_tms
    setattr(xy, name, value)
    assert xy.tms is not active
    assert xy.first_order_tms is not optics
    assert not np.allclose(xy.R(0.13)[0], before)


def test_map_cache_depends_on_energy_and_second_order_is_not_fabricated():
    xy = XYQuadrupole(l=0.9, k1=1.0, x_offs=0.02, y_offs=-0.03)
    ultra = xy.R(0)[0]
    finite = xy.R(0.001)[0]
    assert not np.allclose(ultra, finite)
    np.testing.assert_array_equal(xy.R(0)[0], ultra)
    with pytest.raises(RuntimeError, match="does not declare support"):
        xy.set_tm(SecondTM)
    with pytest.raises(NotImplementedError, match="cubic Hamiltonian"):
        xy.element.create_second_order_main_params(0.13)


@pytest.mark.parametrize("k1", [-1.2, 1.2])
def test_planar_survey_matches_integrated_varying_curvature(k1):
    xy = XYQuadrupole(l=1.1, k1=k1, x_offs=0.02)
    end, frame_end, mid, frame_mid = xy.get_transfer_geometry()

    def angle(s):
        w = np.sqrt(abs(k1))
        sine = np.sin(w*s) if k1 > 0 else np.sinh(w*s)
        return k1 * xy.x_offs * sine / w

    for length, position, frame in [(xy.l, end, frame_end), (xy.l/2, mid, frame_mid)]:
        expected_position = [quad(lambda s: -np.sin(angle(s)), 0, length)[0], 0,
                             quad(lambda s: np.cos(angle(s)), 0, length)[0]]
        c, s = np.cos(angle(length)), np.sin(angle(length))
        expected_frame = [[c, 0, -s], [0, 1, 0], [s, 0, c]]
        np.testing.assert_allclose(position, expected_position, atol=2e-13)
        np.testing.assert_allclose(frame, expected_frame, atol=2e-13)


def test_survey_preserves_bishop_frame_and_rolls_the_reference():
    from ocelot.common.geometry import get_tilt_matrix

    xy = XYQuadrupole(l=1.3, k1=1.2, x_offs=0.02, y_offs=-0.03)
    end, frame_end, mid, frame_mid = xy.get_transfer_geometry()
    for frame in (frame_end, frame_mid):
        np.testing.assert_allclose(frame.T @ frame, np.eye(3), atol=2e-12)
        assert np.linalg.det(frame) == pytest.approx(1.0, abs=2e-12)
    assert end[0] < 0 and end[1] < 0 and end[2] < xy.l
    xy.tilt = 0.37
    roll = get_tilt_matrix(xy.tilt)
    rotated_end, rotated_frame_end, rotated_mid, rotated_frame_mid = xy.get_transfer_geometry()
    np.testing.assert_allclose(rotated_end, roll @ end, atol=2e-13)
    np.testing.assert_allclose(rotated_mid, roll @ mid, atol=2e-13)
    np.testing.assert_allclose(rotated_frame_end, roll @ frame_end @ roll.T, atol=2e-13)
    np.testing.assert_allclose(rotated_frame_mid, roll @ frame_mid @ roll.T, atol=2e-13)


@pytest.mark.parametrize("kwargs", [{"l": 0.0, "k1": 1.0, "x_offs": 0.02},
                                    {"l": 1.0, "k1": 0.0, "x_offs": 0.02},
                                    {"l": 1.0, "k1": 1.0}])
def test_degenerate_survey_is_straight(kwargs):
    xy = XYQuadrupole(**kwargs)
    for actual, expected in zip(xy.get_transfer_geometry(), Drift(l=xy.l).get_transfer_geometry()):
        np.testing.assert_array_equal(actual, expected)


def test_slice_boundaries_and_overlong_final_request():
    xy = XYQuadrupole(l=1.0, k1=1.2, x_offs=0.02, y_offs=-0.03)
    section = xy.get_section_tms(0.5, start_l=0.8)[0]
    assert section.delta_length == pytest.approx(0.2)
    expected = _integrated_map(xy.k1, xy.x_offs, xy.y_offs, 0.2, 0, start=0.8)
    np.testing.assert_allclose(section.get_params(0).R, expected, atol=3e-14)
    np.testing.assert_array_equal(xy.get_section_tms(0, start_l=0.5)[0].get_params(0).R, np.eye(6))
    for start, length in [(-0.1, 0.2), (1.1, 0.2), (0.1, -0.2)]:
        with pytest.raises(ValueError, match="slice requires"):
            xy.get_section_tms(length, start_l=start)
