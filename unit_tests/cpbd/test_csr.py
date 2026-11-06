"""CSR charge conservation, MATLAB port regressions and process boundaries."""

from copy import deepcopy

import numpy as np
import pytest

from ocelot.common.globals import m_e_GeV, speed_of_light
from ocelot.cpbd.beam import ParticleArray
from ocelot.cpbd.csr import CSR, CSRConfigurationError, SaferCSR, Smoothing, SubBinning, csr_convolution
from ocelot.cpbd.elements import Bend, Drift, Marker
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.navi import Navigator


@pytest.fixture(scope="module")
def smoothing():
    return Smoothing()


@pytest.fixture(scope="module")
def sub_binning():
    return SubBinning(0, 20, 3)


@pytest.fixture(scope="module")
def csr():
    process = CSR(n_bin=20, m_bin=3, sigma_min=2e-5,
                  traj_step=0.002, apply_step=0.004)
    lat = MagneticLattice([Drift(l=.1, eid="csr_in"), Bend(l=.5, angle=.1, eid="csr_bend"),
                           Drift(l=.1, eid="csr_out")])
    Navigator(lat).add_physics_proc(process, lat.sequence[0], lat.sequence[-1])
    process.z0 = .4
    return process


@pytest.fixture
def particles():
    rng = np.random.default_rng(711)
    p = ParticleArray(2000)
    p.E = .13
    p.rparticles[4] = rng.normal(0, 2e-4, p.n)
    p.q_array[:] = rng.uniform(.1, 3, p.n)
    p.q_array *= 1e-9 / p.q_array.sum()
    return p


@pytest.mark.parametrize("a,b", [([1.], [2.]), ([1., 2.], [3., 4.]),
                               ([1., -2., 3.], [4., 0., -1., 2.])])
def test_convolution_includes_last_sample(a, b):
    np.testing.assert_allclose(csr_convolution(np.array(a), np.array(b)),
                               np.convolve(a, b), atol=1e-14)


@pytest.mark.parametrize("method", [0, 1, 2])
@pytest.mark.parametrize("m_bin", [1, 5])
@pytest.mark.parametrize("edge_only", [False, True])
def test_smoothing_preserves_charge_at_boundaries(smoothing, method, m_bin, edge_only):
    n_bin = 10
    k = n_bin * m_bin
    edges = np.linspace(-1e-3, 1e-3, k + 1)
    q = np.linspace(1, 3, k) * 1e-12
    if edge_only:
        q[1:-1] = 0
    _, _, _, density = smoothing.Q2EQUI(q, [0, n_bin, m_bin, method, .5, 1e-5],
                                       edges, np.ones(k))
    assert np.all(np.isfinite(density))
    assert np.all(density >= 0)
    np.testing.assert_allclose(density.sum(), q.sum(), rtol=2e-13, atol=0)


def test_gaussian_sampler_includes_upper_endpoint(smoothing):
    args = (1, np.array([2.]), np.array([-.1]), np.array([.1]), 5,
            np.array([.2]), .05, 40, -1.)
    points = -1 + np.arange(40) * .05
    expected = 2 * .05 / (np.sqrt(2 * np.pi) * .2) * np.exp(-.5 * (points / .2)**2)
    expected *= 2 / expected.sum()
    for sampler in (smoothing.q_per_step_ip2_py, smoothing.q_per_step_ip2):
        actual = sampler(*args)
        np.testing.assert_allclose(actual, expected, rtol=3e-14, atol=0)
        assert actual[-1] > 0


@pytest.mark.parametrize("edges", [[0., .25, .5, .75, 1.],
                                   [0., 0., .5, .5, 1., 1.], [0., 1.]])
def test_sorted_bin_counts_include_repeated_boundary_particles(sub_binning, edges):
    s = np.array([0., 0., .2, .5, .5, 1., 1.])
    edges = np.array(edges)
    expected, _ = np.histogram(s, bins=edges)
    actual = sub_binning.p_per_subbins(s, edges, len(edges) - 1)
    np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("method", [0, 1, 2])
@pytest.mark.parametrize("scalar_charge", [False, True])
def test_smoothing_empty_bins_match_aggregated_charges(smoothing, method, scalar_charge):
    counts = np.array([0, 2, 0, 1, 0, 0, 3, 0])
    edges = np.linspace(-1e-3, 1e-3, len(counts) + 1)
    if scalar_charge:
        q = 1e-12
        charges_per_bin = counts * q
    else:
        q = np.array([1., 2., 4., 8., 16., 32.]) * 1e-12
        charges_per_bin = np.array([0., 3., 0., 4., 0., 0., 56., 0.]) * 1e-12
    params = [0, len(counts), 1, method, .5, 1e-5]
    actual = smoothing.Q2EQUI(q, params, edges, counts)[3]
    expected = smoothing.Q2EQUI(charges_per_bin, params, edges, np.ones(len(counts)))[3]
    np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=1e-27)
    np.testing.assert_allclose(actual.sum(), charges_per_bin.sum(), rtol=2e-13)


@pytest.mark.parametrize("x_qbin", [0, .5, 1])
def test_weighted_sub_bins_are_symmetric_under_reversal(sub_binning, x_qbin):
    s = np.array([-3., -1., 0., 2., 4.])
    q = np.array([1., 5., 2., 3., 7.])
    left, _ = sub_binning.subbin_bound(q, s, x_qbin, 4, 3)
    right, _ = sub_binning.subbin_bound(q[::-1], -s[::-1], x_qbin, 4, 3)
    np.testing.assert_allclose(left, -right[::-1], atol=2e-14)


@pytest.mark.parametrize("x_qbin", [0, .5, 1])
def test_weighted_kick_does_not_depend_on_particle_order(csr, particles, monkeypatch, x_qbin):
    monkeypatch.setattr(csr, "charge_binning_fraction", x_qbin)
    perm = np.random.default_rng(13).permutation(particles.n)
    shuffled = deepcopy(particles)
    shuffled.rparticles = shuffled.rparticles[:, perm].copy()
    shuffled.q_array = shuffled.q_array[perm].copy()
    csr.apply(particles, .02)
    csr.apply(shuffled, .02)
    np.testing.assert_allclose(particles.p(), shuffled.p()[np.argsort(perm)], rtol=1e-13, atol=1e-18)
    assert np.linalg.norm(particles.p()) > 0


def test_mesh_step_quantum_controls_actual_wake_mesh(csr, particles, monkeypatch):
    spacings = []
    original = csr.CSR_K1

    def record_mesh(i, traj, ndw, gamma):
        spacings.append(ndw[1])
        return original(i, traj, ndw, gamma)

    monkeypatch.setattr(csr, "CSR_K1", record_mesh)
    monkeypatch.setattr(csr, "mesh_step_quantum", 3e-6)
    csr.apply(particles, .02)
    assert spacings
    np.testing.assert_allclose(np.array(spacings) / csr.mesh_step_quantum,
                               np.round(np.array(spacings) / csr.mesh_step_quantum), atol=1e-13)


@pytest.mark.parametrize("process_type", [CSR, SaferCSR])
@pytest.mark.parametrize("name,alias,value", [("charge_binning_fraction", "x_qbin", .5),
                                             ("mesh_step_quantum", "step_unit", 2e-6)])
def test_renamed_parameters_require_new_names(process_type, name, alias, value):
    process = process_type(**{name: value})
    assert getattr(process, name) == value
    assert not hasattr(process, alias)
    with pytest.raises(TypeError, match=name):
        process_type(**{alias: value})
    with pytest.raises(AttributeError, match=name):
        setattr(process, alias, value)


@pytest.mark.parametrize("extent,sigma,expected_step", [(1e-4, 21.2e-6, 6e-6),
                                                       (1e-5, 2.8e-6, 2e-6)])
def test_mesh_step_quantum_documented_examples(smoothing, extent, sigma, expected_step):
    edges = np.linspace(-extent, extent, 51)
    z1, z2, nz, _ = smoothing.Q2EQUI(1e-12, [0, 10, 5, 2, .5, sigma, 2e-6],
                                     edges, np.ones(50))
    step = (z2 - z1) / nz
    np.testing.assert_allclose(step, expected_step, rtol=1e-14)
    np.testing.assert_allclose(np.array([z1, z2]) / step,
                               np.round(np.array([z1, z2]) / step), atol=1e-13)


def test_kernel_interpolates_sorted_distinct_retarded_coordinates(csr, monkeypatch):
    w = np.array([-.01, -.03, -.02, -.02])
    monkeypatch.setattr(csr.k0_fin_anf, "eval", lambda *args: (w, w*w))
    actual = csr.CSR_K1(1, csr.csr_traj, [4, .005], gamma=200)
    # Piecewise-linear interpolation of w**2 at the five requested positions.
    interpolated = np.array([.00065, .0004, .00025, .0001, .0001])
    expected = np.diff(np.diff(interpolated, append=0), append=0) / .005
    expected *= 1e-7 * speed_of_light**2
    np.testing.assert_allclose(actual, expected, atol=1e-8)


def test_analytic_extension_selects_branch_from_first_sample(csr):
    radius, angle = 5., .1
    traj = np.array([[0, radius*angle], [0, radius*(1-np.cos(angle))], [0, 0],
                     [0, radius*np.sin(angle)], [0, np.sin(angle)], [0, 0], [1, np.cos(angle)]])
    gamma = .13 / m_e_GeV
    beta = np.sqrt(1 - gamma**-2)
    w = np.array([-.0014393, -.0014392])
    winf = traj[3, -1] - traj[0, -1]
    a2 = traj[1, -1]**2
    s = winf + gamma*(gamma*(w-winf) - beta*np.sqrt(gamma**2*(w-winf)**2 + a2))
    r = (w-s) / beta
    assert a2/r[0]**2 < 1e-7 < a2/r[1]**2
    cosine = np.cos(angle)
    expected = (beta*(1-cosine)*np.log(r[0]/r)
                - (beta**2*cosine-1)*np.log((winf-s+r)/(winf-s[0]+r[0]))
                + gamma**-2*np.log(w[0]/w))
    np.testing.assert_allclose(csr.K0_fin_inf(1, traj, w, gamma), expected, rtol=1e-10, atol=1e-20)


@pytest.mark.parametrize("method", [0, 1, 2])
@pytest.mark.parametrize("mesh_step_quantum", [0, 3e-6])
def test_mesh_only_helper_matches_deposition(smoothing, method, mesh_step_quantum):
    edges = np.linspace(-1e-3, 1e-3, 101)
    params = [0, 20, 5, method, .5, 1e-5, mesh_step_quantum]
    result = smoothing.Q2EQUI(1e-12, params, edges, np.ones(100))
    assert smoothing.Q2EQUI_modified(params, edges) == result[:3]


def test_sub_trajectory_kicks_do_not_discard_path_length(csr, particles, monkeypatch):
    monkeypatch.setattr(csr, "CSR_K1", lambda i, traj, ndw, gamma: np.ones(int(ndw[0]) + 1))
    twice = deepcopy(particles)
    ds = csr.traj_step / 4
    csr.apply(particles, 2*ds)
    csr.apply(twice, ds)
    csr.apply(twice, ds)
    assert np.linalg.norm(particles.p()) > 0
    np.testing.assert_allclose(twice.p(), particles.p(), rtol=1e-14, atol=0)


def test_empty_array_and_zero_length_are_noops(csr, particles):
    csr.apply(ParticleArray(), .01)
    before = particles.rparticles.copy()
    csr.apply(particles, 0.)
    np.testing.assert_array_equal(before, particles.rparticles)


@pytest.mark.parametrize("n", [1, 10])
def test_degenerate_bunch_raises_clear_error(csr, n):
    p = ParticleArray(n)
    p.E = .13
    p.q_array[:] = 1e-12
    with pytest.raises(CSRConfigurationError, match="longitudinal extent"):
        csr.apply(p, .01)


@pytest.mark.parametrize("name", ["filter_order", "n_mesh", "debug"])
def test_removed_parameters_are_rejected(name):
    with pytest.raises(TypeError, match=name):
        CSR(**{name: 10})


@pytest.mark.parametrize("method", [0, 1, 2])
def test_coarse_quantized_mesh_preserves_charge(smoothing, method):
    edges = np.linspace(-.00151, .00151, 51)
    result = smoothing.Q2EQUI(1., [0, 10, 5, method, .5, 1e-5, .001], edges, np.ones(50))
    np.testing.assert_allclose(result[3].sum(), 50., rtol=2e-13)


@pytest.mark.parametrize("x,n,m", [(-1, 10, 2), (2, 10, 2), (0, 0, 2), (0, 10, 0), (0, 10.5, 2)])
def test_invalid_binning_is_rejected(sub_binning, x, n, m):
    with pytest.raises(CSRConfigurationError):
        sub_binning.subbin_bound(1., np.array([0., 1.]), x, n, m)


def test_safer_csr_uses_local_origin_and_actual_trajectory_spacing(monkeypatch):
    bend = Bend(l=.503, angle=.1, eid="safe_bend")
    lat = MagneticLattice([Drift(l=2., eid="safe_prefix"), bend, Marker(eid="safe_end")])
    process = SaferCSR(apply_kick=False, traj_step=.007)
    Navigator(lat).add_physics_proc(process, bend, lat.sequence[-1])
    radii = []

    def criterion(x, z, radius):
        radii.append(radius)
        return 0.

    monkeypatch.setattr(process, "_derbenev_criterion", criterion)
    for local_s in (0., .01, .25, bend.l, bend.l + 1e-12):
        process.z0 = process.z_csr_start + local_s
        process.query_derbenev_criterion(1e-4, 1e-4, 1e-4)
    assert len(radii) == 5
    np.testing.assert_allclose(radii, bend.l/bend.angle, rtol=1e-11)


def test_safer_csr_warning_does_not_divide_by_rounded_zero(caplog):
    process = SaferCSR(sigma_min=1e-14, apply_kick=False)
    process.query_bunch_length(1e-15)
    assert "Bunch is too small" in caplog.text
