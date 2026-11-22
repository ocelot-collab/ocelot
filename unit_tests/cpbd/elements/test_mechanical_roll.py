import numpy as np
import pytest

from ocelot.cpbd.elements import Bend, Octupole, Quadrupole, RBend, SBend, Sextupole, XYQuadrupole
from ocelot.cpbd.errors import errors_seed
from ocelot.cpbd.latticeIO import LatticeIO
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.r_matrix import rot_mtx
from ocelot.cpbd.transformations import KickTM, RungeKuttaGlobalTM, RungeKuttaOcelotTM, SecondTM, TransferMap


def track_maps(maps, initial, energy=1.):
    state = initial.copy()
    for tm in maps:
        tm.map_function(state, energy)
    return state


@pytest.mark.parametrize("method", [TransferMap, SecondTM, KickTM, RungeKuttaGlobalTM])
@pytest.mark.parametrize("family,strength", [(Quadrupole, {"k1": 2.}), (Sextupole, {"k2": 30.}),
                                            (Octupole, {"k3": 80.})])
def test_straight_magnet_rolls_about_displaced_axis(method, family, strength):
    common = dict(l=.4, dx=.001, dy=-.002, tm=method, npoints=11, **strength)
    magnet = family(tilt=.2, **common)
    before = magnet.first_order_tms[0]
    before.get_params(1.)
    magnet.tms[0].get_params(1.)
    magnet.dtilt = -.31
    assert magnet.element.dtilt == -.31 and magnet.tilt == .2
    assert before is not magnet.first_order_tms[0]
    assert magnet.B(1.)[0].shape == (6, 1)
    expected = family(tilt=magnet.tilt + magnet.dtilt, **common)
    x = np.array([[.003, -.001, .004, .002, 0., .001]]).T
    np.testing.assert_allclose(track_maps(magnet.tms, x), track_maps(expected.tms, x), atol=2e-15)
    np.testing.assert_allclose(magnet.R(1.), expected.R(1.), atol=2e-15)
    if method is RungeKuttaGlobalTM:
        assert magnet.tms[0].npoints == 11


# MAD-X 5.09.03 TWISS, positron at 1 GeV, L=.4, ANGLE=.1, EALIGN,DPSI=.001.
MADX_EXIT = np.array([9.9417518360e-9, 4.9667328818e-8, -1.9983335458e-5, -9.9833400505e-5])


@pytest.mark.parametrize("method", [TransferMap, SecondTM, KickTM, RungeKuttaOcelotTM])
def test_sector_roll_matches_madx_on_nominal_planes(method):
    bend = SBend(l=.4, angle=.1, tilt=.32, tm=method)
    initial = np.zeros((6, 1))
    survey = bend.get_transfer_geometry()
    assert np.max(np.abs(track_maps(bend.tms, initial))) < 1e-12
    bend.dtilt = .001
    np.testing.assert_array_equal(bend.get_transfer_geometry()[0], survey[0])
    np.testing.assert_array_equal(bend.get_transfer_geometry()[1], survey[1])
    expected = rot_mtx(-bend.tilt)[:4, :4] @ MADX_EXIT
    full = track_maps(bend.tms, initial)
    np.testing.assert_allclose(full[:4, 0], expected, rtol=0, atol=9e-13)
    pieces = [tm for s in (0., .1, .2, .3) for tm in bend.get_section_tms(.1, start_l=s)]
    np.testing.assert_allclose(track_maps(pieces, initial), full, rtol=0, atol=1e-12)
    bend.dtilt = 0.
    assert np.max(np.abs(track_maps(bend.tms, initial))) < 1e-12


@pytest.mark.parametrize("family", [Bend, RBend, SBend])
def test_roll_changes_lattice_map_but_preserves_survey_and_serializes(family):
    bend = family(l=.7, angle=-.23, k1=.3, e1=.05, tilt=.5, dx=.001, dy=-.002, eid="B")
    lattice = MagneticLattice([bend], method={"global": SecondTM})
    original = lattice.transfer_maps(1.)
    survey = lattice.survey()
    bend.dtilt = .012
    rolled = lattice.transfer_maps(1.)
    for before, after in zip(original, rolled):
        assert not np.allclose(before, after)
    for actual, expected in zip(lattice.survey(), survey):
        for row, reference in zip(actual, expected):
            for key in ("W", "r_end", "TILT"):
                np.testing.assert_array_equal(row[key], reference[key])
    bend.name = "b"
    namespace = {family.__name__: family}
    exec(LatticeIO.element_def_string(bend), namespace)
    restored = namespace["b"]
    assert (restored.tilt, restored.dtilt, restored.dx, restored.dy) == (.5, .012, .001, -.002)
    for actual, expected in zip(MagneticLattice([restored], method={"global": SecondTM}).transfer_maps(1.), rolled):
        np.testing.assert_array_equal(actual, expected)


def test_roll_only_random_seed_reaches_cached_bend_maps():
    bend = SBend(l=.4, angle=.1, tilt=.3, eid="B")
    lattice = MagneticLattice([bend])
    lattice.transfer_maps(1.)
    _, errors = errors_seed(lattice, {SBend: {"offset": 0., "dtilt": .001}})
    assert errors["B"]["dtilt"] == bend.element.dtilt
    assert bend.dx == bend.dy == 0.
    assert bend.tilt == .3
    assert lattice.transfer_maps(1.)[0][3, 0] == pytest.approx(
        -np.cos(.3)*np.sin(.1)*np.sin(bend.dtilt), rel=.002, abs=1e-14)


def test_xyquadrupole_survey_and_slices_use_original_reference():
    q = XYQuadrupole(l=.6, k1=1.3, x_offs=.03, y_angle=.01, tilt=.3)
    survey = q.get_transfer_geometry()
    q.dtilt = .004
    for actual, expected in zip(q.get_transfer_geometry(), survey):
        np.testing.assert_array_equal(actual, expected)
    full = track_maps(q.tms, np.zeros((6, 1)))
    assert np.linalg.norm(full[:4]) > 1e-6
    pieces = q.get_section_tms(.25) + q.get_section_tms(.35, start_l=.25)
    np.testing.assert_allclose(track_maps(pieces, np.zeros((6, 1))), full, atol=3e-12)


def test_rolled_bend_quadratic_coefficients_match_tracking_through_physical_faces():
    from ocelot.common.geometry import get_tilt_matrix
    from ocelot.common.globals import m_e_GeV

    energy = .002
    beta = np.sqrt(1. - (m_e_GeV / energy)**2)
    bend = Bend(l=.7, angle=-.25, tilt=.3, k1=.7, k2=2., dx=.001, dy=-.002, tm=SecondTM)
    position, axes, _, _ = bend.get_transfer_geometry()
    native = bend.tms[1].get_params(energy)
    b, r, t = native.B[:, 0], native.get_rotated_R(), native.get_rotated_T()
    bend.dtilt = .08
    roll = get_tilt_matrix(bend.dtilt)
    pivot = np.array([bend.dx, bend.dy, 0.])

    def project(x, origin, basis, dest, dest_basis):
        # Cartesian ray-plane intersection, with direction and path in metres.
        momentum = np.sqrt(1. + 2.*x[5]/beta + x[5]**2)
        direction = basis @ np.array([x[1], x[3], np.sqrt(momentum**2 - x[1]**2 - x[3]**2)]) / momentum
        point = origin + basis[:, :2] @ x[[0, 2]]
        normal = dest_basis[:, 2]
        distance = np.dot(dest - point, normal) / np.dot(direction, normal)
        xy = dest_basis[:, :2].T @ (point + distance*direction - dest)
        pxy = momentum * (dest_basis[:, :2].T @ direction)
        return np.array([xy[0], pxy[0], xy[1], pxy[1], x[4] + distance*(1./beta + x[5])/momentum, x[5]])

    def track(x):
        zero = np.zeros(3)
        u = project(x, zero, np.eye(3), pivot - roll @ pivot, roll)
        v = b + r @ u + np.einsum("ijk,j,k->i", t, u, u)
        return project(v, pivot + roll @ (position-pivot), roll @ axes, position, axes)

    B, R, T = bend.B(energy)[1][:, 0], bend.R(energy)[1], bend.T(energy)[1]
    origin = track(np.zeros(6))
    np.testing.assert_allclose(B, origin, rtol=0., atol=1e-15)
    # Recover every linear column and every symmetric quadratic coefficient.
    basis = np.eye(6)
    directions = [*basis, *(basis[i] + basis[j] for i in range(6) for j in range(i))]
    step = 1.e-4
    for v in directions:
        plus, minus = track(step*v), track(-step*v)
        np.testing.assert_allclose(R @ v, (plus-minus)/(2*step), rtol=0., atol=7e-9)
        np.testing.assert_allclose(np.einsum("ijk,j,k->i", T, v, v),
                                   (plus+minus-2*origin)/(2*step**2), rtol=0., atol=2e-8)


def test_load_magnet_with_legacy_wrapper_roll():
    import pickle

    bend = SBend(l=.4, angle=.1, tm=SecondTM)
    bend.B(1.)
    # Reproduce the state written before the atom owned dtilt.
    del bend.element.dtilt
    object.__setattr__(bend, "dtilt", .001)
    restored = pickle.loads(pickle.dumps(bend))
    assert "dtilt" not in vars(restored)
    assert restored.element.dtilt == .001
    np.testing.assert_allclose(track_maps(restored.tms, np.zeros((6, 1)))[:4, 0],
                               MADX_EXIT, rtol=0, atol=9e-13)
    restored.dtilt = 0.
    np.testing.assert_allclose(restored.B(1.), 0., atol=0.)
