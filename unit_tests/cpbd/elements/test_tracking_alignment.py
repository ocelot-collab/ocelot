import numpy as np
import pytest

from ocelot.common.globals import m_e_GeV
from ocelot.cpbd.elements import Bend, Drift, Multipole, Quadrupole
from ocelot.cpbd.r_matrix import rot_mtx
from ocelot.cpbd.transformations import KickTM, RungeKuttaGlobalTM, RungeKuttaOcelotTM, SecondTM, TransferMap


def probes():
    return np.array([[0., .001, -.002], [.0002, 0., -.0001], [0., -.002, .003],
                     [-.0001, .0003, 0.], [0., .004, -.002], [.001, -.02, 0.]])


@pytest.mark.parametrize("method", [TransferMap, SecondTM, KickTM, "thin"])
def test_straight_map_translates_and_rolls_its_native_kick(method):
    dx, dy, tilt, dtilt = .0013, -.0021, .4, -.13
    shift = np.array([[dx, 0., dy, 0., 0., 0.]]).T
    rotation = rot_mtx(tilt + dtilt)
    if method == "thin":
        family, kwargs = Multipole, dict(kn=[.04, .6, 30., -80.])
    else:
        family, kwargs = Quadrupole, dict(l=.4, k1=.6, k2=30., tm=method)
    native = family(**kwargs)
    magnet = family(dx=dx, dy=dy, tilt=tilt, dtilt=dtilt, **kwargs)
    local = rotation @ (probes() - shift)
    native.apply(local, energy=1.)
    actual = probes()
    magnet.apply(actual, energy=1.)
    np.testing.assert_allclose(actual, shift + rotation.T @ local, atol=2e-16, rtol=0.)


def test_thin_quadrupole_optics_and_tracking_share_alignment_after_assignment():
    thin = Multipole(kn=[0., .6], tilt=.3)
    thin.first_order_tms[0].get_params(1.)
    thin.dx, thin.dy, thin.dtilt = .001, -.002, -.17
    expected = thin.R(1.)[0] @ probes() + thin.B(1.)[0]
    actual = probes()
    thin.apply(actual, energy=1.)
    np.testing.assert_allclose(actual, expected, atol=2e-16, rtol=0.)


@pytest.mark.parametrize("method", [TransferMap, SecondTM])
@pytest.mark.parametrize("index", [0, 2], ids=["entrance", "exit"])
def test_bend_fringe_translates_the_rolled_map(index, method):
    kwargs = dict(l=.7, angle=-.23, tilt=.42, k1=.8, gap=.03, fint=.5, fintx=.6,
                  e1=.12, e2=-.09, h_pole1=.3, h_pole2=-.4, tm=method)
    native = Bend(**kwargs)
    offset = Bend(dx=.001, dy=-.002, **kwargs)
    r, t = native.R(1.)[index], native.T(1.)[index]
    d = np.array([offset.dx, 0., offset.dy, 0., 0., 0.])
    np.testing.assert_allclose(offset.B(1.)[index][:, 0],
                               d - r @ d + np.einsum("ijk,j,k->i", t, d, d), atol=2e-17, rtol=0.)
    np.testing.assert_allclose(offset.R(1.)[index], r - (t + t.swapaxes(1, 2)) @ d, atol=2e-15, rtol=0.)
    np.testing.assert_allclose(offset.T(1.)[index], t, atol=2e-15, rtol=0.)
    x = probes() - d[:, None]
    native.tms[index].map_function(x, 1.)
    np.testing.assert_allclose(offset.tms[index].map_function(probes(), 1.), x + d[:, None], atol=2e-15)


def test_offset_quadrupole_nominal_reference_keeps_orbit_kick():
    kwargs = dict(l=.4, k1=2., dx=.001, dy=-.0003, tilt=.3, npoints=101)
    initial = np.zeros((6, 1))
    expected = RungeKuttaGlobalTM.from_element(Quadrupole(**kwargs).element, npoints=101)
    expected.map_function(initial, 1.)
    assert abs(initial[1, 0]) > .0001
    q = Quadrupole(**kwargs, tm=RungeKuttaOcelotTM, reference="nominal")
    np.testing.assert_allclose(q.tms[0].map_function(np.zeros_like(initial), 1.), initial, atol=2e-13)
    for mode in ({}, {"reference": "field"}):
        following = Quadrupole(**kwargs, tm=RungeKuttaOcelotTM, **mode)
        np.testing.assert_allclose(following.tms[0].map_function(np.zeros_like(initial), 1.), 0., atol=2e-13)


def test_nominal_reference_uses_surveyed_bend_faces_and_original_slice_origin():
    # A field-free surveyed bend has an exactly straight physical ray.
    b = Bend(l=.5, angle=.14, tilt=.37, mag_field=lambda x, y, z: (0., 0., 0.))
    b.set_tm(RungeKuttaOcelotTM, reference="nominal", npoints=31)
    energy = .002
    beta = np.sqrt(1. - (m_e_GeV / energy)**2)
    origin, frame, _, _ = b.get_transfer_geometry()
    normal = frame[:, 2]
    z = origin @ normal / normal[2]
    x, y = (frame.T @ (np.array([0., 0., z]) - origin))[:2]
    expected = [x, frame[2, 0], y, frame[2, 1], (z - b.l)/beta, 0.]
    actual = np.zeros((6, 1))
    b.apply(actual, energy)
    np.testing.assert_allclose(actual[:, 0], expected, atol=2e-13)
    sliced = np.zeros_like(actual)
    for start in [0., .125, .25, .375]:
        for tm in b.get_section_tms(.125, start_l=start):
            tm.map_function(sliced, energy)
    np.testing.assert_allclose(sliced, actual, atol=2e-13)


def test_rk_rejects_unknown_reference():
    with pytest.raises(ValueError, match="reference"):
        Drift(l=.3, tm=RungeKuttaOcelotTM, reference="unknown")


def test_nominal_rk_requires_a_forward_survey_throughout_the_slice():
    bend = Bend(l=1., angle=2*np.pi + .1, tm=RungeKuttaOcelotTM, reference="nominal")
    with pytest.raises(ValueError, match="forward"):
        bend.apply(np.zeros((6, 1)), energy=1.)


@pytest.mark.parametrize("energy", [0., .002])
def test_thin_dipole_dispersion_and_arrival_time_are_conjugate(energy):
    # The path through a thin positive bend is longer for an outward ray.
    # In Ocelot's coordinates R_tau,x = R_px,delta = angle / beta0.
    beta = np.sqrt(1. - (m_e_GeV / energy)**2) if energy else 1.
    bend = Multipole(kn=[.07, .4-.1j], dx=.003, dy=-.002, tilt=.21, dtilt=-.1)
    affine = bend.R(energy)[0]
    direction = np.array([np.cos(.11), np.sin(.11)])
    np.testing.assert_allclose(affine[[1, 3], 5], .07 * direction / beta)
    np.testing.assert_allclose(affine[4, [0, 2]], .07 * direction / beta)
    # (tau, -delta) is the canonical longitudinal pair.
    form = np.kron(np.diag([1., 1., -1.]), [[0., 1.], [-1., 0.]])
    np.testing.assert_allclose(affine.T @ form @ affine, form, atol=5e-16)
    expected = affine @ probes() + bend.B(energy)[0]
    actual = probes()
    bend.apply(actual, energy)
    np.testing.assert_allclose(actual, expected, atol=3e-17)


def test_sector_kicks_converge_to_linear_hamiltonian_with_pole_faces():
    from scipy.linalg import expm

    # Canonical paraxial equations in the horizontal bending plane, with
    # finite velocity, sector weak focusing and geometric time of flight.
    energy, length, angle, k1 = .002, .6, .17, .4
    beta = np.sqrt(1. - (m_e_GeV / energy)**2)
    h = angle / length
    generator = np.zeros((6, 6))
    generator[0, 1] = generator[2, 3] = 1.
    generator[1, 0], generator[3, 2] = -(h*h + k1), k1
    generator[1, 5] = generator[4, 0] = h / beta
    generator[4, 5] = -(m_e_GeV / (energy * beta))**2
    e1, e2 = .05, -.07
    entrance, exit_face = np.eye(6), np.eye(6)
    for face, edge in [(entrance, e1), (exit_face, e2)]:
        face[1, 0] = h * np.tan(edge)
        face[3, 2] = -h * np.tan(edge)
    expected = exit_face @ expm(length * generator) @ entrance @ probes()
    bend = Bend(l=length, angle=angle, k1=k1, e1=e1, e2=e2, tm=KickTM, nkick=256)
    actual = probes()
    bend.apply(actual, energy)
    np.testing.assert_allclose(actual, expected, atol=4e-9, rtol=0.)


def test_nominal_rk_bend_local_offset_against_displaced_circular_orbit():
    # In the native sector an orbit launched at -dx travels on a circle of
    # unchanged radius. Intersect this circle with the exit normal plane.
    energy, angle, length, dx = .002, .1, .4, .001
    radius = length / angle
    beta = np.sqrt(1. - (m_e_GeV / energy)**2)
    dphi = -np.arcsin(dx * np.sin(angle) / radius)
    expected = [dx * (1. - np.cos(angle)) - 2. * radius * np.sin(dphi / 2.)**2,
                -np.sin(dphi), 0., 0., radius * dphi / beta, 0.]
    bend = Bend(l=length, angle=angle, dx=dx, tm=RungeKuttaOcelotTM, reference="nominal", npoints=51)
    actual = np.zeros((6, 1))
    bend.apply(actual, energy)
    np.testing.assert_allclose(actual[:, 0], expected, atol=2e-13, rtol=0.)

    bend.tilt, bend.dtilt, bend.dy = -.21, .003, -.0002
    full, sliced = probes(), probes()
    bend.apply(full, energy)
    for start in [0., .1, .2, .3]:
        for tm in bend.get_section_tms(.1, start_l=start):
            tm.map_function(sliced, energy)
    np.testing.assert_allclose(sliced, full, atol=2e-13, rtol=0.)

    # A field-specified reference keeps its native orbit, including offsets.
    field_reference = Bend(l=length, angle=angle, dx=dx, tm=RungeKuttaOcelotTM)
    origin = np.zeros((6, 1))
    field_reference.apply(origin, energy)
    np.testing.assert_allclose(origin, 0., atol=2e-13)


def test_face_rejects_negative_energy_branch_and_preserves_transverse_only_timing():
    bend = Bend(l=.4, angle=.15, dx=.001, dtilt=.01, tm=RungeKuttaOcelotTM, reference="nominal")
    face = bend.element._rolled_face(.4)
    unphysical = probes()
    unphysical[5] = -3.
    with pytest.raises(ValueError, match="positive particle energy"):
        face.apply(unphysical, 1.)
    with pytest.raises(ValueError, match="rest energy"):
        face.apply(probes(), m_e_GeV)
    points = probes()
    for tm in bend.tms:
        tm.long_dynamics = False
    bend.apply(points, .002)
    np.testing.assert_array_equal(points[4], probes()[4])


def test_lattice_kick_resolution_and_symplecticity_on_a_displaced_orbit():
    from ocelot.cpbd.elements import Hcor
    from ocelot.cpbd.magnetic_lattice import MagneticLattice
    from ocelot.cpbd.optics import MethodTM

    bend = Bend(l=.6, angle=.2, k1=1., k2=11., tilt=.4, dtilt=.08,
                dx=.002, dy=.003, e1=.003, e2=.05, tm=KickTM, nkick=6)
    corrector = Hcor(angle=.001)
    method = MethodTM()
    method.global_method = KickTM
    MagneticLattice([bend, corrector], method=method)
    assert bend.tms[1].nkick == 6
    method.nkick = 12
    MagneticLattice([bend, corrector], method=method)
    assert bend.tms[1].nkick == 12
    assert type(corrector.tms[0]) is TransferMap

    start = np.array([.01, .004, -.002, .002, .0005, -.15])[:, None]
    step = 1e-6
    points = np.column_stack((start + step*np.eye(6), start - step*np.eye(6)))
    bend.apply(points, .002)
    jac = (points[:, :6] - points[:, 6:]) / (2*step)
    form = np.kron(np.diag([1., 1., -1.]), [[0., 1.], [-1., 0.]])
    np.testing.assert_allclose(jac.T @ form @ jac, form, atol=5e-10, rtol=0.)

    # Two half maps of 12 kicks put the kicks at the same 24 positions
    # as the full map. Entrance and exit pole faces are crossed once.
    split = start.copy()
    for offset in [0., .3]:
        for tm in bend.get_section_tms(.3, start_l=offset):
            tm.map_function(split, .002)
    bend.set_tm(KickTM, nkick=24)
    bend.apply(start, .002)
    np.testing.assert_allclose(split, start, atol=1e-14, rtol=0.)
