import numpy as np
import pytest

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
    from ocelot.common.globals import m_e_GeV

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
