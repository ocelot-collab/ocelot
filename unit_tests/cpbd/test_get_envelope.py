import numpy as np
import matplotlib.pyplot as plt
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.elements import Drift, Bend, Marker
from ocelot.cpbd.beam import Twiss, generate_parray
from ocelot.cpbd.track import track
from ocelot.cpbd.optics import twiss
from ocelot.cpbd.navi import Navigator

def test_twiss_dispersion_correction():
    # Construct a simple lattice
    d = Drift(l=1.2)
    b = Bend(l=0.3, angle=5 / 180 * np.pi)
    m2 = Marker()
    lat = MagneticLattice([d, b, d, m2])

    # Initial Twiss parameters
    tws0 = Twiss(beta_x=10, beta_y=10, E=1, emit_xn=0.5e-6, emit_yn=0.5e-6)

    # Calculate reference Twiss parameters along the lattice
    tws = twiss(lat, tws0)
    tws_ref = tws[-1]

    # Generate a particle array with small chirp
    parray = generate_parray(nparticles=10000, chirp=0.01, tws=tws0)
    navi = Navigator(lat)

    # Track with dispersion correction
    tws_track, _ = track(lat, parray, navi, twiss_disp_correction=True)
    tws_track_end = tws_track[-1]

    # Assert that tracked Twiss parameters are close to reference (dispersion-corrected)
    np.testing.assert_allclose(tws_track_end.beta_x, tws_ref.beta_x, rtol=0.01)
    np.testing.assert_allclose(tws_track_end.beta_y, tws_ref.beta_y, rtol=0.01)
    np.testing.assert_allclose(tws_track_end.Dx, tws_ref.Dx, rtol=0.01)
    np.testing.assert_allclose(tws_track_end.Dxp, tws_ref.Dxp, rtol=0.01)
    np.testing.assert_allclose(tws_track_end.emit_xn, tws_ref.emit_xn, rtol=0.02)
    np.testing.assert_allclose(tws_track_end.emit_yn, tws_ref.emit_yn, rtol=0.02)


def get_envelope_weights(p_array, auto_disp=False, **kwargs):
    """
    Custom function to compute Twiss parameters from a ParticleArray.

    This function processes particle data to compute Twiss parameters,
    accounting for the particle weights, and optionally applies dispersion correction.

    Parameters
    ----------
    p_array : ParticleArray
        Input particle array containing phase-space coordinates.
    auto_disp : bool, optional

    Returns
    -------
    Twiss
        Computed Twiss parameters for the (optionally filtered and corrected) particle array.
    """

    tau = p_array.tau()
    p = p_array.p()
    x = p_array.x()
    px = p_array.px()
    y = p_array.y()
    py = p_array.py()
    q = p_array.q_array

    tws = Twiss()
    tws.q = np.sum(q)
    tws.E = np.copy(p_array.E)
    tws.p = np.average(p, weights=q)

    covdx = np.cov(x, p, aweights=q)
    covdy = np.cov(y, p, aweights=q)
    covdpx = np.cov(px, p, aweights=q)
    covdpy = np.cov(py, p, aweights=q)

    tws.Dx = covdx[0, 1] / covdx[1, 1]
    tws.Dy = covdy[0, 1] / covdy[1, 1]
    tws.Dxp = covdpx[0, 1] / covdpx[1, 1]
    tws.Dyp = covdpy[0, 1] / covdpy[1, 1]

    dx = tws.Dx * p
    dy = tws.Dy * p
    dpx = tws.Dxp * p
    dpy = tws.Dyp * p

    if auto_disp:
        x = x - dx
        px = px - dpx
        y = y - dy
        py = py - dpy

    tws.x = np.average(x, weights=q)
    tws.y = np.average(y, weights=q)
    tws.px = np.average(px, weights=q)
    tws.py = np.average(py, weights=q)
    tws.tau = np.average(tau, weights=q)

    covx = np.cov(x, px, aweights=q)
    tws.xx = covx[0, 0]
    tws.xpx = covx[0, 1]
    tws.pxpx = covx[1, 1]

    covy = np.cov(y, py, aweights=q)
    tws.yy = covy[0, 0]
    tws.ypy = covy[0, 1]
    tws.pypy = covy[1, 1]

    tws.emit_x = np.sqrt(tws.xx * tws.pxpx - tws.xpx ** 2)
    tws.emit_y = np.sqrt(tws.yy * tws.pypy - tws.ypy ** 2)

    relgamma = p_array.E / 0.511e-3
    relbeta = np.sqrt(1 - relgamma ** -2) if relgamma != 0 else 1.
    tws.emit_xn = tws.emit_x * relgamma * relbeta
    tws.emit_yn = tws.emit_y * relgamma * relbeta
    
    tws.beta_x = tws.xx / tws.emit_x
    tws.beta_y = tws.yy / tws.emit_y
    tws.alpha_x = -tws.xpx / tws.emit_x
    tws.alpha_y = -tws.ypy / tws.emit_y

    return tws

def test_twiss_with_custom_function():
    # Construct a simple lattice
    d = Drift(l=1.2)
    b = Bend(l=0.3, angle=5 / 180 * np.pi)
    m2 = Marker()
    lat = MagneticLattice([d, b, d, m2])

    # Initial Twiss parameters
    tws0 = Twiss(beta_x=10, beta_y=10, E=1, emit_xn=0.5e-6, emit_yn=0.5e-6)

    # Calculate reference Twiss parameters along the lattice
    tws = twiss(lat, tws0)
    tws_ref = tws[-1]

    # Generate a particle array with small chirp
    parray = generate_parray(nparticles=10000, chirp=0.01, tws=tws0)
    navi = Navigator(lat)

    # Track with dispersion correction
    tws_track, _ = track(lat, parray, navi, twiss_disp_correction=True,
                         get_twiss=get_envelope_weights)
    tws_track_end = tws_track[-1]

    # Assert that tracked Twiss parameters are close to reference (dispersion-corrected)
    np.testing.assert_allclose(tws_track_end.beta_x, tws_ref.beta_x, rtol=0.01)
    np.testing.assert_allclose(tws_track_end.beta_y, tws_ref.beta_y, rtol=0.01)
    np.testing.assert_allclose(tws_track_end.Dx, tws_ref.Dx, rtol=0.01)
    np.testing.assert_allclose(tws_track_end.Dxp, tws_ref.Dxp, rtol=0.01)
    np.testing.assert_allclose(tws_track_end.emit_xn, tws_ref.emit_xn, rtol=0.02)
    np.testing.assert_allclose(tws_track_end.emit_yn, tws_ref.emit_yn, rtol=0.02)


if __name__ == "__main__":
    test_twiss_dispersion_correction()
    print("Twiss dispersion correction test passed.")