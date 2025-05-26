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




if __name__ == "__main__":
    test_twiss_dispersion_correction()
    print("Twiss dispersion correction test passed.")