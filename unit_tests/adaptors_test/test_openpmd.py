import pytest
import numpy as np
import os

from ocelot.common.globals import speed_of_light as c
from ocelot.common.globals import m_e_GeV
from ocelot import MagneticLattice, Twiss, SecondTM, track, generate_parray, Navigator
from ocelot import SBend, Drift, Marker
from ocelot.adaptors.openpmd import write_parray_to_openpmd, read_openpmd_to_parray, \
    SaveBeamOPMD, openpmd_viewer_installed, openpmd_api_installed

# Define decorator to skip tests if both openpmd_viewer and openpmd_api
# are not installed.
openpmd_installed = openpmd_viewer_installed & openpmd_api_installed
only_if_openpmd_installed = pytest.mark.skipif(
    not openpmd_installed, reason='Both openPMD-api and openPMD-viewer required to run tests')


def generate_example_beam():
    # Generate ocelot parray.
    E_ref_GeV = 1.0  # reference energy in GeV
    gamma = E_ref_GeV / 0.511e-3
    sigma_p = 1e-2
    sigma_tau = 1e-6  # in meters
    ex_n = 1e-6
    ey_n = 1e-6
    current = 1000
    charge = current * 2 * sigma_tau / c
    NP = int(1e5)

    # initialization of Twiss object
    tws0 = Twiss()
    tws0.beta_x = 10
    tws0.beta_y = 10
    tws0.alpha_x = 0.25
    tws0.alpha_y = 0.25
    tws0.E = E_ref_GeV
    tws0.emit_x = ex_n / gamma
    tws0.emit_y = ey_n / gamma

    # Generate the beam
    charge = current * np.sqrt(2 * np.pi) * sigma_tau / c
    p_array = generate_parray(sigma_tau=sigma_tau, sigma_p=sigma_p,
                              tws=tws0, chirp=0, charge=charge,
                              nparticles=NP)

    return p_array


@only_if_openpmd_installed
def test_write_parray_to_openpmd():
    """
    Writes an ocelot parray to disk in openpmd format.
    Then it reads the array from disk and compares with the original.

    """

    p_array = generate_example_beam()

    # Beam diagnostics data
    beam_outdir = os.path.join('diags_0', 'hdf5')
    os.makedirs(beam_outdir, exist_ok=True)

    # writes p_array in openpmd format
    write_parray_to_openpmd(p_array, folder_path=beam_outdir, species='beam')

    # reads the beam back
    p_array_read = read_openpmd_to_parray(beam_outdir, species='beam',
                                          gamma_ref=p_array.E / m_e_GeV, z_ref=p_array.s)

    # Check that the read beam is equal to the original.
    assert_equal_beams(p_array, p_array_read)


@only_if_openpmd_installed
def test_write_parray_to_openpmd_with_beamline():
    """
    Tests openPMD writting capabilities in a beamline simulation.
    Then it reads the array from disk and compares with the original.

    """

    # Chicane
    R56_mod = 1e-3  # Target R56
    Ldip = 0.6  # dipole length
    Ld = 0.5  # distance between the central and outer dipoles
    ch_angle = np.sqrt(R56_mod / (2 * (Ld + 2 * Ldip / 3)))

    ch1_B1 = SBend(l=Ldip, angle=ch_angle, e1=0.0, e2=ch_angle, tilt=0.0, fint=0.0, eid='ch1_B1')
    ch1_B2 = SBend(l=Ldip, angle=-ch_angle, e1=-ch_angle, e2=0.0, tilt=0.0, fint=0.0, eid='ch1_B2')
    ch1_B3 = SBend(l=Ldip, angle=-ch_angle, e1=0.0, e2=-ch_angle, tilt=0.0, fint=0.0, eid='ch1_B3')
    ch1_B4 = SBend(l=Ldip, angle=ch_angle, e1=ch_angle, e2=0.0, tilt=0.0, fint=0.0, eid='ch1_B4')
    ch1_D1a = Drift(l=Ld / 2, eid='ch1_D1a')
    ch1_D1b = Drift(l=Ld / 2, eid='ch1_D1b')
    ch1_D2a = Drift(l=Ld / 2, eid='ch1_D2a')
    ch1_D2b = Drift(l=Ld / 2, eid='ch1_D2b')
    ch1_D3a = Drift(l=Ld / 2, eid='ch1_D3a')
    ch1_D3b = Drift(l=Ld / 2, eid='ch1_D3b')

    # Main chicane lattice
    mchic_start = Marker()
    mchic_dip1 = Marker()
    mchic_mid = Marker()
    mchic_dip3 = Marker()
    mchic_end = Marker()

    chicane = [mchic_start, ch1_B1, ch1_D1a, mchic_dip1, ch1_D1b, ch1_B2,
               ch1_D2a, mchic_mid, ch1_D2b, ch1_B3, ch1_D3a, mchic_dip3, ch1_D3b, ch1_B4, mchic_end]

    # Define object Lattice and compute optics functions
    mend = Marker()
    cell = [Drift(l=0.5)] + chicane + [Drift(l=0.5)] + [mend]

    # Initialize Tracking Method
    method = {'global': SecondTM}
    lat = MagneticLattice(cell, method=method)
    navi = Navigator(lat)

    # Beam diagnostics data
    beam_outdir = os.path.join('diags_1', 'hdf5')
    os.makedirs(beam_outdir, exist_ok=True)
    markers = [mchic_start, mchic_dip1, mchic_mid, mchic_dip3, mchic_end, mend]
    for i, mark in enumerate(markers):
        navi.add_physics_proc(SaveBeamOPMD(folder_path=beam_outdir, iteration=i), mark, mark)

    # TRACK BEAM
    p_array_init = generate_example_beam()
    tws_track, p_array = track(lat, p_array_init, navi)

    # reads the beam back
    p_array_read = read_openpmd_to_parray(beam_outdir,
                                          gamma_ref=p_array.E / m_e_GeV, z_ref=p_array.s)

    # Check that the read beam is equal to the original.
    assert_equal_beams(p_array, p_array_read)


def assert_equal_beams(p_array_1, p_array_2):
    """ Check that two p_array are equal. """

    assert np.allclose(p_array_1.rparticles[0], p_array_2.rparticles[0], rtol=1e-5) and \
        np.allclose(p_array_1.rparticles[1], p_array_2.rparticles[1], rtol=1e-5) and \
        np.allclose(p_array_1.rparticles[2], p_array_2.rparticles[2], rtol=1e-5) and \
        np.allclose(p_array_1.rparticles[3], p_array_2.rparticles[3], rtol=1e-5) and \
        np.allclose(p_array_1.rparticles[4], p_array_2.rparticles[4], rtol=1e-5) and \
        np.allclose(p_array_1.rparticles[5], p_array_2.rparticles[5], rtol=1e-5) and \
        np.allclose(p_array_1.q_array, p_array_2.q_array, rtol=1e-5)


if __name__ == '__main__':
    test_write_parray_to_openpmd()
    test_write_parray_to_openpmd_with_beamline()
