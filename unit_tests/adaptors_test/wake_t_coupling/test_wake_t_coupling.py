from time import time
from copy import deepcopy

import pytest
import numpy as np
try:
    from wake_t.driver_witness import LaserPulse
    from wake_t.beamline_elements import (Beamline, PlasmaStage, PlasmaLens,
                                          Drift as WtDrift)
    from wake_t.utilities.bunch_generation import get_gaussian_bunch_from_twiss
    wake_t_installed = True
except ImportError:
    wake_t_installed = False

from ocelot import *
from ocelot.gui.accelerator import *
from ocelot.adaptors.wake_t import wake_t_beam_to_parray, parray_to_wake_t_beam


# Define decorator to skip tests if Wake-T is not installed.
only_if_wake_t_installed = pytest.mark.skipif(
    not wake_t_installed, reason='Wake-T required to run tests')


@only_if_wake_t_installed
def test_conversion():
    """
    Converts a Wake-T beam to Ocelot and back to Wake-T and checks that the
    initial and final beams are the same.

    """
    # Generate Wake-T bunch.
    bunch = get_gaussian_bunch_from_twiss(
        en_x=1e-6, en_y=1e-6, a_x=0, a_y=0, b_x=0.3e-3, b_y=0.3e-3, ene=600,
        ene_sp=0.1, s_t=3, xi_c=20e-6, q_tot=20, n_part=1e4)

    # Convert to Ocelot ParticleArray.
    p_array = wake_t_beam_to_parray(bunch)

    # Convert back to Wake-T.
    bunch_2 = parray_to_wake_t_beam(p_array)

    # Check that the final beam is unchanged by the two-way conversion.
    assert_equal_beams(bunch, bunch_2)


@only_if_wake_t_installed
def test_conversion_with_beamline():
    """
    Tests that the results of a beamline simulation combining Wake-T and
    Ocelot (Plasma stage [Wake-T] -> drift [Ocelot] -> plasma lens [Wake-T])
    are the same than those obtained simulating the same beamline but only
    with Wake-T.

    This therefore checks that the information regarding the propagation
    distance along the beamline and longitudinal position of the particles
    are not affected by the conversion between codes.

    """
    # Generate Wake-T bunch.
    orig_bunch = get_gaussian_bunch_from_twiss(
        en_x=1e-6, en_y=1e-6, a_x=0, a_y=0, b_x=0.3e-3, b_y=0.3e-3, ene=600,
        ene_sp=0.1, s_t=3, xi_c=20e-6, q_tot=20, n_part=1e4)

    # Track beamline only with Wake-T.
    bunch_1 = deepcopy(orig_bunch)
    final_bunch_wt = track_beamline_only_wake_t(bunch_1)

    # Track beamline combining Wake-T and Ocelot.
    bunch_2 = deepcopy(orig_bunch)
    final_bunch_wt_oc = track_beamline_wake_t_and_ocelot(bunch_2)

    # Check that the end results are the same.
    assert_equal_beams(final_bunch_wt, final_bunch_wt_oc)


def track_beamline_only_wake_t(bunch):
    """ Track test beamline only with Wake-T. """
    # Create laser driver.
    laser = LaserPulse(100e-6, l_0=800e-9, w_0=30e-6, a_0=3, tau=30e-15)

    # Create beamline elements.
    plasma = PlasmaStage(
        1e-2, 1e23, laser=laser, wakefield_model='simple_blowout', n_out=10)
    dr = WtDrift(0.1, n_out=10)
    p_lens = PlasmaLens(1e-2, 1000, n_out=5)
    bl = Beamline([plasma, dr, p_lens])

    # Do tracking.
    bl.track(bunch)
    return bunch


def track_beamline_wake_t_and_ocelot(bunch):
    """ Track test beamline combining Wake-T and Ocelot. """
    # 1. Wake-T (simulate plasma cell).
    # ---------------------------------

    # Create laser driver.
    laser = LaserPulse(100e-6, l_0=800e-9, w_0=30e-6, a_0=3, tau=30e-15)

    # Define plasma stage.
    plasma = PlasmaStage(
        1e-2, 1e23, laser=laser, wakefield_model='simple_blowout', n_out=10)

    # Track through plasma stage.
    plasma.track(bunch)

    # 2. Ocelot (simulate drift).
    # ---------------------------

    # Convert beam to particle array.
    p_array_i = wake_t_beam_to_parray(bunch)

    # Define beamline elements.
    d1 = Drift(l=0.1)

    # Initialization of tracking method.
    method = MethodTM()
    method.global_method = SecondTM

    lat = MagneticLattice([d1], method=method)

    navi = Navigator(lat)
    navi.unit_step = 0.01  # m

    # Track lattice.
    p_array_f = deepcopy(p_array_i)
    print("\nTracking with Ocelot ... ")
    start = time()
    tws, p_array_f = track(lat, p_array_f, navi)
    print("\n time exec:", time() - start, "sec")

    # 3. Wake-T (simulate plasma lens).
    # ---------------------------------
    # Convert back to Wake-T.
    bunch = parray_to_wake_t_beam(p_array_f)

    # Define plasma lens.
    p_lens = PlasmaLens(1e-2, 1000, n_out=5)

    # Track through plasma lens.
    p_lens.track(bunch)

    return bunch


def assert_equal_beams(bunch_1, bunch_2):
    """ Check that two Wake-T particle bunches are equal. """
    assert (np.allclose(bunch_1.x, bunch_2.x) and
            np.allclose(bunch_1.px, bunch_2.px) and
            np.allclose(bunch_1.y, bunch_2.y) and
            np.allclose(bunch_1.py, bunch_2.py) and
            np.allclose(bunch_1.xi, bunch_2.xi) and
            np.allclose(bunch_1.pz, bunch_2.pz) and
            np.allclose(bunch_1.q, bunch_2.q) and
            np.allclose(bunch_1.prop_distance, bunch_2.prop_distance))
