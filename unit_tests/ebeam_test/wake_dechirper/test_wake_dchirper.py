#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 11:46:19 2024

@author: weilun
"""
import pytest
import numpy as np
from copy import deepcopy
import sys
import os

from ocelot.cpbd.wake3D import (Wake, 
                                Wake3, 
                                WakeTableDechirperOffAxis, 
                                WakeTableParallelPlate, 
                                WakeTableParallelPlate_origin, 
                                WakeTableParallelPlate3,
                                WakeTableParallelPlate3_origin)
from ocelot import MagneticLattice, SecondTM, track, Navigator, generate_parray
from ocelot import Drift, Marker
from ocelot.gui.accelerator import show_e_beam
from ocelot.cpbd.io import save_particle_array, load_particle_array

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'


def generate_example_beam():
    np.random.seed(10)
    sigma_x = 20e-6
    sigma_xp = 1e-6
    sigma_tau = 7.5e-6
    sigma_delta = 1.5e-4
    p_array = generate_parray(sigma_x=sigma_x, sigma_px=sigma_xp, sigma_y=sigma_x, sigma_py=sigma_xp, 
                             sigma_tau=7.5e-6, sigma_p=sigma_delta, chirp=0., charge=0.25e-09, 
                             nparticles=10000, energy=14, tau_trunc=None, shape='Flat')
    return p_array

def track_with_tilt(p_array_init, model, b, tiltX=0, tiltY=0):

    if model=='02':
        wk_tv_kick = WakeTableParallelPlate_origin(
        b=b,  # distance from the plate in [m]
        a=0.01,  # half gap between plates in [m]
        t=0.25 * 1e-3,  # longitudinal gap in [m]
        p=0.5 * 1e-3,  # period of corrugation in [m]
        length=1,  # length of the corrugated structure in [m]
        sigma=12e-6,  # characteristic (rms) longitudinal beam size in [m]
        orient="horz"  # "horz" or "vert" plate orientation
        )
        wake = Wake()
    elif model == '03':
        wk_tv_kick = WakeTableParallelPlate3_origin(
        b=b,  # distance from the plate in [m]
        a=0.01,  # half gap between plates in [m]
        t=0.25 * 1e-3,  # longitudinal gap in [m]
        p=0.5 * 1e-3,  # period of corrugation in [m]
        length=1,  # length of the corrugated structure in [m]
        sigma=12e-6,  # characteristic (rms) longitudinal beam size in [m]
        orient="horz"  # "horz" or "vert" plate orientation
        )
        wake = Wake3()
    elif model == '12':
        wk_tv_kick = WakeTableParallelPlate(
        b=b,  # distance from the plate in [m]
        a=0.01,  # half gap between plates in [m]
        t=0.25 * 1e-3,  # longitudinal gap in [m]
        p=0.5 * 1e-3,  # period of corrugation in [m]
        length=1,  # length of the corrugated structure in [m]
        sigma=12e-6,  # characteristic (rms) longitudinal beam size in [m]
        orient="horz"  # "horz" or "vert" plate orientation
        )
        wake = Wake()
    elif model == 'N12':
        wk_tv_kick = WakeTableDechirperOffAxis(
        b=b,  # distance from the plate in [m]
        a=0.01,  # half gap between plates in [m]
        t=0.25 * 1e-3,  # longitudinal gap in [m]
        p=0.5 * 1e-3,  # period of corrugation in [m]
        length=1,  # length of the corrugated structure in [m]
        sigma=12e-6,  # characteristic (rms) longitudinal beam size in [m]
        orient="horz"  # "horz" or "vert" plate orientation
        )
        wake = Wake()
    elif model == '13':
        wk_tv_kick = WakeTableParallelPlate3(
        b=b,  # distance from the plate in [m]
        a=0.01,  # half gap between plates in [m]
        t=0.25 * 1e-3,  # longitudinal gap in [m]
        p=0.5 * 1e-3,  # period of corrugation in [m]
        length=1,  # length of the corrugated structure in [m]
        sigma=12e-6,  # characteristic (rms) longitudinal beam size in [m]
        orient="horz"  # "horz" or "vert" plate orientation
        )
        wake = Wake3()
    
    # creation of wake object with parameters
    
    # w_sampling - defines the number of the equidistant sampling points for the one-dimensional
    # wake coefficients in the Taylor expansion of the 3D wake function.
    wake.w_sampling = 500
    wake.wake_table = wk_tv_kick
    wake.step = 1  # step in Navigator.unit_step, dz = Navigator.unit_step * wake.step [m]
    wake.factor = 5
    
    # create a simple lattice MagneticLattice
    m1 = Marker()
    m2 = Marker()
    
    lattice = (m1, Drift(l=1), m2)
    
    lat = MagneticLattice(lattice, method={"global": SecondTM})

    navi = Navigator(lat)
    
    # add physics proccesses
    navi.add_physics_proc(wake, m1, m2)

    p_array = deepcopy(p_array_init)
    p_array.rparticles[0, :] = p_array.x() + p_array.rparticles[4, :] * tiltX
    p_array.rparticles[2, :] = p_array.y() + p_array.rparticles[4, :] * tiltY
    
    p_array_before = deepcopy(p_array)
    
    tws_track, p_array = track(lat, p_array, navi)

    return p_array_before, p_array

def assert_equal_beams(p_array_1, p_array_2):
    """ Check that two p_array are equal. """

    assert np.allclose(p_array_1.rparticles[0], p_array_2.rparticles[0], rtol=1e-5) and \
        np.allclose(p_array_1.rparticles[1], p_array_2.rparticles[1], rtol=1e-5) and \
        np.allclose(p_array_1.rparticles[2], p_array_2.rparticles[2], rtol=1e-5) and \
        np.allclose(p_array_1.rparticles[3], p_array_2.rparticles[3], rtol=1e-5) and \
        np.allclose(p_array_1.rparticles[4], p_array_2.rparticles[4], rtol=1e-5) and \
        np.allclose(p_array_1.rparticles[5], p_array_2.rparticles[5], rtol=1e-5) and \
        np.allclose(p_array_1.q_array, p_array_2.q_array, rtol=1e-5)

def test_wake_02(model='02'):
    p_array_init = generate_example_beam()
    
    b = 500e-6
    
    p1, p2 = track_with_tilt(p_array_init, model, b, tiltY=-15)
    
    #show_e_beam(p2, nparts_in_slice=100)

    p_ref = load_particle_array(REF_RES_DIR + '02.npz')
    assert_equal_beams(p2, p_ref)

def test_wake_03(model='03'):
    p_array_init = generate_example_beam()
    
    b = 500e-6
    
    p1, p2 = track_with_tilt(p_array_init, model, b, tiltY=-15)
    
    #show_e_beam(p2, nparts_in_slice=100)
    
    p_ref = load_particle_array(REF_RES_DIR + '03.npz')
    assert_equal_beams(p2, p_ref)

def test_wake_12(model='12'):
    p_array_init = generate_example_beam()
    
    b = 500e-6
    
    p1, p2 = track_with_tilt(p_array_init, model, b, tiltY=-15)
    
    #show_e_beam(p2, nparts_in_slice=100)

    p_ref = load_particle_array(REF_RES_DIR + '12.npz')
    assert_equal_beams(p2, p_ref)

def test_wake_N12(model='N12'):
    p_array_init = generate_example_beam()
    
    b = 500e-6
    
    p1, p2 = track_with_tilt(p_array_init, model, b, tiltY=-15)
    
    #show_e_beam(p2, nparts_in_slice=100)
 
    p_ref = load_particle_array(REF_RES_DIR + 'N12.npz')
    assert_equal_beams(p2, p_ref)


def test_wake_13(model='13'):
    p_array_init = generate_example_beam()
    
    b = 500e-6
    
    p1, p2 = track_with_tilt(p_array_init, model, b, tiltY=-15)
  
    #show_e_beam(p2, nparts_in_slice=100)
    #save_particle_array('13.npz', p2)
    
    p_ref = load_particle_array(REF_RES_DIR +'13.npz')
    assert_equal_beams(p2, p_ref)


if __name__ == '__main__':
    test_wake_02()
    test_wake_03()
    test_wake_12()
    test_wake_N12()
    test_wake_13()
