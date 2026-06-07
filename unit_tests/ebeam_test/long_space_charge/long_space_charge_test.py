import os
import numpy as np
import pytest
import copy
from long_space_charge_conf import *
from ocelot.cpbd.track import track
from ocelot.cpbd.io import save_particle_array, load_particle_array


REF_DIR = os.path.join(os.path.dirname(__file__), "ref_results")
os.makedirs(REF_DIR, exist_ok=True)


def test_track_with_lsc(lattice, p_array, parameter=None, update_ref_values=False):
    # Track with LSC
    parray_in = p_array.copy()
    lsc = LSC(step=1)
    navi = Navigator(lattice, unit_step=0.1)
    navi.add_physics_proc(lsc, m1, m2)
    _, parray_in = track(lattice, parray_in, navi)
    print(parray_in)

    ref_path = os.path.join(REF_DIR, "test_track_with_lsc.npz")

    if update_ref_values:
        save_particle_array(ref_path, parray_in)
        return

    # Load reference
    parray_ref = load_particle_array(ref_path)
    # Compare mean energy and z-position
    np.testing.assert_allclose(parray_in.tau(), parray_ref.tau(), rtol=1e-5)
    np.testing.assert_allclose(parray_in.p(), parray_ref.p(), rtol=1e-5)


def test_track_wo_lsc(lattice, p_array, parameter=None, update_ref_values=False):
    # Track with LSC
    parray_in = p_array.copy()

    _, parray_in = track(lattice, parray_in)

    ref_path = os.path.join(REF_DIR, "test_track_wo_lsc.npz")

    if update_ref_values:
        save_particle_array(ref_path, parray_in)
        return

    # Load reference
    parray_ref = load_particle_array(ref_path)

    # Compare mean energy and z-position
    np.testing.assert_allclose(parray_in.tau(), parray_ref.tau(), rtol=1e-5)
    np.testing.assert_allclose(parray_in.p(), parray_ref.p(), rtol=1e-5)


def test_K_s_func(lattice, p_array, parameter=None, update_ref_values=False):
    parray_in = p_array.copy().thin_out(10)
    lsc = LSC(step=1)
    navi = Navigator(lattice, unit_step=0.1)
    navi.add_physics_proc(lsc, m1, m2)
    _, parray_out = track(lattice, parray_in, navi)

    s = np.linspace(lsc.s_start, lsc.s_stop, num=10)
    K = lsc.K_s_func(s)
    ref_path = os.path.join(REF_DIR, "test_K_s_func.npz")
    if update_ref_values:
        np.savez(ref_path, s=s, K=K)
        return
    data = np.load(ref_path)
    s_ref = data["s"]
    K_ref = data["K"]
    np.testing.assert_allclose(s, s_ref, rtol=1e-5)
    np.testing.assert_allclose(K, K_ref, rtol=1e-5)

def test_track_with_lsc_K0(lattice, p_array, parameter=None, update_ref_values=False):
    # Track with LSC
    parray_in = p_array.copy()
    lsc = LSC(step=1)
    u.Kx = 0.
    navi = Navigator(lattice, unit_step=0.1)
    navi.add_physics_proc(lsc, m1, m2)
    _, parray_in = track(lattice, parray_in, navi)
    print(parray_in)

    ref_path = os.path.join(REF_DIR, "test_track_with_lsc_K0.npz")

    if update_ref_values:
        save_particle_array(ref_path, parray_in)
        return

    # Load reference
    parray_ref = load_particle_array(ref_path)
    # Compare mean energy and z-position
    np.testing.assert_allclose(parray_in.tau(), parray_ref.tau(), rtol=1e-5)
    np.testing.assert_allclose(parray_in.p(), parray_ref.p(), rtol=1e-5)

@pytest.mark.update
def test_update_ref_values(lattice, p_array, cmdopt):
    update_functions = []
    update_functions.append('test_track_with_lsc')
    update_functions.append('test_track_wo_lsc')
    update_functions.append('test_K_s_func')
    update_functions.append('test_track_with_lsc_K0')
    update_function_parameters = {}

    parameter = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parameter:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p_arr, p, True)

