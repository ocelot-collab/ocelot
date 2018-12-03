"""Test of the demo file demos/ebeam/csr_ex.py"""

import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from csr_ex_conf import *


def test_lattice_transfer_map(lattice, p_array, parameter=None, update_ref_values=False):
    """R matrix calculation test"""

    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)

    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)

@pytest.mark.parametrize('parameter', [0, 1])
def test_lattice_transfer_map_RT(lattice, p_array, parameter, update_ref_values=False):
    """test R56 and T566 of the chicane"""

    r56, t566, u5666, Sref = chicane_RTU(yoke_len=b1.l/b1.angle*np.sin(b1.angle), dip_dist=d1.l*np.cos(b1.angle), r=b1.l/b1.angle, type='c')
    lattice = copy.deepcopy(lattice)

    if parameter == 1:
        for elem in lattice.sequence:
            if elem.__class__ == Bend:
                elem.tilt = np.pi / 2

        lattice.update_transfer_maps()

    r_matrix = lattice_transfer_map(lattice, 0.0)
    result1 = check_value(r_matrix[4, 5], r56, tolerance=1.0e-14, assert_info=" R56 ")
    result2 = check_value(lattice.T[4, 5, 5], t566, tolerance=1.0e-14, assert_info=" T566 ")

    assert check_result([result1, result2])


@pytest.mark.parametrize('parameter', [0, 1])
def test_track_without_csr(lattice, p_array, parameter, update_ref_values=False):
    """
    track function test without CSR

    0 - normal tracking
    1 - tilt bending magnets and tilt back electron beam then untilt beam and compare with ref beam (twiss not checked)
    """

    tilt = 0
    if parameter == 1:
        tilt = np.pi/2.

    for elem in lattice.sequence:
        if elem.__class__ == Bend:
            elem.tilt = tilt

    lattice.update_transfer_maps()

    navi = Navigator(lattice)
    navi.unit_step = 0.05

    p_array.rparticles[:] = np.dot(rot_mtx(-tilt), p_array.rparticles)[:]
    pytest.tws_track_wo, pytest.p_array_wo = track(lattice, p_array, navi)
    p_array.rparticles[:] = np.dot(rot_mtx(tilt), p_array.rparticles)[:]

    pytest.istracked_wo = True

    tws_track = obj2dict(pytest.tws_track_wo)
    p = obj2dict(pytest.p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')


    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    if parameter == 1:
        result1 = [None]
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result1+result2)


def test_track_without_csr_tilted(lattice, p_array, parameter=None, update_ref_values=False):
    """track function test without CSR in vertical plane """
    lattice = copy.deepcopy(lattice)


    for elem in lattice.sequence:
        if elem.__class__ == Bend:
            elem.tilt = np.pi / 2

    lattice.update_transfer_maps()
    navi = Navigator(lattice)
    navi.unit_step = 0.05
    tws_track_wo, p_array_wo = track(lattice, p_array, navi)


    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], TOL, assert_info=' p - ')
    assert check_result(result1 + result2)


def test_track_without_csr_rotated(lattice, p_array, parameter=None, update_ref_values=False):
    """track function test without CSR, comparison between the rotated beam and track in vertical plane"""
    lattice_copy = copy.deepcopy(lattice)
    p_array_copy = copy.deepcopy(p_array)

    navi = Navigator(lattice_copy)
    navi.unit_step = 0.05

    tilt = np.pi / 2
    for elem in lattice_copy.sequence:
        if elem.__class__ == Bend:
            elem.tilt = tilt

    lattice_copy.update_transfer_maps()

    tws_track_tilted, p_array_wo_tilted = track(lattice_copy, p_array_copy, navi)

    tilt = np.pi / 2
    for elem in lattice.sequence:
        if elem.__class__ == Bend:
            elem.tilt = 0

    lattice.update_transfer_maps()

    navi = Navigator(lattice)
    navi.unit_step = 0.05

    p_array.rparticles[:] = np.dot(rot_mtx(tilt), p_array.rparticles)[:]
    tws_track_rot, p_array = track(lattice, p_array, navi)
    p_array.rparticles[:] = np.dot(rot_mtx(-tilt), p_array.rparticles)[:]

    p1 = obj2dict(p_array_wo_tilted)

    p2 = obj2dict(p_array)
    result2 = check_dict(p1, p2, TOL, assert_info=' p - ')
    assert check_result(result2)


@pytest.mark.parametrize('parameter', [0, 1])
def test_track_with_csr(lattice, p_array, parameter, update_ref_values=False):
    """
    track function test with CSR

    0 - normal tracking
    1 - tilt bending magnets and tilt back electron beam then untilt beam and compare with ref beam (twiss not checked)
    """

    tilt = 0
    if parameter == 1:
        tilt = np.pi/2.

    for elem in lattice.sequence:
        if elem.__class__ == Bend:
            elem.tilt = tilt

    lattice.update_transfer_maps()
    
    csr = CSR()
    csr.traj_step = 0.0002
    csr.apply_step = 0.0005

    navi = Navigator(lattice)
    navi.add_physics_proc(csr, lattice.sequence[0], lattice.sequence[-1])
    navi.unit_step = 0.05

    p_array.rparticles[:] = np.dot(rot_mtx(-tilt), p_array.rparticles)[:]
    pytest.tws_track_w, pytest.p_array_w = track(lattice, p_array, navi)
    p_array.rparticles[:] = np.dot(rot_mtx(tilt), p_array.rparticles)[:]

    pytest.istracked_w = True

    tws_track = obj2dict(pytest.tws_track_w)
    p = obj2dict(pytest.p_array_w)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    if parameter == 1:
        result1 = [None]
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], TOL, assert_info=' p_array - ')
    assert check_result(result1+result2)


@pytest.mark.parametrize('parameter', [0, 1])
def test_get_current(lattice, p_array, parameter, update_ref_values=False):
    """Get current function test
    :parametr=0 - tracking was done without CSR
    :parametr=1 - tracking was done with CSR
    """

    if parameter == 0:
        if not hasattr(pytest, 'istracked_wo') or not pytest.istracked_wo:
            test_track_without_csr(lattice, p_array, 0, True)

        p = pytest.p_array_wo        
    else:
        if not hasattr(pytest, 'istracked_w') or not pytest.istracked_w:
            test_track_with_csr(lattice, p_array, 0, True)
        
        p = pytest.p_array_w
    
    sI1, I1 = get_current(p, charge=p.q_array[0], num_bins=200)

    if update_ref_values:
        return numpy2json([sI1, I1])
    
    I_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) +'.json'))
    
    result1 = check_matrix(sI1, I_ref[0], TOL, assert_info=' sI1 - ')
    result2 = check_matrix(I1, I_ref[1], TOL, assert_info=' I1 - ')
    assert check_result(result1+result2)


def test_track_undulator_with_csr(lattice, p_array, parameter=None, update_ref_values=False):
    """track function test without CSR in vertical plane """

    d1 = Drift(l=0.1)
    d2 = Drift(l=1)

    und = Undulator(lperiod=0.4, nperiods=9, Kx=44.81)

    m = MethodTM()
    m.global_method = SecondTM

    lat = MagneticLattice((d1, und, d2), method=m)

    navi = Navigator(lat)
    navi.unit_step = 0.05

    csr = CSR()
    csr.energy = p_array.E
    navi.add_physics_proc(csr, lat.sequence[0], lat.sequence[-1])

    tws_track_wo, p_array_wo = track(lat, p_array, navi)


    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], TOL, assert_info=' p - ')
    assert check_result(result1 + result2)


def test_csr_arcline_rk_traj(lattice, p_array, parameter=None, update_ref_values=False):
    """track function test without CSR, comparison between the rotated beam and track in vertical plane"""
    lattice_copy = copy.deepcopy(lattice)
    p_array_copy = copy.deepcopy(p_array)

    navi = Navigator(lattice_copy)
    navi.unit_step = 0.05

    csr = CSR()
    csr.energy = p_array_copy.E
    navi.add_physics_proc(csr, lattice_copy.sequence[0], lattice_copy.sequence[-1])

    tws_track_arc, p_array_arc = track(lattice_copy, p_array_copy, navi)

    lattice_copy = copy.deepcopy(lattice)
    p_array_copy = copy.deepcopy(p_array)

    navi = Navigator(lattice_copy)
    navi.unit_step = 0.05

    csr = CSR()
    csr.energy = 40
    csr.rk_traj = True
    navi.add_physics_proc(csr, lattice_copy.sequence[0], lattice_copy.sequence[-1])

    tws_track_rk, p_array_rk = track(lattice_copy, p_array_copy, navi)


    p1 = obj2dict(p_array_arc)

    p2 = obj2dict(p_array_rk)

    result2 = check_dict(p1, p2, tolerance=1.0e-8,tolerance_type='absolute', assert_info=' p - ')
    assert check_result(result2)


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### CSR_EX START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### CSR_EX END ###\n\n\n')
    f.close()


def setup_function(function):
    
    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write(function.__name__)
    f.close()

    pytest.t_start = time.time()


def teardown_function(function):
    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write(' execution time is ' + '{:.3f}'.format(time.time() - pytest.t_start) + ' sec\n\n')
    f.close()
    

@pytest.mark.update
def test_update_ref_values(lattice, p_array, cmdopt):
    
    update_functions = []
    update_functions.append('test_lattice_transfer_map')
    update_functions.append('test_track_without_csr')
    update_functions.append('test_track_without_csr_tilted')
    update_functions.append('test_track_with_csr')

    update_functions.append('test_get_current')
    update_function_parameters = {}
    update_function_parameters['test_get_current'] = [0, 1]

    update_functions.append('test_track_undulator_with_csr')
    update_functions.append("test_csr_arcline_rk_traj")

    parametr = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parametr:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p_arr, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
