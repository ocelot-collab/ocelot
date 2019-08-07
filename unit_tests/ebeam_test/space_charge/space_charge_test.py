"""Test of the demo file demos/ebeam/space_charge.py"""

import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from space_charge_conf import *


def test_track_without_sp(lattice, p_array, parameter=None, update_ref_values=False):
    """track function test without space charge"""

    tws_track, p = track_wrapper(lattice, p_array, 0)

    tws_track = obj2dict(tws_track)
    p = obj2dict(p)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, tolerance_type='relative', assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], TOL, assert_info=' p_array - ')
    assert check_result(result1+result2)

def test_track_without_sp_bounds(lattice, p_array, parameter=None, update_ref_values=False):
    """track function test without space charge"""
    navi = Navigator(lattice)
    navi.unit_step = 0.02

    tws_track, p = track(lattice, p_array, navi, bounds=[-0.5, 0.5])

    tws_track = obj2dict(tws_track)
    p = obj2dict(p)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, tolerance_type='relative', assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], TOL, assert_info=' p_array - ')
    assert check_result(result1+result2)

def test_track_with_sp(lattice, p_array, parameter=None, update_ref_values=False):
    """track function test with space charge"""

    tws_track, p = track_wrapper(lattice, p_array, 1)

    tws_track = obj2dict(tws_track)
    p = obj2dict(p)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=1e-12, tolerance_type='absolute', assert_info=' p_array - ')
    assert check_result(result1+result2)


def test_track_with_sp_bounds(lattice, p_array, parameter=None, update_ref_values=False):
    """track function test with space charge"""
    sc1 = SpaceCharge()
    sc1.nmesh_xyz = [63, 63, 63]
    sc1.step = 1

    sc5 = SpaceCharge()
    sc5.nmesh_xyz = [63, 63, 63]
    sc5.step = 5

    navi = Navigator(lattice)
    navi.add_physics_proc(sc1, lattice.sequence[0], C_A1_1_2_I1)
    navi.add_physics_proc(sc5, C_A1_1_2_I1, lattice.sequence[-1])
    navi.unit_step = 0.02
    tws_track, p = track(lattice, p_array, navi, bounds=[-0.5, 0.5])

    tws_track = obj2dict(tws_track)
    p = obj2dict(p)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], tolerance=1e-8, tolerance_type='absolute', assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=1e-12, tolerance_type='absolute',
                         assert_info=' p_array - ')
    assert check_result(result1 + result2)


def test_track_with_lsc(lattice, p_array, parameter=None, update_ref_values=False):
    """track function test with space charge"""
    sc1 = LSC()
    sc1.step = 1

    sc5 = LSC()
    sc5.step = 5

    navi = Navigator(lattice)
    navi.add_physics_proc(sc1, lattice.sequence[0], C_A1_1_2_I1)
    navi.add_physics_proc(sc5, C_A1_1_2_I1, lattice.sequence[-1])
    navi.unit_step = 0.02
    tws_track, p = track(lattice, p_array, navi)

    tws_track = obj2dict(tws_track)
    p = obj2dict(p)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], tolerance=1e-11, tolerance_type='absolute',
                         assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=1e-12, tolerance_type='absolute',
                         assert_info=' p_array - ')
    assert check_result(result1 + result2)

@pytest.mark.parametrize('parameter', [0, 1])
def test_get_current(lattice, p_array, parameter, update_ref_values=False):
    """Get current function test
    :parametr=0 - tracking was done without space charge
    :parametr=1 - tracking was done with space charge
    """

    tws_track, p = track_wrapper(lattice, p_array, parameter)
    
    sI1, I1 = get_current(p, charge=p.q_array[0], num_bins=200)

    if update_ref_values:
        return numpy2json([sI1, I1])
    
    I_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) +'.json'))

    result1 = check_matrix(sI1, I_ref[0], TOL, assert_info=' sI1 - ')
    result2 = check_matrix(I1, I_ref[1], TOL, assert_info=' I1 - ')
    assert check_result(result1+result2)


def track_wrapper(lattice, p_array, param, bounds=None):

    if not hasattr(pytest, 'sp_track_list') or not hasattr(pytest, 'sp_p_array'):
        pytest.sp_track_list = [None, None]
        pytest.sp_p_array = [None, None]

    if pytest.sp_track_list[param] is None or pytest.sp_p_array[param] is None:
        if param == 0:
            navi = Navigator(lattice)
            navi.unit_step = 0.02

        elif param == 1:
            sc1 = SpaceCharge()
            sc1.nmesh_xyz = [63, 63, 63]
            sc1.step = 1

            sc5 = SpaceCharge()
            sc5.nmesh_xyz = [63, 63, 63]
            sc5.step = 5

            navi = Navigator(lattice)
            navi.add_physics_proc(sc1, lattice.sequence[0], C_A1_1_2_I1)
            navi.add_physics_proc(sc5, C_A1_1_2_I1, lattice.sequence[-1])
            navi.unit_step = 0.02

        pytest.sp_track_list[param], pytest.sp_p_array[param] = track(lattice, p_array, navi, bounds=bounds)

    return pytest.sp_track_list[param], pytest.sp_p_array[param]


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### SPACE_CHARGE START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### SPACE_CHARGE END ###\n\n\n')
    f.close()

    if hasattr(pytest, 'sp_track_list'):
        del pytest.sp_track_list

    if hasattr(pytest, 'sp_p_array'):
        del pytest.sp_p_array


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
    update_functions.append('test_track_without_sp')
    update_functions.append("test_track_without_sp_bounds")
    update_functions.append('test_track_with_sp')
    update_functions.append('test_track_with_sp_bounds')
    update_functions.append('test_track_with_lsc')
    update_functions.append('test_get_current')
    
    update_function_parameters = {}
    update_function_parameters['test_get_current'] = [0, 1]
    
    parametr = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parametr:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p_arr, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
