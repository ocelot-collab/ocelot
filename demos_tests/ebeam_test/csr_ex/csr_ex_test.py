"""Test of the demo file demos/ebeam/csr_ex.py"""

import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from demos_tests.params import *
from csr_ex_conf import *


def test_lattice_transfer_map(lattice, p_array, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""

    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)

    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_track_without_csr(lattice, p_array, parametr=None, update_ref_values=False):
    """track function test without CSR"""
    
    navi = Navigator(lattice)
    navi.unit_step = 0.05

    pytest.tws_track_wo, pytest.p_array_wo = track(lattice, p_array, navi)

    pytest.istracked_wo = True

    tws_track = obj2dict(pytest.tws_track_wo)
    p = obj2dict(pytest.p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], TOL, assert_info=' p - ')
    assert check_result(result1+result2)


def test_track_with_csr(lattice, p_array, parametr=None, update_ref_values=False):
    """track function test with CSR"""
    
    csr = CSR()
    csr.traj_step = 0.0002
    csr.apply_step = 0.0005

    navi = Navigator(lattice)
    navi.add_physics_proc(csr, lattice.sequence[0], lattice.sequence[-1])
    navi.unit_step = 0.05

    pytest.tws_track_w, pytest.p_array_w = track(lattice, p_array, navi)

    pytest.istracked_w = True

    tws_track = obj2dict(pytest.tws_track_w)
    p = obj2dict(pytest.p_array_w)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], TOL, assert_info=' p - ')
    assert check_result(result1+result2)


@pytest.mark.parametrize('parametr', [0, 1])
def test_get_current(lattice, p_array, parametr, update_ref_values=False):
    """Get current function test
    :parametr=0 - tracking was done without CSR
    :parametr=1 - tracking was done with CSR
    """

    if parametr == 0:
        if not hasattr(pytest, 'istracked_wo') or not pytest.istracked_wo:
            test_track_without_csr(lattice, p_array, None, True)

        p = pytest.p_array_wo        
    else:
        if not hasattr(pytest, 'istracked_w') or not pytest.istracked_w:
            test_track_with_csr(lattice, p_array, None, True)
        
        p = pytest.p_array_w
    
    sI1, I1 = get_current(p, charge=p.q_array[0], num_bins=200)

    if update_ref_values:
        return numpy2json([sI1, I1])
    
    I_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parametr) +'.json'))
    
    result1 = check_matrix(sI1, I_ref[0], TOL, assert_info=' sI1 - ')
    result2 = check_matrix(I1, I_ref[1], TOL, assert_info=' I1 - ')
    assert check_result(result1+result2)


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
    update_functions.append('test_track_with_csr')
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
