"""Test of the demo file demos/ebeam/csr_ex.py"""

import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from rk_track_conf import *


def test_generate_parray(lattice, p_array, parameter=None, update_ref_values=False):
    """ func generate_parray testing """


    p = obj2dict(p_array)
    
    if update_ref_values:
        return {'p_array': p}

    p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    result = check_dict(p, p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result)


@pytest.mark.parametrize('parameter', [0, 1, 2])
def test_track_undulator_with_diff_chirp(lattice, p_array, parameter, update_ref_values=False):
    """
    test Runge_Kutta transfer map for undulator

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    2 - tracking of the electron beam with zero energy chirp trough undulator
    """

    p_array_track = copy.deepcopy(p_array)
    if parameter == 1:
        p_array_track.p()[:] *= -1

    else:
        p_array_track.p()[:] = 0

    navi = Navigator(lattice)

    tws_track_wo, p_array_after = track(lattice, p_array_track, navi, calc_tws=False)

    p = obj2dict(p_array_after)

    if update_ref_values:
        return {'p_array': p}

    p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')


    result = check_dict(p, p_array_ref['p_array'], tolerance=1.0e-14, tolerance_type='absolute', assert_info=' p - ')
    assert check_result(result)


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### RK TRACK START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### RK TRACK END ###\n\n\n')
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
    update_functions.append('test_generate_parray')
    update_functions.append('test_track_undulator_with_diff_chirp')

    update_function_parameters = {}
    update_function_parameters['test_track_undulator_with_diff_chirp'] = [0, 1, 2]
    
    parameter = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']
    if cmdopt in update_functions:
        for p in parameter:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p_arr, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
