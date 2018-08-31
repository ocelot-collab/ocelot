"""Test of the demo file demos/sr/spatial.py"""

import os
import sys
import time
import copy

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from demos_tests.params import *
from rad_beam_conf import *


@pytest.mark.parametrize('parameter', [0, 1, 2])
def test_calculate_radiation(lattice, screen, beam, parameter, update_ref_values=False):
    """calculate_radiation fucntion test"""
    screen.nullify()
    accuracy = 1
    beam = copy.deepcopy(beam)
    if parameter == 0:
        p_array = beam

    elif parameter == 1:
        p_array = beam
        p_array.tau()[:] *= -1

    elif parameter == 2:
        p_array = beam
        p_array.tau()[:] = 0

    screen = calculate_radiation(lattice, screen, p_array, accuracy=accuracy)
    print(screen.Total[:2])
    if update_ref_values:
        return {'Eph':screen.Eph.tolist(), 'Yph':screen.Yph.tolist(), 'Xph':screen.Xph.tolist(), 'Total':screen.Total.tolist(), 'Sigma':screen.Sigma.tolist(), 'Pi':screen.Pi.tolist()}
    
    screen_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')

    result1 = check_matrix(screen.Eph, screen_ref['Eph'], TOL, assert_info=' Eph - ')
    result2 = check_matrix(screen.Yph, screen_ref['Yph'], TOL, assert_info=' Yph - ')
    result3 = check_matrix(screen.Xph, screen_ref['Xph'], TOL, assert_info=' Xph - ')
    result4 = check_matrix(screen.Total, screen_ref['Total'], TOL, assert_info=' Total - ')
    result5 = check_matrix(screen.Sigma, screen_ref['Sigma'], TOL, assert_info=' Sigma - ')
    result6 = check_matrix(screen.Pi, screen_ref['Pi'], TOL, assert_info=' Pi - ')
    assert check_result(result1 + result2 + result3 + result4 + result5 + result6)


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### SPATIAL START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### SPATIAL END ###\n\n\n')
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
def test_update_ref_values(lattice, screen, beam, cmdopt):
    
    update_functions = []
    update_functions.append('test_calculate_radiation')
    
    update_function_parameters = {}
    update_function_parameters['test_calculate_radiation'] = [0, 1, 2]
    
    parameter = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parameter:
            result = eval(cmdopt)(lattice, screen, beam, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
