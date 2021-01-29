"""Test of the demo file demos/ebeam/csr_ex.py"""

import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from wake_conf import *



def test_get_long_wake(lattice, p_array, parameter=None, update_ref_values=False):
    """ func generate_parray testing """

    wt = WakeTable(wake_file="./unit_tests/ebeam_test/wake/wake_table.dat")
    w = Wake()
    w.wake_table = wt
    w.prepare(None)

    x = np.arange(-30, 30, 0.01) * 1e-6
    sigma = 7e-6
    y2_x = lambda x: 5000 if np.abs(x) < 10e-6 else 0.
    y = np.array([y2_x(xi) for xi in x])


    profile = np.hstack((x.reshape(-1, 1), y.reshape(-1, 1)))
    x, Wz = w.get_long_wake(profile)

    if update_ref_values:
        return {'x': list(x), 'Wz': list(Wz)}

    current_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    result1 = check_matrix(x, current_ref['x'], TOL, assert_info=' x - ')
    result2 = check_matrix(Wz, current_ref['Wz'], TOL, assert_info=' Wz - ')
    assert check_result(result1 + result2)


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### PHYS PROC START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### PHYS PROC END ###\n\n\n')
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
    update_functions.append('test_get_long_wake')

    update_function_parameters = {}


    parametr = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parametr:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p_arr, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
