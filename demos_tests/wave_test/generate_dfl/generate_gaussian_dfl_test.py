"""Test for the function generate_gaussian_dfl(...) from ocelot.optics.wave"""

import os
import sys
import time

from generate_gaussian_dfl_conf import *

import ocelot.optics.wave as wave
from demos_tests.params import *

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'
sys.path.append(os.path.dirname(os.path.dirname(FILE_DIR)))


@pytest.mark.parametrize('parameter', range(4))
def test_generate_gaussian_dfl(args_array_ref, parameter):
    """generate_gaussian_dfl() function testing"""

    dfl_array_ref = np.load(REF_RES_DIR + sys._getframe().f_code.co_name + '.npy')
    (args_ref, kwargs_ref) = args_array_ref[parameter]
    dfl_ref = dfl_array_ref[parameter]
    dfl = wave.generate_gaussian_dfl(*args_ref, **kwargs_ref).fld
    result = check_matrix(dfl, dfl_ref, TOL, tolerance_type='relative',
                          assert_info=f'- generate_gaussian_dfl(){parameter}')
    assert check_result(result)


def setup_module(module):
    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### GENERATE_GAUSSIAN_DFL START ###\n\n')
    f.close()


def teardown_module(module):
    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### GENERATE_GAUSSIAN_DFL END ###\n\n\n')
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

# @pytest.mark.update
# def test_update_ref_values(args_array_ref, cmdopt):
#     update_functions = []
#     update_functions.append('test_generate_gaussian_dfl')
#
#     update_function_parameters = {}
#     update_function_parameters['test_get_current'] = []
#
#     parametr = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']
#
#     if cmdopt in update_functions:
#         for p in parametr:
#             p_arr = copy.deepcopy(p_array)
#             result = eval(cmdopt)(lattice, p_arr, p, True)
#
#             if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.npy'):
#                 os.rename(REF_RES_DIR + cmdopt + str(p) + '.npy', REF_RES_DIR + cmdopt + str(p) + '.old')
#
#             json_save(result, REF_RES_DIR + cmdopt + str(p) + '.npy')
