"""Test of the demo file demos/ebeam/csr_ex.py"""

import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from bump_utils_conf import *


def test_lattice_transfer_map(lattice, parameter=None, update_ref_values=False):
    """R matrix calculation test"""

    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)

    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_bump(lattice, parameter=None, update_ref_values=False):
    """test R56 and T566 of the chicane"""

    cor_list = [c1, c2, c3, c4]
    a = bump_4cors(lattice, cor_list, marker=m, x=0.001, xp=-0.00, energy=1)
    print("corrector, strength: ", a * 1000, " mrad")

    plist = lattice_track(lattice, Particle(E=1))

    plist = obj2dict(plist)

    if update_ref_values:
        return plist

    plist_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result = check_dict(plist, plist_ref, TOL, assert_info=' plist - ')

    assert check_result(result)


def test_bump_disp(lattice, parameter=None, update_ref_values=False):

    cor_list = [c1, c2, c3, c4]

    a = bump_4cors(lattice, cor_list, marker=m, x=0.001, xp=-0.00, energy=1)

    lat_new = convert_cors2dipoles(lattice, cor_list, energy=14)
    tws0 = Twiss()
    tws0.beta_x = 10
    tws0.beta_y = 10

    tws = twiss(lat_new, tws0=tws0)
    tws = obj2dict(tws)

    if update_ref_values:
        return tws

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws - ')

    assert check_result(result)



def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### BUMP_UTILS START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### BUMP_UTILS END ###\n\n\n')
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
def test_update_ref_values(lattice, cmdopt):
    
    update_functions = []
    update_functions.append('test_lattice_transfer_map')
    update_functions.append('test_bump')
    update_functions.append('test_bump_disp')
    update_function_parameters = {}

    parameter = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parameter:
            result = eval(cmdopt)(lattice, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
