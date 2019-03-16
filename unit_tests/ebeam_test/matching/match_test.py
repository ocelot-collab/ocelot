"""Test of the demo file demos/ebeam/dba.py"""

import os
import sys
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from match_conf import *


def test_lattice_transfer_map(lattice, update_ref_values=False):
    """R maxtrix calculation test"""

    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_twiss(lattice, update_ref_values=False):
    """Twiss parameters calculation function test"""

    tws0 = Twiss()
    tws0.beta_x = 8.4
    tws0.beta_y = 8.4
    tws0.alpha_x = -55.8
    tws0.alpha_y = -55.8
    tws0.E = 0.005071

    tws = twiss(lattice, tws0, nPoints=100)
    
    tws = obj2dict(tws)
    
    if update_ref_values:
        return tws

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws - ')
    assert check_result(result)


def test_solenoid_match(lattice, update_ref_values=False):
    """After matching R maxtrix calculcation test"""
    tws0 = Twiss()
    tws0.beta_x = 8.4
    tws0.beta_y = 8.4
    tws0.alpha_x = -55.8
    tws0.alpha_y = -55.8
    tws0.E = 0.005071

    # match beta functions at the end
    constr = {end: {"beta_x": 10, "beta_y": 10}}
    vars = [sol1]


    match(lattice, constr, vars, tws0, verbose=False)

    tws = twiss(lattice, Twiss(), nPoints=100)

    tws = obj2dict(tws)

    if update_ref_values:
        return tws

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws after matching - ')
    assert check_result(result)



def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### DBA START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### DBA END ###\n\n\n')
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
    update_functions.append('test_twiss')
    update_functions.append('test_solenoid_match')

    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, True)
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
