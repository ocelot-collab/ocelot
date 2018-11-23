"""Test of lattice save function in cpbd/io.py file"""

import os
import sys
import time

from ocelot.cpbd.io import *

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from demos_tests.params import *
from io_lattice_conf import *


def test_original_lattice_transfer_map(lattice, tws0, method, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""

    r_matrix = lattice_transfer_map(lattice, tws0.E)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_original_twiss(lattice, tws0, method, parametr=None, update_ref_values=False):
    """Twiss parameters calculation function test"""

    tws = twiss(lattice, tws0, nPoints=None)
    
    tws = obj2dict(tws)
    
    if update_ref_values:
        return tws

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws - ')
    assert check_result(result)


@pytest.mark.parametrize('parametr', [False, True])
def test_lat2input(lattice, tws0, method, parametr, update_ref_values=False):
    """lat2input with tws0 saving function test"""

    lines_arr = lat2input(lattice, tws0=tws0)
    lines = ''.join(lines_arr)
    
    loc_dict = {}
    try:
        exec(lines, globals(), loc_dict)
    except Exception as err:
        assert check_result(['Exception error during the lattice file execution, parametr is ' + str(parametr)])

    if parametr:
        if "tws0" in loc_dict:
            tws0_new = loc_dict['tws0']
        else:
            assert check_result(['No tws0 in the lattice file, parametr is ' + str(parametr)])
    else:
        tws0_new = tws0

    if 'cell' in loc_dict:
        lattice_new = MagneticLattice(loc_dict['cell'], method=method)

        lattice_new_transfer_map_check(lattice_new, tws0_new)
        twiss_new_check(lattice_new, tws0_new)
    else:
        assert check_result(['No cell variable in the lattice file, parametr is ' + str(parametr)])
    
        
def lattice_new_transfer_map_check(lattice, tws0):

    r_matrix = lattice_transfer_map(lattice, tws0.E)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + 'test_original_lattice_transfer_map.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix for new lattice - ')
    assert check_result(result)


def twiss_new_check(lattice, tws0):

    tws = twiss(lattice, tws0, nPoints=None)
    
    tws = obj2dict(tws)

    tws_ref = json_read(REF_RES_DIR + 'test_original_twiss.json')
    
    result = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws for new lattice - ')
    assert check_result(result)
    

def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### IO Lattice START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### IO Lattice END ###\n\n\n')
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
def test_update_ref_values(lattice, tws0, method, cmdopt):
    
    update_functions = []
    update_functions.append('test_original_lattice_transfer_map')
    update_functions.append('test_original_twiss')
    
    # function test_lat2input function need not be added here.
    # It is used reference results from test_original_lattice_transfer_map and test_original_twiss functions

    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, tws0, method, None, True)
        if result is None:
            return
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
