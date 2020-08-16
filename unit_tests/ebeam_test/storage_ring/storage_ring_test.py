"""Test of the demo file demos/ebeam/storage_ring.py"""

import os
import sys
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from storage_ring_conf import *


def test_lattice_transfer_map(lattice, beam=None, update_ref_values=False):
    """R maxtrix test"""

    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_twiss(lattice, beam, update_ref_values=False):
    """Twiss parameters calculation function test"""

    tws = twiss(lattice, Twiss(beam), nPoints=1000)

    tws = obj2dict(tws)
    
    if update_ref_values:
        return tws

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws - ')
    assert check_result(result)


def test_match_tunes(lattice, beam, update_ref_values=False):
    """Twiss parameters calculation function test"""

    match_tolerance = 1.0e-5

    nu_x, nu_y, ncells = match_tunes_wrapper(lattice, beam)
    tws = twiss(lattice, Twiss(beam), nPoints=1000)

    new_nu_x = tws[-1].mux / 2.0 / pi * ncells
    new_nu_y = tws[-1].muy / 2.0 / pi * ncells

    result1 = check_value(new_nu_x, nu_x, match_tolerance, assert_info=' nu_x - \n')
    result2 = check_value(new_nu_y, nu_y, match_tolerance, assert_info=' nu_y - \n')
    assert check_result([result1, result2])


def test_lattice_transfer_map_after_matching(lattice, beam, update_ref_values=False):
    """R maxtrix after matching test"""

    nu_x, nu_y, ncells = match_tunes_wrapper(lattice, beam)
        
    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix after matching - ')
    assert check_result(result)


def test_twiss_after_matching(lattice, beam, update_ref_values=False):
    """Twiss parameters after matching calculation function test"""

    nu_x, nu_y, ncells = match_tunes_wrapper(lattice, beam)

    tws = twiss(lattice, Twiss(beam), nPoints=1000)

    tws = obj2dict(tws)
    
    if update_ref_values:
        return tws

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws after matching - ')
    assert check_result(result)


def test_e_beam_params(lattice, beam, update_ref_values=False):
    """Ebeam parameters calculation function test"""

    nu_x, nu_y, ncells = match_tunes_wrapper(lattice, beam)
    tws = Twiss(beam)
    ebp = EbeamParams(lattice, tws0=tws, nsuperperiod=8)

    ebp = obj2dict([ebp])
     
    if update_ref_values:
        return ebp

    ebp_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result = check_dict(ebp, ebp_ref, TOL, 'absotute', assert_info=' EbeamParams - ')
    assert check_result(result)


def match_tunes_wrapper(lattice, beam):

    nu_x = 1.2
    nu_y = 0.91
    ncells = 1

    if not hasattr(pytest, 'sr_match_tunes'):
        quads = [Q1, Q2, Q3, Q4]
        tws = twiss(lattice, Twiss(beam), nPoints=1000)
        match_tunes(lattice, tws[-1], quads, nu_x, nu_y, ncells=ncells, print_proc=0)
        pytest.sr_match_tunes = True
    
    return nu_x, nu_y, ncells


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### STORAGE_RING START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### STORAGE_RING END ###\n\n\n')
    f.close()

    if hasattr(pytest, 'sr_match_tunes'):
        del pytest.sr_match_tunes


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
def test_update_ref_values(lattice, beam, cmdopt):
    
    update_functions = []
    update_functions.append('test_lattice_transfer_map')
    update_functions.append('test_twiss')
    update_functions.append('test_lattice_transfer_map_after_matching')
    update_functions.append('test_twiss_after_matching')
    update_functions.append('test_e_beam_params')
    
    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, beam, True)
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
