"""Test of the demo file demos/ebeam/ring_orb_correct.py"""

import os
import sys
import time

from ocelot.cpbd.orbit_correction import *
from ocelot.cpbd.response_matrix import *

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from ring_orb_correct_conf import *


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

    beam = Beam()
    beam.E = 2.5
    beam.sigma_E = 0.001
    beam.I = 0.1

    tws = twiss(lattice, Twiss(beam), nPoints=1000)
    
    tws = obj2dict(tws)
    
    if update_ref_values:
        return tws

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws - ')
    assert check_result(result)


def test_responce_matrix(lattice, update_ref_values=False):
    """Responce maxtrix calculation test"""

    orb = NewOrbit(lattice)
    ring_method = RingRM(lattice=orb.lat, hcors=orb.hcors, vcors=orb.vcors, bpms=orb.bpms)

    orb.response_matrix = ResponseMatrix(method=ring_method)
    orb.response_matrix.calculate()

    if update_ref_values:
        return numpy2json(orb.response_matrix.matrix)

    response_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(orb.response_matrix.matrix, response_matrix_ref, TOL, assert_info=' response_matrix - ')
    assert check_result(result)
    

def test_correction(lattice, update_ref_values=False):
    """Orbit correction test"""

    orb = NewOrbit(lattice)
    ring_method = RingRM(lattice=orb.lat, hcors=orb.hcors, vcors=orb.vcors, bpms=orb.bpms)

    orb.response_matrix = ResponseMatrix(method=ring_method)
    orb.response_matrix.calculate()
    
    x_bpm_b, y_bpm_b, x_bpm, y_bpm = correction_wrapper(orb, ring_method)

    if update_ref_values:
        return numpy2json([x_bpm_b, y_bpm_b, x_bpm, y_bpm])

    bpm_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))

    result1 = check_matrix(x_bpm_b, bpm_ref[0], TOL, assert_info=' x_bpm_before_correction - ')
    result2 = check_matrix(y_bpm_b, bpm_ref[1], TOL, assert_info=' y_bpm_before_correction - ')
    result3 = check_matrix(x_bpm, bpm_ref[2], TOL, assert_info=' x_bpm - ')
    result4 = check_matrix(y_bpm, bpm_ref[3], TOL, assert_info=' y_bpm - ')
    assert check_result(result1+result2+result3+result4)


def test_lattice_track(lattice, update_ref_values=False):
    """Lattice track function test"""

    orb = NewOrbit(lattice)
    ring_method = RingRM(lattice=orb.lat, hcors=orb.hcors, vcors=orb.vcors, bpms=orb.bpms)

    orb.response_matrix = ResponseMatrix(method=ring_method)
    orb.response_matrix.calculate()

    correction_wrapper(orb, ring_method)

    p_list = lattice_track(lattice, ring_method.particle0)
    p_list = obj2dict(p_list, unpack=['particle'])

    if update_ref_values:
        return p_list

    p_list_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result = check_dict(p_list, p_list_ref, TOL, assert_info=' p_list - ')
    assert check_result(result)


def correction_wrapper(orb, ring_method):
    
    x_bpm_b, y_bpm_b = ring_method.read_virtual_orbit()
    x_bpm, y_bpm = None, None
     
    if not hasattr(pytest, 'roc_correction'):
        orb.correction(beta=5000)
        x_bpm, y_bpm = ring_method.read_virtual_orbit()
        pytest.roc_correction = True

    return x_bpm_b, y_bpm_b, x_bpm, y_bpm

def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### DBA START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### DBA END ###\n\n\n')
    f.close()
    
    if hasattr(pytest, 'roc_correction'):
        del pytest.roc_correction


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
    update_functions.append('test_responce_matrix')
    update_functions.append('test_correction')
    update_functions.append('test_lattice_track')

    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, True)
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
