"""Test of the demo file demos/ebeam/csr_ex.py"""

import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from acc_utils_conf import *


def test_lattice_transfer_map(lattice, p_array, parameter=None, update_ref_values=False):
    """R matrix calculation test"""

    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)

    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)

@pytest.mark.parametrize('parameter', [0, 1])
def test_lattice_transfer_map_RT(lattice, p_array, parameter, update_ref_values=False):
    """test R56 and T566 of the chicane"""

    r56, t566, u5666, Sref = chicane_RTU(yoke_len=b1.l/b1.angle*np.sin(b1.angle), dip_dist=d1.l*np.cos(b1.angle), r=b1.l/b1.angle, type='c')
    lattice = copy.deepcopy(lattice)

    if parameter == 1:
        for elem in lattice.sequence:
            if elem.__class__ == Bend:
                elem.tilt = np.pi / 2

        lattice.update_transfer_maps()

    r_matrix = lattice_transfer_map(lattice, 0.0)
    result1 = check_value(r_matrix[4, 5], r56, tolerance=1.0e-14, assert_info=" R56 ")
    result2 = check_value(lattice.T[4, 5, 5], t566, tolerance=1.0e-14, assert_info=" T566 ")

    assert check_result([result1, result2])


def test_rf2beam(lattice, p_array, parameter=None, update_ref_values=False):
    """
    track function test without CSR

    0 - normal tracking
    1 - tilt bending magnets and tilt back electron beam then untilt beam and compare with ref beam (twiss not checked)
    """
    v1 = 0.14746291505994155
    phi1 = -11.105280079934298
    vh = 0.030763428944485114
    phih = 132.9179951484828 - 360
    E1, chirp, curvature, skewness = rf2beam(v1, phi1, vh, phih, n=3, freq=1.3e9, E0=0.00675, zeta1=0., zeta2=0.,
                                 zeta3=0.)

    v1_r, phi1_r, vh_r, phih_r = beam2rf(E1, chirp, curvature, skewness, n=3, freq=1.3e9, E0=0.00675, zeta1=0., zeta2=0.,
                                 zeta3=0.)

    r1 = check_value(v1_r, v1, tolerance=1.0e-8, tolerance_type='relative', assert_info='v1')
    r2 = check_value(phi1_r, phi1, tolerance=1.0e-8, tolerance_type='relative', assert_info='phi1')
    r3 = check_value(vh_r, vh, tolerance=1.0e-8, tolerance_type='relative', assert_info='vh')
    r4 = check_value(phih_r, phih, tolerance=1.0e-8, tolerance_type='relative', assert_info='phih')

    assert check_result([r1, r2, r3, r4])



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

    update_function_parameters = {}

    parameter = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parameter:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p_arr, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
