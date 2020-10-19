"""Test of the demo file demos/ebeam/rf_twiss.py"""

import os
import sys
import time
import copy

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from rf_twiss_conf import *


def test_lattice_transfer_map(lattice, p_array, parameter=None, update_ref_values=False):
    """R maxtrix test"""

    beam = Beam()
    beam.E = 2.4

    r_matrix = lattice_transfer_map(lattice, beam.E)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)

@pytest.mark.parametrize('parameter', [0, 1])
def test_lattice_r_t_maps_wo_coupler(lattice, p_array, parameter, update_ref_values=False):
    """R maxtrix test"""
    if parameter == 0:
        for elem in lattice.sequence:
            if elem.__class__ is Cavity:
                elem.vx_up = -5.6813e-5 + 1.0751e-5j
                elem.vy_up = -4.1091e-5 + 5.739e-7j
                elem.vxx_up = 0.00099943 - 0.00081401j
                elem.vxy_up = 0.0034065 - 0.0004146j
                elem.vx_down = -2.4014e-5 + 1.2492e-5j
                elem.vy_down = 3.6481e-5 + 7.9888e-6j
                elem.vxx_down = -0.004057 - 0.0001369j
                elem.vxy_down = 0.0029243 - 1.2891e-5j
    else:
        for elem in lattice.sequence:
            if elem.__class__ is Cavity:
                elem.vx_up = 0
                elem.vy_up = 0
                elem.vxx_up = 0
                elem.vxy_up = 0
                elem.vx_down = 0
                elem.vy_down = 0
                elem.vxx_down = 0
                elem.vxy_down = 0
    lattice.update_transfer_maps()

    beam = Beam()
    beam.E = 2.4

    r_matrix = lattice_transfer_map(lattice, beam.E)


    if update_ref_values:
        return {'r': r_matrix.tolist(), 't': lattice.T.tolist()}

    matrices_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')

    result = check_matrix(r_matrix, matrices_ref["r"], TOL, assert_info=' r_matrix - ')
    result2 = check_matrix(lattice.T, matrices_ref["t"], TOL, assert_info=' t_matrix - ')
    assert check_result(result + result2)


@pytest.mark.parametrize('parameter', [0, 1])
def test_track_wo_coupler(lattice, p_array, parameter, update_ref_values=False):
    """R maxtrix test"""
    if parameter == 0:
        for elem in lattice.sequence:
            if elem.__class__ is Cavity:
                elem.vx_up = -5.6813e-5 + 1.0751e-5j
                elem.vy_up = -4.1091e-5 + 5.739e-7j
                elem.vxx_up = 0.00099943 - 0.00081401j
                elem.vxy_up = 0.0034065 - 0.0004146j
                elem.vx_down = -2.4014e-5 + 1.2492e-5j
                elem.vy_down = 3.6481e-5 + 7.9888e-6j
                elem.vxx_down = -0.004057 - 0.0001369j
                elem.vxy_down = 0.0029243 - 1.2891e-5j
    else:
        for elem in lattice.sequence:
            if elem.__class__ is Cavity:
                elem.vx_up = 0
                elem.vy_up = 0
                elem.vxx_up = 0
                elem.vxy_up = 0
                elem.vx_down = 0
                elem.vy_down = 0
                elem.vxx_down = 0
                elem.vxy_down = 0
    lattice.update_transfer_maps()

    navi = Navigator(lattice)
    navi.unit_step = 0.05
    p_array0 = copy.deepcopy(p_array)
    tws_track, p_array0 = track(lattice, p_array0, navi)

    tws = obj2dict(tws_track)
    p = obj2dict(p_array0)

    if update_ref_values:
        return {'tws_track': tws, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')
    result1 = check_dict(tws, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], TOL, assert_info=' p_array - ')

    assert check_result(result1 + result2)


def test_twiss(lattice, p_array, parameter=None, update_ref_values=False):
    """Twiss parameters calculation function test"""
    
    beam = Beam()
    beam.E = 2.4
    beam.beta_x = 41.1209
    beam.beta_y = 86.3314
    beam.alpha_x = 1.9630
    beam.alpha_y = 4.0972

    tw0 = Twiss(beam)
    tws = twiss(lattice, tw0, nPoints=None)
    
    tws = obj2dict(tws)
    
    if update_ref_values:
        return tws
    
    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws - ')
    assert check_result(result)


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### RF_TWISS START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### RF_TWISS END ###\n\n\n')
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
    update_functions.append('test_twiss')
    #update_functions.append('test_lattice_transfer_map_after_matching')
    #update_functions.append('test_twiss_after_matching')

    update_functions.append('test_lattice_r_t_maps_wo_coupler')
    update_functions.append('test_track_wo_coupler')
    update_function_parameters = {}
    update_function_parameters['test_lattice_r_t_maps_wo_coupler'] = [0, 1]
    update_function_parameters['test_track_wo_coupler'] = [0, 1]


    parameter = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']
    if cmdopt in update_functions:
        for p in parameter:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p_arr, p, True)

            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')

            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
