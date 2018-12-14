"""Test of the demo file demos/ebeam/csr_ex.py"""

import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from laser_heater_conf import *


def test_generate_parray(lattice, p_array, parameter=None, update_ref_values=False):
    """ func generate_parray testing """

    p = obj2dict(p_array)

    if update_ref_values:
        return {'p_array': p}

    p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    result = check_dict(p, p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result)

@pytest.mark.parametrize('parameter', [0, 1])
def test_chicane_trajectory(lattice, p_array, parameter, update_ref_values=False):
    """
    test Runge_Kutta transfer map for undulator

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """

    # p_array_track = copy.deepcopy(p_array)

    for elem in lattice.sequence:
        if elem.__class__ == Undulator:
            if parameter == 1:
                elem.Kx = 0
            else:
                elem.Kx = 1.36 * 1.414213

    navi = Navigator(lattice)
    navi.unit_step = 0.05

    csr = CSR()
    csr.energy = p_array.E

    navi.add_physics_proc(csr, lattice.sequence[0], lattice.sequence[-1])

    if update_ref_values:
        return numpy2json(csr.csr_traj)

    traj_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json'))

    result = check_matrix(csr.csr_traj, traj_ref, tolerance=1.0e-14, tolerance_type='absolute', assert_info=' trajectory - ')

    assert check_result(result)



@pytest.mark.parametrize('parameter', [0, 1])
def test_track_chicane_und_wo_csr(lattice, p_array, parameter, update_ref_values=False):
    """
    test Runge_Kutta transfer map for undulator

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """

    p_array_track = copy.deepcopy(p_array)


    for elem in lattice.sequence:
        if elem.__class__ == Undulator:
            if parameter == 1:
                elem.Kx = 0
            else:
                elem.Kx = 1.36 * 1.414213

    lattice.update_transfer_maps()


    navi = Navigator(lattice)
    navi.unit_step = 0.05

    csr = CSR()
    csr.energy = p_array_track.E

    #navi.add_physics_proc(csr, lattice.sequence[0], lattice.sequence[-1])

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi)

    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}


    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')


    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result1 + result2)


@pytest.mark.parametrize('parameter', [0, 1])
def test_track_chicane_und_csr(lattice, p_array, parameter, update_ref_values=False):
    """
    test Runge_Kutta transfer map for undulator

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """

    # p_array_track = copy.deepcopy(p_array)

    for elem in lattice.sequence:
        if elem.__class__ == Undulator:
            if parameter == 1:
                elem.Kx = 0
            else:
                elem.Kx = 1.36 * 1.414213

    lattice.update_transfer_maps()

    navi = Navigator(lattice)
    navi.unit_step = 0.05

    csr = CSR()
    csr.energy = p_array.E

    navi.add_physics_proc(csr, lattice.sequence[0], lattice.sequence[-1])

    tws_track_wo, p_array_wo = track(lattice, p_array, navi)

    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result1 + result2)


@pytest.mark.parametrize('parameter', [0, 1])
def test_track_with_laser_heater(lattice, p_array, parameter, update_ref_values=False):
    """
    test Runge_Kutta transfer map for undulator

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """

    # p_array_track = copy.deepcopy(p_array)

    for elem in lattice.sequence:
        if elem.__class__ == Undulator:
            if parameter == 1:
                elem.Kx = 0
            else:
                elem.Kx = 1.36 * 1.414213

    lattice.update_transfer_maps()

    navi = Navigator(lattice)
    navi.unit_step = 0.05

    lh = LaserModulator()
    lh.step = 2
    lh.dE = 12500e-9  # GeV
    lh.Ku = 1.294  # undulator parameter
    lh.Lu = 0.8  # [m] - undulator length
    lh.lperiod = 0.074  # [m] - undulator period length
    lh.sigma_l = 30000e-2  # [m]
    lh.sigma_x = 300e-6
    lh.sigma_y = 300e-6
    navi.add_physics_proc(lh, und_start, und_stop)

    tws_track_wo, p_array_wo = track(lattice, p_array, navi)

    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result1 + result2)


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### LASER HEATER START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### LASER HEATER END ###\n\n\n')
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
    update_functions.append('test_generate_parray')
    update_functions.append('test_chicane_trajectory')
    update_functions.append('test_track_chicane_und_wo_csr')
    update_functions.append('test_track_chicane_und_csr')
    update_functions.append('test_track_with_laser_heater')

    update_function_parameters = {}
    update_function_parameters['test_chicane_trajectory'] = [0, 1]
    update_function_parameters['test_track_chicane_und_wo_csr'] = [0, 1]
    update_function_parameters['test_track_chicane_und_csr'] = [0, 1]
    update_function_parameters['test_track_with_laser_heater'] = [0, 1]

    parametr = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parametr:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p_arr, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
