"""Test of the demo file demos/ebeam/csr_ex.py"""

import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from phys_proc_conf import *


def test_generate_parray(lattice, p_array, parameter=None, update_ref_values=False):
    """ func generate_parray testing """

    p = obj2dict(p_array)

    if update_ref_values:
        return {'p_array': p}

    p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    result = check_dict(p, p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result)


@pytest.mark.parametrize('parameter', [0, 1])
def test_track_smooth(lattice, p_array, parameter, update_ref_values=False):
    """
    test Runge_Kutta transfer map for undulator

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.05

    smb = SmoothBeam()
    if parameter == 1:
        navi.add_physics_proc(smb, lattice.sequence[0], lattice.sequence[0])

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi)

    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}


    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')


    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result1 + result2)


def test_track_beam_transform(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test Runge_Kutta transfer map for undulator

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.05

    tws = Twiss()
    tws.beta_x = 20.0
    tws.beta_y = 15.0

    bt = BeamTransform(tws=tws)
    bt.remove_offsets = True
    navi.add_physics_proc(bt, lattice.sequence[-1], lattice.sequence[-1])

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi)

    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}


    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')


    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], tolerance=1.0e-12, tolerance_type='absolute', assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result1 + result2)


def test_track_smooth_csr(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test Runge_Kutta transfer map for undulator

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.05

    smb = SmoothBeam()

    navi.add_physics_proc(smb, lattice.sequence[0], lattice.sequence[0])

    csr = CSR()
    csr.step = 1
    navi.add_physics_proc(csr, lattice.sequence[0], lattice.sequence[-1])

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi)

    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}


    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')


    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result1 + result2)

def test_track_smooth_csr_sc(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test Runge_Kutta transfer map for undulator

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.05

    smb = SmoothBeam()

    navi.add_physics_proc(smb, lattice.sequence[0], lattice.sequence[0])

    csr = CSR()
    csr.step = 1
    navi.add_physics_proc(csr, lattice.sequence[0], lattice.sequence[-1])

    sc = SpaceCharge()
    sc.step = 5
    navi.add_physics_proc(sc, lattice.sequence[0], lattice.sequence[-1])

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi)

    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}


    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')


    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=1.0e-14, tolerance_type='absolute', assert_info=' p - ')
    assert check_result(result1 + result2)


def test_track_laser_modulator(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test PhysicsProc LaserModulator

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """

    p_array_track = copy.deepcopy(p_array)

    navi = Navigator(lattice)
    navi.unit_step = 0.1

    lm = LaserModulator()
    lm.dE = 12500e-9  # GeV
    lm.Ku = 1.294  # undulator parameter
    lm.Lu = 0.74  # [m] - undulator length
    lm.lperiod = 0.074  # [m] - undulator period length
    lm.sigma_l = 300  # [m]
    lm.sigma_x = 300e-6  # [m]
    lm.sigma_y = 300e-6  # [m]
    lm.x_mean = 0
    lm.y_mean = 0

    navi.add_physics_proc(lm, m1, m2)

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi)

    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result1 + result2)

def test_track_lsc(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test PhysicsProc LaserModulator

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """

    p_array_track = copy.deepcopy(p_array)

    navi = Navigator(lattice)
    navi.unit_step = 0.1

    lsc = LSC()
    lsc.smooth_param = 0.1

    navi.add_physics_proc(lsc, lattice.sequence[0], lattice.sequence[-1])

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi)

    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result1 + result2)

def test_track_spontan_rad_effects(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test PhysicsProc LaserModulator

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """
    p_array_track = copy.deepcopy(p_array)

    navi = Navigator(lattice)
    navi.unit_step = 0.1

    rad = SpontanRadEffects(K=4, lperiod=0.05)
    rad.type = "planar"
    rad.energy_loss = True
    rad.quant_diff = True

    navi.add_physics_proc(rad, m1, m2)

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi)

    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
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
    update_functions.append('test_generate_parray')
    update_functions.append('test_track_smooth')
    update_functions.append('test_track_beam_transform')
    update_functions.append('test_track_smooth_csr')
    update_functions.append('test_track_smooth_csr_sc')
    update_functions.append('test_track_laser_modulator')
    update_functions.append('test_track_lsc')
    update_functions.append('test_track_spontan_rad_effects')

    update_function_parameters = {}
    update_function_parameters['test_track_smooth'] = [0, 1]

    parametr = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parametr:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p_arr, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
