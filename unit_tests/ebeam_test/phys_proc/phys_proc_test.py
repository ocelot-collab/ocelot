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

def test_s2current(lattice, p_array, parameter=None, update_ref_values=False):
    """ func generate_parray testing """

    I = s2current(s_array=p_array.tau(), q_array=p_array.q_array, n_points=300, filter_order=5, mean_vel=speed_of_light)

    if update_ref_values:
        return {'s': list(I[:, 0]), 'I': list(I[:, 1])}

    current_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    result1 = check_matrix(I[:, 0], current_ref['s'], TOL, assert_info=' s - ')
    result2 = check_matrix(I[:, 1], current_ref['I'], TOL, assert_info=' I - ')
    assert check_result(result1 + result2)


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

@pytest.mark.parametrize('parameter', [0, 1, 2])
def test_track_aperture(lattice, p_array, parameter, update_ref_values=False):
    """
    test RectAperture

    0 - all aperture activated
    1 - first aperture activated
    2 - second aperture activated
    """

    p_array_track = copy.deepcopy(p_array)

    navi = Navigator(lattice)

    if parameter == 0:
        navi.activate_apertures()
    elif parameter == 1:
        navi.activate_apertures(lattice.sequence[0], lattice.sequence[12])
    else:
        navi.activate_apertures(lattice.sequence[12])

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi)

    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result1 + result2)


@pytest.mark.parametrize('parameter', [0, 1, 2])
def test_track_ellipt_aperture(lattice, p_array, parameter, update_ref_values=False):
    """
    test EllipticalAperture

    0 - all aperture activated
    1 - first aperture activated
    2 - second aperture activated
    """

    p_array_track = copy.deepcopy(p_array)
    apx.y = 0.00005
    apx.type = "ellipt"
    apy.x = 0.00005
    apy.type = "ellipt"
    lattice.update_transfer_maps()
    navi = Navigator(lattice)

    if parameter == 0:
        navi.activate_apertures()
    elif parameter == 1:
        navi.activate_apertures(lattice.sequence[0], lattice.sequence[12])
    else:
        navi.activate_apertures(lattice.sequence[12])

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi)

    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result1 + result2)

def test_track_navi_reset_position(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test Runge_Kutta transfer map for undulator

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """

    p_array_track = copy.deepcopy(p_array)

    navi = Navigator(lattice)
    sc = SpaceCharge()
    sc.step = 5
    navi.add_physics_proc(sc, lattice.sequence[0], lattice.sequence[-1])

    navi.activate_apertures()


    tws_track_1, p_array_1 = track(lattice, p_array_track, navi)

    tws_track_obj_1 = obj2dict(tws_track_1)
    p_obj_1 = obj2dict(p_array_1)

    navi.reset_position()
    p_array_track_2 = copy.deepcopy(p_array)
    tws_track_2, p_array_2 = track(lattice, p_array_track_2, navi)

    tws_track_obj_2 = obj2dict(tws_track_2)
    p_obj_2 = obj2dict(p_array_2)

    result1 = check_dict(tws_track_obj_1, tws_track_obj_2, TOL, assert_info=' tws_track - ')
    result2 = check_dict(p_obj_2, p_obj_2, tolerance=TOL, assert_info=' p - ')
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
    # result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
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

def test_dechirper_offaxis(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test PhysicsProc WakeTable

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """
    p_array_track = copy.deepcopy(p_array)

    navi = Navigator(lattice)
    navi.unit_step = 0.1

    wake_table = WakeTableDechirperOffAxis(b=500 * 1e-6)
    ws = Wake()
    ws.wake_table = wake_table

    navi.add_physics_proc(ws, m1, m1)

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi)

    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result1 + result2)

@pytest.mark.parametrize('parameter', [0, 1, 2])
def test_phase_space_aperture(lattice, p_array, parameter, update_ref_values=False):
    """
    test PhysicsProc WakeTable

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """
    p_array_track = copy.deepcopy(p_array)

    navi = Navigator(lattice)

    ap1 = PhaseSpaceAperture()
    ap1.taumin = -3
    ap1.taumax = 3

    ap2 = PhaseSpaceAperture()
    ap2.horizontal = True
    ap2.xmin = -1
    ap2.xmax = 1
    if parameter == 0:
        navi.add_physics_proc(ap1, m1, m1)
    elif parameter == 1:
        navi.add_physics_proc(ap2, m1, m1)
    else:
        navi.add_physics_proc(ap1, m1, m1)
        navi.add_physics_proc(ap2, m1, m1)

    print(p_array_track)
    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi)
    print(p_array_wo)
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
    update_functions.append("test_s2current")
    update_functions.append('test_track_smooth')
    update_functions.append('test_track_aperture')
    update_functions.append("test_track_ellipt_aperture")
    update_functions.append('test_track_beam_transform')
    update_functions.append('test_track_smooth_csr')
    update_functions.append('test_track_smooth_csr_sc')
    update_functions.append('test_track_laser_modulator')
    update_functions.append('test_track_lsc')
    update_functions.append('test_track_spontan_rad_effects')
    update_functions.append("test_dechirper_offaxis")
    update_functions.append("test_phase_space_aperture")

    update_function_parameters = {}
    update_function_parameters['test_track_smooth'] = [0, 1]
    update_function_parameters['test_track_aperture'] = [0, 1, 2]
    update_function_parameters['test_track_ellipt_aperture'] = [0, 1, 2]
    update_function_parameters['test_phase_space_aperture'] = [0, 1, 2]

    parametr = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parametr:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p_arr, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
