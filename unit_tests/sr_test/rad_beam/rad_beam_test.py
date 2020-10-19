"""Test of the demo file demos/sr/spatial.py"""

import os
import sys
import time
import copy

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from rad_beam_conf import *


@pytest.mark.parametrize('parameter', [0, 1, 2])
def test_coherent_radiation(lattice, screen, beam, parameter, update_ref_values=False):
    """calculate_radiation function test"""
    screen.nullify()
    accuracy = 1
    beam_c = copy.deepcopy(beam)
    if parameter == 0:
        p_array = beam_c

    elif parameter == 1:
        p_array = beam_c
        p_array.tau()[:] *= -1

    elif parameter == 2:
        p_array = beam_c
        p_array.tau()[:] = 0

    screen = coherent_radiation(lattice, screen, p_array, accuracy=accuracy)

    if update_ref_values:
        return {'Eph': screen.Eph.tolist(), 'Yph': screen.Yph.tolist(), 'Xph': screen.Xph.tolist(),
                'Total': screen.Total.tolist(), 'Sigma': screen.Sigma.tolist(), 'Pi': screen.Pi.tolist()}

    screen_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')

    result1 = check_matrix(screen.Eph, screen_ref['Eph'], tolerance=1.0e-7, tolerance_type='relative', assert_info=' Eph - ')
    result2 = check_matrix(screen.Yph, screen_ref['Yph'], tolerance=1.0e-7, tolerance_type='relative', assert_info=' Yph - ')
    result3 = check_matrix(screen.Xph, screen_ref['Xph'], tolerance=1.0e-7, tolerance_type='relative', assert_info=' Xph - ')
    result4 = check_matrix(screen.Total, screen_ref['Total'], tolerance=1.0e-7, tolerance_type='relative', assert_info=' Total - ')
    result5 = check_matrix(screen.Sigma, screen_ref['Sigma'], tolerance=1.0e-7, tolerance_type='relative', assert_info=' Sigma - ')
    result6 = check_matrix(screen.Pi, screen_ref['Pi'], tolerance=1.0e-7, tolerance_type='relative', assert_info=' Pi - ')
    assert check_result(result1 + result2 + result3 + result4 + result5 + result6)

@pytest.mark.parametrize('parameter', [0, 1, 2])
def test_coherent_radiation_traj(lattice, screen, beam, parameter, update_ref_values=False):
    """calculate_radiation function test"""
    screen.nullify()
    accuracy = 1
    beam_c = copy.deepcopy(beam)
    if parameter == 0:
        p_array = beam_c

    elif parameter == 1:
        p_array = beam_c
        p_array.tau()[:] *= -1

    elif parameter == 2:
        p_array = beam_c
        p_array.tau()[:] = 0

    screen = coherent_radiation(lattice, screen, p_array, accuracy=accuracy)

    if update_ref_values:
        return {'x0': list(screen.beam_traj.x(0)[::10]), 'y0': list(screen.beam_traj.y(0)[::10]), 'xp0': list(screen.beam_traj.xp(0)[::10]),
                'yp0': list(screen.beam_traj.yp(0)[::10]), 'z0': list(screen.beam_traj.z(0)[::10]), 's0': list(screen.beam_traj.s(0)[::10]),
                'x1': list(screen.beam_traj.x(1)[::10]), 'x2': list(screen.beam_traj.x(2)[::10])}

    screen_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')

    result1 = check_matrix(screen.beam_traj.x(0)[::10], screen_ref['x0'], TOL, assert_info=' x - ')
    result2 = check_matrix(screen.beam_traj.y(0)[::10], screen_ref['y0'], TOL, assert_info=' y - ')
    result3 = check_matrix(screen.beam_traj.xp(0)[::10], screen_ref['xp0'], TOL, assert_info=' xp - ')
    result4 = check_matrix(screen.beam_traj.yp(0)[::10], screen_ref['yp0'], TOL, assert_info=' yp - ')
    result5 = check_matrix(screen.beam_traj.z(0)[::10], screen_ref['z0'], TOL, assert_info=' z - ')
    result6 = check_matrix(screen.beam_traj.s(0)[::10], screen_ref['s0'], TOL, assert_info=' s - ')
    result7 = check_matrix(screen.beam_traj.x(1)[::10], screen_ref['x1'], TOL, assert_info=' x1 - ')
    result8 = check_matrix(screen.beam_traj.x(2)[::10], screen_ref['x2'], TOL, assert_info=' x2 - ')
    assert check_result(result1 + result2 + result3 + result4 + result5 + result6 + result7 + result8)

@pytest.mark.parametrize('parameter', [0, 1, 2])
def test_coherent_radiation_fields(lattice, screen, beam, parameter, update_ref_values=False):
    """calculate_radiation fucntion test"""
    screen.nullify()
    accuracy = 1
    beam_c = copy.deepcopy(beam)
    if parameter == 0:
        s = Screen()
        s.z = 1000.0
        s.size_x = 15
        s.size_y = 15
        s.nx = 101
        s.ny = 1
        s.start_energy = 0.00850446  # eV
        s.end_energy = 15e-3  # eV
        s.num_energy = 1

    elif parameter == 1:
        s = Screen()
        s.z = 1000.0
        s.size_x = 15
        s.size_y = 15
        s.nx = 11
        s.ny = 11
        s.start_energy = 0.00850446  # eV
        s.end_energy = 15e-3  # eV
        s.num_energy = 1

    elif parameter == 2:
        s = Screen()
        s.z = 1000.0
        s.size_x = 15
        s.size_y = 15
        s.nx = 11
        s.ny = 11
        s.start_energy = 0.00850446  # eV
        s.end_energy = 15e-3  # eV
        s.num_energy = 5

    p_array = beam_c

    screen = coherent_radiation(lattice, screen, p_array, accuracy=accuracy)

    if update_ref_values:
        return {'arReEx': screen.arReEx.tolist(), 'arImEx': screen.arImEx.tolist(), 'arReEy': screen.arReEy.tolist(),
                'arImEy': screen.arImEy.tolist(), 'arPhase': screen.arPhase.tolist()}

    screen_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')

    result1 = check_matrix(screen.arReEx, screen_ref['arReEx'], TOL, assert_info=' arReEx - ')
    result2 = check_matrix(screen.arImEx, screen_ref['arImEx'], TOL, assert_info=' arImEx - ')
    result3 = check_matrix(screen.arReEy, screen_ref['arReEy'], TOL, assert_info=' arReEy - ')
    result4 = check_matrix(screen.arImEy, screen_ref['arImEy'], TOL, assert_info=' arImEy - ')
    result5 = check_matrix(screen.arPhase, screen_ref['arPhase'], TOL, assert_info=' arPhase - ')
    assert check_result(result1 + result2 + result3 + result4 + result5 )


@pytest.mark.parametrize('parameter', [0, 1, 2])
def test_coherent_radiation_parray(lattice, screen, beam, parameter, update_ref_values=False):
    """calculate_radiation function test"""
    screen.nullify()
    accuracy = 1
    beam_c = copy.deepcopy(beam)
    if parameter == 0:
        p_array = beam_c

    elif parameter == 1:
        p_array = beam_c
        p_array.tau()[:] *= -1

    elif parameter == 2:
        p_array = beam_c
        p_array.tau()[:] = 0

    screen = coherent_radiation(lattice, screen, p_array, accuracy=accuracy)

    p = obj2dict(p_array)
    if update_ref_values:
        return {'p_array': p}

    parray_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')
    result1 = check_dict(p, parray_ref['p_array'], tolerance=1.0e-12, tolerance_type='absolute', assert_info=' p - ')

    assert check_result(result1)


def test_energy_loss_spatial(lattice, screen, beam, parameter=None, update_ref_values=False):
    """calculate_radiation fucntion test"""
    beam = Beam()
    beam.E = 17.5
    beam.I = 0.1  # A

    screen = Screen()
    screen.z = 5000.0
    screen.size_x = 0.05
    screen.size_y = 0.05
    screen.nx = 10
    screen.ny = 10

    screen.start_energy = 8078  # eV
    screen.end_energy = 8090  # eV
    screen.num_energy = 1


    U40_short = Undulator(nperiods=5, lperiod=0.040, Kx=4, eid="und")

    seq = (U40_short,)*10
    lat = MagneticLattice(seq)
    screen = calculate_radiation(lat, copy.deepcopy(screen), beam, energy_loss=True, accuracy=1)



    if update_ref_values:
        return {'Eph': screen.Eph.tolist(), 'Yph': screen.Yph.tolist(), 'Xph': screen.Xph.tolist(),
                'Total': screen.Total.tolist(), 'Sigma': screen.Sigma.tolist(), 'Pi': screen.Pi.tolist()}

    screen_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_matrix(screen.Eph, screen_ref['Eph'], TOL, assert_info=' Eph - ')
    result2 = check_matrix(screen.Yph, screen_ref['Yph'], TOL, assert_info=' Yph - ')
    result3 = check_matrix(screen.Xph, screen_ref['Xph'], TOL, assert_info=' Xph - ')
    result4 = check_matrix(screen.Total, screen_ref['Total'], TOL, assert_info=' Total - ')
    result5 = check_matrix(screen.Sigma, screen_ref['Sigma'], TOL, assert_info=' Sigma - ')
    result6 = check_matrix(screen.Pi, screen_ref['Pi'], TOL, assert_info=' Pi - ')
    assert check_result(result1 + result2 + result3 + result4 + result5 + result6)


def test_energy_loss_spectrum(lattice, screen, beam, parameter=None, update_ref_values=False):
    """calculate_radiation fucntion test"""
    beam = Beam()
    beam.E = 17.5
    beam.I = 0.1  # A

    screen = Screen()
    screen.z = 5000.0
    screen.size_x = 0.0
    screen.size_y = 0.0
    screen.nx = 1
    screen.ny = 1

    screen.start_energy = 8030  # eV
    screen.end_energy = 8090  # eV
    screen.num_energy = 100


    U40_short = Undulator(nperiods=5, lperiod=0.040, Kx=4, eid="und")

    seq = (U40_short,)*10
    lat = MagneticLattice(seq)
    screen = calculate_radiation(lat, copy.deepcopy(screen), beam, energy_loss=True, accuracy=1)



    if update_ref_values:
        return {'Eph': screen.Eph.tolist(), 'Yph': screen.Yph.tolist(), 'Xph': screen.Xph.tolist(),
                'Total': screen.Total.tolist(), 'Sigma': screen.Sigma.tolist(), 'Pi': screen.Pi.tolist()}

    screen_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_matrix(screen.Eph, screen_ref['Eph'], TOL, assert_info=' Eph - ')
    result2 = check_matrix(screen.Yph, screen_ref['Yph'], TOL, assert_info=' Yph - ')
    result3 = check_matrix(screen.Xph, screen_ref['Xph'], TOL, assert_info=' Xph - ')
    result4 = check_matrix(screen.Total, screen_ref['Total'], TOL, assert_info=' Total - ')
    result5 = check_matrix(screen.Sigma, screen_ref['Sigma'], TOL, assert_info=' Sigma - ')
    result6 = check_matrix(screen.Pi, screen_ref['Pi'], TOL, assert_info=' Pi - ')
    assert check_result(result1 + result2 + result3 + result4 + result5 + result6)

def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### RAD BEAM START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### RAD BEAM END ###\n\n\n')
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
def test_update_ref_values(lattice, screen, beam, cmdopt):
    
    update_functions = []
    update_functions.append('test_coherent_radiation')
    update_functions.append('test_coherent_radiation_traj')
    update_functions.append('test_coherent_radiation_fields')
    update_functions.append("test_coherent_radiation_parray")
    update_functions.append('test_energy_loss_spatial')
    update_functions.append('test_energy_loss_spectrum')

    update_function_parameters = {}
    update_function_parameters['test_coherent_radiation'] = [0, 1, 2]
    update_function_parameters['test_coherent_radiation_traj'] = [0, 1, 2]
    update_function_parameters['test_coherent_radiation_fields'] = [0, 1, 2]
    update_function_parameters['test_coherent_radiation_parray'] = [0, 1, 2]
    parameter = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parameter:
            result = eval(cmdopt)(lattice, screen, beam, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
