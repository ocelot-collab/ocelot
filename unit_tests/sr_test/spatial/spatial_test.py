"""Test of the demo file demos/sr/spatial.py"""

import os
import sys
import time
import copy

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from spatial_conf import *


@pytest.mark.parametrize('parameter', [0, 1, 2, 3])
def test_calculate_radiation(lattice, screen, beam, parameter, update_ref_values=False):
    """calculate_radiation fucntion test"""

    accuracy = 1

    if parameter == 0:
        accuracy = 2
        
        screen.size_x = 0.002
        screen.size_y = 0.0
        screen.nx = 100
        screen.ny = 1
        screen.start_energy = 7761.2
        screen.end_energy = 7900
        screen.num_energy = 1

    elif parameter == 1:
        screen.size_x = 0.002
        screen.size_y = 0.002
        screen.nx = 51
        screen.ny = 51
        screen.start_energy = 7761.2
        screen.end_energy = 7900
        screen.num_energy = 1

    elif parameter == 2:
        screen.size_x = 0.002
        screen.size_y = 0.002
        screen.nx = 1
        screen.ny = 1
        screen.start_energy = 7700
        screen.end_energy = 7800
        screen.num_energy = 100

    elif parameter == 3:
        screen.size_x = 0.002
        screen.size_y = 0.002
        screen.nx = 5
        screen.ny = 5
        screen.start_energy = 7700
        screen.end_energy = 7800
        screen.num_energy = 5

    screen = calculate_radiation(lattice, screen, beam, accuracy=accuracy)

    if update_ref_values:
        return {'Eph':screen.Eph.tolist(), 'Yph':screen.Yph.tolist(), 'Xph':screen.Xph.tolist(), 'Total':screen.Total.tolist(), 'Sigma':screen.Sigma.tolist(), 'Pi':screen.Pi.tolist()}
    
    screen_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')

    result1 = check_matrix(screen.Eph, screen_ref['Eph'], TOL, assert_info=' Eph - ')
    result2 = check_matrix(screen.Yph, screen_ref['Yph'], TOL, assert_info=' Yph - ')
    result3 = check_matrix(screen.Xph, screen_ref['Xph'], TOL, assert_info=' Xph - ')
    result4 = check_matrix(screen.Total, screen_ref['Total'], TOL, assert_info=' Total - ')
    result5 = check_matrix(screen.Sigma, screen_ref['Sigma'], TOL, assert_info=' Sigma - ')
    result6 = check_matrix(screen.Pi, screen_ref['Pi'], TOL, assert_info=' Pi - ')
    assert check_result(result1+result2+result3+result4+result5+result6)


def test_segments(lattice, screen, beam, update_ref_values=False):
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

    seq = (U40_short,)*5
    lat = MagneticLattice(seq)
    screen_segm = calculate_radiation(lat, copy.deepcopy(screen), beam, accuracy=4)

    U40 = Undulator(nperiods=25, lperiod=0.040, Kx=4, eid="und")

    lat = MagneticLattice((U40,))
    screen_whole = calculate_radiation(lat, copy.deepcopy(screen), beam, accuracy=4)


    #if update_ref_values:
    #    return {'Eph': screen.Eph.tolist(), 'Yph': screen.Yph.tolist(), 'Xph': screen.Xph.tolist(),
    #            'Total': screen.Total.tolist(), 'Sigma': screen.Sigma.tolist(), 'Pi': screen.Pi.tolist()}

    #screen_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_matrix(screen_segm.Eph,   screen_whole.Eph   , TOL, assert_info=' Eph - ')
    result2 = check_matrix(screen_segm.Yph,   screen_whole.Yph   , TOL, assert_info=' Yph - ')
    result3 = check_matrix(screen_segm.Xph,   screen_whole.Xph   , TOL, assert_info=' Xph - ')
    result4 = check_matrix(screen_segm.Total, screen_whole.Total , TOL, assert_info=' Total - ')
    result5 = check_matrix(screen_segm.Sigma, screen_whole.Sigma , TOL, assert_info=' Sigma - ')
    result6 = check_matrix(screen_segm.Pi,    screen_whole.Pi    , tolerance=1.0e-6, assert_info=' Pi - ')
    assert check_result(result1 + result2 + result3 + result4 + result5 + result6)

def test_segments_spectrum(lattice, screen, beam, update_ref_values=False):
    """calculate_radiation fucntion test"""
    beam = Beam()
    beam.E = 17.5
    beam.I = 0.1  # A

    screen = Screen()
    screen.z = 5000.0
    screen.size_x = 0.05
    screen.size_y = 0.05
    screen.nx = 1
    screen.ny = 1

    screen.start_energy = 7900  # eV
    screen.end_energy = 8200  # eV
    screen.num_energy = 100


    U40_short = Undulator(nperiods=5, lperiod=0.040, Kx=4, eid="und")

    seq = (U40_short,)*5
    lat = MagneticLattice(seq)
    screen_segm = calculate_radiation(lat, copy.deepcopy(screen), beam, accuracy=4)

    U40 = Undulator(nperiods=25, lperiod=0.040, Kx=4, eid="und")

    lat = MagneticLattice((U40,))
    screen_whole = calculate_radiation(lat, copy.deepcopy(screen), beam, accuracy=4)


    #if update_ref_values:
    #    return {'Eph': screen.Eph.tolist(), 'Yph': screen.Yph.tolist(), 'Xph': screen.Xph.tolist(),
    #            'Total': screen.Total.tolist(), 'Sigma': screen.Sigma.tolist(), 'Pi': screen.Pi.tolist()}

    #screen_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_matrix(screen_segm.Eph,   screen_whole.Eph   , TOL, assert_info=' Eph - ')
    result2 = check_matrix(screen_segm.Yph,   screen_whole.Yph   , TOL, assert_info=' Yph - ')
    result3 = check_matrix(screen_segm.Xph,   screen_whole.Xph   , TOL, assert_info=' Xph - ')
    result4 = check_matrix(screen_segm.Total, screen_whole.Total , TOL, assert_info=' Total - ')
    result5 = check_matrix(screen_segm.Sigma, screen_whole.Sigma , TOL, assert_info=' Sigma - ')
    result6 = check_matrix(screen_segm.Pi,    screen_whole.Pi    , tolerance=1.0e-6, assert_info=' Pi - ')
    assert check_result(result1 + result2 + result3 + result4 + result5 + result6)


def test_segments_with_3d_screen(lattice, screen, beam, update_ref_values=False):
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

    screen.start_energy = 8000  # eV
    screen.end_energy = 8100  # eV
    screen.num_energy = 10


    U40_short = Undulator(nperiods=5, lperiod=0.040, Kx=4, eid="und")

    seq = (U40_short,)*5
    lat = MagneticLattice(seq)
    screen_segm = calculate_radiation(lat, copy.deepcopy(screen), beam, accuracy=4)

    U40 = Undulator(nperiods=25, lperiod=0.040, Kx=4, eid="und")

    lat = MagneticLattice((U40,))
    screen_whole = calculate_radiation(lat, copy.deepcopy(screen), beam, accuracy=4)


    #if update_ref_values:
    #    return {'Eph': screen.Eph.tolist(), 'Yph': screen.Yph.tolist(), 'Xph': screen.Xph.tolist(),
    #            'Total': screen.Total.tolist(), 'Sigma': screen.Sigma.tolist(), 'Pi': screen.Pi.tolist()}

    #screen_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_matrix(screen_segm.Eph,   screen_whole.Eph   , TOL, assert_info=' Eph - ')
    result2 = check_matrix(screen_segm.Yph,   screen_whole.Yph   , TOL, assert_info=' Yph - ')
    result3 = check_matrix(screen_segm.Xph,   screen_whole.Xph   , TOL, assert_info=' Xph - ')
    result4 = check_matrix(screen_segm.Total, screen_whole.Total , TOL, assert_info=' Total - ')
    result5 = check_matrix(screen_segm.Sigma, screen_whole.Sigma , TOL, assert_info=' Sigma - ')
    result6 = check_matrix(screen_segm.Pi,    screen_whole.Pi    , tolerance=1.0e-6, assert_info=' Pi - ')
    assert check_result(result1 + result2 + result3 + result4 + result5 + result6)


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### SPATIAL START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### SPATIAL END ###\n\n\n')
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
    update_functions.append('test_calculate_radiation')
    update_functions.append('test_segments')
    update_functions.append("test_segments_spectrum")
    update_functions.append('test_segments_with_3d_screen')
    update_function_parameters = {}
    update_function_parameters['test_calculate_radiation'] = [0, 1, 2, 3]
    
    parametr = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parametr:
            result = eval(cmdopt)(lattice, screen, beam, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
