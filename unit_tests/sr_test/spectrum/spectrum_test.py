"""Test of the demo file demos/sr/spectrum.py"""

import os
import sys
import time
import copy
FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from spectrum_conf import *


def test_calculate_radiation(lattice, screen, beam, update_ref_values=False):
    """calculate_radiation fucntion test"""

    screen = calculate_radiation(lattice, screen, beam)

    if update_ref_values:
        return {'Eph':screen.Eph.tolist(), 'Yph':screen.Yph.tolist(), 'Xph':screen.Xph.tolist(), 'Total':screen.Total.tolist(), 'Sigma':screen.Sigma.tolist(), 'Pi':screen.Pi.tolist()}
    
    screen_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result1 = check_matrix(screen.Eph, screen_ref['Eph'], TOL, assert_info=' Eph - ')
    result2 = check_matrix(screen.Yph, screen_ref['Yph'], TOL, assert_info=' Yph - ')
    result3 = check_matrix(screen.Xph, screen_ref['Xph'], TOL, assert_info=' Xph - ')
    result4 = check_matrix(screen.Total, screen_ref['Total'], TOL, assert_info=' Total - ')
    result5 = check_matrix(screen.Sigma, screen_ref['Sigma'], TOL, assert_info=' Sigma - ')
    result6 = check_matrix(screen.Pi, screen_ref['Pi'], TOL, assert_info=' Pi - ')
    assert check_result(result1+result2+result3+result4+result5+result6)


def test_calculate_radiation_endpoles(lattice, screen, beam, update_ref_values=False):
    """calculate_radiation fucntion test"""

    screen = calculate_radiation(lattice, screen, beam, end_poles=True)

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

def test_segments(lattice, screen, beam, update_ref_values=False):
    """calculate_radiation fucntion test"""
    beam = Beam()
    beam.E = 17.5
    beam.I = 0.1  # A

    screen = Screen()
    screen.z = 5000.0
    screen.size_x = 0.00
    screen.size_y = 0.00
    screen.nx = 1
    screen.ny = 1

    screen.start_energy = 8030  # eV
    screen.end_energy = 8090  # eV
    screen.num_energy = 100


    U40_short = Undulator(nperiods=5, lperiod=0.040, Kx=4, eid="und")

    seq = (U40_short,)*5
    lat = MagneticLattice(seq)
    screen_segm = calculate_radiation(lat, copy.deepcopy(screen), beam, accuracy=2)

    U40 = Undulator(nperiods=25, lperiod=0.040, Kx=4, eid="und")

    lat = MagneticLattice((U40,))
    screen_whole = calculate_radiation(lat, copy.deepcopy(screen), beam, accuracy=3)


    #if update_ref_values:
    #    return {'Eph': screen.Eph.tolist(), 'Yph': screen.Yph.tolist(), 'Xph': screen.Xph.tolist(),
    #            'Total': screen.Total.tolist(), 'Sigma': screen.Sigma.tolist(), 'Pi': screen.Pi.tolist()}

    #screen_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_matrix(screen_segm.Eph,   screen_whole.Eph   , TOL, assert_info=' Eph - ')
    result2 = check_matrix(screen_segm.Yph,   screen_whole.Yph   , TOL, assert_info=' Yph - ')
    result3 = check_matrix(screen_segm.Xph,   screen_whole.Xph   , TOL, assert_info=' Xph - ')
    result4 = check_matrix(screen_segm.Total, screen_whole.Total , TOL, assert_info=' Total - ')
    result5 = check_matrix(screen_segm.Sigma, screen_whole.Sigma , TOL, assert_info=' Sigma - ')
    result6 = check_matrix(screen_segm.Pi,    screen_whole.Pi    , TOL, assert_info=' Pi - ')
    assert check_result(result1 + result2 + result3 + result4 + result5 + result6)


def test_calculate_rad_from_lattice(lattice, screen, beam, update_ref_values=False):
    """calculate_radiation fucntion test"""
    und = Undulator(Kx=4., nperiods=125, lperiod=0.04, eid="und")
    D = Drift(l=0.5, eid="D")
    b1 = Hcor(l=0.1, angle=10 * -0.00001, eid="b1")
    b2 = Hcor(l=0.2, angle=10 * 0.00002, eid="b2")
    b3 = Hcor(l=0.1, angle=10 * -0.00001, eid="b3")
    phase_shift = (b1, b2, b3)
    cell = (und, D, phase_shift, D, und)
    lattice = MagneticLattice(cell)

    screen = calculate_radiation(lattice, screen, beam)

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


def test_spectrum_on_axis_3D(lattice, screen, beam, update_ref_values=False):
    """calculate_radiation fucntion test"""

    beam = Beam()
    beam.E = 2.5

    beam.I = 0.1

    beam.beta_x = 12.84
    beam.beta_y = 6.11
    beam.Dx = 0.526

    und = Undulator(Kx=0.43, nperiods=500, lperiod=0.007, eid="und")

    lat = MagneticLattice((und))

    ###
    screen = Screen()
    screen.z = 100.0
    screen.size_x = 0.002  # m
    screen.size_y = 0.002  # m
    screen.nx = 3
    screen.ny = 3

    screen.start_energy = 7700  # eV
    screen.end_energy = 7800  # eV
    screen.num_energy = 100

    screen = calculate_radiation(lat, screen, beam)

    screen2 = Screen()
    screen2.z = 100.0
    screen2.size_x = 0.002  # m
    screen2.size_y = 0.002  # m
    screen2.nx = 1
    screen2.ny = 1

    screen2.start_energy = 7700  # eV
    screen2.end_energy = 7800  # eV
    screen2.num_energy = 100

    screen2 = calculate_radiation(lat, screen2, beam)

    n1 = screen.nx * screen.ny

    total1 = screen.Total[np.floor(n1 / 2).astype(int):-np.floor(n1 / 2).astype(int):n1]
    sigma1 = screen.Sigma[np.floor(n1 / 2).astype(int):-np.floor(n1 / 2).astype(int):n1]
    pi1 = screen.Pi[np.floor(n1 / 2).astype(int):-np.floor(n1 / 2).astype(int):n1]

    result1 = check_matrix(screen2.Eph, screen.Eph, TOL, assert_info=' Eph - ')
    result4 = check_matrix(screen2.Total, total1, TOL, assert_info=' Total - ')
    result5 = check_matrix(screen2.Sigma, sigma1, TOL, assert_info=' Sigma - ')
    result6 = check_matrix(screen2.Pi, pi1, TOL, assert_info=' Pi - ')
    assert check_result(result1 + result4 + result5 + result6)


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### SPECTRUM START ###\n\n')
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
    update_functions.append("test_calculate_radiation_endpoles")
    update_functions.append('test_segments')
    update_functions.append("test_calculate_rad_from_lattice")
    
    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, screen, beam, True)
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
