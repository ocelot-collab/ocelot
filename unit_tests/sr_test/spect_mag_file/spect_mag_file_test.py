"""Test of the demo file demos/sr/spect_mag_file.py"""

import os
import sys
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from spect_mag_file_conf import *


def test_calculate_radiation(lattice, screen, beam, update_ref_values=False):
    """calculate_radiation fucntion test"""

    screen = calculate_radiation(lattice, screen, beam, accuracy=2)

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


def test_radiation_file_and_function(lattice, screen, beam, update_ref_values=False):
    """calculate_radiation fucntion test"""
    lperiod = 0.04  # [m] undulator period
    nperiods = 30  # number of periods
    B0 = 1  # [T] amplitude of the magnetic field

    # longitudinal coordinates from 0 to lperiod*nperiods in [mm]
    z = np.linspace(0, lperiod * nperiods, num=8000) * 1000  # [mm]
    lperiod_mm = lperiod * 1000  # in [mm]
    By = B0 * np.cos(2 * np.pi / lperiod_mm * z)

    filed_map = np.vstack((z, By)).T

    np.savetxt("filed_map.txt", filed_map)
    und_m = Undulator(field_file="filed_map.txt", eid="und")

    lat_m = MagneticLattice((und_m))

    beam = Beam()
    beam.E = 17.5  # beam energy in [GeV]
    beam.I = 0.1  # beam current in [A]

    screen_ref = Screen()
    screen_ref.z = 1000.0  # distance from the begining of lattice to the screen

    screen_ref.start_energy = 7000  # [eV], starting photon energy
    screen_ref.end_energy = 12000  # [eV], ending photon energy
    screen_ref.num_energy = 100  # number of energy points[eV]

    # Calculate radiation
    screen_ref = calculate_radiation(lat_m, screen_ref, beam)


    # function
    und = Undulator(lperiod=lperiod, nperiods=nperiods, Kx=0.0, eid="und")
    und.mag_field = lambda x, y, z: (0., B0*np.cos(2*np.pi/lperiod*z), 0.)

    # next, all the same.

    lat = MagneticLattice((und))

    #beam = Beam()
    #beam.E = 17.5  # beam energy in [GeV]
    #beam.I = 0.1  # beam current in [A]

    screen = Screen()
    screen.z = 1000.0  # distance from the begining of lattice to the screen

    screen.start_energy = 7000  # [eV], starting photon energy
    screen.end_energy = 12000  # [eV], ending photon energy
    screen.num_energy = 100  # number of energy points[eV]

    # Calculate radiation
    screen = calculate_radiation(lat, screen, beam)

    result1 = check_matrix(screen.Eph, screen_ref.Eph, TOL, assert_info=' Eph - ')
    result2 = check_matrix(screen.Yph, screen_ref.Yph, TOL, assert_info=' Yph - ')
    result3 = check_matrix(screen.Xph, screen_ref.Xph, TOL, assert_info=' Xph - ')
    result4 = check_matrix(screen.Total, screen_ref.Total, tolerance=1.0e-4, assert_info=' Total - ')
    result5 = check_matrix(screen.Sigma, screen_ref.Sigma, tolerance=1.0e-4, assert_info=' Sigma - ')
    result6 = check_matrix(screen.Pi, screen_ref.Pi, tolerance=1.0e-4, assert_info=' Pi - ')
    assert check_result(result1 + result2 + result3 + result4 + result5 + result6)


def test_radiation_from_bm(lattice, screen, beam, update_ref_values=False):
    """calculate_radiation function test"""
    phi = 0.15046332005984778
    x_off = 0.018870166315482287

    beam.E = 2  # beam energy in [GeV]
    beam.I = 0.1  # beam current in [A]

    # set new initial coordinates for the beam
    beam.xp = phi / 2  # initial angle x'
    beam.x = -x_off  # initial offset

    By = 1.  # T - amplitude of vertical magnetic field.

    b = Undulator(lperiod=0.10, nperiods=10, eid="und")
    b.mag_field = lambda x, y, z: (0, By, 0)

    # in the Undulator element parameters lperiod and nperiods are needed
    # just for definition of the length of the element
    lattice = MagneticLattice((b,))

    screen = Screen()
    screen.z = 1000.0  # distance from the begining of lattice to the screen

    screen.start_energy = 100  # [eV], starting photon energy
    screen.end_energy = 10000  # [eV], ending photon energy
    screen.num_energy = 100  # number of energy points[eV]


    screen = calculate_radiation(lattice, screen, beam, accuracy=6)

    x = screen.beam_traj.x(0)[::50]
    z = screen.beam_traj.z(0)[::50]

    if update_ref_values:
        return {'Eph': screen.Eph.tolist(), 'Yph': screen.Yph.tolist(), 'Xph': screen.Xph.tolist(),
                'Total': screen.Total.tolist(), 'Sigma': screen.Sigma.tolist(), 'Pi': screen.Pi.tolist(),
                "x": x.tolist(), "z": z.tolist()}

    screen_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_matrix(screen.Eph, screen_ref['Eph'], TOL, assert_info=' Eph - ')
    result2 = check_matrix(screen.Yph, screen_ref['Yph'], TOL, assert_info=' Yph - ')
    result3 = check_matrix(screen.Xph, screen_ref['Xph'], TOL, assert_info=' Xph - ')
    result4 = check_matrix(screen.Total, screen_ref['Total'], TOL, assert_info=' Total - ')
    result5 = check_matrix(screen.Sigma, screen_ref['Sigma'], TOL, assert_info=' Sigma - ')
    result6 = check_matrix(screen.Pi, screen_ref['Pi'], TOL, assert_info=' Pi - ')

    result7 = check_matrix(x, screen_ref['x'], TOL, assert_info=' traj X - ')
    result8 = check_matrix(z, screen_ref['z'], TOL, assert_info=' traj Z - ')


    assert check_result(result1 + result2 + result3 + result4 + result5 + result6 + result7 + result8)


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### SPECT_MAG_FILE START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### SPECT_MAG_FILE END ###\n\n\n')
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
    update_functions.append('test_radiation_file_and_function')
    update_functions.append("test_radiation_from_bm")

    
    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, screen, beam, True)
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
