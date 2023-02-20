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
    
    result = check_matrix(r_matrix, r_matrix_ref, tolerance=1.0e-10, tolerance_type='absolute', assert_info=' r_matrix - ')
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


def test_single_plane_dipole_wake(lattice, p_array, parameter=None, update_ref_values=False):

    wq = single_plane_dipole_wake(p=0.5e-3, t=0.25e-3, b=500e-6, l=5)

    s = np.linspace(0, 50e-6, num=20)
    W = np.array([[wq(si) for si in s]])
    if update_ref_values:
        return numpy2json(W)

    W_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))

    result = check_matrix(W, W_ref, tolerance=1.0e-10, tolerance_type='absolute',
                          assert_info=' dipole wake - ')
    assert check_result(result)

def test_single_plane_quad_wake(lattice, p_array, parameter=None, update_ref_values=False):

    wq = single_plate_quadrupole_wake(p=0.5e-3, t=0.25e-3, b=500e-6, l=5)

    s = np.linspace(0, 50e-6, num=20)
    W = np.array([[wq(si) for si in s]])
    if update_ref_values:
        return numpy2json(W)

    W_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))

    result = check_matrix(W, W_ref, tolerance=1.0e-10, tolerance_type='absolute',
                          assert_info=' quad wake - ')
    assert check_result(result)

def test_resolution(lattice, p_array, parameter=None, update_ref_values=False):
    R = np.array([[0.236135, -6.412734, 0.000000, 0.000000, 0.000000, 0.380000],
                  [0.133821, 0.600669, 0.000000, 0.000000, 0.000000, 0.080015],
                  [0.000000, 0.000000, 0.938983, -35.735219, 0.000000, 0.000000],
                  [0.000000, 0.000000, -0.031794, 2.274971, 0.000000, 0.000000],
                  [0.031958, 0.741370, 0.000000, 0.000000, 1.000000, -0.000199],
                  [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000]])

    distance = 500e-6

    tw1 = Twiss()
    tw1.beta_x = 26.384
    tw1.alpha_x = 0.495
    tw1.beta_y = 42.413
    tw1.alpha_y = -1.195
    parray = generate_parray(sigma_x=1e-4, sigma_px=2e-5, sigma_tau=1e-3 / 170, sigma_p=1e-4, chirp=0.00,
                             charge=250e-12,
                             nparticles=20000, energy=14, tws=tw1, shape="gauss")

    I = parray.I()
    print(np.shape(I))

    wyd = single_plane_dipole_wake(p=0.5e-3, t=0.25e-3, b=distance, l=5)
    wyq = single_plate_quadrupole_wake(p=0.5e-3, t=0.25e-3, b=distance, l=5)

    quad_kick = convolve_beam(I, wyq)
    dipole_kick = convolve_beam(I, wyd)

    r_temp, r_energy, sigma_x2, sigma_y2 = passive_streaker_resolutions(dipole_kick, quad_kick, R, tw1, kick="vert",
                                                                 emittn_x=0.6e-6,
                                                                 emittn_y=0.6e-6, energy=14, sigma_R=30e-6)
    R = r_temp[::10,:]
    if update_ref_values:
        return numpy2json(R)

    R_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))

    result = check_matrix(R, R_ref, tolerance=1.0e-10, tolerance_type='absolute',
                          assert_info=' temp res - ')
    assert check_result(result)


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
    update_functions.append('test_single_plane_dipole_wake')
    update_functions.append('test_single_plane_quad_wake')
    update_functions.append("test_resolution")

    update_function_parameters = {}

    parameter = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parameter:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p_arr, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
