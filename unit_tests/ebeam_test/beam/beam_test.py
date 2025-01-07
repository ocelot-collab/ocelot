"""Test of the demo file demos/ebeam/csr_ex.py"""

import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from beam_conf import *


@pytest.mark.parametrize('parameter', [0, 1, 2])
def test_generate_parray_shape(parameter, update_ref_values=False):
    """R matrix calculation test"""
    if parameter == 0:
        shape = "gauss"
    elif parameter == 1:
        shape = "tri"
    else:
        shape = "rect"

    parray = generate_parray(sigma_x=1e-4, sigma_px=2e-5, sigma_y=None, sigma_py=None,
                    sigma_tau=1e-3, sigma_p=1e-4, chirp=0.0, charge=0.25e-9, nparticles=5000000, energy=0.13,
                    tau_trunc=2.5, tws=None, shape=shape)

    s, I = get_current(parray, num_bins=100)
    B = np.column_stack((s, I))
    if update_ref_values:
        return numpy2json(B)
    B_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json'))
    result = check_matrix(B, B_ref, tolerance=1.0e-1, tolerance_type='relative', assert_info=' current - ')

    assert check_result(result)


def test_generate_parray_arb_shape(parameter=None, update_ref_values=False):
    sigma_tau = 30e-6
    A = np.array([0.2, 0.9, 0.3])

    mu = (np.array([0.1, 0.3, 0.2]) - 0.5) * 2 * sigma_tau

    sigma = np.array([0.9, 0.8, 0.7]) * sigma_tau

    s = np.linspace(-2 * sigma_tau, 2 * sigma_tau, num=500)

    s_reshaped = s[:, np.newaxis]

    f = np.sum(A * np.exp(-(s_reshaped - mu) ** 2 / (2.0 * sigma ** 2)), axis=1)

    shape = [s, f]

    parray = generate_parray(sigma_tau=sigma_tau, tau_trunc=2.5, sigma_p=2.5/14000, chirp=-0.00, charge=250e-12, nparticles=5000000, energy=14, tws=None, shape=shape)

    s, I = get_current(parray, num_bins=100)
    B = np.column_stack((s, I))
    if update_ref_values:
        return numpy2json(B)
    B_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    result = check_matrix(B, B_ref, tolerance=1.0e-1, tolerance_type='relative', assert_info=' current - ')

    assert check_result(result)

def test_get_envelope_slice(parameter=None, update_ref_values=False):
    sigma_x = 1e-4
    sigma_px = 2e-5
    sigma_y = 5e-5
    sigma_py = 4e-5
    sigma_p = 1e-4
    parray = generate_parray(sigma_x=sigma_x, sigma_px=sigma_px, sigma_y=sigma_y, sigma_py=sigma_py,
                    sigma_tau=1e-3, sigma_p=sigma_p, chirp=0.0, charge=0.25e-9, nparticles=5000000, energy=0.13,
                    tau_trunc=2.5, tws=None, shape="gauss")

    tws = get_envelope(parray, bounds=[-0.5, 0.5], slice="Imax")

    assert np.isclose(sigma_x, np.sqrt(tws.xx), rtol=1e-05, atol=1e-07, equal_nan=False)
    assert np.isclose(sigma_px, np.sqrt(tws.pxpx), rtol=1e-05, atol=1e-07, equal_nan=False)
    assert np.isclose(sigma_y, np.sqrt(tws.yy), rtol=1e-05, atol=1e-07, equal_nan=False)
    assert np.isclose(sigma_py, np.sqrt(tws.pypy), rtol=1e-05, atol=1e-07, equal_nan=False)
    assert np.isclose(sigma_p, np.sqrt(tws.pp), rtol=1e-05, atol=1e-07, equal_nan=False)


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
def test_update_ref_values(cmdopt):

    update_functions = []
    update_functions.append('test_generate_parray_shape')
    update_functions.append('test_generate_parray_arb_shape')

    update_function_parameters = {}
    update_function_parameters['test_generate_parray_shape'] = [0, 1, 2]

    parameter = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']
    print(parameter)
    if cmdopt in update_functions:
        for p in parameter:
            print(p)
            result = eval(cmdopt)( p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')

