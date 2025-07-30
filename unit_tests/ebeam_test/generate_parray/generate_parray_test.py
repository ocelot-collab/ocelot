"""Test of the demo file demos/ebeam/dba_tracking.py"""

import os
import sys
from copy import copy
import time

import numpy as np
import pytest

from ocelot import *
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *


def test_default_beam(parameter=None, update_ref_values=False):
    """R maxtrix calculation test"""

    sigma_tau = 50e-6
    parray = generate_parray(sigma_x=1e-4, sigma_px=2e-5, sigma_y=None, sigma_py=None,
                    sigma_tau=sigma_tau, sigma_p=1e-4, chirp=0.01, charge=250e-12, nparticles=500000, energy=0.13,
                    tau_trunc=None, tws=None, shape="gauss")

    cur = parray.I()

    f = lambda x: max(cur[:, 1])*np.exp(-x ** 2 / (2. * sigma_tau ** 2))
    s = np.linspace(-5 * sigma_tau, 5 * sigma_tau, num=500)

    f_g = interp1d(cur[:, 0], cur[:, 1], bounds_error=False, fill_value=(0, 0))
    df = trapezoid(np.abs(f_g(s) - f(s)), s)

    output = f"value is {df} \n shape difference should be less than 0.001"

    result = [check_value(np.max(cur[:, 1]), 594, tolerance=1.0e-2, tolerance_type='relative', assert_info='', numerical_zero=1e-15),
              check_value(np.std(parray.tau()), 50e-6, tolerance=1.0e-2, tolerance_type='relative', assert_info='', numerical_zero=1e-15),
              None if df < 0.001 else output]

    assert check_result(result)

def test_tri_beam(parameter=None, update_ref_values=False):
    """R maxtrix calculation test"""

    sigma_tau = 50e-6
    parray = generate_parray(sigma_x=1e-4, sigma_px=2e-5, sigma_y=None, sigma_py=None,
                             sigma_tau=sigma_tau, sigma_p=1e-4, chirp=0.01, charge=250e-12, nparticles=500000,
                             energy=0.13,
                             tau_trunc=None, tws=None, shape="tri")

    cur = parray.I()

    f = lambda s: max(cur[:, 1]) * np.maximum(sigma_tau - np.abs(s), 0) / sigma_tau
    s = np.linspace(-2 * sigma_tau, 2 * sigma_tau, num=500)

    f_g = interp1d(cur[:, 0], cur[:, 1], bounds_error=False, fill_value=(0, 0))
    df = trapezoid(np.abs(f_g(s) - f(s)), s)

    output = f"value is {df} \n shape difference should be less than 0.005"

    result = [None if df < 0.005 else output]

    assert check_result(result)


def test_arbitrary_shape_beam(parameter=None, update_ref_values=False):
    """R maxtrix calculation test"""
    sigma_tau = 10e-6
    A1, A2, A3 = 0.5, 0.2, 1
    mu1, mu2, mu3 = -7.7e-08, 3.7e-06, -6.3e-06
    sigma1, sigma2, sigma3 = 7.2e-06, 9.2e-06, 3e-06

    f = lambda x: A1 * np.exp(-(x - mu1) ** 2 / (2. * sigma1 ** 2)) + A2 * np.exp(
        -(x - mu2) ** 2 / (2. * sigma2 ** 2)) + A3 * np.exp(-(x - mu3) ** 2 / (2. * sigma3 ** 2))

    s = np.linspace(-5 * sigma_tau, 5 * sigma_tau, num=200)

    shape = [s, f(s)]

    parray = generate_parray(sigma_x=1e-4, sigma_px=2e-5, sigma_y=None, sigma_py=None,
                             sigma_tau=sigma_tau, sigma_p=1e-4, chirp=0.01, charge=250e-12, nparticles=500000,
                             energy=0.13,
                             tau_trunc=None, tws=None, shape=shape)
    cur = parray.I()

    f_g = interp1d(cur[:, 0], cur[:, 1], bounds_error=False, fill_value=(0, 0))
    df = trapezoid(np.abs(f_g(s) - f(s) / max(f(s)) * max(cur[:, 1])), s)

    output = f"value is {df} \n shape difference should be less than 0.005"

    result = [None if df < 0.005 else output]

    assert check_result(result)
    

def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### DBA_TRACKING START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### DBA_TRACKING END ###\n\n\n')
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
def test_update_ref_values( cmdopt):
    
    update_functions = []


    
    update_function_parameters = {}

    parametr = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parametr:
            result = eval(cmdopt)( p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.npz'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.npz', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            np.savez( REF_RES_DIR + cmdopt + str(p) + '.npz', result)