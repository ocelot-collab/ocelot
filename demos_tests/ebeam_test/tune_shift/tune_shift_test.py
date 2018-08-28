"""Test of the demo file demos/ebeam/tune_shift.py"""

import os
import sys
import numpy as np

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from demos_tests.params import *
from tune_shift_conf import *


def test_lattice_transfer_map(lattice, update_ref_values=False):
    """R maxtrix calculation test"""

    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_twiss(lattice, update_ref_values=False):
    """Twiss parameters calculation function test"""
    
    mu_y_ref = 0.304171994243
    
    beam = Beam()
    beam.E = 0.0 #GeV
    beam.I = 0.1 #A

    tw0 = Twiss(beam)
    tws = twiss(lattice, tw0, nPoints=1000)
    
    mu_y_no_u = 1.0 - tws[-1].muy % (2.0*np.pi) / (2.0*np.pi)
    
    tws = obj2dict(tws)
    
    if update_ref_values:
        return tws

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_value(mu_y_no_u, mu_y_ref, TOL, assert_info=' mu_y_no_u - \n')
    result2 = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws - ')
    assert check_result([result1]+result2)


#@pytest.mark.skip(reason='TOO LONG')
def test_track_nturns(lattice, update_ref_values=False):
    """Track N turns function test"""

    pxy_list, nturns = track_nturns_wrapper(lattice)

    pxy_list = obj2dict(pxy_list, unpack=['particle'])
    
    if update_ref_values:
        return pxy_list
    
    pxy_list_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result = check_dict(pxy_list, pxy_list_ref, TOL, 'absolute', assert_info=' pxy_list track_nturns - ')
    assert check_result(result)


#@pytest.mark.skip(reason='TOO LOGN')
def test_freq_analysis(lattice, update_ref_values=False):
    """Frequency analysis function test"""

    mu_y_no_u = 0.304171994243
    mu_y_ref = 0.30322265625
    mu_y_h_ref = 0.3034375
    mu_y_sim_ref = 0.000734494242981
    
    pxy_list, nturns = track_nturns_wrapper(lattice)
    pxy_list = freq_analysis(pxy_list, lattice, nturns, harm=True)

    pxy_list = obj2dict(pxy_list, unpack=['particle'])
    
    if update_ref_values:
        return pxy_list
    
    pxy_list_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(pxy_list, pxy_list_ref, TOL, 'absolute', assert_info=' pxy_list freq_analysis- ')
    result2 = check_value(pxy_list[0]['muy'], mu_y_ref, TOL, assert_info=' muy - \n')

    x = np.linspace(0.3, 0.31, nturns+2)
    y = [p[2] for p in pxy_list[0]['p_list']]
    f = np.abs(dft(y, x))

    mu_y_h = x[np.argmax(f)]
    result3 = check_value(mu_y_h, mu_y_h_ref, TOL, assert_info=' mu_y_h - \n')

    mu_y_sim = mu_y_no_u - mu_y_h
    result4 = check_value(mu_y_sim, mu_y_sim_ref, TOL, assert_info=' mu_y_sim - \n')
    
    assert check_result(result1+[result2, result3, result4])


def dft(sample, freqs):

    n = len(sample)
    x_freq = freqs * n
    transf = np.zeros(len(freqs), dtype=np.complex)

    for i, ai in enumerate(sample):
        transf += ai * np.exp(-2.0 * np.pi * i / n * 1.0j * x_freq)
    return transf


def track_nturns_wrapper(lattice):

    nturns = 2048 - 1

    if not hasattr(pytest, 'ts_pxy_list'):
        pytest.ts_pxy_list = [Track_info(Particle(y=0.0001, E=2.0), 0.00, 0.0001)]
        pytest.ts_pxy_list = track_nturns(lattice, nturns, pytest.ts_pxy_list, nsuperperiods=1, save_track=True, print_progress=False)

    return pytest.ts_pxy_list, nturns


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### TUNE_SHIFT START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### TUNE_SHIFT END ###\n\n\n')
    f.close()

    if hasattr(pytest, 'ts_pxy_list'):
        del pytest.ts_pxy_list


def setup_function(function):
    
    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write(function.__name__)
    f.close()

    pytest.t_start = time()


def teardown_function(function):
    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write(' execution time is ' + '{:.3f}'.format(time() - pytest.t_start) + ' sec\n\n')
    f.close()
    

@pytest.mark.update
def test_update_ref_values(lattice, cmdopt):
    
    update_functions = []
    update_functions.append('test_lattice_transfer_map')
    update_functions.append('test_twiss')
    update_functions.append('test_track_nturns')
    update_functions.append('test_freq_analysis')
    
    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, True)
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
