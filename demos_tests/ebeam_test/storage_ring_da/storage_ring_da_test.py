"""Test of the demo file demos/ebeam/storage_ring_da.py"""

import os
import sys
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from demos_tests.params import *
from storage_ring_da_conf import *
from ocelot.cpbd.chromaticity import *


def test_lattice_transfer_map(lattice, tws=None, update_ref_values=False):
    """R maxtrix test"""

    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_natural_chromaticity(lattice, tws, update_ref_values=False):
    """Natural chromaticity calculation function test"""

    natural_ksi = natural_chromaticity(lattice, tws[0], nsuperperiod=8)
    
    if update_ref_values:
        return numpy2json([natural_ksi])
    
    natural_ksi_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))

    result = check_matrix(natural_ksi, natural_ksi_ref, TOL, assert_info=' natural_chromaticity - ')
    assert check_result(result)


def test_chromaticity(lattice, tws, update_ref_values=False):
    """Chromaticity calculation function test"""
    
    ksi = chromaticity(lattice, tws[0], nsuperperiod=8)
    
    if update_ref_values:
        return numpy2json([ksi])
    
    ksi_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))

    result = check_matrix(ksi, ksi_ref, TOL, assert_info=' chromaticity - ')
    assert check_result(result)


def test_compensate_chromaticity(lattice, tws, update_ref_values=False):
    """Compensate chromaticity function test"""

    match_tolerance = 5.0e-4

    ksi_x, ksi_y, nsuperperiod = compensate_chromaticity_wrapper(lattice)

    ksi = chromaticity(lattice, tws[0], nsuperperiod=nsuperperiod)

    result1 = check_value(ksi[0], ksi_x, match_tolerance, 'absotute', assert_info=' ksi_x - \n')
    result2 = check_value(ksi[1], ksi_y, match_tolerance, 'absotute', assert_info=' ksi_y - \n')
    assert check_result([result1, result2])


def test_create_track_list(lattice, tws=None, update_ref_values=False):
    """Create track list function test"""

    pxy_list = create_track_list_wrapper()
    
    pxy_list = obj2dict(pxy_list, unpack=['particle'])
    
    if update_ref_values:
        return pxy_list
    
    pxy_list_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result = check_dict(pxy_list, pxy_list_ref, TOL, assert_info=' pxy_list - ')
    assert check_result(result)


#@pytest.mark.skip(reason='TOO LONG')
def test_track_nturns(lattice, tws, update_ref_values=False):
    """Track N turns function test"""

    pxy_list, nturns = track_nturns_wrapper(lattice)

    da = [{'turns': pxy.turn} for pxy in pxy_list]
    
    if update_ref_values:
        return da

    da_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result = check_dict(da, da_ref, TOL, assert_info=' da - ')
    assert check_result(result)


#@pytest.mark.skip(reason='TOO LONG')
def test_freq_analysis(lattice, tws, update_ref_values=False):
    """Frequency analysis function test"""

    pxy_list, nturns = track_nturns_wrapper(lattice)    
    pxy_list = freq_analysis(pxy_list, lattice, nturns, harm=True)

    da_mu = [{'turns': [pxy.mux, pxy.muy]} for pxy in pxy_list]

    if update_ref_values:
        return da_mu
    
    da_mu_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result = check_dict(da_mu, da_mu_ref, TOL, assert_info=' da_mu - ')
    assert check_result(result)


def compensate_chromaticity_wrapper(lattice):

    ksi_x = 0.0
    ksi_y = 0.0
    nsuperperiod = 8

    if not hasattr(pytest, 'srda_comp_chromaticity'):
        compensate_chromaticity(lattice, ksi_x_comp=ksi_x, ksi_y_comp=ksi_y, nsuperperiod=nsuperperiod)
        pytest.srda_comp_chromaticity = True

    return ksi_x, ksi_y, nsuperperiod


def create_track_list_wrapper():

    if not hasattr(pytest, 'srda_pxy_list'):
        nx = 100
        ny = 80

        x_array = np.linspace(-0.03, 0.03, nx)
        y_array = np.linspace(0.0001, 0.03, ny)

        pytest.srda_pxy_list = create_track_list(x_array, y_array, p_array=[0.0])
    
    return pytest.srda_pxy_list


def track_nturns_wrapper(lattice):

    nturns = 1000

    if not hasattr(pytest, 'srda_istracked'):
        ksi_x, ksi_y, nsuperperiods = compensate_chromaticity_wrapper(lattice)
        pxy_list = create_track_list_wrapper()
        
        pytest.srda_pxy_list = track_nturns(lattice, nturns, pxy_list, nsuperperiods=nsuperperiods, save_track=True, print_progress=False)
        pytest.srda_istracked = True

    return pytest.srda_pxy_list, nturns


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### STORAGE_RING_DA START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### STORAGE_RING_DA END ###\n\n\n')
    f.close()

    if hasattr(pytest, 'srda_comp_chromaticity'):
        del pytest.srda_comp_chromaticity

    if hasattr(pytest, 'srda_pxy_list'):
        del pytest.srda_pxy_list

    if hasattr(pytest, 'srda_istracked'):
        del pytest.srda_istracked


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
def test_update_ref_values(lattice, tws, cmdopt):
    
    update_functions = []
    update_functions.append('test_lattice_transfer_map')
    update_functions.append('test_natural_chromaticity')
    update_functions.append('test_chromaticity')
    update_functions.append('test_compensate_chromaticity')
    update_functions.append('test_create_track_list')
    update_functions.append('test_track_nturns')
    update_functions.append('test_freq_analysis')
    
    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, tws, True)
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
