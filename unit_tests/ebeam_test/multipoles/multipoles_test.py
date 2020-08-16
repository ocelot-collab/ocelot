"""Test of the demo file demos/ebeam/multipoles.py"""

import os
import sys
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from multipoles_conf import *
from ocelot.cpbd.chromaticity import *


def test_lattice_transfer_map(lattice, update_ref_values=False):
    """R maxtrix test"""

    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_natural_chromaticity(lattice, update_ref_values=False):
    """Natural chromaticity calculation function test"""

    tws = twiss(lattice)
    natural_ksi = natural_chromaticity(lattice, tws[0])
    
    if update_ref_values:
        return numpy2json([natural_ksi])
    
    natural_ksi_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))

    result = check_matrix(natural_ksi, natural_ksi_ref, TOL, assert_info=' natural chromaticity - ')
    assert check_result(result)


def test_tracking_step(lattice, update_ref_values=False):
    """Tracking step function test"""

    compensate_chromaticity_wrapper(lattice)

    p1 = Particle(x=0.001)
    p2 = Particle(x=-0.0001)

    navi = Navigator(lattice)
    dz = 1.0
    P1 = []
    P2 = []
    for i in range(int(lattice.totalLen/dz)):
        tracking_step(lattice, [p1, p2], dz=dz, navi=navi)
        P1.append(copy(p1))
        P2.append(copy(p2))

    P1 = obj2dict(P1)
    P2 = obj2dict(P2)

    if update_ref_values:
        return [P1, P2]

    P_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(P1, P_ref[0], TOL, 'absotute', assert_info=' P1 - ')
    result2 = check_dict(P2, P_ref[1], TOL, 'absotute', assert_info=' P2 - ')
    assert check_result(result1+result2)


def test_create_track_list(lattice, update_ref_values=False):
    """Create track list function test"""

    track_list = create_track_list_wrapper()
    
    track_list = obj2dict(track_list, unpack=['particle'])
    
    if update_ref_values:
        return track_list
    
    track_list_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result = check_dict(track_list, track_list_ref, TOL, 'absotute', assert_info=' track_list - ')
    assert check_result(result)


#@pytest.mark.skip(reason='TOO LONG')
def test_track_nturns(lattice, update_ref_values=False):
    """Track N turns function test"""

    nturns = 1000
    compensate_chromaticity_wrapper(lattice)
        
    track_list = create_track_list_wrapper()
    track_list = track_nturns(lattice, nturns, track_list, save_track=True, print_progress=False)
    track_list_stable = stable_particles(track_list, nturns)

    track_list_stable = obj2dict(track_list_stable, unpack=['particle'])
    
    p_list = []
    for p in track_list_stable:
        tmp = []
        for i in p['p_list']:
            if isinstance(i, np.ndarray):
                tmp.append(i.tolist())
            else:
                tmp.append(i)
        p_list.append(tmp)

    if update_ref_values:
        return numpy2json(p_list)

    p_list_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(p_list, p_list_ref, TOL, assert_info=' p_list - ')
    assert check_result(result)


def compensate_chromaticity_wrapper(lattice):
    
    if not hasattr(pytest, 'multipoles_comp_chromaticity'):
        compensate_chromaticity(lattice, ksi_x_comp=0.0, ksi_y_comp=0.0)
        pytest.multipoles_comp_chromaticity = True


def create_track_list_wrapper():
    
    if not hasattr(pytest, 'multipoles_track_list'):
        x_array = np.linspace(0.4, 0.60, num=100)
        pytest.multipoles_track_list = create_track_list(x_array, [0.0], [0.0], energy=0.0)

    return pytest.multipoles_track_list


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### MULTIPOLES START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### MULTIPOLES END ###\n\n\n')
    f.close()

    if hasattr(pytest, 'multipoles_comp_chromaticity'):
        del pytest.multipoles_comp_chromaticity

    if hasattr(pytest, 'multipoles_track_list'):
        del pytest.multipoles_track_list


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
def test_update_ref_values(lattice, cmdopt):
    
    update_functions = []
    update_functions.append('test_lattice_transfer_map')
    update_functions.append('test_natural_chromaticity')
    update_functions.append('test_tracking_step')
    update_functions.append('test_create_track_list')
    update_functions.append('test_track_nturns')

    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, True)
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
