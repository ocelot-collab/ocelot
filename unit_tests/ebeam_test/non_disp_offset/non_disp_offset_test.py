"""Test of the demo file demos/ebeam/non_disp_offset.py"""

import os
import sys
from copy import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from non_disp_offset_conf import *


def test_lattice_transfer_map(lattice, parametr=None, update_ref_values=False):
    """R maxtrix test"""

    r_matrix = lattice_transfer_map(lattice[0], 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)
    

def test_twiss(lattice, parametr=None, update_ref_values=False):
    """Twiss parameters calculation function test"""

    tw0 = Twiss()
    tw0.beta_x = 5.0
    tw0.alpha_x = 1.4
    tw0.beta_y = 16.0
    tw0.alpha_y = 5.1

    tws = twiss(lattice[0], tw0, nPoints=500)
    
    tws = obj2dict(tws)
    
    if update_ref_values:
        return tws

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws - ')
    assert check_result(result)


@pytest.mark.parametrize('parametr', [0, 1])
def test_tracking_step_energy_offset(lattice, parametr, update_ref_values=False):
    """Tracking step function test
    :parametr=0 - tracking with default MethodTM()
    :parametr=1 - tracking without  MethodTM() - global_method = SecondTM
    """

    p = Particle(x=0.0, p=0.02)

    navi = Navigator(lattice[parametr])
    dz = 0.01
    P = [copy.deepcopy(p)]

    n_end = int(lattice[parametr].totalLen/dz)
    for iii in range(n_end):
        tracking_step(lattice[parametr], [p], dz=dz, navi=navi)
        P.append(copy.deepcopy(p))
    
    P = obj2dict(P)
    
    if update_ref_values:
        return P

    p_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parametr) +'.json')

    result = check_dict(P, p_ref, TOL, 'absotute', assert_info=' P - ')
    assert check_result(result)


@pytest.mark.parametrize('parametr', [0, 1])
def test_tracking_step_phase_space(lattice, parametr, update_ref_values=False):
    """Tracking step function test
    :parametr=0 - tracking with default MethodTM()
    :parametr=1 - tracking without  MethodTM() - global_method = SecondTM
    """
    
    t = np.linspace(0.0, 2.0*pi, num=100)
    x, xp = 0.1*np.cos(t), 0.1*np.sin(t)
    plist = []
    for xi, xpi in zip(x, xp):
        plist.append(Particle(x=xi, px=xpi))

    navi = Navigator(lattice[parametr])
    dz = lattice[parametr].totalLen

    tracking_step(lattice[parametr], plist, dz=dz, navi=navi)

    plist = obj2dict(plist)
    
    if update_ref_values:
        return plist

    plist_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parametr) +'.json')

    result = check_dict(plist, plist_ref, TOL, 'absotute', assert_info=' plist - ')
    assert check_result(result)


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### NON_DISP_OFFCET START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### NON_DISP_OFFCET END ###\n\n\n')
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
def test_update_ref_values(lattice, cmdopt):
    
    update_functions = []
    update_functions.append('test_lattice_transfer_map')
    update_functions.append('test_twiss')
    update_functions.append('test_tracking_step_energy_offset')
    update_functions.append('test_tracking_step_phase_space')
    
    update_function_parameters = {}
    update_function_parameters['test_tracking_step_energy_offset'] = [0, 1]
    update_function_parameters['test_tracking_step_phase_space'] = [0, 1]
    
    parametr = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parametr:
            result = eval(cmdopt)(lattice, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
