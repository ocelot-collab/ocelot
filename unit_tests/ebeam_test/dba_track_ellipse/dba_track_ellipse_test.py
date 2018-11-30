"""Test of the demo file demos/ebeam/dba_track_ellipse.py"""

import os
import sys
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from dba_track_ellipse_conf import *


def test_lattice_transfer_map(lattice, parametr=None, update_ref_values=False):
    """R maxtrix test"""

    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


@pytest.mark.parametrize('parametr', [0, 1])
def test_tracking_step(lattice, parametr, update_ref_values=False):
    """Tracking step function test
    :parametr=0 - tracking with sextupoles
    :parametr=1 - tracking without sextupoles
    """

    t = np.linspace(0.0, 2.0*np.pi, num=100)
    x, xp = 0.1 * np.cos(t), 0.1 * np.sin(t)

    plist = []
    for xi, xpi in zip(x, xp):
        plist.append(Particle(x=xi, px=xpi))

    if parametr == 1:
        for element in lattice.sequence:
            if element.__class__ == Sextupole:
                element.k2 = 0.0
        lattice.update_transfer_maps()

    navi = Navigator(lattice)
    dz = 10.0

    tracking_step(lattice, plist, dz=dz, navi=navi)
    
    plist = obj2dict(plist)
    
    if update_ref_values:
        return plist

    plist_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parametr) +'.json')

    result = check_dict(plist, plist_ref, TOL, assert_info=' plist - ')
    assert check_result(result)


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### DBA_TRACK_ELLIPSE START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### DBA_TRACK_ELLIPSE END ###\n\n\n')
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
    update_functions.append('test_tracking_step')
    
    update_function_parameters = {}
    update_function_parameters['test_tracking_step'] = [0, 1]
    
    parametr = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parametr:
            result = eval(cmdopt)(lattice, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
