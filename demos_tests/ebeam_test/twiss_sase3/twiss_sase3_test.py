"""Test of the demo file demos/ebeam/twiss_sase3.py"""

import os
import sys
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from demos_tests.params import *
from twiss_sase3_conf import *


def test_lattice_transfer_map_before_rematch(lattice, beam=None, update_ref_values=False):
    """R maxtrix test"""
    
    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix before_rematch - ')
    assert check_result(result)


def test_lattice_transfer_map_after_rematch(lattice, beam, update_ref_values=False):
    """R maxtrix test"""

    rematch_wrapper(lattice, beam)

    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix after_rematch - ')
    assert check_result(result)


def test_twiss_after_rematch(lattice, beam, update_ref_values=False):
    """Twiss parameters calculation function test"""

    rematch_wrapper(lattice, beam)
        
    tws = twiss(lattice, Twiss(beam), nPoints=1000)
    
    tws = obj2dict(tws)
    
    if update_ref_values:
        return tws

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    assert check_dict(tws, tws_ref, TOL, 'absotute')
    

def rematch(beta_mean, l_fodo, qdh, lat, extra_fodo, beam, qf, qd):
    """requires l_fodo to be defined in the lattice"""

    k, betaMin, betaMax, __ = fodo_parameters(betaXmean=beta_mean, L=l_fodo, verbose=True)
    k1 = k[0] / qdh.l
    tw0 = Twiss(beam)

    extra = MagneticLattice(extra_fodo)
    tws = twiss(extra, tw0)
    tw2 = tws[-1]

    tw2m = Twiss(tw2)
    tw2m.beta_x = betaMin[0]
    tw2m.beta_y = betaMax[0]
    tw2m.alpha_x = 0.0
    tw2m.alpha_y = 0.0
    tw2m.gamma_x = (1 + tw2m.alpha_x * tw2m.alpha_x) / tw2m.beta_x
    tw2m.gamma_y = (1 + tw2m.alpha_y * tw2m.alpha_y) / tw2m.beta_y

    qf.k1 = k1
    qd.k1 = -k1
    qdh.k1 = -k1

    lat.update_transfer_maps()
    extra.update_transfer_maps()

    R1 = lattice_transfer_map(extra, beam.E)
    Rinv = np.linalg.inv(R1)
    m1 = TransferMap()
    m1.R = lambda energy: Rinv
    tw0m = m1.map_x_twiss(tw2m)
    
    beam.beta_x, beam.alpha_x = tw0m.beta_x, tw0m.alpha_x
    beam.beta_y, beam.alpha_y = tw0m.beta_y, tw0m.alpha_y


def rematch_wrapper(lattice, beam):

    if not hasattr(pytest, 'twiss_sase_isrematch'):
        rematch(19.0, l_fodo, qdh, lattice, extra_fodo, beam, qf, qd)
        lattice.update_transfer_maps()
        pytest.twiss_sase_isrematch = True


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### TWISS_SASE3 START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### TWISS_SASE3 END ###\n\n\n')
    f.close()

    if hasattr(pytest, 'twiss_sase_isrematch'):
        del pytest.twiss_sase_isrematch


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
def test_update_ref_values(lattice, beam, cmdopt):
    
    update_functions = []
    update_functions.append('test_lattice_transfer_map_before_rematch')
    update_functions.append('test_lattice_transfer_map_after_rematch')
    update_functions.append('test_twiss_after_rematch')
    
    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, beam, True)
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
