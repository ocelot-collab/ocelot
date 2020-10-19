"""Test of the demo file demos/ebeam/csr_ex.py"""
from ocelot.adaptors.astra2ocelot import exact_xxstg_2_xp_mad, exact_xp_2_xxstg_mad
import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *

from io_conf import *




def test_npz(p_array, parameter=None, update_ref_values=False):
    """
    testing applying one marker as start ans stop
    """

    p_array_ref = copy.deepcopy(p_array)

    save_particle_array("test.npz", p_array)
    p_array_reload = load_particle_array("test.npz")

    p_rel = obj2dict(p_array_reload)
    p_ref = obj2dict(p_array_ref)

    result2 = check_dict(p_rel, p_ref, tolerance=TOL, assert_info=' p - ')


    assert check_result(result2 )


def test_ast_mad_transf(p_array, parameter=None, update_ref_values=False):
    """
    testing applying one marker as start ans stop
    """


    gamref = p_array.E / m_e_GeV
    P = p_array.rparticles.view()
    xp = exact_xxstg_2_xp_mad(P, gamref)

    xxstg = exact_xp_2_xxstg_mad(xp, gamref)



    result2 =  check_matrix(xxstg.T, P, tolerance=1.0e-10, tolerance_type='relative', assert_info='rp')
    assert check_result(result2 )


def test_ast(p_array, parameter=None, update_ref_values=False):
    """
    testing applying one marker as start ans stop
    """

    p_array_ref = copy.deepcopy(p_array)
    p_array_ref.rparticles[4, 0] *= 0
    p_array_ref.rparticles[5, 0] *= 0
    save_particle_array("test.ast", p_array_ref)
    p_array_reload = load_particle_array("test.ast")

    p_rel = obj2dict(p_array_reload)
    p_ref = obj2dict(p_array_ref)

    result2 = check_dict(p_rel, p_ref, tolerance=1.0e-7, tolerance_type='absolute', assert_info=' p - ')

    assert check_result(result2 )

def test_ast_new(p_array, parameter=None, update_ref_values=False):
    """
    testing applying one marker as start ans stop
    """

    p_array_ref = copy.deepcopy(p_array)
    #p_array_ref.rparticles[4, 0] *= 0
    #p_array_ref.rparticles[5, 0] *= 0
    save_particle_array("test.ast", p_array_ref)
    p_array_reload = load_particle_array("test.ast")

    p_rel = obj2dict(p_array_reload)
    p_ref = obj2dict(p_array_ref)

    result2 = check_dict(p_rel, p_ref, tolerance=1.0e-7, tolerance_type='absolute', assert_info=' p - ')

    assert check_result(result2 )

def test_fmt1(p_array, parameter=None, update_ref_values=False):
    """
    testing applying one marker as start ans stop
    """

    p_array_ref = copy.deepcopy(p_array)
    p_array_ref.rparticles[4, 0] *= 0
    p_array_ref.rparticles[5, 0] *= 0
    save_particle_array("test.fmt1", p_array_ref)
    p_array_reload = load_particle_array("test.fmt1")
    print(np.array_equal(p_array_reload.rparticles, p_array_ref.rparticles))
    print(p_array_reload.rparticles[5, :10], p_array_ref.rparticles[5, :10])
    p_rel = obj2dict(p_array_reload)
    p_ref = obj2dict(p_array_ref)

    result2 = check_dict(p_rel, p_ref, tolerance=1.0e-7, tolerance_type='absolute', assert_info=' p - ')

    assert check_result(result2 )

def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### PHYS PROC START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### PHYS PROC END ###\n\n\n')
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
def test_update_ref_values(p_array, cmdopt):
    
    update_functions = []
    #update_functions.append('test_generate_parray')

    update_function_parameters = {}
    #update_function_parameters['test_track_smooth'] = [0, 1]

    parameter = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parameter:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(p_arr, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
