"""Test of the demo file demos/ebeam/dogleg.py"""

import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from demos_tests.params import *
from dogleg_conf import *

np.random.seed(1)

def test_lattice_transfer_map(lattice, parameter=None, update_ref_values=False):
    """R matrix calculation test"""

    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)

    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)




@pytest.mark.parametrize('parameter', [0, 1, 2])
def test_track_with_energy_shift(lattice, parameter, update_ref_values=False):

    """
    Function generates the particleArray, matches the beam at the entrance of the dogleg and tracks it though lattice
    parameter = 0 - shift -0.5%
    parameter = 1 - shift 0.%
    parameter = 2 - shift +0.5%
    """


    energy_shift = (-1 + parameter) * 0.05
    # generation of the ParticleArray
    emitt = 3.9308e-09
    nparticles = 5000
    np.random.seed(1)
    p_array = generate_parray(sigma_x=np.sqrt(emitt*tws0.beta_x), sigma_px=np.sqrt(emitt*tws0.gamma_x),
                               energy=0.13, nparticles=nparticles)

    # introduce energy shift
    p_array.p()[:] += energy_shift

    # Beam transform for the matching the particleArray at the entrance of the DL
    bt = BeamTransform(tws0)
    navi = Navigator(lattice)
    navi.add_physics_proc(bt, lattice.sequence[0], lattice.sequence[0])

    # tracking
    pytest.tws_track, pytest.p_array = track(lattice, p_array, navi=navi, calc_tws=True, print_progress=False)


    tws_track = obj2dict(pytest.tws_track)
    p = obj2dict(pytest.p_array)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')


    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    if parameter == 1:
        result1 = [None]
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result1+result2)


@pytest.mark.parametrize('parameter', [0, 1, 2])
def test_track_with_energy_shift_tilted(lattice, parameter, update_ref_values=False):
    """
    Function generates the particleArray, matches the beam at the entrance of the dogleg and tracks it though lattice
    parameter = 0 - shift -0.5%
    parameter = 1 - shift 0.%
    parameter = 2 - shift +0.5%
    """

    for elem in lattice.sequence:
        if elem.__class__ in [Sextupole, SBend]:
            elem.tilt = 0
        if elem.__class__ == Quadrupole:
            elem.k1 *= -1
    lattice.update_transfer_maps()


    energy_shift = (-1 + parameter) * 0.05
    # generation of the ParticleArray
    emitt = 3.9308e-09
    nparticles = 5000
    np.random.seed(1)
    p_array = generate_parray(sigma_x=np.sqrt(emitt * tws0.beta_x), sigma_px=np.sqrt(emitt * tws0.gamma_x),
                              energy=0.13, nparticles=nparticles)

    # introduce energy shift
    p_array.p()[:] += energy_shift

    # Beam transform for the matching the particleArray at the entrance of the DL
    bt = BeamTransform(tws0)
    navi = Navigator(lattice)
    navi.add_physics_proc(bt, lattice.sequence[0], lattice.sequence[0])

    # tracking
    pytest.tws_track, pytest.p_array = track(lattice, p_array, navi=navi, calc_tws=True, print_progress=False)

    tws_track = obj2dict(pytest.tws_track)
    p = obj2dict(pytest.p_array)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter)+ '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    if parameter == 1:
        result1 = [None]
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result(result1 + result2)


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### DOGLEG START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### DOGLEG END ###\n\n\n')
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
    update_functions.append('test_track_with_energy_shift')
    update_functions.append('test_track_with_energy_shift_tilted')

    update_function_parameters = {}
    update_function_parameters['test_track_with_energy_shift'] = [0, 1, 2]
    update_function_parameters['test_track_with_energy_shift_tilted'] = [0, 1, 2]
    #update_functions.append('test_track_undulator_with_csr')

    parametr = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parametr:
            #p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
