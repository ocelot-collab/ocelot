"""Test of the demo file ebeam/moga_accelerator_optimization.py"""

import os
import sys
import time

from ocelot.cpbd.optics import *
from ocelot.cpbd.moga.nsga2 import *
from ocelot.cpbd.chromaticity import *

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from params import *
from moga_accelerator_optimization_conf import *

INF = float("inf")


def test_lattice_transfer_map(lattice, update_ref_values=False):
    """R maxtrix test"""

    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_moga_nsga2(lattice, update_ref_values=False):
    """NSGA2 algorithm results test"""

    bounds = ((1.5, 3.5), (-4.0, -2.0), (2.0, 5.0), (-4.0, -2.0), (3.0, 5.0), (-4.0, -1.5))
    weights = (-1.0, -1.0)

    init_pop = []
    #init_pop.append((2.62, -3.1, 2.8, -3.7, 4.0782, -3.534859))
    #init_pop.append((2.62, -3.08, 3.6, -2.5, 3.09, -1.84))
    #init_pop.append((2.865, -3.16, 4.48, -3.74, 4.1, -3.39))

    args = [[Q1, Q2, Q3, Q4, Q5, Q6], lattice]

    res_file = REF_RES_DIR + sys._getframe().f_code.co_name + '.tmp'

    opt = NSGA2(cxpb=0.95, mutpb=0.95)
    opt.set_params(n_gen=100)
    opt.set_params(n_pop=50)
    opt.set_params(bounds=bounds, weights=weights)
    opt.set_params(fit_function=fit_func, fit_function_args=args, init_pop=init_pop)
    opt.set_params(log_print=False, last_data_file=res_file, data_file=False)
    opt.set_params(seed=11)

    opt.start()

    result = json_read(res_file)
    if os.path.exists(res_file):
        os.remove(res_file)

    if update_ref_values:
        return result

    result_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_value(result['iteration'], result_ref['iteration'], TOL, assert_info=' iterations number\n')
    result2 = check_dict(result['data'], result_ref['data'], TOL, assert_info=' data - ')
    assert check_result([result1]+result2)


def fit_func(x0, iter_data, args):

    # update lattice for the new candidate solution
    lattice = args[1]

    vars = args[0]
    for i in range(len(args[0])):
        vars[i].k1 = x0[i]
    lattice.update_transfer_maps()

    beam = Beam()
    beam.E = 2.5

    # additional restriction for a periodic solution existence
    tw = Twiss()
    tw = periodic_twiss(tw, lattice_transfer_map(lattice, beam.E))
    if tw == None:
        return (INF, INF)

    # beam emittance calculation (fitness function parameter 1)
    eb = EbeamParams(lattice, beam, nsuperperiod=6)

    # additional restriction for the beam emittance
    if eb.emittance < 0.0 or eb.emittance > 100.0e-9:
        return (INF, INF)

    # cromaticity compensation
    tws = twiss(lattice)
    ksi = natural_chromaticity(lattice, tws[0], nsuperperiod=6)
    compensate_chromaticity(lattice, ksi_x_comp=0, ksi_y_comp=0, nsuperperiod=6)

    # dynamic aperture calculation (fitness function parameter 2)
    nturns = 5
    nx = 100
    ny = 50

    x_array = np.linspace(-0.1, 0.1, nx)
    y_array = np.linspace(0.0001, 0.05, ny)

    pxy_list = create_track_list(x_array, y_array, p_array=[0.])
    pxy_list = track_nturns(lattice, nturns, pxy_list, nsuperperiods=6, save_track=True)

    da = np.array([pxy.turn for pxy in pxy_list])

    da_num = 1.0
    for i in da:
        if i >= nturns-1:
            da_num += 1.0

    # additional restriction for dynamic aperture
    if da_num < 100.0:
        return (eb.emittance, INF)

    return (eb.emittance, 1.0/da_num)


def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### MOGA_ACCELERATOR_OPTIMIZATION START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### MOGA_ACCELERATOR_OPTIMIZATION END ###\n\n\n')
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
    update_functions.append('test_moga_nsga2')
    
    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, True)
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
