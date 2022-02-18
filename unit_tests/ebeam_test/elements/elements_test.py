"""Test of the demo file demos/ebeam/multipoles.py"""

import os
import sys
import time
import copy

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from elements_conf import *


def test_lattice_transfer_map(lattice, p_array, parameter=None, update_ref_values=False):
    """R maxtrix test for whole line"""

    r_matrix = lattice_transfer_map(lattice, 0.0)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_tracking_step(lattice, p_array, parameter=None, update_ref_values=False):
    """Tracking step function test"""

    p1 = Particle(x=0.001, y=0.0005, p=0.0001)
    p2 = Particle(x=-0.0001, y=0.0005, p=-0.0001)

    navi = Navigator(lattice)
    dz = 1.0
    P1 = []
    P2 = []
    for i in range(int(lattice.totalLen/dz)):
        tracking_step(lattice, [p1, p2], dz=dz, navi=navi)
        P1.append(copy.copy(p1))
        P2.append(copy.copy(p2))

    P1 = obj2dict(P1)
    P2 = obj2dict(P2)

    if update_ref_values:
        return [P1, P2]

    P_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(P1, P_ref[0], TOL, 'absolute', assert_info=' P1 - ')
    result2 = check_dict(P2, P_ref[1], TOL, 'absolute', assert_info=' P2 - ')
    assert check_result(result1+result2)


@pytest.mark.parametrize('parameter', [0, 1])
def test_track(lattice, p_array, parameter, update_ref_values=False):
    """
    test Runge_Kutta transfer map for undulator

    0 - tracking of the electron beam with positive energy chirp trough undulator
    1 - tracking of the electron beam with negative energy chirp trough undulator
    """

    p_array_track = copy.deepcopy(p_array)

    navi = Navigator(lattice)
    navi.unit_step = 0.05

    smb = EmptyProc()
    if parameter == 1:
        navi.add_physics_proc(smb, lattice.sequence[0], lattice.sequence[0])

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi)

    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json')

    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result( result2)

def test_oct_vs_multip(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test Octupole vs Multipole with the same parameters
    """

    oct = Octupole(l=0.5, k3=300, tm=KickTM)
    mult3 = Multipole(kn=[0, 0, 0, 300 * 0.5], tm=KickTM)

    parray = ParticleArray(n=1)
    parray.x()[:] = 0.001
    parray.y()[:] = -0.001
    parray.p()[:] = 0.001

    lat = MagneticLattice((oct), method={Octupole: KickTM, "nkick": 1})
    track(lat, parray, Navigator(lat))
    p1 = obj2dict(parray)
    parray2 = ParticleArray(n=1)
    parray2.x()[:] = 0.001
    parray2.y()[:] = -0.001
    parray2.p()[:] = 0.001

    lat = MagneticLattice((Drift(l=0.250), mult3, Drift(l=0.25)), method={Multipole: KickTM, "nkick": 1})
    track(lat, parray2, Navigator(lat))

    p2 = obj2dict(parray2)

    result2 = check_dict(p1, p2, tolerance=TOL, assert_info=' p - ')
    assert check_result( result2)

def test_xyquad(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test XYQuadrupole
    """
    qxy = XYQuadrupole(l=0.5, x_offs=0.01, y_offs=0.025, k1=-0.5, eid='QK')
    lat = MagneticLattice((qxy), method={"global": TransferMap})
    track(lat, p_array, Navigator(lat))
    p = obj2dict(p_array)

    if update_ref_values:
        return {'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result2 = check_dict(p, tws_track_p_array_ref['p_array'], tolerance=TOL, assert_info=' p - ')
    assert check_result( result2)

def test_self_xyquad(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test Octupole vs Multipole with the same parameters
    """

    parray = ParticleArray(n=1)
    parray.x()[:] = 0.001
    parray.y()[:] = -0.001
    parray.p()[:] = 0.001

    qxy = XYQuadrupole(l=0.5, x_offs=0.01, y_offs=0.025, k1=-0.5, eid='QK')
    lat = MagneticLattice((qxy), method={"global": TransferMap})
    track(lat, parray, Navigator(lat))
    p1 = obj2dict(parray)

    parray2 = ParticleArray(n=1)
    parray2.x()[:] = 0.001
    parray2.y()[:] = -0.001
    parray2.p()[:] = 0.001

    qxy = XYQuadrupole(l=0.5, x_offs=0.01, y_offs=0.025, k1=-0.5, eid='QK')
    lat = MagneticLattice((qxy), method={"global": SecondTM})
    track(lat, parray2, Navigator(lat))

    p2 = obj2dict(parray2)

    result2 = check_dict(p1, p2, tolerance=TOL, assert_info=' p - ')
    assert check_result( result2)

def test_quad_vs_matrix(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test quadrupole with offsets vs Matrix
    """

    qd = Quadrupole(l=0.2, k1=-1, k2=-20, tm=SecondTM)
    qd.dx, qd.dy = -1e-3, 1e-3

    m = Matrix(tm=SecondTM)
    m.r = qd.R(1)[0]
    m.t = qd.T(1)[0]
    m.b = qd.B(1)[0]
    m.l = qd.l

    parray = ParticleArray(n=2)
    parray.E = 1
    parray.rparticles[0, :] = 0.001
    parray.rparticles[1, :] = -0.0002
    parray.rparticles[2, :] = -0.0005
    parray.rparticles[3, :] = 0.0003
    parray.rparticles[4, :] = [-0.0005, 0.0005]
    parray.rparticles[5, :] = [0.0005, -0.0005]

    parray_1 = copy.deepcopy(parray)
    for t in qd.tms:
        t.apply(parray_1)

    p1 = obj2dict(parray_1)

    parray_2 = copy.deepcopy(parray)
    for t in m.tms:
        t.apply(parray_2)

    p2 = obj2dict(parray_2)

    result2 = check_dict(p1, p2, tolerance=TOL, assert_info=' p - ')
    assert check_result( result2)

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
def test_update_ref_values(lattice, p_array, cmdopt):
    
    update_functions = []
    update_functions.append('test_lattice_transfer_map')
    update_functions.append('test_tracking_step')
    update_functions.append('test_track')
    update_functions.append("test_xyquad")

    update_function_parameters = {}
    update_function_parameters['test_track'] = [0, 1]

    parameter = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parameter:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p_arr, p, True)

            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')

            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
