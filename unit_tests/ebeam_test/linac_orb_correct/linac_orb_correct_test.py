"""Test of the demo file demos/ebeam/ring_orb_correct.py"""

import os
import sys
import time

from ocelot.cpbd.orbit_correction import *
from ocelot.cpbd.response_matrix import *

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from linac_orb_correct_conf import *


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

    tws0 = Twiss()
    tws0.beta_x = 3.7650634043210496
    tws0.beta_y = 0.8380337780945593
    tws0.alpha_x = 2.071053894766802
    tws0.alpha_y = -0.61098754078337
    tws0.gamma_x = 1.404827400504972
    tws0.gamma_y = 1.6387236539737111
    tws0.E = 0.13
    tws = twiss(lattice, tws0)
    
    tws = obj2dict(tws)
    
    if update_ref_values:
        return tws

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws - ')
    assert check_result(result)


def test_response_matrix(lattice, update_ref_values=False):
    """Responce maxtrix calculation test"""

    orb = Orbit(lattice)
    linac_method = LinacRmatrixRM(lattice=orb.lat, hcors=orb.hcors, vcors=orb.vcors, bpms=orb.bpms)

    orb.response_matrix = ResponseMatrix(method=linac_method)
    orb.response_matrix.calculate()

    if update_ref_values:
        return numpy2json(orb.response_matrix.matrix)

    response_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(orb.response_matrix.matrix, response_matrix_ref, TOL, assert_info=' response_matrix - ')
    assert check_result(result)


def test_save_load_rm(lattice, update_ref_values=False):
    """Responce maxtrix calculation test"""

    orb = Orbit(lattice)
    linac_method = LinacRmatrixRM(lattice=orb.lat, hcors=orb.hcors, vcors=orb.vcors, bpms=orb.bpms)

    orb.response_matrix = ResponseMatrix(method=linac_method)
    orb.response_matrix.calculate()
    response_matrix_ref = np.copy(orb.response_matrix.get_matrix())
    orb.response_matrix.dump(REF_RES_DIR + "test" + '.p')
    orb.response_matrix.load(REF_RES_DIR + "test" + '.p')
    #if update_ref_values:
    #    return numpy2json(orb.response_matrix.matrix)

    #response_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))

    result = check_matrix(orb.response_matrix.get_matrix(), response_matrix_ref, TOL, assert_info=' response_matrix - ')
    assert check_result(result)

def test_inject_extract_rm(lattice, update_ref_values=False):
    """Responce maxtrix calculation test"""

    orb = Orbit(lattice)
    linac_method = LinacRmatrixRM(lattice=orb.lat, hcors=orb.hcors, vcors=orb.vcors, bpms=orb.bpms)

    orb.response_matrix = ResponseMatrix(method=linac_method)
    orb.response_matrix.calculate()

    response_matrix_ref = np.copy(orb.response_matrix.get_matrix())

    rm = orb.response_matrix.extract(cor_list=['CIX.78.I1', 'CIY.80.I1'], bpm_list=['BPMA.85.I1', 'BPMA.87.I1'])
    rm[0, 0] = 2
    rm[2, 0] = 0.1
    rm[1, 1] = -0.1
    rm[3, 1] = -1.3
    response_matrix_ref = np.copy(rm)
    orb.response_matrix.inject(cor_list=['CIX.78.I1', 'CIY.80.I1'], bpm_list=['BPMA.85.I1', 'BPMA.87.I1'],
                               inj_matrix=rm)

    rm = orb.response_matrix.extract(cor_list=['CIX.78.I1', 'CIY.80.I1'], bpm_list=['BPMA.85.I1', 'BPMA.87.I1'])

    result = check_matrix(rm, response_matrix_ref, TOL, assert_info=' response_matrix - ')
    assert check_result(result)

def test_correction(lattice, update_ref_values=False):
    """Orbit correction test"""

    orb = Orbit(lattice)
    linac_method = LinacRmatrixRM(lattice=orb.lat, hcors=orb.hcors, vcors=orb.vcors, bpms=orb.bpms)

    orb.response_matrix = ResponseMatrix(method=linac_method)
    orb.response_matrix.calculate()
    
    x_bpm_b, y_bpm_b, x_bpm, y_bpm = correction_wrapper(orb, linac_method)
    vcors_angle = np.array([cor.angle for cor in orb.vcors])
    hcors_angle = np.array([cor.angle for cor in orb.hcors])

    if update_ref_values:
        return numpy2json([x_bpm_b, y_bpm_b, x_bpm, y_bpm, hcors_angle, vcors_angle])

    bpm_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))

    result1 = check_matrix(x_bpm_b, bpm_ref[0], TOL, assert_info=' x_bpm_before_correction - ')
    result2 = check_matrix(y_bpm_b, bpm_ref[1], TOL, assert_info=' y_bpm_before_correction - ')
    result3 = check_matrix(x_bpm, bpm_ref[2], TOL, assert_info=' x_bpm - ')
    result4 = check_matrix(y_bpm, bpm_ref[3], TOL, assert_info=' y_bpm - ')
    result5 = check_matrix(hcors_angle, bpm_ref[4], TOL, assert_info=' hor corrector - ')
    result6 = check_matrix(vcors_angle, bpm_ref[5], TOL, assert_info=' ver corrector - ')
    assert check_result(result1+result2+result3+result4 + result5 + result6)


def test_correction_micado(lattice, update_ref_values=False):
    """Orbit correction test"""
    for e in lattice.sequence:
        if e.__class__ in [Hcor, Vcor]:
            e.angle = 0
    lattice.update_transfer_maps()
    orb = Orbit(lattice)
    orb.orbit_solver = MICADO(epsilon_x=0.001, epsilon_y=0.001, epsilon_ksi=1e-5)

    linac_method = LinacRmatrixRM(lattice=orb.lat, hcors=orb.hcors, vcors=orb.vcors, bpms=orb.bpms)

    orb.response_matrix = ResponseMatrix(method=linac_method)
    orb.response_matrix.calculate()


    x_bpm_b, y_bpm_b = linac_method.read_virtual_orbit(p_init=Particle())
    orb.correction(beta=0)
    x_bpm, y_bpm = linac_method.read_virtual_orbit(p_init=Particle())

    print(x_bpm_b, y_bpm_b, x_bpm, y_bpm )
    vcors_angle = np.array([cor.angle for cor in orb.vcors])
    hcors_angle = np.array([cor.angle for cor in orb.hcors])
    if update_ref_values:
        return numpy2json([x_bpm_b, y_bpm_b, x_bpm, y_bpm, hcors_angle, vcors_angle])

    bpm_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))

    result1 = check_matrix(x_bpm_b, bpm_ref[0], TOL, assert_info=' x_bpm_before_correction - ')
    result2 = check_matrix(y_bpm_b, bpm_ref[1], TOL, assert_info=' y_bpm_before_correction - ')
    result3 = check_matrix(x_bpm, bpm_ref[2], TOL, assert_info=' x_bpm - ')
    result4 = check_matrix(y_bpm, bpm_ref[3], TOL, assert_info=' y_bpm - ')
    result5 = check_matrix(hcors_angle, bpm_ref[4], TOL, assert_info=' hor corrector - ')
    result6 = check_matrix(vcors_angle, bpm_ref[5], TOL, assert_info=' ver corrector - ')
    assert check_result(result1+result2+result3+result4 + result5 + result6)


def test_lattice_track(lattice, update_ref_values=False):
    """Lattice track function test"""
    #qi_74_i1.dx = 0.001
    #qi_74_i1.dy = 0.001
    for e in lattice.sequence:
        if e.__class__ in [Hcor, Vcor]:
            e.angle = 0.
        if e.__class__ is Quadrupole:
            #if e.id is 'QI.74.I1':
            #    print("here")
            #    e.dx = 0.001
            #    e.dy = 0.001
            print(e.id, e.dx, e.dy)
    lattice.update_transfer_maps()

    orb = Orbit(lattice)
    linac_method = LinacRmatrixRM(lattice=orb.lat, hcors=orb.hcors, vcors=orb.vcors, bpms=orb.bpms)

    orb.response_matrix = ResponseMatrix(method=linac_method)
    orb.response_matrix.calculate()
    x_bpm_b, y_bpm_b = linac_method.read_virtual_orbit(p_init=Particle())
    #correction_wrapper(orb, linac_method)
    orb.correction(beta=0)
    p_list = lattice_track(lattice, Particle())
    p_list = obj2dict(p_list, unpack=['particle'])

    if update_ref_values:
        return p_list

    p_list_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result = check_dict(p_list, p_list_ref, TOL, assert_info=' p_list - ')
    assert check_result(result)


def correction_wrapper(orb, correction_method):
    
    x_bpm_b, y_bpm_b = correction_method.read_virtual_orbit(p_init=Particle())
    x_bpm, y_bpm = None, None
     
    if not hasattr(pytest, 'roc_correction'):
        orb.correction(beta=0)
        x_bpm, y_bpm = correction_method.read_virtual_orbit(p_init=Particle())
        pytest.roc_correction = True
    print("correction wrapper: ", x_bpm, y_bpm)
    return x_bpm_b, y_bpm_b, x_bpm, y_bpm

def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### DBA START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### DBA END ###\n\n\n')
    f.close()
    
    if hasattr(pytest, 'roc_correction'):
        del pytest.roc_correction


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
    update_functions.append('test_response_matrix')
    update_functions.append('test_correction')
    update_functions.append('test_correction_micado')
    update_functions.append('test_lattice_track')

    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, True)
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
