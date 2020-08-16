"""Test of lattice save function in adaptors/elegant_lattice_converter.py file"""

import os
import sys
import time

from ocelot.adaptors.elegant_lattice_converter import *

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from elegant_lattice_conf import *


def test_elegant2ocelot(tws0, method, update_ref_values=False):
    """elegant2ocelot convertion function test"""
    
    SC = ElegantLatticeConverter()
    cell = SC.elegant2ocelot(REF_RES_DIR + 'XFEL_elegant_TD1_S2E.lte')

    lattice = MagneticLattice(cell, method=method)

    r_matrix = lattice_transfer_map(lattice, tws0.E)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_ocelot2elegant(tws0, method, update_ref_values=False):
    """elegant2ocelot + ocelot2elegant + elegant2ocelot converters test"""
    
    # convert Elegant -> Ocelot
    SC = ElegantLatticeConverter()
    cell = SC.elegant2ocelot(REF_RES_DIR + 'XFEL_elegant_TD1_S2E.lte')

    lattice = MagneticLattice(cell, method=method)

    # convert Ocelot -> Elegant
    tmp_file_name = REF_RES_DIR+'tmp.lte'
    SC = ElegantLatticeConverter()
    SC.ocelot2elegant(lattice, file_name=tmp_file_name)

    # convert Elegant -> Ocelot
    SC = ElegantLatticeConverter()
    cell = SC.elegant2ocelot(tmp_file_name)
    if os.path.isfile(tmp_file_name):
        os.remove(tmp_file_name)

    lattice = MagneticLattice(cell, method=method)

    r_matrix = lattice_transfer_map(lattice, tws0.E)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_elegant2ocelot_flash(tws0, method, update_ref_values=False):
    """elegant2ocelot convertion function test"""

    SC = ElegantLatticeConverter()
    cell = SC.elegant2ocelot(REF_RES_DIR + 'FLASH1.lat')

    lattice = MagneticLattice(cell, method=method)

    r_matrix = lattice_transfer_map(lattice, tws0.E)

    if update_ref_values:
        return numpy2json(r_matrix)

    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))

    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_ocelot2elegant_flash(tws0, method, update_ref_values=False):
    """elegant2ocelot + ocelot2elegant + elegant2ocelot converters test"""

    # convert Elegant -> Ocelot
    SC = ElegantLatticeConverter()
    cell = SC.elegant2ocelot(REF_RES_DIR + 'FLASH1.lat')

    lattice = MagneticLattice(cell, method=method)

    # convert Ocelot -> Elegant
    tmp_file_name = REF_RES_DIR + 'tmp1.lte'
    SC = ElegantLatticeConverter()
    SC.ocelot2elegant(lattice, file_name=tmp_file_name)

    # convert Elegant -> Ocelot
    SC = ElegantLatticeConverter()
    cell = SC.elegant2ocelot(tmp_file_name)
    if os.path.isfile(tmp_file_name):
        os.remove(tmp_file_name)

    lattice = MagneticLattice(cell, method=method)

    r_matrix = lattice_transfer_map(lattice, tws0.E)

    if update_ref_values:
        return numpy2json(r_matrix)

    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))

    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)

def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### Elegant lattice converter Lattice START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### Elegant lattice converter END ###\n\n\n')
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
def test_update_ref_values(tws0, method, cmdopt):
    
    update_functions = []
    update_functions.append('test_elegant2ocelot')
    update_functions.append('test_ocelot2elegant')
    update_functions.append('test_elegant2ocelot_flash')
    update_functions.append('test_ocelot2elegant_flash')

    if cmdopt in update_functions:
        result = eval(cmdopt)(tws0, method, True)
        if result is None:
            return
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
