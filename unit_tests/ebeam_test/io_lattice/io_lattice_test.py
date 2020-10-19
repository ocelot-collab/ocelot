"""Test of lattice save function in cpbd/io.py file"""

import os
import sys
import time
import copy

from ocelot.cpbd.io import *

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from io_lattice_conf import *


def test_original_lattice_transfer_map(lattice, tws0, method, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""

    r_matrix = lattice_transfer_map(lattice, tws0.E)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_write_lattice(lattice, tws0, method, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""

    write_lattice(lattice, file_name="tmp_lattice.py")
    import tmp_lattice as tmp
    new_lat = MagneticLattice(tmp.cell, method=lattice.method)

    res = []
    for i, elem in enumerate(lattice.sequence):
        if (elem.__class__.__name__ == new_lat.sequence[i].__class__.__name__  and elem.id== new_lat.sequence[i].id
                and elem.l == new_lat.sequence[i].l):
            res.append(None)
        else:
            res.append(-1)


    assert check_result(res)


def test_write_lattice_w_coupler(lattice, tws0, method, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""
    lattice0 = copy.deepcopy(lattice)
    for elem in lattice0.sequence:
        if elem.id == 'C.A1.1.1.I1':
            elem.vx_up = -5.6813e-5 + 1.0751e-5j
            elem.vy_up = -4.1091e-5 + 5.739e-7j
            elem.vxx_up = 0.00099943 - 0.00081401j
            elem.vxy_up = 0.0034065 - 0.0004146j
            elem.vx_down = -2.4014e-5 + 1.2492e-5j
            elem.vy_down = 3.6481e-5 + 7.9888e-6j
            elem.vxx_down = -0.004057 - 0.0001369j
            elem.vxy_down = 0.0029243 - 1.2891e-5j
    lattice0.update_transfer_maps()
    write_lattice(lattice0, file_name="tmp_lattice.py")
    import tmp_lattice as tmp
    new_lat = MagneticLattice(tmp.cell, method=lattice0.method)

    res = []
    for i, elem in enumerate(lattice0.sequence):
        if (elem.__class__.__name__ == new_lat.sequence[i].__class__.__name__  and elem.id== new_lat.sequence[i].id
                and elem.l == new_lat.sequence[i].l):
            res.append(None)
        else:
            res.append(-1)
    #for elem in lattice.sequence:
    #    if elem.id == 'C.A1.1.1.I1':
    #        elem.vx_up = 0
    #        elem.vy_up = 0
    #        elem.vxx_up = 0
    #        elem.vxy_up = 0
    #        elem.vx_down = 0
    #        elem.vy_down = 0
    #        elem.vxx_down = 0
    #        elem.vxy_down = 0
    lattice.update_transfer_maps()

    assert check_result(res)

def test_original_twiss(lattice, tws0, method, parametr=None, update_ref_values=False):
    """Twiss parameters calculation function test"""

    tws = twiss(lattice, tws0, nPoints=None)
    
    tws = obj2dict(tws)
    
    if update_ref_values:
        return tws

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws - ')
    assert check_result(result)


@pytest.mark.parametrize('parametr', [False, True])
def test_lat2input(lattice, tws0, method, parametr, update_ref_values=False):
    """lat2input with tws0 saving function test"""

    lines_arr = lat2input(lattice, tws0=tws0)
    lines = ''.join(lines_arr)
    
    loc_dict = {}
    try:
        exec(lines, globals(), loc_dict)
    except Exception as err:
        assert check_result(['Exception error during the lattice file execution, parametr is ' + str(parametr)])

    if parametr:
        if "tws0" in loc_dict:
            tws0_new = loc_dict['tws0']
        else:
            assert check_result(['No tws0 in the lattice file, parametr is ' + str(parametr)])
    else:
        tws0_new = tws0

    if 'cell' in loc_dict:
        lattice_new = MagneticLattice(loc_dict['cell'], method=method)

        lattice_new_transfer_map_check(lattice_new, tws0_new)
        twiss_new_check(lattice_new, tws0_new)
    else:
        assert check_result(['No cell variable in the lattice file, parametr is ' + str(parametr)])
    
        
def lattice_new_transfer_map_check(lattice, tws0):

    r_matrix = lattice_transfer_map(lattice, tws0.E)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + 'test_original_lattice_transfer_map.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix for new lattice - ')
    assert check_result(result)


def twiss_new_check(lattice, tws0):

    tws = twiss(lattice, tws0, nPoints=None)
    
    tws = obj2dict(tws)

    tws_ref = json_read(REF_RES_DIR + 'test_original_twiss.json')
    
    result = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws for new lattice - ')
    assert check_result(result)


def test_merger(lattice, tws0, method, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""
    d = Drift(l=0.5)
    q = Quadrupole(l=0.3, k1=3, k2=3.3, eid="quad")
    b = Bend(l=1.2, angle=0.02, k1=-1, k2=.9, eid="bend")
    s = Sextupole(l=0.23, k2=22, eid="sext")
    c = Cavity(l=2, v=0.02, freq=1.3e9, phi=10, eid="cav")
    cor = Hcor(l=0.1, eid="cor")
    sol = Solenoid(l=0.12, k=0.002)
    tds = TDCavity(l=0.5, freq=1.2e6, phi=90, v=0.02, tilt=0.0, eid="tds")
    m = Monitor(eid="mon")
    mat = Matrix(l=0.5, r11=1.1, r22=1, r33=1, r44=1, r55=1, r66=1, r12=0.1, t111=0.8)
    b2 = RBend(l=1.2, angle=0.01, k1=-1, k2=-.9, eid="bend")
    b3 = SBend(l=1.2, angle=0.02, k1=1, k2=.5, eid="bend")

    init_energy = 0.5

    cell = (d, q, b, s, c, cor, sol, tds, m, mat, b2, b3)

    lat = MagneticLattice(cell, method=MethodTM({'global': SecondTM}))

    R = lattice_transfer_map(lat, energy=init_energy)
    new_lat = merger(lat, remaining_types=[], remaining_elems=[], init_energy=init_energy)
    R_new = lattice_transfer_map(new_lat, energy=init_energy)

    result = check_matrix(R, R_new, TOL, assert_info=' r_matrix - ')
    result2 = check_matrix(lat.T, new_lat.T, TOL, assert_info=' t_matrix - ')
    assert check_result(result + result2)


def test_merger_elem(lattice, tws0, method, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""
    d = Drift(l=0.5)
    q = Quadrupole(l=0.3, k1=3, k2=3.3, eid="quad")
    b = Bend(l=1.2, angle=0.02, k1=-1, k2=.9, eid="bend")
    s = Sextupole(l=0.23, k2=22, eid="sext")
    c = Cavity(l=2, v=0.02, freq=1.3e9, phi=10, eid="cav")
    cor = Hcor(l=0.1, eid="cor")
    sol = Solenoid(l=0.12, k=0.002)
    tds = TDCavity(l=0.5, freq=1.2e6, phi=90, v=0.02, tilt=0.0, eid="tds")
    m = Monitor(eid="mon")
    mat = Matrix(l=0.5, r11=1.1, r22=1, r33=1, r44=1, r55=1, r66=1, r12=0.1, t111=0.8)
    b2 = RBend(l=1.2, angle=0.01, k1=-1, k2=-.9, eid="bend")
    b3 = SBend(l=1.2, angle=0.02, k1=1, k2=.5, eid="bend")

    init_energy = 0.5

    cell = (d, q, b, s, c, cor, sol, tds, m, mat, b2, b3)

    lat = MagneticLattice(cell, method=MethodTM({'global': SecondTM}))

    R = lattice_transfer_map(lat, energy=init_energy)
    new_lat = merger(lat, remaining_types=[], remaining_elems=[sol], init_energy=init_energy)
    R_new = lattice_transfer_map(new_lat, energy=init_energy)

    result = check_matrix(R, R_new, TOL, assert_info=' r_matrix - ')
    result2 = check_matrix(lat.T, new_lat.T, TOL, assert_info=' t_matrix - ')
    assert check_result(result + result2)

def test_merger_elem_w_coupler(lattice, tws0, method, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""
    d = Drift(l=0.5)
    q = Quadrupole(l=0.3, k1=3, k2=3.3, eid="quad")
    b = Bend(l=1.2, angle=0.02, k1=-1, k2=.9, eid="bend")
    s = Sextupole(l=0.23, k2=22, eid="sext")
    c = Cavity(l=2, v=0.02, freq=1.3e9, phi=10, eid="cav")
    c.vx_up = -5.6813e-5 + 1.0751e-5j
    c.vy_up = -4.1091e-5 + 5.739e-7j
    c.vxx_up = 0.00099943 - 0.00081401j
    c.vxy_up = 0.0034065 - 0.0004146j
    c.vx_down = -2.4014e-5 + 1.2492e-5j
    c.vy_down = 3.6481e-5 + 7.9888e-6j
    c.vxx_down = -0.004057 - 0.0001369j
    c.vxy_down = 0.0029243 - 1.2891e-5j
    cor = Hcor(l=0.1, eid="cor")
    sol = Solenoid(l=0.12, k=0.002)
    tds = TDCavity(l=0.5, freq=1.2e6, phi=90, v=0.02, tilt=0.0, eid="tds")
    m = Monitor(eid="mon")
    mat = Matrix(l=0.5, r11=1.1, r22=1, r33=1, r44=1, r55=1, r66=1, r12=0.1, t111=0.8)
    b2 = RBend(l=1.2, angle=0.01, k1=-1, k2=-.9, eid="bend")
    b3 = SBend(l=1.2, angle=0.02, k1=1, k2=.5, eid="bend")

    init_energy = 0.5

    cell = (d, q, b, s, c, cor, sol, tds, m, mat, b2, b3)

    lat = MagneticLattice(cell, method=MethodTM({'global': SecondTM}))

    R = lattice_transfer_map(lat, energy=init_energy)
    new_lat = merger(lat, remaining_types=[], remaining_elems=[sol], init_energy=init_energy)
    R_new = lattice_transfer_map(new_lat, energy=init_energy)

    result = check_matrix(R, R_new, TOL, assert_info=' r_matrix - ')
    result2 = check_matrix(lat.T, new_lat.T, TOL, assert_info=' t_matrix - ')
    assert check_result(result + result2)

def test_merger_type(lattice, tws0, method, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""
    d = Drift(l=0.5)
    q = Quadrupole(l=0.3, k1=3, k2=3.3, eid="quad")
    b = Bend(l=1.2, angle=0.02, k1=-1, k2=.9, eid="bend")
    s = Sextupole(l=0.23, k2=22, eid="sext")
    c = Cavity(l=2, v=0.02, freq=1.3e9, phi=10, eid="cav")
    cor = Hcor(l=0.1, eid="cor")
    ap = Aperture()
    sol = Solenoid(l=0.12, k=0.002)
    d2 = Drift(l=0.5)
    tds = TDCavity(l=0.5, freq=1.2e6, phi=90, v=0.02, tilt=0.0, eid="tds")
    m = Monitor(eid="mon")
    mat = Matrix(l=0.5, r11=1.1, r22=1, r33=1, r44=1, r55=1, r66=1, r12=0.1, t111=0.8)
    b2 = RBend(l=1.2, angle=0.01, k1=-1, k2=-.9, eid="bend")
    b3 = SBend(l=1.2, angle=0.02, k1=1, k2=.5, eid="bend")

    init_energy = 0.5

    cell = (d, q, b, s, c, cor, ap, sol,d2, tds, m, mat, b2, b3)

    lat = MagneticLattice(cell, method=MethodTM({'global': SecondTM}))

    R = lattice_transfer_map(lat, energy=init_energy)
    new_lat = merger(lat, remaining_types=[Drift], remaining_elems=[sol], init_energy=init_energy)
    R_new = lattice_transfer_map(new_lat, energy=init_energy)

    result = check_matrix(R, R_new, tolerance=1.0e-12, tolerance_type='absolute', assert_info=' r_matrix - ')
    result2 = check_matrix(lat.T, new_lat.T, tolerance=1.0e-12, tolerance_type='absolute', assert_info=' t_matrix - ')
    assert check_result(result + result2)


def test_merger_extensive(lattice, tws0, method, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""


    init_energy = 0.05


    R = lattice_transfer_map(lattice, energy=init_energy)
    new_lat = merger(lattice, remaining_types=[Hcor, Vcor, Monitor], remaining_elems=[MPBPMF_47_I1, START_96_I1], init_energy=init_energy)
    R_new = lattice_transfer_map(new_lat, energy=init_energy)

    result = check_matrix(R, R_new, TOL, assert_info=' r_matrix - ')
    result2 = check_matrix(lattice.T, new_lat.T, TOL, assert_info=' t_matrix - ')
    assert check_result(result + result2)


def test_merger_tilt(lattice, tws0, method, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""
    d = Drift(l=0.5)
    q = Quadrupole(l=0.3, k1=3, k2=3.3, eid="quad", tilt=1.)
    b = Bend(l=1.2, angle=0.02, k1=-1, k2=.9,tilt=0.2, eid="bend")
    s = Sextupole(l=0.23, k2=22, eid="sext", tilt=2.)
    c = Cavity(l=2, v=0.02, freq=1.3e9, phi=10, eid="cav")
    cor = Hcor(l=0.1, eid="cor")
    sol = Solenoid(l=0.12, k=0.002)
    d2 = Drift(l=0.5)
    tds = TDCavity(l=0.5, freq=1.2e6, phi=90, v=0.02, tilt=0.0, eid="tds")
    m = Monitor(eid="mon")
    mat = Matrix(l=0.5, r11=1.1, r22=1, r33=1, r44=1, r55=1, r66=1, r12=0.1, t111=0.8, tilt=0.5)
    b2 = RBend(l=1.2, angle=0.01, k1=-1, k2=-.9, eid="bend")
    b3 = SBend(l=1.2, angle=0.02, k1=1, k2=.5, eid="bend")

    init_energy = 0.5

    cell = (d, q, b, s, c, cor, sol,d2, tds, m, mat, b2, b3)

    lat = MagneticLattice(cell, method=MethodTM({'global': SecondTM}))

    R = lattice_transfer_map(lat, energy=init_energy)
    new_lat = merger(lat, remaining_types=[Drift], remaining_elems=[sol], init_energy=init_energy)
    R_new = lattice_transfer_map(new_lat, energy=init_energy)

    result = check_matrix(R, R_new, TOL, assert_info=' r_matrix - ')
    result2 = check_matrix(lat.T, new_lat.T, TOL, assert_info=' t_matrix - ')
    assert check_result(result + result2)


def test_merger_write_read(lattice, tws0, method, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""

    R = lattice_transfer_map(lattice, energy=tws0.E)

    new_lat = merger(lattice, remaining_types=[Hcor, Vcor, Monitor], remaining_elems=[MPBPMF_47_I1, START_96_I1], init_energy=tws0.E)

    write_lattice(new_lat, file_name="tmp_merger_lat.py")
    import tmp_merger_lat as ml
    new_lat2 = MagneticLattice(ml.cell, method=lattice.method)
    R_new = lattice_transfer_map(new_lat2, energy=tws0.E)


    result = check_matrix(R, R_new, tolerance=1.0e-8, tolerance_type='absolute', assert_info=' r_matrix - ')
    result2 = check_matrix(lattice.T, new_lat2.T, tolerance=1.0e-8, tolerance_type='absolute', assert_info=' t_matrix - ')
    assert check_result(result + result2)

def test_matrix_write_read(lattice, tws0, method, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""
    m = Matrix(l=0.3, delta_e=0.1)
    m.r = np.random.random((6, 6))
    m.t = np.random.random((6, 6, 6))
    m2 = Matrix(l=0.4)
    m2.r = np.random.random((6, 6))
    m2.t = np.random.random((6, 6, 6))

    lat = MagneticLattice((m, m2), method=method)
    R = lattice_transfer_map(lat, energy=tws0.E)

    write_lattice(lat, file_name="tmp_mat_lat.py")
    import tmp_mat_lat as mat
    lat2 = MagneticLattice(mat.cell, method=lat.method)
    R2 = lattice_transfer_map(lat2, energy=tws0.E)


    result = check_matrix(R, R2, TOL, assert_info='r_matrix - ')
    result2 = check_matrix(lat.T, lat2.T, TOL, assert_info='t_matrix - ')
    result3 = check_matrix(np.array([lat.E, lat.totalLen]), np.array([lat2.E,lat2.totalLen ]), TOL, assert_info='t_matrix - ')
    assert check_result(result + result2 + result3)


def test_matrix_b_vector(lattice, tws0, method, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""
    d = Drift(l=0.5)
    q = Quadrupole(l=0.3, k1=3, k2=3.3, eid="quad")
    q.dx = 0.0013
    q.dy = -0.0056
    b = Bend(l=1.2, angle=0.02, k1=-1, k2=.9, eid="bend")
    s = Sextupole(l=0.23, k2=22, eid="sext")
    c = Cavity(l=2, v=0.02, freq=1.3e9, phi=10, eid="cav")
    cor = Hcor(l=0.1, eid="cor")
    sol = Solenoid(l=0.12, k=0.002)
    d2 = Drift(l=0.5)
    tds = TDCavity(l=0.5, freq=1.2e6, phi=90, v=0.02, tilt=0.0, eid="tds")
    m = Monitor(eid="mon")
    mat = Matrix(l=0.5, r11=1.1, r22=1, r33=1, r44=1, r55=1, r66=1, r12=0.1, t111=0.8)
    b2 = RBend(l=1.2, angle=0.01, k1=-1, k2=-.9, eid="bend")
    b3 = SBend(l=1.2, angle=0.02, k1=1, k2=.5, eid="bend")
    b3.dx = -0.0004
    b3.dy = 0.0015

    init_energy = 0.5

    cell = (d, q, b, s, c, cor, sol,d2, tds, m, mat, b2, b3)

    lat = MagneticLattice(cell, method=MethodTM({'global': SecondTM}))

    R = lattice_transfer_map(lat, energy=init_energy)
    new_lat = merger(lat, remaining_types=[Drift], remaining_elems=[sol], init_energy=init_energy)
    R_new = lattice_transfer_map(new_lat, energy=init_energy)
    result = check_matrix(R, R_new, tolerance=1.0e-12, tolerance_type='absolute', assert_info=' r_matrix - ')
    result2 = check_matrix(lat.T, new_lat.T, tolerance=1.0e-12, tolerance_type='absolute', assert_info=' t_matrix - ')
    result3 = check_matrix(lat.B, new_lat.B, tolerance=1.0e-12, tolerance_type='absolute', assert_info=' b_vector - ')
    assert check_result(result + result2 + result3)


def test_matrix_b_vector_read_write(lattice, tws0, method, parametr=None, update_ref_values=False):
    """R maxtrix calculation test"""
    d = Drift(l=0.5)
    q = Quadrupole(l=0.3, k1=3, k2=3.3, eid="quad")
    q.dx = 0.0013
    q.dy = -0.0056
    b = Bend(l=1.2, angle=0.02, k1=-1, k2=.9, eid="bend")
    s = Sextupole(l=0.23, k2=22, eid="sext")
    c = Cavity(l=2, v=0.02, freq=1.3e9, phi=10, eid="cav")
    cor = Hcor(l=0.1, eid="cor")
    sol = Solenoid(l=0.12, k=0.002)
    d2 = Drift(l=0.5)
    tds = TDCavity(l=0.5, freq=1.2e6, phi=90, v=0.02, tilt=0.0, eid="tds")
    m = Monitor(eid="mon")
    mat = Matrix(l=0.5, r11=1.1, r22=1, r33=1, r44=1, r55=1, r66=1, r12=0.1, t111=0.8)
    b2 = RBend(l=1.2, angle=0.01, k1=-1, k2=-.9, eid="bend")
    b3 = SBend(l=1.2, angle=0.02, k1=1, k2=.5, eid="bend")
    b3.dx = -0.0004
    b3.dy = 0.0015

    init_energy = 0.5

    cell = (d, q, b, s, c, cor, sol,d2, tds, m, mat, b2, b3)

    lat = MagneticLattice(cell, method=method)

    R = lattice_transfer_map(lat, energy=init_energy)

    new_lat = merger(lat, remaining_types=[Drift], remaining_elems=[sol], init_energy=init_energy)
    write_lattice(new_lat, file_name="tmp_b_vec.py")
    import tmp_b_vec as b_vec
    lat_read = MagneticLattice(b_vec.cell, method=method)

    R_read = lattice_transfer_map(lat_read, energy=init_energy)
    result = check_matrix(R, R_read, tolerance=1.0e-8, tolerance_type='absolute', assert_info=' r_matrix - ')
    result2 = check_matrix(lat.T, lat_read.T, tolerance=1.0e-8, tolerance_type='absolute', assert_info=' t_matrix - ')
    result3 = check_matrix(lat.B, lat_read.B, tolerance=1.0e-8, tolerance_type='absolute', assert_info=' b_vector - ')
    assert check_result(result + result2 + result3)

def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### IO Lattice START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### IO Lattice END ###\n\n\n')
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
def test_update_ref_values(lattice, tws0, method, cmdopt):
    
    update_functions = []
    update_functions.append('test_original_lattice_transfer_map')
    update_functions.append('test_original_twiss')
    
    # function test_lat2input function need not be added here.
    # It is used reference results from test_original_lattice_transfer_map and test_original_twiss functions

    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, tws0, method, None, True)
        if result is None:
            return
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
