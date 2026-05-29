"""Test of the demo file demos/ebeam/csr_ex.py"""

from pathlib import Path
import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from wake_conf import *

HERE = Path(__file__).parent


def test_long_wake_reads_supported_headers(tmp_path):
    """LongWake reads V/pC files with common % and # headers."""

    long_file = tmp_path / "long_cm.dat"
    long_file.write_text(
        "%bunch rms\n"
        "1.000000e-01\t0.000000e+00\n"
        "%monopole wake\n"
        "%s[cm]\tW(s)[V/pC]\n"
        "-1.0\t2.0\n"
        "0.0\t3.0\n"
    )
    transverse_file = tmp_path / "transverse_mm.dat"
    transverse_file.write_text(
        "#Parameters = {sigma=0.25; x=7.025}\n"
        '#"s / mm"\t"Y [Real Part / V/pC]"\n'
        "#-------------------------------\n"
        "-1.0\t1.0\n"
        "0.0\t2.0\n"
    )

    wake = LongWake(long_wake_file=long_file, transverse_wake_file=transverse_file, driver_charge=2.0)
    wake.prepare(None)

    assert np.allclose(wake.long_wake[:, 0], [-0.01, 0.0])
    assert np.allclose(wake.long_wake[:, 1], [4.0, 6.0])
    assert np.allclose(wake.transverse_wake[:, 0], [-0.001, 0.0])
    assert np.allclose(wake.transverse_wake[:, 1], [2.0, 4.0])


def test_long_wake_requires_at_least_one_file():
    wake = LongWake()

    with pytest.raises(ValueError, match="requires"):
        wake.prepare(None)


def test_long_wake_applies_longitudinal_file_only(tmp_path):
    long_file = tmp_path / "long_m.dat"
    long_file.write_text("-1e-3 1.0\n0.0 2.0\n1e-3 3.0\n")
    p_array = ParticleArray(3)
    p_array.E = 1.0
    p_array.rparticles[4] = [-1e-3, 0.0, 1e-3]

    wake = LongWake(long_wake_file=long_file, driver_charge=10.0, wake_file_unit="m")
    wake.prepare(None)
    wake.s_start = 0.0
    wake.s_stop = 2.0
    wake.apply(p_array, dz=1.0)

    expected = np.array([1.0, 2.0, 3.0]) * 10.0 * 0.5 * 1e-9 / p_array.p0c
    assert np.allclose(p_array.p(), expected)
    assert np.allclose(p_array.px(), 0.0)
    assert np.allclose(p_array.py(), 0.0)


def test_long_wake_applies_transverse_file_only(tmp_path):
    transverse_file = tmp_path / "transverse_m.dat"
    transverse_file.write_text("-1e-3 4.0\n1e-3 8.0\n")
    p_array = ParticleArray(2)
    p_array.E = 2.0
    p_array.rparticles[4] = [-1e-3, 1e-3]

    wake = LongWake(
        transverse_wake_file=transverse_file,
        driver_charge=5.0,
        factor=2.0,
        wake_file_unit="m",
        transverse_plane="x",
        transverse_sign=-1.0,
    )
    wake.prepare(None)
    wake.s_start = 0.0
    wake.s_stop = 1.0
    wake.apply(p_array, dz=1.0)

    expected = -np.array([4.0, 8.0]) * 5.0 * 2.0 * 1e-9 / p_array.p0c
    assert np.allclose(p_array.px(), expected)
    assert np.allclose(p_array.py(), 0.0)
    assert np.allclose(p_array.p(), 0.0)


def test_long_wake_sorts_loaded_wake_and_interpolates(tmp_path):
    long_file = tmp_path / "long_unsorted.dat"
    long_file.write_text("1e-3 3.0\n-1e-3 1.0\n0.0 2.0\n")

    wake = LongWake(long_wake_file=long_file, driver_charge=10.0, wake_file_unit="m")
    wake.prepare(None)

    assert np.allclose(wake.long_wake[:, 0], [-1e-3, 0.0, 1e-3])
    assert np.allclose(wake.sample_wake(wake.long_wake, np.array([-2e-3, -0.5e-3, 0.5e-3, 2e-3])),
                       [0.0, 15.0, 25.0, 0.0])


def test_long_wake_centers_witness_at_beam_position(tmp_path):
    long_file = tmp_path / "long_centered.dat"
    long_file.write_text("9e-3 9.0\n10e-3 10.0\n11e-3 11.0\n")
    p_array = ParticleArray(3)
    p_array.E = 1.0
    p_array.rparticles[4] = [1e-3, 2e-3, 3e-3]

    wake = LongWake(
        long_wake_file=long_file,
        driver_charge=2.0,
        beam_position=10e-3,
        wake_file_unit="m",
    )
    wake.prepare(None)
    wake.s_start = 0.0
    wake.s_stop = 1.0
    wake.apply(p_array, dz=1.0)

    expected = np.array([9.0, 10.0, 11.0]) * 2.0 * 1e-9 / p_array.p0c
    assert np.allclose(p_array.p(), expected)


def test_lin_long_wake_applies_linear_voltage():
    p_array = ParticleArray(3)
    p_array.E = 1.0
    p_array.rparticles[4] = [-1e-3, 0.0, 1e-3]

    wake = LinLongWake(voltage=10.0, derivative=1000.0, position=0.0)
    wake.prepare(None)
    wake.s_start = 0.0
    wake.s_stop = 2.0
    wake.apply(p_array, dz=1.0)

    voltage = np.array([9.0, 10.0, 11.0])
    expected = voltage * 0.5 * 1e-9 / p_array.p0c
    assert np.allclose(p_array.p(), expected)


def test_lin_long_wake_uses_position_for_constant_voltage_point():
    p_array = ParticleArray(2)
    p_array.E = 2.0
    p_array.rparticles[4] = [0.0, 2e-3]

    wake = LinLongWake(voltage=5.0, derivative=1000.0, position=1e-3)
    wake.prepare(None)
    wake.apply(p_array, dz=0.0)

    voltage = np.array([4.0, 6.0])
    expected = voltage * 1e-9 / p_array.p0c
    assert np.allclose(p_array.p(), expected)


def test_get_long_wake(lattice, p_array, parameter=None, update_ref_values=False):
    """ func generate_parray testing """

    wt = WakeTable(wake_file=HERE / "wake_table.dat")
    w = Wake()
    w.wake_table = wt
    w.prepare(None)

    x = np.arange(-30, 30, 0.01) * 1e-6
    sigma = 7e-6
    y2_x = lambda x: 5000 if np.abs(x) < 10e-6 else 0.
    y = np.array([y2_x(xi) for xi in x])


    profile = np.hstack((x.reshape(-1, 1), y.reshape(-1, 1)))
    x, Wz = w.get_long_wake(profile)

    if update_ref_values:
        return {'x': list(x), 'Wz': list(Wz)}

    current_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    result1 = check_matrix(x, current_ref['x'], TOL, assert_info=' x - ')
    result2 = check_matrix(Wz, current_ref['Wz'], TOL, assert_info=' Wz - ')
    assert check_result(result1 + result2)

def test_get_dipole_wake(lattice, p_array, parameter=None, update_ref_values=False):
    """ func generate_parray testing """

    wt = WakeTable(wake_file=HERE / "wake_table.dat")
    w = Wake()
    w.wake_table = wt
    w.prepare(None)

    x = np.arange(-30, 30, 0.01) * 1e-6
    sigma = 7e-6
    y2_x = lambda x: 5000 if np.abs(x) < 10e-6 else 0.
    y = np.array([y2_x(xi) for xi in x])


    profile = np.hstack((x.reshape(-1, 1), y.reshape(-1, 1)))
    x, Wd = w.get_dipole_wake(profile)

    if update_ref_values:
        return {'x': list(x), 'Wd': list(Wd)}

    current_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    result1 = check_matrix(x, current_ref['x'], TOL, assert_info=' x - ')
    result2 = check_matrix(Wd, current_ref['Wd'], TOL, assert_info=' Wd - ')
    assert check_result(result1 + result2)


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
def test_update_ref_values(lattice, p_array, cmdopt):
    
    update_functions = []
    update_functions.append('test_get_long_wake')
    update_functions.append("test_get_dipole_wake")

    update_function_parameters = {}


    parametr = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parametr:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p_arr, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
