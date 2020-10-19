"""Test of the demo file demos/ebeam/csr_ex.py"""

import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from navi_conf import *


def test_generate_parray(lattice, p_array, parameter=None, update_ref_values=False):
    """ func generate_parray testing """

    result = check_matrix(p_array.rparticles, np.zeros((6, 10)), tolerance=TOL, assert_info=' rparticles - ')
    assert check_result(result)


def test_navi_wo_procs(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test kick physProc with the same thick element
    """

    p_array_track = copy.deepcopy(p_array)
    lat = MagneticLattice(lattice.sequence)

    navi = Navigator(lat)

    tws_track_wo, p_array_wo = track(lat, p_array_track, navi, calc_tws=False)
    tws_track = obj2dict(tws_track_wo)
    p = obj2dict(p_array_wo)

    if update_ref_values:
        return {'tws_track': tws_track, 'p_array': p}

    tws_track_p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result1 = check_dict(tws_track, tws_track_p_array_ref['tws_track'], TOL, assert_info=' tws_track - ')
    result2 = check_dict(p, tws_track_p_array_ref['p_array'], TOL, assert_info=' p - ')

    assert check_result(result1 + result2)


def test_navi_wo_procs_reset_pos(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test kick physProc with the same thick element
    """

    p_array_track = copy.deepcopy(p_array)
    lat = MagneticLattice(lattice.sequence)

    navi = Navigator(lat)

    tws_1, p_array_1 = track(lat, p_array_track, navi, calc_tws=False)

    tws_o_1 = obj2dict(tws_1)
    p_o_1 = obj2dict(p_array_1)

    navi.reset_position()
    p_array_track_2 = copy.deepcopy(p_array)

    tws_2, p_array_2 = track(lat, p_array_track_2, navi, calc_tws=False)

    tws_o_2 = obj2dict(tws_2)
    p_o_2 = obj2dict(p_array_2)

    result1 = check_dict(tws_o_1, tws_o_2, TOL, assert_info=' tws_track - ')
    result2 = check_dict(p_o_1, p_o_2, TOL, assert_info=' p - ')

    assert check_result(result1 + result2)


def test_kick_marker(lattice, p_array, parameter=None, update_ref_values=False):
    """
    testing applying one marker as start ans stop
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.05

    t = LogProc()

    navi.add_physics_proc(t, m_kick, m_kick)

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi, calc_tws=False)



    result0 = check_value(t.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen - ')
    result1 = check_matrix(np.array(t.dz_list), np.array([0]), tolerance=TOL, assert_info=' dz - ')
    result2 = check_matrix(np.array(t.L_list), np.array([0]), tolerance=TOL, assert_info=' L - ')
    result3 = check_matrix(np.array(t.s_list), np.array([0.3]), tolerance=TOL, assert_info=' p.s - ')
    result4 = check_matrix(np.array(t.s_stop_list), np.array([0.3]), tolerance=TOL, assert_info=' s_stop - ')
    result5 = check_matrix(np.array(t.s_start_list), np.array([0.3]), tolerance=TOL, assert_info=' s_start - ')


    assert check_result([result0] + result1 + result2 + result3 + result4 + result5)


def test_two_markers(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test physProc between two markers - standard case
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.1

    t = LogProc()
    navi.add_physics_proc(t, m1, m2)

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi, calc_tws=False)


    result0 = check_value(t.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen - ')
    result1 = check_matrix(np.array(t.dz_list), np.ones(5)*0.1, tolerance=TOL, assert_info=' dz - ')
    result2 = check_matrix(np.array(t.L_list), np.ones(5)*0.5, tolerance=TOL, assert_info=' L - ')
    result3 = check_matrix(np.array(t.s_list), np.arange(1.05, 1.55, 0.1), tolerance=TOL, assert_info=' p.s - ')
    result4 = check_matrix(np.array(t.s_stop_list), np.ones(5)*1.45, tolerance=TOL, assert_info=' s_stop - ')
    result5 = check_matrix(np.array(t.s_start_list), np.ones(5)*0.95, tolerance=TOL, assert_info=' s_start - ')
    assert check_result([result0] + result1 + result2 + result3 + result4 + result5)


def test_two_elements(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test physProc between two elements with non zero length
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.1

    t = LogProc()
    navi.add_physics_proc(t, B, D1)

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi, calc_tws=False)


    result0 = check_value(t.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen - ')
    result1 = check_matrix(np.array(t.dz_list), np.array([0.1, 0.1, 0.05]), tolerance=TOL, assert_info=' dz - ')
    result2 = check_matrix(np.array(t.L_list), np.ones(3)*0.25, tolerance=TOL, assert_info=' L - ')
    result3 = check_matrix(np.array(t.s_list), np.array([0.4, 0.5, 0.55]), tolerance=TOL, assert_info=' p.s - ')
    result4 = check_matrix(np.array(t.s_stop_list), np.ones(3)*0.55, tolerance=TOL, assert_info=' s_stop - ')
    result5 = check_matrix(np.array(t.s_start_list), np.ones(3)*0.3, tolerance=TOL, assert_info=' s_start - ')
    assert check_result([result0] + result1 + result2 + result3 + result4 + result5)

def test_two_markers2(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test physProc between two markers but non equal steps - standard case
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.2

    t = LogProc()
    navi.add_physics_proc(t, m1, m2)

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi, calc_tws=False)


    result0 = check_value(t.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen - ')
    result1 = check_matrix(np.array(t.dz_list), np.array([0.2, 0.2, 0.1]), tolerance=TOL, assert_info=' dz - ')
    result2 = check_matrix(np.array(t.L_list), np.ones(3)*0.5, tolerance=TOL, assert_info=' L - ')
    result3 = check_matrix(np.array(t.s_list),  np.array([1.15, 1.35, 1.45]), tolerance=TOL, assert_info=' p.s - ')
    result4 = check_matrix(np.array(t.s_stop_list), np.ones(3)*1.45, tolerance=TOL, assert_info=' s_stop - ')
    result5 = check_matrix(np.array(t.s_start_list), np.ones(3)*0.95, tolerance=TOL, assert_info=' s_start - ')
    assert check_result([result0] + result1 + result2 + result3 + result4 + result5)


def test_start_stop(lattice, p_array, parameter=None, update_ref_values=False):
    """
    testing PhysProc between two markers at the ends of lattice
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.2

    t = LogProc()
    navi.add_physics_proc(t, lattice.sequence[0], lattice.sequence[-1])

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi, calc_tws=False)


    result0 = check_value(t.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen - ')
    result1 = check_matrix(np.array(t.dz_list), np.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.05]), tolerance=TOL, assert_info=' dz - ')
    result2 = check_matrix(np.array(t.L_list), np.ones(8)*1.45, tolerance=TOL, assert_info=' L - ')
    result3 = check_matrix(np.array(t.s_list),  np.cumsum(np.array([0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.05])), tolerance=TOL, assert_info=' p.s - ')
    result4 = check_matrix(np.array(t.s_stop_list), np.ones(8)*1.45, tolerance=TOL, assert_info=' s_stop - ')
    result5 = check_matrix(np.array(t.s_start_list), np.ones(8)*0.0, tolerance=TOL, assert_info=' s_start - ')
    assert check_result([result0] + result1 + result2 + result3 + result4 + result5)


def test_ss_and_markers(lattice, p_array, parameter=None, update_ref_values=False):
    """
    Testing two PhysProc one inside another
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.2

    t = LogProc()
    navi.add_physics_proc(t, lattice.sequence[0], lattice.sequence[-1])

    t2 = LogProc()
    navi.add_physics_proc(t2, m1, m2)

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi, calc_tws=False)


    result0 = check_value(t.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen - ')
    result1 = check_matrix(np.array(t.dz_list), np.array([0.2, 0.2, 0.2, 0.2, 0.15, 0.2, 0.2, 0.1]), tolerance=TOL, assert_info=' dz - ')
    result2 = check_matrix(np.array(t.L_list), np.ones(8)*1.45, tolerance=TOL, assert_info=' L - ')
    result3 = check_matrix(np.array(t.s_list),  np.cumsum(np.array([0.2, 0.2, 0.2, 0.2, 0.15, 0.2, 0.2, 0.1])), tolerance=TOL, assert_info=' p.s - ')
    result4 = check_matrix(np.array(t.s_stop_list), np.ones(8)*1.45, tolerance=TOL, assert_info=' s_stop - ')
    result5 = check_matrix(np.array(t.s_start_list), np.ones(8)*0.0, tolerance=TOL, assert_info=' s_start - ')
    assert check_result([result0] + result1 + result2 + result3 + result4 + result5)

    result02 = check_value(t2.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen2 - ')
    result12 = check_matrix(np.array(t2.dz_list), np.array([0.2, 0.2, 0.1]), tolerance=TOL, assert_info=' dz2 - ')
    result22 = check_matrix(np.array(t2.L_list), np.ones(3)*0.5, tolerance=TOL, assert_info=' L2 - ')
    result32 = check_matrix(np.array(t2.s_list),  np.array([1.15, 1.35, 1.45]), tolerance=TOL, assert_info=' p.s2 - ')
    result42 = check_matrix(np.array(t2.s_stop_list), np.ones(3)*1.45, tolerance=TOL, assert_info=' s_stop2 - ')
    result52 = check_matrix(np.array(t2.s_start_list), np.ones(3)*0.95, tolerance=TOL, assert_info=' s_start2 - ')
    assert check_result([result02] + result12 + result22 + result32 + result42 + result52)


def test_ss_and_markers_steps(lattice, p_array, parameter=None, update_ref_values=False):
    """
    Testing two PhysProc one inside another but different step size (1 and 2)
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.1

    t = LogProc()
    t.step = 2

    navi.add_physics_proc(t, lattice.sequence[0], lattice.sequence[-1])


    t2 = LogProc()


    navi.add_physics_proc(t2, m1, m2)

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi, calc_tws=False)


    result0 = check_value(t.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen - ')
    result1 = check_matrix(np.array(t.dz_list), np.array([0.2, 0.2, 0.2, 0.2, 0.15, 0.2, 0.2, 0.1]), tolerance=TOL, assert_info=' dz - ')
    result2 = check_matrix(np.array(t.L_list), np.ones(8)*1.45, tolerance=TOL, assert_info=' L - ')
    result3 = check_matrix(np.array(t.s_list),  np.cumsum(np.array([0.2, 0.2, 0.2, 0.2, 0.15, 0.2, 0.2, 0.1])), tolerance=TOL, assert_info=' p.s - ')
    result4 = check_matrix(np.array(t.s_stop_list), np.ones(8)*1.45, tolerance=TOL, assert_info=' s_stop - ')
    result5 = check_matrix(np.array(t.s_start_list), np.ones(8)*0.0, tolerance=TOL, assert_info=' s_start - ')
    assert check_result([result0] + result1 + result2 + result3 + result4 + result5)

    result02 = check_value(t2.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen2 - ')
    result12 = check_matrix(np.array(t2.dz_list), np.ones(5)*0.1, tolerance=TOL, assert_info=' dz2 - ')
    result22 = check_matrix(np.array(t2.L_list), np.ones(5)*0.5, tolerance=TOL, assert_info=' L2 - ')
    result32 = check_matrix(np.array(t2.s_list), np.arange(1.05, 1.55, 0.1), tolerance=TOL, assert_info=' p.s2 - ')
    result42 = check_matrix(np.array(t2.s_stop_list), np.ones(5)*1.45, tolerance=TOL, assert_info=' s_stop2 - ')
    result52 = check_matrix(np.array(t2.s_start_list), np.ones(5)*0.95, tolerance=TOL, assert_info=' s_start2 - ')
    assert check_result([result02] + result12 + result22 + result32 + result42 + result52)


def test_ends_markers(lattice, p_array, parameter=None, update_ref_values=False):
    """
    Testing Kick PhysProcs at the ends of lattice with other PhhysProces
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.1

    t = LogProc()
    t.step = 2

    navi.add_physics_proc(t, lattice.sequence[0], lattice.sequence[-1])


    t2 = LogProc()

    navi.add_physics_proc(t2, m1, m2)

    t3 = LogProc()

    navi.add_physics_proc(t3, start, m_extra1)

    t4 = LogProc()

    navi.add_physics_proc(t4, m_extra2, stop)

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi, calc_tws=False)


    result0 = check_value(t.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen - ')
    result1 = check_matrix(np.array(t.dz_list), np.array([0, 0.2, 0.2, 0.2, 0.2, 0.15, 0.2, 0.2, 0.1]), tolerance=TOL, assert_info=' dz - ')
    result2 = check_matrix(np.array(t.L_list), np.ones(9)*1.45, tolerance=TOL, assert_info=' L - ')
    result3 = check_matrix(np.array(t.s_list),  np.cumsum(np.array([0, 0.2, 0.2, 0.2, 0.2, 0.15, 0.2, 0.2, 0.1])), tolerance=TOL, assert_info=' p.s - ')
    result4 = check_matrix(np.array(t.s_stop_list), np.ones(9)*1.45, tolerance=TOL, assert_info=' s_stop - ')
    result5 = check_matrix(np.array(t.s_start_list), np.ones(9)*0.0, tolerance=TOL, assert_info=' s_start - ')
    assert check_result([result0] + result1 + result2 + result3 + result4 + result5)

    result02 = check_value(t2.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen2 - ')
    result12 = check_matrix(np.array(t2.dz_list), np.ones(5)*0.1, tolerance=TOL, assert_info=' dz2 - ')
    result22 = check_matrix(np.array(t2.L_list), np.ones(5)*0.5, tolerance=TOL, assert_info=' L2 - ')
    result32 = check_matrix(np.array(t2.s_list), np.arange(1.05, 1.55, 0.1), tolerance=TOL, assert_info=' p.s2 - ')
    result42 = check_matrix(np.array(t2.s_stop_list), np.ones(5)*1.45, tolerance=TOL, assert_info=' s_stop2 - ')
    result52 = check_matrix(np.array(t2.s_start_list), np.ones(5)*0.95, tolerance=TOL, assert_info=' s_start2 - ')
    assert check_result([result02] + result12 + result22 + result32 + result42 + result52)

    result03 = check_value(t3.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen3 - ')
    result13 = check_matrix(np.array(t3.dz_list), np.ones(1)*0., tolerance=TOL, assert_info=' dz3 - ')
    result23 = check_matrix(np.array(t3.L_list), np.ones(1)*0., tolerance=TOL, assert_info=' L3 - ')
    result33 = check_matrix(np.array(t3.s_list), np.array(0.), tolerance=TOL, assert_info=' p.s3 - ')
    result43 = check_matrix(np.array(t3.s_stop_list), np.ones(1)*0., tolerance=TOL, assert_info=' s_stop3 - ')
    result53 = check_matrix(np.array(t3.s_start_list), np.ones(1)*0.0, tolerance=TOL, assert_info=' s_start3 - ')
    assert check_result([result03] + result13 + result23 + result33 + result43 + result53)

    result04 = check_value(t4.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen4 - ')
    result14 = check_matrix(np.array(t4.dz_list), np.ones(1)*0., tolerance=TOL, assert_info=' dz4 - ')
    result24 = check_matrix(np.array(t4.L_list), np.ones(1)*0., tolerance=TOL, assert_info=' L4 - ')
    result34 = check_matrix(np.array(t4.s_list), np.array(1.45), tolerance=TOL, assert_info=' p.s4 - ')
    result44 = check_matrix(np.array(t4.s_stop_list), np.ones(1)*1.45, tolerance=TOL, assert_info=' s_stop4 - ')
    result54 = check_matrix(np.array(t4.s_start_list), np.ones(1)*1.45, tolerance=TOL, assert_info=' s_start4 - ')
    assert check_result([result04] + result14 + result24 + result34 + result44 + result54)


def test_kick_start(lattice, p_array, parameter=None, update_ref_values=False):
    """
    Testing different kick PhysProces applied to the same markers
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.1

    t = LogProc()
    t.step = 2

    navi.add_physics_proc(t, start, start)

    t2 = LogProc()

    navi.add_physics_proc(t2, start, m_extra1)


    t3 = LogProc()
    t3.step = 1
    navi.add_physics_proc(t3, start, m_kick)

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi, calc_tws=False)

    dz_list = np.array([ 0.])
    result0 = check_value(t.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen - ')
    result1 = check_matrix(np.array(t.dz_list), dz_list, tolerance=TOL, assert_info=' dz - ')
    result2 = check_matrix(np.array(t.L_list), np.ones(len(dz_list))*0.0, tolerance=TOL, assert_info=' L - ')
    result3 = check_matrix(np.array(t.s_list),  np.cumsum(dz_list), tolerance=TOL, assert_info=' p.s - ')
    result4 = check_matrix(np.array(t.s_stop_list), np.ones(len(dz_list))*0.0, tolerance=TOL, assert_info=' s_stop - ')
    result5 = check_matrix(np.array(t.s_start_list), np.ones(len(dz_list))*0.0, tolerance=TOL, assert_info=' s_start - ')
    assert check_result([result0] + result1 + result2 + result3 + result4 + result5)

    dz_list = np.array([ 0.])
    result02 = check_value(t2.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen - ')
    result12 = check_matrix(np.array(t2.dz_list), dz_list, tolerance=TOL, assert_info=' dz - ')
    result22 = check_matrix(np.array(t2.L_list), np.ones(len(dz_list))*0.0, tolerance=TOL, assert_info=' L - ')
    result32 = check_matrix(np.array(t2.s_list),  np.cumsum(dz_list), tolerance=TOL, assert_info=' p.s - ')
    result42 = check_matrix(np.array(t2.s_stop_list), np.ones(len(dz_list))*0.0, tolerance=TOL, assert_info=' s_stop - ')
    result52 = check_matrix(np.array(t2.s_start_list), np.ones(len(dz_list))*0.0, tolerance=TOL, assert_info=' s_start - ')
    assert check_result([result02] + result12 + result22 + result32 + result42 + result52)

    dz_list = np.array([ 0, 0.1, 0.1, 0.1])
    result03 = check_value(t3.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen3 - ')
    result13 = check_matrix(np.array(t3.dz_list), dz_list, tolerance=TOL, assert_info=' dz3 - ')
    result23 = check_matrix(np.array(t3.L_list), np.ones(len(dz_list))*0.3, tolerance=TOL, assert_info=' L3 - ')
    result33 = check_matrix(np.array(t3.s_list), np.cumsum(dz_list)+0, tolerance=TOL, assert_info=' p.s3 - ')
    result43 = check_matrix(np.array(t3.s_stop_list), np.ones(len(dz_list))*0.3, tolerance=TOL, assert_info=' s_stop3 - ')
    result53 = check_matrix(np.array(t3.s_start_list), np.ones(len(dz_list))*0., tolerance=TOL, assert_info=' s_start3 - ')
    assert check_result([result03] + result13 + result23 + result33 + result43 + result53)

def test_two_dist_elems(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test physProc between two markers but non equal steps - standard case
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.2

    t = LogProc()
    navi.add_physics_proc(t, B, D2)

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi, calc_tws=False)

    dz_list = np.array([0.2, 0.2, 0.2, 0.05])
    result0 = check_value(t.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen - ')
    result1 = check_matrix(np.array(t.dz_list), dz_list, tolerance=TOL, assert_info=' dz - ')
    result2 = check_matrix(np.array(t.L_list), np.ones(len(dz_list))*0.65, tolerance=TOL, assert_info=' L - ')
    result3 = check_matrix(np.array(t.s_list),  np.cumsum(dz_list) + 0.3, tolerance=TOL, assert_info=' p.s - ')
    result4 = check_matrix(np.array(t.s_stop_list), np.ones(len(dz_list))*0.95, tolerance=TOL, assert_info=' s_stop - ')
    result5 = check_matrix(np.array(t.s_start_list), np.ones(len(dz_list))*0.3, tolerance=TOL, assert_info=' s_start - ')
    assert check_result([result0] + result1 + result2 + result3 + result4 + result5)

def test_kick_with_one_elem(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test kick physProc with the same thick element
    """

    p_array_track = copy.deepcopy(p_array)


    navi = Navigator(lattice)
    navi.unit_step = 0.2

    t = LogProc()
    navi.add_physics_proc(t, B, B)

    tws_track_wo, p_array_wo = track(lattice, p_array_track, navi, calc_tws=False)

    dz_list = np.array([0.0])
    result0 = check_value(t.totalLen, 1.45, tolerance=TOL, assert_info=' totalLen - ')
    result1 = check_matrix(np.array(t.dz_list), dz_list, tolerance=TOL, assert_info=' dz - ')
    result2 = check_matrix(np.array(t.L_list), np.ones(len(dz_list))*0., tolerance=TOL, assert_info=' L - ')
    result3 = check_matrix(np.array(t.s_list),  np.cumsum(dz_list) + 0.3, tolerance=TOL, assert_info=' p.s - ')
    result4 = check_matrix(np.array(t.s_stop_list), np.ones(len(dz_list))*0.3, tolerance=TOL, assert_info=' s_stop - ')
    result5 = check_matrix(np.array(t.s_start_list), np.ones(len(dz_list))*0.3, tolerance=TOL, assert_info=' s_start - ')
    assert check_result([result0] + result1 + result2 + result3 + result4 + result5)


def test_kick_with_thick_elem(lattice, p_array, parameter=None, update_ref_values=False):
    """
    test kick physProc with the same thick element
    """

    p_array_track = copy.deepcopy(p_array)
    lat = MagneticLattice(lattice.sequence, start=B)

    navi = Navigator(lat)
    navi.unit_step = 0.2

    t = LogProc()
    t2 = LogProc()
    navi.add_physics_proc(t, B, B)
    navi.add_physics_proc(t2, B, B)

    tws_track_wo, p_array_wo = track(lat, p_array_track, navi, calc_tws=False)

    dz_list = np.array([0.0])
    result0 = check_value(t.totalLen, 1.15, tolerance=TOL, assert_info=' totalLen - ')
    result1 = check_matrix(np.array(t.dz_list), dz_list, tolerance=TOL, assert_info=' dz t1 - ')
    result2 = check_matrix(np.array(t2.dz_list), dz_list, tolerance=TOL, assert_info=' dz t2 - ')
    #result2 = check_matrix(np.array(t.L_list), np.ones(len(dz_list))*0., tolerance=TOL, assert_info=' L - ')
    result3 = check_matrix(np.array(t.s_list),  np.cumsum(dz_list) , tolerance=TOL, assert_info=' p.s t1- ')
    result4 = check_matrix(np.array(t2.s_list),  np.cumsum(dz_list) , tolerance=TOL, assert_info=' p.s t2- ')

    #result4 = check_matrix(np.array(t.s_stop_list), np.ones(len(dz_list))*0.3, tolerance=TOL, assert_info=' s_stop - ')
    #result5 = check_matrix(np.array(t.s_start_list), np.ones(len(dz_list))*0.3, tolerance=TOL, assert_info=' s_start - ')
    assert check_result([result0] + result1 + result2 + result3 + result4)




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
    update_functions.append('test_navi_wo_procs')
    #update_functions.append('test_kick_marker')

    update_function_parameters = {}
    #update_function_parameters['test_track_smooth'] = [0, 1]

    parameter = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parameter:
            p_arr = copy.deepcopy(p_array)
            result = eval(cmdopt)(lattice, p_arr, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')
