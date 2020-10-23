"""Test of the demo file demos/ebeam/csr_ex.py"""

import os
import sys
import copy
import time

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

from unit_tests.params import *
from section_track_conf import *


@pytest.mark.parametrize('parameter', [0, 1, 2, 3, 4])
def test_lattice_transfer_map(section_lat, p_array, parameter, update_ref_values=False):
    """R matrix calculation test"""
    sec_lat = copy.deepcopy(section_lat)
    sections = [A1, AH1, LH, DL, BC0, L1, BC1]
    sec = sections[parameter]
    s = sec_lat.dict_sections[sec]
    lat = s.lattice
    if parameter == 0:
        energy = 0.005
    elif parameter == 1:
        energy = 0.15
    else:
        energy = 0.13
    r_matrix = lattice_transfer_map(lat, energy)
    
    if update_ref_values:
        return numpy2json(r_matrix)

    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, tolerance=1.0e-12, tolerance_type='absolute', assert_info=' r_matrix - ')
    assert check_result(result)


@pytest.mark.parametrize('parameter', [0, 1, 2, 3, 4, 5, 6])
def test_lattice_transfer_map_update(section_lat, p_array, parameter, update_ref_values=False):
    SC_exec = False
    wake_exec = False
    CSR_exec = False
    match_exec = False
    smooth_exec = False
    v1 = 0.14747
    phi1 = -11.10
    v13 = 0.03079
    phi13 = -227.1
    V21 = 0.65919
    phi21 = 30.17
    r1 = 3.6587343247857467
    r2 = 9.392750737348779

    config = {
        A1: {"phi": phi1, "v": v1 / 8.,
             "SC": SC_exec, "smooth": True, "wake": wake_exec},
        AH1: {"phi": phi13, "v": v13 / 8,
              "match": True, "bounds": [-5, 5], "SC": SC_exec, "wake": wake_exec},
        LH: {"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec, "match": False},
        DL: {"match": False, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        BC0: {"rho": r1,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        L1: {"phi": phi21, "v": V21 / 32, "match": False,
             "SC": SC_exec, "wake": wake_exec, "smooth": smooth_exec},
        BC1: {"rho": r2,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
    }
    sec_lat = copy.deepcopy(section_lat)

    sections = [A1, AH1, LH, DL, BC0, L1, BC1]
    sec_lat.update_sections(sections, config=config, coupler_kick=False)

    sec = sections[parameter]
    s = sec_lat.dict_sections[sec]
    lat = s.lattice
    if parameter == 0:
        energy = 0.005
    elif parameter == 1:
        energy = 0.15
    else:
        energy = 0.13
    r_matrix = lattice_transfer_map(lat, energy)

    if update_ref_values:
        return numpy2json(r_matrix)

    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + str(parameter) + '.json'))

    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_tracking(section_lat, p_array, parameter=None, update_ref_values=False):
    sec_lat = copy.deepcopy(section_lat)

    parray = copy.deepcopy(p_array)
    SC_exec = False
    wake_exec = False
    CSR_exec = False
    match_exec = False
    smooth_exec = False
    v1 = 0.14747
    phi1 = -11.10
    v13 = 0.03079
    phi13 = -227.1
    V21 = 0.65919
    phi21 = 30.17
    r1 = 3.6587343247857467
    r2 = 9.392750737348779
    config = {
        A1: {"phi": phi1, "v": v1 / 8.,
             "SC": SC_exec, "smooth": True, "wake": wake_exec},
        AH1: {"phi": phi13, "v": v13 / 8,
              "match": True, "bounds": [-5, 5], "SC": SC_exec, "wake": wake_exec},
        LH: {"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec, "match": False},
        DL: {"match": False, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        BC0: {"rho": r1,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        L1: {"phi": phi21, "v": V21 / 32, "match": False,
             "SC": SC_exec, "wake": wake_exec, "smooth": smooth_exec},
        BC1: {"rho": r2,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
    }

    sections = [A1, AH1, LH, DL, BC0, L1, BC1]
    sec_lat.update_sections(sections, config=config, coupler_kick=False)

    parray = sec_lat.track_sections(sections=sections, p_array=parray, config=config, force_ext_p_array=False,
                                         coupler_kick=False)
    p = obj2dict(parray)

    if update_ref_values:
        return {'p_array': p}

    p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result2 = check_dict(p, p_array_ref['p_array'], TOL, assert_info=' p - ')
    assert check_result(result2)


def test_tracking_sc(section_lat, p_array, parameter=None, update_ref_values=False):
    sec_lat = copy.deepcopy(section_lat)

    parray = copy.deepcopy(p_array)

    SC_exec = True
    wake_exec = False
    CSR_exec = False
    match_exec = False
    smooth_exec = False
    v1 = 0.14747
    phi1 = -11.10
    v13 = 0.03079
    phi13 = -227.1
    V21 = 0.65919
    phi21 = 30.17
    r1 = 3.6587343247857467
    r2 = 9.392750737348779
    config = {
        A1: {"phi": phi1, "v": v1 / 8.,
             "SC": SC_exec, "smooth": True, "wake": wake_exec},
        AH1: {"phi": phi13, "v": v13 / 8,
              "match": True, "bounds": [-5, 5], "SC": SC_exec, "wake": wake_exec},
        LH: {"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec, "match": False},
        DL: {"match": False, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        BC0: {"rho": r1,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        L1: {"phi": phi21, "v": V21 / 32, "match": False,
             "SC": SC_exec, "wake": wake_exec, "smooth": smooth_exec},
        BC1: {"rho": r2,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
    }

    sections = [A1, AH1, LH, DL, BC0, L1, BC1]
    sec_lat.update_sections(sections, config=config, coupler_kick=False)

    parray = sec_lat.track_sections(sections=sections, p_array=parray, config=config, force_ext_p_array=False,
                                         coupler_kick=False)
    p = obj2dict(parray)

    if update_ref_values:
        return {'p_array': p}

    p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result2 = check_dict(p, p_array_ref['p_array'], tolerance=1.0e-10, tolerance_type='absolute', assert_info=' p - ')
    assert check_result(result2)


def test_tracking_match(section_lat, p_array, parameter=None, update_ref_values=False):
    sec_lat = copy.deepcopy(section_lat)

    parray = copy.deepcopy(p_array)

    SC_exec = False
    wake_exec = False
    CSR_exec = False
    match_exec = True
    smooth_exec = False
    v1 = 0.14747
    phi1 = -11.10
    v13 = 0.03079
    phi13 = -227.1
    V21 = 0.65919
    phi21 = 30.17
    r1 = 3.6587343247857467
    r2 = 9.392750737348779
    config = {
        A1: {"phi": phi1, "v": v1 / 8.,
             "SC": SC_exec, "smooth": True, "wake": wake_exec},
        AH1: {"phi": phi13, "v": v13 / 8,
              "match": True, "bounds": [-5, 5], "SC": SC_exec, "wake": wake_exec},
        LH: {"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec, "match": False},
        DL: {"match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        BC0: {"rho": r1,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        L1: {"phi": phi21, "v": V21 / 32, "match": False,
             "SC": SC_exec, "wake": wake_exec, "smooth": smooth_exec},
        BC1: {"rho": r2,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
    }

    sections = [A1, AH1, LH, DL, BC0, L1, BC1]
    sec_lat.update_sections(sections, config=config, coupler_kick=False)

    parray = sec_lat.track_sections(sections=sections, p_array=parray, config=config, force_ext_p_array=False,
                                         coupler_kick=False)
    p = obj2dict(parray)

    if update_ref_values:
        return {'p_array': p}

    p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result2 = check_dict(p, p_array_ref['p_array'], tolerance=1.0e-10, tolerance_type='absolute', assert_info=' p - ')
    assert check_result(result2)


def test_tracking_csr(section_lat, p_array, parameter=None, update_ref_values=False):
    sec_lat = copy.deepcopy(section_lat)

    parray = copy.deepcopy(p_array)
    SC_exec = False
    wake_exec = False
    CSR_exec = True
    match_exec = False
    smooth_exec = False
    v1 = 0.14747
    phi1 = -11.10
    v13 = 0.03079
    phi13 = -227.1
    V21 = 0.65919
    phi21 = 30.17
    r1 = 3.6587343247857467
    r2 = 9.392750737348779
    config = {
        A1: {"phi": phi1, "v": v1 / 8.,
             "SC": SC_exec, "smooth": smooth_exec, "wake": wake_exec},
        AH1: {"phi": phi13, "v": v13 / 8,
              "match": True, "bounds": [-5, 5], "SC": SC_exec, "wake": wake_exec},
        LH: {"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec, "match": False},
        DL: {"match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        BC0: {"rho": r1,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        L1: {"phi": phi21, "v": V21 / 32, "match": False,
             "SC": SC_exec, "wake": wake_exec, "smooth": smooth_exec},
        BC1: {"rho": r2,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
    }


    sections = [A1, AH1, LH, DL, BC0] #, L1, BC1]
    sec_lat.update_sections(sections, config=config, coupler_kick=False)

    parray = sec_lat.track_sections(sections=sections, p_array=parray, config=config, force_ext_p_array=False,
                                         coupler_kick=False)
    p = obj2dict(parray)

    if update_ref_values:
        return {'p_array': p}

    p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result2 = check_dict(p, p_array_ref['p_array'], tolerance=1.0e-10, tolerance_type='absolute', assert_info=' p - ')
    assert check_result(result2)


def test_tracking_ck(section_lat, p_array, parameter=None, update_ref_values=False):
    sec_lat = copy.deepcopy(section_lat)
    parray = copy.deepcopy(p_array)
    SC_exec = False
    wake_exec = False
    CSR_exec = False
    match_exec = False
    smooth_exec = False
    v1 = 0.14747
    phi1 = -11.10
    v13 = 0.03079
    phi13 = -227.1
    V21 = 0.65919
    phi21 = 30.17
    r1 = 3.6587343247857467
    r2 = 9.392750737348779
    config = {
        A1: {"phi": phi1, "v": v1 / 8.,
             "SC": SC_exec, "smooth": True, "wake": wake_exec},
        AH1: {"phi": phi13, "v": v13 / 8,
              "match": True, "bounds": [-5, 5], "SC": SC_exec, "wake": wake_exec},
        LH: {"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec, "match": False},
        DL: {"match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        BC0: {"rho": r1,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        L1: {"phi": phi21, "v": V21 / 32, "match": False,
             "SC": SC_exec, "wake": wake_exec, "smooth": smooth_exec},
        BC1: {"rho": r2,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
    }

    sections = [A1, AH1, LH, DL, BC0, L1, BC1]
    sec_lat.update_sections(sections, config=config, coupler_kick=True)

    parray = sec_lat.track_sections(sections=sections, p_array=parray, config=config, force_ext_p_array=False,
                                         coupler_kick=True)
    p = obj2dict(parray)

    if update_ref_values:
        return {'p_array': p}

    p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result2 = check_dict(p, p_array_ref['p_array'], tolerance=1.0e-10, tolerance_type='absolute', assert_info=' p - ')
    assert check_result(result2)


def test_tracking_wake(section_lat, p_array, parameter=None, update_ref_values=False):
    sec_lat = copy.deepcopy(section_lat)

    parray = copy.deepcopy(p_array)
    SC_exec = False
    wake_exec = True
    CSR_exec = False
    match_exec = False
    smooth_exec = False
    v1 = 0.14747
    phi1 = -11.10
    v13 = 0.03079
    phi13 = -227.1
    V21 = 0.65919
    phi21 = 30.17
    r1 = 3.6587343247857467
    r2 = 9.392750737348779
    config = {
        A1: {"phi": phi1, "v": v1 / 8.,
             "SC": SC_exec, "smooth": smooth_exec, "wake": wake_exec},
        AH1: {"phi": phi13, "v": v13 / 8,
              "match": True, "bounds": [-5, 5], "SC": SC_exec, "wake": wake_exec},
        LH: {"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec, "match": False},
        DL: {"match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        BC0: {"rho": r1,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        L1: {"phi": phi21, "v": V21 / 32, "match": False,
             "SC": SC_exec, "wake": wake_exec, "smooth": smooth_exec},
        BC1: {"rho": r2,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
    }


    sections = [A1, AH1, LH, DL, BC0, L1, BC1]
    sec_lat.update_sections(sections, config=config, coupler_kick=False)

    parray = sec_lat.track_sections(sections=sections, p_array=parray, config=config, force_ext_p_array=False,
                                         coupler_kick=False)
    p = obj2dict(parray)

    if update_ref_values:
        return {'p_array': p}

    p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result2 = check_dict(p, p_array_ref['p_array'], tolerance=1.0e-10, tolerance_type='absolute', assert_info=' p - ')
    assert check_result(result2)


def test_tracking_smooth(section_lat, p_array, parameter=None, update_ref_values=False):
    sec_lat = copy.deepcopy(section_lat)

    parray = copy.deepcopy(p_array)
    SC_exec = False
    wake_exec = False
    CSR_exec = False
    match_exec = False
    smooth_exec = True
    v1 = 0.14747
    phi1 = -11.10
    v13 = 0.03079
    phi13 = -227.1
    V21 = 0.65919
    phi21 = 30.17
    r1 = 3.6587343247857467
    r2 = 9.392750737348779
    config = {
        A1: {"phi": phi1, "v": v1 / 8.,
             "SC": SC_exec, "smooth": smooth_exec, "wake": wake_exec},
        AH1: {"phi": phi13, "v": v13 / 8,
              "match": True, "bounds": [-5, 5], "SC": SC_exec, "wake": wake_exec},
        LH: {"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec, "match": False},
        DL: {"match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        BC0: {"rho": r1,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        L1: {"phi": phi21, "v": V21 / 32, "match": False,
             "SC": SC_exec, "wake": wake_exec, "smooth": smooth_exec},
        BC1: {"rho": r2,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
    }


    sections = [A1, AH1, LH, DL, BC0, L1, BC1]
    sec_lat.update_sections(sections, config=config, coupler_kick=False)

    parray = sec_lat.track_sections(sections=sections, p_array=parray, config=config, force_ext_p_array=False,
                                         coupler_kick=False)
    p = obj2dict(parray)

    if update_ref_values:
        return {'p_array': p}

    p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result2 = check_dict(p, p_array_ref['p_array'], tolerance=1.0e-10, tolerance_type='absolute', assert_info=' p - ')
    assert check_result(result2)


def test_twiss(section_lat, p_array, parameter=None, update_ref_values=False):
    sec_lat = copy.deepcopy(section_lat)

    parray = copy.deepcopy(p_array)
    SC_exec = False
    wake_exec = False
    CSR_exec = False
    match_exec = False
    smooth_exec = False
    v1 = 0.14747
    phi1 = -11.10
    v13 = 0.03079
    phi13 = -227.1
    V21 = 0.65919
    phi21 = 30.17
    r1 = 3.6587343247857467
    r2 = 9.392750737348779
    config = {
        A1: {"phi": phi1, "v": v1 / 8.,
             "SC": SC_exec, "smooth": smooth_exec, "wake": wake_exec},
        AH1: {"phi": phi13, "v": v13 / 8,
              "match": True, "bounds": [-5, 5], "SC": SC_exec, "wake": wake_exec},
        LH: {"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec, "match": False},
        DL: {"match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        BC0: {"rho": r1,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        L1: {"phi": phi21, "v": V21 / 32, "match": False,
             "SC": SC_exec, "wake": wake_exec, "smooth": smooth_exec},
        BC1: {"rho": r2,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
    }


    sections = [A1, AH1, LH, DL, BC0, L1, BC1]
    sec_lat.update_sections(sections, config=config, coupler_kick=False)

    parray = sec_lat.track_sections(sections=sections, p_array=parray, config=config, force_ext_p_array=False,
                                         coupler_kick=False)

    tws_track_global = []
    L = 0
    for s in sections:
        sec = sec_lat.dict_sections[s]
        for tws in sec.tws_track:
            tws.s += L
        tws_track_global = np.append(tws_track_global, sec.tws_track)
        L += sec.lattice.totalLen

    tws = obj2dict(tws_track_global)

    if update_ref_values:
        return {'tws': tws}

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result2 = check_dict(tws, tws_ref['tws'], tolerance=1.0e-10, tolerance_type='absolute',  assert_info=' tws - ')
    assert check_result(result2)


def test_twiss_section_lat(section_lat, p_array, parameter=None, update_ref_values=False):
    sec_lat = copy.deepcopy(section_lat)

    parray = copy.deepcopy(p_array)
    SC_exec = False
    wake_exec = False
    CSR_exec = False
    match_exec = False
    smooth_exec = False
    v1 = 0.14747
    phi1 = -11.10
    v13 = 0.03079
    phi13 = -227.1
    V21 = 0.65919
    phi21 = 30.17
    r1 = 3.6587343247857467
    r2 = 9.392750737348779
    config = {
        A1: {"phi": phi1, "v": v1 / 8.,
             "SC": SC_exec, "smooth": smooth_exec, "wake": wake_exec},
        AH1: {"phi": phi13, "v": v13 / 8,
              "match": True, "bounds": [-5, 5], "SC": SC_exec, "wake": wake_exec},
        LH: {"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec, "match": False},
        DL: {"match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        BC0: {"rho": r1,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        L1: {"phi": phi21, "v": V21 / 32, "match": False,
             "SC": SC_exec, "wake": wake_exec, "smooth": smooth_exec},
        BC1: {"rho": r2,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
    }


    sections = [A1, AH1, LH, DL, BC0, L1, BC1]
    sec_lat.update_sections(sections, config=config, coupler_kick=False)

    parray = sec_lat.track_sections(sections=sections, p_array=parray, config=config, force_ext_p_array=False,
                                         coupler_kick=False)



    tws = obj2dict(sec_lat.tws_track)

    if update_ref_values:
        return {'tws': tws}

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result2 = check_dict(tws, tws_ref['tws'], tolerance=1.0e-10, tolerance_type='absolute', assert_info=' tws - ')
    assert check_result(result2)


def test_tracking_tds(section_lat, p_array, parameter=None, update_ref_values=False):
    sec_lat = copy.deepcopy(section_lat)

    parray = copy.deepcopy(p_array)
    SC_exec = False
    wake_exec = False
    CSR_exec = False
    match_exec = False
    smooth_exec = False
    v1 = 0.14747
    phi1 = -11.10
    v13 = 0.03079
    phi13 = -227.1
    V21 = 0.65919
    phi21 = 30.17
    r1 = 3.6587343247857467
    r2 = 9.392750737348779
    config = {
        A1: {"phi": phi1, "v": v1 / 8.,
             "SC": SC_exec, "smooth": True, "wake": wake_exec},
        AH1: {"phi": phi13, "v": v13 / 8,
              "match": True, "bounds": [-5, 5], "SC": SC_exec, "wake": wake_exec},
        LH: {"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec, "match": False, "tds.v": 0.001, "tds.phi":0},
        DL: {"match": False, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        BC0: {"rho": r1,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        L1: {"phi": phi21, "v": V21 / 32, "match": False,
             "SC": SC_exec, "wake": wake_exec, "smooth": smooth_exec},
        BC1: {"rho": r2,
              "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
    }

    sections = [A1, AH1, LH]
    sec_lat.update_sections(sections, config=config, coupler_kick=False)

    parray = sec_lat.track_sections(sections=sections, p_array=parray, config=config, force_ext_p_array=False,
                                         coupler_kick=False)
    p = obj2dict(parray)

    if update_ref_values:
        return {'p_array': p}

    p_array_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')

    result2 = check_dict(p, p_array_ref['p_array'], TOL, assert_info=' p - ')
    assert check_result(result2)

def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### CSR_EX START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### CSR_EX END ###\n\n\n')
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
def test_update_ref_values(section_lat, p_array, cmdopt):
    
    update_functions = []
    update_functions.append('test_lattice_transfer_map')
    update_functions.append('test_lattice_transfer_map_update')
    update_functions.append('test_tracking')
    update_functions.append('test_tracking_sc')
    update_functions.append('test_tracking_match')
    update_functions.append('test_tracking_csr')
    update_functions.append('test_tracking_ck')
    update_functions.append('test_tracking_wake')
    update_functions.append('test_tracking_smooth')
    update_functions.append('test_twiss')
    update_functions.append('test_twiss_section_lat')
    update_functions.append('test_tracking_tds')


    update_function_parameters = {}
    update_function_parameters['test_lattice_transfer_map'] = [0, 1, 2, 3, 4]
    update_function_parameters['test_lattice_transfer_map_update'] = [0, 1, 2, 3, 4, 5, 6]

    parameter = update_function_parameters[cmdopt] if cmdopt in update_function_parameters.keys() else ['']

    if cmdopt in update_functions:
        for p in parameter:
            p_arr = copy.deepcopy(p_array)
            sec_lat = copy.deepcopy(section_lat)
            result = eval(cmdopt)(sec_lat, p_arr, p, True)
        
            if os.path.isfile(REF_RES_DIR + cmdopt + str(p) + '.json'):
                os.rename(REF_RES_DIR + cmdopt + str(p) + '.json', REF_RES_DIR + cmdopt + str(p) + '.old')
            
            json_save(result, REF_RES_DIR + cmdopt + str(p) + '.json')


