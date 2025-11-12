"""Test of the demo file demos/ebeam/dba.py"""

import os
import sys
import time

import numpy as np

from ocelot.utils import bba
from unit_tests.params import *
from bba_conf import *

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
REF_RES_DIR = FILE_DIR + '/ref_results/'

def _ref_path(name: str) -> str:
    return os.path.join(REF_RES_DIR, f"{name}.json")


def test_lattice_transfer_map(lattice, update_ref_values=False):
    """R maxtrix calculation test"""

    r_matrix = lattice_transfer_map(lattice, 14)
    
    if update_ref_values:
        return numpy2json(r_matrix)
    
    r_matrix_ref = json2numpy(json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json'))
    
    result = check_matrix(r_matrix, r_matrix_ref, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_twiss(lattice, update_ref_values=False):
    """Twiss parameters calculation function test"""

    tws = twiss(lattice, sase1.tws0)
    
    tws = obj2dict(tws)
    
    if update_ref_values:
        return tws

    tws_ref = json_read(REF_RES_DIR + sys._getframe().f_code.co_name + '.json')
    
    result = check_dict(tws, tws_ref, TOL, 'absotute', assert_info=' tws - ')
    assert check_result(result)


def test_read_orbit(lattice, update_ref_values=False):
    # --- pick elements
    quads, bpms = [], []
    for e in lattice.sequence:
        eid = getattr(e, "id", "")
        if isinstance(e, Quadrupole) and ".SA1" in eid:
            quads.append(e)
        elif isinstance(e, Monitor) and ".SA1" in eid:
            bpms.append(e)
    bx = np.zeros(len(bpms))
    by = np.zeros(len(bpms))
    bx[3] = -0.001
    bx[6] = 0.001
    by[8] = -0.002
    by[20] = 0.003
    mx, my, s = bba.read_orbit(
                lattice,
                bpms,
                pinit= Particle(x=1e-3, y=-1e-3, px=1e-5, py=-1e-4),
                bpm_offset_x = bx,
                bpm_offset_y = by,
                noise_rms = (0.0, 0.0),
                noise_truncated= 3.0)


    # --- serialize for update mode
    payload = {
        "n_bpms": len(bpms),
        "n_quads": len(quads),
        "mx": numpy2json(mx),
        "my": numpy2json(my),
        "s":  numpy2json(s),
    }

    ref_path = _ref_path(sys._getframe().f_code.co_name)
    if update_ref_values:
        return payload

    # --- compare
    ref = json_read(ref_path)
    results = []

    # Optional meta checks
    if payload["n_bpms"] != ref["n_bpms"]:
        results.append(f"n_bpms mismatch: {payload['n_bpms']} vs {ref['n_bpms']}\n")
    else:
        results.append(None)
    if payload["n_quads"] != ref["n_quads"]:
        results.append(f"n_quads mismatch: {payload['n_quads']} vs {ref['n_quads']}\n")
    else:
        results.append(None)

    # Single-array comparisons
    cmp_array("mx", mx, ref["mx"], results, tolerance=1e-10, tolerance_type='absolute')
    cmp_array("my", my, ref["my"], results, tolerance=1e-10, tolerance_type='absolute')
    cmp_array("s",  s,  ref["s"],  results, tolerance=1e-12, tolerance_type='absolute')

    assert check_result(results)


def test_list_quads_bpms(lattice, update_ref_values=False):
    """
    Collect all Quadrupole and BPM elements whose IDs contain '.SA1',
    compare them against a stored reference list.
    Compatible with your check_result() helper.
    """
    quads, bpms = [], []
    for e in lattice.sequence:
        eid = getattr(e, "id", "")
        if isinstance(e, Quadrupole) and ".SA1" in eid:
            quads.append(e)
        elif isinstance(e, Monitor) and ".SA1" in eid:
            bpms.append(e)

    quad_ids = [q.id for q in quads]
    bpm_ids  = [b.id for b in bpms]

    payload = {
        "quad_ids": quad_ids,
        "bpm_ids": bpm_ids,
        "n_quads": len(quad_ids),
        "n_bpms": len(bpm_ids),
    }

    ref_path = os.path.join(REF_RES_DIR, sys._getframe().f_code.co_name + ".json")

    if update_ref_values:
        # Called from the update routine; return JSON-serializable payload
        return payload

    # --- regular comparison mode ---
    ref = json_read(ref_path)
    results = []

    # 1. number of quads
    if payload["n_quads"] != ref["n_quads"]:
        results.append(
            f"Number of quadrupoles mismatch: {payload['n_quads']} vs {ref['n_quads']}\n"
        )
    else:
        results.append(None)

    # 2. number of BPMs
    if payload["n_bpms"] != ref["n_bpms"]:
        results.append(
            f"Number of BPMs mismatch: {payload['n_bpms']} vs {ref['n_bpms']}\n"
        )
    else:
        results.append(None)

    # 3. quad id list equality
    if payload["quad_ids"] != ref["quad_ids"]:
        results.append("Quadrupole ID list differs (order or membership)\n")
    else:
        results.append(None)

    # 4. bpm id list equality
    if payload["bpm_ids"] != ref["bpm_ids"]:
        results.append("BPM ID list differs (order or membership)\n")
    else:
        results.append(None)

    # The helper prints info and returns True/False
    assert check_result(results)


def test_response_matrices(lattice, update_ref_values: bool = False):
    # --- pick elements
    quads, bpms = [], []
    for e in lattice.sequence:
        eid = getattr(e, "id", "")
        if isinstance(e, Quadrupole) and ".SA1" in eid:
            quads.append(e)
        elif isinstance(e, Monitor) and ".SA1" in eid:
            bpms.append(e)

    # --- compute
    Eref = 14.0
    energies = [16.0, 14.0, 10.0]
    tws0 = None  # or import from your context

    Rxs, Rys, Pxs, Pys = bba.generate_response_matrices_for_energies(
        lattice, quads, bpms, energies, Eref, plot=False, tws0=tws0
    )

    # --- serialize for update mode
    payload = {
        "Eref": Eref,
        "energies": energies,
        "n_bpms": len(bpms),
        "n_quads": len(quads),
        "Rxs": [numpy2json(m) for m in Rxs],
        "Rys": [numpy2json(m) for m in Rys],
        "Pxs": [numpy2json(m) for m in Pxs],
        "Pys": [numpy2json(m) for m in Pys],
    }

    ref_path = _ref_path(sys._getframe().f_code.co_name)
    if update_ref_values:
        return payload

    # --- compare
    ref = json_read(ref_path)

    results = []


    # basic metadata checks (optional; remove if not needed)
    if payload["Eref"] != ref["Eref"]:
        results.append(f"Eref mismatch: {payload['Eref']} vs {ref['Eref']}\n")
    else:
        results.append(None)

    if payload["energies"] != ref["energies"]:
        results.append("Energies list mismatch (order or values)\n")
    else:
        results.append(None)

    if payload["n_bpms"] != ref["n_bpms"]:
        results.append(f"n_bpms mismatch: {payload['n_bpms']} vs {ref['n_bpms']}\n")
    else:
        results.append(None)

    if payload["n_quads"] != ref["n_quads"]:
        results.append(f"n_quads mismatch: {payload['n_quads']} vs {ref['n_quads']}\n")
    else:
        results.append(None)

    # matrix comparisons (simple loops, clear failure locations)
    cmp_lists_of_mats("Rx", Rxs, ref["Rxs"], energies, results)
    cmp_lists_of_mats("Ry", Rys, ref["Rys"], energies, results)
    cmp_lists_of_mats("Px", Pxs, ref["Pxs"], energies, results)
    cmp_lists_of_mats("Py", Pys, ref["Pys"], energies, results)

    assert check_result(results)


def test_bpm_read_vs_energy(lattice, update_ref_values: bool = False):

    # --- pick elements
    quads, bpms = [], []
    for e in lattice.sequence:
        eid = getattr(e, "id", "")
        if isinstance(e, Quadrupole) and ".SA1" in eid:
            quads.append(e)
        elif isinstance(e, Monitor) and ".SA1" in eid:
            bpms.append(e)

    Xinit = (1e-5, -3e-6)
    Yinit = (-1e-5, 2e-6)

    # --- compute
    Eref = 14.0
    energies = [16.0, 14.0, 10.0]
    tws0 = None  # or import from your context

    Mx, My = bba.read_bpm_trajectories_vs_energy(
        lattice,
        quads,
        bpms,
        energies,
        Eref,
        Xinit=Xinit,
        Yinit=Yinit,
        bpm_offset_x=None,
        bpm_offset_y=None,
        noise_rms=(0e-6, 0e-6),  # BPM accuracy
        noise_truncated=3,
        plot=False,
    )


    # --- serialize for update mode
    payload = {
        "Eref": Eref,
        "energies": energies,
        "n_bpms": len(bpms),
        "n_quads": len(quads),
        "Mx": [numpy2json(m) for m in Mx],
        "My": [numpy2json(m) for m in My],
    }

    ref_path = _ref_path(sys._getframe().f_code.co_name)
    if update_ref_values:
        return payload

    # --- compare
    ref = json_read(ref_path)

    results = []

    if payload["Eref"] != ref["Eref"]:
        results.append(f"Eref mismatch: {payload['Eref']} vs {ref['Eref']}\n")
    else:
        results.append(None)

    if payload["energies"] != ref["energies"]:
        results.append("Energies list mismatch (order or values)\n")
    else:
        results.append(None)

    if payload["n_bpms"] != ref["n_bpms"]:
        results.append(f"n_bpms mismatch: {payload['n_bpms']} vs {ref['n_bpms']}\n")
    else:
        results.append(None)

    if payload["n_quads"] != ref["n_quads"]:
        results.append(f"n_quads mismatch: {payload['n_quads']} vs {ref['n_quads']}\n")
    else:
        results.append(None)

    # matrix comparisons (simple loops, clear failure locations)
    cmp_lists_of_mats("Mx", Mx, ref["Mx"], energies, results)
    cmp_lists_of_mats("My", My, ref["My"], energies, results)

    assert check_result(results)





def setup_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### DBA START ###\n\n')
    f.close()


def teardown_module(module):

    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### DBA END ###\n\n\n')
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
    update_functions.append('test_twiss')
    update_functions.append("test_read_orbit")
    update_functions.append('test_list_quads_bpms')
    update_functions.append('test_response_matrices')
    update_functions.append("test_bpm_read_vs_energy")
    
    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, True)
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
