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

def test_quad_response_matrix(lattice, update_ref_values: bool = False):
    # --- pick elements (same SA1 logic as in other tests)
    quads, bpms = [], []
    for e in lattice.sequence:
        eid = getattr(e, "id", "")
        if isinstance(e, Quadrupole) and ".SA1" in eid:
            quads.append(e)
        elif isinstance(e, Monitor) and ".SA1" in eid:
            bpms.append(e)

    energy = 14.0  # GeV (or whatever reference you use elsewhere)

    Px, Py = bba.quad_response_matrix(
        lat=lattice,
        quads=quads,
        bpms=bpms,
        energy=energy,
    )

    # --- serialize for update mode
    payload = {
        "energy": energy,
        "n_bpms": len(bpms),
        "n_quads": len(quads),
        "quad_ids": [q.id for q in quads],
        "bpm_ids": [b.id for b in bpms],
        "Px": numpy2json(Px),   # 2D array
        "Py": numpy2json(Py),   # 2D array
    }

    ref_path = _ref_path(sys._getframe().f_code.co_name)
    if update_ref_values:
        return payload

    # --- compare with reference
    ref = json_read(ref_path)
    results = []

    # meta checks (very useful if lattice/topology changes)
    if payload["energy"] != ref["energy"]:
        results.append(
            f"energy mismatch: {payload['energy']} vs {ref['energy']}\n"
        )
    else:
        results.append(None)

    if payload["n_bpms"] != ref["n_bpms"]:
        results.append(
            f"n_bpms mismatch: {payload['n_bpms']} vs {ref['n_bpms']}\n"
        )
    else:
        results.append(None)

    if payload["n_quads"] != ref["n_quads"]:
        results.append(
            f"n_quads mismatch: {payload['n_quads']} vs {ref['n_quads']}\n"
        )
    else:
        results.append(None)

    # element IDs (helps to see if you accidentally picked different devices)
    if payload["quad_ids"] != ref.get("quad_ids", []):
        results.append(
            "quad_ids mismatch (order or membership differs from reference)\n"
        )
    else:
        results.append(None)

    if payload["bpm_ids"] != ref.get("bpm_ids", []):
        results.append(
            "bpm_ids mismatch (order or membership differs from reference)\n"
        )
    else:
        results.append(None)

    # numeric comparisons for Px, Py
    Px_ref = json2numpy(ref["Px"])
    Py_ref = json2numpy(ref["Py"])

    if Px.shape != Px_ref.shape:
        results.append(
            f"Px shape mismatch: {Px.shape} vs {Px_ref.shape}\n"
        )
    else:
        # check_matrix returns list[None | str]
        results.extend(
            check_matrix(
                Px,
                Px_ref,
                tolerance=1.0e-10,
                tolerance_type="absolute",
                assert_info="Px - ",
                numerical_zero=1.0e-15,
            )
        )

    if Py.shape != Py_ref.shape:
        results.append(
            f"Py shape mismatch: {Py.shape} vs {Py_ref.shape}\n"
        )
    else:
        results.extend(
            check_matrix(
                Py,
                Py_ref,
                tolerance=1.0e-10,
                tolerance_type="absolute",
                assert_info="Py - ",
                numerical_zero=1.0e-15,
            )
        )

    assert check_result(results)


def test_launch_response_matrix(lattice, update_ref_values: bool = False):
    # --- pick elements, same SA1 logic as in other tests ---
    quads, bpms = [], []
    for e in lattice.sequence:
        eid = getattr(e, "id", "")
        if isinstance(e, Quadrupole) and ".SA1" in eid:
            quads.append(e)
        elif isinstance(e, Monitor) and ".SA1" in eid:
            bpms.append(e)

    energy = 14.0  # GeV (keep consistent with other tests)

    Rx, Ry = bba.launch_response_matrix(
        lat=lattice,
        bpms=bpms,
        energy=energy,
    )

    # --- serialize for update mode ---
    payload = {
        "energy": energy,
        "n_bpms": len(bpms),
        "n_quads": len(quads),
        "bpm_ids": [b.id for b in bpms],
        # quads not used in this function, but stored for consistency/topology checks
        "quad_ids": [q.id for q in quads],
        "Rx": numpy2json(Rx),
        "Ry": numpy2json(Ry),
    }

    ref_path = _ref_path(sys._getframe().f_code.co_name)
    if update_ref_values:
        return payload

    # --- compare with reference ---
    ref = json_read(ref_path)
    results = []

    # meta checks
    if payload["energy"] != ref["energy"]:
        results.append(
            f"energy mismatch: {payload['energy']} vs {ref['energy']}\n"
        )
    else:
        results.append(None)

    if payload["n_bpms"] != ref["n_bpms"]:
        results.append(
            f"n_bpms mismatch: {payload['n_bpms']} vs {ref['n_bpms']}\n"
        )
    else:
        results.append(None)

    if payload["n_quads"] != ref["n_quads"]:
        results.append(
            f"n_quads mismatch: {payload['n_quads']} vs {ref['n_quads']}\n"
        )
    else:
        results.append(None)

    if payload["bpm_ids"] != ref.get("bpm_ids", []):
        results.append(
            "bpm_ids mismatch (order or membership differs from reference)\n"
        )
    else:
        results.append(None)

    if payload["quad_ids"] != ref.get("quad_ids", []):
        results.append(
            "quad_ids mismatch (order or membership differs from reference)\n"
        )
    else:
        results.append(None)

    # numeric comparisons
    Rx_ref = json2numpy(ref["Rx"])
    Ry_ref = json2numpy(ref["Ry"])

    if Rx.shape != Rx_ref.shape:
        results.append(
            f"Rx shape mismatch: {Rx.shape} vs {Rx_ref.shape}\n"
        )
    else:
        results.extend(
            check_matrix(
                Rx,
                Rx_ref,
                tolerance=1.0e-10,
                tolerance_type="absolute",
                assert_info="Rx - ",
                numerical_zero=1.0e-15,
            )
        )

    if Ry.shape != Ry_ref.shape:
        results.append(
            f"Ry shape mismatch: {Ry.shape} vs {Ry_ref.shape}\n"
        )
    else:
        results.extend(
            check_matrix(
                Ry,
                Ry_ref,
                tolerance=1.0e-10,
                tolerance_type="absolute",
                assert_info="Ry - ",
                numerical_zero=1.0e-15,
            )
        )

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


def test_generate_response_matrices_for_energies_cached(tmp_path, monkeypatch):
    class _Elem:
        def __init__(self, eid: str, k1: float = 0.0):
            self.id = eid
            self.k1 = k1

    q1 = _Elem("Q1", k1=1.0)
    q2 = _Elem("Q2", k1=-2.0)
    b1 = _Elem("B1")
    b2 = _Elem("B2")
    lat = type("Lat", (), {"sequence": [q1, b1, q2, b2]})()
    quads = [q1, q2]
    bpms = [b1, b2]
    energies = [10.0, 12.0]
    Eref = 14.0

    expected = (
        [np.full((2, 2), 1.0), np.full((2, 2), 2.0)],  # Rxs
        [np.full((2, 2), 3.0), np.full((2, 2), 4.0)],  # Rys
        [np.full((2, 2), 5.0), np.full((2, 2), 6.0)],  # Pxs
        [np.full((2, 2), 7.0), np.full((2, 2), 8.0)],  # Pys
    )
    calls = {"n": 0}

    def _fake_generate(
        lat,
        quads,
        bpms,
        energies,
        Eref,
        plot=False,
        tws0=None,
    ):
        calls["n"] += 1
        return expected

    monkeypatch.setattr(bba, "generate_response_matrices_for_energies", _fake_generate)

    out1 = bba.generate_response_matrices_for_energies_cached(
        lat=lat,
        quads=quads,
        bpms=bpms,
        energies=energies,
        Eref=Eref,
        cache_dir=tmp_path,
        force=False,
        verbose=False,
    )
    assert calls["n"] == 1

    # Second call with identical signature should come from cache.
    out2 = bba.generate_response_matrices_for_energies_cached(
        lat=lat,
        quads=quads,
        bpms=bpms,
        energies=energies,
        Eref=Eref,
        cache_dir=tmp_path,
        force=False,
        verbose=False,
    )
    assert calls["n"] == 1

    # Force=True must recompute.
    out3 = bba.generate_response_matrices_for_energies_cached(
        lat=lat,
        quads=quads,
        bpms=bpms,
        energies=energies,
        Eref=Eref,
        cache_dir=tmp_path,
        force=True,
        verbose=False,
    )
    assert calls["n"] == 2

    # Changing signature (energies) must create a different cache key and recompute.
    _ = bba.generate_response_matrices_for_energies_cached(
        lat=lat,
        quads=quads,
        bpms=bpms,
        energies=[10.0, 13.0],
        Eref=Eref,
        cache_dir=tmp_path,
        force=False,
        verbose=False,
    )
    assert calls["n"] == 3

    for out in (out1, out2, out3):
        for got_group, exp_group in zip(out, expected):
            assert len(got_group) == len(exp_group)
            for got, exp in zip(got_group, exp_group):
                np.testing.assert_allclose(got, exp)

    assert any(tmp_path.glob("response_matrices_*.pkl"))


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


def test_build_full_matrix_shape_and_blocks():
    # S: number of energy settings / runs
    S = 3
    Nbpm = 4
    Nqx = 3
    Nqy = 3

    # Construct deterministic matrices so we can check blocks
    Rx_sets = []
    Ry_sets = []
    Px_sets = []
    Py_sets = []

    for s in range(S):
        # Rx, Ry: (Nbpm, 2)
        Rx_s = np.full((Nbpm, 2), fill_value=10 + s, dtype=float)
        Ry_s = np.full((Nbpm, 2), fill_value=20 + s, dtype=float)
        # Px: (Nbpm, Nqx), Py: (Nbpm, Nqy)
        Px_s = np.full((Nbpm, Nqx), fill_value=30 + s, dtype=float)
        Py_s = np.full((Nbpm, Nqy), fill_value=40 + s, dtype=float)

        Rx_sets.append(Rx_s)
        Ry_sets.append(Ry_s)
        Px_sets.append(Px_s)
        Py_sets.append(Py_s)

    A = bba.build_full_matrix(R=(Rx_sets, Ry_sets), P=(Px_sets, Py_sets))

    # Expected overall shape
    m_exp = 2 * S * Nbpm
    n_exp = (2 + Nqx + Nbpm) + (2 + Nqy + Nbpm)
    assert A.shape == (m_exp, n_exp)

    # Rebuild the pieces that build_full_matrix is supposed to form
    Rx_all = np.vstack(Rx_sets)        # (S*Nbpm, 2)
    Ry_all = np.vstack(Ry_sets)        # (S*Nbpm, 2)
    Px_all = np.vstack(Px_sets)        # (S*Nbpm, Nqx)
    Py_all = np.vstack(Py_sets)        # (S*Nbpm, Nqy)

    I_bpm = -np.eye(Nbpm)
    Bx_all = np.vstack([I_bpm for _ in range(S)])   # (S*Nbpm, Nbpm)
    By_all = np.vstack([I_bpm for _ in range(S)])   # (S*Nbpm, Nbpm)

    Ax_ref = np.hstack([Rx_all, Px_all, Bx_all])    # (S*Nbpm, 2+Nqx+Nbpm)
    Ay_ref = np.hstack([Ry_all, Py_all, By_all])    # (S*Nbpm, 2+Nqy+Nbpm)

    n_x = Ax_ref.shape[1]
    n_y = Ay_ref.shape[1]

    # Extract blocks from A
    Ax = A[: S * Nbpm, :n_x]
    Zxy = A[: S * Nbpm, n_x:]
    Zyx = A[S * Nbpm :, :n_x]
    Ay = A[S * Nbpm :, n_x:]

    # Check main blocks match what we expect
    np.testing.assert_allclose(Ax, Ax_ref)
    np.testing.assert_allclose(Ay, Ay_ref)

    # Off-diagonal blocks must be exactly zero
    assert np.allclose(Zxy, 0.0)
    assert np.allclose(Zyx, 0.0)


def test_build_full_matrix_per_measurement_launch_shape_and_blocks():
    # S: number of energy settings / runs
    S = 3
    Nbpm = 4
    Nqx = 3
    Nqy = 3

    # Construct deterministic matrices so we can check block placement
    Rx_sets = []
    Ry_sets = []
    Px_sets = []
    Py_sets = []

    for s in range(S):
        Rx_sets.append(np.full((Nbpm, 2), fill_value=10 + s, dtype=float))
        Ry_sets.append(np.full((Nbpm, 2), fill_value=20 + s, dtype=float))
        Px_sets.append(np.full((Nbpm, Nqx), fill_value=30 + s, dtype=float))
        Py_sets.append(np.full((Nbpm, Nqy), fill_value=40 + s, dtype=float))

    A = bba.build_full_matrix(
        R=(Rx_sets, Ry_sets),
        P=(Px_sets, Py_sets),
        per_measurement_launch=True,
    )

    # Expected overall shape with 2*S launch parameters per plane
    m_exp = 2 * S * Nbpm
    n_exp = (2 * S + Nqx + Nbpm) + (2 * S + Nqy + Nbpm)
    assert A.shape == (m_exp, n_exp)

    # Build expected launch blocks: block-diagonal in measurement index
    Lx_ref = np.zeros((S * Nbpm, 2 * S))
    Ly_ref = np.zeros((S * Nbpm, 2 * S))
    for s in range(S):
        r0 = s * Nbpm
        Lx_ref[r0 : r0 + Nbpm, 2 * s : 2 * s + 2] = Rx_sets[s]
        Ly_ref[r0 : r0 + Nbpm, 2 * s : 2 * s + 2] = Ry_sets[s]

    Px_all = np.vstack(Px_sets)
    Py_all = np.vstack(Py_sets)

    I_bpm = -np.eye(Nbpm)
    Bx_all = np.vstack([I_bpm for _ in range(S)])
    By_all = np.vstack([I_bpm for _ in range(S)])

    Ax_ref = np.hstack([Lx_ref, Px_all, Bx_all])
    Ay_ref = np.hstack([Ly_ref, Py_all, By_all])

    n_x = Ax_ref.shape[1]

    Ax = A[: S * Nbpm, :n_x]
    Zxy = A[: S * Nbpm, n_x:]
    Zyx = A[S * Nbpm :, :n_x]
    Ay = A[S * Nbpm :, n_x:]

    np.testing.assert_allclose(Ax, Ax_ref)
    np.testing.assert_allclose(Ay, Ay_ref)
    assert np.allclose(Zxy, 0.0)
    assert np.allclose(Zyx, 0.0)


def test_build_full_matrix_wrong_shapes_raises():
    S = 2
    Nbpm = 4

    # Correct Rx, Ry shapes
    Rx_sets = [np.zeros((Nbpm, 2)), np.zeros((Nbpm, 2))]
    Ry_sets = [np.zeros((Nbpm, 2)), np.zeros((Nbpm, 2))]

    # Px has wrong row count for one run
    Px_good = np.zeros((Nbpm, 3))
    Px_bad = np.zeros((Nbpm + 1, 3))   # wrong Nbpm
    Px_sets = [Px_good, Px_bad]

    Py_sets = [np.zeros((Nbpm, 2)), np.zeros((Nbpm, 2))]

    with pytest.raises(ValueError, match="All Px in Px_sets must have Nbpm rows"):
        bba.build_full_matrix(R=(Rx_sets, Ry_sets), P=(Px_sets, Py_sets))

    # Fix Px, but break Rx shape
    Px_sets = [Px_good, Px_good]
    Rx_bad = np.zeros((Nbpm + 1, 2))
    Rx_sets = [np.zeros((Nbpm, 2)), Rx_bad]

    with pytest.raises(ValueError, match="All Rx in Rx_sets must have shape"):
        bba.build_full_matrix(R=(Rx_sets, Ry_sets), P=(Px_sets, Py_sets))

def test_solve_svd_vector_rhs_exact_solution():
    """
    For a well-conditioned full-rank system (m > n), solve_svd should recover
    the exact x_true (within numerical precision) for a 1D RHS.
    """
    rng = np.random.default_rng(12345)

    m, n = 10, 4
    A = rng.normal(size=(m, n))
    # make A reasonably well-conditioned
    # (random Gaussian is usually fine for this kind of unit test)
    x_true = rng.normal(size=n)
    M = A @ x_true

    x_est = bba.solve_svd(A, M, rcutoff=1e-12, print_spectrum=False)

    assert x_est.shape == (n,)
    np.testing.assert_allclose(x_est, x_true, rtol=1e-10, atol=1e-12)


def test_solve_svd_multi_rhs_exact_solution():
    """
    For multiple RHS (2D M), solve_svd should solve all systems simultaneously.
    """
    rng = np.random.default_rng(54321)

    m, n, k = 12, 5, 3
    A = rng.normal(size=(m, n))
    X_true = rng.normal(size=(n, k))
    M = A @ X_true   # shape (m, k)

    X_est = bba.solve_svd(A, M, rcutoff=1e-12, print_spectrum=False)

    assert X_est.shape == (n, k)
    np.testing.assert_allclose(X_est, X_true, rtol=1e-10, atol=1e-12)


@pytest.mark.parametrize("rcutoff", [1e-2, 1e-4, 1e-8])
def test_solve_svd_matches_numpy_pinv_rcutoff(rcutoff):
    """
    Check that solve_svd(A, M, rcutoff) matches np.linalg.pinv(A, rcond=rcutoff)
    for both 1D and 2D RHS.
    """
    rng = np.random.default_rng(111)

    m, n, k = 8, 6, 2
    A = rng.normal(size=(m, n))

    # 1D RHS
    M_vec = rng.normal(size=m)
    X_ref_vec = np.linalg.pinv(A, rcond=rcutoff) @ M_vec
    X_est_vec = bba.solve_svd(A, M_vec, rcutoff=rcutoff, print_spectrum=False)

    assert X_est_vec.shape == (n,)
    np.testing.assert_allclose(X_est_vec, X_ref_vec, rtol=1e-10, atol=1e-12)

    # 2D RHS
    M_mat = rng.normal(size=(m, k))
    X_ref_mat = np.linalg.pinv(A, rcond=rcutoff) @ M_mat
    X_est_mat = bba.solve_svd(A, M_mat, rcutoff=rcutoff, print_spectrum=False)

    assert X_est_mat.shape == (n, k)
    np.testing.assert_allclose(X_est_mat, X_ref_mat, rtol=1e-10, atol=1e-12)


def test_extract_solution_splits_correctly():
    """
    Check that extract_solution returns the right slices for a known pattern.
    """
    Nquad = 3
    Nbpm = 4
    expected_len = 2 * (2 + Nquad + Nbpm)

    # Build X_est with an easy-to-check pattern: X_est[i] = i
    X_est = np.arange(expected_len, dtype=float)

    (
        Xinit_est,
        Yinit_est,
        qx_est,
        qy_est,
        bx_est,
        by_est,
    ) = bba.extract_solution(X_est, Nquad=Nquad, Nbpm=Nbpm)

    # Manual reference slices
    # X plane
    Xinit_ref = X_est[0:2]
    qx_ref    = X_est[2 : 2 + Nquad]
    bx_ref    = X_est[2 + Nquad : 2 + Nquad + Nbpm]

    off = 2 + Nquad + Nbpm
    Yinit_ref = X_est[off : off + 2]
    qy_ref    = X_est[off + 2 : off + 2 + Nquad]
    by_ref    = X_est[off + 2 + Nquad : off + 2 + Nquad + Nbpm]

    np.testing.assert_allclose(Xinit_est, Xinit_ref)
    np.testing.assert_allclose(Yinit_est, Yinit_ref)
    np.testing.assert_allclose(qx_est,    qx_ref)
    np.testing.assert_allclose(qy_est,    qy_ref)
    np.testing.assert_allclose(bx_est,    bx_ref)
    np.testing.assert_allclose(by_est,    by_ref)


def test_extract_solution_per_measurement_launch_splits_correctly():
    Nquad = 3
    Nbpm = 4
    S = 3
    n_launch = 2 * S
    expected_len = 2 * (n_launch + Nquad + Nbpm)

    X_est = np.arange(expected_len, dtype=float)

    (
        Xinit_est,
        Yinit_est,
        qx_est,
        qy_est,
        bx_est,
        by_est,
    ) = bba.extract_solution(
        X_est,
        Nquad=Nquad,
        Nbpm=Nbpm,
        n_measurements=S,
        per_measurement_launch=True,
    )

    Xinit_ref = X_est[0:n_launch].reshape(S, 2)
    qx_ref = X_est[n_launch : n_launch + Nquad]
    bx_ref = X_est[n_launch + Nquad : n_launch + Nquad + Nbpm]

    off = n_launch + Nquad + Nbpm
    Yinit_ref = X_est[off : off + n_launch].reshape(S, 2)
    qy_ref = X_est[off + n_launch : off + n_launch + Nquad]
    by_ref = X_est[off + n_launch + Nquad : off + n_launch + Nquad + Nbpm]

    np.testing.assert_allclose(Xinit_est, Xinit_ref)
    np.testing.assert_allclose(Yinit_est, Yinit_ref)
    np.testing.assert_allclose(qx_est, qx_ref)
    np.testing.assert_allclose(qy_est, qy_ref)
    np.testing.assert_allclose(bx_est, bx_ref)
    np.testing.assert_allclose(by_est, by_ref)


@pytest.mark.parametrize("Nquad, Nbpm", [(1, 1), (5, 2), (3, 7)])
def test_extract_solution_raises_on_wrong_length(Nquad, Nbpm):
    """
    If X_est.size != 2*(2+Nquad+Nbpm) the function must raise ValueError.
    """
    expected_len = 2 * (2 + Nquad + Nbpm)

    # too short
    X_short = np.zeros(expected_len - 1, dtype=float)
    with pytest.raises(ValueError, match="X_est length mismatch"):
        bba.extract_solution(X_short, Nquad=Nquad, Nbpm=Nbpm)

    # too long
    X_long = np.zeros(expected_len + 3, dtype=float)
    with pytest.raises(ValueError, match="X_est length mismatch"):
        bba.extract_solution(X_long, Nquad=Nquad, Nbpm=Nbpm)


def test_extract_solution_per_measurement_launch_raises_on_wrong_length():
    Nquad = 3
    Nbpm = 4
    S = 3
    n_launch = 2 * S
    expected_len = 2 * (n_launch + Nquad + Nbpm)

    X_short = np.zeros(expected_len - 1, dtype=float)
    with pytest.raises(ValueError, match="X_est length mismatch"):
        bba.extract_solution(
            X_short,
            Nquad=Nquad,
            Nbpm=Nbpm,
            n_measurements=S,
            per_measurement_launch=True,
        )

    with pytest.raises(ValueError, match="n_measurements"):
        bba.extract_solution(
            np.zeros(expected_len, dtype=float),
            Nquad=Nquad,
            Nbpm=Nbpm,
            n_measurements=0,
            per_measurement_launch=True,
        )


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
    update_functions.append('test_quad_response_matrix')
    update_functions.append("test_launch_response_matrix")
    update_functions.append('test_response_matrices')
    update_functions.append("test_bpm_read_vs_energy")
    
    if cmdopt in update_functions:
        result = eval(cmdopt)(lattice, True)
        
        if os.path.isfile(REF_RES_DIR + cmdopt + '.json'):
            os.rename(REF_RES_DIR + cmdopt + '.json', REF_RES_DIR + cmdopt + '.old')
        
        json_save(result, REF_RES_DIR + cmdopt + '.json')
