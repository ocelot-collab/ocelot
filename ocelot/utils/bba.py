from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union
from pathlib import Path
import hashlib
import json
import pickle
import numpy as np
from ocelot import MagneticLattice, Quadrupole, Monitor, Twiss, twiss
from ocelot.gui import plot_opt_func
from ocelot.cpbd.tm_utils import transfer_maps_mult
from scipy.stats import truncnorm
from ocelot import Particle, lattice_track
import copy
import matplotlib.pyplot as plt
from numpy.typing import NDArray

def truncated_normal(mean: float, sigma: float, clip: float = 3.0) -> float:
    a, b = (-clip, clip)
    return truncnorm.rvs(a, b, loc=mean, scale=sigma)


def quad_response_matrix(
    lat: MagneticLattice,
    quads: List[Quadrupole],
    bpms: List[Monitor],
    energy: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the quadrupole-to-BPM response matrices for beam-based alignment (BBA).

    This function evaluates the linear response of BPM readings to transverse
    orbit kicks induced by quadrupole misalignments in a given magnetic lattice.
    The implementation follows the standard formalism used at LCLS, FLASH,
    SwissFEL, and EuXFEL (see P. Emma et al., *NIM A 429 (1999) 407–413*).

    For each quadrupole *j* and BPM *i*, the response coefficients are
    evaluated using the local thick-lens transport matrices:

        Px[i, j] = (1 − Q₁₁) R₁₁ − Q₂₁ R₁₂
        Py[i, j] = (1 − Q₃₃) R₃₃ − Q₄₃ R₃₄

    where:
        • Q is the 4×4 transverse transfer matrix across quadrupole *j*
        • R is the transfer matrix from the exit of *j* to BPM *i*
        • Indices correspond to (x, px, y, py) coordinates

    Parameters
    ----------
    lat : MagneticLattice
        Full Ocelot lattice containing all beamline elements.
    quads : list[Quadrupole]
        Ordered list of quadrupole magnets included in the alignment section.
    bpms : list[Monitor]
        Ordered list of beam position monitors used for orbit measurements.
     energy : float
        Initial beam energy used for evaluating element transfer maps (typically in GeV).

    Returns
    -------
    Px : np.ndarray, shape (n_bpm, n_quad)
        Horizontal response matrix between quadrupole offsets and BPM readings.
    Py : np.ndarray, shape (n_bpm, n_quad)
        Vertical response matrix between quadrupole offsets and BPM readings.

    Notes
    -----
    • Only BPMs located downstream of a given quadrupole contribute to that
      column of the response matrix.
    • The optics is currently evaluated at a fixed reference beam energy
      (14 GeV by default).
    • These matrices form the building blocks of the full BBA linear system
      used to extract BPM and quadrupole misalignments from multi-energy orbit data.
    """

    n = len(bpms)
    k = len(quads)
    Px = np.zeros((n, k))
    Py = np.zeros((n, k))
    for j, q in enumerate(quads):
        Q = q.R(energy=energy)[0]
        for i, bpm in enumerate(bpms):
            start_indx = lat.sequence.index(q)
            if lat.sequence.index(bpm) < start_indx:
                continue
            temp_lat = MagneticLattice(lat.sequence[start_indx+1:], stop=bpm)
            _, R, _ = temp_lat.transfer_maps(energy=energy)
            Px[i, j] = (1 - Q[0,0]) * R[0, 0] - Q[1, 0] * R[0, 1]
            Py[i, j] = (1 - Q[2, 2]) * R[2, 2] - Q[3, 2] * R[2, 3]
    return Px, Py


def launch_response_matrix(
    lat: MagneticLattice,
    bpms: List[Monitor],
    energy: float,
    ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the launch response matrices between the initial beam parameters
    (x0, x0′, y0, y0′) and BPM readings along the lattice.

    This function evaluates the linear transport of the initial orbit through
    the lattice and constructs the 2×N_BPM matrices `Rx` and `Ry`, which map
    the initial horizontal and vertical coordinates to the measured beam
    positions at each BPM:

        x_i = Rx[i, 0] * x0 + Rx[i, 1] * x0′
        y_i = Ry[i, 0] * y0 + Ry[i, 1] * y0′

    For each BPM, the accumulated transfer maps (`B`, `R`, `T`) from the
    lattice entrance up to the BPM are multiplied and then reduced to the
    relevant linear terms from the 6×6 transfer matrix:
        Rx[i, 0] = R[0, 0],   Rx[i, 1] = R[0, 1]
        Ry[i, 0] = R[2, 2],   Ry[i, 1] = R[2, 3]
    assuming the standard coordinate order (x, px, y, py, z, δ).

    Parameters
    ----------
    lat : MagneticLattice
        Complete Ocelot magnetic lattice containing the beamline elements.
    bpms : list[Monitor]
        Ordered list of BPMs for which the response is computed.
    energy : float
        Initial beam energy used for evaluating element transfer maps (typically in GeV).

    Returns
    -------
    Rx : np.ndarray, shape (n_bpm, 2)
        Horizontal launch response matrix relating (x0, x0′) to BPM readouts.
    Ry : np.ndarray, shape (n_bpm, 2)
        Vertical launch response matrix relating (y0, y0′) to BPM readouts.

    Notes
    -----
    • The accumulated transfer maps are reset after each BPM to store the
      propagation up to that location individually.
    • The beam energy is incremented at each transport step by the value
      returned from `tm.get_delta_e()`.
    • The matrices `Rx` and `Ry` are later used as part of the global
      beam-based alignment system, together with the quadrupole–BPM
      response matrices, to solve for initial orbit, BPM offsets, and
      quadrupole misalignments.
    • The implementation assumes a linear transport model; second-order
      tensors (`T`) are propagated for completeness but not used here.

    Raises
    ------
    Relies on valid Ocelot element interfaces (`.R(E)`, `.B(E)`, `.T(E)`, `.tms`);
    invalid or non-standard lattice elements may raise errors upstream.
    """
    Rx = np.zeros((len(bpms), 2))
    Ry = np.zeros((len(bpms), 2))
    Ra = np.eye(6)
    Ta = np.zeros((6, 6, 6))
    Ba = np.zeros((6, 1))
    E = energy
    i = 0
    for elem in lat.sequence:
        for Rb, Bb, Tb, tm in zip(elem.R(E), elem.B(E), elem.T(E), elem.tms):
            Ba, Ra, Ta = transfer_maps_mult(Ba, Ra, Ta, Bb, Rb, Tb)
            E += tm.get_delta_e()
        if elem in bpms:
            Rx[i, 0], Rx[i, 1], Ry[i, 0], Ry[i, 1] = Ra[0, 0], Ra[0, 1], Ra[2, 2], Ra[2, 3]
            i += 1
    return Rx, Ry


def read_orbit(
    lat: MagneticLattice,
    bpms: List[Monitor],
    pinit: Particle,
    bpm_offset_x: Optional[Sequence[float]] = None,
    bpm_offset_y: Optional[Sequence[float]] = None,
    noise_rms: Tuple[float, float] = (0.0, 0.0),
    noise_truncated: float = 3.0,
    launch_jitter: Tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0),
    launch_jitter_truncated: float = 3.0,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Track a single macroparticle through a lattice and simulate BPM readbacks
    with offsets and truncated Gaussian noise.

    The function propagates a single `Particle` through `lat`, records its (x, y)
    coordinates at each BPM, and applies BPM offsets and measurement noise to
    emulate realistic orbit readouts.

    Physics convention
    ------------------
    Each BPM measures beam position relative to its electrical center:
        m_x = x_true − b_x + n_x,
        m_y = y_true − b_y + n_y,
    where:
        • b_x, b_y — BPM offsets (its center position w.r.t. design orbit, in meters)
        • n_x, n_y — measurement noise, truncated Gaussian with σ = noise_rms

    Thus, a positive BPM offset (center shifted toward +x) produces a *negative*
    reading for a perfectly centered beam.

    Parameters
    ----------
    lat : MagneticLattice
        Complete Ocelot lattice containing all beamline elements.
    bpms : list[Monitor]
        Ordered list of BPM elements where readings are sampled
        (must be present in `lat.sequence`).
    pinit : Particle
        Initial particle; will be deep-copied so that the original is not modified.
    bpm_offset_x, bpm_offset_y : Optional[Sequence[float]]
        Per-BPM offsets (in meters). Length must match `len(bpms)` if provided.
        Offsets are subtracted from the true position (`m = x_true − b`).
    noise_rms : (float, float), default (0.0, 0.0)
        RMS of Gaussian measurement noise in (x, y) [m].
        Use (0, 0) to disable noise.
    noise_truncated : float, default 3.0
        Truncation level (in σ) for the Gaussian noise.
        For example, `3.0` keeps samples within ±3σ.
    launch_jitter : (float, float, float, float), default (0.0, 0.0, 0.0, 0.0)
        RMS of Gaussian launch orbit jitter in (x [m], x'[rad], y[m], y'[rad]).
        Use (0, 0, 0, 0) to disable jitter.
    launch_jitter_truncated : float, default 3.0
        Truncation level (in σ) for the Gaussian launch orbit jitter.
        For example, `3.0` keeps samples within ±3σ.

    Returns
    -------
    mx : np.ndarray
        Array of horizontal BPM readings (in meters), shape (n_bpm,).
    my : np.ndarray
        Array of vertical BPM readings (in meters), shape (n_bpm,).

    Raises
    ------
    ValueError
        If BPM offsets have incorrect length or BPMs are not found in the lattice.

    Notes
    -----
    - Each element's internal transfer maps (`elem.tms`) are applied in sequence
      using `tm.apply([p])`, assuming in-place updates of the particle state.
    - The simulation represents the single-particle (centroid) orbit, ignoring
      collective or statistical effects.
    - The resulting `(mx, my)` arrays can be used directly in BBA or DFS studies,
      e.g., for generating synthetic multi-energy orbit data following
      Emma et al. (1999).

    References
    ----------
    P. Emma, R. Carr, H.-D. Nuhn, “Beam-based alignment for the LCLS FEL undulator,”
    Nucl. Instrum. Methods Phys. Res. A 429 (1999) 407–413.
    """
    # --- Validations
    if bpm_offset_x is not None and len(bpm_offset_x) != len(bpms):
        raise ValueError("Length of `bpm_offset_x` must match length of `bpms`.")
    if bpm_offset_y is not None and len(bpm_offset_y) != len(bpms):
        raise ValueError("Length of `bpm_offset_y` must match length of `bpms`.")

    seq_set = set(lat.sequence)
    missing = [b for b in bpms if b not in seq_set]
    if missing:
        ids = [getattr(b, "id", repr(b)) for b in missing]
        raise ValueError(f"BPMs not found in lattice sequence: {ids}")

    # Map BPM -> index for O(1) lookup
    bpm_idx = {b: i for i, b in enumerate(bpms)}

    # Copy initial particle (avoid modifying caller’s instance)
    p = copy.deepcopy(pinit)
    p.x = p.x + truncated_normal(0.0, launch_jitter[0], clip=launch_jitter_truncated)
    p.px = p.px + truncated_normal(0.0, launch_jitter[1], clip=launch_jitter_truncated)
    p.y = p.y + truncated_normal(0.0, launch_jitter[2], clip=launch_jitter_truncated)
    p.py = p.py + truncated_normal(0.0, launch_jitter[3], clip=launch_jitter_truncated)

    mx, my = [], []
    s = []
    sx, sy = noise_rms

    # --- Tracking

    for elem in lat.sequence:
        for tm in getattr(elem, "tms", ()):
            tm.apply([p])

        if elem in bpm_idx:
            i = bpm_idx[elem]
            p_bpm = copy.copy(p)

            # Apply BPM offsets (subtract)
            if bpm_offset_x is not None:
                p_bpm.x -= bpm_offset_x[i]
            if bpm_offset_y is not None:
                p_bpm.y -= bpm_offset_y[i]

            # Add truncated Gaussian noise
            if sx:
                p_bpm.x += truncated_normal(0.0, sx, clip=noise_truncated)
            if sy:
                p_bpm.y += truncated_normal(0.0, sy, clip=noise_truncated)

            mx.append(p_bpm.x)
            my.append(p_bpm.y)
            s.append(p_bpm.s)

    return np.array(mx), np.array(my), np.array(s)


def get_random_offsets(
    elements: Sequence[Union[Quadrupole, Monitor]],
    sigma_x: float = 100e-6,
    sigma_y: float = 100e-6,
    n_sigma: float = 3.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate random (x, y) offsets for a list of lattice elements
    (BPMs or quadrupoles) using a truncated Gaussian distribution.

    Each element receives independent horizontal and vertical offsets
    drawn from N(0, σ) truncated symmetrically to ±n_sigma·σ.

    Parameters
    ----------
    elements : Sequence[Quadrupole | Monitor]
        List of elements for which offsets should be generated.
        Only the length of the list is used.
    sigma_x, sigma_y : float, optional
        RMS offset in meters for horizontal and vertical planes.
        Default is 100 µm for both.
    n_sigma : float, optional
        Truncation half-width in multiples of sigma (default ±3σ).

    Returns
    -------
    x, y : np.ndarray
        Arrays of truncated Gaussian offsets for each element, in meters.

    Notes
    -----
    - This unified version can be used for both BPM and quadrupole offsets.
      The interpretation depends on context:
          • For Quadrupoles → mechanical misalignment
          • For BPMs → readout/electrical zero offset
    - Example:
        >>> x_off, y_off = get_random_offsets(bpms, sigma_x=50e-6, sigma_y=30e-6)
        >>> x_off.shape
        (len(bpms),)
    """
    n = len(elements)
    a, b = -n_sigma, n_sigma
    x = truncnorm.rvs(a, b, loc=0.0, scale=sigma_x, size=n)
    y = truncnorm.rvs(a, b, loc=0.0, scale=sigma_y, size=n)
    return x, y


def _twiss_signature(tws0: Optional[Twiss]) -> Optional[Dict[str, Optional[float]]]:
    if tws0 is None:
        return None
    keys = ("beta_x", "beta_y", "alpha_x", "alpha_y", "E", "Dx", "Dy", "Dxp", "Dyp")
    out: Dict[str, Optional[float]] = {}
    for key in keys:
        value = getattr(tws0, key, None)
        out[key] = None if value is None else float(value)
    return out


def _response_cache_payload(
    lat: MagneticLattice,
    quads: Sequence[Quadrupole],
    bpms: Sequence[Monitor],
    energies: Sequence[float],
    Eref: float,
    tws0: Optional[Twiss],
) -> Dict[str, Any]:
    lat_sequence_ids = []
    for i, elem in enumerate(lat.sequence):
        elem_id = getattr(elem, "id", "")
        if elem_id in ("", None):
            elem_id = f"{type(elem).__name__}_{i}"
        lat_sequence_ids.append(str(elem_id))
    return {
        "algo": "generate_response_matrices_for_energies_cached_v1",
        "lat_sequence_ids": lat_sequence_ids,
        "quad_ids": [str(getattr(q, "id", f"quad_{i}")) for i, q in enumerate(quads)],
        "quad_k1": [float(q.k1) for q in quads],
        "bpm_ids": [str(getattr(b, "id", f"bpm_{i}")) for i, b in enumerate(bpms)],
        "energies": [float(e) for e in energies],
        "Eref": float(Eref),
        "tws0": _twiss_signature(tws0),
    }


def generate_response_matrices_for_energies(
    lat: MagneticLattice,
    quads: List[Quadrupole],
    bpms: List[Monitor],
    energies: List[float],
    Eref: float,
    plot: bool = False,
    tws0: Twiss = None
) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
    """
    Generate response matrices (P, R) for multiple beam energies.

    Parameters
    ----------
    lat : MagneticLattice
        Lattice object.
    quads : list of Quadrupole
        Quadrupoles whose strengths will be scaled with energy.
    bpms : list of BPM
        BPM elements used for measurement.
    energies : list of float
        Energies (in GeV) for which to compute response matrices.
    Eref : float
        Reference energy (in GeV) corresponding to nominal quadrupole settings.
    plot : bool, optional
        If True, plot optics for each scaled configuration.
    tws0 : Twiss, optional
        Initial Twiss parameters for optics calculation.

    Returns
    -------
    Rx, Ry, Px, Py : tuple of list[np.ndarray]
        Lists of response matrices for each energy:
        - Rx, Ry : launch response matrices
        - Px, Py : quadrupole response matrices
    """
    Rx: List[np.ndarray] = []
    Ry: List[np.ndarray] = []
    Px: List[np.ndarray] = []
    Py: List[np.ndarray] = []

    k1_ref = [q.k1 for q in quads]

    for energy in energies:
        coeff = Eref / energy
        for i, q in enumerate(quads):
            q.k1 = k1_ref[i] * coeff

        if plot:
            tws = twiss(lat, tws0)
            plot_opt_func(lat, tws, title=f"Energy {energy} GeV")

        px, py = quad_response_matrix(lat, quads, bpms, energy=energy)
        rx, ry = launch_response_matrix(lat, bpms, energy=energy)

        Px.append(px)
        Py.append(py)
        Rx.append(rx)
        Ry.append(ry)

    # restore original strengths
    for i, q in enumerate(quads):
        q.k1 = k1_ref[i]

    return Rx, Ry, Px, Py


def generate_response_matrices_for_energies_cached(
    lat: MagneticLattice,
    quads: List[Quadrupole],
    bpms: List[Monitor],
    energies: List[float],
    Eref: float,
    plot: bool = False,
    tws0: Twiss = None,
    cache_dir: Union[str, Path] = ".cache",
    force: bool = False,
    verbose: bool = True,
) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
    """
    Generate response matrices for multiple energies using on-disk caching.

    This helper wraps `generate_response_matrices_for_energies()` and avoids
    recomputing matrices when the effective lattice/input signature is unchanged.

    Parameters
    ----------
    lat, quads, bpms, energies, Eref, plot, tws0
        Same as in `generate_response_matrices_for_energies`.
    cache_dir : str | pathlib.Path, optional
        Directory where cache files are stored. Default: ``.cache``.
    force : bool, optional
        If True, ignore existing cache and recompute.
    verbose : bool, optional
        If True, print cache hit/miss/save messages.

    Returns
    -------
    Rx, Ry, Px, Py : tuple of list[np.ndarray]
        Same output as `generate_response_matrices_for_energies`.

    Notes
    -----
    Cache files are Python pickles and should only be loaded from trusted paths.
    """
    payload = _response_cache_payload(
        lat=lat,
        quads=quads,
        bpms=bpms,
        energies=energies,
        Eref=Eref,
        tws0=tws0,
    )
    payload_str = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    cache_key = hashlib.sha256(payload_str.encode("utf-8")).hexdigest()[:16]
    cache_path = Path(cache_dir) / f"response_matrices_{cache_key}.pkl"

    if cache_path.exists() and not force:
        try:
            with cache_path.open("rb") as f:
                cached = pickle.load(f)
            needed = ("Rxs", "Rys", "Pxs", "Pys")
            if all(k in cached for k in needed):
                if verbose:
                    print(f"Loaded response matrices from cache: {cache_path}")
                return cached["Rxs"], cached["Rys"], cached["Pxs"], cached["Pys"]
            if verbose:
                print(f"Cache file missing required keys, recomputing: {cache_path}")
        except Exception as exc:
            if verbose:
                print(f"Failed to load cache ({cache_path}), recomputing: {exc}")

    if verbose:
        print("Cache miss or forced recompute. Building response matrices...")
    Rxs, Rys, Pxs, Pys = generate_response_matrices_for_energies(
        lat=lat,
        quads=quads,
        bpms=bpms,
        energies=energies,
        Eref=Eref,
        plot=plot,
        tws0=tws0,
    )

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = cache_path.with_suffix(cache_path.suffix + ".tmp")
    with tmp_path.open("wb") as f:
        pickle.dump(
            {
                "payload": payload,
                "Rxs": Rxs,
                "Rys": Rys,
                "Pxs": Pxs,
                "Pys": Pys,
            },
            f,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
    tmp_path.replace(cache_path)
    if verbose:
        print(f"Saved response matrices cache: {cache_path}")
    return Rxs, Rys, Pxs, Pys


def build_full_matrix(
    R,
    P,
    per_measurement_launch: bool = False,
):
    """
    Build the full design matrix A for the linear system:
        A @ params = X
    but this function ONLY returns A.

    Parameters
    ----------
    R : tuple(list[np.ndarray], list[np.ndarray])
        (Rx_sets, Ry_sets)
        Rx_sets: list of S arrays, each shape (Nbpm, 2)
        Ry_sets: list of S arrays, each shape (Nbpm, 2)
        Each row maps initial conditions [x0, x0'] or [y0, y0'] to BPM readouts.
    P : tuple(list[np.ndarray], list[np.ndarray])
        (Px_sets, Py_sets)
        Px_sets: list of S arrays, each shape (Nbpm, Nqx)
        Py_sets: list of S arrays, each shape (Nbpm, Nqy)
        Each row maps quadrupole offsets to BPM readouts.
    per_measurement_launch : bool, optional
        If False (default), use one shared launch vector (2 parameters per plane)
        across all measurements/energies.
        If True, use an independent launch vector for each measurement
        (2*S parameters per plane).

    Returns
    -------
    A : np.ndarray
        Block-diagonal design matrix for both planes:
            [ Ax   0 ]
            [  0  Ay ]
        where
            Ax = [Rx_all  Px_all  -I_bpm_repeated]
            Ay = [Ry_all  Py_all  -I_bpm_repeated]
        Dimensions:
            shared launch (per_measurement_launch=False):
                Ax: (S*Nbpm, 2 + Nqx + Nbpm)
                Ay: (S*Nbpm, 2 + Nqy + Nbpm)
                A : (2*S*Nbpm, 2+Nqx+Nbpm + 2+Nqy+Nbpm)
            per-measurement launch (per_measurement_launch=True):
                Ax: (S*Nbpm, 2*S + Nqx + Nbpm)
                Ay: (S*Nbpm, 2*S + Nqy + Nbpm)
                A : (2*S*Nbpm, 2*S+Nqx+Nbpm + 2*S+Nqy+Nbpm)
    """
    Rx_sets, Ry_sets = R
    Px_sets, Py_sets = P

    # --- basic consistency checks (kept minimal) ---
    Sx = len(Rx_sets)
    Sy = len(Ry_sets)
    if not (Sx == len(Px_sets) == Sy == len(Py_sets)):
        raise ValueError("Inconsistent number of runs between R and P for x/y.")

    Nbpm = Rx_sets[0].shape[0]  # correct if all runs use the same BPM set
    # enforce same Nbpm per run
    if any(r.shape[0] != Nbpm or r.shape[1] != 2 for r in Rx_sets):
        raise ValueError("All Rx in Rx_sets must have shape (Nbpm, 2).")
    if any(r.shape[0] != Nbpm or r.shape[1] != 2 for r in Ry_sets):
        raise ValueError("All Ry in Ry_sets must have shape (Nbpm, 2).")
    if any(p.shape[0] != Nbpm for p in Px_sets):
        raise ValueError("All Px in Px_sets must have Nbpm rows.")
    if any(p.shape[0] != Nbpm for p in Py_sets):
        raise ValueError("All Py in Py_sets must have Nbpm rows.")

    Nqx = Px_sets[0].shape[1]
    Nqy = Py_sets[0].shape[1]
    if any(p.shape[1] != Nqx for p in Px_sets):
        raise ValueError("All Px in Px_sets must have same number of columns.")
    if any(p.shape[1] != Nqy for p in Py_sets):
        raise ValueError("All Py in Py_sets must have same number of columns.")

    def _assemble_launch_block(R_sets, S):
        if not per_measurement_launch:
            return np.vstack(R_sets)  # (S*Nbpm, 2)
        launch_block = np.zeros((S * Nbpm, 2 * S))  # (S*Nbpm, 2*S)
        for i, Ri in enumerate(R_sets):
            row0 = i * Nbpm
            launch_block[row0: row0 + Nbpm, 2 * i: 2 * i + 2] = Ri
        return launch_block

    # --- stack response matrices over runs ---
    Rx_all = _assemble_launch_block(Rx_sets, Sx)
    Ry_all = _assemble_launch_block(Ry_sets, Sy)
    Px_all = np.vstack(Px_sets)          # (S*Nbpm, Nqx)
    Py_all = np.vstack(Py_sets)          # (S*Nbpm, Nqy)

    # --- BPM offset blocks: -I repeated S times ---
    I_bpm = -np.eye(Nbpm)
    Bx_all = np.vstack([I_bpm for _ in range(Sx)])  # (S*Nbpm, Nbpm)
    By_all = np.vstack([I_bpm for _ in range(Sy)])  # (S*Nbpm, Nbpm)

    # --- plane blocks ---
    Ax = np.hstack([Rx_all, Px_all, Bx_all])        # (S*Nbpm, 2 + Nqx + Nbpm)
    Ay = np.hstack([Ry_all, Py_all, By_all])        # (S*Nbpm, 2 + Nqy + Nbpm)

    # --- final block-diagonal matrix ---
    Zxy = np.zeros((Ax.shape[0], Ay.shape[1]))
    Zyx = np.zeros((Ay.shape[0], Ax.shape[1]))
    A = np.block([[Ax, Zxy],
                  [Zyx, Ay]])

    return A


FloatArrayLike = Union[float, Sequence[float], NDArray[np.floating]]

def read_bpm_trajectories_vs_energy(
    lat: MagneticLattice,
    quads: List[Quadrupole],
    bpms: List[Monitor],
    energies: List[float],
    Eref: float,
    Xinit: Tuple[float, float] = (0, 0),
    Yinit: Tuple[float, float] = (0, 0),
    bpm_offset_x: Optional[FloatArrayLike] = None,
    bpm_offset_y: Optional[FloatArrayLike] = None,
    noise_rms: Tuple[float, float] = (0.0, 0.0),
    noise_truncated: Optional[float] = 3,
    launch_jitter: Tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0),
    launch_jitter_truncated: float = 3.0,
    plot: bool = True,
)-> Tuple[List[NDArray[np.floating]], List[NDArray[np.floating]]]:
    """
    Read BPM orbits for multiple energies and (optionally) plot two figures: x(s) and y(s).
    Returns only BPM readings: Mx, My (lists with one array per energy).

    Parameters
    ----------
    lat : MagneticLattice
    quads : list
        Quadrupole elements whose k1 must be scaled with energy.
    bpms : list
        BPM elements (order defines s_bpm).
    energies : iterable of float
        Beam energies [GeV] to scan.
    Eref : float
        Reference energy [GeV] for which current k1 values are defined.
    Xinit, Yinit : (x, px), (y, py)
        Initial conditions at lattice start (meters, radians).
    bpm_offset_x, bpm_offset_y : array-like or None
        BPM offsets passed to bba.read_orbit (meters). If None, zeros.
    noise_rms : (sx, sy)
        BPM readout noise RMS in meters.
    noise_truncated : float or None
        Truncation in sigmas for noise (e.g., 3).
    plot : bool
        If True, create exactly two plots: x(s) and y(s), with BPM markers.
    launch_jitter : (float, float, float, float), default (0.0, 0.0, 0.0, 0.0)
        RMS of Gaussian launch orbit jitter in (x [m], x'[rad], y[m], y'[rad]).
        Use (0, 0, 0, 0) to disable jitter.
    launch_jitter_truncated : float, default 3.0
        Truncation level (in σ) for the Gaussian launch orbit jitter.
        For example, `3.0` keeps samples within ±3σ.
    Returns
    -------
    Mx, My : list[np.ndarray], list[np.ndarray]
        Lists of BPM readings (per energy) in meters.
    """

    # save original strengths
    k1_ref = [q.k1 for q in quads]

    # containers
    S, Xtraj, Ytraj = [], [], []
    Mx_list, My_list = [], []
    s_bpm_ref = None

    energies = list(energies)

    try:
        for energy in energies:
            # scale quadrupoles for this energy
            coeff = Eref / energy
            for i, q in enumerate(quads):
                q.k1 = k1_ref[i] * coeff

            # track a representative particle for continuous trajectory
            pinit = Particle(x=Xinit[0], px=Xinit[1], y=Yinit[0], py=Yinit[1])
            p_list = lattice_track(lat, p=copy.copy(pinit))
            s = np.array([p.s for p in p_list])
            x_traj = np.array([p.x for p in p_list])
            y_traj = np.array([p.y for p in p_list])

            # read BPM "measurements"
            mx, my, s_bpm = read_orbit(
                lat,
                bpms=bpms,
                pinit=copy.copy(pinit),
                bpm_offset_x=bpm_offset_x,
                bpm_offset_y=bpm_offset_y,
                noise_rms=noise_rms,
                noise_truncated=noise_truncated,
                launch_jitter=launch_jitter,
                launch_jitter_truncated=launch_jitter_truncated,
            )

            # store
            S.append(s)
            Xtraj.append(x_traj)
            Ytraj.append(y_traj)
            Mx_list.append(np.asarray(mx))
            My_list.append(np.asarray(my))

            if s_bpm_ref is None:
                s_bpm_ref = np.asarray(s_bpm)
            else:
                if len(s_bpm_ref) != len(s_bpm) or np.max(np.abs(s_bpm_ref - s_bpm)) > 1e-9:
                    raise ValueError("BPM s-positions differ between energy settings.")

        # plotting: EXACTLY TWO FIGURES (mm units)
        if plot:
            scale = 1e3  # plot in mm
            # X-plot
            plt.figure()
            for i, (s, x, E) in enumerate(zip(S, Xtraj, energies)):
                color = f"C{i}"
                plt.plot(s, x * scale, lw=1.5, color=color, label=f"E={E:.3g} GeV")
                plt.plot(s_bpm_ref, Mx_list[i] * scale, "o", color=color, ms=4)
            plt.xlabel("s [m]")
            plt.ylabel("x(s) [mm]")
            plt.legend()
            plt.grid(True)
            plt.tight_layout()
            # Y-plot
            plt.figure()
            for i, (s, y, E) in enumerate(zip(S, Ytraj, energies)):
                color = f"C{i}"
                plt.plot(s, y * scale, lw=1.5, color=color, label=f"E={E:.3g} GeV")
                plt.plot(s_bpm_ref, My_list[i] * scale, "o", color=color, ms=4)
            plt.xlabel("s [m]")
            plt.ylabel("y(s) [mm]")
            plt.legend()
            plt.grid(True)
            plt.tight_layout()

        return Mx_list, My_list

    finally:
        # restore original strengths in any case
        for i, q in enumerate(quads):
            q.k1 = k1_ref[i]


def solve_svd(
    A: NDArray[np.floating],
    M: NDArray[np.floating],
    rcutoff: Optional[float] = None,
    print_spectrum: bool = True,
) -> NDArray[np.floating]:
    """
    Solve A * X = M using SVD with relative cutoff.

    Parameters
    ----------
    A : array, shape (m, n)
        Design matrix.
    M : array, shape (m,) or (m, k)
        Measurement vector(s). If 2D, solves for multiple RHS.
    rcutoff : float or None
        Relative cutoff on singular values (relative to max(S)).
        If None -> 1e-12.
    print_spectrum : bool
        Print singular values and cutoff info.

    Returns
    -------
    X_est : array, shape (n,) or (n, k)
        Estimated parameter vector(s), matching RHS multiplicity.
    """
    # SVD
    U, S, Vt = np.linalg.svd(A, full_matrices=False)  # U:(m,r), S:(r,), Vt:(r,n)

    if print_spectrum:
        print("---- SVD spectrum ----")
        print("Singular values:")
        print(S)
        print()

    # Relative cutoff
    if rcutoff is None:
        rcutoff = 1e-12
    Smax = np.max(S) if S.size else 0.0
    thresh = rcutoff * Smax

    # Mask of kept singular values
    keep = S > thresh
    nkeep = int(np.sum(keep))

    if print_spectrum:
        print(f"Cutoff (relative rcond) = {rcutoff}  -> absolute threshold = {thresh}")
        print(f"Keeping {nkeep} singular values out of {len(S)}")
        # Show inverted values (zeros where cut)
        Sinv_show = np.zeros_like(S)
        Sinv_show[keep] = 1.0 / S[keep]
        print("Inverted singular values after cutoff:")
        print(Sinv_show)
        print()

    # Build pseudoinverse: A_pinv = V * diag(1/S) * U^T (with cutoff)
    # Instead of forming A_pinv explicitly, apply it to M efficiently.
    # Step 1: U^T @ M
    Mt = M if M.ndim == 2 else M[:, None]              # shape (m, k)
    UtM = U.T @ Mt                                      # shape (r, k)

    # Step 2: diag(1/S) @ (U^T M), with cutoff
    Sinv = np.zeros_like(S)
    np.divide(1.0, S, out=Sinv, where=keep)            # zeros where not kept
    Sinv_UtM = (Sinv[:, None] * UtM)                   # shape (r, k)

    # Step 3: V @ (previous)
    X_est = (Vt.T @ Sinv_UtM)                          # shape (n, k)

    # Return 1D if input RHS was 1D
    if M.ndim == 1:
        X_est = X_est[:, 0]

    return X_est

def extract_solution(
    X_est: NDArray[np.floating],
    Nquad: int,
    Nbpm: int,
    n_measurements: int = 1,
    per_measurement_launch: bool = False,
) -> Tuple[
    NDArray[np.floating],  # Xinit_est
    NDArray[np.floating],  # Yinit_est
    NDArray[np.floating],  # qx_est
    NDArray[np.floating],  # qy_est
    NDArray[np.floating],  # bx_est
    NDArray[np.floating],  # by_est
]:
    """
    Split the SVD solution vector X_est into:
        X-launch, Y-launch, qx, qy, bx, by

    Structure of the solution vector:
        per_measurement_launch=False:
            X_est = [ Xinit(2), q_dx(Nq), b_dx(Nb),
                      Yinit(2), q_dy(Nq), b_dy(Nb) ]
        per_measurement_launch=True:
            X_est = [ Xinit(2*S), q_dx(Nq), b_dx(Nb),
                      Yinit(2*S), q_dy(Nq), b_dy(Nb) ]

    Parameters
    ----------
    X_est : np.ndarray
        Flattened solution vector from SVD solver.
    Nquad : int
        Number of quadrupoles (offset parameters).
    Nbpm : int
        Number of BPMs (offset parameters).
    n_measurements : int, optional
        Number of measurements/energies S. Used only when
        per_measurement_launch=True.
    per_measurement_launch : bool, optional
        If True, split launch blocks into S pairs and return arrays with shape
        (S, 2). If False, return launch vectors with shape (2,).

    Returns
    -------
    Xinit_est, Yinit_est, qx_est, qy_est, bx_est, by_est : np.ndarray
        Estimated components of the full solution.
    """
    if n_measurements < 1:
        raise ValueError("`n_measurements` must be >= 1.")

    n_launch = 2 * n_measurements if per_measurement_launch else 2
    expected_len = 2 * (n_launch + Nquad + Nbpm)
    if X_est.size != expected_len:
        msg = f"X_est length mismatch: expected {expected_len}, got {X_est.size}"
        if per_measurement_launch:
            msg += ". For per_measurement_launch=True, pass n_measurements=len(energies)."
            remainder = X_est.size - 2 * (Nquad + Nbpm)
            if remainder > 0 and remainder % 4 == 0:
                inferred = remainder // 4
                msg += f" This vector suggests n_measurements={inferred}."
        raise ValueError(
            msg
        )

    # X-plane
    Xinit_est = X_est[0:n_launch]
    qx_est = X_est[n_launch : n_launch + Nquad]
    bx_est = X_est[n_launch + Nquad : n_launch + Nquad + Nbpm]

    # Y-plane
    off = n_launch + Nquad + Nbpm
    Yinit_est = X_est[off : off + n_launch]
    qy_est = X_est[off + n_launch : off + n_launch + Nquad]
    by_est = X_est[off + n_launch + Nquad : off + n_launch + Nquad + Nbpm]

    if per_measurement_launch:
        Xinit_est = Xinit_est.reshape(n_measurements, 2)
        Yinit_est = Yinit_est.reshape(n_measurements, 2)

    return Xinit_est, Yinit_est, qx_est, qy_est, bx_est, by_est


if __name__ == "__main__":
    print("here")
