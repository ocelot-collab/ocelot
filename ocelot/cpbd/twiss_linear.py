import numpy as np

from ocelot.common.globals import m_e_GeV


_TWISS_STATE_DEFAULTS = {
    "beta_x": 0.0,
    "alpha_x": 0.0,
    "beta_y": 0.0,
    "alpha_y": 0.0,
    "Dx": 0.0,
    "Dxp": 0.0,
    "Dy": 0.0,
    "Dyp": 0.0,
    "mux": 0.0,
    "muy": 0.0,
    "E": 0.0,
    "s": 0.0,
    "p": 0.0,
}

_ZERO_TOL = 1.0e-10


def make_twiss_state(*, beta_x=0.0, alpha_x=0.0, beta_y=0.0, alpha_y=0.0,
                     Dx=0.0, Dxp=0.0, Dy=0.0, Dyp=0.0, mux=0.0, muy=0.0,
                     E=0.0, s=0.0, p=0.0, xp=np):
    """Create a backend-aware linear Twiss state dictionary."""
    return {
        "beta_x": xp.asarray(beta_x),
        "alpha_x": xp.asarray(alpha_x),
        "beta_y": xp.asarray(beta_y),
        "alpha_y": xp.asarray(alpha_y),
        "Dx": xp.asarray(Dx),
        "Dxp": xp.asarray(Dxp),
        "Dy": xp.asarray(Dy),
        "Dyp": xp.asarray(Dyp),
        "mux": xp.asarray(mux),
        "muy": xp.asarray(muy),
        "E": xp.asarray(E),
        "s": xp.asarray(s),
        "p": xp.asarray(p),
    }


def _normalize_twiss_state(state, xp=np):
    return {
        key: xp.asarray(state.get(key, default))
        for key, default in _TWISS_STATE_DEFAULTS.items()
    }


def _twiss_gamma(beta, alpha, xp=np):
    beta = xp.asarray(beta)
    alpha = xp.asarray(alpha)
    nonzero_beta = beta != 0
    beta_safe = xp.where(nonzero_beta, beta, xp.ones_like(beta))
    gamma = (1.0 + alpha * alpha) / beta_safe
    return xp.where(nonzero_beta, gamma, xp.zeros_like(gamma))


def _phase_advance_increment(numerator, denominator, xp=np):
    delta_mu = xp.arctan2(numerator, denominator)
    if xp is np:
        if delta_mu < 0:
            delta_mu += np.pi
        return delta_mu
    return xp.where(delta_mu < 0, delta_mu + xp.pi, delta_mu)


def _twiss_energy_scale(energy, delta_e, xp=np):
    if xp is np:
        if abs(delta_e) > _ZERO_TOL:
            return np.sqrt((energy + delta_e) / energy)
        return 1.0

    energy = xp.asarray(energy)
    delta_e = xp.asarray(delta_e)
    needs_scale = xp.abs(delta_e) > _ZERO_TOL
    energy_safe = xp.where(energy != 0, energy, xp.ones_like(energy))
    scale = xp.sqrt((energy + delta_e) / energy_safe)
    return xp.where(needs_scale, scale, xp.ones_like(scale))


def twiss_step_from_r(R, state, delta_e=0.0, length=0.0, xp=np):
    """
    Propagate a linear Twiss state through one 6x6 first-order matrix.

    Parameters
    ----------
    R
        6x6 linear transport matrix.
    state
        Mapping with the keys created by ``make_twiss_state``.
    delta_e
        Reference-energy change across the element. The first autodiff-focused
        path assumes passive lattices and typically leaves this at zero.
    length
        Path length increment added to ``state["s"]``.
    xp
        Backend namespace such as ``numpy`` or ``jax.numpy``.
    """
    state = _normalize_twiss_state(state, xp=xp)
    scale = _twiss_energy_scale(state["E"], delta_e, xp=xp)

    R00 = R[0, 0] * scale
    R01 = R[0, 1] * scale
    R10 = R[1, 0] * scale
    R11 = R[1, 1] * scale
    R22 = R[2, 2] * scale
    R23 = R[2, 3] * scale
    R32 = R[3, 2] * scale
    R33 = R[3, 3] * scale
    R05 = R[0, 5]
    R15 = R[1, 5]
    R25 = R[2, 5]
    R35 = R[3, 5]

    gamma_x = _twiss_gamma(state["beta_x"], state["alpha_x"], xp=xp)
    gamma_y = _twiss_gamma(state["beta_y"], state["alpha_y"], xp=xp)

    beta_x = R00 * R00 * state["beta_x"] - 2.0 * R00 * R01 * state["alpha_x"] + R01 * R01 * gamma_x
    beta_y = R22 * R22 * state["beta_y"] - 2.0 * R22 * R23 * state["alpha_y"] + R23 * R23 * gamma_y

    alpha_x = -R00 * R10 * state["beta_x"] + (R01 * R10 + R11 * R00) * state["alpha_x"] - R01 * R11 * gamma_x
    alpha_y = -R22 * R32 * state["beta_y"] + (R23 * R32 + R33 * R22) * state["alpha_y"] - R23 * R33 * gamma_y

    Dx = R00 * state["Dx"] + R01 * state["Dxp"] + R05
    Dxp = R10 * state["Dx"] + R11 * state["Dxp"] + R15
    Dy = R22 * state["Dy"] + R23 * state["Dyp"] + R25
    Dyp = R32 * state["Dy"] + R33 * state["Dyp"] + R35

    d_mux = _phase_advance_increment(R01, R00 * state["beta_x"] - R01 * state["alpha_x"], xp=xp)
    d_muy = _phase_advance_increment(R23, R22 * state["beta_y"] - R23 * state["alpha_y"], xp=xp)

    return {
        "beta_x": beta_x,
        "alpha_x": alpha_x,
        "beta_y": beta_y,
        "alpha_y": alpha_y,
        "Dx": Dx,
        "Dxp": Dxp,
        "Dy": Dy,
        "Dyp": Dyp,
        "mux": state["mux"] + d_mux,
        "muy": state["muy"] + d_muy,
        "E": state["E"] + xp.asarray(delta_e),
        "s": state["s"] + xp.asarray(length),
        "p": state["p"],
    }


def periodic_twiss_from_r(R, energy=0.0, xp=np):
    """
    Compute periodic linear Twiss parameters from a one-turn matrix.

    The JAX-oriented path assumes the caller evaluates this in a stable region
    where the periodic solution exists.
    """
    R55 = R[5, 5]
    if xp is np:
        scale = 1.0
        if R55 != 1:
            if energy == 0:
                raise TypeError("Lattice contains a cavity; periodic Twiss requires non-zero reference energy.")
            g0 = energy / m_e_GeV
            g1 = np.sqrt(g0 ** 2 - 1.0 + R55 ** 2) / R55
            scale = np.sqrt(g1 / g0)
    else:
        energy = xp.asarray(energy)
        R55 = xp.asarray(R55)
        has_energy_change = R55 != 1
        scale = xp.ones_like(R55)
        if energy is not None:
            g0 = energy / m_e_GeV
            g1 = xp.sqrt(g0 * g0 - 1.0 + R55 * R55) / R55
            scale = xp.where(has_energy_change, xp.sqrt(g1 / g0), scale)

    R00 = R[0, 0] * scale
    R01 = R[0, 1] * scale
    R10 = R[1, 0] * scale
    R11 = R[1, 1] * scale
    R22 = R[2, 2] * scale
    R23 = R[2, 3] * scale
    R32 = R[3, 2] * scale
    R33 = R[3, 3] * scale
    R05 = R[0, 5]
    R15 = R[1, 5]
    R25 = R[2, 5]
    R35 = R[3, 5]

    cosmx = (R00 + R11) / 2.0
    cosmy = (R22 + R33) / 2.0
    stable_x = xp.abs(cosmx) < 1.0
    stable_y = xp.abs(cosmy) < 1.0

    if xp is np and (not stable_x or not stable_y):
        return {
            "stable_x": False,
            "stable_y": False,
            "beta_x": np.nan,
            "beta_y": np.nan,
            "alpha_x": np.nan,
            "alpha_y": np.nan,
            "Dx": np.nan,
            "Dxp": np.nan,
            "Dy": np.nan,
            "Dyp": np.nan,
            "mux": 0.0,
            "muy": 0.0,
            "E": energy,
            "s": 0.0,
            "p": 0.0,
        }

    sinmx = xp.sign(R01) * xp.sqrt(xp.maximum(1.0 - cosmx * cosmx, 0.0))
    sinmy = xp.sign(R23) * xp.sqrt(xp.maximum(1.0 - cosmy * cosmy, 0.0))

    beta_x = xp.abs(R01 / sinmx)
    beta_y = xp.abs(R23 / sinmy)
    alpha_x = (R00 - R11) / (2.0 * sinmx)
    alpha_y = (R22 - R33) / (2.0 * sinmy)

    one = xp.asarray(1.0)
    Hx = xp.stack([
        xp.stack([R00 - one, R01]),
        xp.stack([R10, R11 - one]),
    ])
    Hhx = xp.stack([R05, R15])
    disp_x = xp.linalg.solve(-Hx, Hhx)

    Hy = xp.stack([
        xp.stack([R22 - one, R23]),
        xp.stack([R32, R33 - one]),
    ])
    Hhy = xp.stack([R25, R35])
    disp_y = xp.linalg.solve(-Hy, Hhy)

    state = make_twiss_state(
        beta_x=beta_x,
        alpha_x=alpha_x,
        beta_y=beta_y,
        alpha_y=alpha_y,
        Dx=disp_x[0],
        Dxp=disp_x[1],
        Dy=disp_y[0],
        Dyp=disp_y[1],
        E=energy,
        xp=xp,
    )
    state["stable_x"] = stable_x
    state["stable_y"] = stable_y
    return state


def _element_linear_overrides(overrides, element):
    if overrides is None:
        return {}
    if element in overrides:
        return overrides[element]
    elem_id = getattr(element, "id", None)
    if elem_id is not None and elem_id in overrides:
        return overrides[elem_id]
    return {}


def trace_linear_twiss(sequence_or_lattice, twiss0, xp=np, include_start=True, overrides=None):
    """
    Trace linear Twiss state at element boundaries using ``elem.linear_r_full``.

    Parameters
    ----------
    sequence_or_lattice
        Either a ``MagneticLattice`` or an iterable of elements.
    twiss0
        Initial state mapping created by ``make_twiss_state``.
    overrides
        Optional mapping from element instances or element ids to keyword
        overrides passed into ``elem.linear_r_full(...)``.
    """
    sequence = getattr(sequence_or_lattice, "sequence", sequence_or_lattice)
    state = _normalize_twiss_state(twiss0, xp=xp)
    states = [state] if include_start else []

    for elem in sequence:
        elem_overrides = _element_linear_overrides(overrides, elem)
        blocks = elem.linear_r_blocks(energy=state["E"], xp=xp, **elem_overrides)
        delta_e_blocks = elem.linear_delta_e_blocks(**elem_overrides)
        length_blocks = elem.linear_length_blocks(**elem_overrides)

        for R, delta_e, length in zip(blocks, delta_e_blocks, length_blocks):
            state = twiss_step_from_r(R, state, delta_e=delta_e, length=length, xp=xp)
        states.append(state)

    return states
