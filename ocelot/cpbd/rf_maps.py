"""Longitudinal RF maps shared by standing- and travelling-wave cavities.

The coordinates are tau = c * (t - t_ref) and delta = (E_particle - E) / (p*c).
The entrance RF phase is phi - k*tau, with k = 2*pi*freq/c. Both cavity models
hold this phase fixed while ramping each particle's energy linearly through
the body. This is an effective cavity model, not integration through an RF
field with evolving phase slip.

The time-of-flight derivatives below are evaluated without differences of
nearly equal energies. They remain finite at zero voltage and at either RF
zero crossing. See demos/ipython_tutorials/R_matrices.ipynb for the derivation.

Units: V [GV], E [GeV], z [m], freq [Hz], phi [rad].
"""

import numpy as np

from ocelot.common.globals import m_e_GeV, speed_of_light


def _tof_derivatives(g0, g1):
    """Momenta and derivatives of F = integral_0^1 du / beta(g0+(g1-g0)*u).

    F = (g0+g1)/(p0+p1), where p_i = sqrt(g_i**2-1). Write
    C = (1+g0*g1-p0*p1)/(p0+p1)**2, rationalizing the numerator to avoid
    relativistic cancellation. Then F_i = -C/p_i and
    F_ij = C*F/(p_i*p_j) + Kronecker_ij*C*g_i/p_i**3.
    Both reference energies must be above the rest energy.
    """
    p0 = np.sqrt((g0 - 1.) * (g0 + 1.))
    p1 = np.sqrt((g1 - 1.) * (g1 + 1.))
    ps = p0 + p1
    f = (g0 + g1) / ps
    c = (1. + (g0 * g0 + g1 * g1 - 1.) / (g0 * g1 + p0 * p1)) / ps ** 2
    f0, f1 = -c / p0, -c / p1
    f00 = c * (f / p0 ** 2 + g0 / p0 ** 3)
    f01 = c * f / (p0 * p1)
    f11 = c * (f / p1 ** 2 + g1 / p1 ** 3)
    return p0, p1, f0, f1, f00, f01, f11


def rf_long_r_terms(z, V, E, freq, phi):
    """Return (R55 - 1, R56, R65, R66) for a cavity body or slice.

    V and z are the slice voltage and length; E is its entrance total energy.
    phi is in radians. The RF phase depends on the time-like coordinate tau,
    without an additional beta factor.
    """
    if E == 0. and V == 0.:
        return 0., 0., 0., 1.
    g0 = E / m_e_GeV
    g1 = (E + V * np.cos(phi)) / m_e_GeV
    p0, p1, f0, f1, _, _, _ = _tof_derivatives(g0, g1)
    a = V / m_e_GeV * (2. * np.pi * freq / speed_of_light) * np.sin(phi)
    return z * a * f1, z * p0 * (f0 + f1), a / p1, p0 / p1


def rf_long_t_terms(z, V, E, freq, phi):
    """Return (T566, T556, T555) in tau1 += T566*delta0**2
    + T556*tau0*delta0 + T555*tau0**2.

    These are polynomial coefficients: T556 includes the full mixed term,
    and is twice the corresponding entry of a symmetric Taylor tensor.
    Arguments use the same units and slice convention as rf_long_r_terms.
    """
    if E == 0. and V == 0.:
        return 0., 0., 0.
    g0 = E / m_e_GeV
    g1 = (E + V * np.cos(phi)) / m_e_GeV
    p0, _, _, f1, f00, f01, f11 = _tof_derivatives(g0, g1)
    k = 2. * np.pi * freq / speed_of_light
    v = V / m_e_GeV
    a = v * k * np.sin(phi)
    T566 = 0.5 * z * p0 ** 2 * (f00 + 2. * f01 + f11)
    T556 = z * p0 * a * (f01 + f11)
    T555 = 0.5 * z * (a ** 2 * f11 - v * k ** 2 * np.cos(phi) * f1)
    return T566, T556, T555


def rf_longitudinal_kick(X, X4, X5, E, V, z, freq, phi):
    """Apply the nonlinear longitudinal map after the linear (R, B) map.

    X is modified in place. X4 and X5 are copies of entrance tau and delta,
    before application of the linear map. Other arguments follow rf_long_r_terms.
    """
    if E == 0. and V == 0.:
        return X
    E1 = E + V * np.cos(phi)
    pc0 = np.sqrt((E - m_e_GeV) * (E + m_e_GeV))
    pc1 = np.sqrt((E1 - m_e_GeV) * (E1 + m_e_GeV))
    phase_shift = X4 * (2. * np.pi * freq / speed_of_light)
    # Stable form of cos(phi - phase_shift) - cos(phi), including tiny tau.
    cos_difference = 2. * np.sin(phi - phase_shift / 2.) * np.sin(phase_shift / 2.)
    X[5] = X5 * pc0 / pc1 + V / pc1 * cos_difference
    T566, T556, T555 = rf_long_t_terms(z, V, E, freq, phi)
    X[4] += T566 * X5 * X5 + T556 * X4 * X5 + T555 * X4 * X4
    return X
