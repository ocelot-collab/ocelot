"""
Longitudinal RF physics shared by the standing-wave and travelling-wave cavities.

Both structures ramp the reference energy linearly along the body,
``gamma(z) = gamma0 + gamma' z``, and both damp the momentum deviation
adiabatically, ``delta1 = delta0 * p0 / p1``. The longitudinal map therefore
depends only on (z, V, E, freq, phi) and is identical for the two families -
only the transverse block differs (ponderomotive focusing for the standing
wave, edge lenses for the travelling wave). These helpers hold that common
part so ``CavityAtom``/``CavityTM`` and ``TWCavityAtom``/``TWCavityTM`` cannot
drift apart.

Units: V in [GV], E in [GeV], z in [m], freq in [Hz], phi in [rad].
"""

import numpy as np

from ocelot.common.globals import m_e_GeV, speed_of_light

# Relative |dgamma| below which the general expressions lose their accuracy to
# cancellation and the analytic small-energy-gain limits are used instead.
SMALL_ENERGY_GAIN = 1e-8


def rf_long_r_terms(z, V, E, freq, phi):
    """
    First-order longitudinal block of an RF structure.

    :param z: length of the (slice of the) structure [m]
    :param V: voltage of the (slice of the) structure [GV]
    :param E: reference energy at the entrance [GeV]
    :param freq: RF frequency [Hz]
    :param phi: RF phase [rad]
    :return: (r55_cor, r56, r65, r66), where R[4, 4] = 1 + r55_cor
    """
    de = V * np.cos(phi)
    g0 = E / m_e_GeV
    g1 = (E + de) / m_e_GeV

    beta0 = 1.
    beta1 = 1.
    r56 = 0.
    r55_cor = 0.
    k = 2. * np.pi * freq / speed_of_light

    if V != 0 and E != 0:
        beta0 = np.sqrt(1. - 1. / (g0 * g0))
        beta1 = np.sqrt(1. - 1. / (g1 * g1))

        r56 = -z / (g1 * g1 * g0 * beta1) * (g1 + g0) / (beta1 + beta0)

        # R55 differs from 1 because the RF kick couples back into the arrival
        # time through the velocity dispersion of the ramp.
        #
        # Two known limitations are kept EXACTLY as they are so that this
        # helper reproduces the standing-wave cavity bit for bit. Do not
        # "tidy" either one without deciding to change Cavity's results:
        #
        #  1. The zero-crossing limit below has no sin(phi) factor, so it
        #     returns the same value at phi = 90 and phi = 270 deg. R55 - 1 is
        #     odd in the phase, so the sign is wrong at 270 deg (by 2.1e-5 for
        #     1 m / 20 MV / 150 MeV). The correct limit is the expression
        #     below times sin(phi).
        #  2. The general branch is taken whenever |cos(phi)| >= 1e-3, even
        #     for a vanishing voltage, where (g0 g1 (b0 b1 - 1) + 1) vanishes
        #     as (g1 - g0)**2 and the ratio is a finite 0/0 that must not be
        #     evaluated numerically. At V = 1e-12 GV it returns -1.8e-5 where
        #     the limit is 0. Widening the guard to |dgamma| < 1e-8 |g0| fixes
        #     it.
        #
        # Neither is reachable in normal on-crest or off-crest operation.
        delta_gamma = g1 - g0
        zero_energy_gain = abs(delta_gamma) < SMALL_ENERGY_GAIN * abs(g0)
        if zero_energy_gain and abs(np.cos(phi)) < 1e-3:
            r55_cor = -k * z * V / (2. * m_e_GeV * g0 ** 3 * beta0 ** 2)
        else:
            r55_cor = (k * z * beta0 * V / m_e_GeV * np.sin(phi)
                       * (g0 * g1 * (beta0 * beta1 - 1.) + 1.)
                       / (beta1 * g1 * (g0 - g1) ** 2))

    r66 = g0 / g1 * beta0 / beta1
    r65 = k * np.sin(phi) * V / (g1 * beta1 * m_e_GeV)
    return r55_cor, r56, r65, r66


def rf_long_t_terms(z, V, E, freq, phi):
    """
    Second-order longitudinal coefficients of an RF structure.

    :param z: length of the (slice of the) structure [m]
    :param V: voltage of the (slice of the) structure [GV]
    :param E: reference energy at the entrance [GeV]
    :param freq: RF frequency [Hz]
    :param phi: RF phase [rad]
    :return: (T566, T556, T555)
    """
    beta0 = 1.
    igamma2 = 0.
    g0 = 1e10  # large gamma if E == 0
    if E != 0.:
        g0 = E / m_e_GeV
        igamma2 = 1. / (g0 * g0)
        beta0 = np.sqrt(1. - igamma2)

    # Default: drift-like second order (T566^drift, T555 = T556 = 0)
    T566 = 1.5 * z * igamma2 / (beta0 ** 3)  # = 3L/(2 beta^3 gamma^2)
    T556 = 0.
    T555 = 0.

    E1 = E + V * np.cos(phi)
    if E1 <= 0.:
        # Non-physical final energy: fall back to the drift map.
        return T566, T556, T555

    k = 2. * np.pi * freq / speed_of_light
    g1 = E1 / m_e_GeV
    beta1 = np.sqrt(1. - 1. / (g1 * g1))

    dgamma = V / m_e_GeV

    if abs(g1 - g0) < SMALL_ENERGY_GAIN * abs(g0):
        # --- Small dgamma limit (g1 ~ g0): T566 keeps its drift value ---
        if abs(np.cos(phi)) < 1e-3:
            # Zero-crossing cavity: analytic limits for T555 and T556
            T556 = 1.5 * z * k * dgamma / (beta0 ** 3 * g0 ** 3)
            T555 = 0.5 * z * k * k * dgamma * dgamma / (beta0 ** 3 * g0 ** 4)
        # Away from the zero crossing a vanishing gradient means the structure
        # behaves like a drift in second order, so T556 = T555 = 0.
    else:
        # --- General case: full RF formulas, valid for any sign of dE ---
        T566 = (z * (beta0 ** 3 * g0 ** 3 - beta1 ** 3 * g1 ** 3)
                / (2. * beta0 * beta1 ** 3 * g0 * (g0 - g1) * g1 ** 3))

        T556 = (beta0 * k * z * dgamma * g0 * (beta1 ** 3 * g1 ** 3 + beta0 * (g0 - g1 ** 3)) * np.sin(phi) /
                (beta1 ** 3 * g1 ** 3 * (g0 - g1) ** 2))

        T555 = (beta0 ** 2 * k ** 2 * z * dgamma / 2.
                * (dgamma * (2. * g0 * g1 ** 3 * (beta0 * beta1 ** 3 - 1.)
                             + g0 ** 2 + 3. * g1 ** 2 - 2.) / (beta1 ** 3 * g1 ** 3 * (g0 - g1) ** 3) * np.sin(phi) ** 2
                   - (g1 * g0 * (beta1 * beta0 - 1.) + 1.)
                   / (beta1 * g1 * (g0 - g1) ** 2)
                   * np.cos(phi)
                   )
                )

    return T566, T556, T555


def rf_longitudinal_kick(X, X4, X5, E, V, z, freq, phi):
    """
    Apply the RF longitudinal map to particles whose linear (R, B) part has
    already been applied.

    :param X: rparticles, modified in place
    :param X4: copy of X[4] (tau) before the linear map
    :param X5: copy of X[5] (delta) before the linear map
    :param E: reference energy at the entrance [GeV]
    :param V: voltage of the (slice of the) structure [GV]
    :param z: length of the (slice of the) structure [m]
    :param freq: RF frequency [Hz]
    :param phi: RF phase [rad]
    :return: X
    """
    beta0 = 1.
    if E != 0.:
        beta0 = np.sqrt(1. - (m_e_GeV / E) ** 2)

    E1 = E + V * np.cos(phi)
    T566, T556, T555 = rf_long_t_terms(z, V, E, freq, phi)

    if E1 > 0.:
        k = 2. * np.pi * freq / speed_of_light
        g1 = E1 / m_e_GeV
        beta1 = np.sqrt(1. - 1. / (g1 * g1))

        # delta1 = E0 b0 / (E1 b1) delta0 + V b0 / (E1 b1) [cos(phi - k b0 tau0) - cos phi]
        X[5] = (X5 * E * beta0 / (E1 * beta1)
                + V * beta0 / (E1 * beta1) * (np.cos(-X4 * beta0 * k + phi) - np.cos(phi)))

    X[4] += T566 * X5 * X5 + T556 * X4 * X5 + T555 * X4 * X4
    return X
