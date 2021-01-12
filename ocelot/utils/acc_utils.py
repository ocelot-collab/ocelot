__author__ = 'Sergey Tomin'

import numpy as np
from scipy.integrate import simps
from ocelot.common.globals import speed_of_light
from ocelot.cpbd.beam import s_to_cur


def RTU_56(LB, LD, r, m):
    LB2 = LB * LB
    r2 = r * r
    K = np.sqrt(1 - LB2 / r2)
    r56 = m * r * np.arcsin(LB / r) - m * LB / K - 2 * LB2 * LD / (K ** 3 * r2) if r != np.inf else m * LB - m * LB / K

    t566 = (r2 * LB2 * (6 * LD + m * LB) - m * LB ** 5) / (2 * K * (LB2 - r2) ** 2)
    p1 = m * LB ** 7 - LB ** 4 * (5 * m * LB - 6 * LD) * r2 + 4 * LB2 * (6 * LD + m * LB) * r ** 4
    p2 = 6 * K * (LB2 - r2) ** 3
    u5666 = p1 / p2
    Sref = m * r * np.arcsin(LB / r) + 2 * LD / np.cos(np.arcsin(LB / r))
    return r56, t566, u5666, Sref


def chicane_RTU(yoke_len, dip_dist, r, type):
    """
    Method calculate R56, T566, U5666 and path length of the reference particle for chicanes 's' and 'c' type

    :param yoke_len: dipole yoke length
    :param dip_dist: distance between 1st and 2nd dipoles on Z-axis (not a particle path but projection on Z-axis)
    :param r: radii of the dipoles
    :param type: type of the chicane "s" or "c"
    :return: R56, T566, U5666, Sref (distance between magnet 2 and 3 is 0)
    """
    dx = lambda r, LB: np.sqrt(1 - LB ** 2 / r ** 2) * r * (r - np.sqrt(r ** 2 - LB ** 2)) / LB
    r56 = 0
    t566 = 0
    u5666 = 0
    Sref = 0

    if type == 'c':
        r56, t566, u5666, Sref = RTU_56(yoke_len, dip_dist, r, 4)

    elif type == 's':
        r56, t566, u5666, Sref = RTU_56(yoke_len, 2 * dip_dist + dx(r, yoke_len), r, 6)
    else:
        print("Unknown chicane type. Use 'c' or 's'. ")
    return r56, t566, u5666, Sref


def bunching(p_array, lambda_mod, smooth_sigma=None):
    """
    Function calculates bunching factor for wavelength lambda_mod

    $b(\lambda) = \frac{1}{N_0}\left| \langle e^{- i\frac{ 2 \pi}{\lambda} s} N(s)\rangle \right|$

    :param p_array: ParticleArray
    :param lambda_mod: wavelength
    :param smooth_sigma: smoothing parameter
    :return: bunching factor
    """
    if smooth_sigma is None:
        smooth_sigma = min(np.std(p_array.tau())*0.01, lambda_mod*0.1)

    B = s_to_cur(p_array.tau(), sigma=smooth_sigma, q0=np.sum(p_array.q_array), v=speed_of_light)

    b = np.abs(simps(B[:, 1] / speed_of_light * np.exp(-1j * 2 * np.pi / lambda_mod * B[:, 0]), B[:, 0])) / np.sum(
        p_array.q_array)
    return b


def slice_bunching(tau, charge, lambda_mod, smooth_sigma=None):
    """
    Function calculates bunching factor for wavelength lambda_mod

    $b(\lambda) = \frac{1}{N_0}\left| \langle e^{- i\frac{ 2 \pi}{\lambda} s} N(s)\rangle \right|$

    :param p_array: ParticleArray
    :param lambda_mod: wavelength
    :param smooth_sigma: smoothing parameter
    :return: bunching factor
    """
    if smooth_sigma is None:
        smooth_sigma = lambda_mod/10.

    B = s_to_cur(tau, sigma=smooth_sigma, q0=charge, v=speed_of_light)

    b = np.abs(simps(B[:, 1] / speed_of_light * np.exp(-1j * 2 * np.pi / lambda_mod * B[:, 0]), B[:, 0])) / charge
    return b


def calculate_mismatch(tws_des, tws_err):
    """
    Function calculates mismatch and mismatch phase using two twiss lists.
    Result is saved in tws_err list as M_x, M_y, psi_x, psi_y

    :param tws_des: list, design twiss parameters
    :param tws_err: list, error twiss parameters
    :return: (Mx, My, phi_x, phi_y) mismatch at the end of lattice
    """
    if len(tws_des) != len(tws_err):
        raise Exception("length of the twiss lists must be the same")

    Sigma = lambda tws: np.array([[tws.beta_x, -tws.alpha_x, 0, 0],
                                  [-tws.alpha_x, tws.gamma_x, 0, 0],
                                  [0, 0, tws.beta_y, -tws.alpha_y],
                                  [0, 0, -tws.alpha_y, tws.gamma_y]])

    S = lambda tws: np.array([[1 / np.sqrt(tws.beta_x), 0, 0, 0],
                              [tws.alpha_x / np.sqrt(tws.beta_x), np.sqrt(tws.beta_x), 0, 0],
                              [0, 0, 1 / np.sqrt(tws.beta_y), 0],
                              [0, 0, tws.alpha_y / np.sqrt(tws.beta_y), np.sqrt(tws.beta_y)]])

    for i in range(len(tws_des)):
        tws_e = tws_err[i]
        tws_m = tws_des[i]

        s = S(tws_m)
        Sigma_b = np.dot(np.dot(s, Sigma(tws_e)), s.T)
        bx, gx, ax = Sigma_b[0,0], Sigma_b[1,1], -Sigma_b[0,1]
        by, gy, ay = Sigma_b[2, 2], Sigma_b[3, 3], -Sigma_b[2, 3]
        if bx + gx < 2:
            tws_err[i].M_x = 1
        else:
            tws_err[i].M_x = 0.5*(bx + gx + np.sqrt((bx + gx)**2 - 4))


        if by + gy < 2:
            tws_err[i].M_y = 1
        else:
            tws_err[i].M_y = 0.5*(by + gy + np.sqrt((by + gy)**2 - 4))

        theta = np.arctan2(-2*ax, bx-gx)/2

        if tws_err[i].M_x > 1.05:
            tws_err[i].psi_x = (theta + tws_m.mux)%(np.pi)*180/np.pi - 90
        else:
            tws_err[i].psi_x = 0

        thetay = np.arctan2(-2*ay, by-gy)/2
        if tws_err[i].M_y > 1.05:
            tws_err[i].psi_y = (thetay + tws_m.muy)%(np.pi)*180/np.pi - 90
        else:
            tws_err[i].psi_y = (0)

    return tws_err[-1].M_x, tws_err[-1].M_y, tws_err[-1].psi_x, tws_err[-1].psi_y


def rf2beam(v1, phi1, vh, phih, n=3, freq=1.3e9, E0=0.00675, zeta1=0., zeta2=0., zeta3=0.):
    """
    Function calculates beam parameters: the final beam energy, chirp, curvature, skewness,
    from the RF parameters: voltages and phases.
    Note: in EuXFEL case, L1: E0 = 130 MeV, to get correct Sum Voltage for L1/L2, E1' = E1 - E0

    :param v1: voltage [GeV] of the first harmonic cavity
    :param phi1: phase [deg] of the first harmonic cavity
    :param vh: voltage [GeV] of the high harmonic cavity
    :param phih: phase [deg] of the high harmonic cavity
    :param n: 3, number of harmonic of the high harmonic cavity. If n = 0 the high harmonic cavity does not exist
    :param freq: frequency [Hz] of the first harmonic cavity
    :param E0: initial beam energy [GeV] (from the gun)
    :param zeta1: initial beam chirp (from the gun)
    :param zeta2: initial beam curvature (from the gun)
    :param zeta3: initial beam skewness (from the gun)
    :return: E1, chirp, curvature, skewness
    """
    k = 2 * np.pi * freq / speed_of_light
    deg2rad = np.pi / 180

    M = np.array([[1, 0, 1, 0],
                  [0, -k, 0, -(n * k)],
                  [-k ** 2, 0, -(n * k) ** 2, 0],
                  [0, k ** 3, 0, (n * k) ** 3]])

    V = np.array([v1 * np.cos(phi1 * deg2rad),
                  v1 * np.sin(phi1 * deg2rad),
                  vh * np.cos(phih * deg2rad),
                  vh * np.sin(phih * deg2rad)]).T

    if n == 0:
        M = np.array([[1, 0],
                      [0, -k]])

        V = np.array([v1 * np.cos(phi1 * deg2rad),
                      v1 * np.sin(phi1 * deg2rad)]).T
    R = np.dot(M, V)
    E1 = E0 + R[0]
    chirp = (R[1] + E0 * zeta1) / E1

    if n == 0:
        return E1, chirp, 0, 0

    curvature = (R[2] + E0 * zeta2 / 2) / E1
    skewness = (R[3] + E0 * zeta3 / 6) / E1

    return E1, chirp, curvature, skewness


def beam2rf(E1, chirp, curvature, skewness, n, freq, E0=0.00675, zeta1=0., zeta2=0., zeta3=0.):
    """
    Function calculates RF parameters: cavities (first and high harmonic) voltage [GV] and phase [deg]
    from the final beam energy, chirp, curvature, skewness.

    :param E1: the beam energy [GeV] after RF system
    :param chirp: the beam chirp
    :param curvature: the beam curvature
    :param skewness: the beam skewness
    :param n: 3, number of harmonic of the high harmonic cavity.  If n = 0 the high harmonic cavity does not exist.
    :param freq: frequency [Hz] of the first harmonic cavity
    :param E0: initial beam energy [GeV] (from the gun)
    :param zeta1: initial beam chirp (from the gun)
    :param zeta2: initial beam curvature (from the gun)
    :param zeta3: initial beam skewness (from the gun)
    :return: v1, phi1, vh, phih
    """

    k = 2 * np.pi * freq / speed_of_light
    M = np.array([[1, 0, 1, 0],
                  [0, -k, 0, -(n * k)],
                  [-k ** 2, 0, -(n * k) ** 2, 0],
                  [0, k ** 3, 0, (n * k) ** 3]])

    r = np.array([E1 - E0, chirp * E1 - E0 * zeta1, curvature * E1 - E0 * zeta2 / 2, skewness * E1 - E0 * zeta3 / 6])

    if n == 0:
        M = np.array([[1, 0],
                      [0, -k]])
        r = np.array(
            [E1 - E0, chirp * E1 - E0 * zeta1])
    rf = np.dot(np.linalg.inv(M), r)
    X1 = rf[0]
    Y1 = rf[1]
    rad2deg = 180 / np.pi
    v1 = np.sqrt(X1 ** 2 + Y1 ** 2)
    phi1 = (np.arctan(Y1 / X1) + np.pi / 2 * (1 - np.sign(X1))) * rad2deg

    if n == 0:
        return v1, phi1, 0, 0

    X13 = rf[2]
    Y13 = rf[3]
    vh = np.sqrt(X13 ** 2 + Y13 ** 2)
    phih = (np.arctan(Y13 / X13) + np.pi / 2 * (1 - np.sign(X13)) - 2 * np.pi) * rad2deg
    return v1, phi1, vh, phih


def beam2rf_xfel_linac(sum_voltage, chirp, init_energy=0.13):
    """
    wrapped up function for EuXFEL linacs

    :param sum_voltage: in control system [GeV]
    :param chirp: in control system [GeV]
    :param init_energy: for L1 it is 0.13 GeV, L2 = 0.7 GeV
    :return: v1, phi1
    """
    v1, phi1, _, _ = beam2rf(sum_voltage + init_energy, chirp, curvature=0, skewness=0, n=0, freq=1.3e9,
                             E0=init_energy, zeta1=0., zeta2=0., zeta3=0.)
    return v1, phi1


def rf2beam_xfel_linac(v, phi, init_energy=0.13):
    """
    wrapped up function for EuXFEL linacs

    :param v:
    :param phi:
    :param init_energy: for L1 it is 0.13 GeV, L2 = 0.7 GeV
    :return:
    """
    E1, chirp, _, _ = rf2beam(v, phi, vh=0, phih=0, n=0, freq=1.3e9, E0=init_energy, zeta1=0., zeta2=0., zeta3=0.)
    sum_voltage = E1 - init_energy
    return sum_voltage, chirp

