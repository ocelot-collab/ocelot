from __future__ import annotations

__author__ = 'Sergey Tomin'

from typing import Any, Callable, Dict, List, Tuple, Type, Union

import numpy as np
from scipy.integrate import simpson

from ocelot.common.globals import Z0, m_e_GeV, speed_of_light
from ocelot.cpbd.beam import SliceParameters, Twiss, s_to_cur


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

    b = np.abs(simpson(B[:, 1] / speed_of_light * np.exp(-1j * 2 * np.pi / lambda_mod * B[:, 0]), x=B[:, 0])) / np.sum(
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

    b = np.abs(simpson(B[:, 1] / speed_of_light * np.exp(-1j * 2 * np.pi / lambda_mod * B[:, 0]), x=B[:, 0])) / charge
    return b


def calculate_BMAG(tws_des: Twiss, tws_err: Twiss) -> tuple[Union[float, Any], Union[float, Any]]:
    """
    Function calculates mismatch and mismatch phase using two twiss lists.
    Result is saved in tws_err list as M_x, M_y, psi_x, psi_y

    :param tws_des: Twiss, design twiss parameters
    :param tws_err: Twiss, error twiss parameters
    :return: (Mx, My, phi_x, phi_y) mismatch at the end of lattice
    """
    gamma_x_des = (1 + tws_des.alpha_x * tws_des.alpha_x) / tws_des.beta_x
    gamma_y_des = (1 + tws_des.alpha_y * tws_des.alpha_y) / tws_des.beta_y
    gamma_x_err = (1 + tws_err.alpha_x * tws_err.alpha_x) / tws_err.beta_x
    gamma_y_err = (1 + tws_err.alpha_y * tws_err.alpha_y) / tws_err.beta_y
    mp_x = 0.5 * (tws_err.beta_x * gamma_x_des - 2 * tws_err.alpha_x * tws_des.alpha_x + tws_des.beta_x * gamma_x_err)
    mp_y = 0.5 * (tws_err.beta_y * gamma_y_des - 2 * tws_err.alpha_y * tws_des.alpha_y + tws_des.beta_y * gamma_y_err)
    lam_x = mp_x + np.sqrt(mp_x * mp_x - 1)
    lam_y = mp_y + np.sqrt(mp_y * mp_y - 1)
    return lam_x, lam_y


def calculate_slice_BMAG(slice_params: SliceParameters, tws_des: Twiss) -> tuple[Union[float, Any], Union[float, Any]]:
    """
    Function calculates mismatch and mismatch phase using two twiss lists.
    Result is saved in tws_err list as M_x, M_y, psi_x, psi_y

    :param tws_des: Twiss, design twiss parameters
    :param tws_err: Twiss, error twiss parameters
    :return: (Mx, My, phi_x, phi_y) mismatch at the end of lattice
    """
    gamma_x_des = (1 + tws_des.alpha_x * tws_des.alpha_x) / tws_des.beta_x
    gamma_y_des = (1 + tws_des.alpha_y * tws_des.alpha_y) / tws_des.beta_y

    mp_x = 0.5 * (slice_params.beta_x * gamma_x_des - 2 * slice_params.alpha_x * tws_des.alpha_x + tws_des.beta_x * slice_params.gamma_x)
    mp_y = 0.5 * (slice_params.beta_y * gamma_y_des - 2 * slice_params.alpha_y * tws_des.alpha_y + tws_des.beta_y * slice_params.gamma_y)
    lam_x = mp_x + np.sqrt(mp_x * mp_x - 1)
    lam_y = mp_y + np.sqrt(mp_y * mp_y - 1)
    return lam_x, lam_y


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


def single_plane_dipole_wake(p=0.5e-3, t=0.25e-3, b=500e-6, l=5):
    """
    Function calculates dipole (monopole) wake for single plane corrugated structure. Default parameters are taken
    for EuXFEL diagnostic streaker https://www.slac.stanford.edu/pubs/slacpubs/16750/slac-pub-16881.pdf
    Coefficient for s0yd is corrected by I.Zagorodnov 8/9 -> 1/2

    :param l: 5, corrugated plane length in m
    :param p: 0.5e-3  # period in m
    :param t: 0.25e-3  # Longitudinal gap im m
    :param b: 250e-6  # Distance from the plate im m
    :return: wake in V/Q
    """

    alpha = 1 - 0.465 * np.sqrt(t / p) - 0.070 * (t / p)
    s0yd = b ** 2 * t / (2 * np.pi * alpha ** 2 * p ** 2)
    wyd = lambda s: l * 2. / b ** 3 * s0yd * (
                1 - (1 + np.sqrt(s / s0yd)) * np.exp(-np.sqrt(s / s0yd))) * Z0 * speed_of_light / (4 * np.pi)

    return wyd


def single_plate_quadrupole_wake(p=0.5e-3, t=0.25e-3, b=500e-6, l=5):
    """
    Function calculates quadrupole wake for single plane corrugated structure. Default parameters are taken for EuXFEL
    diagnostic streaker https://www.slac.stanford.edu/pubs/slacpubs/16750/slac-pub-16881.pdf
    Coefficient for s0yq is corrected by I.Zagorodnov 8/9 -> 1/2

    :param l: corrugated plane length in m
    :param p: 0.5e-3  # period
    :param t: 0.25e-3  # Longitudinal gap
    :param b: 250e-6  # Distance from the plate
    :return: wake in V/Q
    """

    alpha = 1 - 0.465 * np.sqrt(t / p) - 0.070 * (t / p)
    s0yq = b ** 2 * t / (2 * np.pi * alpha ** 2 * p ** 2)
    wyq = lambda s: l * 3. / b ** 4 * s0yq * (
                1 - (1 + np.sqrt(s / s0yq)) * np.exp(-np.sqrt(s / s0yq))) * Z0 * speed_of_light / (4 * np.pi)

    return wyq


def convolve_beam(current, wake):
    """
    Function to convolve wake with beam current

    :param current: current[:, 0] - s in [m], current[:, 1] - current in [A]. The beam head is on the left
    :param wake: wake function in form: wake(s)
    :return: wake_kick[:, 0] - s in [m], wake_kick[:, 1] - V
    """
    s_shift = current[0, 0]
    current[:, 0] -= s_shift
    s = current[:, 0]
    step = (s[-1] - s[0]) / (len(s) - 1)

    q = current[:, 1] / speed_of_light

    w = np.array([wake(si) for si in s]).flatten()

    wake = np.convolve(q, w) * step
    s_new = (np.cumsum(np.ones(len(wake))) - 1.) * step
    wake_kick = np.vstack((s_new, wake))
    return wake_kick.T


def passive_streaker_resolutions(dipole_kick, quad_kick, R, tw, kick="vert", emittn_x=1e-6, emittn_y=1e-6, energy=14, sigma_R=30e-6):
    """
    Function to calculate time and energy resolution
    Example to use:
    I = self.get_current(num=50)
    R = self.R_matrix
    distance = 500e-6  # m
    wyd = recon.chirper_dipole_wake(p=0.5e-3, t=0.25e-3, b=distance, l=6)
    wyq = recon.chirper_quadrupole_wake(p=0.5e-3, t=0.25e-3, b=distance)

    quad_kick = recon.convolve_beam(I, wyq)
    dipole_kick = recon.convolve_beam(I, wyd)
    tw = self.tws_chirper
    r_temp, r_energy, sigma_x2, sigma_y2 = recon.calculate_resolutions(dipole_kick, quad_kick, R, tw, sigma_R, chirper_len,
                                                                   energy, emitt_x, emitt_y)

    :param kick:
    :param dipole_kick: dipole wake kick, example dipole_kick = convolve_beam(current, single_plane_dipole_wake(...))
    :param quad_kick: quad wake kick, example quad_kick = convolve_beam(current, single_plane_quad_wake(...))
    :param R: R-Matrix lattice between a passive streaker and a screen
    :param tw: Twiss parameters at a passive streaker position
    :param emittn_x: normalized emittance in horizontal plane
    :param emittn_y: normalized emittance in vertical plane
    :param energy: beam energy in GeV
    :param sigma_R: Resolution of the screen
    :return:

    Examples
    --------
    from ocelot import *
    from ocelot.gui import *
    from ocelot.utils.acc_utils import *
    from lattice import sase2
    from lattice import t3_bump_fin as t3


    lat = MagneticLattice(sase2.cell + t3.cell, stop=t3.otrb_2560_t3)
    B, R, T = lat.transfer_maps(energy=14, start=t3.ws_center, stop=t3.otrb_2560_t3)

    tws = twiss(lat, sase2.tws, attach2elem=True)

    distance = 500e-6

    tw = t3.ws_center.tws
    parray = generate_parray(sigma_x=1e-4, sigma_px=2e-5, sigma_tau=1e-3/200, sigma_p=1e-4, chirp=0.00, charge=250e-12,
                    nparticles=200000, energy=14, tws=tw, shape="gauss")

    I = parray.I()

    wyd = single_plane_dipole_wake(p=0.5e-3, t=0.25e-3, b=distance, l=5)
    wyq = single_plate_quadrupole_wake(p=0.5e-3, t=0.25e-3, b=distance, l=5)

    quad_kick = convolve_beam(I, wyq)
    dipole_kick = convolve_beam(I, wyd)

    r_temp, r_energy, sigma_x2, sigma_y2 = calculate_resolutions(dipole_kick, quad_kick, R, tw, kick="vert", emittn_x=0.6e-6,
                                                                 emittn_y=0.6e-6, energy=14, sigma_R=30e-6)

    fig1, ax1 = plt.subplots()
    fig1.suptitle("time resolution")
    color = 'tab:red'
    ax1.set_xlabel('s [mm]')
    ax1.set_ylabel('R(s) [fs]', color=color)
    ax1.plot(r_temp[:, 0]* 1000, r_temp[:, 1]/speed_of_light*1e15, "r-")
    ax1.set_ylim([0, 50])
    ax1.set_xlim([I[0, 0] * 1000, I[-1, 0] * 1000])

    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('I [kA]', color=color)  # we already handled the x-label with ax1
    ax2.plot(I[:, 0] * 1000, I[:, 1] / 1000, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax1.set_xlim([I[0, 0] * 1000, I[-1, 0] * 1000])


    fig2, ax1 = plt.subplots()
    fig2.suptitle("energy resolution")
    color = 'tab:red'
    ax1.set_xlabel('s [mm]')
    ax1.set_ylabel('R(s) [MV]', color=color)
    ax1.plot(r_energy[:, 0]* 1000, r_energy[:, 1]*1e-6, "r-")
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_ylim([0, 5])
    ax1.set_xlim([I[0, 0] * 1000, I[-1, 0] * 1000])

    ax2 = ax1.twinx()
    color = 'tab:blue'
    ax2.set_ylabel('I [kA]', color=color)  # we already handled the x-label with ax1
    ax2.plot(I[:, 0] * 1000, I[:, 1] / 1000, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    ax1.set_xlim([I[0, 0] * 1000, I[-1, 0] * 1000])

    plt.show()

    """
    gamma = energy / m_e_GeV
    emitt_x = emittn_x/gamma
    emitt_y = emittn_y/gamma
    energy_eV = energy * 1e9
    wq = quad_kick[:, 1] / energy_eV
    wd = dipole_kick[:, 1] / energy_eV
    ds = dipole_kick[1, 0] - dipole_kick[0, 0]
    wdp = np.gradient(wd, ds)

    if kick == "vert":
        sigma_y2 = emitt_y * (
                R[2, 3] ** 2 * wq ** 2 * tw.beta_y + 2 * R[2, 3] * wq * (R[2, 2] * tw.beta_y - R[2, 3] * tw.alpha_y) +
                (R[2, 2] * tw.beta_y - R[2, 3] * tw.alpha_y) ** 2 / tw.beta_y + R[2, 3] ** 2 / tw.beta_y)

        dy_ds = R[2, 3] * wdp
        r_temp = np.sqrt(sigma_R ** 2 + sigma_y2) / np.abs(dy_ds)

        # ENERGY RESOLUTION
        sigma_x2 = emitt_x * ((R[0, 0] - R[0, 1] * wq) ** 2 * tw.beta_x -
                          2 * R[0, 1] * (R[0, 0] - R[0, 1] * wq) * tw.alpha_y + R[0, 1] ** 2 * tw.gamma_x)

        r_energy = energy_eV / R[0, 5] * np.sqrt(sigma_R ** 2 + np.abs(sigma_x2))
    else:
        sigma_x2 = emitt_x * (R[0, 1] ** 2 * wq ** 2 * tw.beta_x + 2 * R[0, 1] * wq * (
                    R[0, 0] * tw.beta_x - R[0, 1] * tw.alpha_x) +
                              (R[0, 0] * tw.beta_x - R[0, 1] * tw.alpha_x) ** 2 / tw.beta_x + R[0, 1] ** 2 / tw.beta_x)

        dy_ds = R[0, 1] * wdp
        r_temp = np.sqrt(sigma_R ** 2 + sigma_x2) / np.abs(dy_ds)

        # ENERGY RESOLUTION
        sigma_y2 = emitt_y * ((R[2, 2] - R[2, 3] * wq) ** 2 * tw.beta_x -
                          2 * R[2, 3] * (R[2, 2] - R[2, 3] * wq) * tw.alpha_y + R[2, 3] ** 2 * tw.gamma_x)

        r_energy = energy_eV / R[0, 5] * np.sqrt(sigma_R ** 2 + np.abs(sigma_y2))

    r_temp = np.vstack([quad_kick[:, 0], r_temp]).T
    r_energy = np.vstack([quad_kick[:, 0], r_energy]).T
    return r_temp, r_energy, sigma_x2, sigma_y2
