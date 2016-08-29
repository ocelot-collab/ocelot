# from numpy.core.umath import sqrt
from ocelot.common.globals import m_e_eV
from ocelot.cpbd.beam import *
import numpy as np


def exact_xp_2_xxstg_mad(xp, gamref):
    # to mad format
    N = xp.shape[0]
    xxstg = np.zeros((N, 6))
    pref = m_e_eV * np.sqrt(gamref ** 2 - 1)
    u = np.c_[xp[:, 3], xp[:, 4], xp[:, 5] + pref]
    gamma = np.sqrt(1 + np.sum(u * u, 1) / m_e_eV ** 2)
    beta = np.sqrt(1 - gamma ** -2)
    betaref = np.sqrt(1 - gamref ** -2)
    if np.__version__ > "1.8":
        p0 = np.linalg.norm(u, 2, 1).reshape((N, 1))
    else:
        p0 = np.sqrt(u[:, 0] ** 2 + u[:, 1] ** 2 + u[:, 2] ** 2).reshape((N, 1))
    u = u / p0
    cdt = -xp[:, 2] / (beta * u[:, 2])
    xxstg[:, 0] = xp[:, 0] + beta * u[:, 0] * cdt
    xxstg[:, 2] = xp[:, 1] + beta * u[:, 1] * cdt
    xxstg[:, 4] = cdt
    xxstg[:, 1] = xp[:, 3] / pref
    xxstg[:, 3] = xp[:, 4] / pref
    xxstg[:, 5] = (gamma / gamref - 1) / betaref
    return xxstg


def exact_xxstg_2_xp_mad(xxstg, gamref):
    # from mad format
    N = len(xxstg) / 6
    xp = np.zeros((N, 6))
    pref = m_e_eV * np.sqrt(gamref ** 2 - 1)
    betaref = np.sqrt(1 - gamref ** -2)
    gamma = (betaref * xxstg[5::6] + 1) * gamref
    beta = np.sqrt(1 - gamma ** -2)

    pz2pref = np.sqrt(((gamma * beta) / (gamref * betaref)) ** 2 - xxstg[1::6] ** 2 - xxstg[3::6] ** 2)

    u = np.c_[xxstg[1::6] / pz2pref, xxstg[3::6] / pz2pref, np.ones(N)]
    if np.__version__ > "1.8":
        norm = np.linalg.norm(u, 2, 1).reshape((N, 1))
    else:
        norm = np.sqrt(u[:, 0] ** 2 + u[:, 1] ** 2 + u[:, 2] ** 2).reshape((N, 1))
    u = u / norm
    xp[:, 0] = xxstg[0::6] - u[:, 0] * beta * xxstg[4::6]
    xp[:, 1] = xxstg[2::6] - u[:, 1] * beta * xxstg[4::6]
    xp[:, 2] = -u[:, 2] * beta * xxstg[4::6]
    xp[:, 3] = u[:, 0] * gamma * beta * m_e_eV
    xp[:, 4] = u[:, 1] * gamma * beta * m_e_eV
    xp[:, 5] = u[:, 2] * gamma * beta * m_e_eV - pref
    return xp


def exact_xp_2_xxstg_dp(xp, gamref):
    # dp/p0
    N = xp.shape[0]
    xxstg = np.zeros((N, 6))
    pref = m_e_eV * np.sqrt(gamref ** 2 - 1)
    u = np.c_[xp[:, 3], xp[:, 4], xp[:, 5] + pref]
    gamma = np.sqrt(1 + np.sum(u * u, 1) / m_e_eV ** 2)
    beta = np.sqrt(1 - gamma ** -2)
    if np.__version__ > "1.8":
        p0 = np.linalg.norm(u, 2, 1).reshape((N, 1))
    else:
        p0 = np.sqrt(u[:, 0] ** 2 + u[:, 1] ** 2 + u[:, 2] ** 2).reshape((N, 1))
    u = u / p0
    cdt = -xp[:, 2] / (beta * u[:, 2])
    xxstg[:, 0] = xp[:, 0] + beta * u[:, 0] * cdt
    xxstg[:, 2] = xp[:, 1] + beta * u[:, 1] * cdt
    xxstg[:, 4] = cdt
    xxstg[:, 1] = u[:, 0] / u[:, 2]
    xxstg[:, 3] = u[:, 1] / u[:, 2]
    xxstg[:, 5] = p0.reshape(N) / pref - 1
    return xxstg


def exact_xxstg_2_xp_dp(xxstg, gamref):
    # dp/p0
    N = len(xxstg) / 6
    xp = np.zeros((N, 6))
    pref = m_e_eV * np.sqrt(gamref ** 2 - 1)

    p = pref * (1 + xxstg[5::6])
    gamma = np.sqrt((p / m_e_eV) ** 2 + 1)

    beta = np.sqrt(1 - gamma ** -2)
    u = np.c_[xxstg[1::6], xxstg[3::6], np.ones(N)]
    if np.__version__ > "1.8":
        norm = np.linalg.norm(u, 2, 1).reshape((N, 1))
    else:
        norm = np.sqrt(u[:, 0] ** 2 + u[:, 1] ** 2 + u[:, 2] ** 2).reshape((N, 1))
    u = u / norm
    xp[:, 0] = xxstg[0::6] - u[:, 0] * beta * xxstg[4::6]
    xp[:, 1] = xxstg[2::6] - u[:, 1] * beta * xxstg[4::6]
    xp[:, 2] = -u[:, 2] * beta * xxstg[4::6]
    xp[:, 3] = u[:, 0] * gamma * beta * m_e_eV
    xp[:, 4] = u[:, 1] * gamma * beta * m_e_eV
    xp[:, 5] = u[:, 2] * gamma * beta * m_e_eV - pref
    return xp


def exact_xp_2_xxstg_de(xp, gamref):
    # dE/E0
    N = xp.shape[0]
    xxstg = np.zeros((N, 6))
    pref = m_e_eV * np.sqrt(gamref ** 2 - 1)
    u = np.c_[xp[:, 3], xp[:, 4], xp[:, 5] + pref]
    gamma = np.sqrt(1 + np.sum(u * u, 1) / m_e_eV ** 2)
    beta = np.sqrt(1 - gamma ** -2)
    if np.__version__ > "1.8":
        p0 = np.linalg.norm(u, 2, 1).reshape((N, 1))
    else:
        p0 = np.sqrt(u[:, 0] ** 2 + u[:, 1] ** 2 + u[:, 2] ** 2).reshape((N, 1))
    u = u / p0
    cdt = -xp[:, 2] / (beta * u[:, 2])
    xxstg[:, 0] = xp[:, 0] + beta * u[:, 0] * cdt
    xxstg[:, 2] = xp[:, 1] + beta * u[:, 1] * cdt
    xxstg[:, 4] = cdt
    xxstg[:, 1] = u[:, 0] / u[:, 2]
    xxstg[:, 3] = u[:, 1] / u[:, 2]
    xxstg[:, 5] = gamma / gamref - 1
    return xxstg


def exact_xxstg_2_xp_de(xxstg, gamref):
    # dE/E0
    N = len(xxstg) / 6
    xp = np.zeros((N, 6))
    pref = m_e_eV * np.sqrt(gamref ** 2 - 1)
    gamma = gamref * (1 + xxstg[5::6])
    beta = np.sqrt(1 - gamma ** -2)
    u = np.c_[xxstg[1::6], xxstg[3::6], np.ones(N)]
    if np.__version__ > "1.8":
        norm = np.linalg.norm(u, 2, 1).reshape((N, 1))
    else:
        norm = np.sqrt(u[:, 0] ** 2 + u[:, 1] ** 2 + u[:, 2] ** 2).reshape((N, 1))
    u = u / norm
    xp[:, 0] = xxstg[0::6] - u[:, 0] * beta * xxstg[4::6]
    xp[:, 1] = xxstg[2::6] - u[:, 1] * beta * xxstg[4::6]
    xp[:, 2] = -u[:, 2] * beta * xxstg[4::6]
    xp[:, 3] = u[:, 0] * gamma * beta * m_e_eV
    xp[:, 4] = u[:, 1] * gamma * beta * m_e_eV
    xp[:, 5] = u[:, 2] * gamma * beta * m_e_eV - pref
    return xp


# from pylab import *
def astraBeam2particleArray(filename):
    """
    function convert Astra beam distribution to Ocelot format - ParticleArray
    :type filename: str
    :return: ParticleArray
    """
    P0 = np.loadtxt(filename)
    charge_array = -P0[:, 7] * 1e-9  # charge in nC -> in C
    print("Astra to Ocelot: charge = ", sum(charge_array))
    print("Astra to Ocelot: particles number = ", len(charge_array))
    xp = P0[:, :6]
    Pref = xp[0, 5]
    s_ref = xp[0, 2]
    xp[0, 5] = 0
    xp[0, 2] = 0.
    #    print(xp[1:, 2])
    #    plot(xp[1:, 2], xp[1:, 5]/Pref, "b.")
    #    show()
    gamref = np.sqrt((Pref / m_e_eV) ** 2 + 1)
    xxstg = exact_xp_2_xxstg_mad(xp, gamref)

    p_array = ParticleArray(len(charge_array))
    p_array.s = s_ref
    p_array.E = np.sqrt((Pref / m_e_eV) ** 2 + 1) * m_e_GeV
    print("Astra to Ocelot: energy = ", p_array.E)
    print("Astra to Ocelot: s pos = ", p_array.s)
    p_array.particles[0::6] = xxstg[:, 0]
    p_array.particles[1::6] = xxstg[:, 1]
    p_array.particles[2::6] = xxstg[:, 2]
    p_array.particles[3::6] = xxstg[:, 3]
    p_array.particles[4::6] = xxstg[:, 4]
    p_array.particles[5::6] = xxstg[:, 5]
    p_array.q_array = charge_array
    return p_array


def particleArray2astraBeam(p_array, filename="tytest.ast"):
    """
    function convert  Ocelot's ParticleArray to Astra beam distribution and save to "filename".
    :param p_array:
    :param filename:
    :return:
    """
    gamref = p_array.E / m_e_GeV
    s0 = p_array.s
    P = p_array.particles.view()
    Np = len(P) / 6
    xp = exact_xxstg_2_xp_mad(P, gamref)
    Pref = np.sqrt(p_array.E ** 2 / m_e_GeV ** 2 - 1) * m_e_eV
    xp[:, 5] = xp[:, 5] + Pref
    xp[:, 2] = xp[:, 2] + s0
    xp[1:Np, 5] = xp[1:Np, 5] - xp[0, 5]
    xp[1:Np, 2] = xp[1:Np, 2] - xp[0, 2]

    charge_array = -p_array.q_array.reshape(len(p_array.q_array), 1) * 1e+9  # charge in C -> in nC
    flag = np.zeros((len(charge_array), 1))
    astra = np.append(xp, flag, axis=1)
    astra = np.append(astra, charge_array, axis=1)
    astra = np.append(astra, flag, axis=1)
    astra = np.append(astra, flag, axis=1)
    print("SAVE")
    np.savetxt(filename, astra, fmt='%.7e')