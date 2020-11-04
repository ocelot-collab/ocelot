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
    N = int(xxstg.size / 6)
    xp = np.zeros((N, 6))
    pref = m_e_eV * np.sqrt(gamref ** 2 - 1)
    betaref = np.sqrt(1 - gamref ** -2)
    gamma = (betaref * xxstg[5] + 1) * gamref
    beta = np.sqrt(1 - gamma ** -2)

    pz2pref = np.sqrt(((gamma * beta) / (gamref * betaref)) ** 2 - xxstg[1] ** 2 - xxstg[3] ** 2)

    u = np.c_[xxstg[1] / pz2pref, xxstg[3] / pz2pref, np.ones(N)]
    if np.__version__ > "1.8":
        norm = np.linalg.norm(u, 2, 1).reshape((N, 1))
    else:
        norm = np.sqrt(u[:, 0] ** 2 + u[:, 1] ** 2 + u[:, 2] ** 2).reshape((N, 1))
    u = u / norm
    xp[:, 0] = xxstg[0] - u[:, 0] * beta * xxstg[4]
    xp[:, 1] = xxstg[2] - u[:, 1] * beta * xxstg[4]
    xp[:, 2] = -u[:, 2] * beta * xxstg[4]
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


def astraBeam2particleArray(filename, print_params=True):
    """
    function convert Astra beam distribution to Ocelot format - ParticleArray.
    Note that downloading ParticleArray from the astra file and saving it back does not give the same distribution.
    The difference arises because the array of particles does not have a reference particle, and in this case
    the first particle is used as a reference.

    :param print_params:
    :type filename: str
    :return: ParticleArray
    """
    P0 = np.loadtxt(filename)

    # remove particles lost or not injected
    inds = np.argwhere(P0[:, 9] > 0)
    inds = inds.reshape(inds.shape[0])

    P0 = P0[inds,:]

    s_ref = P0[0, 2]
    Pref = P0[0, 5]

    if P0[0, 7] == 0:
        xp = P0[1:, :6]
        charge_array = -P0[1:, 7] * 1e-9  # charge in nC -> in C
    else:
        charge_array = -P0[:, 7] * 1e-9  # charge in nC -> in C
        xp = P0[:, :6]
        xp[0, 2] = 0.
        xp[0, 5] = 0.

    gamref = np.sqrt((Pref / m_e_eV) ** 2 + 1)
    xxstg = exact_xp_2_xxstg_mad(xp, gamref)

    p_array = ParticleArray(len(charge_array))
    p_array.s = s_ref
    p_array.E = np.sqrt((Pref / m_e_eV) ** 2 + 1) * m_e_GeV
    p_array.rparticles[0] = xxstg[:, 0]
    p_array.rparticles[1] = xxstg[:, 1]
    p_array.rparticles[2] = xxstg[:, 2]
    p_array.rparticles[3] = xxstg[:, 3]
    p_array.rparticles[4] = xxstg[:, 4]
    p_array.rparticles[5] = xxstg[:, 5]
    p_array.q_array = charge_array

    if print_params:
        print("Astra to Ocelot: charge = ", sum(charge_array), " C")
        print("Astra to Ocelot: particles number = ", len(charge_array))
        print("Astra to Ocelot: energy = ", p_array.E, " GeV")
        print("Astra to Ocelot: s pos = ", p_array.s, " m")

    return p_array


def particleArray2astraBeam(p_array, filename="tytest.ast"):
    """
    function convert  Ocelot's ParticleArray to Astra beam distribution and save to "filename".

    Note that downloading ParticleArray from the astra file and saving it back does not give the same distribution.
    The difference arises because the array of particles does not have a reference particle, and in this case
    the first particle is used as a reference.

    :param p_array:
    :param filename:
    :return:
    """

    gamref = p_array.E / m_e_GeV
    s0 = p_array.s
    P = p_array.rparticles.view()
    Np = int(P.size / 6)
    xp = exact_xxstg_2_xp_mad(P, gamref)
    Pref = np.sqrt(p_array.E ** 2 / m_e_GeV ** 2 - 1) * m_e_eV

    ref_particle=np.array([0, 0, s0, 0, 0, Pref])
    xp = np.vstack((ref_particle, xp))

    charge_array = -p_array.q_array.reshape(len(p_array.q_array), 1) * 1e+9  # charge in C -> in nC
    charge_array = np.vstack((0, charge_array))
    flag = np.ones((len(charge_array+1), 1))
    astra = np.append(xp, flag * 0, axis=1)  # time in [ns]
    astra = np.append(astra, charge_array, axis=1)
    astra = np.append(astra, flag, axis=1)  # 1 - electron, 2 - positron, 3 - protons and 4 - hydrogen ions.
    astra = np.append(astra, flag * 5, axis=1)  # 5 - standard particle
    np.savetxt(filename, astra, fmt='%.7e')


def emittance_analysis(fileprefix="Exfel", trace_space=True, s_offset=None):
    """
    To calculate emittance in the trace space the flag "Tr_EmitS=.T" is needed.

    The phase space is (x, px, y, py, z, pz)
    The trace space is (x, x’, y, y’, z, z’)
    :param fileprefix: file prefix to read files:   fileprefix + '.Xemit.001',
                                                    fileprefix + '.Yemit.001',
                                                    fileprefix + '.Zemit.001',
                                                    fileprefix + '.TRemit.001'
    :param trace_space: True, to calculate emittance in the trace space
    :param s_offset: None, if None use s coordinates from a file, if not None the s coordinate starts from "s_offset"
    :return: list of the twiss objects
    """

    optx = np.loadtxt(fileprefix + '.Xemit.001')
    opty = np.loadtxt(fileprefix + '.Yemit.001')
    optz = np.loadtxt(fileprefix + '.Zemit.001')

    r_E = optz[:, 2] * 1e6 + m_e_eV
    gamma = r_E / m_e_eV
    emitx = optx[:, 5] * 1e-6
    sigmax = optx[:, 3] * 1e-3
    betax = sigmax ** 2 / emitx * gamma
    emity = opty[:, 5] * 1e-6
    sigmay = opty[:, 3] * 1e-3
    betay = sigmay ** 2 / emity * gamma
    corx = optx[:, 6] * 1e-3 * sigmax
    alphax = -corx / emitx * gamma
    cory = opty[:, 6] * 1e-3 * sigmay
    alphay = -cory / emity * gamma

    if trace_space:
        optTS = np.loadtxt(fileprefix + '.TRemit.001')
        emitxTS = optTS[:, 3 - 1] * 1e-6
        emityTS = optTS[:, 4 - 1] * 1e-6
        betaxTS = sigmax ** 2 / emitxTS * gamma
        betayTS = sigmay ** 2 / emityTS * gamma
        betax = betaxTS
        betay = betayTS
        emitx = emitxTS
        emity = emityTS

    tws = []
    for i in range(len(optx[:, 0])):
        tw = Twiss()
        tw.beta_x = betax[i]
        tw.beta_y = betay[i]
        tw.alpha_x = alphax[i]
        tw.alpha_y = alphay[i]
        tw.gamma_x = (1 + tw.alpha_x * tw.alpha_x) / tw.beta_x
        tw.gamma_y = (1 + tw.alpha_y * tw.alpha_y) / tw.beta_y
        tw.emit_x = emitx[i]
        tw.emit_y = emity[i]
        tw.E = r_E * 1e-9  # in GeV

        if s_offset != None:
            tw.s = optx[i, 0] - optx[0, 0] + s_offset
        else:
            tw.s = optx[i, 0]

        tws.append(tw)

    return tws
