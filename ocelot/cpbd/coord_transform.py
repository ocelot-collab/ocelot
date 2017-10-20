
import numpy as np
from ocelot.common.globals import *
try:
    import numexpr as ne
    ne_flag = True
except:
    print("coord_transform.py: module NUMEXPR is not installed. Install it if you want higher speed calculation.")
    ne_flag = False

def xp_2_xxstg_mad(xp, xxstg, gamref):
    # to mad format
    N = xp.shape[1]
    pref = m_e_eV * np.sqrt(gamref ** 2 - 1)
    betaref = np.sqrt(1 - gamref ** -2)
    u = np.c_[xp[3], xp[4], xp[5]]
    if ne_flag:
        sum_u2 = ne.evaluate('sum(u * u, 1)')
        gamma = ne.evaluate('sqrt(1 + sum_u2 / m_e_eV ** 2)')
        beta = ne.evaluate('sqrt(1 - gamma ** -2)')
    else:
        gamma = np.sqrt(1 + np.sum(u * u, 1) / m_e_eV ** 2)
        beta = np.sqrt(1 - gamma ** -2)

    if np.__version__ > "1.8":
        p0 = np.linalg.norm(u, 2, 1).reshape((N, 1))
    else:
        p0 = np.sqrt(u[:, 0] ** 2 + u[:, 1] ** 2 + u[:, 2] ** 2).reshape((N, 1))
    u = u / p0
    u0 = u[:, 0]
    u1 = u[:, 1]
    u2 = u[:, 2]
    if ne_flag:
        xp0 = xp[0]
        xp1 = xp[1]
        xp2 = xp[2]
        cdt = ne.evaluate('-xp2 / (beta * u2)')
        xxstg[0] = ne.evaluate('xp0 + beta * u0 * cdt')
        xxstg[2] = ne.evaluate('xp1 + beta * u1 * cdt')
        xxstg[5] = ne.evaluate('(gamma / gamref - 1) / betaref')
    else:
        cdt = -xp[2] / (beta * u2)
        xxstg[0] = xp[0] + beta * u0 * cdt
        xxstg[2] = xp[1] + beta * u1 * cdt
        xxstg[5] = (gamma / gamref - 1) / betaref
    xxstg[4] = cdt
    xxstg[1] = xp[3] / pref
    xxstg[3] = xp[4] / pref
    return xxstg


def xxstg_2_xp_mad(xxstg, xp, gamref):
    # from mad format
    N = xxstg.shape[1]
    #pref = m_e_eV * np.sqrt(gamref ** 2 - 1)
    betaref = np.sqrt(1 - gamref ** -2)
    if ne_flag:
        xxstg1 = xxstg[1]
        xxstg3 = xxstg[3]
        xxstg5 = xxstg[5]
        gamma = ne.evaluate('(betaref * xxstg5 + 1) * gamref')
        beta = ne.evaluate('sqrt(1 - gamma ** -2)')
        pz2pref = ne.evaluate('sqrt(((gamma * beta) / (gamref * betaref)) ** 2 - xxstg1 ** 2 - xxstg3 ** 2)')
    else:
        gamma = (betaref * xxstg[5] + 1) * gamref
        beta = np.sqrt(1 - gamma ** -2)
        pz2pref = np.sqrt(((gamma * beta) / (gamref * betaref)) ** 2 - xxstg[1] ** 2 - xxstg[3] ** 2)

    u = np.c_[xxstg[1] / pz2pref, xxstg[3] / pz2pref, np.ones(N)]
    if np.__version__ > "1.8":
        norm = np.linalg.norm(u, 2, 1).reshape((N, 1))
    else:
        norm = np.sqrt(u[:, 0] ** 2 + u[:, 1] ** 2 + u[:, 2] ** 2).reshape((N, 1))
    u = u / norm
    u0 = u[:, 0]
    u1 = u[:, 1]
    u2 = u[:, 2]
    if ne_flag:
        xxstg0 = xxstg[0]
        xxstg2 = xxstg[2]
        xxstg4 = xxstg[4]
        xp[0] = ne.evaluate('xxstg0 - u0 * beta * xxstg4')
        xp[1] = ne.evaluate('xxstg2 - u1 * beta * xxstg4')
        xp[2] = ne.evaluate('-u2 * beta * xxstg4')
        xp[3] = ne.evaluate('u0 * gamma * beta * m_e_eV')
        xp[4] = ne.evaluate('u1 * gamma * beta * m_e_eV')
        xp[5] = ne.evaluate('u2 * gamma * beta * m_e_eV')
    else:
        xp[0] = xxstg[0] - u0 * beta * xxstg[4]
        xp[1] = xxstg[2] - u1 * beta * xxstg[4]
        xp[2] = -u2 * beta * xxstg[4]
        xp[3] = u0 * gamma * beta * m_e_eV
        xp[4] = u1 * gamma * beta * m_e_eV
        xp[5] = u2 * gamma * beta * m_e_eV
    return xp

