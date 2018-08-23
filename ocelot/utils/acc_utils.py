__author__ = 'Sergey Tomin'

import numpy as np

def RTU_56(LB, LD, r, m):
    LB2 = LB * LB
    r2 = r * r
    K = np.sqrt(1 - LB2 / r2)
    r56 = m * r * np.arcsin(LB / r) - m * LB / K - 2 * LB2 * LD / (K ** 3 * r2)

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
    :return: R56, T566, U5666, Sref
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
