__author__ = 'Sergey Tomin'

import logging
import numpy as np

from ocelot.common.globals import m_e_GeV, speed_of_light

logger = logging.getLogger(__name__)


def rot_mtx(angle):
    cs = np.cos(angle)
    sn = np.sin(angle)
    return np.array([[cs, 0., sn, 0., 0., 0.],
                     [0., cs, 0., sn, 0., 0.],
                     [-sn, 0., cs, 0., 0., 0.],
                     [0., -sn, 0., cs, 0., 0.],
                     [0., 0., 0., 0., 1., 0.],
                     [0., 0., 0., 0., 0., 1.]])


def uni_matrix(z, k1, hx, sum_tilts=0., energy=0.):
    """
    universal matrix. The function creates R-matrix from given parameters.
    r = element.l/element.angle
    +K - focusing lens, -K - defoc

    :param z: element length [m]
    :param k1: quadrupole strength [1/m**2]
    :param hx: the curvature (1/r) of the element [1/m]
    :param sum_tilts: rotation relative to longitudinal axis [rad]
    :param energy: the beam energy [GeV]
    :return: R-matrix [6, 6]
    """

    gamma = energy / m_e_GeV

    kx2 = (k1 + hx * hx)
    ky2 = -k1
    kx = np.sqrt(kx2 + 0.j)
    ky = np.sqrt(ky2 + 0.j)
    cx = np.cos(z * kx).real
    cy = np.cos(z * ky).real
    sy = (np.sin(ky * z) / ky).real if ky != 0 else z

    igamma2 = 0.

    if gamma != 0:
        igamma2 = 1. / (gamma * gamma)

    beta = np.sqrt(1. - igamma2)

    if kx != 0:
        sx = (np.sin(kx * z) / kx).real
        dx = hx / kx2 * (1. - cx)
        r56 = hx * hx * (z - sx) / kx2 / beta ** 2
    else:
        sx = z
        dx = z * z * hx / 2.
        r56 = hx * hx * z ** 3 / 6. / beta ** 2

    r56 -= z / (beta * beta) * igamma2

    u_matrix = np.array([[cx, sx, 0., 0., 0., dx / beta],
                         [-kx2 * sx, cx, 0., 0., 0., sx * hx / beta],
                         [0., 0., cy, sy, 0., 0.],
                         [0., 0., -ky2 * sy, cy, 0., 0.],
                         [hx * sx / beta, dx / beta, 0., 0., 1., r56],
                         [0., 0., 0., 0., 0., 1.]])
    if sum_tilts != 0:
        u_matrix = np.dot(np.dot(rot_mtx(-sum_tilts), u_matrix), rot_mtx(sum_tilts))
    return u_matrix

