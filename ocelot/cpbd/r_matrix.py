__author__ = 'Sergey Tomin'

import logging
from ocelot.common.globals import m_e_GeV, speed_of_light
from ocelot.cpbd.elements import *

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


def create_r_matrix(element):
    k1 = element.k1
    if element.l == 0:
        hx = 0.
    else:
        hx = element.angle / element.l

    r_z_e = lambda z, energy: uni_matrix(z, k1, hx=hx, sum_tilts=0, energy=energy)

    if element.__class__ == Edge:
        sec_e = 1. / np.cos(element.edge)
        phi = element.fint * element.h * element.gap * sec_e * (1. + np.sin(element.edge) ** 2)
        # phi = element.fint * element.h * element.gap * sec_e * (1. + np.sin(2*element.edge) )
        r = np.eye(6)
        r[1, 0] = element.h * np.tan(element.edge)
        r[3, 2] = -element.h * np.tan(element.edge - phi)
        r_z_e = lambda z, energy: r

    if element.__class__ in [Hcor, Vcor]:
        r_z_e = lambda z, energy: uni_matrix(z, 0, hx=0, sum_tilts=0, energy=energy)

    elif element.__class__ == Undulator:
        """
        in OCELOT coordinates:
        R56 = - Lu/(gamma**2 * beta**2) * (1 + 0.5 * K**2 * beta**2)
        S.Tomin, Varenna, 2017.
        """

        def undulator_r_z(z, lperiod, Kx, Ky, energy):
            gamma = energy / m_e_GeV
            r = np.eye(6)
            r[0, 1] = z
            if gamma != 0 and lperiod != 0 and Kx != 0:
                beta = 1 / np.sqrt(1.0 - 1.0 / (gamma * gamma))

                omega_x = np.sqrt(2.0) * np.pi * Kx / (lperiod * gamma * beta)
                omega_y = np.sqrt(2.0) * np.pi * Ky / (lperiod * gamma * beta)
                r[2, 2] = np.cos(omega_x * z)
                r[2, 3] = np.sin(omega_x * z) / omega_x
                r[3, 2] = -np.sin(omega_x * z) * omega_x
                r[3, 3] = np.cos(omega_x * z)

                r[4, 5] = - z / (gamma * beta) ** 2 * (1 + 0.5 * (Kx * beta) ** 2)

            else:
                r[2, 3] = z
            return r

        r_z_e = lambda z, energy: undulator_r_z(z, lperiod=element.lperiod, Kx=element.Kx, Ky=element.Ky, energy=energy)
        # b_z = lambda z, energy: dot((eye(6) - R_z(z, energy)), array([dx, 0., dy, 0., 0., 0.]))

    elif element.__class__ == Cavity:

        def cavity_R_z(z, V, E, freq, phi=0.):
            """
            :param z: length
            :param de: delta E
            :param freq: frequency
            :param E: initial energy
            :return: matrix
            """

            phi = phi * np.pi / 180.
            de = V * np.cos(phi)
            # pure pi-standing-wave case
            eta = 1
            # gamma = (E + 0.5 * de) / m_e_GeV
            Ei = E / m_e_GeV
            Ef = (E + de) / m_e_GeV
            Ep = (Ef - Ei) / z  # energy derivative
            if Ei == 0:
                logger.error("CAVITY: Initial energy is 0, check ParticleArray.E or Twiss.E OR cavity.v must be 0")

            cos_phi = np.cos(phi)
            alpha = np.sqrt(eta / 8.) / cos_phi * np.log(Ef / Ei)
            sin_alpha = np.sin(alpha)

            cos_alpha = np.cos(alpha)
            r11 = (cos_alpha - np.sqrt(2. / eta) * cos_phi * sin_alpha)

            if abs(Ep) > 1e-10:
                r12 = np.sqrt(8. / eta) * Ei / Ep * cos_phi * sin_alpha
            else:
                r12 = z
            r21 = -Ep / Ef * (cos_phi / np.sqrt(2. * eta) + np.sqrt(eta / 8.) / cos_phi) * sin_alpha

            r22 = Ei / Ef * (cos_alpha + np.sqrt(2. / eta) * cos_phi * sin_alpha)
            # print(f"z = {z}, V = {V}, E = {E}, phi = {phi}, alpha = {alpha}, r11 = {r11}")

            r56 = 0.
            beta0 = 1
            beta1 = 1

            k = 2. * np.pi * freq / speed_of_light
            r55_cor = 0.
            if V != 0 and E != 0:
                gamma2 = Ei * Ei
                beta0 = np.sqrt(1. - 1 / gamma2)
                gamma2 = Ef * Ef
                beta1 = np.sqrt(1. - 1 / gamma2)

                # r56 = (beta0 / beta1 - 1) * Ei / (Ef - Ei) * z
                r56 = - z / (Ef * Ef * Ei * beta1) * (Ef + Ei) / (beta1 + beta0)
                g0 = Ei
                g1 = Ef
                r55_cor = k * z * beta0 * V / m_e_GeV * np.sin(phi) * (g0 * g1 * (beta0 * beta1 - 1) + 1) / (
                        beta1 * g1 * (g0 - g1) ** 2)

            r66 = Ei / Ef * beta0 / beta1
            r65 = k * np.sin(phi) * V / (Ef * beta1 * m_e_GeV)
            cav_matrix = np.array([[r11, r12, 0., 0., 0., 0.],
                                   [r21, r22, 0., 0., 0., 0.],
                                   [0., 0., r11, r12, 0., 0.],
                                   [0., 0., r21, r22, 0., 0.],
                                   [0., 0., 0., 0., 1. + r55_cor, r56],
                                   [0., 0., 0., 0., r65, r66]]).real

            return cav_matrix

        if element.v == 0.:
            r_z_e = lambda z, energy: uni_matrix(z, 0., hx=0., sum_tilts=element.dtilt + element.tilt, energy=energy)
        else:
            r_z_e = lambda z, energy: cavity_R_z(z, V=element.v * z / element.l, E=energy, freq=element.freq,
                                                 phi=element.phi)

    elif element.__class__ == CouplerKick:

        def ck_matrix(v, phi, vxx, vxy, energy):
            """
            matrix for coupler kick

            :param v: voltage of the cavity in GV
            :param phi: phase [deg] of the cavity
            :param vxx: first order coefficients of the coupler kicks
            :param vxy: first order coefficients of the coupler kicks
            :param energy: beam energy in GeV
            :return:
            """
            phi = phi * np.pi / 180.
            m21 = (vxx * v * np.exp(1j * phi)).real / energy
            m43 = - m21
            m23 = (vxy * v * np.exp(1j * phi)).real / energy

            coupl_kick = np.array([[1, 0., 0., 0., 0., 0.],
                                    [m21, 1, m23, 0., 0., 0.],
                                    [0., 0., 1, 0., 0., 0.],
                                    [m23, 0., m43, 1, 0., 0.],
                                    [0., 0., 0., 0., 1., 0.],
                                    [0., 0., 0., 0., 0., 1]])
            return coupl_kick

        r_z_e = lambda z, energy: ck_matrix(v=element.v, phi=element.phi,
                                            vxx=element.vxx, vxy=element.vxy, energy=energy)

    elif element.__class__ == TWCavity:

        def tw_cavity_R_z(z, V, E, freq, phi=0.):
            """
            :param z: length
            :param de: delta E
            :param f: frequency
            :param E: initial energy
            :return: matrix
            """
            phi = phi * np.pi / 180.
            de = V * np.cos(phi)
            r12 = z * E / de * np.log(1. + de / E) if de != 0 else z
            r22 = E / (E + de)
            r65 = V * np.sin(phi) / (E + de) * (2 * np.pi / (speed_of_light / freq)) if freq != 0 else 0
            r66 = r22
            cav_matrix = np.array([[1, r12, 0., 0., 0., 0.],
                                   [0, r22, 0., 0., 0., 0.],
                                   [0., 0., 1, r12, 0., 0.],
                                   [0., 0., 0, r22, 0., 0.],
                                   [0., 0., 0., 0., 1., 0],
                                   [0., 0., 0., 0., r65, r66]]).real
            return cav_matrix

        def f_entrance(z, V, E, phi=0.):
            phi = phi * np.pi / 180.
            de = V * np.cos(phi)
            r = np.eye(6)
            r[1, 0] = -de / z / 2. / E
            r[3, 2] = r[1, 0]
            return r

        def f_exit(z, V, E, phi=0.):
            phi = phi * np.pi / 180.
            de = V * np.cos(phi)
            r = np.eye(6)
            r[1, 0] = +de / z / 2. / (E + de)
            r[3, 2] = r[1, 0]
            return r

        def cav(z, V, E, freq, phi):
            R_z = np.dot(tw_cavity_R_z(z, V, E, freq, phi), f_entrance(z, V, E, phi))
            R = np.dot(f_exit(z, V, E, phi), R_z)
            return R

        if element.v == 0.:
            r_z_e = lambda z, energy: uni_matrix(z, 0., hx=0., sum_tilts=element.dtilt + element.tilt, energy=energy)
        else:
            r_z_e = lambda z, energy: cav(z, V=element.v * z / element.l, E=energy, freq=element.freq,
                                          phi=element.phi)

    elif element.__class__ == Solenoid:
        def sol(l, k, energy):
            """
            K.Brown, A.Chao.
            :param l: efective length of solenoid
            :param k: B0/(2*Brho), B0 is field inside the solenoid, Brho is momentum of central trajectory
            :return: matrix
            """
            gamma = energy / m_e_GeV
            c = np.cos(l * k)
            s = np.sin(l * k)
            if k == 0:
                s_k = l
            else:
                s_k = s / k
            r56 = 0.
            if gamma != 0:
                gamma2 = gamma * gamma
                beta = np.sqrt(1. - 1. / gamma2)
                r56 -= l / (beta * beta * gamma2)
            sol_matrix = np.array([[c * c, c * s_k, s * c, s * s_k, 0., 0.],
                                   [-k * s * c, c * c, -k * s * s, s * c, 0., 0.],
                                   [-s * c, -s * s_k, c * c, c * s_k, 0., 0.],
                                   [k * s * s, -s * c, -k * s * c, c * c, 0., 0.],
                                   [0., 0., 0., 0., 1., r56],
                                   [0., 0., 0., 0., 0., 1.]]).real
            return sol_matrix

        r_z_e = lambda z, energy: sol(z, k=element.k, energy=energy)

    elif element.__class__ == TDCavity:
        """
        R - matrix for TDS - NOT TESTED
        """

        def tds_R_z(z, energy, freq, v, phi):
            """

            :param z:  length [m]
            :param freq: freq [Hz]
            :param v: voltage in [GeV]
            :param phi: phase [deg]
            :param energy: Energy in [GeV]
            :return:
            """
            phi = phi * np.pi / 180.

            gamma = energy / m_e_GeV
            igamma2 = 0.
            k0 = 2 * np.pi * freq / speed_of_light
            if gamma != 0:
                igamma2 = 1. / (gamma * gamma)
            if gamma > 1:
                pref = m_e_GeV * np.sqrt(gamma ** 2 - 1)
                K = v * k0 / pref
            else:
                K = 0.
            cos_phi = np.cos(phi)
            cos2_phi = np.cos(2 * phi)

            rm = np.eye(6)

            rm[0, 1] = z
            rm[0, 4] = -z * K * cos_phi / 2.
            rm[1, 4] = -K * cos_phi
            rm[2, 3] = z
            rm[4, 5] = - z * igamma2 / (1. - igamma2)
            rm[5, 0] = -rm[1, 4]
            rm[5, 1] = -rm[0, 4]
            rm[5, 4] = -z * K ** 2 * cos2_phi / 6
            return rm

        r_z_e = lambda z, energy: tds_R_z(z, energy, freq=element.freq, v=element.v * z / element.l, phi=element.phi)

    elif element.__class__ == Matrix:
        rm = np.eye(6)
        rm = element.r

        def r_matrix(z, l, rm):
            if z < l:
                r_z = uni_matrix(z, 0, hx=0)
            else:
                r_z = rm
            return r_z

        r_z_e = lambda z, energy: r_matrix(z, element.l, rm)

    elif element.__class__ == Multipole:
        r = np.eye(6)
        r[1, 0] = -element.kn[1]
        r[3, 2] = element.kn[1]
        r[1, 5] = element.kn[0]
        r_z_e = lambda z, energy: r

    elif element.__class__ == XYQuadrupole:
        k1 = element.k1

        if element.l == 0:
            hx = 0.
            hy = 0.
        else:
            hx = k1 * element.x_offs
            hy = -k1 * element.y_offs

        def r_mtx(z, k1, hx, hy, sum_tilts=0., energy=0.):
            # r = element.l/element.angle
            #  +K - focusing lens , -K - defoc
            gamma = energy / m_e_GeV

            kx2 = (k1 + hx * hx)
            ky2 = hy * hy - k1
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
                dy = hy / ky2 * (1. - cy)
                r56 = hx * hx * (z - sx) / kx2 / beta ** 2 + hy * hy * (z - sy) / ky2 / beta ** 2
            else:
                sx = z
                dx = z * z * hx / 2.
                dy = z * z * hy / 2.
                r56 = hx * hx * z ** 3 / 6. / beta ** 2 + hy * hy * z ** 3 / 6. / beta ** 2

            r56 -= z / (beta * beta) * igamma2

            u_matrix = np.array([[cx, sx, 0., 0., 0., dx / beta],
                                 [-kx2 * sx, cx, 0., 0., 0., sx * hx / beta],
                                 [0., 0., cy, sy, 0., dy / beta],
                                 [0., 0., -ky2 * sy, cy, 0., sy * hy / beta],
                                 [hx * sx / beta, dx / beta, hy * sy / beta, dy / beta, 1., r56],
                                 [0., 0., 0., 0., 0., 1.]])
            if sum_tilts != 0:
                u_matrix = np.dot(np.dot(rot_mtx(-sum_tilts), u_matrix), rot_mtx(sum_tilts))
            return u_matrix

        r_z_e = lambda z, energy: r_mtx(z, k1, hx=hx, hy=hy, sum_tilts=0, energy=energy)

    return r_z_e
