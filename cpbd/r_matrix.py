__author__ = 'Sergey Tomin'

import numpy as np
from ocelot.common.globals import m_e_GeV, speed_of_light
from ocelot.cpbd.elements import *


def rot_mtx(angle):
    cs = np.cos(angle)
    sn = np.sin(angle)
    return np.array([[cs, 0., sn, 0., 0., 0.],
                    [0.,  cs, 0., sn, 0., 0.],
                    [-sn, 0., cs, 0., 0., 0.],
                    [0., -sn, 0., cs, 0., 0.],
                    [0.,  0., 0., 0., 1., 0.],
                    [0.,  0., 0., 0., 0., 1.]])


def uni_matrix(z, k1, hx, sum_tilts=0., energy=0.):
    # r = element.l/element.angle
    # - K - focusing lens , +K - defoc
    gamma = energy/m_e_GeV

    kx2 = (k1 + hx*hx)
    ky2 = -k1
    kx = np.sqrt(kx2 + 0.j)
    ky = np.sqrt(ky2 + 0.j)
    cx = np.cos(z*kx).real
    cy = np.cos(z*ky).real
    sy = (np.sin(ky*z)/ky).real if ky != 0 else z

    igamma2 = 0.

    if gamma != 0:
        igamma2 = 1./(gamma*gamma)

    beta = np.sqrt(1. - igamma2)

    if kx != 0:
        sx = (np.sin(kx*z)/kx).real
        dx = hx/kx2*(1. - cx)
        r56 = hx*hx*(z - sx)/kx2/beta**2
    else:
        sx = z
        dx = z*z*hx/2.
        r56 = hx*hx*z**3/6./beta**2

    r56 -= z/(beta*beta)*igamma2
    u_matrix = np.array([[cx, sx, 0., 0., 0., dx/beta],
                        [-kx2*sx, cx, 0., 0., 0., sx*hx/beta],
                        [0., 0., cy, sy, 0., 0.],
                        [0., 0., -ky2*sy, cy, 0., 0.],
                        [hx*sx, dx, 0., 0., 1., r56],
                        [0., 0., 0., 0., 0., 1.]])
    if sum_tilts != 0:
        u_matrix = np.dot(np.dot(rot_mtx(-sum_tilts), u_matrix), rot_mtx(sum_tilts))
    return u_matrix


def create_r_matrix(element):

    #dx = element.dx
    #dy = element.dy
    #tilt = element.dtilt + element.tilt
    k1 = element.k1
    if element.l == 0:
        hx = 0.
    else:
        hx = element.angle / element.l

    r_z_e = lambda z, energy: uni_matrix(z, k1, hx=hx, sum_tilts=0, energy=energy)

    if element.__class__ == Edge:
        sec_e = 1. / np.cos(element.edge)
        phi = element.fint * element.h * element.gap * sec_e * (1. + np.sin(element.edge) ** 2)
        #phi = element.fint * element.h * element.gap * sec_e * (1. + np.sin(2*element.edge) )
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

                #print("here", r[2, 2], r[2, 3], r[3, 2], r[3, 3])
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
            :param f: frequency
            :param E: initial energy
            :return: matrix
            """
            phi = phi * np.pi / 180.
            de = V * np.cos(phi)
            # pure pi-standing-wave case
            eta = 1.0
            #gamma = (E + 0.5 * de) / m_e_GeV
            Ei = E / m_e_GeV
            Ef = (E + de) / m_e_GeV
            Ep = (Ef - Ei) / z  # energy derivative
            if Ei == 0:
                print("Warning! Initial energy is zero and cavity.delta_e != 0! Change Ei or cavity.delta_e must be 0")

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

            r56 = 0.
            if V != 0 and E != 0:
                #gamma2 = gamma * gamma
                #beta = np.sqrt(1. - 1 / gamma2)
                #r56 = -z / (beta * beta * gamma2)
                #gs = (Ef-Ei)/z
                #r56 = -(1.0/Ei-1.0/Ef)/gs

                gamma2 = Ei * Ei
                beta0 = np.sqrt(1. - 1 / gamma2)
                gamma2 = Ef * Ef
                beta1 = np.sqrt(1. - 1 / gamma2)
                r56 = (beta0 / beta1 - 1) * Ei / (Ef - Ei) * z
                # the same equation:
                # gamma = E/m_e_GeV
                # beta = np.sqrt(1. - 1 / gamma**2)
                # r56 = -1./(beta*gamma)**2 * z*(1 - 1.5*de/m_e_GeV/beta**2/gamma)

            k = 2.*np.pi*freq/speed_of_light
            r66 = Ei/Ef
            r65 = k*np.sin(phi)*V/(Ef*m_e_GeV)
            cav_matrix = np.array([[r11, r12, 0., 0., 0., 0.],
                                [r21, r22, 0., 0., 0., 0.],
                                [0., 0., r11, r12, 0., 0.],
                                [0., 0., r21, r22, 0., 0.],
                                [0., 0., 0., 0., 1., r56],
                                [0., 0., 0., 0., r65, r66]]).real
            if element.coupler_kick:
                #element.vxx_up = 1.0003 - 0.8132j
                #element.vxy_up = (3.4075 - 0.41223j)
                m21 = (element.vxx_up * V * np.exp(1j*phi)).real*1e-3 /E
                m43 = - m21
                m23 = (element.vxy_up* V * np.exp(1j*phi)).real*1e-3 /E

                coupl_kick_up = np.array([[1, 0., 0., 0., 0., 0.],
                                      [m21, 1, m23, 0., 0., 0.],
                                      [0., 0., 1, 0., 0., 0.],
                                      [m23, 0., m43, 1, 0., 0.],
                                      [0., 0., 0., 0., 1., 0.],
                                      [0., 0., 0., 0., 0., 1]]).real

                #vxx = ((-4.9278 - 2.2112j) * V * np.exp(1j*phi)).real*1e-3 /(E + de)
                #vyy = - vxx
                #vxy = ((2.9224 - 0.027228j) * V * np.exp(1j*phi)).real *1e-3 /(E + de)

                #element.vxx_down = (-4.9278 - 2.2112j)
                #element.vxy_down = (2.9224 - 0.027228j)
                m21 = (element.vxx_down * V * np.exp(1j*phi)).real*1e-3 /(E + de)
                m43 = - m21
                m23 = (element.vxy_down* V * np.exp(1j*phi)).real*1e-3 /(E + de)
                coupl_kick_down = np.array([[1, 0., 0., 0., 0., 0.],
                                      [m21, 1, m23, 0., 0., 0.],
                                      [0., 0., 1, 0., 0., 0.],
                                      [m23, 0., m43, 1, 0., 0.],
                                      [0., 0., 0., 0., 1., 0.],
                                      [0., 0., 0., 0., 0., 1]]).real
                return np.dot(np.dot(coupl_kick_up, cav_matrix), coupl_kick_down)
            return cav_matrix

        if element.delta_e == 0. and element.v == 0.:
            r_z_e = lambda z, energy: uni_matrix(z, 0., hx=0., sum_tilts=element.dtilt + element.tilt, energy=energy)
        else:
            r_z_e = lambda z, energy: cavity_R_z(z, V=element.v * z / element.l, E=energy, freq=element.f,
                                               phi=element.phi)

        #delta_e_z = lambda z: element.v * np.cos(element.phi * np.pi / 180.) * z / element.l
        #delta_e = element.v * np.cos(element.phi * np.pi / 180.)



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
                gamma2 = gamma*gamma
                beta = np.sqrt(1. - 1./gamma2)
                r56 -= l/(beta*beta*gamma2)
            sol_matrix = np.array([[c * c, c * s_k, s * c, s * s_k, 0., 0.],
                                [-k * s * c, c * c, -k * s * s, s * c, 0., 0.],
                                [-s * c, -s * s_k, c * c, c * s_k, 0., 0.],
                                [k * s * s, -s * c, -k * s * c, c * c, 0., 0.],
                                [0., 0., 0., 0., 1., r56],
                                [0., 0., 0., 0., 0., 1.]]).real
            return sol_matrix

        r_z_e = lambda z, energy: sol(z, k=element.k, energy=energy)

    elif element.__class__ == Matrix:
        rm = np.eye(6)
        rm[0, 0] = element.rm11
        rm[0, 1] = element.rm12
        rm[0, 2] = element.rm13
        rm[0, 3] = element.rm14
        rm[0, 4] = element.rm15
        rm[0, 5] = element.rm16

        rm[1, 0] = element.rm21
        rm[1, 1] = element.rm22
        rm[1, 2] = element.rm23
        rm[1, 3] = element.rm24
        rm[1, 4] = element.rm25
        rm[1, 5] = element.rm26

        rm[2, 0] = element.rm31
        rm[2, 1] = element.rm32
        rm[2, 2] = element.rm33
        rm[2, 3] = element.rm34
        rm[2, 4] = element.rm35
        rm[2, 5] = element.rm36

        rm[3, 0] = element.rm41
        rm[3, 1] = element.rm42
        rm[3, 2] = element.rm43
        rm[3, 3] = element.rm44
        rm[3, 4] = element.rm45
        rm[3, 5] = element.rm46

        rm[4, 0] = element.rm51
        rm[4, 1] = element.rm52
        rm[4, 2] = element.rm53
        rm[4, 3] = element.rm54
        rm[4, 4] = element.rm55
        rm[4, 5] = element.rm56

        rm[5, 0] = element.rm61
        rm[5, 1] = element.rm62
        rm[5, 2] = element.rm63
        rm[5, 3] = element.rm64
        rm[5, 4] = element.rm65
        rm[5, 5] = element.rm66

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

    # else:
    #    print (element.__class__, " : unknown type of magnetic element. Cannot create transfer map ")

    #b_z = lambda z, energy: dot((eye(6) - R_z(z, energy)), array([dx, 0., dy, 0., 0., 0.]))
    return r_z_e