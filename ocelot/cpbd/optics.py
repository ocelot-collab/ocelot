__author__ = 'Sergey'

from numpy.linalg import inv
from math import factorial
from ocelot.cpbd.beam import Particle, Twiss, ParticleArray
from ocelot.cpbd.physics_proc import RectAperture, EllipticalAperture
from ocelot.cpbd.high_order import *
from ocelot.cpbd.r_matrix import *
from copy import deepcopy
import logging
import numpy as np

_logger = logging.getLogger(__name__)
_logger_navi = logging.getLogger(__name__ + ".navi")

try:
    import numba as nb
    nb_flag = True
except ImportError as error:
    _logger.debug(" optics.py: module NUMBA is not installed. Install it to speed up calculation")
    nb_flag = False


class SecondOrderMult:
    """
    The class includes three different methods for transforming the particle
    coordinates:

    1. NUMBA module - DEACTIVATED, because new numpy implementation shows higher performance.
        Slightly faster than NUMPY for simulations with a large number of time
        steps. Uses full matrix multiplication.

    2. NUMPY module
        Base method to be used.
        Uses full matrix multiplication.

    """
    def __init__(self):
        self.full_matrix = False

        if nb_flag and False:
            # print("SecondTM: NUMBA")
            self.tmat_multip = nb.njit()(self.numba_apply)
        else:
            # print("SecondTM: Numpy")
            self.tmat_multip = self.numpy_apply

    @staticmethod
    def numba_apply(X, R, T):
        Xcopy = np.copy(X)
        for n in range(X.shape[1]):
            for i in range(6):
                x_new = 0
                for j in range(6):
                    x_new += R[i, j]*Xcopy[j, n]
                    for k in range(6):
                        x_new += T[i, j, k] * Xcopy[j, n] * Xcopy[k, n]
                X[i, n] = x_new

    def numpy_apply(self, X, R, T):
        X[:] = np.matmul(R, X) + np.einsum('ijk,j...,k...->i...', T, X, X)


def transform_vec_ent(X, dx, dy, tilt):
    rotmat = rot_mtx(tilt)
    x_add = np.add(X, np.array([[-dx], [0.], [-dy], [0.], [0.], [0.]]))
    X[:] = np.dot(rotmat, x_add)[:]
    return X


def transform_vec_ext(X, dx, dy, tilt):
    rotmat = rot_mtx(-tilt)
    x_tilt = np.dot(rotmat, X)
    X[:] = np.add(x_tilt, np.array([[dx], [0.], [dy], [0.], [0.], [0.]]))[:]
    return X


def transfer_maps_mult_py(Ra, Ta, Rb, Tb):
    """
    cell = [A, B]
    Rc = Rb * Ra
    :param Ra:
    :param Ta:
    :param Rb:
    :param Tb:
    :param sym_flag:
    :return:
    """
    Rc = np.dot(Rb, Ra)
    Tc = np.zeros((6, 6, 6))
    for i in range(6):
        for j in range(6):
            for k in range(6):
                t1 = 0.
                t2 = 0.
                for l in range(6):
                    t1 += Rb[i, l] * Ta[l, j, k]
                    for m in range(6):
                        t2 += Tb[i, l, m] * Ra[l, j] * Ra[m, k]
                Tc[i, j, k] = t1 + t2
    return Rc, Tc


transfer_maps_mult = transfer_maps_mult_py if nb_flag is not True else nb.jit(transfer_maps_mult_py)


def transfer_map_rotation(R, T, tilt):
    rotmat = rot_mtx(tilt)
    Rc, Tc = transfer_maps_mult(Ra=rotmat, Ta=np.zeros((6, 6, 6)), Rb=R, Tb=T)
    R, T = transfer_maps_mult(Ra=Rc, Ta=Tc, Rb=rot_mtx(-tilt), Tb=np.zeros((6, 6, 6)))
    return R, T


class TransferMap:
    """
    TransferMap is a basic linear transfer map for all elements.
    """
    def __init__(self):
        self.dx = 0.
        self.dy = 0.
        self.tilt = 0.
        self.length = 0
        self.hx = 0.
        self.delta_e = 0.0
        self.delta_e_z = lambda z: 0.0
        # 6x6 linear transfer matrix
        self.R = lambda energy: np.eye(6)
        self.R_z = lambda z, energy: np.zeros((6, 6))
        self.B_z = lambda z, energy: np.dot((np.eye(6) - self.R_z(z, energy)), np.array([[self.dx], [0.], [self.dy], [0.], [0.], [0.]]))
        self.B = lambda energy: self.B_z(self.length, energy)
        self.map = lambda u, energy: self.mul_p_array(u, energy=energy)

    def map_x_twiss(self, tws0):
        E = tws0.E
        M = self.R(E)
        zero_tol = 1.e-10
        if abs(self.delta_e) > zero_tol:
            Ei = tws0.E
            Ef = tws0.E + self.delta_e
            k = np.sqrt(Ef / Ei)
            M[0, 0] = M[0, 0] * k
            M[0, 1] = M[0, 1] * k
            M[1, 0] = M[1, 0] * k
            M[1, 1] = M[1, 1] * k
            M[2, 2] = M[2, 2] * k
            M[2, 3] = M[2, 3] * k
            M[3, 2] = M[3, 2] * k
            M[3, 3] = M[3, 3] * k
            E = Ef

        m = tws0
        tws = Twiss(tws0)
        tws.E = E
        tws.p = m.p
        tws.beta_x = M[0, 0] * M[0, 0] * m.beta_x - 2 * M[0, 1] * M[0, 0] * m.alpha_x + M[0, 1] * M[0, 1] * m.gamma_x
        # tws.beta_x = ((M[0,0]*tws.beta_x - M[0,1]*m.alpha_x)**2 + M[0,1]*M[0,1])/m.beta_x
        tws.beta_y = M[2, 2] * M[2, 2] * m.beta_y - 2 * M[2, 3] * M[2, 2] * m.alpha_y + M[2, 3] * M[2, 3] * m.gamma_y
        # tws.beta_y = ((M[2,2]*tws.beta_y - M[2,3]*m.alpha_y)**2 + M[2,3]*M[2,3])/m.beta_y
        tws.alpha_x = -M[0, 0] * M[1, 0] * m.beta_x + (M[0, 1] * M[1, 0] + M[1, 1] * M[0, 0]) * m.alpha_x - M[0, 1] * M[
            1, 1] * m.gamma_x
        tws.alpha_y = -M[2, 2] * M[3, 2] * m.beta_y + (M[2, 3] * M[3, 2] + M[3, 3] * M[2, 2]) * m.alpha_y - M[2, 3] * M[
            3, 3] * m.gamma_y

        tws.gamma_x = (1. + tws.alpha_x * tws.alpha_x) / tws.beta_x
        tws.gamma_y = (1. + tws.alpha_y * tws.alpha_y) / tws.beta_y

        tws.Dx = M[0, 0] * m.Dx + M[0, 1] * m.Dxp + M[0, 5]
        tws.Dy = M[2, 2] * m.Dy + M[2, 3] * m.Dyp + M[2, 5]

        tws.Dxp = M[1, 0] * m.Dx + M[1, 1] * m.Dxp + M[1, 5]
        tws.Dyp = M[3, 2] * m.Dy + M[3, 3] * m.Dyp + M[3, 5]
        denom_x = M[0, 0] * m.beta_x - M[0, 1] * m.alpha_x
        if denom_x == 0.:
            d_mux = np.pi / 2. * M[0, 1] / np.abs(M[0, 1])
        else:
            d_mux = np.arctan(M[0, 1] / denom_x)

        if d_mux < 0:
            d_mux += np.pi
        tws.mux = m.mux + d_mux
        denom_y = M[2, 2] * m.beta_y - M[2, 3] * m.alpha_y
        if denom_y == 0.:
            d_muy = np.pi / 2. * M[2, 3] / np.abs(M[2, 3])
        else:
            d_muy = np.arctan(M[2, 3] / denom_y)
        if d_muy < 0:
            d_muy += np.pi
        tws.muy = m.muy + d_muy

        return tws

    def mul_p_array(self, rparticles, energy=0.):
        a = np.add(np.dot(self.R(energy), rparticles), self.B(energy))
        rparticles[:] = a[:]
        return rparticles

    def __mul__(self, m):
        """
        :param m: TransferMap, Particle or Twiss
        :return: TransferMap, Particle or Twiss
        Ma = {Ba, Ra, Ta}
        Mb = {Bb, Rb, Tb}
        X1 = R*(X0 - dX) + dX = R*X0 + B
        B = (E - R)*dX
        """

        if m.__class__ in [TransferMap]:
            m2 = TransferMap()
            m2.R = lambda energy: np.dot(self.R(energy), m.R(energy))
            m2.B = lambda energy: np.dot(self.R(energy), m.B(energy)) + self.B(energy)
            m2.length = m.length + self.length

            return m2

        elif m.__class__ == Particle:
            self.apply(m)
            return deepcopy(m)

        elif m.__class__ == Twiss:

            tws = self.map_x_twiss(m)
            # trajectory
            # X0 = array([m.x, m.xp, m.y, m.yp, m.tau, m.p])
            # tws.x, tws.xp, tws.y, tws.yp, tws.tau, tws.dE = self.mul_p_array(X0, energy=tws.E, order=1)
            tws.s = m.s + self.length
            return tws

        else:
            _logger.error(" TransferMap.__mul__: unknown object in transfer map multiplication: " + str(m.__class__.__name__))
            raise Exception(" TransferMap.__mul__: unknown object in transfer map multiplication: " + str(m.__class__.__name__))

    def apply(self, prcl_series):
        """
        :param prcl_series: can be list of Particles [Particle_1, Particle_2, ... ] or ParticleArray
        :return: None
        """
        if prcl_series.__class__ == ParticleArray:
            self.map(prcl_series.rparticles, energy=prcl_series.E)
            prcl_series.E += self.delta_e
            prcl_series.s += self.length

        elif prcl_series.__class__ == Particle:
            p = prcl_series
            p.x, p.px, p.y, p.py, p.tau, p.p = self.map(np.array([[p.x], [p.px], [p.y], [p.py], [p.tau], [p.p]]), p.E)[:,0]
            p.s += self.length
            p.E += self.delta_e

        elif prcl_series.__class__ == list and prcl_series[0].__class__ == Particle:
            # If the energy is not the same (p.E) for all Particles in the list of Particles
            # in that case cycle is applied. For particles with the same energy p.E
            list_e = np.array([p.E for p in prcl_series])
            if False in (list_e[:] == list_e[0]):
                for p in prcl_series:
                    self.map(np.array([[p.x], [p.px], [p.y], [p.py], [p.tau], [p.p]]), energy=p.E)
                    p.E += self.delta_e
                    p.s += self.length
            else:
                pa = ParticleArray()
                pa.list2array(prcl_series)
                pa.E = prcl_series[0].E
                self.map(pa.rparticles, energy=pa.E)
                pa.E += self.delta_e
                pa.s += self.length
                pa.array2ex_list(prcl_series)

        else:
            _logger.error(" TransferMap.apply(): Unknown type of Particle_series: " + str(prcl_series.__class__.__name))
            raise Exception(" TransferMap.apply(): Unknown type of Particle_series: " + str(prcl_series.__class__.__name))

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.delta_e = m.delta_e_z(s)
        m.map = lambda u, energy: m.mul_p_array(u, energy=energy)
        return m


class SecondTM(TransferMap):
    def __init__(self, r_z_no_tilt, t_mat_z_e):
        TransferMap.__init__(self)
        self.r_z_no_tilt = r_z_no_tilt
        self.t_mat_z_e = t_mat_z_e

        self.multiplication = None
        self.map = lambda X, energy: self.t_apply(self.r_z_no_tilt(self.length, energy),
                                                  self.t_mat_z_e(self.length, energy), X, self.dx, self.dy, self.tilt)

        self.R_tilt = lambda energy: np.dot(np.dot(rot_mtx(-self.tilt), self.r_z_no_tilt(self.length, energy)), rot_mtx(self.tilt))

        self.T_tilt = lambda energy: transfer_map_rotation(self.r_z_no_tilt(self.length, energy),
                                                             self.t_mat_z_e(self.length, energy), self.tilt)[1]

    def t_apply(self, R, T, X, dx, dy, tilt, U5666=0.):
        if dx != 0 or dy != 0 or tilt != 0:
            X = transform_vec_ent(X, dx, dy, tilt)
        self.multiplication(X, R, T)
        if dx != 0 or dy != 0 or tilt != 0:
            X = transform_vec_ext(X, dx, dy, tilt)

        return X

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.T = lambda s, energy: m.t_mat_z_e(s, energy)
        m.delta_e = m.delta_e_z(s)
        m.map = lambda X, energy: m.t_apply(m.r_z_no_tilt(s, energy), m.t_mat_z_e(s, energy), X, m.dx, m.dy, m.tilt)
        return m


class CorrectorTM(SecondTM):
    def __init__(self, angle_x=0., angle_y=0., r_z_no_tilt=None, t_mat_z_e=None):
        SecondTM.__init__(self, r_z_no_tilt, t_mat_z_e)
        self.angle_x = angle_x
        self.angle_y = angle_y
        self.map = lambda X, energy: self.kick(X, self.length, self.length, self.angle_x, self.angle_y, energy)
        self.B_z = lambda z, energy: self.kick_b(z, self.length, angle_x, angle_y)

    def kick_b(self, z, l, angle_x, angle_y):
        if l == 0:
            hx = 0.
            hy = 0.
        else:
            hx = angle_x / l
            hy = angle_y / l

        dx = hx * z * z / 2.
        dy = hy * z * z / 2.
        dx1 = hx * z if l != 0 else angle_x
        dy1 = hy * z if l != 0 else angle_y
        b = np.array([[dx], [dx1], [dy], [dy1], [0.], [0.]])
        return b

    def kick(self, X, z, l, angle_x, angle_y, energy):
        _logger.debug('invoking kick_b')
        b = self.kick_b(z, l, angle_x, angle_y)
        if self.multiplication is not None and self.t_mat_z_e is not None:
            X1 = self.t_apply(R=self.R(energy), T=self.t_mat_z_e(z, energy), X=X, dx=0, dy=0, tilt=0, U5666=0.)
        else:
            X1 = np.dot(self.R(energy), X)
        X1 = np.add(X1, b)
        X[:] = X1[:]
        return X

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.delta_e = m.delta_e_z(s)
        m.map = lambda X, energy: m.kick(X, s, self.length, m.angle_x, m.angle_y, energy)
        return m


class PulseTM(TransferMap):
    def __init__(self, kn):
        TransferMap.__init__(self)

    def mul_parray(self, rparticles , energy=0.):
        n = len(rparticles)
        #if 'pulse' in self.__dict__:
        _logger.debug('TD transfer map')
        if n > 6: _logger.debug(
                'warning: time-dependent transfer maps not implemented for an array. Using 1st particle value')
        if n > 6: _logger.debug('warning: time-dependent transfer maps not implemented for steps inside element')
        tau = rparticles[4]
        dxp = self.pulse.kick_x(tau)
        dyp = self.pulse.kick_y(tau)
        _logger.debug('kick ' + str(dxp) + ' ' + str(dyp))
        b = np.array([0.0, dxp, 0.0, dyp, 0., 0.])
        #a = np.add(np.transpose(dot(self.R(energy), np.transpose(particles.reshape(int(n / 6), 6)))), b).reshape(n)
        a = np.add(np.dot(self.R(energy), rparticles), b)
        rparticles[:] = a[:]
        _logger.debug('return trajectory, array ' + str(len(rparticles)))
        return rparticles


class MultipoleTM(TransferMap):
    def __init__(self, kn):
        TransferMap.__init__(self)
        self.kn = kn
        self.map = lambda X, energy: self.kick(X, self.kn)

    def kick(self, X, kn):
        p = -kn[0] * X[5] + 0j
        for n in range(1, len(kn)):
            p += kn[n] * (X[0] + 1j * X[2]) ** n / factorial(n)
        X[1] = X[1] - np.real(p)
        X[3] = X[3] + np.imag(p)
        X[4] = X[4] - kn[0] * X[0]
        return X

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.delta_e = m.delta_e_z(s)
        m.map = lambda X, energy: m.kick(X, m.kn)
        return m


class CavityTM(TransferMap):
    def __init__(self, v=0, freq=0., phi=0.):
        TransferMap.__init__(self)
        self.v = v
        self.freq = freq
        self.phi = phi
        self.vx_up = 0.
        self.vy_up = 0.
        self.vx_down = 0.
        self.vy_down = 0.
        phase_term = np.cos(self.phi * np.pi / 180.)
        self.delta_e_z = lambda z: phase_term * self.v * z / self.length if self.length != 0 else phase_term * self.v
        self.delta_e = self.v * np.cos(self.phi * np.pi / 180.)
        self.map = lambda X, energy: self.map4cav(X, energy, self.v, self.freq, self.phi, self.length)

    def map4cav(self, X, E, V, freq, phi, z=0):
        beta0 = 1
        igamma2 = 0
        g0 = 1e10
        if E != 0:
            g0 = E / m_e_GeV
            igamma2 = 1. / (g0 * g0)
            beta0 = np.sqrt(1. - igamma2)

        phi = phi * np.pi / 180.

        X4 = np.copy(X[4])
        X5 = np.copy(X[5])
        X = self.mul_p_array(X, energy=E)  # t_apply(R, T, X, dx, dy, tilt)
        delta_e = V * np.cos(phi)

        T566 = 1.5 * z*igamma2/(beta0**3)
        T556 = 0.
        T555 = 0.
        if E + delta_e > 0:
            k = 2. * np.pi * freq / speed_of_light
            E1 = E + delta_e
            g1 = E1 / m_e_GeV
            beta1 = np.sqrt(1. - 1. / (g1 * g1))

            X[5] = X5 * E*beta0/(E1*beta1) + V*beta0 / (E1*beta1) * (np.cos(-X4*beta0 * k + phi) - np.cos(phi))

            dgamma = V / m_e_GeV
            if delta_e > 0:
                T566 = z * (beta0**3*g0**3 - beta1**3*g1**3)/(2*beta0*beta1**3*g0*(g0 - g1)*g1**3)
                T556 = beta0 * k * z * dgamma *g0 * (beta1**3*g1**3 + beta0 * (g0 - g1**3)) * np.sin(phi)/ (beta1**3 * g1**3 * (g0 - g1)**2)
                T555 = beta0**2 * k**2 * z * dgamma/2.*(dgamma*(2*g0*g1**3*(beta0*beta1**3 - 1) + g0**2 + 3*g1**2 - 2)/(beta1**3*g1**3*(g0 - g1)**3)*np.sin(phi)**2 -
                                                    (g1*g0*(beta1*beta0 - 1) + 1)/(beta1*g1*(g0 - g1)**2)*np.cos(phi))
        X[4] += T566 * X5*X5 + T556*X4*X5 + T555 * X4*X4

        return X

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.delta_e = m.delta_e_z(s)
        v = m.v * s / self.length if self.length != 0 else m.v
        m.map = lambda X, energy: m.map4cav(X, energy, v, m.freq, m.phi, s)
        return m


class CouplerKickTM(TransferMap):
    def __init__(self, v=0, freq=0., phi=0., vx=0., vy=0.):
        TransferMap.__init__(self)
        self.v = v
        self.freq = freq
        self.phi = phi
        self.vx = vx
        self.vy = vy
        self.map = lambda X, energy: self.kick(X, self.v, self.phi, energy)
        self.B_z = lambda z, energy: self.kick_b(self.v, self.phi, energy)

    def kick_b(self, v, phi, energy):
        phi = phi * np.pi / 180.
        dxp = (self.vx * v * np.exp(1j * phi)).real / energy
        dyp = (self.vy * v * np.exp(1j * phi)).real / energy
        b = np.array([[0.], [dxp], [0.], [dyp], [0.], [0.]])
        return b

    def kick(self, X, v, phi, energy):
        _logger.debug('invoking coupler kick zero order')
        b = self.kick_b(v, phi, energy)
        X1 = np.dot(self.R(energy), X)
        X1 = np.add(X1, b)
        X[:] = X1[:]
        return X

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.delta_e = m.delta_e_z(s)
        m.map = lambda X, energy: m.kick(X, self.v, self.phi, energy)
        return m


class KickTM(TransferMap):
    def __init__(self, angle=0., k1=0., k2=0., k3=0., nkick=1):
        TransferMap.__init__(self)
        self.angle = angle
        self.k1 = k1
        self.k2 = k2
        self.k3 = k3
        self.nkick = nkick

    def kick(self, X, l, angle, k1, k2, k3, energy, nkick=1):
        """
        does not work for dipole
        """
        gamma = energy / m_e_GeV
        coef = 0.
        beta = 1.
        if gamma != 0:
            gamma2 = gamma * gamma
            beta = 1. - 0.5 / gamma2
            coef = 1. / (beta * beta * gamma2)
        l = l / nkick
        angle = angle / nkick

        dl = l / 2.
        k1 = k1 * l
        k2 = k2 * l/2.
        k3 = k3 * l/6.

        for i in range(nkick):
            x = X[0] + X[1] * dl
            y = X[2] + X[3] * dl

            p = -angle * X[5] + 0j
            xy1 = x + 1j * y
            xy2 = xy1 * xy1
            xy3 = xy2 * xy1
            p += k1 * xy1 + k2 * xy2 + k3 * xy3
            X[1] = X[1] - np.real(p)
            X[3] = X[3] + np.imag(p)
            X[4] = X[4] + np.real(angle * xy1)/beta - X[5] * dl * coef

            X[0] = x + X[1] * dl
            X[2] = y + X[3] * dl
            #X[4] -= X[5] * dl * coef
        return X

    def kick_apply(self, X, l, angle, k1, k2, k3, energy, nkick, dx, dy, tilt):
        if dx != 0 or dy != 0 or tilt != 0:
            X = transform_vec_ent(X, dx, dy, tilt)
        self.kick(X, l, angle, k1, k2, k3, energy, nkick=nkick)
        if dx != 0 or dy != 0 or tilt != 0:
            X = transform_vec_ext(X, dx, dy, tilt)

        return X

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.delta_e = m.delta_e_z(s)
        m.map = lambda X, energy: m.kick_apply(X, s, m.angle, m.k1, m.k2, m.k3, energy, m.nkick, m.dx, m.dy, m.tilt)
        return m


class UndulatorTestTM(TransferMap):
    def __init__(self, lperiod, Kx, ax=0, ndiv=10):
        TransferMap.__init__(self)
        self.lperiod = lperiod
        self.Kx = Kx
        self.ax = ax
        self.ndiv = ndiv
        self.map = lambda X, energy: self.map4undulator(X, self.length, self.lperiod, self.Kx, self.ax, energy,
                                                        self.ndiv)

    def map4undulator(self, u, z, lperiod, Kx, ax, energy, ndiv):
        kz = 2. * np.pi / lperiod
        if ax == 0:
            kx = 0
        else:
            kx = 2. * np.pi / ax
        zi = np.linspace(0., z, num=ndiv)
        h = zi[1] - zi[0]
        kx2 = kx * kx
        kz2 = kz * kz
        ky2 = kz * kz + kx * kx
        ky = np.sqrt(ky2)
        gamma = energy / m_e_GeV
        h0 = 0.
        if gamma != 0:
            h0 = 1. / (gamma / Kx / kz)
        h02 = h0 * h0
        h = h / (1. + u[5])
        x = u[0]
        y = u[2]
        for z in range(len(zi) - 1):
            chx = np.cosh(kx * x)
            chy = np.cosh(ky * y)
            shx = np.sinh(kx * x)
            shy = np.sinh(ky * y)
            u[1] -= h / 2. * chx * shx * (kx * ky2 * chy * chy + kx2 * kx * shy * shy) / (ky2 * kz2) * h02
            u[3] -= h / 2. * chy * shy * (ky2 * chx * chx + kx2 * shx * shx) / (ky * kz2) * h02
            u[4] -= h / 2. / (1. + u[5]) * ((u[1] * u[1] + u[3] * u[3]) + chx * chx * chy * chy / (
                2. * kz2) * h02 + shx * shx * shy * shy * kx2 / (2. * ky2 * kz2) * h02)
            u[0] = x + h * u[1]
            u[2] = y + h * u[3]
        return u

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.delta_e = m.delta_e_z(s)
        m.map = lambda X, energy: m.map4undulator(X, m.length, m.lperiod, m.Kx, m.ax, energy, m.ndiv)
        return m


class RungeKuttaTM(TransferMap):
    def __init__(self, s_start=0, npoints=200):
        TransferMap.__init__(self)
        self.s_start = s_start
        self.npoints = npoints
        self.long_dynamics = True
        self.mag_field = lambda x, y, z: (0, 0, 0)
        self.map = lambda X, energy: rk_field(X, self.s_start, self.length, self.npoints, energy, self.mag_field,
                                              self.long_dynamics)

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.delta_e = m.delta_e_z(s)
        m.map = lambda X, energy: rk_field(X, m.s_start, s, m.npoints, energy, m.mag_field, m.long_dynamics)
        return m


class RungeKuttaTrTM(RungeKuttaTM):
    """
    THe same method as RungeKuttaTM but only transverse dynamics is included, longitudinal dynamics is skipped
    """
    def __init__(self, s_start=0, npoints=200):
        RungeKuttaTM.__init__(self, s_start=s_start, npoints=npoints)
        self.long_dynamics = False


class TWCavityTM(TransferMap):
    def __init__(self, l=0, v=0, phi=0, freq=0):
        TransferMap.__init__(self)
        self.length = l
        self.dx = 0
        self.dy = 0
        self.tilt = 0
        self.v = v
        self.phi = phi
        self.freq = freq
        self.delta_e_z = lambda z: self.v * np.cos(self.phi * np.pi / 180.) * z / self.length
        self.delta_e = self.v * np.cos(self.phi * np.pi / 180.)
        self.R_z = lambda z, energy: np.dot(
            self.tw_cavity_R_z(z, self.v * z / self.length, energy, self.freq, self.phi),
            self.f_entrance(z, self.v * z / self.length, energy, self.phi))
        self.R = lambda energy: np.dot(self.f_exit(self.length, self.v, energy, self.phi),
                                       self.R_z(self.length, energy))

    def tw_cavity_R_z(self, z, V, E, freq, phi=0.):
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

    def f_entrance(self, z, V, E, phi=0.):
        phi = phi * np.pi / 180.
        de = V * np.cos(phi)
        r = np.eye(6)
        r[1, 0] = -de / z / 2. / E
        r[3, 2] = r[1, 0]
        return r

    def f_exit(self, z, V, E, phi=0.):
        phi = phi * np.pi / 180.
        de = V * np.cos(phi)
        r = np.eye(6)
        r[1, 0] = +de / z / 2. / (E + de)
        r[3, 2] = r[1, 0]
        return r

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.delta_e = m.delta_e_z(s)
        m.map = lambda u, energy: m.mul_p_array(u, energy=energy)
        return m


class MethodTM:
    """
    The class creates a transfer map for elements that depend on user-defined parameters ("parameters").
    By default, the parameters = {"global": TransferMap}, which means that all elements will have linear transfer maps.
    You can also specify different transfer maps for any type of element.

    Example:
    --------
    # use linear matrices for all elements except Sextupole which will have nonlinear kick map (KickTM)
    method = MethodTM()
    method.global_method = TransferMap
    method.params[Sextupole] = KickTM

    # All elements are assigned matrices of the second order.
    # For elements for which there are no matrices of the second order are assigned default matrices, e.g. linear matrices.
    method2 = MethodTM()
    method2.global_method = SecondTM

    """
    def __init__(self, params=None):
        if params is None:
            self.params = {'global': TransferMap}
        else:
            self.params = params

        if "global" in self.params:
            self.global_method = self.params['global']
        else:
            self.global_method = TransferMap
        self.sec_order_mult = SecondOrderMult()
        self.nkick = self.params['nkick'] if 'nkick' in self.params else 1

    def create_tm(self, element):

        if element.__class__ in self.params:
            transfer_map = self.set_tm(element, self.params[element.__class__])
        else:
            transfer_map = self.set_tm(element, self.global_method)
        return transfer_map

    def set_tm(self, element, method):
        dx = element.dx
        dy = element.dy
        tilt = element.dtilt + element.tilt
        if element.l == 0:
            hx = 0.
        else:
            hx = element.angle / element.l

        r_z_e = create_r_matrix(element)

        # global method
        if method == KickTM:
            try:
                k3 = element.k3
            except:
                k3 = 0.
            tm = KickTM(angle=element.angle, k1=element.k1, k2=element.k2, k3=k3, nkick=self.nkick)

        elif method == SecondTM:

            T_z_e = lambda z, energy: t_nnn(z, hx, element.k1, element.k2, energy)

            if element.__class__ == Edge:
                if element.pos == 1:

                    _, T = fringe_ent(h=element.h, k1=element.k1, e=element.edge, h_pole=element.h_pole,
                                      gap=element.gap, fint=element.fint)
                else:

                    _, T = fringe_ext(h=element.h, k1=element.k1, e=element.edge, h_pole=element.h_pole,
                                      gap=element.gap, fint=element.fint)
                T_z_e = lambda z, energy: T
            if element.__class__ == XYQuadrupole:
                T = np.zeros((6, 6, 6))
            if element.__class__ == Matrix:
                T_z_e = lambda z, energy: element.t

            tm = SecondTM(r_z_no_tilt=r_z_e, t_mat_z_e=T_z_e)
            tm.multiplication = self.sec_order_mult.tmat_multip

        elif method == TWCavityTM:
            tm = TWCavityTM(l=element.l, v=element.v, phi=element.phi, freq=element.freq)
            return tm

        else:
            tm = TransferMap()

        if element.__class__ == Undulator and method == UndulatorTestTM:
            try:
                ndiv = element.ndiv
            except:
                ndiv = 5
            tm = UndulatorTestTM(lperiod=element.lperiod, Kx=element.Kx, ax=element.ax, ndiv=ndiv)

        if method in [RungeKuttaTM, RungeKuttaTrTM]:
            try:
                s_start = element.s_start
            except:
                s_start = 0.
            try:
                npoints = element.npoints
            except:
                npoints = 200
            tm = method(s_start=s_start, npoints=npoints)
            tm.mag_field = element.mag_field

        if element.__class__ == Cavity:
            tm = CavityTM(v=element.v, freq=element.freq, phi=element.phi)

        if element.__class__ == CouplerKick:
            tm = CouplerKickTM(v=element.v, freq=element.freq, phi=element.phi, vx=element.vx, vy=element.vy)

        if element.__class__ == TWCavity:
            tm = TWCavityTM(v=element.v, freq=element.freq, phi=element.phi)

        if element.__class__ == Matrix:
            tm.delta_e = element.delta_e
            tm.B_z = lambda z, energy: element.b
            tm.B = lambda energy: element.b

        if element.__class__ == Multipole:
            tm = MultipoleTM(kn=element.kn)

        if element.__class__ == Hcor:
            t_mat_z_e = lambda z, energy: t_nnn(z, 0, 0, 0, energy)
            tm = CorrectorTM(angle_x=element.angle, angle_y=0., r_z_no_tilt=r_z_e, t_mat_z_e=t_mat_z_e)
            tm.multiplication = self.sec_order_mult.tmat_multip

        if element.__class__ == Vcor:
            t_mat_z_e = lambda z, energy: t_nnn(z, 0, 0, 0, energy)
            tm = CorrectorTM(angle_x=0, angle_y=element.angle, r_z_no_tilt=r_z_e, t_mat_z_e=t_mat_z_e)
            tm.multiplication = self.sec_order_mult.tmat_multip

        tm.length = element.l
        tm.dx = dx
        tm.dy = dy
        tm.tilt = tilt
        tm.R_z = lambda z, energy: np.dot(np.dot(rot_mtx(-tilt), r_z_e(z, energy)), rot_mtx(tilt))
        tm.R = lambda energy: tm.R_z(element.l, energy)
        # tm.B_z = lambda z, energy: dot((eye(6) - tm.R_z(z, energy)), array([dx, 0., dy, 0., 0., 0.]))
        # tm.B = lambda energy: tm.B_z(element.l, energy)

        return tm


def sym_matrix(T):
    for i in range(6):
        for j in range(6):
            for k in range(j, 6):
                if j != k:
                    a = T[i, j, k] / 2.
                    T[i, k, j] = a
                    T[i, j, k] = a
    return T


def unsym_matrix(T):
    for i in range(6):
        for j in range(6):
            for k in range(j, 6):
                if j != k:
                    a = T[i, j, k] * 2.
                    T[i, k, j] = 0
                    T[i, j, k] = a
    return T


def lattice_transfer_map(lattice, energy):
    """
    Function calculates transfer maps, the first and second orders (R, T), for the whole lattice.
    Second order matrices are attached to lattice object:
    lattice.T_sym - symmetric second order matrix
    lattice.T - second order matrix
    lattice.R - linear R matrix

    :param lattice: MagneticLattice
    :param energy: the initial electron beam energy [GeV]
    :return: R - matrix
    """

    Ra = np.eye(6)
    Ta = np.zeros((6, 6, 6))
    Ba = np.zeros((6, 1))
    E = energy
    for i, elem in enumerate(lattice.sequence):
        Rb = elem.transfer_map.R(E)
        Bb = elem.transfer_map.B(E)

        if isinstance(elem.transfer_map, SecondTM):  # elem.transfer_map.__class__ == SecondTM:
            Tb = np.copy(elem.transfer_map.T_tilt(E))
            Tb = sym_matrix(Tb)
            Ra, Ta = transfer_maps_mult(Ra, Ta, Rb, Tb)
        else:
            Ra, Ta = transfer_maps_mult(Ra, Ta, Rb, Tb=np.zeros((6, 6, 6)))
        Ba = np.dot(Rb, Ba) + Bb
        E += elem.transfer_map.delta_e
    lattice.E = E
    lattice.T_sym = Ta
    lattice.T = unsym_matrix(deepcopy(Ta))
    lattice.R = Ra
    lattice.B = Ba
    return Ra


def trace_z(lattice, obj0, z_array):
    """
    Z-dependent tracer (twiss(z) and particle(z))
    usage: twiss = trace_z(lattice,twiss_0, [1.23, 2.56, ...]) ,
    to calculate Twiss params at 1.23m, 2.56m etc.
    """
    obj_list = []
    i = 0
    elem = lattice.sequence[i]
    L = elem.l
    obj_elem = obj0
    for z in z_array:
        while z > L:
            obj_elem = lattice.sequence[i].transfer_map * obj_elem
            i += 1
            elem = lattice.sequence[i]
            L += elem.l

        obj_z = elem.transfer_map(z - (L - elem.l)) * obj_elem

        obj_list.append(obj_z)
    return obj_list


def trace_obj(lattice, obj, nPoints=None):
    """
    track object through the lattice
    obj must be Twiss or Particle
    """

    if nPoints is None:
        obj_list = [obj]
        for e in lattice.sequence:
            obj = e.transfer_map * obj
            obj.id = e.id
            obj_list.append(obj)
    else:
        z_array = np.linspace(0, lattice.totalLen, nPoints, endpoint=True)
        obj_list = trace_z(lattice, obj, z_array)
    return obj_list


def periodic_twiss(tws, R):
    """
    initial conditions for a periodic Twiss solution
    """
    tws = Twiss(tws)

    cosmx = (R[0, 0] + R[1, 1]) / 2.
    cosmy = (R[2, 2] + R[3, 3]) / 2.

    if abs(cosmx) >= 1 or abs(cosmy) >= 1:
        _logger.warning(" ************ periodic solution does not exist. return None ***********")
        return None
    sinmx = np.sign(R[0, 1]) * np.sqrt(1. - cosmx * cosmx)
    sinmy = np.sign(R[2, 3]) * np.sqrt(1. - cosmy * cosmy)

    tws.beta_x = abs(R[0, 1] / sinmx)
    tws.beta_y = abs(R[2, 3] / sinmy)

    tws.alpha_x = (R[0, 0] - R[1, 1]) / (2. * sinmx)  # X[0,0]

    tws.gamma_x = (1. + tws.alpha_x * tws.alpha_x) / tws.beta_x  # X[1,0]

    tws.alpha_y = (R[2, 2] - R[3, 3]) / (2 * sinmy)  # Y[0,0]
    tws.gamma_y = (1. + tws.alpha_y * tws.alpha_y) / tws.beta_y  # Y[1,0]

    Hx = np.array([[R[0, 0] - 1, R[0, 1]], [R[1, 0], R[1, 1] - 1]])
    Hhx = np.array([[R[0, 5]], [R[1, 5]]])
    hh = np.dot(inv(-Hx), Hhx)
    tws.Dx = hh[0, 0]
    tws.Dxp = hh[1, 0]
    Hy = np.array([[R[2, 2] - 1, R[2, 3]], [R[3, 2], R[3, 3] - 1]])
    Hhy = np.array([[R[2, 5]], [R[3, 5]]])
    hhy = np.dot(inv(-Hy), Hhy)
    tws.Dy = hhy[0, 0]
    tws.Dyp = hhy[1, 0]
    return tws


def twiss(lattice, tws0=None, nPoints=None):
    """
    twiss parameters calculation

    :param lattice: lattice, MagneticLattice() object
    :param tws0: initial twiss parameters, Twiss() object. If None, try to find periodic solution.
    :param nPoints: number of points per cell. If None, then twiss parameters are calculated at the end of each element.
    :return: list of Twiss() objects
    """
    if tws0 is None:
        tws0 = periodic_twiss(tws0, lattice_transfer_map(lattice, energy=0.))

    if tws0.__class__ == Twiss:
        if tws0.beta_x == 0 or tws0.beta_y == 0:
            R = lattice_transfer_map(lattice, tws0.E)
            tws0 = periodic_twiss(tws0, R)
            if tws0 is None:
                _logger.info(' twiss: Twiss: no periodic solution')
                return None
        else:
            tws0.gamma_x = (1. + tws0.alpha_x ** 2) / tws0.beta_x
            tws0.gamma_y = (1. + tws0.alpha_y ** 2) / tws0.beta_y

        twiss_list = trace_obj(lattice, tws0, nPoints)
        return twiss_list
    else:
        _logger.warning(' Twiss: no periodic solution. return None')
        return None


def twiss_fast(lattice, tws0=None):
    """
    twiss parameters calculation

    :param lattice: lattice, MagneticLattice() object
    :param tws0: initial twiss parameters, Twiss() object. If None, try to find periodic solution.
    :param nPoints: number of points per cell. If None, then twiss parameters are calculated at the end of each element.
    :return: list of Twiss() objects
    """
    if tws0 is None:
        tws0 = periodic_twiss(tws0, lattice_transfer_map(lattice, energy=0.))
    if tws0.__class__ == Twiss:
        if tws0.beta_x == 0 or tws0.beta_y == 0:
            R = lattice_transfer_map(lattice, tws0.E)
            tws0 = periodic_twiss(tws0, R)
            if tws0 is None:
                _logger.warning(' twiss_fast: Twiss: no periodic solution')
                return None
        else:
            tws0.gamma_x = (1. + tws0.alpha_x ** 2) / tws0.beta_x
            tws0.gamma_y = (1. + tws0.alpha_y ** 2) / tws0.beta_y

        obj_list = [tws0]
        for e in lattice.fast_seq:
            e.transfer_map.R = lambda x: e.transfer_map._r
            tws0 = e.transfer_map * tws0
            tws0.id = e.id
            obj_list.append(tws0)
        return obj_list
    else:
        _logger.warning(' twiss_fast: Twiss: no periodic solution')
        return None


class ProcessTable:
    def __init__(self, lattice):
        self.proc_list = []
        self.kick_proc_list = []
        self.lat = lattice

    def searching_kick_proc(self, physics_proc, elem1):
        """
        function finds kick physics process. Kick physics process applies kick only once between two elements
        with zero length (e.g. Marker) or at the beginning of the element if it is the same element,
        others physics processes are applied during finite lengths.
        :return:
        """

        if (physics_proc.indx0 == physics_proc.indx1 or
            (physics_proc.indx0 + 1 == physics_proc.indx1 and elem1.l == 0)):
            physics_proc.indx1 = physics_proc.indx0
            physics_proc.s_stop = physics_proc.s_start
            self.kick_proc_list.append(physics_proc)
            if len(self.kick_proc_list) > 1:
                pos = np.array([proc.s_start for proc in self.kick_proc_list])
                indx = np.argsort(pos)
                self.kick_proc_list = [self.kick_proc_list[i] for i in indx]
        _logger_navi.debug(" searching_kick_proc: self.kick_proc_list.append(): " + str([p.__class__.__name__ for p in self.kick_proc_list]))

    def add_physics_proc(self, physics_proc, elem1, elem2):
        physics_proc.start_elem = elem1
        physics_proc.end_elem = elem2
        physics_proc.indx0 = self.lat.sequence.index(elem1)
        physics_proc.indx1 = self.lat.sequence.index(elem2)
        physics_proc.s_start = np.sum(np.array([elem.l for elem in self.lat.sequence[:physics_proc.indx0]]))
        physics_proc.s_stop = np.sum(np.array([elem.l for elem in self.lat.sequence[:physics_proc.indx1]]))
        self.searching_kick_proc(physics_proc, elem1)
        physics_proc.counter = physics_proc.step
        physics_proc.prepare(self.lat)

        _logger_navi.debug(" add_physics_proc: self.proc_list = " + str([p.__class__.__name__ for p in self.proc_list]) + ".append(" + physics_proc.__class__.__name__ + ")" +
                          "; start: " + str(physics_proc.indx0 ) + " stop: " + str(physics_proc.indx1))

        self.proc_list.append(physics_proc)


class Navigator:
    """
    Navigator defines step (dz) of tracking and which physical process will be applied during each step.
    lattice - MagneticLattice
    Attributes:
        unit_step = 1 [m] - unit step for all physics processes
    Methods:
        add_physics_proc(physics_proc, elem1, elem2)
            physics_proc - physics process, can be CSR, SpaceCharge or Wake,
            elem1 and elem2 - first and last elements between which the physics process will be applied.
    """

    def __init__(self, lattice):

        self.lat = lattice
        self.process_table = ProcessTable(self.lat)
        self.ref_process_table = deepcopy(self.process_table)

        self.z0 = 0.  # current position of navigator
        self.n_elem = 0  # current index of the element in lattice
        self.sum_lengths = 0.  # sum_lengths = Sum[lat.sequence[i].l, {i, 0, n_elem-1}]
        self.unit_step = 1  # unit step for physics processes
        self.proc_kick_elems = []
        self.kill_process = False # for case when calculations are needed to terminated e.g. from gui

    def reset_position(self):
        """
        method to reset Navigator position.
        :return:
        """
        _logger_navi.debug(" reset position")
        self.z0 = 0.  # current position of navigator
        self.n_elem = 0  # current index of the element in lattice
        self.sum_lengths = 0.  # sum_lengths = Sum[lat.sequence[i].l, {i, 0, n_elem-1}]
        self.process_table = deepcopy(self.ref_process_table)

    def go_to_start(self):
        self.reset_position()

    def get_phys_procs(self):
        """
        method return list of all physics processes which were added

        :return: list, list of PhysProc(s)
        """
        return self.process_table.proc_list

    def add_physics_proc(self, physics_proc, elem1, elem2):
        """
        Method adds Physics Process.

        :param physics_proc: PhysicsProc, e.g. SpaceCharge, CSR, Wake ...
        :param elem1: the element in the lattice where to start applying the physical process.
        :param elem2: the element in the lattice where to stop applying the physical process,
                        can be the same as starting element.
        :return:
        """
        _logger_navi.debug(" add_physics_proc: phys proc: " + physics_proc.__class__.__name__)
        self.process_table.add_physics_proc(physics_proc, elem1, elem2)
        self.ref_process_table = deepcopy(self.process_table)

    def activate_apertures(self, start=None, stop=None):
        """
        activate apertures if thea exist in the lattice from

        :param start: element,  activate apertures starting form element 'start' element
        :param stop: element, activate apertures up to 'stop' element
        :return:
        """
        id1 = self.lat.sequence.index(start) if start is not None else None
        id2 = self.lat.sequence.index(stop) if stop is not None else None
        for elem in self.lat.sequence[id1:id2]:
            if elem.__class__ is Aperture:
                if elem.type == "rect":
                    ap = RectAperture(xmin=-elem.xmax + elem.dx, xmax=elem.xmax + elem.dx,
                                      ymin=-elem.ymax + elem.dy, ymax=elem.ymax + elem.dy)
                    self.add_physics_proc(ap, elem, elem)
                elif elem.type == "ellipt":
                    ap = EllipticalAperture(xmax=elem.xmax, ymax=elem.ymax,
                                      dx=elem.dx, dy=elem.dy)
                    self.add_physics_proc(ap, elem, elem)

    def check_overjump(self, dz, processes, phys_steps):
        phys_steps_red = phys_steps - dz
        if len(processes) != 0:
            nearest_stop_elem = min([proc.indx1 for proc in processes])
            L_stop = np.sum(np.array([elem.l for elem in self.lat.sequence[:nearest_stop_elem]]))
            if self.z0 + dz > L_stop:
               dz = L_stop - self.z0

            # check if inside step dz there is another phys process

            proc_list_rest = [p for p in self.process_table.proc_list if p not in processes]
            start_pos_of_rest = np.array([proc.s_start for proc in proc_list_rest])
            start_pos = [pos for pos in start_pos_of_rest if self.z0 < pos < self.z0 + dz]
            if len(start_pos) > 0:
                start_pos.sort()
                dz = start_pos[0] - self.z0
                _logger_navi.debug(" check_overjump: there is phys proc inside step -> dz was decreased: dz = " + str(dz))

        phys_steps = phys_steps_red + dz

        # check kick processes
        kick_list = self.process_table.kick_proc_list
        kick_pos = np.array([proc.s_start for proc in kick_list])
        indx = np.argwhere(self.z0 < kick_pos)

        if 0 in kick_pos and self.z0 == 0 and self.n_elem == 0:
            indx0 = np.argwhere(self.z0 == kick_pos)
            indx = np.append(indx0, indx)
        if len(indx) != 0:
            kick_process = np.array([kick_list[i] for i in indx.flatten()])
            for i, proc in enumerate(kick_process):
                L_kick_stop = proc.s_start
                if self.z0 + dz > L_kick_stop:
                    dz = L_kick_stop - self.z0
                    phys_steps = phys_steps_red + dz
                    processes.append(proc)
                    phys_steps = np.append(phys_steps, 0)
                    continue
                elif self.z0 + dz == L_kick_stop:
                    processes.append(proc)
                    phys_steps = np.append(phys_steps, 0)
                else:
                    pass

        return dz, processes, phys_steps

    def get_proc_list(self):
        _logger_navi.debug(" get_proc_list: all phys proc = " + str([p.__class__.__name__ for p in self.process_table.proc_list]))
        proc_list = []
        for p in self.process_table.proc_list:
            if p.indx0 <= self.n_elem < p.indx1:
                proc_list.append(p)
        return proc_list

    def hard_edge_step(self, dz):
        # self.sum_lengths
        elem1 = self.lat.sequence[self.n_elem]
        L = self.sum_lengths + elem1.l
        if self.z0 + dz > L:
            dz = L - self.z0
        return dz

    def check_proc_bounds(self, dz, proc_list, phys_steps, active_process):
        for p in proc_list:
            if dz + self.z0 >= p.s_stop and p not in active_process:
                active_process.append(p)
                phys_steps = np.append(phys_steps, dz)

        return active_process, phys_steps

    def remove_used_processes(self, processes):
        """
        in case physics processes are applied and do not more needed they are removed from table

        :param processes: list of processes are about to apply
        :return: None
        """
        for p in processes:
            if p in self.process_table.kick_proc_list:
                _logger_navi.debug(" Navigator.remove_used_processes: " + p.__class__.__name__)
                self.process_table.kick_proc_list.remove(p)
                self.process_table.proc_list.remove(p)

    def get_next(self):

        proc_list = self.get_proc_list()

        if len(proc_list) > 0:

            counters = np.array([p.counter for p in proc_list])
            step = counters.min()

            inxs = np.where(counters == step)

            processes = [proc_list[i] for i in inxs[0]]

            phys_steps = np.array([p.step for p in processes])*self.unit_step

            for p in proc_list:
                p.counter -= step
                if p.counter == 0:
                    p.counter = p.step

            dz = np.min(phys_steps)
        else:
            processes = proc_list
            n_elems = len(self.lat.sequence)
            if n_elems >= self.n_elem + 1:
                L = np.sum(np.array([elem.l for elem in self.lat.sequence[:self.n_elem + 1]]))
            else:
                L = self.lat.totalLen
            dz = L - self.z0
            phys_steps = np.array([])
        # check if dz overjumps the stop element
        dz, processes, phys_steps = self.check_overjump(dz, processes, phys_steps)
        processes, phys_steps = self.check_proc_bounds(dz, proc_list, phys_steps, processes)

        _logger_navi.debug(" Navigator.get_next: process: " + " ".join([proc.__class__.__name__ for proc in processes]))

        _logger_navi.debug(" Navigator.get_next: navi.z0=" + str(self.z0) + " navi.n_elem=" + str(self.n_elem) + " navi.sum_lengths="
                     + str(self.sum_lengths) + " dz=" + str(dz))

        _logger_navi.debug(" Navigator.get_next: element type=" + self.lat.sequence[self.n_elem].__class__.__name__ + " element name=" +
                     str(self.lat.sequence[self.n_elem].id))

        self.remove_used_processes(processes)

        return dz, processes, phys_steps


def get_map(lattice, dz, navi):
    nelems = len(lattice.sequence)
    TM = []
    i = navi.n_elem
    z1 = navi.z0 + dz
    elem = lattice.sequence[i]
    # navi.sum_lengths = np.sum([elem.l for elem in lattice.sequence[:i]])
    L = navi.sum_lengths + elem.l
    while z1 + 1e-10 > L:

        dl = L - navi.z0
        TM.append(elem.transfer_map(dl))

        navi.z0 = L
        dz -= dl
        if i >= nelems - 1:
            break

        i += 1
        elem = lattice.sequence[i]
        L += elem.l
        #if i in navi.proc_kick_elems:
        #    break
    if abs(dz) > 1e-10:
        TM.append(elem.transfer_map(dz))
    navi.z0 += dz
    navi.sum_lengths = L - elem.l
    navi.n_elem = i
    return TM


def merge_maps(t_maps):
    tm0 = TransferMap()
    t_maps_new = []
    for tm in t_maps:
        if tm.__class__ == TransferMap:
            tm0 = tm * tm0
        else:
            t_maps_new.append(tm0)
            t_maps_new.append(tm)
            tm0 = TransferMap()
    t_maps_new.append(tm0)
    return t_maps_new


'''
returns two solutions for a periodic fodo, given the mean beta
initial betas are at the center of the focusing quad
'''


def fodo_parameters(betaXmean=36.0, L=10.0, verbose=False):
    lquad = 0.001

    kap1 = np.sqrt(1.0 / 2.0 * (
    (betaXmean / L) * (betaXmean / L) + (betaXmean / L) * np.sqrt(-4.0 + (betaXmean / L) * (betaXmean / L))))
    kap2 = np.sqrt(1.0 / 2.0 * (
    (betaXmean / L) * (betaXmean / L) - (betaXmean / L) * np.sqrt(-4.0 + (betaXmean / L) * (betaXmean / L))))

    k = 1.0 / (lquad * L * kap2)

    f = 1.0 / (k * lquad)

    kappa = f / L
    betaMax = np.array(
        (L * kap1 * (kap1 + 1) / np.sqrt(kap1 * kap1 - 1), L * kap2 * (kap2 + 1) / np.sqrt(kap2 * kap2 - 1)))
    betaMin = np.array(
        (L * kap1 * (kap1 - 1) / np.sqrt(kap1 * kap1 - 1), L * kap2 * (kap2 - 1) / np.sqrt(kap2 * kap2 - 1)))
    betaMean = np.array(
        (L * kap2 * kap2 / (np.sqrt(kap2 * kap2 - 1.0)), L * kap1 * kap1 / (np.sqrt(kap1 * kap1 - 1.0))))
    k = np.array((1.0 / (lquad * L * kap1), 1.0 / (lquad * L * kap2)))

    if verbose:
        print('********* calculating fodo parameters *********')
        print('fodo parameters:')
        print('k*l=', k * lquad)
        print('f=', L * kap1, L * kap2)
        print('kap1=', kap1)
        print('kap2=', kap2)
        print('betaMax=', betaMax)
        print('betaMin=', betaMin)
        print('betaMean=', betaMean)
        print('*********                             *********')

    return k * lquad, betaMin, betaMax, betaMean
