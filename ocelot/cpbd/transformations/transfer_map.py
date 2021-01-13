from copy import deepcopy, copy
import logging

import numpy as np

from ocelot.cpbd.beam import Twiss, Particle, ParticleArray

_logger = logging.getLogger(__name__)


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
        self.B_z = lambda z, energy: np.dot((np.eye(6) - self.R_z(z, energy)),
                                            np.array([[self.dx], [0.], [self.dy], [0.], [0.], [0.]]))
        self.B = lambda energy: self.B_z(self.length, energy)
        self._map = None

    @property
    def map(self):
        if not self._map:
            self._map = self.map_function()
        return self._map

    @map.setter
    def map(self, func):
        self._map = func

    def map_function(self, delta_length=None, length=None):
        """
        This function calculate the map function which can be overload if the map function is different to the first order mapping  
        @param delta_length: delta length of the element
        @param length: total length of the element
        """
        return lambda u, energy: self.mul_p_array(u, energy=energy)

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
        """
        Calculates new rpaticles with the first order transformation, overrides the old rpaticles and returns the new rparticles. 
        :param rparticles: Can be a ParticleArray, Particle or a list of Particle object. 
        :param engery:
        :return: Returns the modified rparticles 
        """
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
            _logger.error(
                " TransferMap.__mul__: unknown object in transfer map multiplication: " + str(m.__class__.__name__))
            raise Exception(
                " TransferMap.__mul__: unknown object in transfer map multiplication: " + str(m.__class__.__name__))

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
            p.x, p.px, p.y, p.py, p.tau, p.p = self.map(np.array([[p.x], [p.px], [p.y], [p.py], [p.tau], [p.p]]), p.E)[
                :, 0]
            p.s += self.length
            p.E += self.delta_e

        elif prcl_series.__class__ == list and prcl_series[0].__class__ == Particle:
            # If the energy is not the same (p.E) for all Particles in the list of Particles
            # in that case cycle is applied. For particles with the same energy p.E
            list_e = np.array([p.E for p in prcl_series])
            if False in (list_e[:] == list_e[0]):
                for p in prcl_series:
                    # TODO: Remove this line of code because its doing nothing. What was the idea?
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
            raise Exception(
                " TransferMap.apply(): Unknown type of Particle_series: " + str(prcl_series.__class__.__name))

    def __call__(self, delta_length):
        m = copy(self)
        m.R = lambda energy: m.R_z(delta_length, energy)
        m.B = lambda energy: m.B_z(delta_length, energy)
        m.delta_e = m.delta_e_z(delta_length)
        m.map = m.map_function(delta_length=delta_length, length=self.length)
        m.length = delta_length
        return m

    @classmethod
    def create_from_element(cls, element, params=None):
        return cls()
