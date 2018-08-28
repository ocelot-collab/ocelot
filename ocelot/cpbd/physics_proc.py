from ocelot.cpbd.io import save_particle_array
from ocelot.common.globals import *
import numpy as np
import logging
logger = logging.getLogger(__name__)


class PhysProc:
    """
    Parent class for all Physics processes

    :method prepare(self, lat): - the method is called at the moment of Physics Process addition to Navigator class.
    :method apply(self, p_array, dz): - the method is called on every step.
    :attribute step: - number of steps in [Navigator.unit_step] self.step*Navigator.unit_step = [m]
    :attribute indx0: - number of start element in lattice.sequence
    :attribute indx1: - number of stop element in lattice.sequence
    """
    def __init__(self, step=1):
        self.step = step
        self.energy = None
        self.indx0 = None
        self.indx1 = None

    def prepare(self, lat):
        """
        method is called at the moment of Physics Process addition to Navigator class.

        :param lat:
        :return:
        """
        pass

    def apply(self, p_array, dz):
        """
        the method is called on every step.

        :param p_array:
        :param dz:
        :return:
        """
        pass


class EmptyProc(PhysProc):
    def __init__(self, step=1):
        PhysProc.__init__(self, step)
        self.energy = None
        self.pict_debug = True
        self.traj_step = 0.0002


class SaveBeam(PhysProc):
    def __init__(self, filename):
        PhysProc.__init__(self)
        self.energy = None
        self.filename = filename

    def apply(self, p_array, dz):
        logger.debug(" SaveBeam applied, dz =", dz)
        save_particle_array(filename=self.filename, p_array=p_array)


class SmoothBeam(PhysProc):
    """
    Physics Process for the beam smoothing. Can be applied when number of particles is not enough.

    :atribuet mslice: number of particles in the slice

    Examples
    --------
    # lat is the MagneticLattice
    navi = Navigator(lat)

    smooth = SmoothBeam()
    smooth.mslice = 10000

    navi.add_physics_process(smooth, start=elem, stop=elem)
    # elem is the lattice element where you want to apply smoothing

    """
    def __init__(self):
        PhysProc.__init__(self)
        self.mslice = 1000

    def apply(self, p_array, dz):
        """
        the method is called on every step.

        :param p_array:
        :param dz:
        :return:
        """

        logger.debug(" SmoothBeam applied, dz =" + str(dz))
        def myfunc(x, A):
            if x < 2 * A:
                y = x - x * x / (4 * A)
            else:
                y = A
            return y

        #Zin = np.copy(p_array.tau())
        Zin = p_array.tau()
        inds = np.argsort(Zin, axis=0)
        Zout = np.copy(Zin[inds])
        N = Zin.shape[0]
        S = np.zeros(N + 1)
        S[N] = 0
        S[0] = 0
        for i in range(N):
            S[i + 1] = S[i] + Zout[i]
        Zout2 = np.zeros(N)
        Zout2[N - 1] = Zout[N - 1]
        Zout2[0] = Zout[0]
        for i in range(1, N - 1):
            m = min(i, N - i + 1)
            m = np.int(np.floor(myfunc(0.5 * m, 0.5 * self.mslice) + 0.500001))
            #print(m)
            Zout2[i] = (S[i + m + 1] - S[i - m]) / (2 * m + 1)
        #Zout[inds] = Zout2
        p_array.tau()[inds] = Zout2



class LaserHeater(PhysProc):
    def __init__(self, step=1):
        PhysProc.__init__(self, step)
        # amplitude of energy modulation on axis
        self.dE = 12500e-9  # GeV
        self.Ku = 1.294  # undulator parameter
        self.Lu = 0.8  # [m] - undulator length
        self.lperiod = 0.074  # [m] - undulator period length
        self.sigma_l = 300e-6  # [m]
        self.sigma_x = self.sigma_l
        self.sigma_y = self.sigma_l
        self.x_mean = 0
        self.y_mean = 0

    def apply(self, p_array, dz):
        logger.debug(" LH applied, dz =", dz)
        gamma = p_array.E / m_e_GeV
        lbda_ph = self.lperiod / (2 * gamma ** 2) * (1 + self.Ku ** 2 / 2)
        k_ph = 2 * np.pi / lbda_ph
        pc = np.sqrt(p_array.E ** 2 - m_e_GeV ** 2)

        A = self.dE / (pc) * dz / self.Lu
        dx = p_array.x()[:] - self.x_mean
        dy = p_array.y()[:] - self.y_mean
        p_array.p()[:] += A * np.cos(k_ph * p_array.tau()[:]) * np.exp(
            -0.25 * dx ** 2 / self.sigma_x ** 2
            - 0.25 * dy ** 2 / self.sigma_y ** 2)

class Apperture(PhysProc):
    def __init__(self, step=1):
        PhysProc.__init__(self, step)
        self.zmin = -5
        self.zmax = 5

    def apply(self, p_array, dz):
        logger.debug(" Apperture applied")
        tau = p_array.tau()[:]
        tau0 = np.mean(tau)
        tau = tau - tau0
        sig = np.std(tau)
        inds = np.argwhere(np.logical_or(tau < sig * self.zmin, tau > sig * self.zmax))
        inds = inds.reshape(inds.shape[0])
        p_array.rparticles = np.delete(p_array.rparticles, inds, 1)
        p_array.q_array = np.delete(p_array.q_array, inds, 0)