from ocelot.cpbd.io import save_particle_array
from ocelot.common.globals import *
import numpy as np
from ocelot.cpbd.beam import Twiss
from scipy import optimize
from ocelot.utils.acc_utils import *
from ocelot.common.logging import *
_logger = logging.getLogger(__name__)


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

    def finalize(self, *args, **kwargs):
        """
        the method is called at the end of tracking

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
        _logger.debug(" SaveBeam applied, dz =", dz)
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

        _logger.debug(" SmoothBeam applied, dz =" + str(dz))
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


class LaserModulator(PhysProc):
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
        self.laser_peak_pos = 0 # relative to the beam center, if 0 laser_peak_pos == mean(p_array.tau()) - laser_peak_pos

    def lambda_ph(self, energy):
        """
        Wavelength of the laser pulse

        :param energy: in [GeV] - beam energy
        :return: wavelength in [m]
        """
        gamma = energy / m_e_GeV
        return self.lperiod / (2 * gamma ** 2) * (1 + self.Ku ** 2 / 2)

    def apply(self, p_array, dz):
        _logger.debug(" LH applied, dz =" + str(dz))
        lbda_ph = self.lambda_ph(p_array.E)
        k_ph = 2 * np.pi / lbda_ph
        pc = np.sqrt(p_array.E ** 2 - m_e_GeV ** 2)

        A = self.dE / (pc) * dz / self.Lu
        dx = p_array.x()[:] - self.x_mean
        dy = p_array.y()[:] - self.y_mean

        tau_mean = np.mean(p_array.tau()) + self.laser_peak_pos
        dtau = p_array.tau()[:] - tau_mean

        p_array.p()[:] += A * np.exp(-dtau**2/(2*self.sigma_l**2))*np.cos(k_ph * p_array.tau()[:]) * np.exp(
            -0.25 * dx ** 2 / self.sigma_x ** 2
            - 0.25 * dy ** 2 / self.sigma_y ** 2)


class LaserHeater(LaserModulator):
    def __init__(self, step=1):
        LaserModulator.__init__(self, step)
        _logger.info("LaserHeater physics process is obsolete. Use 'LaserModulator' instead.")


class Aperture(PhysProc):
    """
    Method to cut beam in longitudinal (by default), horizontal or/and vertical direction
    :param longitudinal: True, cutting in longitudinal direction
    :param vertical: False, cutting in vertical direction
    :param horizontal: False, cutting in horizontal direction

    """
    def __init__(self, step=1):
        PhysProc.__init__(self, step)
        self.longitudinal = True
        self.vertical = False
        self.horizontal = False

        self.zmin = -5   # in simgas
        self.zmax = 5    # in simgas

        self.xmin = -5   # in simgas
        self.xmax = 5    # in simgas

        self.ymin = -5   # in simgas
        self.ymax = 5    # in simgas

    def apply(self, p_array, dz):
        _logger.debug(" Apperture applied")
        if self.longitudinal:
            tau = p_array.tau()[:]
            tau0 = np.mean(tau)
            tau = tau - tau0
            sig = np.std(tau)
            inds = np.argwhere(np.logical_or(tau < sig * self.zmin, tau > sig * self.zmax))
            inds = inds.reshape(inds.shape[0])
            p_array.rparticles = np.delete(p_array.rparticles, inds, 1)
            p_array.q_array = np.delete(p_array.q_array, inds, 0)

        if self.horizontal:
            x = p_array.x()
            x0 = np.mean(x)
            x = x - x0
            sigx = np.std(x)
            inds = np.argwhere(np.logical_or(x < sigx * self.xmin, x > sigx * self.xmax))
            inds = inds.reshape(inds.shape[0])
            p_array.rparticles = np.delete(p_array.rparticles, inds, 1)
            p_array.q_array = np.delete(p_array.q_array, inds, 0)

        if self.vertical:
            y = p_array.y()
            y0 = np.mean(y)
            y = y - y0
            sigy = np.std(y)
            inds = np.argwhere(np.logical_or(y < sigy * self.ymin, y > sigy * self.ymax))
            inds = inds.reshape(inds.shape[0])
            p_array.rparticles = np.delete(p_array.rparticles, inds, 1)
            p_array.q_array = np.delete(p_array.q_array, inds, 0)


class BeamTransform(PhysProc):
    """
    Beam matching
    """
    def __init__(self, tws=None, x_opt=None, y_opt=None):
        """
        :param tws : Twiss object
        :param x_opt (obsolete): [alpha, beta, mu (phase advance)]
        :param y_opt (obsolete): [alpha, beta, mu (phase advance)]
        """
        PhysProc.__init__(self)
        self.bounds = [-5, 5]  # [start, stop] in sigmas
        self.tws = tws       # Twiss
        self.x_opt = x_opt   # [alpha, beta, mu (phase advance)]
        self.y_opt = y_opt   # [alpha, beta, mu (phase advance)]
        self.step = 1
        self.remove_offsets = True

    @property
    def twiss(self):
        if self.tws == None:
            _logger.warning("BeamTransform: x_opt and y_opt are obsolete, use Twiss")
            tws = Twiss()
            tws.alpha_x, tws.beta_x, tws.mux = self.x_opt
            tws.alpha_y, tws.beta_y, tws.muy = self.y_opt
        else:
            tws = self.tws
        return tws

    def apply(self, p_array, dz):
        _logger.debug("BeamTransform: apply")
        self.x_opt = [self.twiss.alpha_x, self.twiss.beta_x, self.twiss.mux]
        self.y_opt = [self.twiss.alpha_y, self.twiss.beta_y, self.twiss.muy]
        self.beam_matching(p_array.rparticles, self.bounds, self.x_opt, self.y_opt)

    def beam_matching(self, particles, bounds, x_opt, y_opt):
        #the beam is centered in the phase space
        pd = np.zeros((int(particles.size / 6), 6))
        dx = 0
        dxp = 0
        dy = 0
        dyp = 0
        if self.remove_offsets:
            dx = np.mean(particles[0])
            dxp = np.mean(particles[1])
            dy = np.mean(particles[2])
            dyp = np.mean(particles[3])

        pd[:, 0] = particles[0] - dx
        pd[:, 1] = particles[1] - dxp
        pd[:, 2] = particles[2] - dy
        pd[:, 3] = particles[3] - dyp
        pd[:, 4] = particles[4]
        pd[:, 5] = particles[5]

        z0 = np.mean(pd[:, 4])
        sig0 = np.std(pd[:, 4])
        inds = np.argwhere((z0 + sig0 * bounds[0] <= pd[:, 4]) * (pd[:, 4] <= z0 + sig0 * bounds[1]))

        mx, mxs, mxx, mxxs, mxsxs, emitx0 = self.moments(pd[inds, 0], pd[inds, 1])
        beta = mxx / emitx0
        alpha = -mxxs / emitx0
        M = self.m_from_twiss([alpha, beta, 0], x_opt)



        particles[0] = M[0, 0] * pd[:, 0] + M[0, 1] * pd[:, 1]
        particles[1] = M[1, 0] * pd[:, 0] + M[1, 1] * pd[:, 1]
        [mx, mxs, mxx, mxxs, mxsxs, emitx0] = self.moments(pd[inds, 2], pd[inds, 3])
        beta = mxx / emitx0
        alpha = -mxxs / emitx0
        M = self.m_from_twiss([alpha, beta, 0], y_opt)
        particles[2] = M[0, 0] * pd[:, 2] + M[0, 1] * pd[:, 3]
        particles[3] = M[1, 0] * pd[:, 2] + M[1, 1] * pd[:, 3]
        return particles

    def moments(self, x, y, cut=0):
        n = len(x)
        inds = np.arange(n)
        mx = np.mean(x)
        my = np.mean(y)
        x = x - mx
        y = y - my
        x2 = x * x
        mxx = np.sum(x2) / n
        y2 = y * y
        myy = np.sum(y2) / n
        xy = x * y
        mxy = np.sum(xy) / n

        emitt = np.sqrt(mxx * myy - mxy * mxy)

        if cut > 0:
            inds = []
            beta = mxx / emitt
            gamma = myy / emitt
            alpha = mxy / emitt
            emittp = gamma * x2 + 2. * alpha * xy + beta * y2
            inds0 = np.argsort(emittp)
            n1 = np.round(n * (100 - cut) / 100)
            inds = inds0[0:n1]
            mx = np.mean(x[inds])
            my = np.mean(y[inds])
            x1 = x[inds] - mx
            y1 = y[inds] - my
            mxx = np.sum(x1 * x1) / n1
            myy = np.sum(y1 * y1) / n1
            mxy = np.sum(x1 * y1) / n1
            emitt = np.sqrt(mxx * myy - mxy * mxy)
        return mx, my, mxx, mxy, myy, emitt

    def m_from_twiss(self, Tw1, Tw2):
        # Transport matrix M for two sets of Twiss parameters (alpha,beta,psi)
        b1 = Tw1[1]
        a1 = Tw1[0]
        psi1 = Tw1[2]
        b2 = Tw2[1]
        a2 = Tw2[0]
        psi2 = Tw2[2]

        psi = psi2-psi1
        cosp = np.cos(psi)
        sinp = np.sin(psi)
        M = np.zeros((2, 2))
        M[0, 0] = np.sqrt(b2/b1)*(cosp+a1*sinp)
        M[0, 1] = np.sqrt(b2*b1)*sinp
        M[1, 0] = ((a1-a2)*cosp-(1+a1*a2)*sinp)/np.sqrt(b2*b1)
        M[1, 1] = np.sqrt(b1/b2)*(cosp-a2*sinp)
        return M


class SpontanRadEffects(PhysProc):
    """
    Effects of the spontaneous radiation:
    energy loss and quantum diffusion
    """
    def __init__(self, K=0.0, lperiod=0.0, type="planar"):
        """

        :param Kx: Undulator deflection parameter
        :param lperiod: undulator period in [m]
        :param type:
        """
        PhysProc.__init__(self)
        self.K = K
        self.lperiod = lperiod
        self.type = type
        self.energy_loss = True
        self.quant_diff = True
        self.filling_coeff = 1.0

    def apply(self, p_array, dz):
        _logger.debug("BeamTransform: apply")
        mean_p = np.mean(p_array.p())
        energy = p_array.E*(1 + mean_p)

        if self.quant_diff:
            sigma_Eq = self.sigma_gamma_quant(energy, dz)
            p_array.p()[:] += sigma_Eq * np.random.randn(p_array.n)*self.filling_coeff

        if self.energy_loss:
            dE = self.energy_loss_und(energy, dz)
            p_array.p()[:] -= dE/energy*self.filling_coeff

    def energy_loss_und(self, energy, dz):
        k = 4. * np.pi * np.pi / 3. * ro_e / m_e_GeV
        U = k * energy ** 2 * self.K ** 2 * dz / self.lperiod ** 2
        return U

    def sigma_gamma_quant(self, energy, dz):
        """
        rate of energy diffusion

        :param energy: electron beam energy
        :param Kx: undulator parameter
        :param lperiod: undulator period
        :param dz: length of the
        :return: sigma_gamma/gamma
        """
        gamma = energy / m_e_GeV
        k = 2*np.pi/self.lperiod

        lambda_compt = h_eV_s/m_e_eV*speed_of_light # m
        lambda_compt_r = lambda_compt / 2. / pi
        if self.type == "helical":
            f = lambda K: 1.42*K + 1./(1 + 1.5*K + 0.95*K*K)
        else:
            f = lambda K: 0.6*K + 1. / (2 + 2.66 * K + 0.8 * K ** 2)

        delta_Eq2 = 14/15.*lambda_compt_r*ro_e*gamma**4*k**3*self.K**2*f(self.K)*dz
        sigma_Eq = np.sqrt(delta_Eq2 / (gamma * gamma))
        return sigma_Eq


class BeamAnalysis(PhysProc):
    def __init__(self, filename):
        PhysProc.__init__(self)
        self.filename = filename
        self.lambda_mod = 1e-6
        self.nlambdas = 4 # +- nlambda for analysis
        self.p = []
        self.phi = []
        self.s = []
        self.energy = []
        self.bunching = []

    def apply(self, p_array, dz):
        _logger.debug(" BeamAnalysis applied, dz =", dz)

        def test_func(x, a, phi, delta, g):
            return a * np.sin(2*np.pi/self.lambda_mod * x + phi) + delta + g*x

        bunch_c = np.mean(p_array.tau())

        slice_min = bunch_c - self.lambda_mod*self.nlambdas
        slice_max = bunch_c + self.lambda_mod*self.nlambdas

        indx = np.where(np.logical_and(np.greater_equal(p_array.tau(), slice_min),
                                       np.less(p_array.tau(), slice_max)))[0]
        p = p_array.p()[indx]
        tau = p_array.tau()[indx]


        params, params_covariance = optimize.curve_fit(test_func, tau, p, p0=[0.001, 0, 0, 0])

        b = slice_bunching(tau, charge=np.sum(p_array.q_array[indx]), lambda_mod=self.lambda_mod,
                           smooth_sigma=self.lambda_mod/5)

        self.p.append(params[0])
        self.phi.append(params[1])
        self.s.append(np.copy(p_array.s))
        self.energy.append(np.copy(p_array.E))
        self.bunching.append(b)

    def finalize(self):
        data = np.array([np.array(self.s), np.array(self.p), np.array(self.phi), np.array(self.energy),
                        np.array(self.bunching)])
        np.savetxt(self.filename, data)
