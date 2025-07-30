from ocelot.cpbd.io import save_particle_array
from ocelot.common.globals import h_eV_s, m_e_eV, m_e_GeV, ro_e, speed_of_light, q_e
from ocelot.cpbd.beam import Twiss, beam_matching, global_slice_analysis, s_to_cur, get_envelope
from ocelot.utils.acc_utils import slice_bunching
from ocelot.common.ocelog import *
from ocelot.cpbd.beam import ParticleArray

import numpy as np
from scipy import optimize

_logger = logging.getLogger(__name__)


class PhysProc:
    """
    Parent class for all Physics processes

    :method prepare(self, lat): - the method is called at the moment of Physics Process addition to Navigator class.
    :method apply(self, p_array, dz): - the method is called on every step.
    :attribute step: - number of steps in [Navigator.unit_step] self.step*Navigator.unit_step = [m]
    :attribute indx0: - number of start element in lattice.sequence - assigned in navigator.add_physics_proc()
    :attribute indx1: - number of stop element in lattice.sequence - assigned in navigator.add_physics_proc()
    :attribute s_start: - position of start element in lattice - assigned in navigator.add_physics_proc()
    :attribute s_stop: - position of stop element in lattice.sequence - assigned in navigator.add_physics_proc()
    :attribute start_elem: -  start element in lattice - assigned in navigator.add_physics_proc()
    :attribute end_elem: -  stop element in lattice.sequence - assigned in navigator.add_physics_proc()
    :attribute z0: - current position of navigator - assigned in track.track() before p.apply()
    """

    def __init__(self, step=1):
        self.step = step
        self.energy = None
        self.indx0 = None
        self.indx1 = None
        self.s_start = None
        self.s_stop = None
        self.start_elem = None
        self.end_elem = None
        self.z0 = None

    def check_step(self):
        if not isinstance(self.step, (int, float)) and float(self.step).is_integer():
            raise ValueError(f'step must be an integer number, instead {self.step}')

    def prepare(self, lat):
        """
        method is called at the moment of Physics Process addition to Navigator class.

        :param lat:
        :return:
        """
        self.check_step()

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
        _logger.debug(" SaveBeam applied, dz =" + str(dz))
        save_particle_array(filename=self.filename, p_array=p_array)


class CopyBeam(PhysProc):
    """Physics process that copies the ParticleArray instance when applied.  Makes
    most sense to be attached to zero-length elements (e.g. Marker instances)."""
    def __init__(self, name: str = ""):
        super().__init__()
        self.name = name
        self.parray = None

    def apply(self, parray: ParticleArray, dz: float) -> None:
        """Copy the given particle array to self"""
        self.parray = parray.copy()

    def __repr__(self) -> str:
        return f"<CopyBeam: {self.name}, at={hex(id(self))}>"


class SmoothBeam(PhysProc):
    """
    Physics Process for the beam smoothing. Can be applied when number of particles is not enough.

    :attribute mslice: number of particles in the slice

    Examples
    --------
    # lat is the MagneticLattice
    navi = Navigator(lat)

    smooth = SmoothBeam()
    smooth.mslice = 10000

    navi.add_physics_process(smooth, start=elem, stop=elem)
    # elem is the lattice element where you want to apply smoothing

    """

    def __init__(self, mslice=1000):
        PhysProc.__init__(self)
        self.mslice = mslice

    def apply(self, p_array, dz):
        """
        the method is called on every step.

        :param p_array:
        :param dz:
        :return:
        """

        _logger.debug(" SmoothBeam applied, dz =" + str(dz))

        def myfunc(x, A):
            y = np.where(x < 2 * A, x - x * x / (4 * A), A)
            return y
        Zin = p_array.tau()
        inds = np.argsort(Zin, axis=0)
        Zout = np.sort(Zin, axis=0)
        N = Zin.shape[0]
        S = np.zeros(N + 1)
        S[1:] = np.cumsum(Zout)
        Zout2 = np.zeros(N)
        Zout2[N - 1] = Zout[N - 1]
        Zout2[0] = Zout[0]

        i = np.arange(1, N - 1)
        m = np.minimum(i, N - i + 1)
        m = np.floor(myfunc(0.5 * m, 0.5 * self.mslice) + 0.500001).astype(int)
        Zout2[i] = (S[i + m + 1] - S[i - m]) / (2 * m + 1)
        #Zout[inds] = Zout2
        p_array.tau()[inds] = Zout2

    def __repr__(self) -> str:
        cname = type(self).__name__
        mslice = self.mslice
        return f"<{cname}: {mslice=}>"


class LaserModulator(PhysProc):
    def __init__(self, **kwargs):
        # Extract 'step' if provided, otherwise default to 1
        step = kwargs.pop('step', 1)
        super().__init__(step)

        # Pull out known parameters with defaults
        self.dE           = kwargs.pop('dE',           12500e-9)  # GeV
        self.Ku           = kwargs.pop('Ku',           1.294)     # Undulator parameter
        self.Lu           = kwargs.pop('Lu',           0.8)       # [m] - Undulator length
        self.lperiod      = kwargs.pop('lperiod',      0.074)     # [m] - Undulator period length
        self.sigma_l      = kwargs.pop('sigma_l',      300e-6)    # [m]
        self.sigma_x      = kwargs.pop('sigma_x',      self.sigma_l)
        self.sigma_y      = kwargs.pop('sigma_y',      self.sigma_l)
        self.x_mean       = kwargs.pop('x_mean',       0)
        self.y_mean       = kwargs.pop('y_mean',       0)
        self.z_waist      = kwargs.pop('z_waist',      None)       # Center of the undulator
        self.include_r56  = kwargs.pop('include_r56',  False)
        self.laser_peak_pos = kwargs.pop('laser_peak_pos', 0)
        # relative to the beam center; if 0, laser_peak_pos == mean(p_array.tau()) - laser_peak_pos

        # If there are any additional kwargs left, you can set them as attributes:
        for key, value in kwargs.items():
            setattr(self, key, value)

    def lambda_ph(self, energy):
        """
        Wavelength of the laser pulse

        :param energy: in [GeV] - beam energy
        :return: wavelength in [m]
        """
        gamma = energy / m_e_GeV
        return self.lperiod / (2 * gamma ** 2) * (1 + self.Ku ** 2 / 2)

    def r56(self, energy):
        """
        Method calculate R56 of the undulator

        :param energy: in [GeV] - beam energy
        :return: R56 in [m]
        """
        gamma = energy / m_e_GeV
        beta = 1 / np.sqrt(1.0 - 1.0 / (gamma * gamma))
        r56 = - self.Lu / (gamma * beta) ** 2 * (1 + 0.5 * (self.Ku * beta) ** 2)
        return r56

    def apply(self, p_array, dz):
        _logger.debug(" LaserModulator applied, dz =" + str(dz))

        L = self.s_stop - self.s_start
        if L == 0:
            _logger.warning(" LaserModulator is not applied, undulator length =" + str(L))
            return
        else:
            if np.abs(self.Lu - L) > 1e-5:
                _logger.warning(
                    " LaserModulator: undulator length ({}) is not equal the distance between Markers ({}). "
                    "Distance between Markers is used".format(self.Lu, L))

        lbda_ph = self.lambda_ph(p_array.E)
        k_ph = 2 * np.pi / lbda_ph
        pc = np.sqrt(p_array.E ** 2 - m_e_GeV ** 2)

        A = self.dE / (pc) * dz / L
        dx = p_array.x()[:] - self.x_mean
        dy = p_array.y()[:] - self.y_mean

        tau_mean = np.mean(p_array.tau()) + self.laser_peak_pos
        dtau = p_array.tau()[:] - tau_mean

        if self.z_waist is None:
            p_array.p()[:] += A * np.exp(-dtau ** 2 / (2 * self.sigma_l ** 2)) * np.cos(
                k_ph * p_array.tau()[:]) * np.exp(
                -0.25 * dx ** 2 / self.sigma_x ** 2 - 0.25 * dy ** 2 / self.sigma_y ** 2)
            if self.include_r56:
                p_array.tau()[:] += self.r56(p_array.E) * p_array.p()[:] * dz / L

        else:

            z_waist_abs = self.s_start + L / 2 + self.z_waist
            z = self.z0 - p_array.tau()[:] - z_waist_abs
            zRx = 2 * k_ph * self.sigma_x ** 2
            zRy = 2 * k_ph * self.sigma_y ** 2
            p = 2 * complex(0, 1) * z / k_ph
            N = np.exp(-dx ** 2 / (4 * self.sigma_x ** 2 - p) - dy ** 2 / (4 * self.sigma_y ** 2 - p))
            D = np.sqrt((1 - complex(0, 1) * z / zRx) * (1 - complex(0, 1) * z / zRy))
            V = np.real(N / D * np.exp(- complex(0, 1) * k_ph * p_array.tau()[:]))
            p_array.p()[:] += A * V * np.exp(-dtau ** 2 / (2 * self.sigma_l ** 2))


class LaserHeater(LaserModulator):
    def __init__(self, step=1):
        LaserModulator.__init__(self, step)
        _logger.info("LaserHeater physics process is obsolete. Use 'LaserModulator' instead.")


class PhaseSpaceAperture(PhysProc):
    """
    Method to cut beam in longitudinal (by default), horizontal or/and vertical direction

    :param longitudinal: True, cutting in longitudinal direction
    :param vertical: False, cutting in vertical direction
    :param horizontal: False, cutting in horizontal direction
    :param taumin: -5 longitudinal plane in [rms] from center of mass
    :param taumax: 5 longitudinal plane in [rms] from center of mass

    :param xmin: -5 horizontal plane in [rms] from center of mass
    :param xmax: 5 horizontal plane in [rms] from center of mass

    :param ymin: -5 vertical plane in [rms] from center of mass
    :param ymax: 5 vertical plane in [rms] from center of mass
    """

    def __init__(self, step=1, **kwargs):
        PhysProc.__init__(self, step)
        self.longitudinal = kwargs.get("longitudinal", True)
        self.vertical = kwargs.get("vertical", False)
        self.horizontal = kwargs.get("horizontal", False)

        self.taumin = kwargs.get("taumin", -5)  # in simgas
        self.taumax = kwargs.get("taumax", 5)  # in simgas

        self.xmin = kwargs.get("xmin", -5)  # in simgas
        self.xmax = kwargs.get("xmax", 5)  # in simgas

        self.ymin = kwargs.get("ymin", -5)  # in simgas
        self.ymax = kwargs.get("ymax", 5)  # in simgas

    def apply(self, p_array, dz):
        _logger.debug(" Aperture applied")
        if self.longitudinal:
            tau = p_array.tau()[:]
            tau0 = np.mean(tau)
            tau = tau - tau0
            sig = np.std(tau)
            inds = np.argwhere(np.logical_or(tau < sig * self.taumin, tau > sig * self.taumax))
            inds = inds.reshape(inds.shape[0])
            p_array.delete_particles(inds)

        if self.horizontal:
            x = p_array.x()
            x0 = np.mean(x)
            x = x - x0
            sigx = np.std(x)
            inds = np.argwhere(np.logical_or(x < sigx * self.xmin, x > sigx * self.xmax))
            inds = inds.reshape(inds.shape[0])
            p_array.delete_particles(inds)

        if self.vertical:
            y = p_array.y()
            y0 = np.mean(y)
            y = y - y0
            sigy = np.std(y)
            inds = np.argwhere(np.logical_or(y < sigy * self.ymin, y > sigy * self.ymax))
            inds = inds.reshape(inds.shape[0])
            p_array.delete_particles(inds)


class RectAperture(PhysProc):
    """
    Method to cut beam in horizontal or/and vertical direction

    :param xmin: -np.inf horizontal plane in [m]
    :param xmax: np.inf horizontal plane in [m]

    :param ymin: -np.inf vertical plane in [m]
    :param ymax: np.inf vertical plane in [m]
    """

    def __init__(self, xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf, step=1):
        PhysProc.__init__(self, step)
        self.xmin = xmin  # in m
        self.xmax = xmax  # in m

        self.ymin = ymin  # in m
        self.ymax = ymax  # in m

    def apply(self, p_array, dz):
        _logger.debug(" RectAperture applied")

        x = p_array.x()
        inds = np.argwhere(np.logical_or(x < self.xmin, x > self.xmax))
        inds = inds.reshape(inds.shape[0])
        p_array.delete_particles(inds)

        y = p_array.y()
        inds = np.argwhere(np.logical_or(y < self.ymin, y > self.ymax))
        inds = inds.reshape(inds.shape[0])
        p_array.delete_particles(inds)


class EllipticalAperture(PhysProc):
    """
    Method to delete particles outside an ellipse with semi-axes xmax, ymax centered at dx, dy.

    :param xmax: horizontal semi-axis in [m]. Default np.inf
    :param ymax: vertical semi-axis in [m]. Default None, then ymax equals xmax (i.e. circular aperture).
    :param dx: offset in the horizontal axis in [m]. Default 0.0.
    :param dy: offset in the vertical axis in [m]. Default 0.0.
    """

    def __init__(self, xmax=np.inf, ymax=None, dx=0.0, dy=0.0, step=1):
        PhysProc.__init__(self, step)
        self.xmax = xmax
        self.ymax = (ymax if not ymax is None else xmax)
        self.dx = dx
        self.dy = dy

    def apply(self, p_array, dz):
        x = p_array.x()
        y = p_array.y()
        inds = np.argwhere((x - self.dx) ** 2 / self.xmax ** 2 + (y - self.dy) ** 2 / self.ymax ** 2 > 1.0)
        inds = inds.reshape(inds.shape[0])
        p_array.delete_particles(inds)

        
class BeamTransform(PhysProc):
    """
    Beam matching
    """

    def __init__(self, tws=None, x_opt=None, y_opt=None, **kw):
        """
        :param tws : Twiss object
        :param x_opt (obsolete): [alpha, beta, mu (phase advance)]
        :param y_opt (obsolete): [alpha, beta, mu (phase advance)]
        :param remove_offsets: True, before apply matching remove offsets from the beam in all planes
        :param bounds: [-5, 5] in tau-sigmas. Twiss parameters will be calculated for that part of the beam
        :param slice: None, if "Imax" or "Emax" beam matched to that slice and 'bound' param is ignored
        """
        PhysProc.__init__(self)
        self.tws = tws      # Twiss
        self.x_opt = x_opt  # [alpha, beta, mu (phase advance)] - obsolete
        self.y_opt = y_opt  # [alpha, beta, mu (phase advance)] - obsolete
        self.step = 1
        self.remove_offsets = kw.get("remove_offsets", True)
        self.bounds = kw.get("bounds", [-5, 5])  # [start, stop] in sigmas
        self.slice = kw.get("slice", None)

    @property
    def twiss(self):
        if self.tws is None:
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
        beam_matching(p_array, self.bounds, self.x_opt, self.y_opt, self.remove_offsets, self.slice)


class SlottedFoil(PhysProc):
    """
    Class to simulate a slotted foil

    :param dx: thickness of foil [um]
    :param X0: radiation length of the foil material in [cm]

    :param xmin: -np.inf left position of the foil slot [m]
    :param xmax: np.inf right position of the foil slot [m]

    :param ymin: -np.inf lower position of the foil slot [m]
    :param ymax: np.inf upper position of the foil slot [m]
    """
    def __init__(self, dx, X0, xmin=-np.inf, xmax=np.inf, ymin=-np.inf, ymax=np.inf, step=1):
        PhysProc.__init__(self, step)
        self.xmin = xmin  # in m
        self.xmax = xmax  # in m

        self.ymin = ymin  # in m
        self.ymax = ymax  # in m
        self.z = 1        # charge number of the incident particle
        self.dx = dx      # thickness of the foil in [um]
        self.X0 = X0      # radiation length in [cm]

    def scattered_particles(self, p_array):
        x = p_array.x()
        inds = np.argwhere(np.logical_or(x < self.xmin, x > self.xmax))
        inds_x = inds.reshape(inds.shape[0])

        y = p_array.y()
        inds = np.argwhere(np.logical_or(y < self.ymin, y > self.ymax))
        inds_y = inds.reshape(inds.shape[0])
        indeces = np.append(inds_x, inds_y)
        return indeces

    def get_scattering_angle(self, p_array):
        """
        formula from The Review of Particle Physics https://pdg.lbl.gov

        :param p_array: ParticleArray
        :return: theta - rms scattering angle
        """
        p = np.sqrt(p_array.E**2 - m_e_GeV**2) * 1000  # MeV/c momentum of the particles
        gamma = p_array.E / m_e_GeV
        igamma2 = 0.
        if gamma != 0:
            igamma2 = 1. / (gamma * gamma)
        beta = np.sqrt(1. - igamma2)
        theta_rms = 13.6 / (beta * p) * self.z * np.sqrt(self.dx * 1e-4 / self.X0) * (1 + 0.038 * np.log(self.dx*1e-4 / self.X0))
        return theta_rms

    def apply(self, p_array, dz):
        _logger.debug(" SlottedFoil applied")

        theta_rms = self.get_scattering_angle(p_array)
        indeces = self.scattered_particles(p_array)
        n = len(indeces)
        thetas = np.random.normal(0, theta_rms, n)
        alphas = np.random.uniform(0, 2*np.pi, n)
        p_array.px()[indeces] += np.cos(alphas) * thetas
        p_array.py()[indeces] += np.sin(alphas) * thetas


class SpontanRadEffects(PhysProc):
    """
    Effects of the spontaneous radiation:
    energy loss and quantum diffusion
    """

    def __init__(self, K=0.0, lperiod=0.0, type="planar", **kwargs):
        """
        Initialize spontaneous radiation effects.

        :param K: Undulator deflection parameter
        :param lperiod: Undulator period in meters
        :param type: "planar"/"helical" undulator or "dipole"
        :param kwargs: Additional keyword arguments for customization
        """
        super().__init__(**kwargs)  # Pass any extra arguments to the parent class

        # Explicitly defined parameters
        self.K = K
        self.lperiod = lperiod
        self.type = type

        # Optional parameters with default values, can be overridden by kwargs
        self.energy_loss = kwargs.get("energy_loss", True)
        self.quant_diff = kwargs.get("quant_diff", True)
        self.filling_coeff = kwargs.get("filling_coeff", 1.0)
        self.radius = kwargs.get("radius", np.inf)

    def apply(self, p_array, dz):
        _logger.debug("SpontanRadEffects: apply")
        mean_p = np.mean(p_array.p())
        energy = p_array.E * (1 + mean_p)
        gamma = energy / m_e_GeV

        if self.type == "dipole":
            self.K = 100.0  # awake the asymptotic for K>>1
            self.lperiod = 2.0 * np.pi * np.abs(self.radius) / gamma * self.K

        if self.quant_diff:
            sigma_Eq = self.sigma_gamma_quant(energy, dz, self.K, self.lperiod, self.type)
            p_array.p()[:] += sigma_Eq * np.random.randn(p_array.n) * self.filling_coeff

        if self.energy_loss:
            dE = self.energy_loss_und(energy, dz)
            p_array.p()[:] -= dE / energy * self.filling_coeff

    def energy_loss_und(self, energy, dz):
        k = 4. * np.pi * np.pi / 3. * ro_e / m_e_GeV
        U = k * energy ** 2 * self.K ** 2 * dz / self.lperiod ** 2
        return U

    @staticmethod
    def sigma_gamma_quant(energy, dz, K, lperiod, type="planar"):
        """
        rate of energy diffusion

        :param energy: electron beam energy [GeV]
        :param dz: length of the undulator [m]
        :param Kx: undulator parameter
        :param lperiod: undulator period [m]
        :param type: str, undulator type "planar" or "helical"
        :return: sigma_gamma/gamma
        """
        gamma = energy / m_e_GeV
        k = 2 * np.pi / lperiod

        lambda_compt = h_eV_s / m_e_eV * speed_of_light  # m
        lambda_compt_r = lambda_compt / 2. / np.pi
        if type == "helical":
            f = lambda K: 1.42 * K + 1. / (1 + 1.5 * K + 0.95 * K * K)
        else:
            f = lambda K: 0.6 * K + 1. / (2 + 2.66 * K + 0.8 * K ** 2)

        delta_Eq2 = 14 / 15. * lambda_compt_r * ro_e * gamma ** 4 * k ** 3 * K ** 2 * f(K) * dz
        sigma_Eq = np.sqrt(delta_Eq2 / (gamma * gamma))
        return sigma_Eq


class BeamAnalysis(PhysProc):
    def __init__(self, filename):
        PhysProc.__init__(self)
        self.filename = filename
        self.lambda_mod = 1e-6
        self.nlambdas = 4  # +- nlambda for analysis
        self.p = []
        self.phi = []
        self.s = []
        self.energy = []
        self.bunching = []
        self.sigma_x = []
        self.sigma_y = []
        # self.sigma_px = []
        # self.sigma_py = []

    def apply(self, p_array, dz):
        _logger.debug(" BeamAnalysis applied, dz =" + str(dz))

        def test_func(x, a, phi, delta, g):
            return a * np.sin(2 * np.pi / self.lambda_mod * x + phi) + delta + g * x

        bunch_c = np.mean(p_array.tau())

        slice_min = bunch_c - self.lambda_mod * self.nlambdas
        slice_max = bunch_c + self.lambda_mod * self.nlambdas

        indx = np.where(np.logical_and(np.greater_equal(p_array.tau(), slice_min),
                                       np.less(p_array.tau(), slice_max)))[0]
        p = p_array.p()[indx]
        tau = p_array.tau()[indx]

        sigma_x, sigma_y = np.std(p_array.x()[indx]), np.std(p_array.y()[indx])
        # sigma_px, sigma_py = np.std(p_array.px()[indx]), np.std(p_array.py()[indx])
        params, params_covariance = optimize.curve_fit(test_func, tau, p, p0=[0.001, 0, 0, 0])

        b = slice_bunching(tau, charge=np.sum(p_array.q_array[indx]), lambda_mod=self.lambda_mod,
                           smooth_sigma=self.lambda_mod / 5)

        self.p.append(params[0])
        self.phi.append(params[1])
        self.s.append(np.copy(p_array.s))
        self.energy.append(np.copy(p_array.E))
        self.bunching.append(b)
        self.sigma_x.append(sigma_x)
        self.sigma_y.append(sigma_y)
        # print("center = ", (slice_max + slice_min)/2.)
        # self.sigma_px.append(sigma_px)
        # self.sigma_py.append(sigma_py)

    def finalize(self):
        data = np.array([np.array(self.s), np.array(self.p), np.array(self.phi), np.array(self.energy),
                         np.array(self.bunching), np.array(self.sigma_x), np.array(self.sigma_y)])
        np.savetxt(self.filename, data)


class Chicane(PhysProc):
    """
    simple physics process to simulate longitudinal dynamics in chicane
    """

    def __init__(self, r56, t566=0.):
        PhysProc.__init__(self)
        self.r56 = r56
        self.t566 = t566

    def apply(self, p_array, dz):
        _logger.debug(" Chicane applied, r56 =" + str(self.r56))

        p_array.rparticles[4] += (
                    self.r56 * p_array.rparticles[5] + self.t566 * p_array.rparticles[5] * p_array.rparticles[5])


class LatticeEnergyProfile(PhysProc):
    """
    The PhysProcess shifts the canonical momentum according to new reference energy Eref
    """
    def __init__(self, Eref):
        PhysProc.__init__(self)
        self.Eref = Eref

    def apply(self, p_array, dz=0):
        Eref_old = p_array.E
        p0c_old = np.sqrt(Eref_old ** 2 - m_e_GeV ** 2)
        p0c_new = np.sqrt(self.Eref ** 2 - m_e_GeV ** 2)
        p_old = p_array.p()[:]
        p_new = (p_old * p0c_old + Eref_old - self.Eref) / p0c_new
        p_array.E = self.Eref
        p_array.p()[:] = p_new[:]


class IBS(PhysProc):
    """
    Intrabeam Scattering (IBS) Physics Process.

    This class models the intrabeam scattering process in particle beams. Two methods are implemented based on different formulations:

    1. **Huang Method:** Based on Z. Huang's work: *Intrabeam Scattering in an X-ray FEL Driver* (LCLS-TN-02-8, 2002).
       URL: https://www-ssrl.slac.stanford.edu/lcls/technotes/LCLS-TN-02-8.pdf
    2. **Nagaitsev Method:** Based on S. Nagaitsev's work: *Intrabeam scattering formulas for fast numerical evaluation*
       (PRAB 8, 064403, 2005). URL: https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.8.064403

    The methods can be selected using `self.method`, which accepts "Huang" or "Nagaitsev".

    Key assumptions and features:
    - A **round beam approximation** is used, where `sigma_xy = (sigma_x + sigma_y) / 2`.
    - Beam parameters are calculated within a slice of width ±0.5 `sigma_tau` relative to the slice with maximum current
      (`self.slice = "Imax"` by default).
    - The Coulomb logarithm (`self.Clog`) is a configurable parameter, defaulting to 8 (based on Z. Huang's work).
    - A factor of 2 is applied to the Nagaitsev formula to align high-energy results with the Huang formula.

    Attributes:
    - `method` (str): Selected method for IBS calculation ("Huang" or "Nagaitsev").
    - `Clog` (float): Coulomb logarithm (default: 8) is used constant if update_Clog = False
    - `update_Clog` (bool): recalculate Clog on each step using Huang's formula without cut off.
    - `bounds` (list): Range of bounds for slice analysis in units of `sigma_tau` (default: `[-0.5, 0.5]`).
    - `slice` (str): Reference slice ("Imax" for maximum current by default).

    Methods:
    - `get_beam_params(p_array)`: Computes beam parameters such as `sigma_xy`, `sigma_z`, and normalized emittance.
    - `apply(p_array, dz)`: Applies the IBS process to a particle array over a specified distance.

    Raises:
    - `ValueError`: If an invalid method is specified.
    """

    def __init__(self, step=1, **kwargs):
        PhysProc.__init__(self)
        self.step = step  # in unit step
        self.method = kwargs.get("method", "Huang")
        self.Clog = kwargs.get("Clog", 8)
        self.update_Clog = kwargs.get("update_Clog", True)
        self.bounds = kwargs.get("bounds", [-0.5, 0.5])
        self.slice = kwargs.get("slice", "Imax")

        self.emit_n = None
        self.sigma_xy = None
        self.sigma_z = None
        self.sigma_dgamma0 = None

    def get_beam_params(self, p_array):
        """
        Compute beam parameters such as sigma_xy, sigma_z, and normalized emittance.

        :param p_array: ParticleArray
            Input particle array containing particle properties.
        """
        tws = get_envelope(p_array, bounds=self.bounds, slice=self.slice)
        self.sigma_z = np.std(p_array.tau())
        self.sigma_xy = (np.sqrt(tws.xx) + np.sqrt(tws.yy)) / 2.
        self.emit_n = (tws.emit_xn + tws.emit_yn) / 2.
        pc = np.sqrt(p_array.E ** 2 - m_e_GeV ** 2)
        self.sigma_dgamma0 = np.sqrt(tws.pp)*pc / m_e_GeV

    def estimate_Clog(self):
        """
        Used formula from Hunag's paper without cut offs

        :return: float - Coulomb logarithm
        """
        Clog = np.log(self.emit_n * self.emit_n / (ro_e * self.sigma_xy))
        return Clog

    def apply(self, p_array, dz):
        """
        Apply the IBS process to a particle array over a given path length.

        :param p_array: ParticleArray
            The particle array to be processed.
        :param dz: float
            The path length over which the IBS process is applied.

        Raises:
        - ValueError: If an invalid method is specified.
        """
        # Number of particles
        Nb = np.sum(p_array.q_array) / q_e

        # Update beam parameters
        self.get_beam_params(p_array)

        # estimate C log
        if self.update_Clog:
            Clog = self.estimate_Clog()
        else:
            Clog = self.Clog

        # Particle energy and momentum
        gamma = p_array.E / m_e_GeV

        igamma2 = 1. / (gamma * gamma)

        beta = np.sqrt(1. - igamma2)

        pc = np.sqrt(p_array.E**2 - m_e_GeV**2)

        # Calculate sigma_gamma based on the selected method
        if self.method == "Huang":
            sigma_gamma = np.sqrt(Clog * ro_e**2 * Nb / (beta**3 * self.sigma_xy * self.emit_n * self.sigma_z) * dz / 4)

        elif self.method == "Nagaitsev":
            xi = (self.sigma_dgamma0 * self.sigma_xy / (gamma * self.emit_n))**2
            F = (1 - xi**0.25) * np.log(xi + 1) / xi
            sigma_gamma = np.sqrt(Clog / 4 * ro_e**2 * Nb / (beta**3 * self.sigma_xy * self.emit_n * self.sigma_z) * F * dz)
        else:
            raise ValueError(f"Invalid method '{self.method}'. Choose 'Huang' or 'Nagaitsev'.")

        # Energy spread and update particle momenta
        sigma_e = sigma_gamma * m_e_GeV / pc
        p_array.p()[:] += np.random.randn(p_array.n) * sigma_e

