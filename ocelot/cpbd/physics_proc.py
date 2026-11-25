import logging
from ocelot.cpbd.io import save_particle_array
from ocelot.common.globals import lambda_C_r, m_e_GeV, ro_e, speed_of_light, q_e
from ocelot.cpbd.beam import Twiss, global_slice_analysis, s_to_cur, get_envelope
from ocelot.utils.acc_utils import slice_bunching
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
        if not (isinstance(self.step, (int, float, np.integer)) and float(self.step).is_integer()):
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
        super().__init__(step=step)
        self.energy = None
        self.pict_debug = True
        self.traj_step = 0.0002


class SaveBeam(PhysProc):
    def __init__(self, filename, step=1):
        super().__init__(step=step)
        self.energy = None
        self.filename = filename

    def apply(self, p_array, dz):
        _logger.debug(" SaveBeam applied, dz =" + str(dz))
        save_particle_array(filename=self.filename, p_array=p_array)


class CopyBeam(PhysProc):
    """Physics process that copies the ParticleArray instance when applied.  Makes
    most sense to be attached to zero-length elements (e.g. Marker instances)."""
    def __init__(self, name: str = "", step=1):
        super().__init__(step=step)
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
        super().__init__()
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
        super().__init__(step=step)

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
        super().__init__(step=step)
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
        super().__init__(step=step)
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
        super().__init__(step=step)
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
        super().__init__(step=step)
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
        super().__init__()
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
        p_array.match_to_twiss(self.twiss, bounds=self.bounds, slice=self.slice, remove_offsets=self.remove_offsets)


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
        super().__init__(step=step)
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
    """Spontaneous radiation loss and energy diffusion of ultrarelativistic electrons.

    For an ideal planar or helical undulator, `K` is the peak deflection
    parameter of one field component and `lperiod` is the period in metres.
    For a uniform dipole, `radius` is the signed bend radius in metres at
    the reference momentum. An infinite radius represents zero field.

    `filling_coeff` is the magnetic fraction of each step. `energy_loss`
    and `quant_diff` independently enable the mean loss and Gaussian energy
    fluctuations. Kicks use each particle's entrance energy and preserve its
    direction at fixed reference momentum.
    """

    def __init__(self, K=0.0, lperiod=0.0, type="planar", **kwargs):
        """
        Initialize spontaneous radiation effects.

        :param K: Undulator deflection parameter
        :param lperiod: Undulator period in meters
        :param type: "planar"/"helical" undulator or "dipole"
        :param kwargs: Additional keyword arguments for customization
        """
        step = kwargs.pop("step", 1)
        super().__init__(step=step)

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
        if not np.isfinite(dz) or dz < 0:
            raise ValueError("Radiation step length must be finite and nonnegative.")
        filling = np.asarray(self.filling_coeff)
        if not np.all((0 <= filling) & (filling <= 1)):
            raise ValueError("Magnetic filling fraction must be between zero and one.")
        length = dz * filling
        if np.all(length == 0) or p_array.n == 0 or not (self.energy_loss or self.quant_diff):
            return

        if self.type == "dipole":
            radius = np.abs(self.radius)
            if np.any(np.isnan(radius) | (radius == 0)):
                raise ValueError("Dipole radius must be nonzero.")
            if np.all(np.isinf(radius)):
                return
        elif self.type in ("planar", "helical") and np.all(np.equal(self.K, 0)):
            return

        if not np.isfinite(p_array.E) or p_array.E <= m_e_GeV:
            raise ValueError("Radiation requires a reference energy above the electron rest energy.")
        p0c = p_array.p0c
        energy = p_array.energies
        if not np.all(np.isfinite(energy) & (energy > m_e_GeV)):
            raise ValueError("Radiation requires particle energies above the electron rest energy.")
        momentum = np.sqrt(energy**2 - m_e_GeV**2)
        delta_energy = 0.
        if self.type == "dipole":
            # Fixed transverse field. Yampolsky & Carlsten, Eqs. (17)-(20).
            curvature = (p0c / momentum) / radius
            gamma = energy / m_e_GeV
            if self.energy_loss:
                delta_energy = -(2 / 3) * m_e_GeV * ro_e * gamma**4 * curvature**2 * length
            if self.quant_diff:
                variance_gamma = (55 / (24 * np.sqrt(3)) * ro_e * lambda_C_r
                                  * gamma**7 * curvature**3 * length)
                sigma_energy = m_e_GeV * np.sqrt(variance_gamma)
        else:
            if self.energy_loss:
                delta_energy = -self.energy_loss_und(energy, length)
            if self.quant_diff:
                sigma_energy = energy * self.sigma_gamma_quant(energy, length, self.K, self.lperiod, self.type)

        if self.quant_diff:
            delta_energy += sigma_energy * np.random.randn(p_array.n)

        energy_after = energy + delta_energy
        if not np.all(np.isfinite(energy_after) & (energy_after > m_e_GeV)):
            raise ValueError("Radiation kick leaves a particle with invalid energy. Reduce the step length.")
        momentum_ratio = np.sqrt(energy_after**2 - m_e_GeV**2) / momentum
        p_array.px()[:] *= momentum_ratio
        p_array.py()[:] *= momentum_ratio
        p_array.p()[:] += delta_energy / p0c

    def energy_loss_und(self, energy, dz):
        """Mean energy radiated in GeV over magnetic length `dz` in metres."""
        k_u = self._undulator_wavenumber(self.K, self.lperiod, self.type)
        mean_k_squared = np.square(self.K) * (1.0 if self.type == "helical" else 0.5)
        return (2 / 3) * ro_e / m_e_GeV * energy**2 * k_u**2 * mean_k_squared * dz

    @staticmethod
    def sigma_gamma_quant(energy, dz, K, lperiod, type="planar"):
        """Induced rms sigma_E/E for energy [GeV] and undulator magnetic length dz [m]."""
        k_u = SpontanRadEffects._undulator_wavenumber(K, lperiod, type)
        K = np.abs(K)
        # Saldin, Schneidmiller & Yurkov, NIM A 381 (1996), 545-547 (peak K).
        if type == "helical":
            fit = 1.42 * K + 1 / (1 + 1.50 * K + 0.95 * K**2)
            coefficient = 14 / 15
        else:
            fit = 1.20 * K + 1 / (1 + 1.33 * K + 0.40 * K**2)
            coefficient = 7 / 15
        return (energy / m_e_GeV) * np.sqrt(coefficient * lambda_C_r * ro_e * k_u**3 * K**2 * fit * dz)

    @staticmethod
    def _undulator_wavenumber(K, lperiod, type):
        if type not in ("planar", "helical"):
            raise ValueError(f"Unknown undulator type {type!r}. Expected 'planar' or 'helical'.")
        if np.all(np.equal(K, 0)):
            return np.zeros_like(lperiod, dtype=float)
        if not np.all(np.isfinite(K)) or not np.all(np.isfinite(lperiod) & (lperiod > 0)):
            raise ValueError("An active undulator requires finite K and a positive finite period.")
        return 2 * np.pi / lperiod


class BeamAnalysis(PhysProc):
    def __init__(self, filename):
        super().__init__()
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
        super().__init__()
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
        super().__init__()
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
    Longitudinal intrabeam-scattering heating for a round electron beam in a linac.

    ``method="Huang"`` uses the longitudinally cold limit of Z. Huang,
    *Intrabeam Scattering in an X-ray FEL Driver*, LCLS-TN-02-8 (2002),
    Eqs. (5)-(9): https://www-ssrl.slac.stanford.edu/lcls/technotes/LCLS-TN-02-8.pdf.
    ``method="Nagaitsev"`` uses the round-beam rate with the approximate
    G function in S. Nagaitsev, PRSTAB 8, 064403 (2005), Eqs. (45)-(49):
    https://doi.org/10.1103/PhysRevSTAB.8.064403. The fit has approximately
    20% accuracy; this is not the general elliptic-integral IBS model.

    The Nagaitsev longitudinal variance rate is multiplied by two because
    there are no synchrotron oscillations in a linac (Xiao and Borland,
    PAC09, TH5PFP038, Eq. (4):
    https://accelconf.web.cern.ch/pac2009/papers/th5pfp038.pdf).
    This doubles the kick variance, not its standard deviation.

    Parameters
    ----------
    step : int
        Number of Navigator unit steps between applications (default: 1).
    method : str
        "Huang" (default) or "Nagaitsev".
    Clog : float
        Fixed Coulomb logarithm, used only when ``update_Clog=False``.
        Default: 8, approximating Huang's cutoff for the Gaussian core.
    update_Clog : bool
        Default: True. Recalculate the uncut, high-energy Coulomb-log
        estimate each step instead of using ``Clog``. ``Clog_used`` records
        the effective value from the last nonzero-charge, positive-length step.
    bounds : sequence of float or None
        Analysis window in units of the whole-bunch rms length, relative
        to ``slice``. Default: [-0.5, 0.5]; None selects the whole bunch.
    slice : str or None
        "Imax" (default) centers the window on the current maximum;
        otherwise it is centered on the mean longitudinal position.

    Notes
    -----
    ``sigma_dgamma0`` measures energy spread after subtracting a fitted
    linear chirp within the analysis window, without modifying the particles.
    Curvature and variation between slices require a narrower window.
    Transverse sizes and normalized emittances are averaged between planes.
    One kick strength is applied to the entire bunch using its total charge
    and rms length: ``sigma_z`` stores std(tau), with tau = c*delta_t; the
    physical bunch length in the rate is beta*sigma_z. Local density
    variation, transverse IBS, and dispersion
    effects are not modelled. Equal macro-particle weights and approximately
    Gaussian, transversely round distributions are assumed.

    Nagaitsev heating requires a longitudinal/transverse rest-frame
    temperature ratio xi <= 1. A larger ratio raises ValueError before any
    kick: longitudinal cooling and the accompanying transverse heating
    require a coupled model. Huang requires xi << 1. Both methods use the
    finite-beta conversion of Nagaitsev's rate to Ocelot's energy coordinate;
    Huang is recovered as the cold limit and agrees with its reference at
    beta approximately 1. This does not extend the Coulomb-log estimate's
    validity to nonrelativistic beams.
    """

    def __init__(self, step=1, **kwargs):
        super().__init__(step=step)
        self.method = kwargs.get("method", "Huang")
        self.Clog = kwargs.get("Clog", 8)
        self.update_Clog = kwargs.get("update_Clog", True)
        self.bounds = kwargs.get("bounds", [-0.5, 0.5])
        self.slice = kwargs.get("slice", "Imax")

        self.emit_n = None
        self.sigma_xy = None
        self.sigma_z = None
        self.sigma_dgamma0 = None
        self.Clog_used = None

    def get_beam_params(self, p_array: ParticleArray) -> None:
        """
        Compute transverse moments and the energy spread after removing linear chirp.

        :param p_array: ParticleArray
            Input particle array containing particle properties.
        """
        if not np.isfinite(p_array.E) or p_array.E <= m_e_GeV:
            raise ValueError("IBS requires a finite reference energy above the electron rest energy.")
        if p_array.n < 3 or not np.all(np.isfinite(p_array.rparticles)):
            raise ValueError("IBS requires at least three particles with finite coordinates.")
        self.sigma_z = np.std(p_array.tau())
        if self.sigma_z <= 0 or not np.isfinite(self.sigma_z):
            raise ValueError("IBS requires a positive, finite rms bunch length.")

        # Select once so the transverse moments and chirp fit use the same particles.
        selected = p_array
        if self.bounds is not None:
            bounds = np.asarray(self.bounds, dtype=float)
            if bounds.shape != (2,) or not np.all(np.isfinite(bounds)) or bounds[0] >= bounds[1]:
                raise ValueError("IBS bounds must contain two finite, increasing values.")
            if self.slice == "Imax":
                current = s_to_cur(p_array.tau(), 0.01 * self.sigma_z,
                                   np.sum(p_array.q_array), speed_of_light)
                center = current[np.argmax(current[:, 1]), 0]
            else:
                center = np.mean(p_array.tau())
            mask = ((center + self.sigma_z * bounds[0] <= p_array.tau()) &
                    (p_array.tau() <= center + self.sigma_z * bounds[1]))
            selected = ParticleArray()
            selected.E = p_array.E
            selected.rparticles = p_array.rparticles[:, mask]
            selected.q_array = p_array.q_array[mask]
        if selected.n < 3:
            raise ValueError("IBS analysis window must contain at least three particles.")

        with np.errstate(divide="ignore", invalid="ignore"):
            tws = get_envelope(selected)
        if not np.all(np.isfinite([tws.emit_xn, tws.emit_yn])) or min(tws.emit_xn, tws.emit_yn) <= 0:
            raise ValueError("IBS requires positive, finite emittances in both transverse planes.")
        self.sigma_xy = (np.sqrt(tws.xx) + np.sqrt(tws.yy)) / 2.
        self.emit_n = (tws.emit_xn + tws.emit_yn) / 2.

        tau = selected.tau() - np.mean(selected.tau())
        residual = selected.p() - np.mean(selected.p())
        tau_variance = np.mean(tau * tau)
        if tau_variance > 0:
            residual = residual - np.mean(tau * residual) / tau_variance * tau
        pc = np.sqrt(p_array.E ** 2 - m_e_GeV ** 2)
        self.sigma_dgamma0 = np.std(residual) * pc / m_e_GeV

    def estimate_Clog(self) -> float:
        """
        Huang's high-energy estimate without the Gaussian-core cutoff.

        :return: float - Coulomb logarithm
        """
        Clog = np.log(self.emit_n * self.emit_n / (ro_e * self.sigma_xy))
        return Clog

    def apply(self, p_array: ParticleArray, dz: float) -> None:
        """
        Apply the IBS process to a particle array over a given path length.

        :param p_array: ParticleArray
            The particle array to be processed.
        :param dz: float
            The path length over which the IBS process is applied.

        Raises:
        - ValueError: For invalid inputs or unsupported Nagaitsev cooling.
        """
        if self.method not in ("Huang", "Nagaitsev"):
            raise ValueError(f"Invalid method '{self.method}'. Choose 'Huang' or 'Nagaitsev'.")
        if not np.isfinite(dz) or dz < 0:
            raise ValueError("IBS path length must be finite and nonnegative.")
        if dz == 0:
            return

        Nb = np.sum(p_array.q_array) / q_e
        if not np.isfinite(Nb) or Nb < 0:
            raise ValueError("IBS requires a finite, nonnegative bunch charge.")
        if Nb == 0:
            return

        # Update beam parameters
        self.get_beam_params(p_array)

        # estimate C log
        if self.update_Clog:
            Clog = self.estimate_Clog()
        else:
            Clog = self.Clog
        if not np.isfinite(Clog) or Clog < 0:
            raise ValueError("IBS Coulomb logarithm must be finite and nonnegative.")
        self.Clog_used = Clog
        if Clog == 0:
            return

        # Particle energy and momentum
        gamma = p_array.E / m_e_GeV

        igamma2 = 1. / (gamma * gamma)

        beta = np.sqrt(1. - igamma2)

        pc = np.sqrt(p_array.E**2 - m_e_GeV**2)

        # Nagaitsev Eq. (46) uses delta = dP/P0, whereas Ocelot uses dE/(P0*c).
        # dgamma = beta**2 * gamma * delta; emit_n = beta*gamma*sigma_xy*theta.
        # Convert dt to ds/(beta*c) and double the ring variance rate for a linac.
        # The cold-limit rate has a beta prefactor for a physical bunch length.
        # Ocelot tau = c*delta_t, so that length is beta*std(tau).
        sigma_s = beta * self.sigma_z
        G = 1.
        if self.method == "Nagaitsev":
            xi = (self.sigma_dgamma0 * self.sigma_xy / (beta * gamma * self.emit_n))**2
            if not np.isfinite(xi) or xi > 1:
                raise ValueError(
                    f"Nagaitsev IBS requires xi <= 1 (got {xi:g}). Longitudinal cooling "
                    "and accompanying transverse heating are not supported by this longitudinal-kick model."
                )
            # Eq. (49), with its finite cold-beam limit and stable log near zero.
            G = (1 - xi**0.25) * np.log1p(xi) / xi if xi > 0 else 1.

        variance_gamma = Clog * ro_e**2 * Nb * beta * G * dz / (
            4 * self.sigma_xy * self.emit_n * sigma_s)
        if not np.isfinite(variance_gamma) or variance_gamma < 0:
            raise ValueError("IBS produced an invalid energy-variance increment.")
        if variance_gamma == 0:
            return

        # Energy spread and update particle momenta
        sigma_e = np.sqrt(variance_gamma) * m_e_GeV / pc
        p_array.p()[:] += np.random.randn(p_array.n) * sigma_e
