"""
@author: Igor Zagorodnov @ Martin Dohlus
Created on 27.03.2015
Revision on 01.06.2017: coordinate transform to the velocity direction
2019: Added LSC: S. Tomin and I. Zagorodnov
"""
import numpy as np
import scipy.ndimage as ndimage
import time
from ocelot.common.globals import *
from ocelot.cpbd.coord_transform import *
from scipy import interpolate
import multiprocessing
from scipy.special import exp1, k1
from ocelot.cpbd.physics_proc import PhysProc
from ocelot.common.math_op import conj_sym
from ocelot.cpbd.beam import s_to_cur
from ocelot.common import conf
from ocelot.cpbd.elements import Undulator
from scipy.interpolate import interp1d
import logging

try:
    from scipy.special import factorial
except:
    from scipy.misc import factorial    # legacy support

logger = logging.getLogger(__name__)

try:
    pyfftw_flag = True
    from pyfftw.interfaces.numpy_fft import fftn
    from pyfftw.interfaces.numpy_fft import ifftn
    import pyfftw
except:
    pyfftw_flag = False
    logger.debug("cs.py: module PYFFTW is not installed. Install it to speed up calculation")
    from numpy.fft import ifftn
    from numpy.fft import fftn

try:
    import numexpr as ne
    ne_flag = True
except:
    logger.debug("sc.py: module NUMEXPR is not installed. Install it to speed up calculation")
    ne_flag = False

def smooth_z(Zin, mslice):

    def myfunc(x, A):
        if x < 2*A:
            y = x - x*x/(4*A)
        else:
            y = A
        return y

    inds = np.argsort(Zin, axis=0)
    Zout = np.sort(Zin, axis=0)
    N = Zin.shape[0]
    S = np.zeros(N+1)
    S[N] = 0
    S[0] = 0
    for i in range(N):
        S[i+1] = S[i] + Zout[i]
    Zout2 = np.zeros(N)
    Zout2[N-1] = Zout[N-1]
    Zout2[0] = Zout[0]
    for i in range(1, N-1):
        m = min(i, N-i+1)
        m = np.floor(myfunc(0.5*m, 0.5*mslice) + 0.500001).astype(int)
        Zout2[i] = (S[i+m+1] - S[i-m])/(2*m + 1)
    Zout[inds] = Zout2
    return Zout


class SpaceCharge(PhysProc):
    """
    Space Charge physics process

    Attributes:
        self.step = 1 [in Navigator.unit_step] - step of the Space Charge kick applying
        self.nmesh_xyz = [63, 63, 63] - 3D mesh

    Description:
        The space charge forces are calculated by solving the Poisson equation in the bunch frame.
    Then the Lorentz transformed electromagnetic field is applied as a kick in the laboratory frame.
    For the solution of the Poisson equation we use an integral representation of the electrostatic potential
    by convolution of the free-space Green's function with the charge distribution.
    The convolution equation is solved with the help of the Fast Fourier Transform (FFT). The same algorithm for
    solution of the 3D Poisson equation is used, for example, in ASTRA
    """
    def __init__(self, step=1, **kwargs):
        PhysProc.__init__(self)
        self.step = step # in unit step
        self.nmesh_xyz = kwargs.get("nmesh_xyz", [63, 63, 63])
        self.low_order_kick = kwargs.get("low_order_kick", True)

        self.start_elem = None
        self.end_elem = None
        self.debug = False
        self.random_mesh = kwargs.get("random_mesh", False)  # if True mesh is shifted slightly on each step in order to reduce numerical noise
        self.random_seed = 10     # random seeding number. if None seeding is random

    def prepare(self, lat):
        self.check_step()
        if self.random_seed is not None:
            np.random.seed(self.random_seed)

    def sym_kernel(self, ijk2, hxyz):
        i2 = ijk2[0]
        j2 = ijk2[1]
        k2 = ijk2[2]
        hx = hxyz[0]
        hy = hxyz[1]
        hz = hxyz[2]
        x = hx*np.r_[0:i2+1] - hx/2
        y = hy*np.r_[0:j2+1] - hy/2
        z = hz*np.r_[0:k2+1] - hz/2
        x, y, z = np.ix_(x, y, z)
        r = np.sqrt(x*x + y*y + z*z)
        if ne_flag:
            IG = ne.evaluate("(-x*x*0.5*arctan(y*z/(x*r)) + y*z*log(x+r) - y*y*0.5*arctan(z*x/(y*r)) + z*x*log(y+r) - z*z*0.5*arctan(x*y/(z*r)) + x*y*log(z+r))")
        else:
            IG = (-x*x*0.5*np.arctan(y*z/(x*r)) + y*z*np.log(x+r)
                  -y*y*0.5*np.arctan(z*x/(y*r)) + z*x*np.log(y+r)
                  -z*z*0.5*np.arctan(x*y/(z*r)) + x*y*np.log(z+r))

        kern = (IG[1:i2+1, 1:j2+1, 1:k2+1] - IG[0:i2, 1:j2+1, 1:k2+1]
               -IG[1:i2+1, 0:j2, 1:k2+1] + IG[0:i2, 0:j2, 1:k2+1]
               -IG[1:i2+1, 1:j2+1, 0:k2] + IG[0:i2, 1:j2+1, 0:k2]
               +IG[1:i2+1, 0:j2, 0:k2] - IG[0:i2, 0:j2, 0:k2])

        return kern

    def potential(self, q, steps):
        hx = steps[0]
        hy = steps[1]
        hz = steps[2]
        Nx = q.shape[0]
        Ny = q.shape[1]
        Nz = q.shape[2]
        out = np.zeros((2*Nx-1, 2*Ny-1, 2*Nz-1))
        out[:Nx, :Ny, :Nz] = q
        K1 = self.sym_kernel(q.shape, steps)
        K2 = np.zeros((2*Nx-1, 2*Ny-1, 2*Nz-1))
        K2[0:Nx, 0:Ny, 0:Nz] = K1
        K2[0:Nx, 0:Ny, Nz:2*Nz-1] = K2[0:Nx, 0:Ny, Nz-1:0:-1] #z-mirror
        K2[0:Nx, Ny:2*Ny-1,:] = K2[0:Nx, Ny-1:0:-1, :]        #y-mirror
        K2[Nx:2*Nx-1, :, :] = K2[Nx-1:0:-1, :, :]             #x-mirror
        t0 = time.time()
        if pyfftw_flag:
            nthreads = int(conf.OCELOT_NUM_THREADS)
            if nthreads < 1:
                nthreads = 1
            K2_fft = pyfftw.builders.fftn(K2, axes=None, overwrite_input=False, planner_effort='FFTW_ESTIMATE',
                                       threads=nthreads, auto_align_input=False, auto_contiguous=False, avoid_copy=True)
            out_fft = pyfftw.builders.fftn(out, axes=None, overwrite_input=False, planner_effort='FFTW_ESTIMATE',
                                          threads=nthreads, auto_align_input=False, auto_contiguous=False, avoid_copy=True)
            out_ifft = pyfftw.builders.ifftn(out_fft()*K2_fft(), axes=None, overwrite_input=False, planner_effort='FFTW_ESTIMATE',
                                          threads=nthreads, auto_align_input=False, auto_contiguous=False, avoid_copy=True)
            out = np.real(out_ifft())

        else:
            out = np.real(ifftn(fftn(out)*fftn(K2)))
        t1 = time.time()
        logger.debug('fft time:' + str(t1-t0) + ' sec')
        out[:Nx, :Ny, :Nz] = out[:Nx,:Ny,:Nz]/(4*pi*epsilon_0*hx*hy*hz)
        return out[:Nx, :Ny, :Nz]

    def el_field(self, X, Q, gamma, nxyz):
        N = X.shape[0]
        X[:, 2] = X[:, 2] * gamma
        XX = np.max(X, axis=0) - np.min(X, axis=0)
        if self.random_mesh:
            XX = XX * np.random.uniform(low=1, high=1.1)
        logger.debug('mesh steps:' + str(XX))
        # here we use a fast 3D "near-point" interpolation
        # we need a stand-alone module with 1D,2D,3D parricles-to-grid functions
        steps = XX / (nxyz - 3)
        X = X / steps
        X_min = np.min(X, axis=0)
        X_mid = np.dot(Q, X) / np.sum(Q)
        X_off = np.floor(X_min - X_mid) + X_mid
        if self.random_mesh:
            X_off = X_off + np.random.uniform(low=-0.5, high=0.5)
        X = X - X_off
        nx = nxyz[0]
        ny = nxyz[1]
        nz = nxyz[2]
        nzny = nz * ny
        Xi = np.int_(np.floor(X) + 1)
        inds = np.int_(Xi[:, 0] * nzny + Xi[:, 1] * nz + Xi[:, 2])  # 3d -> 1d
        q = np.bincount(inds, Q, nzny * nx).reshape(nxyz)
        p = self.potential(q, steps)
        Ex = np.zeros(p.shape)
        Ey = np.zeros(p.shape)
        Ez = np.zeros(p.shape)
        Ex[:nx - 1, :, :] = (p[:nx - 1, :, :] - p[1:nx, :, :]) / steps[0]
        Ey[:, :ny - 1, :] = (p[:, :ny - 1, :] - p[:, 1:ny, :]) / steps[1]
        Ez[:, :, :nz - 1] = (p[:, :, :nz - 1] - p[:, :, 1:nz]) / steps[2]
        Exyz = np.zeros((N, 3))
        Exyz[:, 0] = ndimage.map_coordinates(Ex, np.c_[X[:, 0], X[:, 1] + 0.5, X[:, 2] + 0.5].T, order=1) * gamma
        Exyz[:, 1] = ndimage.map_coordinates(Ey, np.c_[X[:, 0] + 0.5, X[:, 1], X[:, 2] + 0.5].T, order=1) * gamma
        Exyz[:, 2] = ndimage.map_coordinates(Ez, np.c_[X[:, 0] + 0.5, X[:, 1] + 0.5, X[:, 2]].T, order=1)
        return Exyz


    def apply(self, p_array, zstep):
        logger.debug(" apply: zstep = " + str(zstep))
        if zstep == 0:
            logger.debug(" apply: zstep = 0 -> return ")
            return
        nmesh_xyz = np.array(self.nmesh_xyz)
        gamref = p_array.E / m_e_GeV
        betref2 = 1 - gamref ** -2
        betref = np.sqrt(betref2)

        # MAD coordinates!!!
        # Lorentz transformation with V-axis and gamma_av
        xp = np.zeros(p_array.rparticles.shape)
        xp = xxstg_2_xp_mad(p_array.rparticles, xp, gamref)

        # coordinate transformation to the velocity direction
        t3 = np.mean(xp[3:6], axis=1)
        Pav = np.linalg.norm(t3)
        t3 = t3 / Pav
        ey = np.array([0, 1, 0])
        t1 = np.cross(ey, t3)
        t1 = t1 / np.linalg.norm(t1)
        t2 = np.cross(t3, t1)
        T = np.c_[t1, t2, t3]

        xyz = np.dot(xp[0:3].T, T)
        xp[3:6] = np.dot(xp[3:6].T, T).T

        # electric field in the rest frame of bunch
        gamma0 = np.sqrt((Pav / m_e_eV) ** 2 + 1)
        beta02 = 1 - gamma0 ** -2
        beta0 = np.sqrt(beta02)

        Exyz = self.el_field(xyz, p_array.q_array, gamma0, nmesh_xyz)

        # equations of motion in the lab system
        cdT = zstep / betref

        xp[3] = xp[3] + cdT * (1 - beta0 * beta0) * Exyz[:, 0]
        xp[4] = xp[4] + cdT * (1 - beta0 * beta0) * Exyz[:, 1]
        xp[5] = xp[5] + cdT *  Exyz[:, 2]
        T = np.transpose(T)
        xp[3:6] = np.dot(xp[3:6].T, T).T
        xp_2_xxstg_mad(xp, p_array.rparticles, gamref)

    def __repr__(self) -> str:
        cname = type(self).__name__
        step = self.step
        nmesh_xyz = self.nmesh_xyz
        random_mesh = self.random_mesh
        return f"<{cname}: {step=}, {nmesh_xyz=}, {random_mesh=}>"


class LSC(PhysProc):
    """
    Longitudinal Space Charge (LSC) impedance model.

    This class simulates the LSC effect using an impedance-based approach.
    It supports both Gaussian and uniform (step-profile) transverse beam distributions,
    and accounts for undulator sections via a position-dependent longitudinal Lorentz factor γ_z.

    Parameters:
    -----------
    step : int
        Number of lattice elements between applications of the LSC kick (default: 1).
    step_profile : bool
        If True, uses the uniform transverse beam model (step profile). Default is False (Gaussian).
    smooth_param : float
        Smoothing parameter for the current profile (KDE bandwidth).
        Effective resolution is std_tau * smooth_param. Default: 0.1.
    bounds : list of float
        Fractional bounds [min, max] (in units of σ_τ) to define the central slice of the bunch
        used for calculating the transverse beam size. Default: [-0.4, 0.4].
    slice : optional, not yet implemented
        Custom slice index or mask for selecting a portion of the bunch for analysis (default: None).

    Notes:
    ------
    The impedance model for a round Gaussian beam was originally taken from [1] and implemented in Matlab
    by I. Zagorodnov. This Python version was ported and integrated into Ocelot.

    References:
    -----------
    [1] Geloni et al., NIM A 578 (2007) 34-46. https://arxiv.org/abs/physics/0612077
    """
    def __init__(self, step=1, **kwargs):
        PhysProc.__init__(self, step)
        self.K_s_func = None
        self.step_profile = kwargs.get("step_profile", False)
        self.smooth_param = kwargs.get("smooth_param", 0.1)
        self.bounds = kwargs.get("bounds", [-0.4, 0.4])
        self.slice = kwargs.get("slice", None)
        self._is_undul_in_beam_line = False

    def imp_lsc(self, gamma, sigma, w, dz):
        """
        Calculates LSC impedance for a beam with a round Gaussian profile.

        Parameters
        ----------
        gamma : float
            The relativistic Lorentz factor of the beam.
        sigma : float
            The transverse RMS size of the beam in meters.
        w : ndarray
            An array of angular frequencies (omega = 2*pi*f) in rad/s.
        dz : float
            The length of the beamline element in meters.

        Returns
        -------
        ndarray
            The complex longitudinal impedance Z(w) in Ohms. Note: The raw
            impedance is per unit length [Ohm/m], which is then multiplied by `dz`.
        """
        eps = 1e-16
        ass = 40.0

        alpha = w * sigma / (gamma * speed_of_light)
        alpha2 = alpha * alpha

        inda = np.where(alpha2 > ass)[0]
        ind = np.where((alpha2 <= ass) & (alpha2 >= eps))[0]

        T = np.zeros(w.shape)
        T[ind] = np.exp(alpha2[ind]) * exp1(alpha2[ind])

        x = alpha2[inda]
        k = 0
        for i in range(10):
            k += (-1) ** i * factorial(i) / (x ** (i + 1))
        T[inda] = k
        Z = 1j * Z0 / (4 * pi * speed_of_light*gamma**2) * w * T * dz
        return Z # --> Omm/m

    def imp_step_lsc(self, gamma, rb, w, dz):
        """
        Calculates LSC impedance for a beam with a uniform transverse profile.

        Parameters
        ----------
        gamma : float
            The relativistic Lorentz factor of the beam.
        rb : float
            The transverse radius of the beam in meters.
        w : ndarray
            An array of angular frequencies (omega = 2*pi*f) in rad/s.
        dz : float
            The length of the beamline element in meters.

        Returns
        -------
        ndarray
            The complex longitudinal impedance Z(w) in Ohms. Note: The raw
            impedance is per unit length [Ohm/m], which is then multiplied by `dz`.
        """
        indx = np.where(w < 1e-7)[0]
        w[indx] = 1e-7
        x = w*rb/(speed_of_light*gamma)
        Z = 1j*Z0 * speed_of_light / (4 * w * rb*rb) * dz * (1 - x*k1(x))
        Z[indx] = 0
        return Z # --> Omm/m

    def wake2impedance(self, s, w):
        """
        Fourier transform with exp(iwt)
        s - Meter
        w - V/C
        f - Hz
        y - Om
        """
        ds = s[1] - s[0]
        dt = ds / speed_of_light
        n = len(s)
        f = 1 / dt * np.arange(0, n) / n
        shift = 1#np.exp(1j * f * t0 * 2 * np.pi)
        y = dt * np.fft.fft(w, n) * shift
        return f, y

    def impedance2wake(self, f, y):
        """
        Fourier transform with exp(-iwt)
        f - Hz
        y - Om
        s - Meter
        w - V/C
        """
        df = f[1] - f[0]
        n = len(f)
        s = 1 / df * np.arange(0, n) / n * speed_of_light
        # general case
        # w1 = n * df * np.fft.ifft(conj_sym(y), n).real
        w = n * df * np.fft.irfft(y, n)
        return s, w

    def wake_lsc(self, s, bunch, gamma, sigma, dz, K_max=0, fill_factor=0):
        """
        Computes the LSC wakefield for a given bunch profile.

        This function calculates the impedance (for either Gaussian or step-profile
        beams), multiplies it by the bunch spectrum in the frequency domain, and
        transforms the result back to the time domain to get the wake potential.

        Parameters
        ----------
        s : ndarray
            Longitudinal coordinates array [m].
        bunch : ndarray
            The normalized longitudinal current profile (dimensionless).
        gamma : float
            The relativistic Lorentz factor.
        sigma : float
            The characteristic transverse beam size (RMS or radius) [m].
        dz : float
            The length of the interaction region [m].
        K_max : float, optional
            The maximum undulator K-parameter in the section. Defaults to 0.
        fill_factor : float, optional
            The fraction of the length `dz` that is occupied by an undulator
            field. Defaults to 0.

        Returns
        -------
        ndarray
            The wake potential W(s) in Volts.
        """
        ds = s[1] - s[0]
        dt = ds / speed_of_light
        nb = len(s)
        n = nb * 2
        f = 1 / dt * np.arange(0, n) /n
        if self.step_profile:
            # space charge impedance of transverse step profile
            Za = self.imp_step_lsc(gamma, rb=sigma, w=f[0:nb] * 2 * np.pi, dz=dz) * (1 + 0.5 * K_max*K_max * fill_factor)
        else:
            Za = self.imp_lsc(gamma, sigma, w=f[0:nb] * 2 * np.pi, dz=dz) * (1 + 0.5 * K_max*K_max * fill_factor)

        sb1 = s[0] + np.cumsum(np.ones(n)) * ds

        bunch1 = np.append(bunch, np.zeros(nb))

        f, Zb = self.wake2impedance(sb1, bunch1 * speed_of_light)

        Z = np.zeros(n, dtype=complex)
        Z[0:nb] = Za * Zb[0:nb]
        Z[nb:n] = np.flipud(np.conj(Z[0:nb]))

        xa, wa = self.impedance2wake(f, Z)
        res = -wa[0:nb]
        return res

    def prepare(self, lat):
        """
        Scans the lattice to create an undulator strength K(s) profile.

        This setup method is called once before the tracking. It builds an
        interpolation function for the undulator K-value as a function of
        longitudinal position 's' along the beamline. This is used later by
        the `apply` method to account for LSC effects within undulators.

        Parameters
        ----------
        lat : Lattice
            The accelerator lattice object containing the sequence of beamline
            elements.

        Returns
        -------
        None
            Sets the internal `self.K_s_func` and `self._is_undul_in_beam_line`.
        """
        self.check_step()
        seq = lat.get_sequence_part(self.start_elem, self.end_elem)

        s = []
        k = []

        s_current = self.s_start
        for elem in seq:
            s_next = s_current + elem.l
            s.extend([s_current, s_next])
            if isinstance(elem, Undulator) and not (elem.Kx != 0 and elem.Ky != 0):
                k.extend([elem.Kx + elem.Ky, elem.Kx + elem.Ky])
                self._is_undul_in_beam_line = True
            else:
                k.extend([0.0, 0.0])
            s_current = s_next
        self.K_s_func = interp1d(s, k, kind='linear', bounds_error=False, fill_value=(k[0], k[-1]))

    def compute_filling_factor(self, x0, x1, num_points=100):
        """Calculates the undulator filling factor over a longitudinal segment.

        The filling factor is the ratio of the length occupied by an undulator
        field to the total length of the segment. This is used to correctly
        scale the LSC effect in undulators.

        Parameters
        ----------
        x0 : float
            The start position of the segment [m].
        x1 : float
            The end position of the segment [m].
        num_points : int, optional
            The number of points to use for the numerical integration.
            Defaults to 100.

        Returns
        -------
        float
            The filling factor (dimensionless, between 0 and 1).
        """
        x = np.linspace(x0, x1, num_points)
        k_vals = self.K_s_func(x)

        # For each interval between x[i] and x[i+1], check if K is nonzero at least at one end
        nonzero_intervals = (k_vals[:-1] != 0) | (k_vals[1:] != 0)
        dx = (x1 - x0) / (num_points - 1)
        filled_length = np.sum(nonzero_intervals) * dx
        total_length = x1 - x0

        return filled_length / total_length

    def apply(self, p_array, dz):
        """Applies the LSC kick to a particle array for a given integration step.

        This method calculates the LSC-induced wake and applies the corresponding
        energy kick to the particles. It handles both drift spaces and undulators
        by calculating an effective filling factor for the undulator fields.

        Parameters
        ----------
        p_array : ParticleArray
            The particle bunch object. It is modified in-place.
        dz : float
            The length of the integration step in meters.

        Returns
        -------
        None
            The `p_array` object is modified in-place.
        """
        if dz < 1e-10:
            logger.debug(" LSC applied, dz < 1e-10, dz = " + str(dz))
            return
        if self._is_undul_in_beam_line:
            fill_factor = self.compute_filling_factor(self.z0-dz, self.z0, num_points=200)
            K_max = np.max(self.K_s_func(np.linspace(self.z0-dz, self.z0, num=200)))
        else:
            fill_factor = 0
            K_max = 0

        logger.debug(" LSC applied, dz =" + str(dz))
        tau = p_array.tau()
        mean_tau = np.mean(tau)
        sigma_tau = np.std(tau)

        slice_min = mean_tau + sigma_tau * self.bounds[0]
        slice_max = mean_tau + sigma_tau * self.bounds[1]

        indx = np.where((tau >= slice_min) & (tau < slice_max))

        if self.step_profile:
            rb = min(np.max(p_array.x()[indx]) - np.min(p_array.x()[indx]),
                     np.max(p_array.y()[indx]) - np.min(p_array.y()[indx]))/2
            sigma = rb
        else:
            sigma = (np.std(p_array.x()[indx]) + np.std(p_array.y()[indx]))/2.
            # sigma = min(np.std(p_array.x()[indx]), np.std(p_array.y()[indx]))
        q = np.sum(p_array.q_array)
        gamma = p_array.E / m_e_GeV
        v = np.sqrt(1 - 1 / gamma ** 2) * speed_of_light
        B = s_to_cur(p_array.tau(), sigma_tau * self.smooth_param, q, v)
        bunch = B[:, 1] / (q * speed_of_light)
        x = B[:, 0]

        W = - self.wake_lsc(x, bunch, gamma, sigma, dz, K_max, fill_factor) * q

        indx = np.argsort(p_array.tau(), kind="quicksort")
        tau_sort = p_array.tau()[indx]
        dE = np.interp(tau_sort, x, W)

        pc_ref = np.sqrt(p_array.E ** 2 / m_e_GeV ** 2 - 1) * m_e_GeV
        delta_p = dE * 1e-9 / pc_ref
        p_array.rparticles[5][indx] += delta_p


        #fig, axs = plt.subplots(3, 1, sharex=True)
        ##ax = plt.subplot(211)
        #axs[0].plot(-B[:, 0]*1e3, B[:, 1])
        #axs[0].set_ylabel("I, [A]")
        #axs[1].plot(-p_array.tau()[::10]*1e3, p_array.p()[::10], ".")
        #axs[1].set_ylim(-0.01, 0.01)
        #axs[1].set_ylabel("dE/E")
        ##plt.subplot(212)
        #axs[2].plot(-x*1e3, W, label="s = "+str(p_array.s) + "  m")
        #axs[2].set_ylabel("W, [V]")
        #axs[2].set_xlabel("s, [mm]")
        #plt.legend()
        #plt.show()
        ##plt.ylim(-0.5e6, 0.5e6)
        #dig = str(self.napply)
        #name = "0"*(4 - len(dig)) + dig
        #plt.savefig(name)
        #plt.clf()
        #self.napply += 1
        ##plt.show()