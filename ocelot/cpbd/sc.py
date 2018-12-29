"""
@author: Igor Zagorodnov @ Martin Dohlus
Created on 27.03.2015
Revision on 01.06.2017: coordinate transform to the velocity direction
"""

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
import logging

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
    def __init__(self, step=1):
        PhysProc.__init__(self)
        self.step = step # in unit step
        self.nmesh_xyz = [63, 63, 63]
        self.low_order_kick = True

        self.start_elem = None
        self.end_elem = None
        self.debug = False
        self.random_mesh = False  # random mesh if True
        self.random_seed = 10     # random seeding number. if None seeding is random

    def prepare(self, lat):
        if self.random_seed != None:
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
            nthread = multiprocessing.cpu_count()
            K2_fft = pyfftw.builders.fftn(K2, axes=None, overwrite_input=False, planner_effort='FFTW_ESTIMATE',
                                       threads=nthread, auto_align_input=False, auto_contiguous=False, avoid_copy=True)
            out_fft = pyfftw.builders.fftn(out, axes=None, overwrite_input=False, planner_effort='FFTW_ESTIMATE',
                                          threads=nthread, auto_align_input=False, auto_contiguous=False, avoid_copy=True)
            out_ifft = pyfftw.builders.ifftn(out_fft()*K2_fft(), axes=None, overwrite_input=False, planner_effort='FFTW_ESTIMATE',
                                          threads=nthread, auto_align_input=False, auto_contiguous=False, avoid_copy=True)
            out = np.real(out_ifft())
        else:
            out = np.real(ifftn(fftn(out)*fftn(K2)))
        t1 = time.time()
        logger.debug('fft time:' + str(t1-t0) + ' sec')
        out[:Nx, :Ny, :Nz] = out[:Nx,:Ny,:Nz]/(4*pi*epsilon_0*hx*hy*hz)
        return out[:Nx, :Ny, :Nz]

    def el_field(self, X, Q, gamma, nxyz):
        N = X.shape[0]
        X[:, 2] = X[:, 2]*gamma
        X_min = np.amin(X, axis=0)
        XX = np.amax(X, axis=0) - X_min
        if self.random_mesh:
            XX = XX*np.random.uniform(low=1, high=1.1)
        logger.debug('mesh steps:' + str(XX))
        # here we use a fast 3D "near-point" interpolation
        # we need a stand-alone module with 1D, 2D, 3D particles-to-grid functions
        steps = XX/(nxyz-3)
        X = X/steps
        X_min = X_min/steps
        X_mid = np.dot(Q, X)/np.sum(Q)
        X_off = np.floor(X_min - X_mid) + X_mid
        X = X - X_off
        nx = nxyz[0]
        ny = nxyz[1]
        nz = nxyz[2]
        nzny = nz*ny
        Xi = np.int_(np.floor(X)+1)
        inds = np.int_(Xi[:, 0]*nzny + Xi[:, 1]*nz + Xi[:, 2])  # 3d -> 1d
        q = np.bincount(inds, Q, nzny*nx).reshape(nxyz)
        p = self.potential(q, steps)
        Ex = np.zeros(p.shape)
        Ey = np.zeros(p.shape)
        Ez = np.zeros(p.shape)
        Ex[:nx-1, :, :] = (p[:nx-1, :, :] - p[1:nx, :, :])/steps[0]
        Ey[:, :ny-1, :] = (p[:, :ny-1, :] - p[:, 1:ny, :])/steps[1]
        Ez[:, :, :nz-1] = (p[:, :, :nz-1] - p[:, :, 1:nz])/steps[2]
        Exyz = np.zeros((N, 3))
        Exyz[:, 0] = ndimage.map_coordinates(Ex, np.c_[X[:, 0], X[:, 1]+0.5, X[:, 2]+0.5].T, order=1)*gamma
        Exyz[:, 1] = ndimage.map_coordinates(Ey, np.c_[X[:, 0]+0.5, X[:, 1], X[:, 2]+0.5].T, order=1)*gamma
        Exyz[:, 2] = ndimage.map_coordinates(Ez, np.c_[X[:, 0]+0.5, X[:, 1]+0.5, X[:, 2]].T, order=1)
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


class LSC(PhysProc):
    def __init__(self, step=1):
        PhysProc.__init__(self, step)
        self.smooth_param = 0.1
        self.step_profile = False
        self.napply = 0

    def imp_lsc(self, gamma, sigma, w, dz):
        """
        gamma - energy
        sigma - transverse RMS size of the beam
        w - omega = 2*pi*f
        """
        indx = np.where(w < 1e-7)[0]
        w[indx] = 1e-7
        alpha = w * sigma / (gamma * speed_of_light)
        alpha2 = alpha * alpha
        ksi = np.exp(alpha2 + np.log(exp1(alpha2)) - 2*np.log(gamma))
        Z = 1j*Z0 / (4 * pi * speed_of_light) * w * ksi * dz

        Z[indx] = 0
        return Z * 1e-12  # --> V/pC

    def imp_step_lsc(self, gamma, rb, w, dz):
        """
        longitudinal space-charge impedance in case of a stepped profile bunch

        gamma - energy
        rb - transverse radius of the beam
        w - omega = 2*pi*f
        """
        indx = np.where(w < 1e-7)[0]
        w[indx] = 1e-7
        x = w*rb/(speed_of_light*gamma)
        Z = 1j*Z0 * speed_of_light / (4 * w * rb*rb) * dz * (1 - x*k1(x))
        Z[indx] = 0
        return Z * 1e-12  # --> V/pC

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

    def wake_lsc(self, s, bunch, gamma, sigma, dz):
        ds = s[1] - s[0]
        dt = ds / speed_of_light
        nb = len(s)
        n = nb * 2
        f = 1 / dt * np.arange(0, n) /n
        if self.step_profile:
            # space charge impedance of transverse step profile
            Za = self.imp_step_lsc(gamma, rb=sigma, w=f[0:nb] * 2 * np.pi, dz=dz)
        else:
            Za = self.imp_lsc(gamma, sigma, w=f[0:nb] * 2 * np.pi, dz=dz)

        sb1 = s[0] + np.cumsum(np.ones(n)) * ds

        bunch1 = np.append(bunch, np.zeros(nb))

        f, Zb = self.wake2impedance(sb1, bunch1 * speed_of_light)

        Z = np.zeros(n, dtype=np.complex)
        Z[0:nb] = Za * Zb[0:nb]
        Z[nb:n] = np.flipud(np.conj(Z[0:nb]))

        xa, wa = self.impedance2wake(f, Z)
        res = -wa[0:nb]
        return res

    def apply(self, p_array, dz):
        """
        wakes in V/pC

        :param p_array:
        :param dz:
        :return:
        """
        if dz < 1e-10:
            logger.debug(" LSC applied, dz < 1e-10, dz = " + str(dz))
            return
        logger.debug(" LSC applied, dz =" + str(dz))
        mean_b = np.mean(p_array.tau())
        sigma_tau = np.std(p_array.tau())
        slice_min = mean_b - sigma_tau / 10
        slice_max = mean_b + sigma_tau / 10
        indx = np.where(np.logical_and(np.greater_equal(p_array.tau(), slice_min), np.less(p_array.tau(), slice_max)))

        if self.step_profile:
            rb = min(np.max(p_array.x()[indx]) - np.min(p_array.x()[indx]),
                     np.max(p_array.y()[indx]) - np.min(p_array.y()[indx]))/2
            sigma = rb
        else:
            sigma = min(np.std(p_array.x()[indx]), np.std(p_array.y()[indx]))
        q = np.sum(p_array.q_array)
        gamma = p_array.E / m_e_GeV
        v = np.sqrt(1 - 1 / gamma ** 2) * speed_of_light
        B = s_to_cur(p_array.tau(), sigma_tau * self.smooth_param, q, v)
        bunch = B[:, 1] / (q * speed_of_light)
        x = B[:, 0]

        W = - self.wake_lsc(x, bunch, gamma, sigma, dz) * 1e12 * q

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