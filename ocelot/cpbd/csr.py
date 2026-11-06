"""
@ authors Martin Dohlus DESY, 2015, Sergey Tomin XFEL, 2016
"""

import copy
import importlib
import logging
import time

import numpy as np

from ocelot.common import math_op
from ocelot.common.globals import m_e_GeV, pi, speed_of_light
from ocelot.cpbd.beam import s_to_cur
from ocelot.cpbd.elements.bend import Bend
from ocelot.cpbd.elements.rbend import RBend
from ocelot.cpbd.elements.sbend import SBend
from ocelot.cpbd.elements.undulator import Undulator
from ocelot.cpbd.elements.xyquadruple import XYQuadrupole
from ocelot.cpbd.high_order import arcline, rk_track_in_field, _rk_reference
from ocelot.cpbd.physics_proc import PhysProc


# matplotlib may or may not be on the HPC nodes at DESY.
#try:
#    import matplotlib.pyplot as plt
#except ImportError:
#    pass

# Try to import numba, pyfftw and numexpr for improved performance
logger = logging.getLogger(__name__)

try:
    import numba as nb
    nb_flag = True
except:
    logger.info("csr.py: module NUMBA is not installed."
                " Install it to speed up calculation")
    nb_flag = False

try:

    from pyfftw.interfaces.numpy_fft import fft, ifft
except:
    logger.info("csr.py: module PYFFTW is not installed."
                " Install it to speed up calculation.")
    from numpy.fft import fft, ifft

try:
    import numexpr as ne
    ne_flag = True
except:
    logger.info("csr.py: module NUMEXPR is not installed."
                " Install it to speed up calculation")
    ne_flag = False


def nextpow2(p):
    """
    Mimic matlab function nextpow2

    :param p:
    :return:
    """
    return int(np.ceil(np.log2(p)))


def csr_convolution(a, b):
    P = len(a)
    Q = len(b)
    L = P + Q - 1
    K = 2 ** nextpow2(L)
    a_pad = np.pad(a, (0, K - P), 'constant', constant_values=(0))
    b_pad = np.pad(b, (0, K - Q), 'constant', constant_values=(0))
    c = ifft(fft(a_pad)*fft(b_pad))
    c = c[:L].real
    return c


def sample_0(i, a, b):
    y = max(0, min(i + 0.5, b) - max(i - 0.5, a))
    return y


def sample_1(i, a, b, c):
    y = 0
    x1 = max(i - 0.5, a) - a
    x2 = min(i + 0.5, b) - a
    if x2 > x1:
        y = (x2 - x1) * (x1 + x2) / (2 * (b - a))

    x2 = c - max(i - 0.5, b)
    x1 = c - min(i + 0.5, c)
    if x2 > x1:
        y = y + (x2 - x1) * (x1 + x2) / (2 * (b - a))
    return y


class CSRConfigurationError(RuntimeError):
    pass


def _validate_binning(charge_binning_fraction, n_bin, m_bin):
    if not np.isfinite(charge_binning_fraction) or not 0 <= charge_binning_fraction <= 1:
        raise CSRConfigurationError("charge_binning_fraction must be between 0 and 1")
    for name, value in (("n_bin", n_bin), ("m_bin", m_bin)):
        if not isinstance(value, (int, np.integer)) or value < 1:
            raise CSRConfigurationError(f"{name} must be a positive integer")


class Smoothing:
    def __init__(self):
        self.print_log = False
        if nb_flag:
            logger.debug("Smoothing: NUMBA")
            self.q_per_step_ip2 = nb.jit(nopython=True)(self.q_per_step_ip2_py)
        else:
            logger.debug("Smoothing: Python")
            self.q_per_step_ip2 = self.q_per_step_ip2_py

    @staticmethod
    def q_per_step_ip2_py(N_BIN, Q_BIN, BIN0, BIN1, NSIG, RMS, step, Nz, z1):
        Nz = int(Nz)
        charge_per_step = np.zeros(Nz)
        weights = np.empty(Nz)

        for nb in range(N_BIN):
            qbin = Q_BIN[nb]
            aa = BIN0[nb] - z1
            bb = BIN1[nb] - z1
            mitte = 0.5 * (aa + bb)
            sigma = RMS[nb]
            aa = mitte - NSIG * sigma
            bb = mitte + NSIG * sigma
            a = aa / step + 1
            k1 = int(min(Nz, max(1, np.floor(a))))
            b = bb / step + 1
            k2 = int(min(Nz, max(1, np.ceil(b))))
            # Normalize the truncated, sampled Gaussian so that each window
            # deposits its full charge, including on a quantized mesh.
            normalization = 0.
            for k in range(k1, k2 + 1):
                xx = (k - 1) * step
                weights[k - 1] = np.exp(-0.5 * ((xx - mitte) / sigma) ** 2)
                normalization += weights[k - 1]
            if normalization == 0:
                raise ValueError("CSR mesh does not resolve a Gaussian bin; reduce mesh_step_quantum")
            for k in range(k1, k2 + 1):
                charge_per_step[k - 1] += weights[k - 1] * (qbin / normalization)
        return charge_per_step

    def Q2EQUI(self, q, BS_params, SBINB, NBIN):
        """
        input
        BIN = bin boundaries BIN(N_BIN, 2), in time or space
        Q_BIN = charges per bin Q_BIN(N_BIN)
        BS_params = binning and smoothing parameters
        binning.......................
        X_QBIN = length or charge binning
        0... 1 = length...charge
        N_BIN = number of bins
        M_BIN = multiple binning(with shifted bins)
        smoothing.....................
        IP_method = 0 / 1 / 2 for rectangular / triangular / gauss
        SP = Gaussian RMS width divided by the full bin width
        sigma_min = minimal sigma, if IP_method == 2
        mesh_step_quantum = optional seventh entry, in metres. Positive values
            round the automatic mesh spacing to the nearest positive integer
            multiple of this value and align the mesh bounds to that spacing.
            Nonpositive values leave the spacing automatic.
        output
        z1, z2, Nz = equidistant mesh(Nz meshlines)
        charge_per_step = charge per step, charge_per_step(1:Nz)
        bins might overlapp!
        """

        _validate_binning(*BS_params[:3])
        N_BIN = BS_params[1]
        M_BIN = BS_params[2]
        K_BIN = N_BIN * M_BIN  # number of sub - bins
        I_BIN = K_BIN - (M_BIN - 1)  # number of bin intervalls
        # put charges to sub - bins
        counts = np.asarray(NBIN, dtype=np.intp)
        if np.size(q) == 1:
            Q_BIN = np.multiply(q, counts, dtype=float)
        else:
            ends = np.cumsum(counts)
            occupied = counts > 0
            starts = ends - counts
            Q_BIN = np.zeros(K_BIN)
            # Exclude empty bins: reduceat treats repeated indices as singletons.
            Q_BIN[occupied] = np.add.reduceat(q[:ends[-1]], starts[occupied])

        # Share each sub-bin's charge between the windows that contain it.
        # Interior sub-bins occur in M_BIN windows; edge sub-bins occur in fewer.
        coverage = np.minimum(np.arange(1, K_BIN + 1), np.arange(K_BIN, 0, -1))
        coverage = np.minimum(coverage, min(M_BIN, I_BIN))
        qsum = np.append([0], np.cumsum(Q_BIN / coverage))
        BIN = [SBINB[0:I_BIN], SBINB[M_BIN + np.arange(I_BIN)]]
        Q_BIN = qsum[M_BIN + np.arange(I_BIN)] - qsum[0:I_BIN]

        # ......................................................................

        # interpolation parameters
        IP_method = BS_params[3]
        if IP_method not in (0, 1, 2):
            raise CSRConfigurationError("ip_method must be 0, 1 or 2")
        SP = BS_params[4]
        if SP <= 0:
            SP = 0.5

        sigma_min = max(0, BS_params[5])
        if len(BS_params) < 7:
            mesh_step_quantum = 0
        else:
            mesh_step_quantum = max(0, BS_params[6])

        # define mesh
        N_BIN = len(BIN[0])
        NSIG = 5
        RMS = []
        if IP_method == 1:
            z1 = np.min(2 * BIN[0][:] - BIN[1][:])
            z2 = np.max(2 * BIN[1][:] - BIN[0][:])
            step = 0.5 * min(BIN[1][:] - BIN[:][0])
        elif IP_method == 2:
            NSIG = 5
            #MITTE = 0.5 * (BIN[0][:] + BIN[1][:])
            #RMS = SP * (BIN[1][:] - BIN[0][:])
            MITTE = 0.5 * (BIN[0] + BIN[1])
            RMS = np.maximum(SP * (BIN[1] - BIN[0]), sigma_min)
            z1 = np.min(MITTE - NSIG * RMS)
            z2 = np.max(MITTE + NSIG * RMS)
            step = 0.25 * np.min(RMS)
        else:
            z1 = np.min(BIN[0][:])
            z2 = np.max(BIN[1][:])
            step = 0.5 * min(BIN[1][:] - BIN[0][:])
            # Include complete sampling cells at both rectangular edges.
            z1 -= step
            z2 += step

        if not np.isfinite(step) or step <= 0:
            raise CSRConfigurationError(
                "CSR smoothing requires positive bin widths or a positive Gaussian sigma_min"
            )

        if mesh_step_quantum > 0:
            step = mesh_step_quantum * max(1, np.round(step / mesh_step_quantum))
            if IP_method != 2:
                z1 -= step
                z2 += step
            z1 = step * np.floor(z1 / step)
            z2 = step * np.ceil(z2 / step)
            Nz = max(1, int(np.round((z2 - z1) / step)))
        else:
            Nz = max(1, int(np.round((z2 - z1) / step)))
            step = (z2 - z1) / Nz
        charge_per_step = np.zeros(int(Nz))
        if IP_method == 1:
            for nb in range(N_BIN):
                aa = BIN[0][nb] - z1
                bb = BIN[1][nb] - z1
                qps = step * Q_BIN[nb] / (bb - aa)
                a = ((3. * aa - bb) / 2.) / step + 1.
                k1 = int(min(Nz, max(1, np.floor(a))))
                b = ((aa + bb) / 2.) / step + 1.
                c = ((3. * bb - aa) / 2.) / step + 1.
                k2 = int(min(Nz, max(1, np.ceil(c))))
                for k in range(k1, k2 + 1):
                    w = sample_1(k, a, b, c)
                    charge_per_step[k-1] += w * qps

        elif IP_method == 2:
            charge_per_step = self.q_per_step_ip2(
                N_BIN, Q_BIN, BIN[0], BIN[1], NSIG, RMS, step, Nz, z1)
        else:
            for nb in range(N_BIN):
                aa = BIN[0][nb] - z1
                bb = BIN[1][nb] - z1
                qps = step * Q_BIN[nb] / (bb - aa)
                a = aa / step + 1
                k1 = int(min(Nz, max(1, np.floor(a))))
                b = bb / step + 1
                k2 = int(min(Nz, max(1, np.ceil(b))))
                for k in range(k1, k2 + 1):
                    w = sample_0(k, a, b)
                    charge_per_step[k-1] += w * qps
        return z1, z2, Nz, charge_per_step

    def Q2EQUI_modified(self, BS_params, SBINB):
        """
        input
        BIN = bin boundaries BIN(N_BIN, 2), in time or space
        Q_BIN = charges per bin Q_BIN(N_BIN)
        BS_params = binning and smoothing parameters
        binning.......................
        X_QBIN = length or charge binning
        0... 1 = length...charge
        N_BIN = number of bins
        M_BIN = multiple binning(with shifted bins)
        smoothing.....................
        IP_method = 0 / 1 / 2 for rectangular / triangular / gauss
        SP = Gaussian RMS width divided by the full bin width
        sigma_min = minimal sigma, if IP_method == 2
        mesh_step_quantum = optional seventh entry, in metres. Positive values
            round the automatic mesh spacing to the nearest positive integer
            multiple of this value and align the mesh bounds to that spacing.
            Nonpositive values leave the spacing automatic.
        output
        z1, z2, Nz = equidistant mesh(Nz meshlines)
        bins might overlap!
        """

        _validate_binning(*BS_params[:3])
        N_BIN = BS_params[1]
        M_BIN = BS_params[2]
        K_BIN = N_BIN * M_BIN  # number of sub - bins
        I_BIN = K_BIN - (M_BIN - 1)  # number of bin intervalls
        # put charges to sub - bins

        BIN = [SBINB[0:I_BIN], SBINB[M_BIN + np.arange(I_BIN)]]

        # interpolation parameters
        IP_method = BS_params[3]
        if IP_method not in (0, 1, 2):
            raise CSRConfigurationError("ip_method must be 0, 1 or 2")
        SP = BS_params[4]
        if SP <= 0:
            SP = 0.5

        sigma_min = max(0, BS_params[5])
        if len(BS_params) < 7:
            mesh_step_quantum = 0
        else:
            mesh_step_quantum = max(0, BS_params[6])

        # define mesh
        N_BIN = len(BIN[0])

        if IP_method == 1:
            z1 = np.min(2 * BIN[0][:] - BIN[1][:])
            z2 = np.max(2 * BIN[1][:] - BIN[0][:])
            step = 0.5 * min(BIN[1][:] - BIN[:][0])
        elif IP_method == 2:
            NSIG = 5
            MITTE = 0.5 * (BIN[0] + BIN[1])
            RMS = np.maximum(SP * (BIN[1] - BIN[0]), sigma_min)
            z1 = np.min(MITTE - NSIG * RMS)
            z2 = np.max(MITTE + NSIG * RMS)
            step = 0.25 * np.min(RMS)
        else:
            z1 = np.min(BIN[0][:])
            z2 = np.max(BIN[1][:])
            step = 0.5 * min(BIN[1][:] - BIN[0][:])
            z1 -= step
            z2 += step

        if not np.isfinite(step) or step <= 0:
            raise CSRConfigurationError(
                "CSR smoothing requires positive bin widths or a positive Gaussian sigma_min"
            )

        if mesh_step_quantum > 0:
            step = mesh_step_quantum * max(1, np.round(step / mesh_step_quantum))
            if IP_method != 2:
                z1 -= step
                z2 += step
            z1 = step * np.floor(z1 / step)
            z2 = step * np.ceil(z2 / step)
            Nz = max(1, int(np.round((z2 - z1) / step)))
        else:
            Nz = max(1, int(np.round((z2 - z1) / step)))

        return z1, z2, Nz


class SubBinning:
    def __init__(self, charge_binning_fraction, n_bin, m_bin):
        self.charge_binning_fraction = charge_binning_fraction
        self.n_bin = n_bin
        self.m_bin = m_bin
        self.print_log = False
        self.p_per_subbins = self.p_per_subbins_np

    @staticmethod
    def p_per_subbins_py(s, SBINB, K_BIN):
        NBIN = np.zeros(K_BIN)
        ib = 0
        Ns = len(s)
        for n in range(Ns):
            while s[n] >= SBINB[ib + 1] and ib < K_BIN - 1:
                ib = ib + 1
            NBIN[ib] = NBIN[ib] + 1
        return NBIN

    @staticmethod
    def p_per_subbins_np(s, SBINB, K_BIN):
        # Particles are sorted, so find each boundary in the particle array.
        # Internal edges belong to the bin on their right; the last is closed.
        ends = np.empty(K_BIN + 1, dtype=np.intp)
        ends[0] = 0
        ends[-1] = len(s)
        ends[1:-1] = np.searchsorted(s, SBINB[1:K_BIN], side='left')
        return np.diff(ends).astype(float)

    def subbin_bound(self, q, s, charge_binning_fraction, n_bin, m_bin):
        """
        input
            q         = array/scalar with charges of macro particles (sorted)
            s         = longitudinal vector (time or space) (sorted)
            B_params  = binning parameters
                        [X_QBIN,N_BIN,M_BIN]
                        X_QBIN = length or charge binning
                                 0 ... 1 = length ... charge
                        N_BIN  = number of bins
                        M_BIN  = multiple binning (with shifted bins)
        output
            SBINB     = array with subbin boundaries
            NBIN      = particles per subbin
        particles are sorted!
        all particles are valid
        """
        _validate_binning(charge_binning_fraction, n_bin, m_bin)
        if len(s) < 2 or not np.all(np.isfinite(s)) or s[-1] <= s[0]:
            raise CSRConfigurationError(
                "CSR requires at least two particles with a finite, nonzero longitudinal extent"
            )
        X_QBIN = charge_binning_fraction
        N_BIN = n_bin
        M_BIN = m_bin
        K_BIN = N_BIN * M_BIN  # number of sub-bins

        Ns = len(s)
        # binning intervalls
        # aa=monoton "charge" vector

        if X_QBIN == 0:
            aa = (s - s[0]) / (s[-1] - s[0])
        else:
            if np.size(q) != 1:
                # Locate each particle at the midpoint of its cumulative charge.
                # This also makes the binning symmetric under bunch reversal.
                aa = np.cumsum(q) - 0.5 * q
            else:
                aa = q * (np.arange(1, Ns + 1) - 0.5)
            if not np.isfinite(aa[-1] - aa[0]) or aa[-1] == aa[0]:
                raise CSRConfigurationError("Charge binning requires a nonzero finite charge extent")
            if ne_flag:
                aa0 = aa[0]
                aaNs = aa[Ns - 1]
                s0 = s[0]
                sNs = s[Ns - 1]
                aa = ne.evaluate(
                    '(aa - aa0) / (aaNs - aa0) * X_QBIN +'
                    '(s - s0) / (sNs - s0) * (1 - X_QBIN)')
            else:
                aa = (aa - aa[0]) / (aa[Ns - 1] - aa[0])
                bb = (s - s[0]) / (s[Ns - 1] - s[0])
                aa = aa * X_QBIN + bb * (1 - X_QBIN)
        if np.min(np.diff(aa)) == 0:
            aa = 0.999 * aa + 0.001 * (np.arange(0, Ns)) / (Ns - 1)

        # vector with bin boundaries
        SBINB = np.interp(np.arange(K_BIN + 1.) / K_BIN, aa, s)
        SBINB[0] = s[0]
        SBINB[K_BIN] = s[Ns - 1]
        # particles per subbins
        NBIN = self.p_per_subbins(s, SBINB, K_BIN)

        return SBINB, NBIN


class K0_fin_anf:
    def __init__(self):
        self.print_log = False
        if nb_flag:
            logger.debug("K0_fin_anf: NUMBA")
            self.K0_1 = nb.jit(nopython=True)(self.K0_1_jit)
            self.K0_0 = nb.jit(nopython=True)(self.K0_0_jit)
            self.integrate_kernel = nb.njit(self.integrate_kernel_jit)
            self.eval = self.K0_fin_anf_opt
        elif ne_flag:
            logger.debug("K0_fin_anf: NumExpr")
            self.eval = self.K0_fin_anf_numexpr
        else:
            logger.debug("K0_fin_anf: Python")
            self.eval = self.K0_fin_anf_np

    @staticmethod
    def integrate_kernel_np(K, s):
        N = K.shape[0]
        a = np.empty(N, dtype=K.dtype)
        a[:-1] = 0.5 * (K[:-1] + K[1:]) * np.diff(s)
        a[-1] = 0.5 * K[-1] * s[-1]
        KS = np.cumsum(a[::-1])[::-1]
        return KS

    @staticmethod
    def integrate_kernel_jit(K, s):
        N = K.shape[0]
        a = np.empty(N, dtype=K.dtype)
        for k in range(N-1):
            a[k] = 0.5 * (K[k] + K[k+1]) * (s[k+1] - s[k])
        a[N-1] = 0.5 * K[N-1] * s[N-1]
        # Cumulative sum in reverse order
        KS = np.empty(N, dtype=K.dtype)
        total = 0.0
        for k in range(N-1, -1, -1):
            total += a[k]
            KS[k] = total
        return KS

    @staticmethod
    def K0_1_jit(indx, j, R, n, traj4, traj5, traj6, w, gamma):
        g2i = 1. / gamma ** 2
        b2 = 1. - g2i
        beta = np.sqrt(b2)
        K = np.zeros(indx - j)
        t4i = traj4[indx]
        t5i = traj5[indx]
        t6i = traj6[indx]
        for i in range(indx - j):
            Ri_inv = 1/R[i]
            n0i = n[i, 0] * Ri_inv
            n1i = n[i, 1] * Ri_inv
            n2i = n[i, 2] * Ri_inv
            # kernel
            t4 = traj4[i+j]
            t5 = traj5[i+j]
            t6 = traj6[i+j]
            x = n0i * t4 + n1i * t5 + n2i * t6
            K[i] = ((beta * (x - n0i * t4i - n1i * t5i - n2i * t6i) -
                     b2 * (1. - t4 * t4i - t5 * t5i - t6 * t6i) -
                     g2i) * Ri_inv - (1. - beta * x) / w[i] * g2i)
        return K

    @staticmethod
    def K0_0_jit(i, traj0, traj1, traj2, traj3, gamma, s, n, R, w, wmin):
        g2i = 1. / gamma ** 2
        b2 = 1. - g2i
        beta = np.sqrt(b2)
        i_0 = 0
        traj0i = traj0[i]
        traj1i = traj1[i]
        traj2i = traj2[i]
        traj3i = traj3[i]
        for j_i in range(i):
            # Start looping from last element
            j = i - j_i - 1
            s_j = traj0[j] - traj0i
            n1 = traj1i - traj1[j]
            n2 = traj2i - traj2[j]
            n3 = traj3i - traj3[j]
            R_j = np.sqrt(n1 * n1 + n2 * n2 + n3 * n3)
            w_j = s_j + beta * R_j
            s[j] = s_j
            R[j] = R_j
            w[j] = w_j
            n[j, 0] = n1
            n[j, 1] = n2
            n[j, 2] = n3
            if w_j <= wmin:
                i_0 = j
                break
        # return last index where w <= wmin
        return i_0

    def K0_fin_anf_opt(self, i, traj, wmin, gamma):
        s = np.zeros(i)
        n = np.zeros((i, 3))
        R = np.zeros(i)
        w = np.zeros(i)
        j = self.K0_0(
            i, traj[0], traj[1], traj[2], traj[3], gamma, s, n, R, w, wmin)

        s = s[j:]
        n = n[j:]
        R = R[j:]
        w = w[j:]

        K = self.K0_1(i, j, R, n, traj[4], traj[5], traj[6], w, gamma)
        if len(K) > 1:
            KS = self.integrate_kernel(K, s)
        else:
            KS = 0.5 * K[-1] * s[-1]
        return w, KS

    def K0_fin_anf_np(self, i, traj, wmin, gamma):
        g2i = 1. / gamma ** 2
        b2 = 1. - g2i
        beta = np.sqrt(b2)
        i_0 = self.estimate_start_index(i, traj, wmin, beta)
        idx = np.arange(i_0, i)
        s = traj[0, idx] - traj[0, i]
        n = traj[1:4, i, None] - traj[1:4, idx]
        R = np.linalg.norm(n, axis=0)
        w = s + beta * R
        j = np.where(w <= wmin)[0]
        if len(j) > 0:
            j = j[-1]
            idx = idx[j:]
            s = s[j:]
            n = n[:, j:]
            R = R[j:]
            w = w[j:]
        else:
            j = i_0
        R_inv = 1 / R
        n *= R_inv
        t4 = traj[4, idx]
        t5 = traj[5, idx]
        t6 = traj[6, idx]
        t_i = traj[4:7, i]
        x = n[0] * t4 + n[1] * t5 + n[2] * t6
        K = ((beta * (x - (n[0] * t_i[0] + n[1] * t_i[1] + n[2] * t_i[2]))
              - b2 * (1. - (t4 * t_i[0] + t5 * t_i[1] + t6 * t_i[2]))
              - g2i) * R_inv
             - (1. - beta * x) / w * g2i)
        if len(K) > 1:
            KS = self.integrate_kernel_np(K, s)
        else:
            KS = 0.5 * K[-1] * s[-1]
        return w, KS

    def K0_fin_anf_numexpr(self, i, traj, wmin, gamma):
        g2i = 1. / gamma ** 2
        b2 = 1. - g2i
        beta = np.sqrt(b2)
        i_0 = self.estimate_start_index(i, traj, wmin, beta)
        idx = np.arange(i_0, i)
        s = traj[0, idx] - traj[0, i]
        n0 = traj[1, i] - traj[1, idx]
        n1 = traj[2, i] - traj[2, idx]
        n2 = traj[3, i] - traj[3, idx]
        R = ne.evaluate("sqrt(n0**2 + n1**2 + n2**2)")
        w = ne.evaluate('s + beta*R')
        j = np.where(w <= wmin)[0]
        if len(j) > 0:
            j = j[-1]
            idx = idx[j:]
            s = s[j:]
            n0 = n0[j:]
            n1 = n1[j:]
            n2 = n2[j:]
            R = R[j:]
            w = w[j:]
        else:
            j = i_0
        R_inv = 1 / R
        n0 = n0 * R_inv
        n1 = n1 * R_inv
        n2 = n2 * R_inv
        t4 = traj[4, idx]
        t5 = traj[5, idx]
        t6 = traj[6, idx]
        t4i = traj[4, i]
        t5i = traj[5, i]
        t6i = traj[6, i]
        x = ne.evaluate('n0*t4 + n1*t5 + n2*t6')
        K = ne.evaluate(
            '((beta*(x - n0*t4i - n1*t5i - n2*t6i) -'
            '  b2*(1. - t4*t4i - t5*t5i - t6*t6i) - g2i)/R -'
            ' (1. - beta*x)/w*g2i)')
        if len(K) > 1:
            KS = self.integrate_kernel_np(K, s)
        else:
            KS = 0.5 * K[-1] * s[-1]
        return w, KS

    def estimate_start_index(self, i, traj, w_min, beta, i_min=1000, n_test=10):
        """
        This method estimates the index of the first trajectory point from
        which CSR effects should be computed.

        This method can significantly reduce the computing time of CSR effects
        by pre-discaring regions of the reference trajectory which do not
        influence the CSR calculation (i.e. the points where w <= wmin). This
        is performed by testing the w <= wmin condition for a subset of n_test
        equally-spaced points along the reference trajectory. The index of the
        last tested trajectory point in which w <= wmin is returned.

        Parameters
        ----------

        i : int
            Iteration index

        traj : ndarray
            Reference trajectory along which CSR forces are calculated

        w_min : float
            Leftmost edge of the longitudinal bunch binning.

        beta : float
            Relativistic factor.

        i_min : int
            Minimum iteration index. When i<i_min, no estimation of the starting
            index is performed (0 is returned).

        n_test : int
            Number of points along the trajectory in which to test whether they
            should be taken into account for the CSR calculation.

        Returns
        -------
        The estimated start index, which is always <= than the real one.

        """
        i_0 = 0
        if i > i_min:
            idx = np.linspace(0, i, n_test, dtype=np.int32)
            s = traj[0, idx] - traj[0, i]
            n0 = traj[1, i] - traj[1, idx]
            n1 = traj[2, i] - traj[2, idx]
            n2 = traj[3, i] - traj[3, idx]
            R = np.sqrt(n0*n0 + n1*n1 + n2*n2)
            w = s + beta * R
            j = np.where(w <= w_min)[0]
            if len(j) > 0:
                i_0 = idx[j[-1]]
        return i_0


class CSR(PhysProc):
    """
    This class simulates the CSR wakefield, its interaction with the beam, and applies the corresponding CSR kick to the particle array.

    Attributes:
        step (int): Step size for applying the CSR kick to the beam (ParticleArray). Default is 1 [in Navigator.unit_step].
        sigma_min (float): Minimal sigma used in Gaussian filtering. Default is 1.e-4.
        traj_step (float): Trajectory step, or integration step for CSR-wake calculation. Default is 0.0002 [m].
        apply_step (float): Spacing of kernel evaluations within a kick, in metres.
            Zero samples every trajectory point. Default is 0.0005.

        charge_binning_fraction (float): Blend between length and charge binning,
            from 0 (equal longitudinal widths) to 1 (approximately equal bin
            charges). Intermediate values mix normalized length and cumulative
            charge coordinates. Default is 0.
        n_bin (int): Number of bins used for binning. Default is 100.
        m_bin (int): Multiple binning factor with shifted bins. Default is 5.

        ip_method (int): Method for smoothing (0: rectangular, 1: triangular, 2: Gaussian). Default is 2 (Gaussian).
        sp (float): Gaussian RMS width divided by the full bin width, subject to
            sigma_min. Nonpositive values use 0.5. Default is 0.5.
        mesh_step_quantum (float): Spacing increment for the longitudinal density
            and wake mesh, in metres. If positive, replace the automatically
            estimated spacing h by u * max(1, round(h / u)), where u is this
            value, and align the mesh bounds to the resulting spacing. For
            example, u=2e-6 rounds h=5.3e-6 to 6e-6 m. Nonpositive values keep
            automatic spacing. This does not set the trajectory or kick step.
            Default is 0.

        energy (float or None): Reference energy in GeV, required for an RK
            trajectory. CSR kicks always use the finite energy p_array.E.

        z_csr_start (float): Position of the start element in the trajectory in [m]. Default is 0.
        z0 (float): Position in the navigator for the starting point in the track. Default is 0.

        end_poles (bool): If True, magnetic field configuration will be applied (1/4, -3/4, 1, ...). Default is False.
        rk_traj (bool): If True, the trajectory of the reference particle is calculated using the Runge-Kutta method. Default is False.

        pict_debug (bool): If True, the trajectory of the reference particle and CSR wakes will be saved for each step in the working folder. Default is False.

        sub_bin (SubBinning): An instance of the SubBinning class for particle binning.
        bin_smoth (Smoothing): An instance of the Smoothing class for applying smoothing to the CSR wakefield.
        k0_fin_anf (K0_fin_anf): An instance of the K0_fin_anf class for further CSR-related calculations.
    """

    _renamed_parameters = {"x_qbin": "charge_binning_fraction", "step_unit": "mesh_step_quantum"}

    def __init__(self, **kw):
        removed = sorted({"filter_order", "n_mesh", "debug"}.intersection(kw))
        if removed:
            raise TypeError("Unused CSR parameters have been removed: " + ", ".join(removed))
        for old, new in self._renamed_parameters.items():
            if old in kw:
                raise TypeError(f"CSR parameter {old!r} was renamed to {new!r}; use the new name")
        PhysProc.__init__(self)
        # binning parameters
        self.charge_binning_fraction = kw.get("charge_binning_fraction", 0)
        self.n_bin = kw.get("n_bin", 100) # number of bins
        self.m_bin = kw.get("m_bin", 5)   # multiple binning(with shifted bins)

        # smoothing
        self.ip_method = kw.get("ip_method", 2)     # = 0 / 1 / 2 for rectangular / triangular / gauss
        self.sp = kw.get("sp", 0.5)                 # Gaussian RMS / bin width
        self.sigma_min = kw.get("sigma_min", 1.e-4) # minimal sigma, if ip_method == 2
        self.mesh_step_quantum = kw.get("mesh_step_quantum", 0)

        # trajectory
        self.traj_step = kw.get("traj_step", 0.0002)  # [m] step of the trajectory
        self.energy = kw.get("energy", None)          # [GeV], required for RK trajectory

        # CSR kick
        self.apply_step = kw.get("apply_step", 0.0005) # [m] step of the calculation CSR kick: csr_kick += csr(apply_step)
        self.step = kw.get("step", 1)                  # step in the unit steps, step_in_[m] = self.step * navigator.unit_step [m].
                                                       # The CSR kick is applied at the end of the each step

        self.z_csr_start = kw.get("z_csr_start", 0.) # z [m] position of the start_elem
        self.z0 = kw.get("z0", 0.)                   # self.z0 = navigator.z0 in track.track()

        self.end_poles = kw.get("end_poles", False) # if True magnetic field configuration 1/4, -3/4, 1, ...
        self.rk_traj = kw.get("rk_traj", False)     # calculate trajectory of the reference particle with RK method

        # if True trajectory of the reference particle will be produced
        # and CSR wakes will be saved in the working folder on each step
        self.pict_debug = kw.get("pict_debug", False)

        self.sub_bin = SubBinning(charge_binning_fraction=self.charge_binning_fraction, n_bin=self.n_bin, m_bin=self.m_bin)
        self.bin_smoth = Smoothing()
        self.k0_fin_anf = K0_fin_anf()

    def __setattr__(self, name, value):
        replacement = self._renamed_parameters.get(name)
        if replacement is not None:
            raise AttributeError(f"CSR parameter {name!r} was renamed to {replacement!r}; use the new name")
        super().__setattr__(name, value)

    def K0_inf_anf(self, i, traj, wmin):
        """

        :param i: index of the trajectories points for the convolution kernel is calculated;
        :param traj: trajectory. traj[0,:] - longitudinal coordinate,
                                 traj[1,:], traj[2,:], traj[3,:] - rectangular coordinates, \
                                 traj[4,:], traj[5,:], traj[6,:] - tangential unit vectors
        :param wmin: the first coordinate of the mash
        :return:
        """

        i1 = i-1  # ignore points i1+1:i on linear path to observer
        ra = np.arange(0, i1+1)
        s = traj[0, ra] - traj[0, i]
        n = np.array([traj[1, i] - traj[1, ra],
                      traj[2, i] - traj[2, ra],
                      traj[3, i] - traj[3, ra]])

        R = np.sqrt(np.sum(n**2, axis=0))
        n = np.array([n[0, :]/R, n[1, :]/R, n[2, :]/R])

        w = s + R
        j = np.where(w <= wmin)[0]
        if len(j) > 0:
            j = j[-1]
            ra = np.arange(j, i1+1)
            w = w[ra]
            s = s[ra]

        # kernel
        K = (n[0, ra]*(traj[4, ra] - traj[4, i]) +
             n[1, ra]*(traj[5, ra] - traj[5, i]) +
             n[2, ra]*(traj[6, ra] - traj[6, i]) -
             (1. - traj[4, ra]*traj[4, i] -
              traj[5, ra]*traj[5, i] -
              traj[6, ra]*traj[6, i])) / R[ra]

        # integrated kernel: KS=int_s^0{K(u)*du}=int_0^{-s}{K(-u)*du}
        if np.shape(K)[0] > 1:
            a = np.append(0.5*(K[0:-1] + K[1:])*np.diff(s), 0.5*K[-1]*s[-1])
            KS = np.cumsum(a[::-1])[::-1]
            #KS = np.fliplr(np.cumsum(np.fliplr([0.5*(K[1:-1] + K[2:])*np.diff(s), 0.5*K[-1]*s[-1]])))
        else:
            KS = 0.5*K[-1]*s[-1]

        return w, KS

    def K0_fin_inf(self, i, traj, w_range, gamma):
        """
        Radiative interaction is coming from the infinite straight line
        that is assumed before the specified CSR region and it is calculated analytically

        :param i:
        :param traj:
        :param w_range:
        :param gamma:
        :return:
        """

        g2 = gamma**2
        g2i = 1./g2
        b2 = 1. - g2i
        beta = np.sqrt(b2)
        # winf
        Rv1 = traj[1:4, i] - traj[1:4, 0]
        s1 = traj[0, 0] - traj[0, i]
        ev1 = traj[4:, 0]
        evo = traj[4:, i]
        winfms1 = np.dot(Rv1, ev1)

        aup = -Rv1 + winfms1*ev1
        a2 = np.dot(aup, aup)
        a = np.sqrt(a2)

        uup = aup/a if a != 0 else None

        winf = s1 + winfms1
        s = winf + gamma*(gamma*(w_range - winf)-beta*np.sqrt(g2*(w_range-winf)**2+a2))
        R = (w_range-s)/beta

        if a2/R[0]**2 > 1e-7:
            KS = (beta*(1. - np.dot(ev1, evo))*np.log(R[0]/R) - beta*np.dot(uup, evo)*(np.arctan((s[0] - winf)/a) - np.arctan((s-winf)/a))
                  - (b2*np.dot(ev1, evo) - 1)*np.log((winf - s + R)/(winf-s[0] + R[0]))
                  + g2i*np.log(w_range[0]/w_range))

        else:
            KS = (beta*(1. - np.dot(ev1, evo))*np.log(R[0]/R)
                  - (b2*np.dot(ev1, evo) - 1.)*np.log((winf - s + R)/(winf - s[0] + R[0]))
                  + g2i*np.log(w_range[0]/w_range))
        return KS

    def K0_inf_inf(self, i, traj, w_range):
        """
        Radiative interaction is coming from the infinite straight line
        that is assumed before the specified CSR region and it is calculated analytically

        :param i:
        :param traj:
        :param w_range:
        :return:
        """
        Rv1 = traj[1:4, i] - traj[1:4, 0]
        s1 = traj[0, 0] - traj[0, i]
        ev1 = traj[4:, 0]
        evo = traj[4:, i]
        winfms1 = np.dot(Rv1, ev1)
        aup = -Rv1 + winfms1*ev1
        a2 = np.dot(aup, aup)
        a = np.sqrt(a2)
        uup = aup/a
        winf = s1 + winfms1

        Nvalid = np.where(winf < w_range)[0]
        if len(Nvalid) > 0:
            Nvalid = Nvalid[0]
            w = w_range[Nvalid:]
            s = (winf+w)/2. + a2/2./(winf-w)
            KS = (1. - np.dot(ev1, evo))*np.log(0.5*a2/(w-winf)/(w-s)) + np.dot(uup, evo)*(np.arctan((s-winf)/a) + np.pi/2.)
            KS = np.append(np.zeros(Nvalid), KS)
        else:
            KS = np.zeros(len(w_range))
        return KS

    def CSR_K1(self, i, traj, NdW, gamma=None):
        """
        :param i: index of the trajectories points for the convolution kernel is calculated;
        :param traj: trajectory. traj[0,:] - longitudinal coordinate,
                                 traj[1,:], traj[2,:], traj[3,:] - rectangular coordinates, \
                                 traj[4,:], traj[5,:], traj[6,:] - tangential unit vectors
        :param NdW: list [N, dW], NdW[0] is number of mesh points, NdW[1] = dW > 0 is increment.
                                  Mesh = Mesh = (-N:0) * dW
        :param gamma:
        :return:
        """

        if i == 0:
            return np.zeros(int(NdW[0]) + 1)
        if gamma is not None:
            L_fin = True
        else:
            L_fin = False

        w_range = np.arange(-NdW[0]-1, 0)*NdW[1]
        if L_fin:
            w, KS = self.k0_fin_anf.eval(i, traj, w_range[0], gamma)
        else:
            w, KS = self.K0_inf_anf(i, traj, w_range[0])

        # Cancellation in s + beta*R can make the retarded grid repeat or
        # decrease. Interpolation requires a strictly increasing grid.
        KS = np.atleast_1d(KS)
        if np.any(w[1:] <= w[:-1]):
            w, indices = np.unique(w, return_index=True)
            KS = KS[indices]
        KS1 = KS[0]
        w_min = w[0]
        if w_range[0] < w_min:
            m = np.where(w_range < w_min)[0][-1]
            if L_fin:
                KS2 = self.K0_fin_inf(i, traj, np.append(w_range[0:m+1], w_min), gamma)
            else:
                KS2 = self.K0_inf_inf(i, traj, np.append(w_range[0:m+1], w_min))

            KS2 = (KS2[-1] - KS2) + KS1
            KS = np.append(KS2[0:-1], np.interp(w_range[m+1:], w, KS))

        else:
            KS = np.interp(w_range, w, KS)
        four_pi_eps0 = 1./(1e-7*speed_of_light**2)
        K1 = np.diff(KS, n=2, append=[0., 0.]) / NdW[1] / four_pi_eps0

        return K1

    def _undulator_trajectory(self, elem):
        """Sample a field orbit to the local Cartesian exit plane Z=elem.l."""
        from scipy.interpolate import CubicHermiteSpline

        field = elem.element.get_csr_mag_field(self.energy, use_end_poles=self.end_poles)
        count = max(1, int(np.ceil(elem.l / self.traj_step)))
        tracked, paths = rk_track_in_field(np.zeros((6, 1)), s_stop=elem.l, N=count + 1,
                                           energy=self.energy, mag_field=field, return_path=True)
        tracked = tracked.reshape(-1, 9)
        positions = tracked[:, [0, 2, 4]]
        tangents = np.column_stack((tracked[:, 1], tracked[:, 3], np.ones(count + 1)))
        tangents /= np.linalg.norm(tangents, axis=1)[:, None]
        path = paths[:, 0]
        count = max(1, int(np.ceil(path[-1] / self.traj_step)))
        arc = np.linspace(0., path[-1], count + 1)
        # Positions and their derivatives describe one smooth curve. Linear
        # interpolation would introduce chord errors into the retarded kernel.
        curve = CubicHermiteSpline(path, positions, tangents, axis=0)
        positions = curve(arc)
        tangents = curve(arc, 1)
        tangents /= np.linalg.norm(tangents, axis=1)[:, None]
        return arc, positions, tangents, positions[:, 2]

    def _lattice_to_trajectory_s(self, lattice_s):
        """Convert distance from the CSR entrance to physical orbit arc length."""
        return np.interp(lattice_s, self._trajectory_lattice_s, self.csr_traj[0])

    def prepare(self, lat):
        """
        Build the reference orbit used by the CSR kernel in one Cartesian frame.

        The origin and +Z tangent are at the CSR entrance. X/Y are the initial
        horizontal/vertical axes; subsequent elements use transported local
        axes. Positive untilted bend angles turn towards -X. Row 0 is physical
        arc length in metres, not Cartesian Z or Navigator lattice distance.
        The latter differ inside undulators and are mapped separately.

        Args:
            lat (MagneticLattice): The magnetic lattice for which the calculations are performed.

        Returns:
            ndarray: The CSR trajectory with the following components:
                - `traj[0, :]`: Physical arc length from the CSR entrance.
                - `traj[1, :]`, `traj[2, :]`, `traj[3, :]`: Rectangular coordinates.
                - `traj[4, :]`, `traj[5, :]`, `traj[6, :]`: Tangential unit vectors.
        """
        self.check_step()
        _validate_binning(self.charge_binning_fraction, self.n_bin, self.m_bin)
        if not np.isfinite(self.traj_step) or self.traj_step <= 0:
            raise CSRConfigurationError("traj_step must be positive and finite")
        if not np.isfinite(self.apply_step) or self.apply_step < 0:
            raise CSRConfigurationError("apply_step must be nonnegative and finite")

        # if pict_debug = True import matplotlib
        if self.pict_debug:
            self.plt = importlib.import_module("matplotlib.pyplot")
            self.napply = 0
            self.total_wake = 0


        self.z_csr_start = sum(p.l for p in lat.sequence[:self.indx0])
        # Navigator stops at the entrance of indx1. Preserve the historical
        # inclusive range for direct prepare() calls that only set indices.
        stop = self.indx1 if self.s_stop is not None else self.indx1 + 1
        elements = lat.sequence[self.indx0:stop]
        if any(isinstance(elem, Undulator) for elem in elements) and not self.rk_traj:
            self.rk_traj = True
            logger.warning("CSR: Undulator element is in CSR section --> rk_traj = True")
        if self.rk_traj and (self.energy is None or not np.isfinite(self.energy)
                             or self.energy <= m_e_GeV):
            raise CSRConfigurationError("RK trajectory requires CSR.energy above the electron rest energy")

        entrance = np.array([0., 0., 0., 0., 0., 0., 1.])
        segments = [entrance[:, None]]
        lattice_samples = [np.array([0.])]
        lattice_s = 0.
        trajectory_frame = np.eye(3)
        for elem in elements:
            if not np.isfinite(elem.l) or elem.l < 0:
                raise CSRConfigurationError("CSR element lengths must be nonnegative and finite")
            if elem.l == 0:
                continue
            count = max(1, int(np.round(elem.l / self.traj_step)))
            arc = np.linspace(0., elem.l, count + 1)
            local_lattice_s = arc
            exit_frame = np.eye(3)
            if isinstance(elem, Undulator):
                arc, positions, tangents, local_lattice_s = self._undulator_trajectory(elem)
                # Parallel transport the transverse axes without imposing a
                # rotation about the tangent (the same Bishop convention as RK).
                for direction in tangents[1:]:
                    axis = np.cross(exit_frame[:, 2], direction)
                    cosine = np.dot(exit_frame[:, 2], direction)
                    cross = np.array([[0., -axis[2], axis[1]],
                                      [axis[2], 0., -axis[0]],
                                      [-axis[1], axis[0], 0.]])
                    exit_frame = (np.eye(3) + cross + cross @ cross / (1. + cosine)) @ exit_frame
            elif self.rk_traj and hasattr(elem.element, "get_mag_field"):
                # Fields are evaluated in this element's entrance frame. Arc
                # integration reaches exactly elem.l, independent of bend angle.
                field = elem.element.get_mag_field(self.energy)
                positions, frames = _rk_reference(arc, self.traj_step, self.energy, field,
                                                  require_forward_z=False)
                tangents = frames[:, :, 2]
                exit_frame = frames[-1]
            elif isinstance(elem, XYQuadrupole):
                positions, frames = elem.element._reference_geometry(arc)
                tangents = frames[:, :, 2]
                exit_frame = frames[-1]
            else:
                radius_vector = np.zeros(3)
                if isinstance(elem, (Bend, RBend, SBend)) and elem.angle != 0:
                    radius = -elem.l / elem.angle
                    radius_vector = radius * np.array([-np.sin(elem.tilt), np.cos(elem.tilt), 0.])
                    axis = radius_vector / np.linalg.norm(radius_vector)
                    cross = np.array([[0., -axis[2], axis[1]],
                                      [axis[2], 0., -axis[0]],
                                      [-axis[1], axis[0], 0.]])
                    angle = abs(elem.angle)
                    exit_frame = np.eye(3) + np.sin(angle) * cross + 2 * np.sin(angle / 2)**2 * (cross @ cross)
                local = arcline(np.array([[0., 0., 0., 0., 0., 0., 1.]]).T,
                                elem.l, self.traj_step, radius_vector)
                arc, positions, tangents = local[0], local[1:4].T, local[4:7].T
                local_lattice_s = arc

            xyz = entrance[1:4, None] + trajectory_frame @ positions.T
            directions = trajectory_frame @ tangents.T
            directions /= np.linalg.norm(directions, axis=0)
            segment = np.vstack((entrance[0] + arc, xyz, directions))
            segments.append(segment[:, 1:])
            lattice_samples.append(lattice_s + local_lattice_s[1:])
            lattice_s += elem.l
            entrance = segment[:, -1]
            trajectory_frame = trajectory_frame @ exit_frame
        self.csr_traj = np.column_stack(segments)
        # Same samples as csr_traj, but measured in Navigator lattice distance.
        self._trajectory_lattice_s = np.concatenate(lattice_samples)

        # plot trajectory of the refernece particle
        if self.pict_debug:
            fig = self.plt.figure(figsize=(10, 8))
            ax1 = self.plt.subplot(211)
            #ax1.plot(self.csr_traj[0, :], self.csr_traj[0, :] -self.csr_traj[3, :], "r", label="X")
            ax1.plot(self.csr_traj[0, :], self.csr_traj[1, :]*1000, "r", label="X")
            ax1.plot(self.csr_traj[0, :], self.csr_traj[2, :]*1000, "b", label="Y")
            self.plt.legend()
            self.plt.ylabel("X/Y [mm]")
            self.plt.setp(ax1.get_xticklabels(), visible=False)
            ax3 = self.plt.subplot(212, sharex=ax1)
            ax3.plot(self.csr_traj[0, :], self.csr_traj[1 + 3, :]*1000, "r", label=r"$X'$")
            ax3.plot(self.csr_traj[0, :], self.csr_traj[2 + 3, :]*1000, "b", label=r"$Y'$")
            self.plt.legend()
            self.plt.ylabel(r"$X'/Y'$ [mrad]")
            self.plt.xlabel("s [m]")
            self.plt.show()
            # data = np.array([np.array(self.s), np.array(self.total_wake)])
            # np.savetxt("trajectory_cos.txt", self.csr_traj)
        return self.csr_traj

    def apply(self, p_array, delta_s):
        if p_array.n == 0 or delta_s == 0:
            return
        if not np.isfinite(delta_s) or delta_s < 0:
            raise CSRConfigurationError("CSR kick length must be nonnegative and finite")
        if not np.isfinite(p_array.E) or p_array.E <= m_e_GeV:
            raise CSRConfigurationError("CSR requires a beam energy above the electron rest energy")
        lattice_s = self.z0 - self.z_csr_start
        s_prev, s_cur = self._lattice_to_trajectory_s([lattice_s - delta_s, lattice_s])
        path_length = s_cur - s_prev
        z = -p_array.tau()

        ind_z_sort = np.argsort(z)
        z_sort = z[ind_z_sort]
        q_sort = p_array.q_array[ind_z_sort]
        SBINB, NBIN = self.sub_bin.subbin_bound(q_sort, z_sort, self.charge_binning_fraction, self.n_bin, self.m_bin)

        B_params = [self.charge_binning_fraction, self.n_bin, self.m_bin,
                    self.ip_method, self.sp, self.sigma_min, self.mesh_step_quantum]
        s1, s2, Ns, lam_ds = self.bin_smoth.Q2EQUI(q_sort, B_params, SBINB, NBIN)

        st = (s2 - s1) / Ns
        sa = s1 + st / 2.
        Ndw = [Ns - 1, st]

        #print(Ndw, z[ind_z_sort][-1] - z[ind_z_sort][0], len(lam_ds), z[ind_z_sort].max() - z[ind_z_sort].min(), s1, s2)
        #plt.plot(lam_ds)
        #plt.show()
        s_array = self.csr_traj[0, :]
        indx = (np.abs(s_array - s_cur)).argmin()
        indx_prev = (np.abs(s_array - s_prev)).argmin()
        gamma = p_array.E/m_e_GeV
        h = max(1., self.apply_step/self.traj_step)

        itr_ra = np.unique(-np.round(np.arange(-indx, -indx_prev, h))).astype(int)
        if len(itr_ra) == 0:
            # A navigator step can be shorter than the trajectory spacing.
            # Still integrate its full length using the nearest kernel sample.
            itr_ra = np.array([indx])
        #print(s1, s2, Ns, Ndw, itr_ra)
        K1 = 0
        for it in itr_ra:
            K1 += self.CSR_K1(it, self.csr_traj, Ndw, gamma=gamma)
            #plt.plot(np.linspace(s1, s2, num=int(Ns)), K1*st)
            #plt.show()


        #x = np.arange(Ndw[0] + 1) * Ndw[1]
        #xp = np.arange(self.kernel_Ndw[0] + 1) * self.kernel_Ndw[1]
        #for it in itr_ra:
        #    fp = self.kernels[it]
        #    #print(len(fp), len(xp), self.kernel_Ndw)
        #    K1 = np.interp(x, xp, fp)
        #    K1 += K1 #self.CSR_K1(it, self.csr_traj, Ndw, gamma=gamma)
        K1 /= len(itr_ra)

        lam_K1 = csr_convolution(lam_ds, K1[::-1]) / st * path_length

        dE = np.interp(z_sort*(1./st)+(0. - sa/st), np.arange(len(lam_K1)), lam_K1)

        pc_ref = np.sqrt(p_array.E ** 2 / m_e_GeV ** 2 - 1) * m_e_GeV
        delta_p = dE * 1e-9 / pc_ref
        p_array.rparticles[5][ind_z_sort] += delta_p

        if self.pict_debug:
            self.plot_wake(p_array, lam_K1, itr_ra, s1, st)

    def finalize(self, *args, **kwargs):
        """
        the method is called at the end of tracking

        :return:
        """
        # if self.pict_debug:
        #     data = np.array([ np.array(self.total_wake)])
        #     np.savetxt("total_wake_test.txt", data)

    def calcualte_csr_wakes(self):
        #B_params = [self.charge_binning_fraction, self.n_bin, self.m_bin, self.ip_method, self.sp, self.sigma_min]
        #s1, s2, Ns = self.bin_smoth.Q2EQUI_modified(B_params, SBINB)
        #st = (s2 - s1) / Ns
        #Ndw = [Ns - 1, st]
        gamma = 14 / m_e_GeV
        s_array = self.csr_traj[0, :]
        self.kernel_Ndw = [500, 25e-6]
        self.kernels = {}
        self.kern_indx = []
        for it in range(len(s_array))[1:]:
            kern = self.CSR_K1(it, self.csr_traj, self.kernel_Ndw, gamma=gamma)
            #print(len(kern))
            self.kernels[it] = kern
            self.kern_indx.append(it)

        return self.kernels

    def plot_wake(self, p_array, lam_K1, itr_ra, s1, st):
        """
        Method to plot CSR wakes on each step and save pictures in the working folder. Might be time-consuming.

        :param p_array:
        :param lam_K1:
        :param itr_ra:
        :param s1:
        :param st:
        :return:
        """
        self.napply += 1
        fig = self.plt.figure(figsize=(10, 8))

        ax1 = self.plt.subplot(311)
        ax1.plot(self.csr_traj[0, :], self.csr_traj[1, :]*1000, "r",  self.csr_traj[0, itr_ra], self.csr_traj[1, itr_ra]*1000, "bo")
        self.plt.ylabel("X [mm]")

        ax2 = self.plt.subplot(312)
        # # CSR wake in keV/m - preferable and not depend on step
        # ax2.plot(np.linspace(s1, s1+st*len(lam_K1), len(lam_K1))*1000, lam_K1/delta_s/1000.)
        # plt.ylabel("dE, keV/m")

        # CSR wake in keV - amplitude changes with step
        ax2.plot(np.linspace(s1, s1 + st * len(lam_K1), len(lam_K1)) * 1000, lam_K1 *1e-3)
        self.plt.setp(ax2.get_xticklabels(), visible=False)
        self.plt.ylim((-50, 50))
        self.plt.ylabel("Wake [keV]")

        # ax3 = plt.subplot(413)
        # n_points = len(lam_K1)
        # #wake = np.interp(np.linspace(s1, s1+st*n_points, n_points), np.linspace(s1, s1+st*len(lam_K1), len(lam_K1)), lam_K1)
        # self.total_wake += lam_K1
        # #plt.xlim(s1 * 1000, (s1 + st * len(lam_K1)) * 1000)
        # ax3.plot(np.linspace(s1, s1+st*n_points, n_points)*1000, self.total_wake*1e-6)
        # plt.ylabel("Total Wake [MeV]")
        # plt.setp(ax3.get_xticklabels(), visible=False)
        # plt.ylim((-10, 10))

        ax3 = self.plt.subplot(313, sharex=ax2)
        self.B = s_to_cur(p_array.tau(), sigma=np.std(p_array.tau())*0.05, q0=np.sum(p_array.q_array), v=speed_of_light)
        ax3.plot(-self.B[:, 0]*1000, self.B[:, 1], lw=2)
        ax3.set_ylabel("I [A]")
        ax3.set_xlabel("s [mm]")
        self.plt.subplots_adjust(hspace=0.2)
        dig = str(self.napply)
        name = "0" * (4 - len(dig)) + dig
        self.plt.savefig( name + '.png')

    def __repr__(self) -> str:
        cname = type(self).__name__
        sigma_min = self.sigma_min
        traj_step = self.traj_step
        apply_step = self.apply_step
        n_bin = self.n_bin
        step = self.step
        return f"<{cname}: {sigma_min=}, {traj_step=}, {apply_step=}, {n_bin=}, {step=}>"


class SaferCSR(CSR):
    r"""Run the CSR model but with checks on the sigma_min and the Derbenev criterion
    at every step.  All options applicable to CSR are also applicable to SaferCSR
    instances.

    The Derbenev criterion coefficient is given by

    .. math:: \kappa = \sigma_x \left( \frac{1}{R \sigma_z^2} \right) ^ \frac{1}{3}
                                \ll 1,

    where the symbols have their usual meanings.

    :param derbenev_criterion: The Derbenev criterion coefficient above which a
    warning should be raised for a ParticleArray instance.  By default 0.5.
    :param sigma_min_factor: The ratio of the bunch length relative to the
    sigma_min below which a warning should be raised.  By default 9.9.
    :apply_kick: Whether or not calculate and apply the CSR kick, or just do the
    checking without any physics process applied.  By default True.

    """

    def __init__(self,
                 derbenev_criterion=0.5,
                 sigma_min_factor=9.9,
                 apply_kick=True,
                 **kw
                 ):
        super().__init__(**kw)
        self.derbenev_criterion = derbenev_criterion
        # If bunch become shorter than sigma_min_factor*sigma_min, raise a warning.
        self.sigma_min_factor = sigma_min_factor
        self.apply_kick = apply_kick


    def query_bunch_length(self, bunch_length):
        minimum_permissable_bunch_length = self.sigma_min_factor * self.sigma_min
        if bunch_length < minimum_permissable_bunch_length:
            ratio = bunch_length / self.sigma_min
            bunch_length = round(bunch_length*1e6, 6)
            sigma_min = round(self.sigma_min*1e6, 6)
            logger.warning(
                f"Bunch is too small relative to sigma_min at s={self.z0}m: "
                f" bunch_length = {bunch_length}um;"
                f" sigma_min = {sigma_min}um;"
                f" ratio: {ratio};"
                f" minimum allowed ratio: {self.sigma_min_factor}"
            )

    @staticmethod
    def _derbenev_criterion(sigma_transverse, sigma_z, bending_radius):
        return sigma_transverse * (bending_radius * sigma_z**2) ** (-1/3)

    def query_derbenev_criterion(self, bunch_x, bunch_y, bunch_length):
        s = self.csr_traj[0]
        if len(s) < 2:
            return
        csr_traj_index = np.clip(np.searchsorted(s, self._lattice_to_trajectory_s(self.z0 - self.z_csr_start)), 1, len(s) - 1)

        # The points and tangents before and after the step takes place
        # pre_step = self.csr_traj[..., csr_traj_index-1][1:4]
        # post_step = self.csr_traj[..., csr_traj_index][1:4]
        tangent_pre_step = self.csr_traj[..., csr_traj_index-1][4:]
        tangent_post_step = self.csr_traj[..., csr_traj_index][4:]
        tangent_pre_step = tangent_pre_step / np.linalg.norm(tangent_pre_step)
        tangent_post_step = tangent_post_step / np.linalg.norm(tangent_post_step)

        angle = np.arctan2(np.linalg.norm(np.cross(tangent_pre_step, tangent_post_step)),
                           np.dot(tangent_post_step, tangent_pre_step))
        if angle == 0:
            return

        bending_radius = (s[csr_traj_index] - s[csr_traj_index-1]) / angle

        # Get index of coordinate that changes the most and basically just assume
        # we're bending in that plane.
        bending_plane_index = np.argmax(abs(tangent_post_step[:2]
                                            - tangent_pre_step[:2]))

        size = bunch_x
        if bending_plane_index == 1:
            size = bunch_y


        crit = self._derbenev_criterion(size, bunch_length, bending_radius)
        # Criterion should be must less than 1 for 1D to be valid.  As it's "much
        # less than", we say <0.5 by default...
        if crit > self.derbenev_criterion:
            s = self.z0
            logger.warning(
                f" Derbenev criterion may be violated at s={s}. "
                f" bunch length = {bunch_length*1e6}um;"
                f" tranverse size = {size*1e6}um;"
                f" bending radius = {bending_radius}m"
                f"  Derbenev coeff. = {crit}"
            )

    def apply(self, parray, ds):
        if parray.n == 0:
            return
        bunch_length = parray.tau().std()
        bunch_x = parray.x().std()
        bunch_y = parray.y().std()

        self.query_bunch_length(bunch_length)
        self.query_derbenev_criterion(bunch_x, bunch_y, bunch_length)

        if self.apply_kick:
            return super().apply(parray, ds)
