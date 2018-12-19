"""
@ authors Martin Dohlus DESY, 2015, Sergey Tomin XFEL, 2016
"""

from scipy import interpolate
from ocelot.common.globals import *
from ocelot.common import math_op
from scipy.ndimage.filters import gaussian_filter
from scipy.optimize import curve_fit
from ocelot.cpbd.beam import *
from ocelot.cpbd.high_order import *
from ocelot.cpbd.magnetic_lattice import *
import time
from scipy.integrate import cumtrapz
from ocelot.cpbd.physics_proc import PhysProc
import copy
import logging
logger = logging.getLogger(__name__)

try:
    import numba as nb
    nb_flag = True
except:
    logger.info("csr.py: module NUMBA is not installed. Install it to speed up calculation")
    nb_flag = False

try:
    from pyfftw.interfaces.numpy_fft import fft
    from pyfftw.interfaces.numpy_fft import ifft
except:
    logger.info("csr.py: module PYFFTW is not installed. Install it to speed up calculation.")
    from numpy.fft import ifft
    from numpy.fft import fft

try:
    import numexpr as ne
    ne_flag = True
except:
    logger.info("csr.py: module NUMEXPR is not installed. Install it to speed up calculation")
    ne_flag = False


def nextpow2(p):
    i = 0
    pow2 = 2**i
    while p > pow2:
        i += 1
        pow2 = 2**i
    return i


def csr_convolution(a, b):
    P = len(a)
    Q = len(b)
    L = P + Q - 1
    K = 2 ** nextpow2(L)
    a_pad = np.pad(a, (0, K - P), 'constant', constant_values=(0))
    b_pad = np.pad(b, (0, K - Q), 'constant', constant_values=(0))
    c = ifft(fft(a_pad)*fft(b_pad))
    c = c[0:L-1].real
    return c


def interp1(x, y, xnew, k=1):
    if len(xnew) > 0:
        if k == 1:
            ynew = np.interp(xnew, x, y)
        else:
            tck = interpolate.splrep(x, y, k=k)
            ynew = interpolate.splev(xnew, tck, der=0)
    else:
        ynew = []
    return ynew


def convolution(xu, u, xw, w):
    """
    convolution of equally spaced functions
    """
    hx = xu[1] - xu[0]
    wc = np.convolve(u, w)*hx
    nw = w.shape[0]
    nu = u.shape[0]
    x0 = xu[0] + xw[0]
    xc = x0 + np.arange(nw + nu)*hx
    return xc, wc


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


class Smoothing:
    def __init__(self):
        self.print_log = False
        if nb_flag:
            logger.debug("Smoothing: NUMBA")
            self.q_per_step_ip2 = nb.jit()(self.q_per_step_ip2_py)
        else:
            logger.debug("Smoothing: Python")
            self.q_per_step_ip2 = self.q_per_step_ip2_py

    def q_per_step_ip2_py(self, N_BIN, Q_BIN, BIN0, BIN1, NSIG, RMS, step, Nz, z1):
        Nz = int(Nz)
        charge_per_step = np.zeros(Nz)

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
            fact = step / (np.sqrt(2 * pi) * sigma)
            for k in range(k1, k2):
                xx = (k - 1) * step
                yy = fact * np.exp(-0.5 * ((xx - mitte) / sigma) ** 2)
                charge_per_step[int(k - 1)] += yy * qbin
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
        SP = ? parameter for gauss
        sigma_min = minimal sigma, if IP_method == 2
        step_unit = if positive --> step=integer * step_unit
        output
        z1, z2, Nz = equidistant mesh(Nz meshlines)
        charge_per_step = charge per step, charge_per_step(1:Nz)
        bins might overlapp!
        """

        N_BIN = BS_params[1]
        M_BIN = BS_params[2]
        K_BIN = N_BIN * M_BIN # number of sub - bins
        I_BIN = K_BIN - (M_BIN - 1) # number of bin intervalls
        # put charges to sub - bins
        Q_BIN = np.zeros(K_BIN)
        n2 = 0
        if np.size(q) == 1:
            for k in range(K_BIN):
                if NBIN[k] > 0:
                    n1 = n2 + 1
                    n2 = n2 + NBIN[k]
                    Q_BIN[k] = q * (n2 - n1 + 1)
        else:
            for k in range(K_BIN):
                if NBIN[k] > 0:
                    n1 = n2 + 1
                    n2 = n2 + NBIN[k]
                    Q_BIN[k] = np.sum(q[int(n1-1):int(n2)])

        # put sub - bins to bins
        qsum = np.append([0], np.cumsum(Q_BIN))
        BIN = [SBINB[0:I_BIN], SBINB[M_BIN + np.arange(I_BIN)]]
        Q_BIN = (qsum[M_BIN + np.arange(I_BIN)] - qsum[0:I_BIN])/M_BIN

        # ......................................................................

        # interpolation parameters
        IP_method = BS_params[3]
        SP = BS_params[4]
        if SP <= 0:
            SP = 0.5

        sigma_min = max(0, BS_params[5])
        if len(BS_params) < 7:
            step_unit = 0
        else:
            step_unit = max(0, BS_params[6])


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
            RMS = SP * (BIN[1] - BIN[0])
            for nb in range(N_BIN):
                RMS[nb] = max(RMS[nb], sigma_min)
            z1 = np.min(MITTE - NSIG * RMS)
            z2 = np.max(MITTE + NSIG * RMS)
            step = 0.25 * min(RMS)
        else:
            z1 = np.min(BIN[0][:])
            z2 = np.max(BIN[1][:])
            step = 0.5 * min(BIN[1][:] - BIN[0][:])

        if step_unit > 0:
            step = step_unit * max(1, np.round(step / step_unit))
            z1 = step * np.floor(z1 / step)
            z2 = step * np.ceil(z2 / step)
            Nz = np.round((z2 - z1) / step)
        else:
            Nz = np.round((z2 - z1) / step)
            step = (z2 - z1) / Nz
        charge_per_step = np.zeros(int(Nz))
        if IP_method == 1:
            for nb in range(N_BIN):
                aa = BIN[0][nb] - z1
                bb = BIN[1][nb] - z1
                qps = step * Q_BIN(nb) / (bb - aa)
                a = ((3. * aa - bb) / 2.) / step + 1.
                k1 = min(Nz, max(1, np.floor(a)))
                b = ((aa + bb) / 2.) / step + 1.
                c = ((3. * bb - aa) / 2.) / step + 1.
                k2 = min(Nz, max(1, np.ceil(c)))
                for k in range(k1, k2):
                    w = sample_1(k, a, b, c)
                    charge_per_step[k-1] += w * qps

        elif IP_method == 2:
            charge_per_step = self.q_per_step_ip2(N_BIN, Q_BIN, BIN[0], BIN[1], NSIG, RMS, step, Nz, z1)
        else:
            for nb in range(N_BIN):
                aa = BIN[0][nb] - z1
                bb = BIN[1][nb] - z1
                qps = step * Q_BIN[nb] / (bb - aa)
                a = aa / step + 1
                k1 = min(Nz, max(1, np.floor(a)))
                b = bb / step + 1
                k2 = min(Nz, max(1, np.ceil(b)))
                for k in range(k1, k2):
                    w = sample_0(k, a, b)
                    charge_per_step[k-1] += w * qps
        return z1, z2, Nz, charge_per_step


class SubBinning:
    def __init__(self, x_qbin, n_bin, m_bin):
        self.x_qbin = x_qbin
        self.n_bin = n_bin
        self.m_bin = m_bin
        self.print_log = False
        if nb_flag:
            logger.debug("SubBinning: NUMBA")
            self.p_per_subbins = nb.jit(nb.double[:](nb.double[:], nb.double[:], nb.int64))(self.p_per_subbins_py)
        else:
            logger.debug("SubBinning: Python")
            self.p_per_subbins = self.p_per_subbins_py

    def p_per_subbins_py(self, s, SBINB, K_BIN):
        NBIN = np.zeros(K_BIN)
        ib = 0
        Ns = len(s)
        for n in range(Ns):
            while s[n] >= SBINB[ib + 1] and ib < K_BIN - 1:
                ib = ib + 1
            NBIN[ib] = NBIN[ib] + 1
        return NBIN


    def subbin_bound(self, q, s, x_qbin, n_bin, m_bin):
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
        X_QBIN = x_qbin
        N_BIN = n_bin
        M_BIN = m_bin
        K_BIN = N_BIN * M_BIN  # number of sub-bins

        Ns = len(s)
        # binning intervalls
        # aa=monoton "charge" vector

        if np.size(q) != 1:
            q_cumsum = np.cumsum(q)
            q_cumsum_r = q_cumsum - q_cumsum[0]
            aa = 0.5 * q_cumsum + 0.5 * np.append([0], q_cumsum_r[1:])
            # aa = 0.5 * np.cumsum(q) + 0.5 * np.append([0], np.cumsum(q[1:]))
        else:
            # % aa(1)=q/2; for n=2:Ns, aa(n)=aa(n-1)+q; end
            aa = q * (np.arange(1, Ns + 1) - 0.5)

        if ne_flag:
            aa0 = aa[0]
            aaNs = aa[Ns - 1]
            s0 = s[0]
            sNs = s[Ns - 1]
            # aa = ne.evaluate('(aa - aa0) / (aaNs - aa0)')
            # bb=monoton "length" vector
            # bb = ne.evaluate('(s - s0) / (sNs - s0)')
            # aa=LK of "charge" and "length" vector; avoid zero stepwidth
            aa = ne.evaluate('(aa - aa0) / (aaNs - aa0) * X_QBIN + (s - s0) / (sNs - s0) * (1 - X_QBIN)')
        else:
            aa = (aa - aa[0]) / (aa[Ns - 1] - aa[0])
            # bb=monoton "length" vector
            bb = (s - s[0]) / (s[Ns - 1] - s[0])
            # aa=LK of "charge" and "length" vector; avoid zero stepwidth
            aa = aa * X_QBIN + bb * (1 - X_QBIN)
        if np.min(np.diff(aa)) == 0:
            aa = 0.999 * aa + 0.001 * (np.arange(0, Ns)) / (Ns - 1)

        # vector with bin boundaries
        SBINB = interp1(aa, s, np.arange(K_BIN + 1.) / K_BIN)
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
            self.K0_1 = nb.jit()(self.K0_1_jit)
            self.K0_0 = nb.jit()(self.K0_0_jit)
            self.eval = self.K0_fin_anf_opt
        elif ne_flag:
            logger.debug("K0_fin_anf: NumExpr")
            self.eval = self.K0_fin_anf_numexpr
        else:
            logger.debug("K0_fin_anf: Python")
            self.eval = self.K0_fin_anf_np

    def K0_1_jit(self, indx, j, R, n, traj4, traj5, traj6, w, gamma):
        g2i = 1. / gamma ** 2
        b2 = 1. - g2i
        beta = np.sqrt(b2)
        K = np.zeros(indx - j)
        for i in range(j, indx):
            Ri = R[i]
            n0i = n[i, 0] / Ri
            n1i = n[i, 1] / Ri
            n2i = n[i, 2] / Ri
            # kernel
            t4 = traj4[i]
            t5 = traj5[i]
            t6 = traj6[i]
            x = n0i * t4 + n1i * t5 + n2i * t6
            K[i - j] = ((beta * (x - n0i * traj4[indx] - n1i * traj5[indx] - n2i * traj6[indx]) -
                         b2 * (1. - t4 * traj4[indx] - t5 * traj5[indx] - t6 * traj6[indx]) - g2i) / Ri - (
                        1. - beta * x) / w[i - j] * g2i)
        return K

    def K0_0_jit(self, i, traj0, traj1, traj2, traj3, gamma, s, n, R, w):
        g2i = 1. / gamma ** 2
        b2 = 1. - g2i
        beta = np.sqrt(b2)
        # i1 = i - 1  # ignore points i1+1:i on linear path to observer

        traj0i = traj0[i]
        traj1i = traj1[i]
        traj2i = traj2[i]
        traj3i = traj3[i]
        for j in range(i):
            s[j] = traj0[j] - traj0i
            n1 = traj1i - traj1[j]
            n2 = traj2i - traj2[j]
            n3 = traj3i - traj3[j]
            R[j] = np.sqrt(n1 * n1 + n2 * n2 + n3 * n3)
            w[j] = s[j] + beta * R[j]
            n[j, 0] = n1
            n[j, 1] = n2
            n[j, 2] = n3

    def K0_fin_anf_opt(self, i, traj, wmin, gamma):
        s = np.zeros(i)
        n = np.zeros((i, 3))
        R = np.zeros(i)
        w = np.zeros(i)
        self.K0_0(i, traj[0], traj[1], traj[2], traj[3], gamma, s, n, R, w)
        j = np.where(w <= wmin)[0]

        if len(j) > 0:
            j = j[-1]
            w = w[j:i]
            s = s[j:i]
        else:
            j = 0
        K = self.K0_1(i, j, R, n, traj[4], traj[5], traj[6], w, gamma)

        if len(K) > 1:
            a = np.append(0.5 * (K[0:-1] + K[1:]) * np.diff(s), 0.5 * K[-1] * s[-1])
            KS = np.cumsum(a[::-1])[::-1]
            # KS = cumsum_inv_jit(a)
            # KS = cumtrapz(K[::-1], -s[::-1], initial=0)[::-1] + 0.5*K[-1]*s[-1]
        else:
            KS = 0.5 * K[-1] * s[-1]
        return w, KS

    def K0_fin_anf_np(self, i, traj, wmin, gamma):
        # function [ w,KS ] = K0_inf_anf( i,traj,wmin,gamma )

        g2i = 1. / gamma ** 2
        b2 = 1. - g2i
        beta = np.sqrt(b2)
        i1 = i - 1  # ignore points i1+1:i on linear path to observer
        ind1 = i1 + 1
        s = traj[0, 0:ind1] - traj[0, i]
        n = np.array([traj[1, i] - traj[1, 0:ind1],
                      traj[2, i] - traj[2, 0:ind1],
                      traj[3, i] - traj[3, 0:ind1]])
        R = np.sqrt(np.sum(n ** 2, axis=0))

        w = s + beta * R
        j = np.where(w <= wmin)[0]

        if len(j) > 0:
            j = j[-1]
            w = w[j:ind1]
            s = s[j:ind1]
        else:
            j = 0
        R = R[j:ind1]
        n0 = n[0, j:ind1] / R
        n1 = n[1, j:ind1] / R
        n2 = n[2, j:ind1] / R

        # kernel
        t4 = traj[4, j:i1 + 1]
        t5 = traj[5, j:i1 + 1]
        t6 = traj[6, j:i1 + 1]

        x = n0 * t4 + n1 * t5 + n2 * t6
        K = ((beta * (x - n0 * traj[4, i] - n1 * traj[5, i] - n2 * traj[6, i]) -
              b2 * (1. - t4 * traj[4, i] - t5 * traj[5, i] - t6 * traj[6, i]) - g2i) / R - (1. - beta * x) / w * g2i)

        # K = ((beta*(n0*(t4 - traj[4, i]) +
        #            n1*(t5 - traj[5, i]) +
        #            n2*(t6 - traj[6, i])) -
        #    b2*(1. - t4*traj[4, i] - t5*traj[5, i] - t6*traj[6, i]) - g2i)/R[ra] -
        #    (1. - beta*(n0*t4 + n1*t5 + n2*t6))/w*g2i)

        # integrated kernel: KS=int_s^0{K(u)*du}=int_0^{-s}{K(-u)*du}

        if len(K) > 1:
            a = np.append(0.5 * (K[0:-1] + K[1:]) * np.diff(s), 0.5 * K[-1] * s[-1])
            KS = np.cumsum(a[::-1])[::-1]
            # KS = cumtrapz(K[::-1], -s[::-1], initial=0)[::-1] + 0.5*K[-1]*s[-1]
        else:
            KS = 0.5 * K[-1] * s[-1]

        return w, KS

    def K0_fin_anf_numexpr(self, i, traj, wmin, gamma):
        # function [ w,KS ] = K0_inf_anf( i,traj,wmin,gamma )

        g2i = 1. / gamma ** 2
        b2 = 1. - g2i
        beta = np.sqrt(b2)
        i1 = i - 1  # ignore points i1+1:i on linear path to observer
        # ra = np.arange(0, i1+1)
        ind1 = i1 + 1
        s = traj[0, 0:ind1] - traj[0, i]
        n0 = traj[1, i] - traj[1, 0:ind1]
        n1 = traj[2, i] - traj[2, 0:ind1]
        n2 = traj[3, i] - traj[3, 0:ind1]
        R = ne.evaluate("sqrt(n0**2 + n1**2 + n2**2)")

        w = ne.evaluate('s + beta*R')
        j = np.where(w <= wmin)[0]

        if len(j) > 0:
            j = j[-1]
            w = w[j:ind1]
            s = s[j:ind1]
        else:
            j = 0
        R = R[j:ind1]
        n0 = n0[j:ind1] / R
        n1 = n1[j:ind1] / R
        n2 = n2[j:ind1] / R

        # kernel
        t4 = traj[4, j:i1 + 1]
        t5 = traj[5, j:i1 + 1]
        t6 = traj[6, j:i1 + 1]

        x = ne.evaluate('n0*t4 + n1*t5 + n2*t6')

        t4i = traj[4, i]
        t5i = traj[5, i]
        t6i = traj[6, i]
        K = ne.evaluate(
            '((beta*(x - n0*t4i- n1*t5i - n2*t6i) - b2*(1. - t4*t4i - t5*t5i - t6*t6i) - g2i)/R - (1. - beta*x)/w*g2i)')

        if len(K) > 1:
            a = np.append(0.5 * (K[0:-1] + K[1:]) * np.diff(s), 0.5 * K[-1] * s[-1])
            KS = np.cumsum(a[::-1])[::-1]
            # KS = cumtrapz(K[::-1], -s[::-1], initial=0)[::-1] + 0.5*K[-1]*s[-1]
        else:
            KS = 0.5 * K[-1] * s[-1]

        return w, KS


class CSR(PhysProc):
    """
    coherent synchrotron radiation
    Attributes:
        self.step = 1 [in Navigator.unit_step] - step of the CSR kick applying for beam (ParticleArray)
        self.sigma_min = 1.e-4  - minimal sigma if gauss filtering applied
        self.traj_step = 0.0002 [m] - trajectory step or, other words, integration step for calculation of the CSR-wake
        self.apply_step = 0.0005 [m] - step of the calculation CSR kick, to calculate average CSR kick
    """
    def __init__(self):
        PhysProc.__init__(self)
        # binning parameters
        self.x_qbin = 0             # length or charge binning; 0... 1 = length...charge
        self.n_bin = 100            # number of bins
        self.m_bin = 5              # multiple binning(with shifted bins)

        # smoothing
        self.ip_method = 2          # = 0 / 1 / 2 for rectangular / triangular / gauss
        self.sp = 0.5               # ? parameter for gauss
        self.sigma_min = 1.e-4      # minimal sigma, if ip_method == 2
        self.step_unit = 0          # if positive --> step=integer * step_unit

        # trajectory
        self.traj_step = 0.0002     # [m] step of the trajectory
        self.energy = None          # [GeV], if None, beta = 1 and calculation of the trajectory with RK is not possible

        # CSR kick
        self.apply_step = 0.0005    # [m] step of the calculation CSR kick: csr_kick += csr(apply_step)
        self.step = 1               # step in the unit steps, step_in_[m] = self.step * navigator.unit_step [m].
                                    # The CSR kick is applied at the end of the each step

        self.z_csr_start = 0.       # z [m] position of the start_elem
        self.z0 = 0.                # self.z0 = navigator.z0 in track.track()

        self.debug = False
        # another filter
        self.filter_order = 10
        self.n_mesh = 345
        self.pict_debug = True
        #self.print_log = False

        self.sub_bin = SubBinning(x_qbin=self.x_qbin, n_bin=self.n_bin, m_bin=self.m_bin)
        self.bin_smoth = Smoothing()
        self.k0_fin_anf = K0_fin_anf()
        self.rk_traj = False
        #self.napply = 0
        #if self.pict_debug:
        #    self.napply = 0
        #    #self.f = plt.figure(figsize=(12, 9))
        #    #plt.ion()
        #    #plt.hold(False)

    def K0_inf_anf(self, i, traj, wmin):
        # function [ w,KS ] = K0_inf_anf( i,traj,wmin )

        i1 = i-1 # ignore points i1+1:i on linear path to observer
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
           (1. - traj[4, ra]*traj[4, i] - traj[5, ra]*traj[5, i] - traj[6, ra]*traj[6, i]))/R[ra]

        # integrated kernel: KS=int_s^0{K(u)*du}=int_0^{-s}{K(-u)*du}
        if np.shape(K)[0] > 1:
            a = np.append(0.5*(K[0:-1] + K[1:])*np.diff(s), 0.5*K[-1]*s[-1])
            KS = np.cumsum(a[::-1])[::-1]
            #KS = np.fliplr(np.cumsum(np.fliplr([0.5*(K[1:-1] + K[2:])*np.diff(s), 0.5*K[-1]*s[-1]])))
        else:
            KS = 0.5*K[-1]*s[-1]

        return w, KS

    def K0_fin_inf(self, i, traj, w_range, gamma):
        # function [ KS ] = K0_inf_inf( i,traj,w_range,gamma )

        g2 = gamma**2
        g2i = 1./g2
        b2 = 1. - g2i
        beta = np.sqrt(b2)
        # winf
        Rv1 = traj[1:4, i] - traj[1:4, 0]
        s1 =  traj[0, 0] - traj[0, i]
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

        if a2/R[1]**2 > 1e-7:
            KS = (beta*(1. - np.dot(ev1, evo))*np.log(R[0]/R) - beta*np.dot(uup, evo)*(np.arctan((s[0] - winf)/a) - np.arctan((s-winf)/a))
               - (b2*np.dot(ev1, evo) - 1)*np.log((winf - s + R)/(winf-s[0] + R[0]))
               + g2i*np.log(w_range[0]/w_range))

        else:
            KS = (beta*(1. - np.dot(ev1, evo))*np.log(R[0]/R)
               - (b2*np.dot(ev1, evo) - 1.)*np.log((winf - s + R)/(winf - s[0] + R[0]))
               + g2i*np.log(w_range[0]/w_range))
        return KS

    def K0_inf_inf(self, i, traj, w_range):
        # winf
        Rv1 = traj[1:4, i] - traj[1:4, 0]
        s1 =  traj[0, 0] - traj[0, i]
        ev1 = traj[4:, 0]
        evo = traj[4:, i]
        winfms1 = np.dot(Rv1, ev1)
        aup = -Rv1 + winfms1*ev1
        a2 = np.dot(aup, aup)
        a = np.sqrt(a2)
        uup = aup/a
        winf = s1 + winfms1

        Nvalid = np.where(winf<w_range)[0]
        if len(Nvalid) >0:
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
        :param NdW: list N[0] 0 number of mesh points, N[1] = dW> 0 - increment, Mesh = Mesh = (N: 0) * dW
        :param gamma:
        :return:
        """
        # function [ K1 ] = CSR_K1( i,traj,NdW,gamma )

        # L_fin=nargin==4 && ~isempty(gamma) && gamma>1
        if gamma != None:
            L_fin = True
        else:
            L_fin = False

        w_range = np.arange(-NdW[0]-1, 0)*NdW[1]

        if L_fin:
            w, KS = self.k0_fin_anf.eval(i, traj, w_range[0], gamma)
        else:
            w, KS = self.K0_inf_anf(i, traj, w_range[0])

        KS1 = KS[0]

        w, idx = np.unique(w, return_index=True)
        KS = KS[idx]

        # sort and unique takes time, but is required to avoid numerical trouble
        if w_range[0] < w[0]:
            m = np.where(w_range < w[0])[0][-1]
            if L_fin:
                KS2 = self.K0_fin_inf(i, traj, np.append(w_range[0:m+1], w[0]), gamma)
            else:
                KS2 = self.K0_inf_inf(i, traj, np.append(w_range[0:m+1], w[0]))

            KS2 = (KS2[-1] - KS2) + KS1
            KS = np.append(KS2[0:-1], interp1(w, KS, w_range[m+1:]))

        else:
            KS = interp1(w, KS, w_range)
        four_pi_eps0 = 1./(1e-7*speed_of_light**2)
        K1 = np.diff(np.append(np.diff(np.append(KS, 0)), 0))/NdW[1]/four_pi_eps0

        return K1

    def prepare(self, lat):
        """
        calculation of trajectory in rectangular coordinates
        calculation of the z_csr_start
        :param lat: Magnetic Lattice
        :return: self.csr_traj: trajectory. traj[0,:] - longitudinal coordinate,
                                 traj[1,:], traj[2,:], traj[3,:] - rectangular coordinates, \
                                 traj[4,:], traj[5,:], traj[6,:] - tangential unit vectors
        """
        self.z_csr_start = sum([p.l for p in lat.sequence[:self.indx0]])
        p = Particle()
        beta = 1. if self.energy == None else np.sqrt(1. - 1./(self.energy/m_e_GeV)**2)
        self.csr_traj = np.transpose([[0, p.x, p.y, p.s, p.px, p.py, 1.]])
        if Undulator in [elem.__class__ for elem in lat.sequence[self.indx0:self.indx1+1]]:
            self.rk_traj = True
        for elem in lat.sequence[self.indx0:self.indx1+1]:

            if elem.l == 0:
                continue
            delta_s = elem.l
            step = self.traj_step
            if elem.__class__ in [Bend, RBend, SBend] and not self.rk_traj:
                if elem.angle != 0:
                    R = -elem.l/elem.angle
                else:
                    R = 0
                Rx = R * np.cos(elem.tilt)
                Ry = R * np.sin(elem.tilt)
                #B = energy*1e9*beta/(R*speed_of_light)
                R_vect = [-Ry, Rx, 0]
                self.csr_traj = arcline(self.csr_traj, delta_s, step, R_vect)

            elif elem.__class__ in [Bend, RBend, SBend, XYQuadrupole, Undulator] and self.energy is not None and self.rk_traj:
                """
                rk_track_in_field accepts initial conditions (initial coordinates) in the ParticleArray.rparticles format
                from another hand the csr module use trajectory in another format (see csr.py)
                """

                if elem.l < 1e-10:
                    continue
                delta_z = delta_s
                if elem.__class__ is XYQuadrupole:
                    hx = elem.k1 * elem.x_offs
                    hy = -elem.k1 * elem.y_offs
                    By = self.energy * 1e9 * beta * hx / speed_of_light
                    Bx = -self.energy * 1e9 * beta * hy / speed_of_light
                    mag_field = lambda x, y, z: (Bx, By, 0)

                elif elem.__class__ == Undulator:
                    gamma = self.energy/m_e_GeV
                    ku = 2 * np.pi / elem.lperiod
                    delta_z = elem.lperiod * elem.nperiods

                    delta_s += delta_z / (gamma ) ** 2 * (1 + 0.5 * (elem.Kx ) ** 2)

                    By = elem.Kx * m_e_eV * 2. * pi / (elem.lperiod * speed_of_light)
                    Bx = elem.Ky * m_e_eV * 2. * pi / (elem.lperiod * speed_of_light)
                    mag_field = lambda x, y, z: (0, By * np.cos(ku * z), 0)

                else:
                    delta_z = delta_s * np.sin(elem.angle) / elem.angle if elem.angle != 0 else delta_s
                    hx = elem.angle / elem.l * np.cos(elem.tilt)
                    hy = elem.angle / elem.l * np.sin(elem.tilt)
                    By = self.energy * 1e9 * beta * hx / speed_of_light
                    Bx = -self.energy * 1e9 * beta * hy / speed_of_light
                    mag_field = lambda x, y, z: (Bx, By , 0)

                sre0 = self.csr_traj[:, -1]
                N = int(max(1, np.round(delta_s / step)))
                SRE2 = np.zeros((7, N))

                rparticle0 = np.array([[self.csr_traj[1, -1]], [self.csr_traj[4, -1] / self.csr_traj[6, -1]],
                                       [self.csr_traj[2, -1]], [self.csr_traj[5, -1] / self.csr_traj[6, -1]], [0], [0]])
                traj = rk_track_in_field(rparticle0, s_stop=delta_z, N=N + 1, energy=self.energy, mag_field=mag_field,
                                         s_start=0)

                x = traj[0::9].flatten()
                y = traj[2::9].flatten()
                xp = traj[1::9].flatten()
                yp = traj[3::9].flatten()
                z = traj[4::9].flatten()
                print(xp[-1])
                xp2 = xp * xp
                yp2 = yp * yp
                zp = np.sqrt(1./(1. + xp2 + yp2))

                dS = float(delta_s) / N

                s = cumtrapz(np.sqrt(1. + xp2 + yp2), z, initial=0)

                s_unif = np.arange(0, N + 1) * dS

                z_s = np.interp(s_unif, s, z)

                x = np.interp(z_s, z, x)
                y = np.interp(z_s, z, y)
                zp_s = np.interp(z_s, z, zp)

                xp = np.interp(z_s, z, xp*zp)
                yp = np.interp(z_s, z, yp*zp)

                SRE2[0, :] = sre0[0] + s_unif[1:]
                SRE2[1, :] = x[1:]
                SRE2[2, :] = y[1:]
                SRE2[3, :] = sre0[3] + z_s[1:]
                SRE2[4, :] = xp[1:]
                SRE2[5, :] = yp[1:]
                SRE2[6, :] = zp_s[1:]
                self.csr_traj = np.append(self.csr_traj, SRE2, axis=1)
            else:
                R_vect = [0, 0, 0.]
                self.csr_traj = arcline(self.csr_traj, delta_s, step, R_vect )

        # import matplotlib.pyplot as plt
        # plt.figure(10)
        # plt.plot(self.csr_traj[0, :], self.csr_traj[1+3, :], "r")
        # plt.plot(self.csr_traj[0, :], self.csr_traj[2+3, :], "b")
        # plt.legend(["X", "Y"])
        # #plt.plot( self.csr_traj[3, :] - csr_traj_arc[3, :], "r.")
        # #plt.plot(csr_traj_arc[0, :], csr_traj_arc[4, :]  - self.csr_traj[4, :], "b.")
        # plt.show()
        return self.csr_traj

    def apply(self, p_array, delta_s):
        if delta_s < self.traj_step:
            logger.debug("CSR delta_s < self.traj_step")
            return
        s_cur = self.z0 - self.z_csr_start
        z = -p_array.tau()
        ind_z_sort = np.argsort(z)
        #SBINB, NBIN = subbin_bound(p_array.q_array, z[ind_z_sort], self.x_qbin, self.n_bin, self.m_bin)
        #B_params = [self.x_qbin, self.n_bin, self.m_bin, self.ip_method, self.sp, self.sigma_min]
        #s1, s2, Ns, lam_ds = Q2EQUI(p_array.q_array[ind_z_sort], B_params, SBINB, NBIN)
        SBINB, NBIN = self.sub_bin.subbin_bound(p_array.q_array, z[ind_z_sort], self.x_qbin, self.n_bin, self.m_bin)
        B_params = [self.x_qbin, self.n_bin, self.m_bin, self.ip_method, self.sp, self.sigma_min]
        s1, s2, Ns, lam_ds = self.bin_smoth.Q2EQUI(p_array.q_array[ind_z_sort], B_params, SBINB, NBIN)
        st = (s2 - s1) / Ns
        sa = s1 + st / 2.
        Ndw = [Ns - 1, st]

        s_array = self.csr_traj[0, :]
        indx = (np.abs(s_array-s_cur)).argmin()
        indx_prev = (np.abs(s_array - (s_cur - delta_s))).argmin()
        gamma = p_array.E/m_e_GeV
        h = max(1., self.apply_step/self.traj_step)

        itr_ra = np.unique(-np.round(np.arange(-indx, -indx_prev, h))).astype(np.int)

        nit = 0
        n_iter = len(itr_ra)
        #start = time.time()
        K1 = self.CSR_K1(itr_ra[nit], self.csr_traj, Ndw, gamma)
        for nit in range(1, n_iter):
            K1 += self.CSR_K1(itr_ra[nit], self.csr_traj, Ndw, gamma=gamma)
        K1 = K1/n_iter


        lam_K1 = csr_convolution(lam_ds, K1[::-1]) / st * delta_s

        z_sort = z[ind_z_sort]
        dE = np.interp(z_sort*(1./st)+(0.-sa/st), np.arange(len(lam_K1)), lam_K1)

        pc_ref = np.sqrt(p_array.E ** 2 / m_e_GeV ** 2 - 1) * m_e_GeV
        delta_p = dE * 1e-9 / pc_ref
        p_array.rparticles[5][ind_z_sort] += delta_p

        # self.napply += 1
        # if self.pict_debug:# and self.napply%2 == 0:
        # from matplotlib import pyplot as plt
        # #from ocelot.gui.accelerator import show_density
        # fig = plt.figure(figsize=(10, 7))
        # ax1 = plt.subplot(311)
        # ax1.plot(self.csr_traj[3, :], self.csr_traj[1, :], "r",  self.csr_traj[3,itr_ra], self.csr_traj[1, itr_ra], "bo")
        # ax2 = plt.subplot(312)
        # # CSR wake
        # # plt.xlim(s1 * 1000, (s1 + st * len(lam_K1)) * 1000)
        # ax2.plot(np.linspace(s1, s1+st*len(lam_K1), len(lam_K1))*1000, lam_K1/delta_s/1000.)
        # # #plt.ylim(-2500, 2500)
        # # plt.ylabel("dE, keV/m")
        # # Energy profil
        # #plt.xlim(s1 * 1000, (s1 + st * len(lam_K1)) * 1000)
        # #show_density(-p_array.tau() * 1e3, p_array.p() * 1e2, ax=ax2, nbins_x=400, nbins_y=400,
        # #             interpolation="bilinear", ylabel='$\delta_E$ [%]',
        # #             title="Longitudinal phase space", grid=False, show_xtick_label=False
        # #self.f.add_subplot(313)
        # ax3 = plt.subplot(313, sharex=ax2)
        # ax3.set_title("Beam current")
        # B = s_to_cur(p_array.tau(), sigma=np.std(p_array.tau())*0.05, q0=np.sum(p_array.q_array), v=speed_of_light)
        # ax3.plot(-B[:, 0]*1000, B[:, 1], lw=2)
        # ax3.set_ylabel("I, [A]")
        # ax3.set_xlabel("s, [mm]")
        # ax3.set_ylim(0, 8000)
        # plt.subplots_adjust(hspace=0.3)

        # dig = str(self.napply)
        # name = "0" * (4 - len(dig)) + dig
        # plt.show()
        # plt.savefig( name + '.png')
        # plt.draw()
        # plt.pause(0.01)
        #
        # if self.debug:
        #    fig, ax1 = plt.subplots()
        #    #ax1.plot(-I[::-1, 0]*(1/st)+(1. - sa/st), I[:, 1], "r")
        #    ax1.plot(lam_ds, "r")
        #    #ax1.plot(K1, "r")
        #    ax1.set_ylabel("I, A", color='r')
        #    ax1.set_xlabel('s, mm')
        #    ax2 = ax1.twinx()
        #    ax2.plot(lam_K1, "b")
        #    ax2.set_ylabel("kernel", color='b')
        #    plt.show()


    def apply_i(self, p_array, delta_s):

        s_cur = self.z0 - self.z_csr_start
        z = p_array.tau()
        s1 = min(z)
        s2 = max(z)
        bunch_size = s2 - s1
        st = bunch_size / (self.n_mesh + 1)
        I = s2current(z, p_array.q_array, n_points=self.n_mesh, filter_order=self.filter_order, mean_vel=speed_of_light)
        Ns = len(I[:, 0])
        sa = s1 + st / 2.
        Ndw = [self.n_mesh, st]
        lam_ds = I[:, 1] / speed_of_light * st

        s_array = self.csr_traj[0, :]
        indx = (np.abs(s_array - s_cur)).argmin()
        indx_prev = (np.abs(s_array - (s_cur - delta_s))).argmin()
        gamma = p_array.E / m_e_GeV
        h = max(1., self.apply_step / self.traj_step)
        itr_ra = np.unique(-np.round(np.arange(-indx, -indx_prev, h))).astype(np.int)

        nit = 0
        n_iter = len(itr_ra)
        K1 = self.CSR_K1(itr_ra[nit], self.csr_traj, Ndw, gamma)
        for nit in range(1, n_iter):
            K1 += self.CSR_K1(itr_ra[nit], self.csr_traj, Ndw, gamma=gamma)
        K1 = K1 / n_iter


        lam_K1 = csr_convolution(lam_ds, K1) / st * delta_s
        Nend = len(lam_K1)
        N = int(abs(self.n_mesh + 1 - Nend) + self.n_mesh + 1)
        x = np.linspace(I[-1, 0], I[-1, 0] - N * Ndw[1], num=N, endpoint=False)[::-1]
        tck = interpolate.splrep(x, lam_K1, k=1)
        dE = interpolate.splev(z, tck, der=0)
        pc_ref = np.sqrt(p_array.E ** 2 / m_e_GeV ** 2 - 1) * m_e_GeV
        delta_p = dE * 1e-9 / pc_ref
        p_array.rparticles[5] += delta_p














