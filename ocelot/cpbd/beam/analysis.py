import numpy as np
import ocelot.common.globals as glb
from numba import njit

from . import beam_utils
from . import core
from ocelot.common.ocelog import *

_logger = logging.getLogger(__name__)


try:
    import numexpr as ne

    ne_flag = True
except:
    _logger.debug("beam.py: module NUMEXPR is not installed. Install it to speed up calculation")
    ne_flag = False

try:
    import numba as nb

    nb_flag = True
except:
    _logger.info("beam.py: module NUMBA is not installed. Install it to speed up calculation")
    nb_flag = False


def get_envelope(p_array, tws_i=None, bounds=None, slice=None, auto_disp=False):
    """
    Calculate Twiss parameters from a ParticleArray.

    This function processes particle data to compute Twiss parameters, with optional dispersion correction
    and selection of a particle subset.

    Parameters
    ----------
    p_array : ParticleArray
        Input particle array containing phase-space coordinates.
    tws_i : Twiss, optional
        Reference Twiss parameters for dispersion correction. If None and auto_disp is True,
        dispersion is estimated from the particle data. Default is None.
    bounds : list, optional
        Bounds as [left, right] in units of std(p_array.tau()) for selecting particles. Default is None.
    slice : str or None, optional
        Reference slice when bounds is set. If None, uses mean(tau). If 'Imax', uses maximum current slice.
    auto_disp : bool, optional
        If True and tws_i is None, estimate and subtract linear dispersion from the statistics of the particle array.
        Default is False.

    Returns
    -------
    Twiss
        Computed Twiss parameters for the (optionally filtered and corrected) particle array.
    """

    tau = p_array.tau()
    if bounds is not None:
        sig0 = np.std(tau)
        if slice == "Imax":
            charge = np.sum(p_array.q_array)
            B = s_to_cur(tau, 0.01 * sig0, charge, glb.speed_of_light)
            z0 = B[np.argmax(B[:, 1]), 0]
        else:
            z0 = np.mean(tau)
        inds = np.argwhere((z0 + sig0 * bounds[0] <= tau) * (tau <= z0 + sig0 * bounds[1]))
        p = p_array.p()[inds]
        x = p_array.x()[inds]
        px = p_array.px()[inds]
        y = p_array.y()[inds]
        py = p_array.py()[inds]
        tau = p_array.tau()[inds]
    else:
        p = p_array.p()
        x = p_array.x()
        px = p_array.px()
        y = p_array.y()
        py = p_array.py()

    tws = core.Twiss()
    tws.E = np.copy(p_array.E)
    tws.q = np.sum(p_array.q_array)
    tws.p = np.mean(p)

    # if less than 3 particles are left in the ParticleArray - return default (zero) Twiss()
    if len(x) < 3:
        _logger.warning("ParticleArray contains less than 3 particles. Moments are not calculated")
        return tws

    if tws_i is None:
        if auto_disp:
            mean_x, mean_px = np.mean(x), np.mean(px)
            mean_y, mean_py = np.mean(y), np.mean(py)
            mean_p = tws.p
            var_p = np.var(p)
            tws.Dx = np.mean((p - mean_p) * (x - mean_x)) / var_p
            tws.Dxp = np.mean((p - mean_p) * (px - mean_px)) / var_p
            tws.Dy = np.mean((p - mean_p) * (y - mean_y)) / var_p
            tws.Dyp = np.mean((p - mean_p) * (py - mean_py)) / var_p
    else:
        tws.Dx = tws_i.Dx
        tws.Dxp = tws_i.Dxp
        tws.Dy = tws_i.Dy
        tws.Dyp = tws_i.Dyp

    dx = tws.Dx * p
    dy = tws.Dy * p
    dpx = tws.Dxp * p
    dpy = tws.Dyp * p
    x = x - dx

    px = px - dpx

    y = y - dy
    py = py - dpy

    if ne_flag:
        local = {'px': px, 'py': py, 'p': p}
        factor = ne.evaluate('1. - p - 0.5 * p**2 + 0.5 * px**2 + 0.5 * py**2', local_dict=local)
        px = ne.evaluate('px * factor', local_dict={'px': px, 'factor': factor})
        py = ne.evaluate('py * factor', local_dict={'py': py, 'factor': factor})
    else:
        factor = 1. - p - 0.5 * p * p + 0.5 * px * px + 0.5 * py * py
        px = px * factor
        py = py * factor

    tws.x = np.mean(x)
    tws.y = np.mean(y)
    tws.px = np.mean(px)
    tws.py = np.mean(py)
    tws.tau = np.mean(tau)


    if ne_flag:
        tw_x = tws.x
        tw_y = tws.y
        tw_px = tws.px
        tw_py = tws.py
        tw_tau = tws.tau
        tw_p = tws.p
        tws.xx = np.mean(ne.evaluate('(x - tw_x)**2'))
        tws.xpx = np.mean(ne.evaluate('(x - tw_x) * (px - tw_px)'))
        tws.pxpx = np.mean(ne.evaluate('(px - tw_px) * (px - tw_px)'))
        tws.yy = np.mean(ne.evaluate('(y - tw_y)**2'))
        tws.ypy = np.mean(ne.evaluate('(y - tw_y) * (py - tw_py)'))
        tws.pypy = np.mean(ne.evaluate('(py - tw_py)**2'))
        tws.tautau = np.mean(ne.evaluate('(tau - tw_tau)**2'))

        tws.xy = np.mean(ne.evaluate('(x - tw_x) * (y - tw_y)'))
        tws.pxpy = np.mean(ne.evaluate('(px - tw_px) * (py - tw_py)'))
        tws.xpy = np.mean(ne.evaluate('(x - tw_x) * (py - tw_py)'))
        tws.ypx = np.mean(ne.evaluate('(y - tw_y) * (px - tw_px)'))
        tws.pp = np.mean(ne.evaluate('(p - tw_p)**2'))

    else:
        tws.xx = np.mean((x - tws.x) ** 2)
        tws.xpx = np.mean((x - tws.x) * (px - tws.px))
        tws.pxpx = np.mean((px - tws.px) ** 2)
        tws.yy = np.mean((y - tws.y) ** 2)
        tws.ypy = np.mean((y - tws.y) * (py - tws.py))
        tws.pypy = np.mean((py - tws.py) ** 2)
        tws.tautau = np.mean((tau - tws.tau) * (tau - tws.tau))

        tws.xy = np.mean((x - tws.x) * (y - tws.y))
        tws.pxpy = np.mean((px - tws.px) * (py - tws.py))
        tws.xpy = np.mean((x - tws.x) * (py - tws.py))
        tws.ypx = np.mean((y - tws.y) * (px - tws.px))
        tws.pp = np.mean((p - tws.p) ** 2)

    Sigma = np.array([[tws.xx, tws.xy, tws.xpx, tws.xpy],
                      [tws.xy, tws.yy, tws.ypx, tws.ypy],
                      [tws.xpx, tws.ypx, tws.pxpx, tws.pxpy],
                      [tws.xpy, tws.ypy, tws.pxpy, tws.pypy]])

    S = np.array([[0, 0, 1, 0],
                  [0, 0, 0, 1],
                  [-1, 0, 0, 0],
                  [0, -1, 0, 0]])
    # w, v = np.linalg.eig(np.dot(Sigma, S))

    tws.emit_x = np.sqrt(tws.xx * tws.pxpx - tws.xpx ** 2)
    tws.emit_y = np.sqrt(tws.yy * tws.pypy - tws.ypy ** 2)
    relgamma = p_array.E / glb.m_e_GeV
    relbeta = np.sqrt(1 - relgamma ** -2) if relgamma != 0 else 1.
    tws.emit_xn = tws.emit_x * relgamma * relbeta
    tws.emit_yn = tws.emit_y * relgamma * relbeta

    xx = tws.xx
    xpx = tws.xpx
    pxpx = tws.pxpx
    yy = tws.yy
    ypy = tws.ypy
    pypy = tws.pypy
    xy = tws.xy
    pxpy = tws.pxpy
    xpy = tws.xpy
    ypx = tws.ypx

    eigemit1 = np.sqrt(xpx ** 2 / 2 - (pxpx * xx) / 2 - pxpy * xy + xpy * ypx + ypy ** 2 / 2 - (pypy * yy) / 2
                       - 1 / 2 * np.sqrt(
        (0j - xpx ** 2 + pxpx * xx + 2 * pxpy * xy - 2 * xpy * ypx - ypy ** 2 + pypy * yy) ** 2
        - 4 * (pxpy ** 2 * xy ** 2 - pxpx * pypy * xy ** 2 + 2 * pypy * xpx * xy * ypx - 2 * pxpy * xpy * xy * ypx
               + xpy ** 2 * ypx ** 2 - pypy * xx * ypx ** 2 - 2 * pxpy * xpx * xy * ypy + 2 * pxpx * xpy * xy * ypy
               - 2 * xpx * xpy * ypx * ypy + 2 * pxpy * xx * ypx * ypy + xpx ** 2 * ypy ** 2 - pxpx * xx * ypy ** 2
               - pypy * xpx ** 2 * yy + 2 * pxpy * xpx * xpy * yy - pxpx * xpy ** 2 * yy - pxpy ** 2 * xx * yy
               + pxpx * pypy * xx * yy)))

    eigemit2 = (1 / np.sqrt(2)) * (np.sqrt(xpx ** 2 - pxpx * xx - 2 * pxpy * xy + 2 * xpy * ypx + ypy ** 2 - pypy * yy
                                           + np.sqrt(
        (xpx ** 2 - pxpx * xx - 2 * pxpy * xy + 2 * xpy * ypx + ypy ** 2 - pypy * yy) ** 2
        + 4 * (-2 * pypy * xpx * xy * ypx - xpy ** 2 * ypx ** 2 + pypy * xx * ypx ** 2 + 2 * xpx * xpy * ypx * ypy
               - xpx ** 2 * ypy ** 2 + pypy * xpx ** 2 * yy
               + 2 * pxpy * (xpy * xy * ypx + xpx * xy * ypy - xx * ypx * ypy - xpx * xpy * yy)
               + pxpy ** 2 * (-xy ** 2 + xx * yy)
               + pxpx * (pypy * xy ** 2 - 2 * xpy * xy * ypy + xx * ypy ** 2 + xpy ** 2 * yy - pypy * xx * yy + 0j)))))

    tws.eigemit_1 = eigemit1.imag  # w[0].imag
    tws.eigemit_2 = eigemit2.imag  # w[2].imag
    tws.beta_x = tws.xx / tws.emit_x
    tws.beta_y = tws.yy / tws.emit_y
    tws.alpha_x = -tws.xpx / tws.emit_x
    tws.alpha_y = -tws.ypy / tws.emit_y

    return tws


def get_current(p_array, num_bins=200, **kwargs):
    """
    Function calculates beam current from particleArray.

    :param p_array: particleArray
    :param charge: - None, OBSOLETE, charge of the one macro-particle.
                    If None, charge of the first macro-particle is used
    :param num_bins: number of bins
    :return s, I -  (np.array, np.array) - beam positions [m] and currents in [A]
    """
    if "charge" in kwargs:
        _logger.warning("argument 'charge' is obsolete use 'get_current(p_array, num_bins)' instead")
        charge = kwargs["charge"]
    else:
        charge = None
    weights = None
    if charge is None:
        weights = p_array.q_array
        charge = 1

    z = p_array.tau()
    hist, bin_edges = np.histogram(z, bins=num_bins, weights=weights)
    bin_edges = (bin_edges[:-1] + bin_edges[1:]) / 2.
    delta_Z = max(z) - min(z)
    delta_z = delta_Z / num_bins
    t_bins = delta_z / glb.speed_of_light
    return bin_edges, hist * charge / t_bins


def s2cur_auxil_py(A, xiA, C, N, I):
    for k in range(len(A)):
        i = I[k]
        if i > N - 1:
            i = N - 1
        C[i] = C[i] + xiA[k]
        C[i + 1] = C[i + 1] + (1 - xiA[k])


s2cur_auxil = s2cur_auxil_py if not nb_flag else nb.jit(nopython=True)(s2cur_auxil_py)


def s_to_cur(A, sigma, q0, v=glb.speed_of_light, ds=None, N=None):
    """
    Calculate beam current profile using first-order (CIC) deposition
    with optional Gaussian smoothing.

    Parameters
    ----------
    A : ndarray
        s-coordinates of particles [m].
    q0 : float
        Total bunch charge [C].
    v : float
        Mean velocity [m/s].
    sigma : float, optional
        Gaussian smoothing width [m]. If None, no smoothing is applied.
    ds : float, optional
        Bin size [m]. If given, overrides default binning.
    N : int, optional
        Number of bins. If given, overrides default binning.
        Mutually exclusive with `ds`.

    Returns
    -------
    B : ndarray, shape (N,2)
        Columns: [s [m], I(s) [A]].
    """

    Nsigma = 3

    # --- grid extent ---
    a = np.min(A)
    b = np.max(A)
    if sigma is not None:
        a -= Nsigma * sigma
        b += Nsigma * sigma

    # --- grid definition ---
    if ds is not None and N is not None:
        raise ValueError("Specify either ds or N, not both.")
    if ds is not None:
        N = int(np.ceil((b - a) / ds))
    elif N is not None:
        ds = (b - a) / N
    else:
        # default resolution if neither ds nor N is given
        if sigma is not None and sigma > 0:
            ds = 0.25 * sigma
        else:
            ds = (b - a) / 1000.0
        N = int(np.ceil((b - a) / ds))

    ds = (b - a) / N  # ensure consistency
    B = np.zeros((N + 1, 2))
    C = np.zeros(N + 1)

    B[:, 0] = np.arange(0, (N + 0.5) * ds, ds) + a
    N = N + 1
    cA = (A - a) / ds
    I = np.int_(np.floor(cA))
    xiA = 1 + I - cA
    s2cur_auxil(A, xiA, C, N, I)

    # --- Gaussian smoothing (optional) ---
    if sigma is not None and sigma > 0:
        K = int(np.floor(Nsigma * sigma / ds + 0.5))
        G = np.exp(-0.5 * (np.arange(-K, K + 1) * ds / sigma) ** 2)
        G = G / np.sum(G)
        B[:, 1] = beam_utils.convmode(C, G, 1)
    else:
        B[:, 1] = C

    # --- normalization ---
    koef = q0 * v / (ds * np.sum(B[:, 1]))
    B[:, 1] = koef * B[:, 1]
    return B


def signal_to_spectrum(s, w):
    """
    Compute the Fourier spectrum of a longitudinal signal using the
    convention exp(+i ω t).

    This is a generic utility: depending on the physical meaning of `w`,
    the output spectrum can represent an impedance, a current spectrum,
    or a bunching spectrum.

    Parameters
    ----------
    s : ndarray
        Longitudinal coordinate array [m]. Uniformly spaced.
    w : ndarray
        Signal defined on `s`. Examples:
          - wake potential W(s) [V/C] → spectrum gives impedance Z(ω) [Ω]
          - current profile I(s) [A]  → spectrum gives current spectrum I(ω)
          - normalized charge density → spectrum gives bunching factor b(k)

    Returns
    -------
    f : ndarray
        Frequency array [Hz] corresponding to the spectrum bins.
    y : ndarray (complex)
        Fourier spectrum of the input signal. Physical units depend on `w`.

    Notes
    -----
    - Uses FFT with exp(+i ω t) convention.
    - The relation between spatial coordinate s and time is t = s / c.
    - The FFT normalization includes a factor `dt` so that Parseval's theorem
      is satisfied and the units of `y` are consistent with the input `w`.
    """
    ds = s[1] - s[0]
    dt = ds / glb.speed_of_light
    n = len(s)
    f = 1 / dt * np.arange(0, n) / n
    #f = np.fft.fftfreq(n, d=dt)
    shift = 1  # optional phase shift, e.g. exp(1j* f * t0 * 2π)
    y = dt * np.fft.fft(w, n) * shift
    return f, y




def slice_analysis_py(x, xp, m_slice):
    """
    Function calculates moments and emittance - <x>, <xs>, <x^2>, <x*xs>, <xs^2>, np.sqrt(<x^2> * <xs^2> - <x*xs>^2)
    based on m_slice particles in moving window.
    NOTE: the coordinate must be sorted with respect to longitudinal coordinate.

    :param x: ndarray, 1st coordinate
    :param xp: ndarray, 2nd coordinate
    :param m_slice: M particles in moving window
    :return: list, [<x>, <xs>, <x^2>, <x*xs>, <xs^2>, np.sqrt(<x^2> * <xs^2> - <x*xs>^2)]
    """

    N = len(x)
    mx = np.zeros(N)
    mxs = np.zeros(N)
    mxx = np.zeros(N)
    mxxs = np.zeros(N)
    mxsxs = np.zeros(N)

    m = np.max(np.array([np.round(m_slice / 2), 1]))
    xc = np.cumsum(x)
    xsc = np.cumsum(xp)
    for i in range(N):
        n1 = int(max(0, i - m))
        n2 = int(min(N - 1, i + m))
        dq = n2 - n1  # window size
        mx[i] = (xc[n2] - xc[n1]) / dq  # average for over window per particle
        mxs[i] = (xsc[n2] - xsc[n1]) / dq

    x = x - mx
    xp = xp - mxs
    x2c = np.cumsum(x * x)
    xs2c = np.cumsum(xp * xp)
    xxsc = np.cumsum(x * xp)
    for i in range(N):
        n1 = int(max(0, i - m))
        n2 = int(min(N - 1, i + m))
        dq = n2 - n1
        mxx[i] = (x2c[n2] - x2c[n1]) / dq
        mxsxs[i] = (xs2c[n2] - xs2c[n1]) / dq
        mxxs[i] = (xxsc[n2] - xxsc[n1]) / dq

    emittx = np.sqrt(mxx * mxsxs - mxxs * mxxs)
    return [mx, mxs, mxx, mxxs, mxsxs, emittx]


slice_analysis = slice_analysis_py if not nb_flag else nb.jit(slice_analysis_py, nopython=True)



def slice_analysis_transverse(parray, Mslice, Mcur, p, iter):
    q1 = np.sum(parray.q_array)
    _logger.debug("slice_analysis_transverse: charge = " + str(q1))
    n = np.int_(parray.rparticles.size / 6)
    PD = parray.rparticles
    PD = beam_utils.sortcols(PD, row=4)

    z = np.copy(PD[4])
    mx, mxs, mxx, mxxs, mxsxs, emittx = slice_analysis(PD[0], PD[1], Mslice)

    my, mys, myy, myys, mysys, emitty = slice_analysis(PD[2], PD[3], Mslice)

    mm, mm, mm, mm, mm, emitty0 = beam_utils.moments(PD[2], PD[3])
    gamma0 = parray.E / glb.m_e_GeV
    emityn = emitty0 * gamma0
    mm, mm, mm, mm, mm, emitt0 = beam_utils.moments(PD[0], PD[1])
    emitxn = emitt0 * gamma0

    z, ind = np.unique(z, return_index=True)
    emittx = emittx[ind]
    emitty = emitty[ind]
    smin = min(z)
    smax = max(z)
    n = 1000
    hs = (smax - smin) / (n - 1)
    s = np.arange(smin, smax + hs, hs)
    ex = beam_utils.interp1(z, emittx, s)
    ey = beam_utils.interp1(z, emitty, s)

    ex = beam_utils.simple_filter(ex, p, iter) * gamma0 * 1e6
    ey = beam_utils.simple_filter(ey, p, iter) * gamma0 * 1e6

    sig0 = np.std(parray.tau())
    B = s_to_cur(z, Mcur * sig0, q1, glb.speed_of_light)
    I = beam_utils.interp1(B[:, 0], B[:, 1], s)
    return [s, I, ex, ey, gamma0, emitxn, emityn]


class SliceParameters:
    SP_TO_TWISS_NAMES: dict[str, str] = {"ex": "emit_x",
                                         "ey": "emit_y",
                                         "exn": "emit_xn",
                                         "eyn": "emit_yn",
                                         "mx": "x",
                                         "mxp": "xp",
                                         "my": "y",
                                         "myp": "yp",
                                         "me": "E",
                                         "beta_x": "beta_x",
                                         "beta_y": "beta_y",
                                         "alpha_x": "alpha_x",
                                         "alpha_y": "alpha_y",
                                         "gamma_x": "gamma_x",
                                         "gamma_y": "gamma_y"}

    VARIANCE_SP_NAMES: dict[str, str] = {"se": "pp",
                                         "sig_x": "xx",
                                         "sig_y": "yy",
                                         "sig_xp": "pxpx",
                                         "sig_yp": "pypy"}
    # SliceParameter energy is in units of eV whereas in Twiss
    # instances it should be units of GeV.  Maybe I missed some here.
    TWISS_UNITS_CONVERSION = {"E": 1e-6}

    def __init__(self):
        self.s = None
        self.I = None
        self.ex = None
        self.ey = None
        self.exn = None
        self.eyn = None
        self.me = None
        self.se = None
        self.gamma0 = None
        self.emitxn = None
        self.emityn = None

        # additional moments <x>, <xp>, <y>, <yp>, <p>
        self.mx = None
        self.mxp = None
        self.my = None
        self.myp = None
        self.mp = None

        self.sig_x = None
        self.sig_y = None
        self.sig_yp = None
        self.sig_xp = None

        # twiss slice parameters
        self.beta_x = None
        self.beta_y = None
        self.alpha_x = None
        self.alpha_y = None
        self.gamma_x = None
        self.gamma_y = None

    def extract_slice(self, index: int) -> core.Twiss:
        # Not all are added if there is no appropriate and unambiguous
        # analogue in the Twiss class
        rtwiss = core.Twiss()

        for slice_parameters_name, twiss_name in self.SP_TO_TWISS_NAMES.items():
            chosen_slice_value = getattr(self, slice_parameters_name)[index]
            setattr(rtwiss, twiss_name, chosen_slice_value)

        for slice_parameters_name, twiss_name in self.VARIANCE_SP_NAMES.items():
            chosen_slice_value = getattr(self, slice_parameters_name)[index] ** 2
            setattr(rtwiss, twiss_name, chosen_slice_value)

        for attr, factor in self.TWISS_UNITS_CONVERSION.items():
            value = getattr(rtwiss, attr)
            new_value = value * factor
            setattr(rtwiss, attr, new_value)

        return rtwiss


def global_slice_analysis_extended(parray, Mslice, Mcur, p, iter):
    """
    Function to calculate slice parameters

    :param parray: ParticleArray
    :param Mslice: 5000, nparticles in the slice
    :param Mcur: 0.01, smoothing parameters to calculate the beam current: smooth_param = m_std * np.std(p_array.tau())
    :param p: 2, filter parameter in the func: simple_filter
    :param iter: 2, filter parameter in the func: simple_filter
    :return: s, I, ex, ey, me, se, gamma0, emitxn, emityn
    """

    q1 = np.sum(parray.q_array)
    n = np.int_(parray.rparticles.size / 6)
    PD = parray.rparticles
    PD = beam_utils.sortcols(PD, row=4)

    z = np.copy(PD[4])
    mx, mxs, mxx, mxxs, mxsxs, emittx = slice_analysis(PD[0], PD[1], Mslice)

    my, mys, myy, myys, mysys, emitty = slice_analysis(PD[2], PD[3], Mslice)

    pc_0 = np.sqrt(parray.E ** 2 - glb.m_e_GeV ** 2)
    E1 = PD[5] * pc_0 + parray.E
    pc_1 = np.sqrt(E1 ** 2 - glb.m_e_GeV ** 2)
    mE, mEs, mEE, mEEs, mEsEs, emittE = slice_analysis(PD[4], pc_1 * 1e9, Mslice)

    mE = mEs  # mean energy
    sE = np.sqrt(mEsEs)  # energy spread
    sig0 = np.std(parray.tau())  # std pulse duration
    B = s_to_cur(z, Mcur * sig0, q1, glb.speed_of_light)
    gamma0 = parray.E / glb.m_e_GeV
    _, _, _, _, _, emitty0 = beam_utils.moments(PD[2], PD[3])
    emityn = emitty0 * gamma0
    _, _, _, _, _, emitt0 = beam_utils.moments(PD[0], PD[1])
    emitxn = emitt0 * gamma0

    z, ind = np.unique(z, return_index=True)
    emittx = emittx[ind]
    emitty = emitty[ind]
    sE = sE[ind]
    mE = mE[ind]
    smin = min(z)
    smax = max(z)
    n = 1000
    hs = (smax - smin) / (n - 1)
    s = np.arange(smin, smax + hs, hs)
    ex = beam_utils.interp1(z, emittx, s)
    ey = beam_utils.interp1(z, emitty, s)
    se = beam_utils.interp1(z, sE, s)
    me = beam_utils.interp1(z, mE, s)
    ex = beam_utils.simple_filter(ex, p, iter) * gamma0 * 1e6
    ey = beam_utils.simple_filter(ey, p, iter) * gamma0 * 1e6
    se = beam_utils.simple_filter(se, p, iter)
    me = beam_utils.simple_filter(me, p, iter)

    I = beam_utils.interp1(B[:, 0], B[:, 1], s)

    return [s, I, ex, ey, me, se, gamma0, emitxn, emityn]


def global_slice_analysis(parray, nparts_in_slice=5000, smooth_param=0.01, filter_base=2, filter_iter=2):
    """
    Function to calculate slice parameters

    :param parray: ParticleArray
    :param nparts_in_slice: 5000, nparticles in the slice (in moving window)
    :param smooth_param: 0.01, smoothing parameters to calculate the beam current: smooth_param = m_std * np.std(p_array.tau())
    :param filter_base: 2, filter parameter in the func: simple_filter
    :param filter_iter: 2, filter parameter in the func: simple_filter
    :return: SliceParameters,
    """
    n = 1000  # number of points

    slc = SliceParameters()

    q1 = np.sum(parray.q_array)

    PD = parray.rparticles
    PD = beam_utils.sortcols(PD, row=4)

    z = np.copy(PD[4])
    mx, mxs, mxx, mxxs, mxsxs, emittx = slice_analysis(PD[0], PD[1], nparts_in_slice)

    my, mys, myy, myys, mysys, emitty = slice_analysis(PD[2], PD[3], nparts_in_slice)
    pc_0 = np.sqrt(parray.E ** 2 - glb.m_e_GeV ** 2)
    E1 = PD[5] * pc_0 + parray.E
    pc_1 = np.sqrt(E1 ** 2 - glb.m_e_GeV ** 2)

    mE, mEs, mEE, mEEs, mEsEs, emittE = slice_analysis(PD[4], pc_1 * 1e9, nparts_in_slice)

    mE = mEs  # mean energy
    sE = np.sqrt(mEsEs)  # energy spread
    sig0 = np.std(parray.tau())  # std pulse duration
    B = s_to_cur(z, smooth_param * sig0, q1, glb.speed_of_light)
    gamma0 = parray.E / glb.m_e_GeV
    _, _, _, _, _, emitty0 = beam_utils.moments(PD[2], PD[3])
    slc.emityn = emitty0 * gamma0
    _, _, _, _, _, emitt0 = beam_utils.moments(PD[0], PD[1])
    slc.emitxn = emitt0 * gamma0

    _, mp, _, _, _, _ = slice_analysis(PD[4], PD[5], nparts_in_slice)

    z, ind = np.unique(z, return_index=True)

    emittx = emittx[ind]
    emitty = emitty[ind]
    sE = sE[ind]
    mE = mE[ind]
    sig_x = np.sqrt(mxx[ind])
    sig_y = np.sqrt(myy[ind])

    sig_xp = np.sqrt(mxsxs[ind])
    sig_yp = np.sqrt(mysys[ind])

    smin = min(z)
    smax = max(z)

    hs = (smax - smin) / (n - 1)
    s = np.linspace(smin, smax, num=n)
    ex = beam_utils.interp1(z, emittx, s)
    ey = beam_utils.interp1(z, emitty, s)
    se = beam_utils.interp1(z, sE, s)
    me = beam_utils.interp1(z, mE, s)
    slc.ex = beam_utils.simple_filter(ex, filter_base, filter_iter)
    slc.ey = beam_utils.simple_filter(ey, filter_base, filter_iter)
    slc.exn = slc.ex * gamma0
    slc.eyn = slc.ey * gamma0
    slc.se = beam_utils.simple_filter(se, filter_base, filter_iter)
    slc.me = beam_utils.simple_filter(me, filter_base, filter_iter)

    slc.I = beam_utils.interp1(B[:, 0], B[:, 1], s)

    mxpx = mxxs[ind]
    mypy = myys[ind]
    xpx_m = beam_utils.interp1(z, mxpx, s)
    ypy_m = beam_utils.interp1(z, mypy, s)
    x_px = beam_utils.simple_filter(xpx_m, filter_base, filter_iter)
    y_py = beam_utils.simple_filter(ypy_m, filter_base, filter_iter)
    # additional moments <x>, <xp>, <y>, <yp>, <p>
    mx = mx[ind]
    mxs = mxs[ind]
    my = my[ind]
    mys = mys[ind]

    xm = beam_utils.interp1(z, mx, s)
    xpm = beam_utils.interp1(z, mxs, s)
    ym = beam_utils.interp1(z, my, s)
    ypm = beam_utils.interp1(z, mys, s)

    sig_x = beam_utils.interp1(z, sig_x, s)
    sig_y = beam_utils.interp1(z, sig_y, s)

    sig_xp = beam_utils.interp1(z, sig_xp, s)
    sig_yp = beam_utils.interp1(z, sig_yp, s)

    slc.mx = beam_utils.simple_filter(xm, filter_base, filter_iter)
    slc.mxp = beam_utils.simple_filter(xpm, filter_base, filter_iter)
    slc.my = beam_utils.simple_filter(ym, filter_base, filter_iter)
    slc.myp = beam_utils.simple_filter(ypm, filter_base, filter_iter)

    slc.sig_x = beam_utils.simple_filter(sig_x, filter_base, filter_iter)
    slc.sig_y = beam_utils.simple_filter(sig_y, filter_base, filter_iter)

    slc.sig_xp = beam_utils.simple_filter(sig_xp, filter_base, filter_iter)
    slc.sig_yp = beam_utils.simple_filter(sig_yp, filter_base, filter_iter)

    # twiss
    # np.full(n, np.nan)
    slc.beta_x = np.divide(slc.sig_x ** 2, slc.ex, out=np.zeros_like(slc.sig_x), where=slc.ex != 0)
    slc.beta_y = np.divide(slc.sig_y ** 2, slc.ey, out=np.zeros_like(slc.sig_y), where=slc.ey != 0)
    slc.alpha_x = -x_px / slc.ex
    slc.alpha_y = -y_py / slc.ey
    slc.gamma_x = np.divide(1 + slc.alpha_x ** 2, slc.beta_x, out=np.zeros_like(slc.alpha_x), where=slc.beta_x != 0)
    slc.gamma_y = np.divide(1 + slc.alpha_y ** 2, slc.beta_y, out=np.zeros_like(slc.alpha_y), where=slc.beta_y != 0)
    mp = mp[ind]
    mp = beam_utils.interp1(z, mp, s)
    slc.mp = beam_utils.simple_filter(mp, filter_base, filter_iter)

    slc.s = s
    slc.gamma0 = gamma0
    return slc


def bunching_spectrum(z, Lwin=None, M=4096):
    """
    Compute |b(k)|^2 spectrum from 1D particle positions using FFT histogram.

    Parameters
    ----------
    z : array_like
        Particle positions (e.g. tau).
    Lwin : float or None
        Total window length. If None, use span of z (max - min).
        Larger L = higher k-resolution but more zero-padding.
    M : int
        Number of grid points for FFT.

    Returns
    -------
    k : ndarray
        Wavenumbers [1/m].
    bk2 : ndarray
        |b(k)|^2 spectrum.
    """
    N = len(z)
    zmin, zmax = z.min(), z.max()
    if Lwin is None:
        Lwin = zmax - zmin
    # center window instead of forcing start at z.min()
    zmid = 0.5*(zmin + zmax)
    edges = np.linspace(zmid - Lwin/2, zmid + Lwin/2, M+1)
    n, _ = np.histogram(z, bins=edges)
    n_hat = np.fft.rfft(n)
    b_k = n_hat / N
    k = 2*np.pi*np.arange(len(b_k)) / Lwin
    # power bunching np.abs(b_k)**2
    return k, b_k


def bunching_at_klist_2(z, k_list):
    z = np.asarray(z)
    b = []
    for k in k_list:
        b.append(np.exp(-1j * k * z).mean())
    return np.array(b)


def bunching_at_klist_py(z, k_list):
    N = z.size
    out = np.empty(k_list.size, dtype=np.complex128)
    for i in range(k_list.size):
        s = 0.0 + 0.0j
        k = k_list[i]
        for j in range(N):
            s += np.exp(-1j * k * z[j])
        out[i] = s / N
    return out


def bunching_at_klist_np(z, k_list):
    """
    Pure NumPy implementation of
    b(k) = <exp(-i k z)>

    Parameters
    ----------
    z : array, shape (N,)
        Longitudinal coordinates
    k_list : array, shape (M,)
        Wavenumbers

    Returns
    -------
    out : array, shape (M,)
        Complex bunching factor
    """
    z = np.asarray(z)
    k_list = np.asarray(k_list)
    return np.exp(-1j * np.outer(k_list, z)).mean(axis=1)


# final dispatch
if nb_flag:
    bunching_at_klist = nb.jit(
        nopython=True, cache=True, fastmath=True
    )(bunching_at_klist_py)
else:
    bunching_at_klist = bunching_at_klist_np


def spectrum_to_z(spectrum_k, M, Lwin):
    """
    Transform e.g. Δδ(k) spectrum back to real-space Δδ(z).

    Parameters
    ----------
    spectrum_k : ndarray
        Complex spectrum (length M//2+1, from rfft).
    M : int
        Number of bins used in FFT (same as in bunching_spectrum).
    Lwin : float
        Window length [m] (same as in bunching_spectrum).

    Returns
    -------
    z : ndarray
        Bin centers [m].
    signal_z : ndarray
        Real-space Δδ(z) on the same grid.
    """
    signal_z = np.fft.irfft(spectrum_k, n=M)
    dz = Lwin / M
    z = np.linspace(-Lwin/2, Lwin/2 - dz, M)  # centered window
    return z, signal_z.real


def compute_bunching(z, sigma, q0, v, ds=None, N=None):
    """
    Compute bunching factor spectrum from particle coordinates.

    Parameters
    ----------
    z : ndarray
        Particle longitudinal positions [m].
    q0 : float
        Bunch charge [C].
    v : float
        Mean velocity [m/s].
    sigma : float
        Smoothing parameter [m].
    ds : float, optional
        Bin size [m].
    N : int, optional
        Number of bins.

    Returns
    -------
    k : ndarray
        Wavenumbers [1/m].
    b_k : ndarray (complex)
        Bunching factor spectrum.
    """
    B = s_to_cur(z, sigma=sigma, q0=q0, v=v, ds=ds, N=N)

    s, I_s = B[:,0], B[:,1]
    s_pad, I_s_pad = zero_pad_signal(s, I_s)

    f, I_f = signal_to_spectrum(s_pad, I_s_pad)

    k = 2*np.pi*f / v
    b_k = I_f / q0

    return k, b_k


def zero_pad_signal(s, I_s):
    """
    Zero-pad signal to double its length (append zeros of same length).

    Parameters
    ----------
    s : ndarray
        Coordinate array [m], uniformly spaced.
    I_s : ndarray
        Signal samples on grid s (e.g. current profile).

    Returns
    -------
    s_pad : ndarray
        Extended coordinate array [m].
    I_s_pad : ndarray
        Zero-padded signal.
    """
    nb = len(s)
    ds = s[1] - s[0]

    # New extended arrays
    s_pad = np.arange(0, 2*nb) * ds + s[0]
    I_s_pad = np.zeros_like(s_pad)
    I_s_pad[:nb] = I_s

    return s_pad, I_s_pad