import numpy as np
import pandas as pd
from typing import Iterable
import ocelot.common.globals as glb
from scipy.special import factorial
from scipy.signal import savgol_filter
from ocelot.common.math_op import find_nearest_idx
from ocelot.common.ocelog import *
from ocelot.cpbd.reswake import pipe_wake
import copy

_logger = logging.getLogger(__name__)


class Twiss:
    """
    class - container for twiss parameters
    """

    def __init__(self, beam=None, **kwargs):

        self._emit_xn = kwargs.get("emit_xn", 0.)
        self._emit_yn = kwargs.get("emit_yn", 0.)
        self._pending_emit_x = None
        self._pending_emit_y = None
        self._E = kwargs.get("E", 0.0)  # ref the beam energy in [GeV]

        # Apply emit_x and emit_y properly (only if present!)
        if "emit_x" in kwargs:
            self.emit_x = kwargs["emit_x"]
        if "emit_y" in kwargs:
            self.emit_y = kwargs["emit_y"]

        if "emit_xn" in kwargs:
            self._emit_xn = kwargs["emit_xn"]
        if "emit_yn" in kwargs:
            self._emit_yn = kwargs["emit_yn"]

        self.eigemit_1 = 0.
        self.eigemit_2 = 0.
        self._beta_x = kwargs.get("beta_x", 0.)
        self._beta_y = kwargs.get("beta_y", 0.)
        self._alpha_x = kwargs.get("alpha_x", 0.)
        self._alpha_y = kwargs.get("alpha_y", 0.)
        self.Dx = kwargs.get("Dx", 0.)
        self.Dy = kwargs.get("Dy", 0.)
        self.Dxp = kwargs.get("Dxp", 0.)
        self.Dyp = kwargs.get("Dyp", 0.)
        self.mux = kwargs.get("mux", 0.)  # phase advance
        self.muy = kwargs.get("muy", 0.)  # phase advance

        # parameters below in the most cases are calculated from the ParticleArray object
        # during tracking (see func 'get_envelop()')


        self.s = kwargs.get("s", 0.0)  # position along the reference trajectory [m]
        self.q = kwargs.get("q", 0.0)  # charge of the whole beam [C]

        # moments
        self.x = kwargs.get("x", 0.0)
        self.y = kwargs.get("y", 0.0)
        self.p = kwargs.get("p", 0.0)
        self.tau = kwargs.get("tau", 0.0)
        self.xp = kwargs.get("xp", 0.0)
        self.yp = kwargs.get("yp", 0.0)
        self.xx = kwargs.get("xx", 0.)
        self.xpx = kwargs.get("xpx", 0.)
        self.pxpx = kwargs.get("pxpx", 0.)
        self.yy = kwargs.get("yy", 0.)
        self.ypy = kwargs.get("ypy", 0.)
        self.pypy = kwargs.get("pypy", 0.)
        self.tautau = kwargs.get("tautau", 0.)
        self.xy = kwargs.get("xy", 0.)
        self.pxpy = kwargs.get("pxpy", 0.)
        self.xpy = kwargs.get("xpy", 0.)
        self.ypx = kwargs.get("ypx", 0.)
        self.pp = kwargs.get("pp", 0)

        self.id = kwargs.get("id", "")

        if isinstance(beam, (Twiss, Beam)):

            self._emit_xn = beam.emit_xn
            self._emit_yn = beam.emit_yn

            self._beta_x = beam.beta_x
            self._beta_y = beam.beta_y
            self._alpha_x = beam.alpha_x
            self._alpha_y = beam.alpha_y
            self.Dx = beam.Dx
            self.Dy = beam.Dy
            self.Dxp = beam.Dxp
            self.Dyp = beam.Dyp
            self.x = beam.x
            self.y = beam.y
            self.xp = beam.xp
            self.yp = beam.yp
            self.E = beam.E

    @property
    def beta_x(self):
        return self._beta_x

    @beta_x.setter
    def beta_x(self, value):
        self._beta_x = value

    @property
    def alpha_x(self):
        return self._alpha_x

    @alpha_x.setter
    def alpha_x(self, value):
        self._alpha_x = value

    @property
    def gamma_x(self):
        if self._beta_x != 0:
            return (1 + self._alpha_x ** 2) / self._beta_x
        return 0

    @property
    def beta_y(self):
        return self._beta_y

    @beta_y.setter
    def beta_y(self, value):
        self._beta_y = value

    @property
    def alpha_y(self):
        return self._alpha_y

    @alpha_y.setter
    def alpha_y(self, value):
        self._alpha_y = value

    @property
    def gamma_y(self):
        if self._beta_y != 0:
            return (1 + self._alpha_y ** 2) / self._beta_y
        return 0

    @property
    def E(self):
        return self._E

    @E.setter
    def E(self, value):
        self._E = value
        # If emit_x or emit_y were set before E, now recompute their xn versions
        if self._pending_emit_x is not None:
            self.emit_x = self._pending_emit_x  # triggers recompute
            self._pending_emit_x = None
        if self._pending_emit_y is not None:
            self.emit_y = self._pending_emit_y
            self._pending_emit_y = None

    @property
    def relgamma(self):
        return self.E / glb.m_e_GeV if self.E != 0 else 0.

    @property
    def relbeta(self):
        g = self.relgamma
        return np.sqrt(1 - g ** -2) if g > 0 else 1.

    # --- X Plane ---
    @property
    def emit_xn(self):
        return self._emit_xn

    @emit_xn.setter
    def emit_xn(self, value):
        self._emit_xn = value

    @property
    def emit_x(self):
        rb = self.relbeta * self.relgamma
        return self._emit_xn / rb if rb != 0 else 0.

    @emit_x.setter
    def emit_x(self, value):
        rb = self.relbeta * self.relgamma
        if rb == 0:
            self._pending_emit_x = value
        else:
            self._emit_xn = value * rb

    # --- Y Plane ---
    @property
    def emit_yn(self):
        return self._emit_yn

    @emit_yn.setter
    def emit_yn(self, value):
        self._emit_yn = value

    @property
    def emit_y(self):
        rb = self.relbeta * self.relgamma
        return self._emit_yn / rb if rb != 0 else 0.

    @emit_y.setter
    def emit_y(self, value):
        rb = self.relbeta * self.relgamma
        if rb == 0:
            self._pending_emit_y = value
        else:
            self._emit_yn = value * rb

    @property
    def sigma_x(self):
        if self.emit_x <= 0 or self.beta_x <= 0:
            return 0
        return float(np.sqrt(self.emit_x * self.beta_x))

    @property
    def sigma_y(self):
        if self.emit_y <= 0 or self.beta_y <= 0:
            return 0
        return float(np.sqrt(self.emit_y * self.beta_y))


    def multiply_with_tm(self, tm: 'TransferMap', length):
        tws = self.map_x_twiss(tm)
        tws.s = self.s + length
        return tws

    @staticmethod
    def track(R, tws0):
        R00, R01, R10, R11 = R[0, 0], R[0, 1], R[1, 0], R[1, 1]
        R22, R23, R32, R33 = R[2, 2], R[2, 3], R[3, 2], R[3, 3]
        R05, R15, R25, R35 = R[0, 5], R[1, 5], R[2, 5], R[3, 5]

        tws = Twiss(tws0)
        tws.p = tws0.p

        tws.beta_x = R00 ** 2 * tws0.beta_x - 2 * R00 * R01 * tws0.alpha_x + R01 ** 2 * tws0.gamma_x
        tws.beta_y = R22 ** 2 * tws0.beta_y - 2 * R22 * R23 * tws0.alpha_y + R23 ** 2 * tws0.gamma_y

        tws.alpha_x = -R00 * R10 * tws0.beta_x + (R01 * R10 + R11 * R00) * tws0.alpha_x - R01 * R11 * tws0.gamma_x
        tws.alpha_y = -R22 * R32 * tws0.beta_y + (R23 * R32 + R33 * R22) * tws0.alpha_y - R23 * R33 * tws0.gamma_y

        tws.Dx = R00 * tws0.Dx + R01 * tws0.Dxp + R05
        tws.Dxp = R10 * tws0.Dx + R11 * tws0.Dxp + R15
        tws.Dy = R22 * tws0.Dy + R23 * tws0.Dyp + R25
        tws.Dyp = R32 * tws0.Dy + R33 * tws0.Dyp + R35

        d_mux = np.arctan2(R01, R00 * tws0.beta_x - R01 * tws0.alpha_x)
        if d_mux < 0: d_mux += np.pi
        tws.mux = tws0.mux + d_mux

        d_muy = np.arctan2(R23, R22 * tws0.beta_y - R23 * tws0.alpha_y)
        if d_muy < 0: d_muy += np.pi
        tws.muy = tws0.muy + d_muy

        return tws

    def map_x_twiss(self, tm):
        E = self.E
        M = tm.get_params(energy=E).get_rotated_R()
        zero_tol = 1.e-10
        if abs(tm.get_delta_e()) > zero_tol:
            Ei = self.E
            Ef = self.E + tm.get_delta_e()
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

        tws = self.track(M, self)
        tws.E = E
        return tws

    def __str__(self):
        val = ""
        val += "emit_x  = " + str(self.emit_x) + "\n"
        val += "emit_y  = " + str(self.emit_y) + "\n"
        val += "emit_xn  = " + str(self.emit_xn) + "\n"
        val += "emit_yn  = " + str(self.emit_yn) + "\n"
        val += "beta_x  = " + str(self.beta_x) + "\n"
        val += "beta_y  = " + str(self.beta_y) + "\n"
        val += "alpha_x = " + str(self.alpha_x) + "\n"
        val += "alpha_y = " + str(self.alpha_y) + "\n"
        val += "Dx      = " + str(self.Dx) + "\n"
        val += "Dy      = " + str(self.Dy) + "\n"
        val += "Dxp     = " + str(self.Dxp) + "\n"
        val += "Dyp     = " + str(self.Dyp) + "\n"
        val += "mux     = " + str(self.mux) + "\n"
        val += "muy     = " + str(self.muy) + "\n"
        val += "nu_x    = " + str(self.mux / 2. / np.pi) + "\n"
        val += "nu_y    = " + str(self.muy / 2. / np.pi) + "\n"
        val += "E       = " + str(self.E) + "\n"
        val += "s        = " + str(self.s) + "\n"
        return val

    def to_series(self) -> pd.Series:
        """Return this Twiss instance as an equivalent Pandas Series instance."""
        keys = [
            'emit_x', 'emit_y', 'emit_xn', 'emit_yn', 'eigemit_1', 'eigemit_2',
            'beta_x', 'beta_y', 'alpha_x', 'alpha_y',
            'Dx', 'Dy', 'Dxp', 'Dyp',
            'mux', 'muy',
            'E', 's', 'q',
            'x', 'y', 'xp', 'yp', 'p', 'tau',
            'xx', 'xpx', 'pxpx',
            'yy', 'ypy', 'pypy',
            'tautau', 'xy', 'xpy', 'ypx', 'pxpy', 'pp',
            'id'
        ]
        data = {k: getattr(self, k) for k in keys if hasattr(self, k)}
        return pd.Series(data)

    @classmethod
    def from_series(cls, series: pd.Series):
        result = cls()
        for key, value in series.items():
            if hasattr(result, key):
                setattr(result, key, np.squeeze(value))
        return result


def twiss_iterable_to_df(twisses: Iterable[Twiss]) -> pd.DataFrame:
    """Convert an iterable of Twiss instances to a single DataFrame with the columns as
    keys

    :param twisses: iterable of twisses to be converted to a pandas DataFrame.
    """
    return pd.DataFrame(data=(t.to_series() for t in twisses))


class Beam:
    def __init__(self, x=0, xp=0, y=0, yp=0):
        # initial conditions
        self.x = x  # [m]
        self.y = y  # [m]
        self.xp = xp  # xprime [rad]
        self.yp = yp  # yprime [rad]

        self.E = 0.0  # electron energy [GeV]
        self.sigma_E = 0.0  # Energy spread [GeV]
        self.I = 0.0  # beam current [A]
        self.emit_x = 0.0  # horizontal emittance [m rad]
        self.emit_y = 0.0  # horizontal emittance [m rad]
        self.tlen = 0.0  # bunch length (rms) in fsec

        # twiss parameters
        self.beta_x = 0.0
        self.beta_y = 0.0
        self.alpha_x = 0.0
        self.alpha_y = 0.0
        self.Dx = 0.0
        self.Dy = 0.0
        self.Dxp = 0.0
        self.Dyp = 0.0

        self.shape = 'gaussian'  # of 'flattop'
        self.filePath = ''

    properties = ['g', 'dg', 'emit_xn', 'emit_yn', 'p', 'pz', 'px', 'py']

    @property
    def g(self):
        return self.E / glb.m_e_GeV

    @g.setter
    def g(self, value):
        self.E = value * glb.m_e_GeV

    @property
    def dg(self):
        return self.sigma_E / glb.m_e_GeV

    @dg.setter
    def dg(self, value):
        self.sigma_E = value * glb.m_e_GeV

    @property
    def emit_xn(self):
        return self.emit_x * self.g

    @emit_xn.setter
    def emit_xn(self, value):
        self.emit_x = value / self.g

    @property
    def emit_yn(self):
        return self.emit_y * self.g

    @emit_yn.setter
    def emit_yn(self, value):
        self.emit_y = value / self.g

    @property
    def p(self):
        return np.sqrt(self.g ** 2 - 1)

    @p.setter
    def p(self, value):
        self.g = np.sqrt(value ** 2 + 1)

    @property
    def pz(self):
        return self.p / (self.xp ** 2 + self.yp ** 2 + 1)

    @pz.setter
    def pz(self, value):
        self.p = value * (self.xp ** 2 + self.yp ** 2 + 1)

    @property
    def px(self):
        return self.p * self.xp

    @px.setter
    def px(self, value):
        self.xp = value / self.p

    @property
    def py(self):
        return self.p * self.yp

    @py.setter
    def py(self, value):
        self.yp = value / self.p

    def len(self):
        return 1

    def to_array(self, nslice=100, window_len=None):

        if self.tlen in [None, 0] and window_len is None:
            raise ValueError('both self.tlen and window_len are not set')
        if window_len is None:
            if self.shape == 'gaussian':
                window_len = self.tlen * 1e-15 * glb.speed_of_light * 6  # sigmas
            elif self.shape == 'flattop':
                window_len = self.tlen * 1e-15 * glb.speed_of_light * 2  # fwhm
            else:
                raise ValueError('Beam() shape can be either "gaussian" or "flattop"')

        beam_arr = BeamArray(nslice)
        for param in beam_arr.params():
            if hasattr(self, param) and len(getattr(beam_arr, param)) == nslice:
                setattr(beam_arr, param, np.ones(nslice) * getattr(self, param))
        beam_arr.s = np.linspace(0, window_len, nslice)

        if self.tlen not in [None, 0, np.inf]:
            beam_slen = self.tlen * 1e-15 * glb.speed_of_light
            Ipeak_pos = (np.amax(beam_arr.s) - np.amin(beam_arr.s)) / 2
            if self.shape == 'gaussian':
                beam_arr.I = self.I * np.exp(-(beam_arr.s - Ipeak_pos) ** 2 / (2 * beam_slen ** 2))
            elif self.shape == 'flattop':
                beam_arr.I = np.ones_like(beam_arr.s) * self.I
                beam_arr.I[abs(beam_arr.s - Ipeak_pos) > beam_slen / 2] = 0
        return beam_arr

    def sizes(self):
        if self.beta_x != 0:
            self.gamma_x = (1. + self.alpha_x ** 2) / self.beta_x
        else:
            self.gamma_x = 0.

        if self.beta_y != 0:
            self.gamma_y = (1. + self.alpha_y ** 2) / self.beta_y
        else:
            self.gamma_y = 0.

        self.sigma_x = np.sqrt((self.sigma_E / self.E * self.Dx) ** 2 + self.emit_x * self.beta_x)
        self.sigma_y = np.sqrt((self.sigma_E / self.E * self.Dy) ** 2 + self.emit_y * self.beta_y)
        self.sigma_xp = np.sqrt((self.sigma_E / self.E * self.Dxp) ** 2 + self.emit_x * self.gamma_x)
        self.sigma_yp = np.sqrt((self.sigma_E / self.E * self.Dyp) ** 2 + self.emit_y * self.gamma_y)

    # def print_sizes(self):
    # self.sizes()
    # print("sigma_E/E and Dx/Dy : ", self.sigma_E/self.E, "  and ", self.Dx, "/",self.Dy, " m")
    # print("emit_x/emit_y     : ",  self.emit_x*1e9, "/",self.emit_y*1e9, " nm-rad")
    # print("sigma_x/y         : ", self.sigma_x*1e6, "/", self.sigma_y*1e6, " um")
    # print("sigma_xp/yp       : ", self.sigma_xp*1e6, "/", self.sigma_yp*1e6, " urad")

class BeamArray(Beam):

        def __init__(self, nslice=0):
            super().__init__()
            # initial conditions
            self.s = np.zeros(nslice)
            self.x = np.zeros(nslice)  # [m]
            self.y = np.zeros(nslice)  # [m]
            self.xp = np.zeros(nslice)  # xprime [rad]
            self.yp = np.zeros(nslice)  # yprime [rad]

            self.E = np.zeros(nslice)  # electron energy [GeV]
            self.sigma_E = np.zeros(nslice)  # Energy spread [GeV]
            self.I = np.zeros(nslice)  # beam current [A]
            self.emit_x = np.zeros(nslice)  # horizontal emittance [m rad]
            self.emit_y = np.zeros(nslice)  # horizontal emittance [m rad]
            # self.tlen = 0.0         # bunch length (rms) in fsec

            # twiss parameters
            self.beta_x = np.zeros(nslice)
            self.beta_y = np.zeros(nslice)
            self.alpha_x = np.zeros(nslice)
            self.alpha_y = np.zeros(nslice)
            self.Dx = np.zeros(nslice)
            self.Dy = np.zeros(nslice)
            self.Dxp = np.zeros(nslice)
            self.Dyp = np.zeros(nslice)

            self.eloss = np.zeros(nslice)

            del self.shape  # inherited, not applicable

        def idx_max(self):
            idx = np.where(self.I == np.nanmax(self.I))[0]
            if len(idx) == 1:
                return idx[0]
            else:
                return idx[int(len(idx) / 2)]

        def len(self):
            return np.size(self.s)

        def charge(self):
            return np.trapz(self.I, self.s / glb.speed_of_light)  # C

        def params(self):
            l = self.len()
            attrs = []
            for attr in dir(self):
                if attr.startswith('__') or attr in self.properties:
                    continue
                # if callable(getattr(self,attr)):
                # continue
                if np.size(getattr(self, attr)) == l:
                    attrs.append(attr)
            return attrs

        def sort(self):
            _logger.debug('sorting beam slices')
            inds = self.s.argsort()
            for attr in self.params():
                _logger.log(5, ind_str + 'sorting {:}'.format(str(attr)))
                values = getattr(self, attr)
                _logger.log(5, ind_str + 'size {:}'.format(values.size))
                setattr(self, attr, values[inds])

        def equidist(self, ds=None):
            dsarr = (self.s - np.roll(self.s, 1))[1:]
            dsm = np.mean(dsarr)
            if (np.abs(dsarr - dsm) / dsm > 1 / 1000).any():
                if ds is None:
                    s_new = np.linspace(np.amin(self.s), np.amax(self.s), self.len())
                else:
                    s_new = np.arange(np.amin(self.s), np.amax(self.s), ds)
                for attr in self.params():
                    if attr == 's':
                        continue
                    val = getattr(self, attr)
                    val = np.interp(s_new, self.s, val)
                    setattr(self, attr, val)
                self.s = s_new
            self.ds = dsm

        def smear(self, sw):
            _logger.debug('smearing the beam by {:.2e} m'.format(sw))
            self.equidist()
            sn = (sw / self.ds).astype(int)
            if sn < 2:
                return
            if not sn % 2:
                sn += 1

            for attr in self.params():
                if attr == 's':
                    continue
                val = getattr(self, attr)
                val = savgol_filter(val, sn, 2, mode='nearest')

                if attr in ['E', 'I', 'beta_x', 'beta_y', 'emit_x', 'emit_y', 'sigma_E']:
                    val[val < 0] = 0
                #    val = convolve(val,spike,mode='same')
                setattr(self, attr, val)

        def get_s(self, s):
            idx = find_nearest_idx(self.s, s)
            return self[idx]

        def get_E(self):
            return self.E[self.idx_max()]

        def set_E(self, E_GeV):
            idx = self.idx_max()
            self.E += (E_GeV - self.get_E())

        def center(self, s):
            beam_s = self.get_s(s)
            self.x -= beam_s.x
            self.xp -= beam_s.xp
            self.y -= beam_s.y
            self.yp -= beam_s.yp

        def cut_empty_I(self, thresh=0.01):
            idx = np.where(self.I <= self.I.max() * thresh)
            del self[idx]

        def start_at_0(self):
            self.s -= self.s.min()

        def cleanup(self):
            self.cut_empty_I()
            self.sort()
            self.equidist()
            self.start_at_0()

        def __getitem__(self, index):
            l = self.len()
            if index.__class__ is not slice:
                if index > l:
                    raise IndexError('slice index out of range')
                beam_slice = Beam()
            else:
                beam_slice = copy.deepcopy(self)

            for attr in dir(self):
                if attr.startswith('__') or callable(getattr(self, attr)):
                    continue
                value = getattr(self, attr)
                if np.size(value) == l:
                    setattr(beam_slice, attr, value[index])
                else:
                    setattr(beam_slice, attr, value)
            return beam_slice

        def __delitem__(self, index):
            l = self.len()
            for attr in self.params():
                if attr.startswith('__') or callable(getattr(self, attr)):
                    continue
                value = getattr(self, attr)
                if np.size(value) == l:
                    setattr(self, attr, np.delete(value, index))

        def pk(self):
            return self[self.idx_max()]

        def add_chirp(self, chirp=0, order=1, s0=None):
            '''
            adds energy chirp of given "order" to the beam around "center" position [m]
            chirp = dE/E0/ds[um]**order = dg/g0/s[um]**order
            so that
            g_new = g0 + dg = g0 + g0 (s-s0)**order * chirp
            '''
            if chirp not in [None, 0]:
                if s0 is None:
                    s0 = (np.amax(self.s) - np.amin(self.s)) / 2
                E_center = self.E[find_nearest_idx(self.s, s0)]
                self.E += (self.s - s0) ** order * chirp * E_center * 1e6

        def add_chirp_poly(self, coeff, s0=None):
            '''
            The method adds a polynomial energy chirp to the beam object.

            coeff   --- coefficients for the chirp
            s0      --- the point with respect to which the chirp will be introduced

            The expression for the chirp:

            E = E0((g0 + coeff[0])/g0 +

                + coeff[1]*(s - s0))**1 / 1! / ((speed_of_light * 1e-15)**1  * g0) +

                + coeff[2]*(s - s0))**2 / 2! / ((speed_of_light * 1e-15)**2  * g0) +

                + coeff[3]*(s - s0))**3 / 3! / ((speed_of_light * 1e-15)**3  * g0) + ...

            ... + coeff[n]*(s - s0))**n / n! / ((speed_of_light * 1e-15)**n  * g0))

            where coeff[n] is represented in [1/fs**n]
            The convention for the coeff is introduced for convenient treatment this
            with respect to a radiation chirp in order to easily satisfy the resonant
            condition along the whole bunch in the case of linear electron bunch chirp.
            Here is the expression:

                2*dw/dt = (w0/g0) * dg/dt

            @author: Andrei Trebushinin

            '''
            _logger.debug('introducing a chirp to the ebeam')
            s = self.s

            if s0 is None:
                s0 = (np.amax(self.s) - np.amin(self.s)) / 2
            elif isinstance(s0, str) is not True:
                s0 = s0  # / 1e6
            else:
                raise ValueError("s0 must be None or some value")

            delta_s = s - s0
            E0 = self.E[find_nearest_idx(self.s, s0)]
            g0 = self.g[find_nearest_idx(self.s, s0)]
            #    coeff[0] += g0
            _logger.debug(ind_str + 'coeffs for chirp = {}'.format(coeff))
            coeff_norm = [ci / ((glb.speed_of_light * 1e-15) ** i * factorial(i) * g0) for i, ci in enumerate(coeff)]
            coeff_norm = list(np.flip(coeff_norm, axis=0))
            _logger.debug(ind_str + 'coeffs_norm = {}'.format(coeff_norm))
            coeff_norm = np.asarray(coeff_norm) * E0
            self.E += np.polyval(coeff_norm, delta_s)

        def add_wake(self, tube_radius=5e-3, tube_len=1, conductivity=3.66e+7, tau=7.1e-15, roughness=600e-9,
                     d_oxid=5e-9):
            self.eloss = pipe_wake(self.s, self.I, tube_radius, tube_len, conductivity, tau, roughness, d_oxid)[1][1][
                ::-1]

        def to_array(self, *args, **kwargs):
            raise NotImplementedError('Method inherited from Beam() class, not applicable for BeamArray objects')

