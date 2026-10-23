"""Synchrotron-radiation trajectory preparation and field integration.

The module tracks particles through a magnetic lattice, interpolates each
trajectory onto three-point Gauss-Legendre integration nodes, and accumulates
the complex radiation field on :class:`ocelot.rad.screen.Screen`.

Beamline coordinates enter in metres. ``traj2motion`` converts positions to
millimetres because the radiation integrator and internal screen coordinates
use millimetres. Photon energies are expressed in eV and electron energies in
GeV.
"""

from __future__ import annotations

import copy
import logging
import numbers
import sys
from collections.abc import Sequence
from math import pi

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.interpolate import splrep, splev

from ocelot.rad.spline_py import cspline_coef, moment_numba
from ocelot.rad.screen import Screen
from ocelot.common.globals import m_e_GeV, h_eV_s, q_e, speed_of_light, ro_e
from ocelot.cpbd.elements import Undulator
from ocelot.cpbd.field_map import field_map2field_func as _field_map2field_func
from ocelot.cpbd.high_order import rk_track_in_field, _rk_energy, _rk_slopes, _rk_momenta
from ocelot.cpbd import track
from ocelot.cpbd.navi import Navigator
from ocelot.cpbd import beam
import ocelot.cpbd.magnetic_lattice as mlattice

__author__ = 'Sergey Tomin'
_logger = logging.getLogger(__name__)

try:
    import numba as nb

    nb_flag = True
except ImportError:
    _logger.info("radiation_py.py: module NUMBA is not installed. Install it to speed up calculation")
    nb_flag = False


FloatArray = NDArray[np.float64]
# Parallel dispatch is reserved for sufficiently long spectra to avoid paying
# its compilation/startup cost for small one-point calculations.
_PARALLEL_SPECTRUM_MIN_ENERGIES = 64


def _screen_field_workspace(screen: Screen) -> Screen:
    """Create an empty field workspace with the same screen geometry.

    A shallow copy is sufficient because radiation calculations only mutate
    the field buffer and scalar geometry metadata on the workspace. This
    avoids copying potentially large result, motion, and trajectory arrays
    attached to a previously used screen.
    """
    workspace = copy.copy(screen)
    workspace.nullify()
    return workspace


def _trajectory_start_mm(
    trajectory: FloatArray,
    particle_index: int = 0,
) -> tuple[float, float, float]:
    """Return the first trajectory position in millimetres."""
    return (
        float(trajectory[0, particle_index] * 1000.0),
        float(trajectory[2, particle_index] * 1000.0),
        float(trajectory[4, particle_index] * 1000.0),
    )


class Motion:
    """Interpolated particle trajectory used by the radiation integrator.

    Position arrays ``x``, ``y``, and ``z`` are stored in millimetres.
    ``bx`` and ``by`` are trajectory slopes, ``Bx`` and ``By`` are magnetic
    field components in tesla, and ``XbetaI2``/``YbetaI2`` contain cumulative
    slope-squared path integrals.
    """

    def __init__(self):
        self.x = []
        self.y = []
        self.z = []
        self.bx = []
        self.by = []
        self.bz = []
        self.Bx = []
        self.By = []
        self.Bz = []
        self.XbetaI2 = []
        self.YbetaI2 = []


class BeamTraject:
    """Access particle coordinates stored in segmented trajectory arrays.

    Each trajectory segment has shape ``(9 * n_steps, n_particles)``. Rows in
    every nine-row block contain ``x, x', y, y', z, p, Bx, By, Bz``. Accessor
    methods concatenate the selected particle coordinate over all segments in
    beamline order.

    Parameters
    ----------
    beam_trajectories
        Sequence of trajectory arrays returned by :func:`track4rad_beam`.
    end_particle_array
        Optional tracked endpoint to retain for exact, repeatable restoration
        by :meth:`p_array_end`, including timing and reference energy.
    """

    def __init__(self, beam_trajectories: Sequence[FloatArray], *, end_particle_array=None):
        self.U = beam_trajectories
        self._end_particle_array = copy.deepcopy(end_particle_array)

    def n(self) -> int:
        """Return the number of particles in each trajectory segment."""
        return np.shape(self.U[0])[1]

    def check(self, n: int) -> None:
        """Validate an upper particle-index bound."""
        if n > self.n() - 1:
            raise Exception('n > number of particles')

    def _coordinate(self, row: int, n: int) -> FloatArray:
        self.check(n)
        return np.concatenate(
            [np.asarray(segment[row::9, n], dtype=float) for segment in self.U]
        )

    def x(self, n: int = 0) -> FloatArray:
        """Return horizontal positions in metres for particle ``n``."""
        return self._coordinate(0, n)

    def y(self, n: int = 0) -> FloatArray:
        """Return vertical positions in metres for particle ``n``."""
        return self._coordinate(2, n)

    def z(self, n: int = 0) -> FloatArray:
        """Return longitudinal Cartesian positions in metres for particle ``n``."""
        return self._coordinate(4, n)

    def xp(self, n: int = 0) -> FloatArray:
        """Return horizontal slopes ``dx/dz`` for particle ``n``."""
        return self._coordinate(1, n)

    def yp(self, n: int = 0) -> FloatArray:
        """Return vertical slopes ``dy/dz`` for particle ``n``."""
        return self._coordinate(3, n)

    def p(self, n: int = 0) -> FloatArray:
        """Return energy deviations (E-E0)/(p0*c) for particle ``n``."""
        return self._coordinate(5, n)

    def s(self, n: int = 0) -> FloatArray:
        """Return cumulative path length in metres for particle ``n``."""
        self.check(n)

        xp2 = self.xp(n) ** 2
        yp2 = self.yp(n) ** 2
        # zp = np.sqrt(1. / (1. + xp2 + yp2))
        s = cumulative_trapezoid(np.sqrt(1. + xp2 + yp2), self.z(n), initial=0)
        return s

    def p_array_end(self, p_array: beam.ParticleArray) -> None:
        """Restore the tracked endpoint, including its reference and arrival offsets.

        Radiation calculations retain the endpoint ParticleArray explicitly.
        For standalone legacy trajectories without that endpoint, the supplied
        array must contain the entrance state; conversion assumes a fixed
        Cartesian frame, constant reference energy and tau=ct-Z/beta0.
        """
        if self._end_particle_array is not None:
            endpoint = self._end_particle_array
            p_array.rparticles[:] = endpoint.rparticles
            p_array.E, p_array.s = endpoint.E, endpoint.s
            return

        tau = p_array.tau().copy()
        for trajectory in self.U:
            beta0, beta, _ = _rk_energy(p_array.E, trajectory[5])
            z = trajectory[4::9]
            path = trapezoid(np.sqrt(1 + trajectory[1::9]**2 + trajectory[3::9]**2), z, axis=0)
            tau += path / beta - (z[-1] - z[0]) / beta0
        final = _rk_momenta(self.U[-1][-9:-3], p_array.E)
        p_array.rparticles[:4] = final[:4]
        p_array.rparticles[4] = tau
        p_array.rparticles[5] = final[5]


def bspline(x: FloatArray, y: FloatArray, x_new: FloatArray) -> FloatArray:
    """Interpolate samples with an exact cubic B-spline.

    Parameters
    ----------
    x, y
        One-dimensional sample coordinates and values.
    x_new
        Coordinates at which to evaluate the spline.

    Returns
    -------
    numpy.ndarray
        Interpolated values at ``x_new``.
    """
    tck = splrep(x, y, s=0)
    return np.asarray(splev(x_new, tck, der=0))


def _integ_beta2_numpy(x: FloatArray, y: FloatArray) -> FloatArray:
    """Evaluate the spline-square integral with NumPy as a fallback."""
    a, b, c, d, _ = cspline_coef(x, y)
    h = np.diff(x)
    increments = h * (
        d * d
        + h
        * (
            c * d
            + h
            * (
                1. / 3. * (c * c + 2 * b * d)
                + h
                * (
                    0.5 * (b * c + a * d)
                    + h
                    * (
                        0.2 * (b * b + 2 * a * c)
                        + h * (1. / 3. * a * b + (a * a * h) / 7.)
                    )
                )
            )
        )
    )
    return np.concatenate(([0.0], np.cumsum(increments)))


def _integ_beta2_compiled_impl(x: FloatArray, y: FloatArray) -> FloatArray:
    """Compiled implementation of the legacy spline-square integral."""
    n = len(x)
    moments = moment_numba(x, y)
    result = np.zeros(n)
    cumulative = 0.0

    for i in range(n - 1):
        h = x[i + 1] - x[i]
        moment = moments[i]
        a = (moments[i + 1] - moment) / (6.0 * h)
        b = moment / 2.0
        c = (
            (y[i + 1] - y[i]) / h
            - moments[i + 1] * h / 6.0
            - moment * h / 3.0
        )
        d = y[i]
        cumulative += h * (
            d * d
            + h
            * (
                c * d
                + h
                * (
                    (c * c + 2.0 * b * d) / 3.0
                    + h
                    * (
                        (b * c + a * d) / 2.0
                        + h
                        * (
                            (b * b + 2.0 * a * c) / 5.0
                            + h * (a * b / 3.0 + a * a * h / 7.0)
                        )
                    )
                )
            )
        )
        result[i + 1] = cumulative

    return result


_integ_beta2_compiled = (
    nb.njit(cache=True)(_integ_beta2_compiled_impl)
    if nb_flag
    else None
)


def integ_beta2(x: FloatArray, y: FloatArray) -> FloatArray:
    """Integrate the square of a cubic-spline representation cumulatively.

    The spline coefficients and analytical interval integrals are evaluated
    by a Numba-compiled implementation when Numba is available. A numerically
    equivalent NumPy implementation is retained as the fallback.

    Parameters
    ----------
    x
        Strictly increasing spline coordinates.
    y
        Values to square and integrate.

    Returns
    -------
    numpy.ndarray
        Cumulative integral with ``result[0] == 0`` and the same length as
        ``x``.
    """
    x_array = np.asarray(x, dtype=np.float64)
    y_array = np.asarray(y, dtype=np.float64)
    if _integ_beta2_compiled is not None:
        return _integ_beta2_compiled(x_array, y_array)
    return _integ_beta2_numpy(x_array, y_array)


def x2xgaus(x: FloatArray) -> FloatArray:
    """Insert three-point Gauss-Legendre nodes into a uniform grid.

    The output contains the initial coordinate, three quadrature nodes for
    every interval, and the final coordinate. The input is assumed to contain
    at least two uniformly spaced coordinates; this preserves the historical
    behavior used by radiation trajectory generation.

    Parameters
    ----------
    x
        One-dimensional, uniformly spaced coordinates.

    Returns
    -------
    numpy.ndarray
        Coordinates with length ``3 * (len(x) - 1) + 2``.
    """
    sqrt35 = 0.5 * np.sqrt(3. / 5.)
    h = x[1] - x[0]
    gauss_offsets = np.array(
        [
            0.5 - sqrt35,
            0.5 - sqrt35 + sqrt35,
            0.5 - sqrt35 + sqrt35 + sqrt35,
        ]
    ) * h
    gauss_nodes = x[:-1, np.newaxis] + gauss_offsets
    return np.concatenate((x[:1], gauss_nodes.reshape(-1), x[-1:]))


def traj2motion(traj: FloatArray) -> Motion:
    """Convert a raw trajectory to interpolated radiation-integration data.

    Parameters
    ----------
    traj
        One particle trajectory as a flat array with repeated
        ``x, x', y, y', z, p, Bx, By, Bz`` blocks. Positions are in metres and
        magnetic-field components are in tesla.

    Returns
    -------
    Motion
        Interpolated trajectory at three-point Gauss-Legendre nodes. Position
        arrays are converted to millimetres.
    """
    motion = Motion()
    motion.x = traj[0::9]
    motion.y = traj[2::9]
    motion.z = traj[4::9]
    motion.bx = traj[1::9]
    motion.by = traj[3::9]
    motion.bz = traj[5::9]
    motion.Bx = traj[6::9]
    motion.By = traj[7::9]
    motion.Bz = traj[8::9]
    # new_motion = Motion()

    motion.z = motion.z.flatten()

    Z = x2xgaus(motion.z)

    motion.x = bspline(motion.z, motion.x, Z) * 1000.
    motion.y = bspline(motion.z, motion.y, Z) * 1000.
    # print "inegr = ", simps(motion.bx.flatten()**2, motion.z*1000)
    Ibx2 = integ_beta2(motion.z * 1000., motion.bx)
    Iby2 = integ_beta2(motion.z * 1000., motion.by)

    motion.XbetaI2 = bspline(motion.z, Ibx2, Z)
    motion.YbetaI2 = bspline(motion.z, Iby2, Z)

    motion.bx = bspline(motion.z, motion.bx, Z)
    motion.by = bspline(motion.z, motion.by, Z)

    motion.Bx = bspline(motion.z, motion.Bx, Z)
    motion.By = bspline(motion.z, motion.By, Z)
    motion.z = Z * 1000.

    return motion


def energy_loss_und(energy, Kx, lperiod, L, energy_loss=False):
    if energy_loss:
        k = 4. * pi * pi / 3. * ro_e / m_e_GeV
        U = k * energy ** 2 * Kx ** 2 * L / lperiod ** 2
    else:
        U = 0.
    return U


def sigma_gamma_quat(energy, Kx, lperiod, L):
    """
    rate of energy diffusion

    :param energy: electron beam energy
    :param Kx: undulator parameter
    :param lperiod: undulator period
    :param L: length
    :return: sigma_gamma/gamma
    """
    gamma = energy / m_e_GeV
    lambda_compt = 2.4263102389e-12  # m
    lambda_compt_r = lambda_compt / 2. / pi
    def f(K): return 1.2 + 1. / (K + 1.33 * K * K + 0.4 * K ** 3)
    delta_Eq2 = 56. * pi ** 3 / 15. * lambda_compt_r * ro_e * gamma ** 4 / lperiod ** 3 * Kx ** 3 * f(Kx) * L
    sigma_Eq = np.sqrt(delta_Eq2 / (gamma * gamma))
    return sigma_Eq


def quantum_diffusion(energy, Kx, lperiod, L, quantum_diff=False):
    if quantum_diff:
        # gamma = energy/m_e_GeV
        # lambda_compt = 2.4263102389e-12 # h_eV_s/m_e_eV*speed_of_light
        # lambda_compt_r = lambda_compt/2./pi
        # f = lambda K: 1.2 + 1./(K + 1.33*K*K + 0.4*K**3)
        # delta_Eq2 = 56.*pi**3/15.*lambda_compt_r*ro_e*gamma**4/lperiod**3*Kx**3*f(Kx)*L
        sigma_Eq = sigma_gamma_quat(energy, Kx, lperiod, L)  # sqrt(delta_Eq2/(gamma*gamma))
        U = sigma_Eq * np.random.randn() * energy
    else:
        U = 0.
    return U


def field_map2field_func(z, By):
    return _field_map2field_func(z, By)


def gintegrator(Xscr, Yscr, Erad, motion, screen, n, n_end, gamma, half_step):
    """

    :param Xscr:
    :param Yscr:
    :param Erad:
    :param motion:
    :param screen:
    :param n:
    :param n_end:
    :param gamma:
    :param half_step:
    :return:
    """
    Q = 0.5866740802042227  # speed_of_light/m_e_eV/1000  // e/mc = (mm*T)^-1
    hc = 1.239841874330e-3  # h_eV_s*speed_of_light*1000  // mm
    k2q3 = 1.1547005383792517  # ;//  = 2./sqrt(3)
    gamma2 = gamma * gamma
    w = [0.5555555555555556 * half_step, 0.8888888888888889 * half_step, 0.5555555555555556 * half_step]
    LenPntrConst = screen.Distance - motion.z[0]  # ; // I have to pay attention to this
    phaseConst = np.pi * Erad / (gamma2 * hc)
    for p in range(3):  # // Gauss integration
        i = n * 3 + p + 1
        # radConstAdd = w[p]*Q*k2q3*(screen.Distance - motion.z[0] - 0*screen.Zstart)
        XX = motion.x[i]
        YY = motion.y[i]
        ZZ = motion.z[i]
        BetX = motion.bx[i]
        BetY = motion.by[i]
        IbetX2 = motion.XbetaI2[i]
        IbetY2 = motion.YbetaI2[i]
        Bx = motion.Bx[i]
        By = motion.By[i]
        LenPntrZ = screen.Distance - ZZ

        prX = Xscr - XX  # //for pointer nx(z)
        prY = Yscr - YY  # //for pointer ny(z)
        nx = prX / LenPntrZ
        ny = prY / LenPntrZ
        tx = gamma * (nx - BetX)
        ty = gamma * (ny - BetY)
        tx2 = tx * tx
        ty2 = ty * ty
        tyx = 2. * tx * ty

        radConst = w[p] * Q * k2q3 * (screen.Distance) / LenPntrZ / ((1. + tx2 + ty2) * (1. + tx2 + ty2))

        radX = radConst * (By * (1. - tx2 + ty2) + Bx * tyx - 2. * tx / Q / LenPntrZ)  # /*sigma*/
        radY = -radConst * (Bx * (1. + tx2 - ty2) + By * tyx + 2. * ty / Q / LenPntrZ)  # ;/*pi*/

        prXconst = Xscr - motion.x[0]
        prYconst = Yscr - motion.y[0]
        phaseConstIn = (prXconst * prXconst + prYconst * prYconst) / LenPntrConst
        phaseConstCur = (prX * prX + prY * prY) / LenPntrZ
        # // string below is for case direct accumulation
        # //double phase = screen->Phase[ypoint*xpoint*je + xpoint*jy + jx] + faseConst*(ZZ - motion->Z[0]  + gamma2*(IbetX2 + IbetY2 + phaseConstCur - phaseConstIn));

        phase = phaseConst * (
            ZZ - motion.z[0] + gamma2 * (IbetX2 + IbetY2 + phaseConstCur - phaseConstIn)) + screen.arPhase

        cosf = np.cos(phase)
        sinf = np.sin(phase)
        EreX = radX * cosf  # //(cosf *cos(fase0) - sinf*sin(fase0));
        EimX = radX * sinf  # //(sinf *cos(fase0) + cosf*sin(fase0));
        EreY = radY * cosf
        EimY = radY * sinf

        screen.arReEx += EreX
        screen.arImEx += EimX
        screen.arReEy += EreY
        screen.arImEy += EimY
        if i == n_end:  # //(n == 5000 && p == 2)
            LenPntrZ = screen.Distance - motion.z[-1]
            prX = Xscr - motion.x[-1]  # //for pointer nx(z)
            prY = Yscr - motion.y[-1]  # //for pointer ny(z)
            IbetX2 = motion.XbetaI2[-1]
            IbetY2 = motion.YbetaI2[-1]
            phase = phaseConst * (motion.z[-1] - motion.z[0] + gamma2 * (
                IbetX2 + IbetY2 + prX * prX / LenPntrZ + prY * prY / LenPntrZ - phaseConstIn))
            screen.arPhase = screen.arPhase + phase
    return screen


def gintegrator_over_traj_py(Nmotion, Xscr, Yscr, Erad, n_end, gamma, half_step, Distance, x, y, z, bx, by,
                             XbetaI2, YbetaI2, Bx, By, arReEx, arImEx, arReEy, arImEy, arPhase):
    """Accumulate radiation fields over one interpolated trajectory segment.

    This is the numerical hot loop compiled by Numba when available. Screen
    coordinate arrays and output field arrays may be one- or two-dimensional,
    depending on whether the caller requests a spectrum, line, or spatial
    screen. All position-like inputs use millimetres.

    The four field arrays and ``arPhase`` are updated in place.
    """

    q = 0.5866740802042227  # speed_of_light/m_e_eV/1000  // e/mc = (mm*T)^-1
    hc = 1.239841874330e-3  # h_eV_s*speed_of_light*1000  // mm
    k2q3 = 1.1547005383792517  # 2./sqrt(3)
    gamma2 = gamma * gamma
    w = np.array([0.5555555555555556 * half_step, 0.8888888888888889 * half_step, 0.5555555555555556 * half_step])
    LenPntrConst = Distance - z[0]  # I have to pay attention to this
    phaseConst = np.pi * Erad / (gamma2 * hc)
    prXconst = Xscr - x[0]
    prYconst = Yscr - y[0]
    phaseConstIn = (prXconst * prXconst + prYconst * prYconst) / LenPntrConst

    for n in range(Nmotion - 1):
        for p in range(3):  # // Gauss integration
            i = n * 3 + p + 1
            # radConstAdd = w[p]*Q*k2q3*(screen.Distance - motion.z[0] - 0*screen.Zstart)

            LenPntrZ = Distance - z[i]

            prX = Xscr - x[i]  # //for pointer nx(z)
            prY = Yscr - y[i]  # //for pointer ny(z)
            nx = prX / LenPntrZ
            ny = prY / LenPntrZ
            tx = gamma * (nx - bx[i])
            ty = gamma * (ny - by[i])
            tx2 = tx * tx
            ty2 = ty * ty
            tyx = 2. * tx * ty

            denominator = 1. + tx2 + ty2
            radConst = w[p] * q * k2q3 * Distance / LenPntrZ / (denominator * denominator)

            radX = radConst * (By[i] * (1. - tx2 + ty2) + Bx[i] * tyx - 2. * tx / q / LenPntrZ)  # /*sigma*/
            radY = -radConst * (Bx[i] * (1. + tx2 - ty2) + By[i] * tyx + 2. * ty / q / LenPntrZ)  # ;/*pi*/

            phaseConstCur = (prX * prX + prY * prY) / LenPntrZ
            # // string below is for case direct accumulation
            # //double phase = screen->Phase[ypoint*xpoint*je + xpoint*jy + jx] + faseConst*(ZZ - motion->Z[0]  + gamma2*(IbetX2 + IbetY2 + phaseConstCur - phaseConstIn));

            # here the constant phase shift was subtracted
            phase = phaseConst * (
                z[i] - z[0] + gamma2 * (XbetaI2[i] + YbetaI2[i] + phaseConstCur - phaseConstIn)) + arPhase

            # phase = phaseConst * (z[i] - z[0] + gamma2 * (XbetaI2[i] + YbetaI2[i] + phaseConstCur)) + arPhase # + (LenPntrConst *2*np.pi * Erad/hc)%(2*np.pi)
            cosf = np.cos(phase)
            sinf = np.sin(phase)
            EreX = radX * cosf  # (cosf *cos(fase0) - sinf*sin(fase0));
            EimX = radX * sinf  # (sinf *cos(fase0) + cosf*sin(fase0));
            EreY = radY * cosf
            EimY = radY * sinf

            arReEx += EreX
            arImEx += EimX
            arReEy += EreY
            arImEy += EimY
            if i == n_end:
                LenPntrZ = Distance - z[-1]
                prX = Xscr - x[-1]  # //for pointer nx(z)
                prY = Yscr - y[-1]  # //for pointer ny(z)
                phase = phaseConst * (z[-1] - z[0] + gamma2 * (
                    XbetaI2[-1] + YbetaI2[-1] + prX * prX / LenPntrZ + prY * prY / LenPntrZ - phaseConstIn))
                # phase = phaseConst * (z[-1] - z[0] + gamma2 * (
                #        XbetaI2[-1] + YbetaI2[-1] + prX * prX / LenPntrZ + prY * prY / LenPntrZ )) # + (LenPntrConst *2*np.pi * Erad/hc)%(2*np.pi)
                arPhase += phase


gintegrator_over_traj = gintegrator_over_traj_py if not nb_flag else nb.jit(gintegrator_over_traj_py, nopython=True)


def gintegrator_over_spectrum_py(
    Nmotion,
    Xscr,
    Yscr,
    Erad,
    n_end,
    gamma,
    half_step,
    Distance,
    x,
    y,
    z,
    bx,
    by,
    XbetaI2,
    YbetaI2,
    Bx,
    By,
    arReEx,
    arImEx,
    arReEy,
    arImEy,
    arPhase,
):
    """Parallel field integration for an on-axis or off-axis spectrum.

    The screen contains one transverse point and multiple independent photon
    energies. Numba parallelizes the outer energy loop; each iteration writes
    to a distinct element of the field and phase arrays.
    """
    q = 0.5866740802042227
    hc = 1.239841874330e-3
    k2q3 = 1.1547005383792517
    gamma2 = gamma * gamma
    weights = np.array(
        [
            0.5555555555555556 * half_step,
            0.8888888888888889 * half_step,
            0.5555555555555556 * half_step,
        ]
    )
    len_pointer_const = Distance - z[0]
    pr_x_const = Xscr - x[0]
    pr_y_const = Yscr - y[0]
    phase_const_in = (
        pr_x_const * pr_x_const + pr_y_const * pr_y_const
    ) / len_pointer_const

    for energy_index in nb.prange(len(Erad)):
        phase_const = np.pi * Erad[energy_index] / (gamma2 * hc)
        re_ex = arReEx[energy_index]
        im_ex = arImEx[energy_index]
        re_ey = arReEy[energy_index]
        im_ey = arImEy[energy_index]
        initial_phase = arPhase[energy_index]
        final_phase = 0.0

        for n in range(Nmotion - 1):
            for quadrature_index in range(3):
                i = n * 3 + quadrature_index + 1
                len_pointer_z = Distance - z[i]
                pr_x = Xscr - x[i]
                pr_y = Yscr - y[i]
                tx = gamma * (pr_x / len_pointer_z - bx[i])
                ty = gamma * (pr_y / len_pointer_z - by[i])
                tx2 = tx * tx
                ty2 = ty * ty
                tyx = 2.0 * tx * ty

                denominator = 1.0 + tx2 + ty2
                radiation_const = (
                    weights[quadrature_index]
                    * q
                    * k2q3
                    * Distance
                    / len_pointer_z
                    / (denominator * denominator)
                )
                radiation_x = radiation_const * (
                    By[i] * (1.0 - tx2 + ty2)
                    + Bx[i] * tyx
                    - 2.0 * tx / q / len_pointer_z
                )
                radiation_y = -radiation_const * (
                    Bx[i] * (1.0 + tx2 - ty2)
                    + By[i] * tyx
                    + 2.0 * ty / q / len_pointer_z
                )

                phase_const_current = (
                    pr_x * pr_x + pr_y * pr_y
                ) / len_pointer_z
                phase = phase_const * (
                    z[i]
                    - z[0]
                    + gamma2
                    * (
                        XbetaI2[i]
                        + YbetaI2[i]
                        + phase_const_current
                        - phase_const_in
                    )
                ) + initial_phase
                cos_phase = np.cos(phase)
                sin_phase = np.sin(phase)
                re_ex += radiation_x * cos_phase
                im_ex += radiation_x * sin_phase
                re_ey += radiation_y * cos_phase
                im_ey += radiation_y * sin_phase

                if i == n_end:
                    len_pointer_end = Distance - z[-1]
                    pr_x_end = Xscr - x[-1]
                    pr_y_end = Yscr - y[-1]
                    final_phase = phase_const * (
                        z[-1]
                        - z[0]
                        + gamma2
                        * (
                            XbetaI2[-1]
                            + YbetaI2[-1]
                            + pr_x_end * pr_x_end / len_pointer_end
                            + pr_y_end * pr_y_end / len_pointer_end
                            - phase_const_in
                        )
                    )

        arReEx[energy_index] = re_ex
        arImEx[energy_index] = im_ex
        arReEy[energy_index] = re_ey
        arImEy[energy_index] = im_ey
        arPhase[energy_index] = initial_phase + final_phase


gintegrator_over_spectrum = (
    nb.njit(parallel=True, cache=True)(gintegrator_over_spectrum_py)
    if nb_flag
    else None
)


def wrap_gintegrator(Nmotion, Xscr, Yscr, Erad, motion, screen, n_end, gamma, half_step):
    Distance = screen.Distance
    x, y, z = motion.x, motion.y, motion.z
    bx, by = motion.bx, motion.by
    XbetaI2, YbetaI2 = motion.XbetaI2, motion.YbetaI2
    Bx, By = motion.Bx, motion.By

    gintegrator_over_traj(Nmotion, Xscr, Yscr, Erad, n_end, gamma, half_step, Distance,
                          x, y, z, bx, by, XbetaI2,
                          YbetaI2, Bx, By,
                          screen.arReEx, screen.arImEx, screen.arReEy, screen.arImEy, screen.arPhase)

    return screen


def wrap_gintegrator_spectrum(
    Nmotion,
    Xscr,
    Yscr,
    Erad,
    motion,
    screen,
    n_end,
    gamma,
    half_step,
):
    """Run the parallel one-point spectrum integration kernel."""
    gintegrator_over_spectrum(
        Nmotion,
        Xscr,
        Yscr,
        Erad,
        n_end,
        gamma,
        half_step,
        screen.Distance,
        motion.x,
        motion.y,
        motion.z,
        motion.bx,
        motion.by,
        motion.XbetaI2,
        motion.YbetaI2,
        motion.Bx,
        motion.By,
        screen.arReEx,
        screen.arImEx,
        screen.arReEy,
        screen.arImEy,
        screen.arPhase,
    )
    return screen


def radiation_py(gamma: float, traj: FloatArray, screen: Screen) -> int:
    """Accumulate one trajectory segment's radiation field on a screen.

    Parameters
    ----------
    gamma
        Relativistic Lorentz factor of the particle in this trajectory
        segment.
    traj
        Flat trajectory array containing repeated
        ``x, x', y, y', z, p, Bx, By, Bz`` blocks.
    screen
        :class:`ocelot.rad.screen.Screen` whose complex field components and
        phase are updated in place.

    Notes
    -----
    The function accumulates rather than replaces existing screen fields.
    Callers that need a segment-local field must clear ``arReEx``, ``arImEx``,
    ``arReEy``, and ``arImEy`` before calling. ``arPhase`` is intentionally
    propagated between consecutive trajectory segments.

    Flattened screen fields use logical order ``(energy, y, x)``.
    One-point spectra with at least 64 energy samples use a dedicated Numba
    kernel that parallelizes independent photon energies.

    Returns
    -------
    int
        The legacy success value ``1``. Radiation data are returned through
        the mutated ``screen`` object.
    """

    motion = traj2motion(traj)

    size = len(motion.z)
    Nmotion = int((size + 1) / 3)
    half_step = (motion.z[-1] - motion.z[0]) / 2. / (Nmotion - 1)

    n_end = len(motion.z) - 2
    Xscr = np.linspace(screen.x_start, screen.x_start + screen.x_step * (screen.nx - 1), num=screen.nx)
    Yscr = np.linspace(screen.y_start, screen.y_start + screen.y_step * (screen.ny - 1), num=screen.ny)
    Yscr = Yscr.reshape((screen.ny, 1))
    Erad = np.linspace(screen.e_start, screen.e_start + screen.e_step * (screen.ne - 1), num=screen.ne)

    shape_array = [screen.ne, screen.ny, screen.nx]
    use_parallel_spectrum = (
        gintegrator_over_spectrum is not None
        and screen.nx == 1
        and screen.ny == 1
        and screen.ne >= _PARALLEL_SPECTRUM_MIN_ENERGIES
    )
    if use_parallel_spectrum:
        screen.arReEx = screen.arReEx.reshape(-1)
        screen.arImEx = screen.arImEx.reshape(-1)
        screen.arReEy = screen.arReEy.reshape(-1)
        screen.arImEy = screen.arImEy.reshape(-1)
        screen.arPhase = screen.arPhase.reshape(-1)
        wrap_gintegrator_spectrum(
            Nmotion,
            Xscr[0],
            Yscr[0, 0],
            Erad,
            motion,
            screen,
            n_end,
            gamma,
            half_step,
        )
    elif 1 in shape_array:
        Erad = Erad.reshape((screen.ne, 1))
        if screen.ny > 1 and screen.ne > 1:
            Yscr = Yscr.reshape((1, screen.ny))
        shape_array.remove(1)
        screen.arReEx = screen.arReEx.reshape(shape_array)
        screen.arImEx = screen.arImEx.reshape(shape_array)
        screen.arReEy = screen.arReEy.reshape(shape_array)
        screen.arImEy = screen.arImEy.reshape(shape_array)
        screen.arPhase = screen.arPhase.reshape(shape_array)

        wrap_gintegrator(Nmotion, Xscr, Yscr, Erad, motion, screen, n_end, gamma, half_step)
        # print("phase", screen.arPhase )
        # for n in range(Nmotion-1):
        #     screen = gintegrator(Xscr, Yscr, Erad, motion, screen, n, n_end, gamma, half_step)
        screen.arReEx = screen.arReEx.reshape(-1)
        screen.arImEx = screen.arImEx.reshape(-1)
        screen.arReEy = screen.arReEy.reshape(-1)
        screen.arImEy = screen.arImEy.reshape(-1)
        screen.arPhase = screen.arPhase.reshape(-1)
    else:
        Erad = Erad.reshape((screen.ne, 1))
        print("SR 3D calculation")
        arReEx = np.empty_like(screen.arReEx)
        arImEx = np.empty_like(screen.arImEx)
        arReEy = np.empty_like(screen.arReEy)
        arImEy = np.empty_like(screen.arImEy)
        screen_segment = copy.copy(screen)

        n_pl = screen.ny * screen.nx
        screen_segment.arReEx = np.zeros((screen.ny, screen.nx))
        screen_segment.arImEx = np.zeros((screen.ny, screen.nx))
        screen_segment.arReEy = np.zeros((screen.ny, screen.nx))
        screen_segment.arImEy = np.zeros((screen.ny, screen.nx))

        for i, erad in enumerate(Erad):
            start = i * n_pl
            stop = start + n_pl
            arPhase = screen.arPhase[start:stop]
            screen_segment.arReEx.fill(0.0)
            screen_segment.arImEx.fill(0.0)
            screen_segment.arReEy.fill(0.0)
            screen_segment.arImEy.fill(0.0)

            screen_segment.arPhase = arPhase.reshape((screen.ny, screen.nx))

            # for n in range(Nmotion-1):
            #    screen_segment = gintegrator(Xscr, Yscr, erad, motion, screen_segment, n, n_end, gamma, half_step)
            wrap_gintegrator(Nmotion, Xscr, Yscr, erad, motion, screen_segment, n_end, gamma, half_step)

            arReEx[start:stop] = screen_segment.arReEx.reshape(-1)
            arImEx[start:stop] = screen_segment.arImEx.reshape(-1)
            arReEy[start:stop] = screen_segment.arReEy.reshape(-1)
            arImEy[start:stop] = screen_segment.arImEy.reshape(-1)
            screen.arPhase[start:stop] = screen_segment.arPhase.reshape(-1)
        screen.arReEx[:] += arReEx[:]
        screen.arImEx[:] += arImEx[:]
        screen.arReEy[:] += arReEy[:]
        screen.arImEy[:] += arImEy[:]
    return 1


def calculate_radiation(
    lat,
    screen,
    ebeam,
    energy_loss=False,
    quantum_diff=False,
    accuracy=1,
    **kwargs,
):
    """Calculate incoherent synchrotron radiation from a beam description.

    Parameters
    ----------
    lat
        Magnetic lattice to track.
    screen
        Observation screen, mutated and returned.
    ebeam
        :class:`~ocelot.cpbd.beam.Beam` supplying electron coordinates,
        reference energy in GeV, and current in amperes.
    energy_loss
        Apply one aggregate classical energy correction per undulator.
    quantum_diff
        Apply one stochastic energy correction per undulator.
    accuracy
        Scale the automatically estimated trajectory-point count. Explicit
        undulator ``npoints`` values override this scale.
    **kwargs
        Compatibility arguments. ``end_poles`` is obsolete and must instead
        be configured on each :class:`Undulator`.

    Returns
    -------
    Screen
        The supplied screen containing photon flux in photons per second, per
        square millimetre, per ``10**-3`` relative bandwidth.

    Raises
    ------
    TypeError
        If ``ebeam`` is not a :class:`~ocelot.cpbd.beam.Beam`.
    ValueError
        If the lattice has no trackable nonzero-length elements.
    """
    if "end_poles" in kwargs:
        _logger.warning("The argument 'end_poles' is obsolete. It has been moved to the Undulator element.")
    screen.update()

    if isinstance(ebeam, beam.Beam):
        p = beam.Particle(x=ebeam.x, y=ebeam.y, px=ebeam.xp, py=ebeam.yp, E=ebeam.E)
        p_array = beam.ParticleArray()
        p_array.list2array([p])

    # elif beam.__class__ is beam.ParticleArray:
    #    b_current = beam.q_array[0] * 1000.
    #    p_array = beam

    else:
        raise TypeError("'beam' object must be Beam class")

    if ebeam.I == 0:
        print("Beam charge or beam current is 0. Default current I=0.1 A is used")
        ebeam.I = 0.1  # A

    p_array.tau()[:] = 0

    U, E = track4rad_beam(p_array, lat, energy_loss=energy_loss, quantum_diff=quantum_diff, accuracy=accuracy)
    # plt.plot(U[0][4::9, :], U[0][::9, :])
    # plt.show()
    screen_copy = _screen_field_workspace(screen)
    gamma_normalization = ebeam.E / m_e_GeV
    for i in range(p_array.n):
        for u, e in zip(U, E):
            # Each segment's kernel returns fields normalized by its gamma.
            # Accumulate in one common normalization when energy loss is enabled.
            for component in (screen_copy.arReEx, screen_copy.arImEx,
                              screen_copy.arReEy, screen_copy.arImEy):
                component.fill(0.)
            gamma = (e + u[5, i] * np.sqrt(e**2 - m_e_GeV**2)) / m_e_GeV
            radiation_py(gamma, u[:, i], screen_copy)
            weight = gamma / gamma_normalization
            screen.arReEx += weight * screen_copy.arReEx
            screen.arImEx += weight * screen_copy.arImEx
            screen.arReEy += weight * screen_copy.arReEy
            screen.arImEy += weight * screen_copy.arImEy
        screen.arPhase += screen_copy.arPhase
    screen.distPhoton(gamma_normalization, current=ebeam.I)
    screen.Ef_electron = p_array.E
    screen.motion = U
    beam_traj = BeamTraject(beam_trajectories=U, end_particle_array=p_array)
    screen.beam_traj = beam_traj

    # adding fast oscillating term to the phase
    screen.rebuild_efields(*_trajectory_start_mm(U[0]))

    return screen


def coherent_radiation(
    lat: mlattice.MagneticLattice,
    screen: Screen,
    p_array: beam.ParticleArray,
    energy_loss: bool = False,
    quantum_diff: bool = False,
    accuracy: float = 1,
    verbose: bool = True,
    **kwargs: object,
) -> Screen:
    """Calculate coherently summed radiation from a macroparticle ensemble.

    Each particle is tracked through ``lat``. Its initial longitudinal
    coordinate supplies the radiation phase, and each trajectory segment is
    weighted by the segment Lorentz factor and the number of electrons
    represented by that macroparticle. Complex fields are summed before the
    photon distribution is calculated.

    Parameters
    ----------
    lat
        Magnetic lattice containing the radiating elements.
    screen
        Observation screen. Coordinates are configured in metres and photon
        energies in eV. Field arrays and photon distributions are replaced in
        place.
    p_array
        Macroparticle ensemble. Reference energy is in GeV, charge is in
        coulombs, and longitudinal coordinates are in metres.
        :class:`ParticleArray` subclasses are accepted.
    energy_loss
        Include the aggregate classical energy loss for each undulator
        element.
    quantum_diff
        Apply a stochastic energy kick for each undulator element.
    accuracy
        Scale factor applied to the default number of trajectory samples.
    verbose
        Write particle progress to standard output.
    **kwargs
        Compatibility arguments. ``end_poles`` is obsolete and must instead
        be configured on each :class:`Undulator`.

    Returns
    -------
    Screen
        The supplied screen after coherent fields and photon distributions
        have been calculated.

    Raises
    ------
    TypeError
        If ``p_array`` is not a :class:`ParticleArray` instance or subclass.
    ValueError
        If the lattice has no trackable nonzero-length elements.

    Notes
    -----
    ``screen`` and ``p_array`` are mutated. The initial particle ``tau`` values
    are used for radiation phase, while final particle coordinates are written
    directly by tracking, retaining their entrance arrival offsets. The endpoint
    is also saved by :class:`BeamTraject`.

    The photon distributions ``screen.Total``, ``screen.Sigma``, and
    ``screen.Pi`` are normalized per passage of the supplied bunch, per
    square millimetre, and per ``10**-3`` relative bandwidth. No bunch
    repetition rate is an input to this function. The shared
    :func:`ocelot.gui.sr_plot.show_flux` plot label says ``ph/sec``; for this
    coherent calculation that label is numerically equivalent to assuming one
    bunch per second. For a machine repetition rate ``f_rep`` in hertz,
    multiply the photon distributions by ``f_rep`` to obtain photons per
    second.

    ``10**-3 BW`` means a relative photon-energy bandwidth
    ``delta_E / E = 10**-3`` (equivalently
    ``abs(delta_lambda) / lambda = 10**-3``). For a one-point spectrum, the
    photon density integrated over an energy interval is therefore

    ``np.trapezoid(screen.Total / (1e-3 * screen.Eph), screen.Eph)``

    after selecting the required energy range. The result remains per bunch
    and per square millimetre. With equally spaced samples, the corresponding
    bin sum is approximately
    ``sum(Total[i] * delta_E / (1e-3 * Eph[i]))``. A transverse integration is
    additionally required to obtain photons per bunch rather than photon
    density at one observation point.

    Segment field arrays are cleared between calls to :func:`radiation_py`,
    but ``arPhase`` is retained so consecutive trajectory segments interfere
    coherently without earlier segment fields being counted again.
    """
    if "end_poles" in kwargs:
        _logger.warning("The argument 'end_poles' is obsolete. It has been moved to the Undulator element.")
    screen.update()

    if not isinstance(p_array, beam.ParticleArray):
        raise TypeError("'p_array' must be a ParticleArray instance")

    tau0 = np.copy(p_array.tau())
    p_array.tau()[:] = 0

    U, E = track4rad_beam(p_array, lat, energy_loss=energy_loss, quantum_diff=quantum_diff, accuracy=accuracy)
    p_array.tau()[:] += tau0
    # plt.plot(U[0][4::9, :], U[0][::9, :])
    # plt.show()
    screen_copy = _screen_field_workspace(screen)
    wavelengths = h_eV_s * speed_of_light / screen.Eph
    electron_counts = p_array.q_array / q_e
    for i in range(p_array.n):
        # print("%i/%i" % (i, p_array.n))
        screen_copy.arPhase[:] = tau0[i] / wavelengths * 2 * np.pi
        # Number of electrons represented by this macroparticle.
        n_e = electron_counts[i]

        for u, e in zip(U, E):
            # radiation_py() accumulates fields in the supplied screen. Reset
            # the field components for each trajectory segment so that a
            # previously calculated segment is not added to the output again.
            # arPhase must be preserved to retain the phase relation between
            # consecutive segments.
            screen_copy.arReEx.fill(0.0)
            screen_copy.arImEx.fill(0.0)
            screen_copy.arReEy.fill(0.0)
            screen_copy.arImEy.fill(0.0)

            gamma = (e + u[5, i] * np.sqrt(e**2 - m_e_GeV**2)) / m_e_GeV

            radiation_py(gamma, u[:, i], screen_copy)

            screen.arReEx += screen_copy.arReEx * n_e * gamma
            screen.arImEx += screen_copy.arImEx * n_e * gamma
            screen.arReEy += screen_copy.arReEy * n_e * gamma
            screen.arImEy += screen_copy.arImEy * n_e * gamma
        if verbose:
            sys.stdout.write("\r" + "n: " + str(i) + " / " + str(p_array.n - 1))
            sys.stdout.flush()
    screen.coherent_photon_dist()
    screen.rebuild_efields(*_trajectory_start_mm(U[0]))

    screen.Ef_electron = p_array.E
    screen.motion = U
    beam_traj = BeamTraject(beam_trajectories=U, end_particle_array=p_array)
    screen.beam_traj = beam_traj

    return screen


def _undulator_trajectory_points(elem: Undulator, accuracy: float) -> int:
    """Return the number of Runge-Kutta trajectory points for radiation.

    An explicit ``Undulator(..., npoints=N)`` uses exactly ``N`` points,
    matching the semantics of Ocelot's Runge-Kutta transformations. Otherwise
    the historical length-based estimate is scaled by ``accuracy``.
    """
    # Normal lattice RK tracking gets npoints from the active transformation.
    # Radiation tracking calls rk_track_in_field() directly, so it must recover
    # the same configuration itself. Prefer an active RK-style map because it
    # includes set_tm(..., npoints=...) overrides, then fall back to constructor
    # transformation configuration stored by OpticElement. This private lookup
    # is intentionally isolated here because the radiation path bypasses the
    # transformation object that normally consumes these parameters.
    npoints = next(
        (
            tm.npoints
            for tm in elem.tms
            if getattr(tm, "npoints", None) is not None
        ),
        None,
    )
    if npoints is None:
        npoints = getattr(elem, "_kwargs", {}).get("npoints")
    if npoints is None:
        return int((elem.l * 1500 + 100) * accuracy)
    if isinstance(npoints, bool) or not isinstance(npoints, numbers.Integral):
        raise TypeError("Undulator npoints must be an integer")
    if npoints < 4:
        raise ValueError(
            "Undulator npoints must be at least 4 for cubic trajectory interpolation"
        )
    return int(npoints)


def _track_non_undulator_segment(
    p_array: beam.ParticleArray,
    elements: Sequence[object],
    z_start: float,
    energy: float,
    accuracy: float,
) -> tuple[FloatArray | None, float]:
    """Track a buffered non-undulator section and sample its trajectory.

    Samples include both endpoints, with slopes converted from mechanical
    momenta. This path assumes a straight Cartesian layout; general placement
    of curved local-map samples needs the corresponding reference geometry.
    """
    if not elements:
        return None, z_start

    section = mlattice.MagneticLattice(elements)
    if section.totalLen == 0:
        return None, z_start

    navigator = Navigator(section)
    n_points = int((section.totalLen * 2000 + 150) * accuracy)
    trajectory = np.zeros(
        (n_points * 9, np.shape(p_array.rparticles)[1])
    )
    step = section.totalLen / (n_points - 1)

    for index, z in enumerate(
        np.linspace(z_start, z_start + section.totalLen, num=n_points)
    ):
        if index:
            track.tracking_step(section, p_array, step, navigator)
        slopes = _rk_slopes(p_array.rparticles, p_array.E)
        trajectory[index * 9 + 0, :] = slopes[0]
        trajectory[index * 9 + 1, :] = slopes[1]
        trajectory[index * 9 + 2, :] = slopes[2]
        trajectory[index * 9 + 3, :] = slopes[3]
        trajectory[index * 9 + 4, :] = z
        trajectory[index * 9 + 5, :] = p_array.rparticles[5]

    return trajectory, z_start + section.totalLen


def track4rad_beam(
    p_array: beam.ParticleArray,
    lat: mlattice.MagneticLattice,
    energy_loss: bool = False,
    quantum_diff: bool = False,
    accuracy: float = 1,
    **kwargs: object,
) -> tuple[list[FloatArray], list[float]]:
    """Track a particle array and collect radiation trajectory segments.

    Parameters
    ----------
    p_array
        Particle ensemble to track. Coordinates are mutated to the end of each
        processed segment, using normalized mechanical momenta and retaining
        entrance tau. Undulator segments use tau=ct-Z/beta0 in the fixed frame.
        General curved laboratory assembly is not implemented by this routine.
    lat
        Magnetic lattice containing the radiating elements. Undulator
        subclasses are recognized as radiating elements.
    energy_loss
        Apply one aggregate classical energy correction per undulator.
    quantum_diff
        Apply one stochastic energy correction per undulator.
    accuracy
        Scale the automatically estimated trajectory-point count. For an
        undulator with an explicit ``npoints`` transformation parameter,
        ``npoints`` is used exactly and overrides this scale factor.
    **kwargs
        Compatibility arguments. ``end_poles`` is obsolete and must instead
        be configured on each :class:`Undulator`.

    Returns
    -------
    trajectories
        List of arrays with shape ``(9 * n_points, n_particles)``. Every
        nine-row block stores ``x, x', y, y', z, relative_momentum,
        Bx, By, Bz``. Consecutive non-undulator elements form one segment;
        this includes a non-undulator section after the final undulator.
    energies
        Reference energy in GeV for each trajectory segment.

    Raises
    ------
    TypeError
        If an explicit undulator ``npoints`` value is not an integer.
    ValueError
        If explicit ``npoints`` is less than four, which is insufficient for
        the cubic interpolation used by :func:`traj2motion`, or if the lattice
        has no trackable nonzero-length elements.
    """
    if "end_poles" in kwargs:
        _logger.warning("The argument 'end_poles' is obsolete. It has been moved to the Undulator element.")
    energy = p_array.E
    # Y0 = [beam.x, beam.xp, beam.y, beam.yp, 0, 0]
    # p = Particle(x=beam.x, px=beam.xp, y=beam.yp, py=beam.yp, E=beam.E)
    L = 0.
    U = []
    E = []
    non_u = []
    for elem in lat.sequence:
        if elem.l == 0:
            continue
        if not isinstance(elem, Undulator):
            non_u.append(elem)
            continue

        u, L = _track_non_undulator_segment(
            p_array,
            non_u,
            L,
            energy,
            accuracy,
        )
        if u is not None:
            U.append(u)
            E.append(energy)
        non_u = []

        energy_correction = energy_loss_und(
            energy,
            elem.Kx,
            elem.lperiod,
            elem.l,
            energy_loss,
        )
        energy_correction += quantum_diffusion(
            energy,
            elem.Kx,
            elem.lperiod,
            elem.l,
            quantum_diff,
        )

        mag_length = elem.l
        mag_field = elem.element.create_runge_kutta_main_params(energy).mag_field

        N = _undulator_trajectory_points(elem, accuracy)
        beta0, beta, _ = _rk_energy(energy, p_array.p())
        u, paths = rk_track_in_field(
            _rk_slopes(p_array.rparticles, energy),
            mag_length,
            N,
            energy,
            mag_field,
            return_path=True,
        )

        p_array.rparticles[:4] = _rk_momenta(u[-9:-3], energy)[:4]
        p_array.tau()[:] += paths[-1] / beta - mag_length / beta0
        p_array.s += mag_length
        s = u[-5, 0]
        u[4::9] += L
        L += s
        U.append(u)
        E.append(energy)
        if energy_correction != 0:
            # The existing radiation model removes the same energy from each
            # particle along its direction, once at the undulator exit.
            p0_before = np.sqrt(energy**2 - m_e_GeV**2)
            particle_energy = energy + p_array.p() * p0_before
            momentum_before = np.sqrt(particle_energy**2 - m_e_GeV**2)
            new_energy = energy - energy_correction
            if new_energy <= m_e_GeV or np.any(particle_energy - energy_correction <= m_e_GeV):
                raise ValueError("Radiation energy loss leaves a particle below rest energy")
            p0_after = np.sqrt(new_energy**2 - m_e_GeV**2)
            momentum_after = np.sqrt((particle_energy - energy_correction)**2 - m_e_GeV**2)
            factor = momentum_after / momentum_before * p0_before / p0_after
            p_array.px()[:] *= factor
            p_array.py()[:] *= factor
            p_array.p()[:] *= p0_before / p0_after
            energy = new_energy
        p_array.E = energy

    u, L = _track_non_undulator_segment(
        p_array,
        non_u,
        L,
        energy,
        accuracy,
    )
    if u is not None:
        U.append(u)
        E.append(energy)

    if not U:
        raise ValueError(
            "Radiation lattice contains no trackable nonzero-length elements"
        )

    # for u in U:
    #     print("here", len(u[4::9, 0]))
    #     plt.plot(u[4::9, :], u[0::9, :])
    # plt.show()
    return U, E


if __name__ == "__main__":
    quantum_diffusion(17.5, 4., 0.04, 200., quantum_diff=True)
    x = np.linspace(0, 1, 4)
    xnew = x2xgaus(x)
    print(x)
    print(xnew)
