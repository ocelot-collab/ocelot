"""
definition of particles, beams and trajectories
"""
import os
from copy import deepcopy
import pandas as pd
import numpy as np
from typing import Iterable

import ocelot.common.globals as glb
from ocelot.common.math_op import find_nearest_idx, invert_cdf
from ocelot.common.ocelog import *
from . import analysis
from . import particle
from . import core
from .beam_utils import beam_matching, m_from_twiss


#TypeParticleArray = TypeVar("TypeParticleArray", bound="ParticleArray")

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


class BeamFormFactor:
    '''
    contains and calculates electron beam form-factor (fourier transform of currenta profile)
    from the electron beam "BeamArray" object
    '''

    def __init__(self, beam_array=None):
        self.cfactor = None  #modulus of form-factor
        self.frequency = None  #frequency in Hz
        self.beam_array = beam_array  #original beam file

        if self.beam_array is not None:
            self.beam_array.sort()
            self.beam_array.equidist()
            self.calc()

    def __len__(self):
        return np.size(self.cfactor)

    # def current_profile(self):
    #     return self.beam_array.s, self.beam_array.I

    def calc(self):
        '''
        calculates the form-factor and populates self.modulus and self.frequency
        '''
        I_norm = self.beam_array.I / self.beam_array.charge()
        self.cfactor = np.fft.fft(I_norm)
        ds = np.abs(self.beam_array.s[1] - self.beam_array.s[0])
        self.frequency = np.fft.fftfreq(len(self), ds / glb.speed_of_light)
        idx = len(self) // 2

        self.cfactor = self.cfactor[:idx]
        self.frequency = self.frequency[:idx]


def recalculate_ref_particle(p_array):
    pref = np.sqrt(p_array.E ** 2 / glb.m_e_GeV ** 2 - 1) * glb.m_e_GeV
    Enew = p_array.p()[0] * pref + p_array.E
    s_new = p_array.s - p_array.tau()[0]
    p_array.rparticles[5, :] -= p_array.p()[0]
    p_array.rparticles[4, :] -= p_array.tau()[0]
    p_array.E = Enew
    p_array.s = s_new
    return p_array


def twiss_parray_slice(parray, slice="Imax", nparts_in_slice=5000, smooth_param=0.05, filter_base=2, filter_iter=2):
    """
    Function calculates twiss parameters in a beam slice

    :param parray: ParticleArray
    :param slice: "Imax" or "Emax" or center of bunch
    :param nparts_in_slice: 5000, nparticles in the slice (in moving window)
    :param smooth_param: 0.01, smoothing parameters to calculate the beam current: smooth_param = m_std * np.std(p_array.tau())
    :param filter_base: 2, filter parameter in the func: simple_filter
    :param filter_iter: 2, filter parameter in the func: simple_filter
    :return: Twiss
    """
    slice_params = analysis.global_slice_analysis(parray, nparts_in_slice=nparts_in_slice, smooth_param=smooth_param,
                                         filter_base=filter_base, filter_iter=filter_iter)

    if slice == "Imax":
        ind0 = np.argmax(slice_params.I)
    elif slice == "Emax":
        ind0 = np.argmax(slice_params.me)
    else:
        ind0 = np.argsort(np.abs(slice_params.s))[0]
    tws = slice_params.extract_slice(ind0)
    return tws


'''
beam funcions proposed
'''


def parray2beam(parray, step=1e-7):
    '''
    reads ParticleArray()
    returns BeamArray()
    step [m] - long. size ob bin to calculate distribution parameters
    '''
    _logger.info('calculating electron beam distribution from particle array')

    part_c = parray.q_array[0]  # fix for general case  # charge per particle
    t_step = step / glb.speed_of_light
    t = parray.tau() / glb.speed_of_light
    t_min = min(t)
    t_max = max(t)
    t_window = t_max - t_min
    npoints = int(t_window / t_step)
    t_step = t_window / npoints
    beam = core.BeamArray()
    for parm in ['I',
                 's',
                 'emit_x',
                 'emit_y',
                 'beta_x',
                 'beta_y',
                 'alpha_x',
                 'alpha_y',
                 'x',
                 'y',
                 'xp',
                 'yp',
                 'E',
                 'sigma_E',
                 ]:
        setattr(beam, parm, np.zeros((npoints - 1)))

    for i in range(npoints - 1):
        indices = (t > t_min + t_step * i) * (t < t_min + t_step * (i + 1))
        beam.s[i] = (t_min + t_step * (i + 0.5)) * glb.speed_of_light

        if np.sum(indices) > 2:
            e0 = parray.E * 1e9
            p0 = np.sqrt((e0 ** 2 - glb.m_e_eV ** 2) / glb.speed_of_light ** 2)
            p = parray.rparticles[5][indices]  # deltaE / average_impulse / speed_of_light
            dist_e = (p * p0 * glb.speed_of_light + e0)
            dist_x = parray.rparticles[0][indices]
            dist_y = parray.rparticles[2][indices]
            dist_xp = parray.rparticles[1][indices]
            dist_yp = parray.rparticles[3][indices]

            beam.I[i] = np.sum(indices) * part_c / t_step
            beam.E[i] = np.mean(dist_e) * 1e-9
            beam.sigma_E[i] = np.std(dist_e) * 1e-9

            dist_x_m = np.mean(dist_x)
            dist_y_m = np.mean(dist_y)
            dist_xp_m = np.mean(dist_xp)
            dist_yp_m = np.mean(dist_yp)

            beam.x[i] = dist_x_m
            beam.y[i] = dist_y_m
            # g = beam.E[i] / m_e_GeV
            # p = np.sqrt(g**2 - 1)
            beam.xp[i] = dist_xp_m
            beam.yp[i] = dist_yp_m

            dist_x -= dist_x_m
            dist_y -= dist_y_m
            dist_xp -= dist_xp_m
            dist_yp -= dist_yp_m

            beam.emit_x[i] = np.sqrt(np.mean(dist_x ** 2) * np.mean(dist_xp ** 2) - np.mean(dist_x * dist_xp) ** 2)
            beam.emit_y[i] = np.sqrt(np.mean(dist_y ** 2) * np.mean(dist_yp ** 2) - np.mean(dist_y * dist_yp) ** 2)
            beam.beta_x[i] = np.mean(dist_x ** 2) / beam.emit_x[i]
            beam.beta_y[i] = np.mean(dist_y ** 2) / beam.emit_y[i]
            beam.alpha_x[i] = -np.mean(dist_x * dist_xp) / beam.emit_x[i]
            beam.alpha_y[i] = -np.mean(dist_y * dist_yp) / beam.emit_y[i]

    idx = np.where(np.logical_or.reduce(
        (beam.I == 0, beam.E == 0, beam.beta_x > np.mean(beam.beta_x) * 100, beam.beta_y > np.mean(beam.beta_y) * 100)))
    del beam[idx]

    if hasattr(parray, 'filePath'):
        beam.filePath = parray.filePath + '.beam'
    return (beam)


def cov_matrix_from_twiss(ex, ey, sigma_tau, sigma_p, **twiss):
    """Generate a covariance matrix from Twiss parameters, dispersions and
    emittances (horizontal) and standard deviations (longitudinal).  No
    correlations between tau and the other coordinates are present in this
    parametrisation.

    :param ex: Geometric emittance in x-plane
    :param ey: Geometric emittance in y-plane
    :param sigma_tau: Standard deviation of tau (c*t)
    :param sigma_p: Standard deviation of cannonical coordinate p=dE/(c*p0).
    :param twiss: Horizontal Twiss parameters.  Required: alpha_x, beta_x,
        alpha_y, beta_y.  Optional (set to 0 if missing): dispersions dx,
        dpx, dy, dpy.
    :return: 6x6 correlation matrix.


    """
    alpha_x = twiss["alpha_x"]
    beta_x = twiss["beta_x"]
    alpha_y = twiss["alpha_y"]
    beta_y = twiss["beta_y"]
    dx = twiss.get("dx", 0)
    dpx = twiss.get("dpx", 0)
    dy = twiss.get("dy", 0)
    dpy = twiss.get("dpy", 0)
    # X block, y block, xy upper block, xy lower block
    xb = _horizontal_2x2_elements(ex, alpha_x, beta_x, dx, dpx, sigma_p)
    yb = _horizontal_2x2_elements(ey, alpha_y, beta_y, dy, dpy, sigma_p)
    xyu = _horizontal_coupling_elements(dx, dy, dpx, dpy, sigma_p)
    xyl = np.array(
        _horizontal_coupling_elements(dx, dy, dpx, dpy, sigma_p)
    ).T
    sp2 = sigma_p ** 2
    return np.array([[xb[0, 0], xb[0, 1], xyu[0, 0], xyu[0, 1], 0., dx * sp2],
                     [xb[1, 0], xb[1, 1], xyu[1, 0], xyu[1, 1], 0., dpx * sp2],
                     [xyl[0, 0], xyl[0, 1], yb[0, 0], yb[0, 1], 0., dy * sp2],
                     [xyl[1, 0], xyl[1, 1], yb[1, 0], yb[1, 1], 0., dpy * sp2],
                     [0., 0, 0, 0, sigma_tau ** 2, 0.0],
                     [dx * sp2, dpx * sp2, dy * sp2, dpy * sp2, 0, sp2]])


def cov_matrix_to_parray(mean, cov, energy, charge, nparticles):
    """Generate a ParticleArray instance using a covariance matrix.

    :param mean: 1-D list of 6 means of the particle distributions.
    :param cov: 6x6 covariance matrix.
    :param energy: Beam energy in GeV.
    :param charge: Total beam charge in Coulombs.
    :param nparticles: Number of particles to populate the ParticleArray
        instance with.
    :return: ParticleArray with given charge and energy populated with
        nparticles and the particle distribution having the correct means
        and covariances.
    :rtype: ParticleArray

    """

    p_array = particle.ParticleArray()
    p_array.E = energy
    p_array.rparticles = np.random.multivariate_normal(mean, cov, nparticles).T
    p_array.q_array = np.ones(nparticles) * charge / nparticles
    return p_array


def _horizontal_2x2_elements(emit, alpha, beta, disp, disp_p, sigma_p):
    """2x2 correlation matrix between x/py  and y/py"""
    gamma = (1 + alpha ** 2) / beta
    offdiag = -emit * alpha + disp * disp_p * sigma_p ** 2
    return np.array([[emit * beta + (disp * sigma_p) ** 2, offdiag],
                     [offdiag, emit * gamma + (disp_p * sigma_p) ** 2]])


def _horizontal_coupling_elements(disp_x, disp_y, disp_px, disp_py, sigma_p):
    """generate cov matrix elements for correlations between horz. and
    vertical."""
    sigp2 = sigma_p ** 2
    return np.array([[disp_x * disp_y * sigp2, disp_x * disp_py * sigp2],
                     [disp_px * disp_y * sigp2, disp_px * disp_py * sigp2]])


def optics_from_moments(mean, cov_matrix, energy=None):
    """Calculate the beam optics from the mean and covariance matrix.

    :param mean: 1x6 array of means
    :param cov_matrix: 6x6 matrix of covariances between the particle vectors.
    :param energy: Energy to additionally calculate the normalised
        emittances from the geometric emittances, optional.

    """
    r = core.Twiss()
    r.x, r.xp, r.y, r.yp, r.tau, r.p = mean
    r.Dx, r.Dxp, r.Dy, r.Dyp = _dispersions_from_cov_matrix(cov_matrix)
    sigp2 = cov_matrix[5, 5]
    r.emit_x, r.alpha_x, r.beta_x, _ = _dispersionless_twiss_parameters(
        cov_matrix[0:2, 0:2],
        r.Dx,
        r.Dxp,
        sigp2
    )
    r.emit_y, r.alpha_y, r.beta_y, _ = _dispersionless_twiss_parameters(
        cov_matrix[2:4, 2:4],
        r.Dy,
        r.Dyp,
        sigp2
    )

    if energy is not None:
        r.E = energy

    return r


def _dispersionless_twiss_parameters(submatrix, dx, dpx, sigp2):
    """Calculate the emittance and twiss parameters from the 2x2 covariance
    matrix accounting for the increase in the spot size due to the
    dispersion.

    """
    x = submatrix[0, 0] - dx ** 2 * sigp2
    px = submatrix[1, 1] - dpx ** 2 * sigp2
    xpx = submatrix[0, 1] - dx * dpx * sigp2
    emittance = np.sqrt(x * px - xpx * xpx)
    beta = x / emittance
    alpha = -xpx / emittance
    gamma = px / emittance

    return emittance, alpha, beta, gamma


def _dispersions_from_cov_matrix(cov_matrix):
    """Calculate the dispersions from the provided 6x6 covariance matrix."""
    sigp2 = cov_matrix[5, 5]
    if sigp2 == 0:
        return 0, 0, 0, 0
    dx = cov_matrix[0, 5] / sigp2
    dpx = cov_matrix[1, 5] / sigp2
    dy = cov_matrix[2, 5] / sigp2
    dpy = cov_matrix[3, 5] / sigp2
    return dx, dpx, dy, dpy


def moments_from_parray(parray, dispersions=None):
    """Calculate the central moments and covariances from the given
    ParticleArray.  Correct for dispersion by providing either a Twiss
    instance, a list of four dispersion [dx, dxp, dy, dyp].
    """
    rpart = parray.rparticles
    try:
        dx = dispersions.Dx
        dxp = dispersions.Dxp
        dy = dispersions.Dy
        dyp = dispersions.Dyp
    except AttributeError:
        try:
            dx, dxp, dy, dyp = dispersions
        except TypeError:
            dx, dxp, dy, dyp = 0, 0, 0, 0
    # TODO: ACTUALLY CORRECT FOR DISPERSION AND DO IT CLEANLY SOMEWHERE.

    return np.mean(rpart, axis=1, dtype=np.float64), np.cov(rpart)



def gauss_from_twiss(emit, beta, alpha):
    phi = 2 * np.pi * np.random.rand()
    u = np.random.rand()
    a = np.sqrt(-2 * np.log((1 - u)) * emit)
    x = a * np.sqrt(beta) * np.cos(phi)
    xp = -a / np.sqrt(beta) * (np.sin(phi) + alpha * np.cos(phi))
    return (x, xp)


def waterbag_from_twiss(emit, beta, alpha):
    phi = 2 * np.pi * np.random.rand()
    a = np.sqrt(emit) * np.random.rand()
    x = a * np.sqrt(beta) * np.cos(phi)
    xp = -a / np.sqrt(beta) * (np.sin(phi) + alpha * np.cos(phi))
    return (x, xp)


def ellipse_from_twiss(emit, beta, alpha):
    phi = 2 * np.pi * np.random.rand()
    # u = np.random.rand()
    # a = np.sqrt(-2*np.log( (1-u)) * emit)
    a = np.sqrt(emit)
    x = a * np.sqrt(beta) * np.cos(phi)
    xp = -a / np.sqrt(beta) * (np.sin(phi) + alpha * np.cos(phi))
    return (x, xp)

