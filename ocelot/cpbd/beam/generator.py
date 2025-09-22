import numpy as np
import ocelot.common.globals as glb
from ocelot.common.math_op import invert_cdf
from . import particle
from . import core
from . import beam

from ocelot.common.ocelog import *

_logger = logging.getLogger(__name__)


def generate_parray(sigma_x=1e-4, sigma_px=2e-5, sigma_y=None, sigma_py=None,
                    sigma_tau=1e-3, sigma_p=1e-4, chirp=0.01, charge=5e-9, nparticles=200000, energy=0.13,
                    tau_trunc=None, tws=None, shape="gauss",
                    quiet=None, quiet_method="stratified", quiet_alpha=1.0, quiet_jitter_fraction=0.2,
                    noise_Ne=None, kmin=2 * np.pi / 1e-5, kmax=2 * np.pi / 1e-6, rng=None):
    """
    Method to generate ParticleArray with gaussian distribution.

    Note: in ParticleArray, {x, px} and {y, py} are canonical coordinates. {tau, p} is not, to make it canonical
            the sign of "tau" should be flipped.

    Example:
    --------

    sigma_tau= 10e-6

    A1, A2, A3 = 0.5, 0.2, 1
    mu1, mu2, mu3 = -7.7e-08, 3.7e-06, -6.3e-06
    sigma1, sigma2, sigma3 = 7.2e-06, 9.2e-06, 3e-06

    f = lambda x:  A1*np.exp(-(x - mu1) ** 2 / (2. * sigma1 ** 2)) + A2*np.exp(-(x - mu2) ** 2 / (2. * sigma2 ** 2)) + A3*np.exp(-(x - mu3) ** 2 / (2. * sigma3 ** 2))

    tau = np.linspace(-5*sigma_tau, 5*sigma_tau, num=200)

    shape = [tau, f(tau)]

    parray = generate_parray(sigma_x=1e-4, sigma_px=2e-5, sigma_y=None, sigma_py=None,
                        sigma_tau=sigma_tau, sigma_p=1e-4, chirp=0.01, charge=250e-12, nparticles=500000, energy=0.13,
                        tau_trunc=None, tws=None, shape=shape)

    ---------


    :param sigma_x: std(x), x is horizontal cartesian coordinate.
    :param sigma_px: std(px), 'px' is conjugate momentum canonical momentum px/p0.
    :param sigma_y: std(y), y is vertical cartesian coordinate.
    :param sigma_py: std(py), 'py' is canonical momentum py/p0.
    :param sigma_tau: std(tau), "tau" = c*t
    :param sigma_p: std(p), 'p' is canonical momentum E/(c*p0)
    :param chirp: energy chirp [unitless], linear correlation - p_i += chirp * tau_i/sigma_tau
    :param charge: beam charge in [C], 5e-9 by default
    :param nparticles: number of particles, 200k by default
    :param energy: beam energy in [GeV], 0.13 [GeV]
    :param tau_trunc: None, if not [float] - truncated gauss distribution in "tau" direction.
    :param tws: None, if Twiss obj - the beam is matched to twiss params.
    :param shape: "gauss", "tri" - triangular, "rect" - shape of the beam current profile. Gaussian distribution by default.
                  Current profile also can be arbitrary with shap=[tau, f(tau)] - shape
    :param quiet : {'tau', 'x', 'y', 'px', 'py', 'p', 'all'} or None
    :param quiet_method : str, default="stratified"
                        Quiet-start method for tau:
                          - "stratified": use uniform quantile grid.
                          - "jittered": stratified + small random jitter.
                          - "sobol": Sobol low-discrepancy sequence (rounded to nearest power-of-two).
    :param quiet_alpha : float, default=1.0
                        Blend factor for tau quiet-start:
                          - 0.0 → keep original random tau.
                          - 1.0 → full quiet tau.
                          - Between 0–1 → partial suppression of noise.
    :param quiet_jitter_fraction : float, default=0.2
                                Jitter width as fraction of quantile spacing (only used if quiet_method="jittered").
    :param  noise_Ne : int or None
                        If set, add synthetic shot-noise to match a bunch of N_e electrons
                        (default: None = no noise injection)
    :param kmin : float. Band for spectral noise injection (in 1/m).
    :param kmax : float. Band for spectral noise injection (in 1/m).

    :return: ParticleArray
    """

    if isinstance(tws, core.Twiss) and np.all(
            np.array([tws.emit_x, tws.emit_y, tws.beta_x, tws.beta_y, tws.gamma_x, tws.gamma_y]) != 0):
        _logger.info("Twiss parameters have priority. sigma_{x, px, y, py} will be redefined")

        cov_x_px = tws.emit_x * np.array([[tws.beta_x, -tws.alpha_x],
                                          [-tws.alpha_x, tws.gamma_x]])
        hor_dist = np.random.multivariate_normal((0, 0), cov_x_px, nparticles)
        x = hor_dist[:, 0]
        px = hor_dist[:, 1]

        cov_y_py = tws.emit_y * np.array([[tws.beta_y, -tws.alpha_y],
                                          [-tws.alpha_y, tws.gamma_y]])
        vert_dist = np.random.multivariate_normal((0, 0), cov_y_py, nparticles)
        y = vert_dist[:, 0]
        py = vert_dist[:, 1]

        if tws.pp != 0:
            sigma_p = np.sqrt(tws.pp)
        if tws.E != 0:
            energy = tws.E
    else:
        if sigma_y is None:
            sigma_y = sigma_x
        if sigma_py is None:
            sigma_py = sigma_px

        x = np.random.randn(nparticles) * sigma_x
        px = np.random.randn(nparticles) * sigma_px
        y = np.random.randn(nparticles) * sigma_y
        py = np.random.randn(nparticles) * sigma_py

    if isinstance(shape, str):
        if shape == "gauss":
            f = lambda x: np.exp(-x ** 2 / (2. * sigma_tau ** 2))

            tau_trunc = 5 if tau_trunc is None else tau_trunc
            s = np.linspace(-tau_trunc * sigma_tau, tau_trunc * sigma_tau, num=500)
            shape = [s, f(s)]
        elif shape == "tri":
            s = np.linspace(-1 * sigma_tau, 1 * sigma_tau, num=500)
            f = lambda s: np.maximum(sigma_tau - np.abs(s), 0)
            shape = [s, f(s)]
        else:
            s = np.linspace(-1 * sigma_tau, 1 * sigma_tau, num=500)
            f = lambda s: np.where((s >= -sigma_tau) & (s <= sigma_tau), 1, 0)
            shape = [s, f(s)]

    inv_cdf = invert_cdf(y=shape[1], x=shape[0])
    tau = inv_cdf(np.random.rand(nparticles))

    #tau = np.random.randn(nparticles) * sigma_tau
    dp = np.random.randn(nparticles) * sigma_p

    # covariance matrix for [tau, p] for beam compression in BC
    # cov_t_p = [[1.30190131e-06, 2.00819771e-05],
    #           [2.00819771e-05, 3.09815718e-04]]
    # k = tau_p_cor*sigma_tau*sigma_p
    # cov_t_p = [[sigma_tau**2, k],
    #           [k, sigma_p**2]]
    # long_dist = np.random.multivariate_normal((0, 0), cov_t_p, nparticles)
    # tau = long_dist[:, 0]
    # dp = long_dist[:, 1]

    p_array = particle.ParticleArray(n=nparticles)
    p_array.E = energy  # GeV
    p_array.rparticles[0][:] = x[:]
    p_array.rparticles[1][:] = px[:]
    p_array.rparticles[2][:] = y[:]
    p_array.rparticles[3][:] = py[:]
    p_array.rparticles[4][:] = tau[:]
    p_array.rparticles[5][:] = dp[:]

    p_array.q_array = np.ones(nparticles) * charge / nparticles

    if quiet is not None:
        if quiet in ["tau", "all"]:
            p_array.quietify_tau(inv_cdf=inv_cdf, quiet_method=quiet_method, quiet_alpha=quiet_alpha,
                         quiet_jitter_fraction=quiet_jitter_fraction, noise_Ne=noise_Ne,
                         kmin=kmin, kmax=kmax, rng=rng)
        if quiet in ["x", "y", "px", "py", "p", "all"]:
            coords = [quiet] if quiet != "all" else ("x", "y", "px", "py", "dp")
            p_array.quietify_transverse(coords=coords, quiet_method=quiet_method, quiet_alpha=quiet_alpha,
                         quiet_jitter_fraction=quiet_jitter_fraction, rng=rng)
    if sigma_tau != 0:
        dp += chirp * p_array.tau() / sigma_tau
    p_array.rparticles[5][:] = dp[:]
    if isinstance(tws, core.Twiss):
        x_opt = [tws.alpha_x, tws.beta_x, tws.mux]
        y_opt = [tws.alpha_y, tws.beta_y, tws.muy]
        bounds = [-5, 5]
        beam.beam_matching(p_array, bounds, x_opt, y_opt, remove_offsets=True)

    return p_array


def generate_beam(E, I=5000, l_beam=3e-6, **kwargs):
    """
    Generates BeamArray object
    accepts arguments with the same names as BeamArray().parameters()
    I - current in Amps
    E - beam ebergy in GeV

    dE - rms energy spread in GeV
    emit_x, emit_n(both normalized), emit_xn, etc.
    shape - beam shape ('gaussian' of 'flattop')
    l_beam [m] - beam length in meters
    l_window [m] - window length in um
        by default: l_beam * 2 if flattop,
                    l_beam * 6 if gaussian,
    nslice - number of slices in the beam
    """

    _logger.info('generating electron beam distribution')

    beam = core.Beam()
    beam.E = E
    beam.tlen = l_beam / glb.speed_of_light * 1e15
    beam.I = I
    nslice = 100

    for key, value in kwargs.items():
        if (key in beam.__dict__ or key in beam.properties) and (key not in ['s', 'E', 'tlen', 'I']):
            setattr(beam, key, value)
        if key == 'emit':
            beam.emit_x = value
            beam.emit_y = value
        if key == 'emit_n':
            beam.emit_xn = value
            beam.emit_yn = value
        if key == 'beta':
            beam.beta_x = value
            beam.beta_y = value
        if key == 'nslice':
            nslice = value
        if key == 'dE':
            beam.dg = value / glb.m_e_GeV

    if 'l_window' not in kwargs:
        if beam.shape is ['gaussian', 'gauss', 'g']:
            l_window = l_beam * 6
        elif beam.shape in ['flattop', 'ft']:
            l_window = l_beam * 2
        else:
            raise ValueError('Beam() shape can be either "gaussian" or "flattop"')
    else:
        l_window = kwargs['l_window']

    beam_arr = beam.to_array(nslice, l_window)

    if 'chirp' in kwargs:
        if type(kwargs['chirp']) == list:
            beam_arr.add_chirp_poly(kwargs['chirp'])
        else:
            beam_arr.add_chirp(kwargs['chirp'])
    return beam_arr
