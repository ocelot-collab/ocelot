import numpy as np
from scipy import interpolate
from ocelot.common.ocelog import *

_logger = logging.getLogger(__name__)
try:
    import numba as nb

    nb_flag = True
except:
    _logger.info("beam.py: module NUMBA is not installed. Install it to speed up calculation")
    nb_flag = False


def simple_filter(x, p, iter):
    n = len(x)
    if iter == 0:
        y = x
        return y
    for k in range(iter):

        y = np.zeros(n)
        for i in range(n):
            i0 = i - p
            if i0 < 0:
                i0 = 0
            i1 = i + p
            if i1 > n - 1:
                i1 = n - 1
            s = 0
            for j in range(i0, i1 + 1):
                s = s + x[j]
            y[i] = s / (i1 - i0 + 1)

        x = y
    return y


def interp1(x, y, xnew, k=1):
    if len(xnew) > 0:
        tck = interpolate.splrep(x, y, k=k)
        ynew = interpolate.splev(xnew, tck, der=0)
    else:
        ynew = []
    return ynew

def sortcols(array2d, row):
    """
    function sorts an 2D array columns based on specific row
    Example: ParticleArray.rparticles is [6, N] array where 5th row is tau coordinate. In case one needs to sort
             particles with respect to tau coordinate the code is sortcols(ParticleArray.rparticles, row=4)
    :param array:
    :param row:
    :return:
    """
    return array2d[:, array2d[row].argsort()]


def convmode_py(A, B, mode):
    if mode == 2:
        C = np.convolve(A, B)
    else:  # if mode == 1:
        i = np.int_(np.floor(len(B) * 0.5))
        n = len(A)
        C = np.zeros(n)
        C1 = np.convolve(A, B)
        C[:n] = C1[i:n + i]
    return C


convmode = convmode_py if not nb_flag else nb.jit(nopython=True)(convmode_py)

def moments(x, y, cut=0):
    n = len(x)
    # inds = np.arange(n)
    mx = np.mean(x)
    my = np.mean(y)
    x = x - mx
    y = y - my
    x2 = x * x
    mxx = np.sum(x2) / n
    y2 = y * y
    myy = np.sum(y2) / n
    xy = x * y
    mxy = np.sum(xy) / n

    emitt = np.sqrt(mxx * myy - mxy * mxy)

    if cut > 0:
        # inds=[]
        beta = mxx / emitt
        gamma = myy / emitt
        alpha = mxy / emitt
        emittp = gamma * x2 + 2. * alpha * xy + beta * y2
        inds0 = np.argsort(emittp)
        n1 = np.round(n * (100 - cut) / 100)
        inds = inds0[0:n1]
        mx = np.mean(x[inds])
        my = np.mean(y[inds])
        x1 = x[inds] - mx
        y1 = y[inds] - my
        mxx = np.sum(x1 * x1) / n1
        myy = np.sum(y1 * y1) / n1
        mxy = np.sum(x1 * y1) / n1
        emitt = np.sqrt(mxx * myy - mxy * mxy)
    return mx, my, mxx, mxy, myy, emitt


def m_from_twiss(Tw1, Tw2):
    """Transport matrix between two uncoupled Twiss parameter sets."""
    b1 = Tw1[1]
    a1 = Tw1[0]
    psi1 = Tw1[2]
    b2 = Tw2[1]
    a2 = Tw2[0]
    psi2 = Tw2[2]

    psi = psi2 - psi1
    cosp = np.cos(psi)
    sinp = np.sin(psi)
    M = np.zeros((2, 2))
    M[0, 0] = np.sqrt(b2 / b1) * (cosp + a1 * sinp)
    M[0, 1] = np.sqrt(b2 * b1) * sinp
    M[1, 0] = ((a1 - a2) * cosp - (1 + a1 * a2) * sinp) / np.sqrt(b2 * b1)
    M[1, 1] = np.sqrt(b1 / b2) * (cosp - a2 * sinp)
    return M


def beam_matching(parray, bounds, x_opt, y_opt, remove_offsets=True, slice=None):
    """
    Match a ``ParticleArray`` to the target transverse Twiss parameters.

    Parameters
    ----------
    parray : ParticleArray
        Particle array to transform in place.
    bounds : sequence of float
        ``[start, stop]`` window in ``tau`` rms units used to estimate the
        source Twiss parameters when ``slice`` is ``None``.
    x_opt : sequence of float
        Target horizontal optics ``[alpha, beta, mu]``.
    y_opt : sequence of float
        Target vertical optics ``[alpha, beta, mu]``.
    remove_offsets : bool, optional
        If ``True``, subtract the mean transverse offsets before matching.
    slice : str or None, optional
        If set, match using the selected beam slice (for example ``"Imax"`` or
        ``"Emax"``) and ignore ``bounds``.

    Returns
    -------
    numpy.ndarray
        The transformed ``parray.rparticles`` array.
    """
    particles = parray.rparticles
    pd = np.zeros((int(particles.size / 6), 6))
    dx = 0.0
    dxp = 0.0
    dy = 0.0
    dyp = 0.0
    if remove_offsets:
        dx = np.mean(particles[0])
        dxp = np.mean(particles[1])
        dy = np.mean(particles[2])
        dyp = np.mean(particles[3])

    pd[:, 0] = particles[0] - dx
    pd[:, 1] = particles[1] - dxp
    pd[:, 2] = particles[2] - dy
    pd[:, 3] = particles[3] - dyp
    pd[:, 4] = particles[4]
    pd[:, 5] = particles[5]

    z0 = np.mean(pd[:, 4])
    sig0 = np.std(pd[:, 4])
    inds = np.argwhere((z0 + sig0 * bounds[0] <= pd[:, 4]) * (pd[:, 4] <= z0 + sig0 * bounds[1]))

    mx, mxs, mxx, mxxs, mxsxs, emitx0 = moments(pd[inds, 0], pd[inds, 1])
    beta_x = mxx / emitx0
    alpha_x = -mxxs / emitx0

    my, mys, myy, myys, mysys, emity0 = moments(pd[inds, 2], pd[inds, 3])
    beta_y = myy / emity0
    alpha_y = -myys / emity0

    if slice is not None:
        tws = parray.get_twiss_from_slice(
            slice=slice,
            nparts_in_slice=5000,
            smooth_param=0.05,
            filter_base=2,
            filter_iter=2,
        )
        beta_x = tws.beta_x
        alpha_x = tws.alpha_x
        beta_y = tws.beta_y
        alpha_y = tws.alpha_y

    Mx = m_from_twiss([alpha_x, beta_x, 0], x_opt)
    particles[0] = Mx[0, 0] * pd[:, 0] + Mx[0, 1] * pd[:, 1]
    particles[1] = Mx[1, 0] * pd[:, 0] + Mx[1, 1] * pd[:, 1]

    My = m_from_twiss([alpha_y, beta_y, 0], y_opt)
    particles[2] = My[0, 0] * pd[:, 2] + My[0, 1] * pd[:, 3]
    particles[3] = My[1, 0] * pd[:, 2] + My[1, 1] * pd[:, 3]
    return particles
