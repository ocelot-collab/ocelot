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

