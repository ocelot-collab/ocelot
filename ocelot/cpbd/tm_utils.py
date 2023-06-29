import logging

import numpy as np

from ocelot.cpbd.r_matrix import rot_mtx

_logger = logging.getLogger(__name__)

try:
    import numba as nb

    nb_flag = True
except ImportError as error:
    _logger.debug(" optics.py: module NUMBA is not installed. Install it to speed up calculation")
    nb_flag = False


class SecondOrderMult:
    """
    The class includes three different methods for transforming the particle
    coordinates:

    1. NUMBA module - DEACTIVATED, because new numpy implementation shows higher performance.
        Slightly faster than NUMPY for simulations with a large number of time
        steps. Uses full matrix multiplication.

    2. NUMPY module
        Base method to be used.
        Uses full matrix multiplication.

    """

    def __init__(self):
        self.full_matrix = False

        if nb_flag and False:
            self.tmat_multip = nb.njit()(SecondOrderMult.numba_apply)
        else:
            self.tmat_multip = SecondOrderMult.numpy_apply

    @staticmethod
    def numba_apply(X, R, T):
        Xcopy = np.copy(X)
        for n in range(X.shape[1]):
            for i in range(6):
                x_new = 0
                for j in range(6):
                    x_new += R[i, j] * Xcopy[j, n]
                    for k in range(6):
                        x_new += T[i, j, k] * Xcopy[j, n] * Xcopy[k, n]
                X[i, n] = x_new

    @staticmethod
    def numpy_apply(X, R, T):
        X[:] = np.matmul(R, X) + np.einsum('ijk,j...,k...->i...', T, X, X)


def transform_vec_ent(X, dx, dy, tilt):
    rotmat = rot_mtx(tilt)
    x_add = np.add(X, np.array([[-dx], [0.], [-dy], [0.], [0.], [0.]]))
    X[:] = np.dot(rotmat, x_add)[:]
    return X


def transform_vec_ext(X, dx, dy, tilt):
    rotmat = rot_mtx(-tilt)
    x_tilt = np.dot(rotmat, X)
    X[:] = np.add(x_tilt, np.array([[dx], [0.], [dy], [0.], [0.], [0.]]))[:]
    return X


def transfer_maps_mult_py(Ba, Ra, Ta, Bb, Rb, Tb):
    """
    cell = [A, B]
    Rc = Rb * Ra
    :param Ra:
    :param Ta:
    :param Rb:
    :param Tb:
    :param sym_flag:
    :return:
    """
    Bc = np.dot(Rb, Ba) + Bb
    Rc = np.dot(Rb, Ra)
    Tc = np.zeros((6, 6, 6))
    for i in range(6):
        for j in range(6):
            for k in range(6):
                t1 = 0.
                t2 = 0.
                for l in range(6):
                    t1 += Rb[i, l] * Ta[l, j, k]
                    for m in range(6):
                        t2 += Tb[i, l, m] * Ra[l, j] * Ra[m, k]
                Tc[i, j, k] = t1 + t2
    return Bc, Rc, Tc


def transfer_maps_mult_full_py(Ba, Ra, Ta, Bb, Rb, Tb):
    """
    cell = [A, B]
    Rc = Rb * Ra
    :param Ba:
    :param Ra:
    :param Ta:
    :param Bb:
    :param Rb:
    :param Tb:
    :return:
    """
    # zero order, multiplication was rewritten due to limitations of Numba
    Bc = np.zeros((6, 1))  # Bb + np.dot(Rb, Ba) + np.dot(Ba.T, np.dot(Tb, Ba))[0]
    for i in range(6):
        b1 = Bb[i, 0]
        b2 = 0.
        for l in range(6):
            b1 += Rb[i, l] * Ba[l, 0]
            for m in range(6):
                b2 += Ba[m, 0] * Ba[l, 0] * Tb[i, m, l]
        Bc[i, 0] += b1 + b2

    # first order
    Rc = np.dot(Rb, Ra)
    for i in range(6):
        for j in range(6):
            r1 = 0.
            r2 = 0.
            for l in range(6):
                for m in range(6):
                    r1 += Ra[l, j] * Ba[m, 0] * Tb[i, m, l]
                    r2 += Ba[l, 0] * Ra[m, j] * Tb[i, m, l]
            Rc[i, j] += r1 + r2

    # second order
    Tc = np.zeros((6, 6, 6))
    for i in range(6):
        for j in range(6):
            for k in range(6):
                t1 = 0.
                t2 = 0.
                t3 = 0.
                t4 = 0.
                for l in range(6):
                    t1 += Rb[i, l] * Ta[l, j, k]
                    for m in range(6):
                        t2 += Tb[i, l, m] * Ra[l, j] * Ra[m, k]
                        t3 += Ta[l, j, k] * Ba[m, 0] * Tb[i, m, l]
                        t4 += Ba[l, 0] * Ta[m, j, k] * Tb[i, m, l]
                Tc[i, j, k] = t1 + t2 + t3 + t4
    return Bc, Rc, Tc


transfer_maps_mult = transfer_maps_mult_py if nb_flag is not True else nb.jit(transfer_maps_mult_py, nopython=True)


def transfer_map_rotation(R, T, tilt):
    rotmat = rot_mtx(tilt)
    Bc, Rc, Tc = transfer_maps_mult(Ba=np.zeros((6, 1)), Ra=rotmat, Ta=np.zeros((6, 6, 6)), Bb=np.zeros((6, 1)), Rb=R,
                                    Tb=T)
    B, R, T = transfer_maps_mult(Ba=Bc, Ra=Rc, Ta=Tc, Bb=np.zeros((6, 1)), Rb=rot_mtx(-tilt), Tb=np.zeros((6, 6, 6)))
    return R, T


def sym_matrix(T):
    for i in range(6):
        for j in range(6):
            for k in range(j, 6):
                if j != k:
                    a = (T[i, j, k] + T[i, k, j]) / 2.
                    T[i, k, j] = a
                    T[i, j, k] = a
    return T


def unsym_matrix(T):
    for i in range(6):
        for j in range(6):
            for k in range(j, 6):
                if j != k:
                    a = T[i, j, k] * 2.
                    T[i, k, j] = 0
                    T[i, j, k] = a
    return T


def map_transform_with_offsets(M, R, T):
    """
    Calculate new R matrix and B vector for an element with misalignment. Second order matrix is unchanged.

    :param M: array(6, 1), element offset, e.g. for dx and dy offsets - [[-dx], [0], [-dy], [0], [0], [0]]
    :param R: array(6, 6), First order matrix
    :param T: array(6, 6, 6), Second Order Matrix
    :return: B, R - zero and first order matrices. T matrix does not change.
    """
    # first order matrix transformation
    R_off = np.copy(R)
    for i in range(6):
        for j in range(6):
            for l in range(6):
                R_off[i, j] += M[l, 0] * (T[i, j, l] + T[i, l, j])

    # zero order matrix transformation
    B_off = np.dot(R, M) - M + np.dot(np.dot(M.T, T), M)[0]

    return B_off, R_off
