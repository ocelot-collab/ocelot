import logging

import numba as nb
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
            # print("SecondTM: NUMBA")
            self.tmat_multip = nb.njit()(SecondOrderMult.numba_apply)
        else:
            # print("SecondTM: Numpy")
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


def transfer_maps_mult_py(Ra, Ta, Rb, Tb):
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
    return Rc, Tc


transfer_maps_mult = transfer_maps_mult_py if nb_flag is not True else nb.jit(transfer_maps_mult_py)


def transfer_map_rotation(R, T, tilt):
    rotmat = rot_mtx(tilt)
    Rc, Tc = transfer_maps_mult(Ra=rotmat, Ta=np.zeros((6, 6, 6)), Rb=R, Tb=T)
    R, T = transfer_maps_mult(Ra=Rc, Ta=Tc, Rb=rot_mtx(-tilt), Tb=np.zeros((6, 6, 6)))
    return R, T


def sym_matrix(T):
    for i in range(6):
        for j in range(6):
            for k in range(j, 6):
                if j != k:
                    a = T[i, j, k] / 2.
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