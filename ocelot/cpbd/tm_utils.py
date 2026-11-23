import numpy as np

from ocelot.cpbd.r_matrix import rot_mtx


class SecondOrderMult:
    """Apply a second-order transfer map to a particle array."""

    def __init__(self):
        self.full_matrix = False
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
    """Compose A followed by B when Ba or Tb is zero."""
    Tc = np.einsum("il,ljk->ijk", Rb, Ta) + Ra.T @ Tb @ Ra
    return Bb + Rb @ Ba, Rb @ Ra, Tc


def transfer_maps_mult_full_py(Ba, Ra, Ta, Bb, Rb, Tb):
    """Compose A followed by B through second order about zero input to A."""
    # Evaluate B(x) = Bb + Rb x + Tb[x,x] and its Jacobian at x = Ba.
    # Accept symmetric tensors and the one-sided mixed terms used by some maps.
    tensor_offset = Tb @ Ba[:, 0]
    orbit = Bb + Rb @ Ba + tensor_offset @ Ba
    jacobian = Rb + tensor_offset + Tb.swapaxes(1, 2) @ Ba[:, 0]
    quadratic = np.einsum("il,ljk->ijk", jacobian, Ta) + Ra.T @ Tb @ Ra
    return orbit, jacobian @ Ra, quadratic


def transfer_maps_mult(Ba, Ra, Ta, Bb, Rb, Tb):
    """Compose the quadratic maps in beam order, evaluating the second map at Ba."""
    if Ba.any() and Tb.any():
        return transfer_maps_mult_full_py(Ba, Ra, Ta, Bb, Rb, Tb)
    return transfer_maps_mult_py(Ba, Ra, Ta, Bb, Rb, Tb)


def transfer_map_rotation(R, T, tilt):
    rotation = rot_mtx(tilt)
    # Transform both tensor inputs to magnet axes and its output to beam axes.
    return (rotation.T @ R @ rotation,
            np.einsum("il,ljk->ijk", rotation.T, rotation.T @ T @ rotation))


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
