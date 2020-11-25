from copy import copy

import numpy as np

from ocelot.cpbd.transformations.first_order import TransferMap
from ocelot.cpbd.transformations.optics import SecondOrderMult, transform_vec_ext, transform_vec_ent, \
    transfer_map_rotation, sym_matrix
from ocelot.cpbd.r_matrix import rot_mtx


class SecondTM(TransferMap):
    def __init__(self, r_z_no_tilt, t_mat_z_e):
        TransferMap.__init__(self)
        self.r_z_no_tilt = r_z_no_tilt
        self.t_mat_z_e = t_mat_z_e

        self.multiplication = None
        self.map = lambda X, energy: self.t_apply(self.r_z_no_tilt(self.length, energy),
                                                  self.t_mat_z_e(self.length, energy), X, self.dx, self.dy, self.tilt)

        self.R_tilt = lambda energy: np.dot(np.dot(rot_mtx(-self.tilt), self.r_z_no_tilt(self.length, energy)),
                                            rot_mtx(self.tilt))

        self.T_tilt = lambda energy: transfer_map_rotation(self.r_z_no_tilt(self.length, energy),
                                                           self.t_mat_z_e(self.length, energy), self.tilt)[1]

    @classmethod
    def create_from_element(cls, element, params=None):
        T_z_e = element.get_T_z_e_func()
        tm = cls(r_z_no_tilt=element.create_r_matrix(), t_mat_z_e=T_z_e)
        tm.multiplication = SecondOrderMult().tmat_multip
        return tm

    def calculate_Tb(self, energy) -> np.ndarray:
        """
        Calculates the Tb matrix which is needed to claculate the transfromation matrix.
        Note: The calculation of the Tb matrix is different between first order TM and second order TM.
        @return: Tb matrix
        """
        Tb = np.copy(self.T_tilt(energy))
        Tb = sym_matrix(Tb)
        return Tb

    def t_apply(self, R, T, X, dx, dy, tilt, U5666=0.):
        if dx != 0 or dy != 0 or tilt != 0:
            X = transform_vec_ent(X, dx, dy, tilt)
        self.multiplication(X, R, T)
        if dx != 0 or dy != 0 or tilt != 0:
            X = transform_vec_ext(X, dx, dy, tilt)

        return X

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.T = lambda s, energy: m.t_mat_z_e(s, energy)
        m.delta_e = m.delta_e_z(s)
        m.map = lambda X, energy: m.t_apply(m.r_z_no_tilt(s, energy), m.t_mat_z_e(s, energy), X, m.dx, m.dy, m.tilt)
        return m
