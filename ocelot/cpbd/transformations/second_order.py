import numpy as np

from ocelot.cpbd.tm_utils import SecondOrderMult, transform_vec_ent, transform_vec_ext, \
    transfer_map_rotation, sym_matrix
from ocelot.cpbd.transformations.transformation import Transformation, TMTypes
from ocelot.cpbd.elements.element import Element


class SecondTM(Transformation):
    """[summary]
    Implementation of the second order transformation.
    The concrete element atom have to implement: 
    create_second_order_main_params(self, energy: float, delta_length: float) -> FirstOrderPrams
    If the element has edges is also have to implement:
    create_second_order_entrance_params(self, energy: float, delta_length: float) -> FirstOrderPrams
    create_second_order_exit_params(self, energy: float, delta_length: float) -> FirstOrderPrams
    """

    def __init__(self, create_tm_param_func, delta_e_func, tm_type: TMTypes, length: float, delta_length: float = 0.0) -> None:
        self.multiplication = SecondOrderMult().tmat_multip
        super().__init__(create_tm_param_func, delta_e_func, tm_type, length, delta_length)

    @classmethod
    def from_element(cls, element: Element, tm_type: TMTypes = TMTypes.MAIN, delta_l=None,  **params):
        return cls.create(entrance_tm_params_func=element.create_second_order_entrance_params if element.has_edge else None,
                          delta_e_func=element.create_delta_e,
                          main_tm_params_func=element.create_second_order_main_params,
                          exit_tm_params_func=element.create_second_order_exit_params if element.has_edge else None,
                          tm_type=tm_type, length=element.l, delta_length=delta_l)

    def t_apply(self, energy, X, U5666=0.):
        params = self.get_params(energy)
        R, T = params.R, params.T
        if params.tilt != 0:
            R = params.get_rotated_R()
            T = params.get_rotated_T()
        self.multiplication(X, R, T)
        X[:] = np.add(X, params.B)
        return X

    def map_function(self, X, energy: float):
        return self.t_apply(energy, X)

    def calculate_Tb(self, energy) -> np.ndarray:
        """
        Calculates the Tb matrix which is needed to calculate the transformation matrix.
        @return: Tb matrix
        """
        raise NotImplementedError("Not implemented yet")
        T_tilt = transfer_map_rotation(self.r_z_no_tilt(self.length, energy),
                                       self.t_mat_z_e(self.length, energy), self.tilt)[1]
        return sym_matrix(T_tilt)
