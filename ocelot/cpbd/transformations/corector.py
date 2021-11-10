import numpy as np

from ocelot.cpbd.transformations.second_order import SecondTM, TMTypes
from ocelot.cpbd.elements.element import Element


class CorrectorTM(SecondTM):
    """[summary]
    Implementation of the Corrector Transformation.
    The concrete element atom have to implement:
    create_cor_second_order_main_params(self, energy: float, delta_length: float) -> CorrectorSecondOrderParams
    """

    def __init__(self, create_tm_param_func, delta_e_func, tm_type: TMTypes, length: float,
                 delta_length: float) -> None:
        super().__init__(create_tm_param_func, delta_e_func, tm_type, length, delta_length=delta_length)

    @classmethod
    def from_element(cls, element: Element, tm_type: TMTypes = TMTypes.MAIN, delta_l=None):
        return cls.create(entrance_tm_params_func=None,
                          delta_e_func=element.create_delta_e,
                          main_tm_params_func=element.create_cor_second_order_main_params,
                          exit_tm_params_func=None,
                          tm_type=tm_type, length=element.l, delta_length=delta_l)

    def kick(self, X, energy):
        params = self.get_params(energy)

        X1 = self.t_apply(energy, X=X)

        X1 = np.add(X1, params.B)
        X[:] = X1[:]
        return X

    def map_function(self, X, energy: float):
        return self.kick(X, energy=energy)
