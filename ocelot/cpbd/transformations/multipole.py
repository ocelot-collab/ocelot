from math import factorial

import numpy as np

from ocelot.cpbd.transformations.transfer_map import TransferMap, TMTypes
from ocelot.cpbd.elements.element import Element


class MultipoleTM(TransferMap):
    """[summary]
    Implementation of the Multipole Transforamtion.
    The concrete element atom have to implement: 
    create_multipole_tm_main_params(self) -> MultipoleParams
    """
        
    def __init__(self, create_tm_param_func, delta_e_func, tm_type: TMTypes, length: float, delta_length: float = None, **params) -> None:
        super().__init__(create_tm_param_func, delta_e_func, tm_type, length, delta_length)

    @classmethod
    def from_element(cls, element: Element, tm_type: TMTypes = TMTypes.MAIN, delta_l=None, **params):
        return cls.create(entrance_tm_params_func=None,
                          delta_e_func=element.create_delta_e,
                          main_tm_params_func=element.create_multipole_tm_main_params,
                          exit_tm_params_func=None,
                          has_params=False,
                          tm_type=tm_type, length=element.l, delta_length=delta_l, params=params)

    def get_params(self, energy: float = 0.):
        return self.create_tm_param_func()

    def kick(self, X, kn):
        p = -kn[0] * X[5] + 0j
        for n in range(1, len(kn)):
            p += kn[n] * (X[0] + 1j * X[2]) ** n / factorial(n)
        X[1] = X[1] - np.real(p)
        X[3] = X[3] + np.imag(p)
        X[4] = X[4] - kn[0] * X[0]
        return X

    def map_function(self, X, energy: float):
        return self.kick(X, self.get_params().kn)
