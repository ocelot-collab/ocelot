from math import factorial

import numpy as np

from ocelot.cpbd.transformations.transfer_map import TransferMap, TMTypes
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.tm_params.multipole_params import thin_dipole_dispersion
from ocelot.cpbd.tm_utils import transform_vec_ent, transform_vec_ext


class MultipoleTM(TransferMap):
    """[summary]
    Implementation of the Multipole Transformation.
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

    def kick(self, X, kn, energy: float = 0.):
        dispersion = thin_dipole_dispersion(kn[0], energy)
        p = -dispersion * X[5] + 0j
        for n in range(1, len(kn)):
            p += kn[n] * (X[0] + 1j * X[2]) ** n / factorial(n)
        X[1] = X[1] - np.real(p)
        X[3] = X[3] + np.imag(p)
        X[4] += dispersion.real * X[0] - dispersion.imag * X[2]
        return X

    def map_function(self, X, energy: float):
        params = self.get_params()
        dx = getattr(params, "dx", 0.)
        dy = getattr(params, "dy", 0.)
        tilt = getattr(params, "tilt", 0.)
        if not (dx or dy or tilt):
            return self.kick(X, params.kn, energy)
        transform_vec_ent(X, dx, dy, tilt)
        self.kick(X, params.kn, energy)
        return transform_vec_ext(X, dx, dy, tilt)
