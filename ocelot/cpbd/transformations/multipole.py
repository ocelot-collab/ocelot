from math import factorial

import numpy as np

from ocelot.cpbd.transformations.transfer_map import TransferMap, TMTypes
from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.tm_params.multipole_params import thin_dipole_dispersion
from ocelot.cpbd.tm_utils import transform_vec_ent, transform_vec_ext


class MultipoleTM(TransferMap):
    """Apply the zero-length kick of a :class:`~ocelot.cpbd.elements.multipole.Multipole`.

    In magnet coordinates let ``z = x + i*y``. An integrated coefficient of
    order n ≥ 1 contributes ``kn[n] * z**n / n!`` to the kick polynomial.
    Its real part is subtracted from px and its imaginary part is added to py.
    For a real ``kn[0] = theta``, the bend adds ``theta*delta/beta0`` to px
    and ``theta*x/beta0`` to tau, where ``delta = (E-E0)/(p0*c)`` and
    ``beta0 = v0/c``. Positions and kicks use the displaced, rolled magnet
    axes before the particle returns to the beam coordinates.
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
