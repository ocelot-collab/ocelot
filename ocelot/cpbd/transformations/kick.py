import numpy as np

from ocelot.cpbd.high_order import m_e_GeV
from ocelot.cpbd.tm_utils import transform_vec_ent, transform_vec_ext
from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.transformations.transformation import Transformation, TMTypes
from ocelot.cpbd.elements.element import Element


class KickTM(Transformation):
    """[summary]
    Implementation of the Kick Transforamtion.
    The concrete element atom have to implement at least 'create_kick_main_params(self) -> KickParams'. If the element has edges 
        'create_kick_entrance_params(self) -> KickParams' and 'create_kick_exit_params(self) -> KickParams' have to be implemented as well.
    """

    def __init__(self, create_tm_param_func, delta_e_func, tm_type: TMTypes, length: float, delta_length: float = 0.0, **params) -> None:    
        super().__init__(create_tm_param_func, delta_e_func, tm_type, length, delta_length)
        nkick = params.get('nkick')
        self.nkick = nkick if nkick else 1

    @classmethod
    def from_element(cls, element: Element, tm_type: TMTypes = TMTypes.MAIN, delta_l=None, **params):
        """[summary]
        :param element: An Element which implements at least 'create_kick_main_params(self) -> KickParams'. If the element has edges 
        'create_kick_entrance_params(self) -> KickParams' and 'create_kick_exit_params(self) -> KickParams' have to be implemented as well.
        :type element: Element
        :param tm_type: Type of Transformation can be TMTypes.ENTRANCE, TMTypes.MAIN or TMTypes.EXIT, defaults to TMTypes.MAIN
        :type tm_type: TMTypes, optional
        :param delta_l: If the parameter is set, just a section of the element is taken into account for the tm calculation, defaults to None
        :type delta_l: [type], optional
        :return: [description]
        :rtype: [type]
        """
        return cls.create(entrance_tm_params_func=element.create_kick_entrance_params if element.has_edge else None,
                          delta_e_func=element.create_delta_e,
                          main_tm_params_func=element.create_kick_main_params,
                          exit_tm_params_func=element.create_kick_exit_params if element.has_edge else None,
                          tm_type=tm_type, length=element.l, delta_length=delta_l, **params)

    def get_params(self, energy: float = 0.):
        return self.create_tm_param_func()

    def kick(self, X, l, angle, k1, k2, k3, energy, nkick=1):
        """
        does not work for dipole
        """
        gamma = energy / m_e_GeV
        coef = 0.
        beta = 1.
        if gamma != 0:
            gamma2 = gamma * gamma
            beta = 1. - 0.5 / gamma2
            coef = 1. / (beta * beta * gamma2)
        l = l / nkick
        angle = angle / nkick

        dl = l / 2.
        k1 = k1 * l
        k2 = k2 * l / 2.
        k3 = k3 * l / 6.

        for i in range(nkick):
            x = X[0] + X[1] * dl
            y = X[2] + X[3] * dl

            p = -angle * X[5] + 0j
            xy1 = x + 1j * y
            xy2 = xy1 * xy1
            xy3 = xy2 * xy1
            p += k1 * xy1 + k2 * xy2 + k3 * xy3
            X[1] = X[1] - np.real(p)
            X[3] = X[3] + np.imag(p)
            X[4] = X[4] + np.real(angle * xy1) / beta - X[5] * dl * coef

            X[0] = x + X[1] * dl
            X[2] = y + X[3] * dl
            # X[4] -= X[5] * dl * coef
        return X

    def kick_apply(self, X, energy):
        params = self.get_params()
        angle = params.angle
        k1 = params.k1
        k2 = params.k2
        k3 = params.k3
        tilt = params.tilt
        
        dx = params.dx
        dy = params.dy

        nkick = self.nkick
        l = self.delta_length if self.delta_length is not None else self.length

        if dx != 0 or dy != 0 or tilt != 0:
            X = transform_vec_ent(X, dx, dy, tilt)
        self.kick(X, l, angle, k1, k2, k3, energy, nkick=nkick)
        if dx != 0 or dy != 0 or tilt != 0:
            X = transform_vec_ext(X, dx, dy, tilt)

        return X

    def map_function(self, X, energy: float):
        return self.kick_apply(X, energy)
