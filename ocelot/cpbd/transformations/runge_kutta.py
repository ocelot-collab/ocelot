from copy import copy

from ocelot.cpbd.high_order import rk_field
from ocelot.cpbd.transformations.transfer_map import TransferMap, TMTypes
from ocelot.cpbd.elements.element import Element


class RungeKuttaTM(TransferMap):
    """[summary]
    Implementation of the Runge Kutta transformation.
    The concrete element atom have to implement: 
    create_runga_kutta_main_params(self) -> RungeKuttaParams
    If the element has edges is also have to implement:
    create_runga_kutta_entrance_params(self) -> RungeKuttaParams
    create_runga_kutta_exit_params(self) -> RungeKuttaParams
    """

    def __init__(self, create_tm_param_func, delta_e_func, tm_type: TMTypes, length: float, delta_length: float = 0.0, **params) -> None:
        super().__init__(create_tm_param_func, delta_e_func, tm_type, length, delta_length)
        s_start = params.get('s_start')
        self.s_start = s_start if s_start else 0
        npoints = params.get('npoints')
        self.npoints = npoints if npoints else 200
        self.long_dynamics = True

    @classmethod
    def from_element(cls, element: Element, tm_type: TMTypes = TMTypes.MAIN, delta_l=None, **params):
        """[summary]
        :param element: An Element which implement at least 'create_runga_kutta_main_params(self) -> RungeKuttaParams'. If the element has edges 
        'create_runga_kutta_entrance_params(self) -> RungeKuttaParams' and 'create_runga_kutta_exit_params(self) -> RungeKuttaParams' have to be implemented as well.
        :type element: Element
        :param tm_type: Type of Transformation can be TMTypes.ENTRANCE, TMTypes.MAIN or TMTypes.EXIT, defaults to TMTypes.MAIN
        :type tm_type: TMTypes, optional
        :param delta_l: If the parameter is set, just a section of the element is taken into account for the tm calculation, defaults to None
        :type delta_l: [type], optional
        :return: [description]
        :rtype: [type]
        """
        return cls.create(entrance_tm_params_func=element.create_runge_kutta_entrance_params if element.has_edge else None,
                          delta_e_func=element.create_delta_e,
                          main_tm_params_func=element.create_runge_kutta_main_params,
                          exit_tm_params_func=element.create_runge_kutta_exit_params if element.has_edge else None,
                          tm_type=tm_type, length=element.l, delta_length=delta_l, **params)

    def get_params(self, energy):
        return self.create_tm_param_func(energy)

    def map_function(self, X, energy: float):
        params = self.get_params(energy)
        return rk_field(X, self.s_start, self.delta_length if self.delta_length != None else self.length, self.npoints, energy, params.mag_field, self.long_dynamics)
