import logging

import numpy as np

from ocelot.cpbd.transformations.transformation import Transformation, TMTypes
from ocelot.cpbd.elements.element import Element

_logger = logging.getLogger(__name__)


class TransferMap(Transformation):
    """[summary]
    Implementation of the first oder transformation.
    The concrete element atom have to implement: 
    create_first_order_main_params(self, energy: float, delta_length: float) -> FirstOrderPrams
    If the element has edges is also have to implement:
    create_first_order_entrance_params(self, energy: float, delta_length: float) -> FirstOrderPrams
    create_first_order_exit_params(self, energy: float, delta_length: float) -> FirstOrderPrams
    """
    
    def __init__(self, create_tm_param_func, delta_e_func, tm_type: TMTypes, length: float, delta_length: float = 0.0) -> None:
        super().__init__(create_tm_param_func, delta_e_func, tm_type, length, delta_length)

    @classmethod
    def from_element(cls, element: Element, tm_type: TMTypes = TMTypes.MAIN, delta_l=None):
        return cls.create(entrance_tm_params_func=element.create_first_order_entrance_params if element.has_edge else None,
                          delta_e_func=element.create_delta_e,
                          main_tm_params_func=element.create_first_order_main_params,
                          exit_tm_params_func=element.create_first_order_exit_params if element.has_edge else None,
                          tm_type=tm_type, length=element.l, delta_length=delta_l)

    def map_function(self, X, energy: float):
        """
        This function calculate the map function which can be overload if the map function is different to the first order mapping  
        @param X: Input Particle 
        @param energy: Initial energy
        """
        return self.mul_p_array(X, energy=energy)

    def mul_p_array(self, rparticles, energy=0.):
        """
        Calculates new rpaticles with the first order transformation, overrides the old rpaticles and returns the new rparticles. 
        :param rparticles: Can be a ParticleArray, Particle or a list of Particle object. 
        :param engery:
        :return: Returns the modified rparticles 
        """
        params = self.get_params(energy)
        a = np.add(np.dot(params.get_rotated_R(), rparticles), params.B)
        rparticles[:] = a[:]
        return rparticles

    def multiply_with_tm(self, tm: 'TransferMap', length) -> 'TransferMap':
        return TransferMap(create_tm_param_func=lambda energy: self.get_params(energy) * tm.get_params(energy),
                           length=length + tm.length)

    def __mul__(self, m):
        """
        :param m: TransferMap, Particle or Twiss
        :return: TransferMap, Particle or Twiss
        Ma = {Ba, Ra, Ta}
        Mb = {Bb, Rb, Tb}
        X1 = R*(X0 - dX) + dX = R*X0 + B
        B = (E - R)*dX
        """
        try:
            return m.multiply_with_tm(self, self.delta_length if self.delta_length is not None else self.length)
        except AttributeError as e:
            _logger.error(
                " TransferMap.__mul__: unknown object in transfer map multiplication: " + str(m.__class__.__name__))
            raise Exception(
                " TransferMap.__mul__: unknown object in transfer map multiplication: " + str(m.__class__.__name__))
