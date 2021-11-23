import enum
import logging
from sys import flags
from abc import ABC, abstractmethod

import numpy as np

from ocelot.cpbd.beam import Particle, ParticleArray
from ocelot.cpbd.elements.element import Element

_logger = logging.getLogger(__name__)


class TMTypes(enum.Enum):
    ROT_ENTRANCE = 1
    ROT_EXIT = 2
    ENTRANCE = 3
    EXIT = 4
    MAIN = 5


class Transformation(ABC):
    """
    Base class for all Transforamtions.
    """

    def __init__(self, create_tm_param_func, delta_e_func, tm_type: TMTypes, length: float, delta_length: float = None) -> None:
        """[summary]

        Args:
            create_tm_param_func ([type]): [description]
            length (float): total length to which the transformation is applied.
        """
        self.create_tm_param_func = create_tm_param_func
        self.length = length if tm_type == TMTypes.MAIN else 0.0  # entrance/exit functions (e.g Edges or Couplerkicks) do not have a length.
        self.delta_length = delta_length
        self._delta_e_func = delta_e_func if tm_type == TMTypes.MAIN else None  # entrance/exit functions (e.g Edges or Couplerkicks) do not change beam energy.
        self.tm_type = tm_type

        self._current_energy = None  # used for caching tm params
        self._params = None

        self._map = None

    def _clean_cashed_values(self):
        self._current_energy = None  # used for caching tm params
        self._params = None

    def get_delta_e(self):
        if self._delta_e_func and self.tm_type == TMTypes.MAIN:  # Entrance and Exit tms (Edge TMs) don't have delta_e
            return self._delta_e_func(delta_length=self.delta_length, total_length=self.length)
        else:
            return 0.0

    @classmethod
    @abstractmethod
    def from_element(cls, element: Element, tm_type: TMTypes = TMTypes.MAIN, delta_l: float = None, **params):
        """[summary]
        This classmethod is used to create a new transforamtion from a element. 
        With 'params' custom parameter can be added to the Transformation. 
        :param element: Element for which the transformation is calculated. The Element have to implement the hook of concrete transformation.
        :type element: Element
        :param tm_type: Type of Transformation can be TMTypes.ENTRANCE, TMTypes.MAIN or TMTypes.EXIT, defaults to TMTypes.MAIN
        :type tm_type: TMTypes, optional
        :param delta_l: If the parameter is set, just a section of the element is taken into account for the tm calculation, defaults to None
        :type delta_l: float, optional
        :raises NotImplementedError: If not implemented
        """
        raise NotImplementedError

    @classmethod
    def create(cls, main_tm_params_func, delta_e_func, length, delta_length=None, entrance_tm_params_func=None,
               exit_tm_params_func=None, tm_type: TMTypes = TMTypes.MAIN, **params):
        """[summary]
        Factory method the concrete transforamtion. 
        :param main_tm_params_func: Function that is called on the element to calculate the transformation parameter for the main transforamtion.
        :type main_tm_params_func: [type]
        :param delta_e_func: Function that is called on the element to calculate the delta energy.
        :type delta_e_func: [type]
        :param length: Length of the element.
        :type length: [type]
        :param delta_length: define the delta length for which the transforamtion of the element will be calculated if delta_length is None length will be used, default is None.
        :type delta_length: [type], optional
        :param entrance_tm_params_func: : Function that is called on the element to calculate the transformation parameter for the entrance transforamtion, defaults to None
        :type entrance_tm_params_func: [type], optional
        :param exit_tm_params_func: Function that is called on the element to calculate the transformation parameter for the exit transforamtion, defaults to None
        :type exit_tm_params_func: [type], optional
        :param tm_type: [description], defaults to TMTypes.MAIN
        :type tm_type: TMTypes, optional
        :raises NotImplementedError: If the element doesn't implement the call back functions.
        :return: Instance of the concrete transfroamt
        :rtype: [type]
        """
        try:
            if tm_type == TMTypes.ENTRANCE:
                tm_params_func = entrance_tm_params_func
            elif tm_type == TMTypes.MAIN:
                tm_params_func = main_tm_params_func
            elif tm_type == TMTypes.EXIT:
                tm_params_func = exit_tm_params_func
            if not tm_params_func:
                raise NotImplementedError(f"{'entrance' if tm_type == TMTypes.ENTRANCE else 'exit'} function is not set in {cls.__class__.__name__}'s __init__")
        except AttributeError:
            raise NotImplementedError(f"The specific element have to implement the function {tm_params_func.__name__}.")
        if delta_length is not None and (delta_length > length):
            _logger.warning("delta_l > length of element. Set delta_l == length of element.")
            delta_length = length
        return cls(create_tm_param_func=tm_params_func, delta_e_func=delta_e_func, tm_type=tm_type, length=length, delta_length=delta_length, **params)

    def get_params(self, energy: float):
        """[summary]
        Calculates the parameters for the transformation
        :param energy: [description]
        :type energy: float
        :return: Depends on the type of transformation
        :rtype: [type]
        """
        if not self._params or self._current_energy != energy:
            self._params = self.create_tm_param_func(energy, self.delta_length)
            self._current_energy = energy
        return self._params

    def apply(self, prcl_series):
        """
        :param prcl_series: can be list of Particles [Particle_1, Particle_2, ... ] or ParticleArray
        :return: None
        """
        if prcl_series.__class__ == ParticleArray:
            self.map_function(prcl_series.rparticles, energy=prcl_series.E)
            #self.map(prcl_series.rparticles, energy=prcl_series.E)
            prcl_series.E += self.get_delta_e()
            #prcl_series.E += self.delta_e
            prcl_series.s += self.delta_length if self.delta_length != None else self.length

        elif prcl_series.__class__ == Particle:
            p = prcl_series
            p.x, p.px, p.y, p.py, p.tau, p.p = self.map_function(np.array([[p.x], [p.px], [p.y], [p.py], [p.tau], [p.p]]), p.E)[
                :, 0]
            p.s += self.delta_length if self.delta_length != None else self.length
            p.E += self.get_delta_e()

        elif prcl_series.__class__ == list and prcl_series[0].__class__ == Particle:
            # If the energy is not the same (p.E) for all Particles in the list of Particles
            # in that case cycle is applied. For particles with the same energy p.E
            list_e = np.array([p.E for p in prcl_series])
            if False in (list_e[:] == list_e[0]):
                for p in prcl_series:
                    self.map(np.array([[p.x], [p.px], [p.y], [p.py], [p.tau], [p.p]]), energy=p.E)

                    p.E += self.get_delta_e()
                    p.s += self.delta_length if self.delta_length != None else self.length
            else:
                pa = ParticleArray()
                pa.list2array(prcl_series)
                pa.E = prcl_series[0].E
                #self.map(pa.rparticles, energy=pa.E)
                self.map_function(pa.rparticles, energy=pa.E)
                pa.E += self.get_delta_e()
                pa.s += self.delta_length if self.delta_length != None else self.length
                pa.array2ex_list(prcl_series)

        else:
            _logger.error(" TransferMap.apply(): Unknown type of Particle_series: " + str(prcl_series.__class__.__name))
            raise Exception(
                " TransferMap.apply(): Unknown type of Particle_series: " + str(prcl_series.__class__.__name))

    @abstractmethod
    def map_function(self, X: np.ndarray, energy: float) -> np.ndarray:
        """[summary]
        This function calculate the transformation. It has to be overloaded by each transformation class.
        :param X: Particle array
        :type X: np.ndarray
        :param energy: Energy for which the transformation is calculated. 
        :type energy: float
        """
        raise NotImplementedError()
