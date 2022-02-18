from copy import copy
from ocelot.cpbd.transformations.transfer_map import TransferMap
from typing import List, Dict, Type, Callable


import numpy as np

from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.transformations import second_order
from ocelot.cpbd.transformations.transformation import Transformation, TMTypes
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


class OpticElement:
    """[summary]
    Facade between old interface an new interface.
    """

    __is_init = False  # needed to disable __getattr__ and __setattr__ until __init__ is executed

    def __init__(self, element: Element, tm: Type[Transformation], default_tm: Type[Transformation], **params) -> None:
        """[summary]
        Creates a optic element which holds element atom and its transformation.
        Each concrete optic element have to implement its own __init__
        it the concrete element parameters.

        :param element: Concrete element atom.
        :type element: Element
        :param tm: Transformation that is used by the element.
        :type tm: Type[Transformation]
        :param default_tm: If tm is not supported by the specific element the default_tm is used.
        :type default_tm: Type[Transformation]
        """
        self.element = element
        self.default_tm = default_tm
        # every optics element has a first_order tm to calculate e.g Twiss Parameters
        self._first_order_tms = self._create_tms(self.element, TransferMap)
        self._kwargs = params  # Storing transforamtion sp
        if tm == TransferMap:
            self._tms = self._first_order_tms
            self._tm_class_type = tm
        else:
            self.__init_tms(tm, **params)
        self.__is_init = True  # needed to disable __getattr__ and __setattr__ in __init__ phase. Do not add new attributes after.

    # To access all getter attributes of element to fulfill old iface

    def __getattr__(self, name):
        if self.__is_init and name in self.element.__dict__:
            return getattr(self.element, name)
        raise AttributeError(f"{self.__class__.__name__} object has no attribute {name}")

    # To access all setter attributes of element to fulfill old iface
    def __setattr__(self, name, value):
        if self.__is_init and name in self.element.__dict__:

            if self._tms is not None:
                for tm in self._tms:
                    tm._clean_cashed_values()
            if self._first_order_tms is not None:
                for tm in self._first_order_tms:
                    tm._clean_cashed_values()
            self._first_order_tms = None
            self._tms = None
            return setattr(self.element, name, value)
        return object.__setattr__(self, name, value)

    @property
    def tms(self) -> List[Transformation]:
        """[summary]
        Returns a list of transformation function [f(X, energy),...].
        :return: List of first order transformations
        :rtype: List[Callable[[np.ndarray, float], np.ndarray]]
        """
        if self._tms is None:
            if self._tm_class_type == TransferMap:
                if self._first_order_tms is None:
                    self._first_order_tms = self._create_tms(self.element, TransferMap)
                self._tms = self._first_order_tms
            else:
                self._tms = self._create_tms(self.element, self._tm_class_type)
        return self._tms

    @property
    def first_order_tms(self) -> List[Transformation]:
        """[summary]
        Returns a list of first order transformation function [f(X, energy),...].
        :return: List of first order transformations
        :rtype: List[Callable[[np.ndarray, float], np.ndarray]]
        """
        if self._first_order_tms is None:
            self._first_order_tms = self._create_tms(self.element, TransferMap)
        return self._first_order_tms

    def B(self, energy: float) -> List[np.ndarray]:
        """[summary]
        Calculates B matrices for second order transformation if second order transformation is set otherwise
        B matrices are calculated for first order transformation.

        :param energy: [description]
        :type energy: float
        :return: List of B matrices
        :rtype: List[np.ndarray]
        """
        tms = self._tms if self._tm_class_type == SecondTM else self.first_order_tms
        res = []
        E = energy
        for tm in tms:
            res.append(tm.get_params(E).B)
            E += tm.get_delta_e()
        return res

    def R(self, energy: float) -> List[np.ndarray]:
        """[summary]
        Calculates R matrices for second order transformation if second order transformation is set otherwise
        the B matrices are calculated for first order transformation.

        :param energy: [description]
        :type energy: float
        :return: [description]
        :rtype: List[np.ndarray]
        """
        tms = self.tms if self._tm_class_type == SecondTM else self.first_order_tms
        res = []
        E = energy
        for tm in tms:
            res.append(tm.get_params(E).get_rotated_R())
            E += tm.get_delta_e()
        return res

    def T(self, energy: float) -> List[np.ndarray]:
        """[summary]
        Calculates the T matrices for second order transformation if second order transformation is set
        otherwise a zero matrices will be returned.

        :param energy: [description]
        :type energy: float
        :return: List of B matrices
        :rtype: List[np.ndarray]
        """
        if self._tm_class_type != SecondTM:
            return [np.zeros((6, 6, 6)) for _ in self.first_order_tms]
        res = []
        E = energy
        for tm in self.tms:
            res.append(tm.get_params(E).get_rotated_T())
            E += tm.get_delta_e()
        return res

    def __init_tms(self, tm: Transformation, **params):
        """[summary]
        Initialize the transformations. If the transformation is not supported default_tm will be used.

        :param tm: Transformation to be set.
        :type tm: Transformation
        """
        try:
            self._tms = self._create_tms(self.element, tm, **params)
            self._tm_class_type = tm
        except AttributeError as e:
            print(f"Can't set {tm.__name__} for {self.__class__.__name__} fall back to default tm which is {self.default_tm.__name__}.")
            self._tms = self._create_tms(self.element, self.default_tm, **params)
            self._tm_class_type = self.default_tm

    def apply(self, X: np.ndarray, energy: float):
        """[summary]
        Apply all transformations on particle array.

        :param X: Array of particles
        :type X: np.ndarray
        :param energy: Given energy
        :type energy: float
        """
        for tm in self.tms:
            tm.map_function(X, energy)

    def set_tm(self, tm: Transformation, **params):
        """[summary]
        Sets a new transformation for the element. Transformation specific parameter can be added via **params
        e.g elem.set_tm(tm=RungaKuttraTM, npoints=2000),
        Note: Transformation specific parameters can be also set via __init__.

        :param tm: 
        :type tm: Transformation
        """
        new_kwargs = params if params and params != self._kwargs else None
        if tm != self._tm_class_type or new_kwargs:
            if tm == self.default_tm:
                self._tms = self._first_order_tms
            else:
                if new_kwargs:
                    self.__init_tms(tm, **new_kwargs)
                    self._kwargs = new_kwargs
                else:
                    self.__init_tms(tm, **self._kwargs)

            self._tm_class_type = tm
            for tm in self.tms:
                tm._clean_cashed_values()

    def get_section_tms(self, delta_l: float, start_l: float = 0.0, ignore_edges=False, first_order_only=False) -> List[Transformation]:
        """[summary]
        Calculates transformations for a section of the Element. The section is defined by start_l and delta_l.

        :param delta_l: The length of the section.
        :type delta_l: float
        :param start_l: Start position in the element, defaults to 0.0
        :type start_l: float, optional
        :param ignore_edges: Ingore the Entrance- and Exist-Transformations if True, defaults to False
        :type ignore_edges: bool, optional
        :param first_order_only: Returns the first order transforamtions if True, defaults to False
        :type first_order_only: bool, optional
        :return: A List of Transformation for the section. 
        :rtype: List[Transformation]
        """
        #tms = [TMTypes.ROT_ENTRANCE]
        tm_list = []
        total_length = self.element.l
        if start_l < 1e-10:
            if self.element.has_edge and not ignore_edges:
                tm = self.get_tm(TMTypes.ENTRANCE, first_order_only)
                tm_list.append(copy(tm))
            if np.isclose(delta_l, total_length):
                tm = self.get_tm(TMTypes.MAIN, first_order_only)
                tm_list.append(copy(tm))
                if self.element.has_edge and not ignore_edges:
                    tm = self.get_tm(TMTypes.EXIT, first_order_only)
                    tm_list.append(copy(tm))
                return tm_list

        if (start_l + delta_l > total_length or np.isclose(start_l + delta_l, total_length)):
            delta_l_red = total_length - start_l
            TM_Class = self.get_tm(TMTypes.MAIN).__class__
            tm_list.append(TM_Class.from_element(element=self.element, tm_type=TMTypes.MAIN, delta_l=delta_l_red))
            if self.element.has_edge and not ignore_edges:
                tm = self.get_tm(TMTypes.EXIT, first_order_only)
                tm_list.append(copy(tm))
        else:
            TM_Class = self.get_tm(TMTypes.MAIN, first_order_only).__class__
            tm_list.append(TM_Class.from_element(element=self.element, tm_type=TMTypes.MAIN, delta_l=delta_l))
        # tms.append(TMTypes.ROT_EXIT)
        return tm_list

    def get_tm(self, tm_type: TMTypes, first_order_only=False):
        tms = self.first_order_tms if first_order_only else self.tms
        for tm in tms:
            if tm.tm_type == tm_type:
                return tm

    @staticmethod
    def _create_tms(element: Element, tm: Type[Transformation], **params) -> List[Transformation]:
        tms = []
        if element.has_edge:
            tms.append(tm.from_element(element, TMTypes.ENTRANCE, **params))
            tms.append(tm.from_element(element, TMTypes.MAIN, **params))
            tms.append(tm.from_element(element, TMTypes.EXIT, **params))
        else:
            tms.append(tm.from_element(element, TMTypes.MAIN, **params))
        return tms

    def __str__(self):
        return self.element.__str__()

    def __repr__(self):
        return f"<{type(self).__name__}: name={self.id} at {hex(id(self))}>"
