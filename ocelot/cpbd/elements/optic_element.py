from copy import copy
from typing import List, Type
import warnings

import numpy as np

from ocelot.cpbd.elements.element import Element
from ocelot.cpbd.transformations.transformation import Transformation, TMTypes
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


class OpticElement:
    """
    Public facade around an element atom.

    Current CPBD contract:
    - the wrapper owns framework state such as cached transformations, the
      selected tracking method, and section slicing helpers
    - the wrapped ``element`` object owns the physics state and the
      ``create_*_params(...)`` hook methods used by transformations
    - most public element attributes are forwarded to ``self.element`` to
      preserve the historical user-facing API

    Two TM roles must be kept separate:

    - ``first_order_tms`` is the always-available ``TransferMap`` cache used by
      linear optics / Twiss code.
    - ``tms`` is the active tracking method selected on the wrapper.

    ``default_tm`` remains the family fallback for active tracking. If a
    wrapper also declares ``supported_tms``, those entries should be read as
    wrapper-selectable active tracking methods, not as a list of every
    internal first-order path that exists for optics. Generic wrappers still
    keep the legacy hook-based fallback path for undeclared requests, while
    pinned families may override ``set_tm`` more strictly.
    """

    default_tm = TransferMap
    supported_tms = None

    __is_init = False  # needed to disable __getattr__ and __setattr__ until __init__ is executed

    def __init__(self, element: Element, tm: Type[Transformation], default_tm: Type[Transformation], **params) -> None:
        """
        Create a wrapper around an element atom.

        Parameters
        ----------
        element
            Physics implementation object.
        tm
            Active transformation requested for tracking.
        default_tm
            Family default transformation used when the requested one is not
            supported or when a wrapper intentionally pins the family to one TM.
        """
        # Physics state lives on `self.element`; framework/cache state stays on
        # the wrapper. This split lets us invalidate maps only when a physics
        # parameter changes.
        self.element = element
        self.default_tm = default_tm
        self._validate_tm_declarations()
        # First-order maps are always available because optics/Twiss code uses
        # them even when the active tracking method is something else.
        self._first_order_tms = self._create_tms(self.element, TransferMap)
        self._kwargs = params  # Storing transforamtion sp
        self._activate_tm(tm, **params)
        self.__is_init = True  # needed to disable __getattr__ and __setattr__ in __init__ phase. Do not add new attributes after.

    def __getattr__(self, name):
        # Forward reads to the wrapped atom to preserve the historical public
        # element API without storing duplicate physics state on the wrapper.
        if self.__is_init and name in self.element.__dict__:
            return getattr(self.element, name)
        raise AttributeError(f"{self.__class__.__name__} object has no attribute {name}")

    def __setattr__(self, name, value):
        if self.__is_init and name in self.element.__dict__:
            # Only changes to physics parameters trigger TM cache invalidation.
            # Internal wrapper bookkeeping (e.g. `_tms`) is stored on the
            # wrapper itself and should not recursively rebuild maps.
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
        """
        Active transformations for tracking.

        These may differ from ``first_order_tms`` when a family uses a custom
        tracking method such as ``CavityTM`` or ``RungeKuttaTM``.
        """
        if self._tms is None:
            if self._tm_class_type == TransferMap:
                if self._first_order_tms is None:
                    self._first_order_tms = self._create_tms(self.element, TransferMap)
                self._tms = self._first_order_tms
            else:
                # Rebuild the active tracking maps lazily from the current
                # physics state on the wrapped element.
                self._tms = self._create_tms(self.element, self._tm_class_type)
        return self._tms

    @property
    def first_order_tms(self) -> List[Transformation]:
        """
        First-order transformations used by optics/Twiss code.

        These are always kept available even when the active tracking method is
        something else. This is why a family such as ``Multipole`` can keep a
        ``TransferMap``-based optics path while still pinning active tracking
        to ``MultipoleTM``.
        """
        if self._first_order_tms is None:
            self._first_order_tms = self._create_tms(self.element, TransferMap)
        return self._first_order_tms

    def B(self, energy: float) -> List[np.ndarray]:
        """
        Return the sequence of B vectors at the requested energy.

        For ``SecondTM`` the active second-order maps are used. Otherwise CPBD
        continues to use the first-order maps for this accessor.
        """
        tms = self._tms if self._tm_class_type == SecondTM else self.first_order_tms
        res = []
        E = energy
        for tm in tms:
            res.append(tm.get_params(E).B)
            E += tm.get_delta_e()
        return res

    def R(self, energy: float) -> List[np.ndarray]:
        """
        Return the sequence of rotated R matrices at the requested energy.

        For ``SecondTM`` the active second-order maps are used. Otherwise CPBD
        continues to expose the first-order matrices.
        """
        tms = self.tms if self._tm_class_type == SecondTM else self.first_order_tms
        res = []
        E = energy
        for tm in tms:
            res.append(tm.get_params(E).get_rotated_R())
            E += tm.get_delta_e()
        return res

    def T(self, energy: float) -> List[np.ndarray]:
        """
        Return the sequence of rotated T tensors.

        Families not using ``SecondTM`` expose zero tensors here by current
        CPBD convention.
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
        """
        Initialize active transformations.

        Requests are normalized through the declared wrapper contract first.
        If a wrapper declares support for a TM, the underlying atom is expected
        to provide the required hook methods.
        """
        try:
            self._tms = self._create_tms(self.element, tm, **params)
            self._tm_class_type = tm
        except AttributeError as exc:
            if self._is_declared_tm(tm) or tm == self.default_tm:
                raise RuntimeError(
                    f"{self.__class__.__name__} declares support for {tm.__name__}, "
                    "but its atom does not implement the required create_*_params hooks."
                ) from exc
            warnings.warn(
                f"Can't set {tm.__name__} for {self.__class__.__name__}; "
                f"falling back to default {self.default_tm.__name__}.",
                stacklevel=3,
            )
            self._activate_tm(self.default_tm, **params)

    def _validate_tm_declarations(self) -> None:
        """Check that the family default is part of the declared TM contract."""
        if self.supported_tms is not None and self.default_tm not in self.supported_tms:
            raise RuntimeError(
                f"{self.__class__.__name__} declares default_tm={self.default_tm.__name__}, "
                "but default_tm is missing from supported_tms."
            )

    def _is_declared_tm(self, tm: Type[Transformation]) -> bool:
        """Return True when the wrapper explicitly declares the TM as supported."""
        return self.supported_tms is not None and tm in self.supported_tms

    def _warn_if_undeclared_tm(self, tm: Type[Transformation]) -> None:
        """Warn when a generic wrapper is asked for an undeclared TM."""
        if not self._is_declared_tm(tm):
            warnings.warn(
                f"{self.__class__.__name__} does not declare support for {tm.__name__}; "
                "trying legacy hook-based fallback.",
                stacklevel=3,
            )

    def _activate_tm(self, tm: Type[Transformation], **params) -> None:
        """Build and store the active TM sequence for the normalized request."""
        if tm == TransferMap:
            # Reuse the cached first-order maps instead of rebuilding them.
            self._tms = self.first_order_tms
            self._tm_class_type = TransferMap
            return
        self._warn_if_undeclared_tm(tm)
        self.__init_tms(tm, **params)

    def apply(self, X: np.ndarray, energy: float):
        """Apply all active transformations to a particle coordinate array."""
        for tm in self.tms:
            tm.map_function(X, energy)

    def set_tm(self, tm: Transformation, **params):
        """
        Set the active tracking method for the wrapper.

        Declared support is expected to be buildable. Undeclared requests keep
        the legacy generic fallback path unless a concrete wrapper pins the
        family more strictly.
        """
        new_kwargs = params if params and params != self._kwargs else None
        if tm != self._tm_class_type or new_kwargs:
            if new_kwargs:
                self._activate_tm(tm, **new_kwargs)
                self._kwargs = new_kwargs
            else:
                self._activate_tm(tm, **self._kwargs)
            for tm in self.tms:
                tm._clean_cashed_values()

    def get_section_tms(self, delta_l: float, start_l: float = 0.0, ignore_edges=False, first_order_only=False) -> List[Transformation]:
        """
        Build transformations for a section of the element.

        Current edge semantics are intentionally preserved:
        - a slice starting at 0 includes ``ENTRANCE``
        - a slice ending at ``l`` includes ``EXIT``
        - the main body is rebuilt for the requested ``delta_l``
        - copied edge maps are not rescaled
        """
        # Build only the maps needed for a slice through the element. For edge
        # elements this preserves the current CPBD convention:
        # ENTRANCE -> MAIN(slice) -> EXIT.
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

    def get_transfer_geometry(self):
        return self.element.get_transfer_geometry()
