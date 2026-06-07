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
    Public wrapper for beamline elements with cached transformations and TM selection.

    Four-Layer Architecture
    =======================
    OpticElement is the public wrapper in a four-layer design:

    1. **Wrapper** (OpticElement subclass): Public API, caching, framework behavior
    2. **Atom** (Element or Magnet): Physics state and hook methods
    3. **TMParams**: Typed data passed from atom to transformation
    4. **Transformation**: Tracking algorithm (TransferMap, SecondTM, etc.)

    This class owns framework state and forwards physics attributes to the atom.

    Two Transformation Paths
    ========================
    **first_order_tms** (always available)
        Linear R-matrix path built with TransferMap.
        Used by Twiss, optics, and linear analysis.
        Always present even if active method is different.

    **tms** (active tracking path)
        May use TransferMap, SecondTM, CavityTM, or family-specific transformation.
        Selected via default_tm and supported_tms.

    Tracking Method Selection
    =========================
    - ``default_tm``: Family fallback when tm=None
    - ``supported_tms``: Wrapper-selectable active methods (if declared)
    - Explicit requests must be in supported_tms (raises otherwise)
    - Global lattice requests may warn and fall back to default_tm

    See Also
    --------
    https://ocelot-collab.github.io/docs/docu/elements/architecture/
    https://ocelot-collab.github.io/docs/docu/elements/optical-element/

    Examples
    --------
    >>> from ocelot.cpbd.elements import Quadrupole
    >>> from ocelot.cpbd.transformations.second_order import SecondTM
    >>> q = Quadrupole(l=0.4, k1=1.2)
    >>> q.set_tm(SecondTM)  # Switch to second-order tracking
    >>> q.R(energy=1.0)  # Get R-matrix
    """

    default_tm = TransferMap
    supported_tms = None

    __is_init = False  # needed to disable __getattr__ and __setattr__ until __init__ is executed

    def __init__(self, element: Element, tm: Type[Transformation] | None = None,
                 default_tm: Type[Transformation] | None = None, **params) -> None:
        """
        Create a wrapper around an element atom.

        Parameters
        ----------
        element
            Physics implementation object.
        tm
            Active transformation requested for tracking. ``None`` means "use
            the wrapper class default".
        default_tm
            Family default transformation used when a broad global lattice
            request asks for an undeclared TM. If omitted, the wrapper class
            attribute ``default_tm`` is used.
        """
        # Physics state lives on `self.element`; framework/cache state stays on
        # the wrapper. This split lets us invalidate maps only when a physics
        # parameter changes.
        self.element = element
        self.default_tm = default_tm if default_tm is not None else type(self).default_tm
        self._validate_tm_declarations()
        # First-order maps are always available because optics/Twiss code uses
        # them even when the active tracking method is something else. For
        # ``has_edge=True`` families this also means the atom must always
        # provide the first-order entrance/exit hooks.
        try:
            self._first_order_tms = self._create_tms(self.element, TransferMap)
        except (AttributeError, NotImplementedError) as exc:
            raise self._tm_contract_error(TransferMap, first_order_only=True) from exc
        atom_tm_params = getattr(self.element, "params", {})
        params = {**atom_tm_params, **params}
        self._kwargs = params  # Storing transforamtion sp
        requested_tm = self.default_tm if tm is None else self._normalize_tm_request(
            tm, request_source="explicit", stacklevel=5
        )
        self._activate_tm(requested_tm, **params)
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

        Important:
        This property returns transformation objects, not already materialized
        matrices. The real energy-dependent params are created later via
        ``tm.get_params(energy)``.

        For example, to get the rotated R matrix for each active TM in a
        sequence, code must do the same energy propagation that ``R(energy)``
        does:

        ``params = tm.get_params(E)``
        ``R = params.get_rotated_R()``
        ``E += tm.get_delta_e()``

        In most cases callers should prefer the wrapper helpers ``R(energy)``,
        ``B(energy)``, and ``T(energy)`` instead of reimplementing that loop.
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
        ``TransferMap``-based optics path while still exposing only
        ``MultipoleTM`` as an active tracking method.

        Important:
        This property returns first-order ``Transformation`` objects, not ready
        R matrices. The actual first-order params are created on demand via
        ``tm.get_params(energy)``.

        For one TM object the materialized linear matrix is typically read as
        ``tm.get_params(E).get_rotated_R()``. For a whole element sequence the
        incoming reference energy must be advanced after each map using
        ``tm.get_delta_e()``. The wrapper helpers ``R(energy)`` and
        ``B(energy)`` already perform that bookkeeping for the optics path.
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
        except (AttributeError, NotImplementedError) as exc:
            if self._is_declared_tm(tm) or tm == self.default_tm:
                raise self._tm_contract_error(tm) from exc
            raise RuntimeError(
                f"{self.__class__.__name__} does not declare support for {tm.__name__}. "
                "The request should have been normalized before TM construction."
            ) from exc

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

    def _warn_global_tm_fallback(self, tm: Type[Transformation], stacklevel: int) -> None:
        """Warn when a global lattice TM request falls back to the family default."""
        warnings.warn(
            f"{self.__class__.__name__} does not declare support for {tm.__name__}; "
            f"global lattice request falls back to default {self.default_tm.__name__}.",
            stacklevel=stacklevel,
        )

    @staticmethod
    def _validate_request_source(request_source: str) -> None:
        """Validate the public request source contract."""
        if request_source not in {"explicit", "global"}:
            raise RuntimeError(
                f"Unsupported request_source={request_source!r}. Expected 'explicit' or 'global'."
            )

    def _normalize_tm_request(self, tm: Type[Transformation], request_source: str, stacklevel: int) -> Type[Transformation]:
        """Normalize TM requests according to the wrapper policy and request source."""
        self._validate_request_source(request_source)
        if tm == self.default_tm or self._is_declared_tm(tm):
            return tm
        if request_source == "global":
            self._warn_global_tm_fallback(tm, stacklevel=stacklevel)
            return self.default_tm
        raise RuntimeError(
            f"{self.__class__.__name__} does not declare support for {tm.__name__}. "
            f"Explicit requests must use supported_tms={self.supported_tms} or default_tm={self.default_tm.__name__}."
        )

    @staticmethod
    def _tm_hook_family(tm: Type[Transformation]) -> str:
        """Return the hook family the TM expects from the atom."""
        if tm is TransferMap:
            return "create_first_order_*_params(...)"
        if tm is SecondTM:
            return "create_second_order_*_params(...)"
        tm_name = tm.__name__
        if "RungeKutta" in tm_name:
            return "create_runge_kutta_*_params(...)"
        if "Kick" in tm_name:
            return "create_kick_*_params(...)"
        if "Cavity" in tm_name:
            return "create_cavity_tm_*_params(...)"
        if "Multipole" in tm_name:
            return "create_multipole_tm_main_params(...)"
        if "UndulatorTest" in tm_name:
            return "create_undulator_test_tm_main_params(...)"
        return "create_*_params(...)"

    def _tm_contract_error(self, tm: Type[Transformation], first_order_only: bool = False) -> RuntimeError:
        """Explain which hook family is missing for the requested TM build."""
        hook_family = self._tm_hook_family(tm)
        if first_order_only and self.element.has_edge:
            return RuntimeError(
                f"{self.__class__.__name__} has has_edge=True, so its atom must implement "
                f"{hook_family} for ENTRANCE, MAIN, and EXIT. "
                "CPBD could not build the required first_order_tms optics path."
            )
        if first_order_only:
            return RuntimeError(
                f"{self.__class__.__name__} could not build its required first_order_tms optics path. "
                f"Its atom must implement {hook_family}."
            )
        if self.element.has_edge:
            return RuntimeError(
                f"{self.__class__.__name__} declares support for {tm.__name__} and has has_edge=True, "
                f"so its atom must implement {hook_family} for ENTRANCE, MAIN, and EXIT."
            )
        return RuntimeError(
            f"{self.__class__.__name__} declares support for {tm.__name__}, "
            f"but its atom does not implement the required {hook_family} hooks."
        )

    def _activate_tm(self, tm: Type[Transformation], **params) -> None:
        """Build and store the active TM sequence for the normalized request."""
        if tm == TransferMap:
            # Reuse the cached first-order maps instead of rebuilding them.
            self._tms = self.first_order_tms
            self._tm_class_type = TransferMap
            return
        self.__init_tms(tm, **params)

    def apply(self, X: np.ndarray, energy: float):
        """Apply all active transformations to a particle coordinate array."""
        for tm in self.tms:
            tm.map_function(X, energy)

    def set_tm(self, tm: Transformation, request_source: str = "explicit", **params):
        """
        Set the active tracking method for the wrapper.

        Declared support is expected to be buildable. Explicit undeclared
        requests are treated as errors, while ``request_source='global'``
        allows a warning and fallback to ``default_tm``.
        """
        requested_tm = self._normalize_tm_request(tm, request_source=request_source, stacklevel=4)
        new_kwargs = params if params and params != self._kwargs else None
        if requested_tm != self._tm_class_type or new_kwargs:
            if new_kwargs:
                self._activate_tm(requested_tm, **new_kwargs)
                self._kwargs = new_kwargs
            else:
                self._activate_tm(requested_tm, **self._kwargs)
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

        main_tm_class = self.get_tm(TMTypes.MAIN, first_order_only).__class__

        if (start_l + delta_l > total_length or np.isclose(start_l + delta_l, total_length)):
            delta_l_red = total_length - start_l
            tm_list.append(main_tm_class.from_element(element=self.element, tm_type=TMTypes.MAIN, delta_l=delta_l_red))
            if self.element.has_edge and not ignore_edges:
                tm = self.get_tm(TMTypes.EXIT, first_order_only)
                tm_list.append(copy(tm))
        else:
            tm_list.append(main_tm_class.from_element(element=self.element, tm_type=TMTypes.MAIN, delta_l=delta_l))
        # tms.append(TMTypes.ROT_EXIT)
        return tm_list

    def get_tm(self, tm_type: TMTypes, first_order_only=False):
        tms = self.first_order_tms if first_order_only else self.tms
        for tm in tms:
            if tm.tm_type == tm_type:
                return tm

    @staticmethod
    def _create_tms(element: Element, tm: Type[Transformation], **params) -> List[Transformation]:
        # ``has_edge`` means the TM sequence is a three-map contract:
        # ENTRANCE -> MAIN -> EXIT. Any TM used on such a family must therefore
        # be able to bind both edge hooks and the main hook.
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
