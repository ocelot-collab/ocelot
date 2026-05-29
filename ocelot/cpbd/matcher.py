"""Object-oriented optics matching for Ocelot.

The API follows a "vary + targets" style similar to Xsuite, but keeps Ocelot
physics objects and naming.

Example
-------
```python
problem = MatchProblem(lat, tw0, periodic=False)
problem.vary_element(q1, quantity="k1", limits=(-10, 10), name="Q1")
problem.vary_linked_elements([q2, q3], scales=[1.0, -1.0], quantity="k1", name="PS_Q23")
problem.target_twiss(end, "beta_x", value=12.0, tol=1e-3)
problem.target_twiss_delta(m1, m2, "muy", value=3*np.pi/2, wrap_phase=True)
problem.target_rmatrix(m1, m2, i=0, j=1, value=8.0)
result = problem.solve(solver="ls_trf", max_iter=300)
print(result.success, result.merit)
```
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import OptimizeResult, least_squares, minimize

from ocelot.cpbd.beam import Twiss
from ocelot.cpbd.beam_params import radiation_integrals
from ocelot.cpbd.elements.optic_element import OpticElement


def _as_float(value: Any, name: str) -> float:
    try:
        return float(value)
    except Exception as exc:
        raise ValueError(f"Expected numeric value for '{name}', got {value!r}") from exc


def _residual_scalar(value: float, target: float, relation: str, tol: float) -> float:
    relation = relation.strip().lower()
    if relation in {"==", "eq"}:
        diff = value - target
        adiff = abs(diff)
        if adiff <= tol:
            return 0.0
        return np.sign(diff) * (adiff - tol)
    if relation in {"<", "<=", "le"}:
        return max(0.0, value - (target + tol))
    if relation in {">", ">=", "ge"}:
        return min(0.0, value - (target - tol))
    if relation in {"abs<", "a<"}:
        return max(0.0, abs(value) - (target + tol))
    if relation in {"abs>", "a>"}:
        return min(0.0, abs(value) - (target - tol))
    raise ValueError(f"Unsupported relation '{relation}'")


def _wrap_to_pi(value: float) -> float:
    return (value + np.pi) % (2.0 * np.pi) - np.pi


@dataclass
class Vary:
    """One optimization variable."""

    name: str
    getter: Callable[[], float]
    setter: Callable[[float], None]
    limits: Optional[Tuple[Optional[float], Optional[float]]] = None
    step: Optional[float] = None
    weight: float = 1.0
    max_step: Optional[float] = None
    active: bool = True
    tag: str = ""

    def get(self) -> float:
        return float(self.getter())

    def set(self, value: float) -> None:
        self.setter(float(value))

    @property
    def lower(self) -> float:
        if self.limits is None or self.limits[0] is None:
            return -np.inf
        return float(self.limits[0])

    @property
    def upper(self) -> float:
        if self.limits is None or self.limits[1] is None:
            return np.inf
        return float(self.limits[1])


class PowerSupplyVary(Vary):
    """Single knob applied to many elements (e.g. one power supply)."""

    def __init__(
        self,
        name: str,
        elements: Sequence[Any],
        quantity: str = "k1",
        scales: Optional[Sequence[float]] = None,
        limits: Optional[Tuple[Optional[float], Optional[float]]] = None,
        step: Optional[float] = None,
        weight: float = 1.0,
        max_step: Optional[float] = None,
        active: bool = True,
        tag: str = "",
    ):
        if len(elements) == 0:
            raise ValueError("PowerSupplyVary requires at least one element")
        if scales is None:
            scales = [1.0] * len(elements)
        if len(scales) != len(elements):
            raise ValueError("Length of 'scales' must match 'elements'")
        if any(s == 0 for s in scales):
            raise ValueError("PowerSupplyVary scale cannot be zero")

        channels: List[Tuple[Any, str, float]] = [
            (elements[i], quantity, float(scales[i])) for i in range(len(elements))
        ]

        def _getter() -> float:
            elem0, attr0, scale0 = channels[0]
            return float(getattr(elem0, attr0)) / scale0

        def _setter(x: float) -> None:
            for elem, a, s in channels:
                setattr(elem, a, float(x) * s)

        super().__init__(
            name=name,
            getter=_getter,
            setter=_setter,
            limits=limits,
            step=step,
            weight=weight,
            max_step=max_step,
            active=active,
            tag=tag,
        )
        self.channels = channels


@dataclass
class MatchState:
    """Snapshot of optics state for one matcher evaluation.

    A new ``MatchState`` is built on every call to :meth:`MatchProblem.evaluate`
    and on every solver residual evaluation during :meth:`MatchProblem.solve`.
    It is not persistent across iterations.

    Attributes
    ----------
    lat:
        Lattice instance used for this evaluation.
    twiss_start:
        Effective initial Twiss used in this evaluation.
        - If ``periodic=False``, this is a copy of ``problem.twiss0``.
        - If ``periodic=True``, this is the periodic Twiss found from the
          current lattice/variables.
    twiss_by_element:
        Dictionary mapping element object -> Twiss at that element (after
        applying its first-order transfer maps).
    twiss_sequence:
        Ordered list of Twiss objects along ``lat.sequence``.
        Item ``i`` corresponds to ``lat.sequence[i]``.
    twiss_end:
        Twiss at end of matched line for this evaluation
        (normally the last entry in ``twiss_sequence``).
    failed, failure_reason:
        Evaluation status and error message if optics propagation failed.
    _r_cache:
        Internal lazy cache for transfer matrices between element pairs.
        Key: ``(start_elem, end_elem)``.
        Value: corresponding 6x6 ``R`` matrix.
    """

    lat: Any
    twiss_start: Twiss
    twiss_by_element: Dict[OpticElement, Twiss]
    twiss_sequence: List[Twiss]
    twiss_end: Twiss
    failed: bool = False
    failure_reason: str = ""
    _r_cache: Dict[Tuple[OpticElement, OpticElement], np.ndarray] = field(default_factory=dict)

    def twiss_at(self, element: OpticElement) -> Twiss:
        """Return Twiss at a given element for this evaluation state."""

        if element not in self.twiss_by_element:
            raise ValueError(f"Element '{getattr(element, 'id', element)}' not found in twiss map")
        return self.twiss_by_element[element]

    def r_matrix(self, start: OpticElement, end: OpticElement) -> np.ndarray:
        """Return cached (or lazily computed) transfer matrix between elements.

        This avoids repeated ``lat.transfer_maps`` calls when several targets
        or objectives request the same ``R`` block in one evaluation.
        """

        key = (start, end)
        if key not in self._r_cache:
            _b, r_mat, _t = self.lat.transfer_maps(energy=self.twiss_start.E, start=start, stop=end)
            self._r_cache[key] = r_mat
        return self._r_cache[key]


@dataclass
class TargetReport:
    name: str
    residual_norm: float
    met: bool
    weight: float
    details: Dict[str, Any] = field(default_factory=dict)


@dataclass
class ObjectiveReport:
    name: str
    residual_norm: float
    weight: float
    details: Dict[str, Any] = field(default_factory=dict)


class Objective:
    """Base objective class.

    Unlike ``Target``, objectives are generally used for minimization terms
    (e.g. radiation integrals, custom lattice metrics).
    """

    def __init__(
        self,
        name: Optional[str] = None,
        weight: float = 1.0,
        active: bool = True,
        tag: str = "",
    ):
        self.name = name or self.__class__.__name__
        self.weight = float(weight)
        self.active = bool(active)
        self.tag = tag

    def residuals(self, state: MatchState) -> np.ndarray:
        raise NotImplementedError

    def weighted_residuals(self, state: MatchState) -> np.ndarray:
        return np.sqrt(self.weight) * np.atleast_1d(self.residuals(state)).astype(float)

    def report(self, state: MatchState) -> ObjectiveReport:
        residuals = np.atleast_1d(self.residuals(state)).astype(float)
        return ObjectiveReport(
            name=self.name,
            residual_norm=float(np.linalg.norm(residuals)),
            weight=self.weight,
            details={},
        )


class FunctionObjective(Objective):
    """Generic objective based on a user-supplied callable.

    Parameters
    ----------
    func:
        Callable ``func(state)`` returning scalar or vector.
    mode:
        - ``"minimize"``: values are treated as residual-like quantities to
          minimize directly (for positive functions this means minimization
          towards zero).
        - ``"target"``: values are constrained to ``target`` using ``relation``.
        - ``"residual"``: returned values are interpreted as residuals directly.
    """

    def __init__(
        self,
        func: Callable[[MatchState], Any],
        mode: str = "minimize",
        target: float = 0.0,
        relation: str = "==",
        tol: float = 0.0,
        **kwargs: Any,
    ):
        super().__init__(**kwargs)
        self.func = func
        self.mode = mode.strip().lower()
        self.target = float(target)
        self.relation = relation
        self.tol = float(tol)
        if self.mode not in {"minimize", "target", "residual"}:
            raise ValueError("FunctionObjective mode must be one of: minimize, target, residual")

    def _values(self, state: MatchState) -> np.ndarray:
        val = self.func(state)
        arr = np.atleast_1d(np.asarray(val, dtype=float))
        if not np.all(np.isfinite(arr)):
            return np.array([1.0e6], dtype=float)
        return arr

    def residuals(self, state: MatchState) -> np.ndarray:
        arr = self._values(state)
        if self.mode == "residual":
            return arr
        if self.mode == "minimize":
            return arr - self.target
        # mode == target
        return np.asarray(
            [_residual_scalar(float(v), self.target, self.relation, self.tol) for v in arr],
            dtype=float,
        )

    def report(self, state: MatchState) -> ObjectiveReport:
        values = self._values(state)
        residuals = self.residuals(state)
        return ObjectiveReport(
            name=self.name,
            residual_norm=float(np.linalg.norm(residuals)),
            weight=self.weight,
            details={
                "type": "function_objective",
                "mode": self.mode,
                "values": values.tolist(),
                "target": self.target,
                "relation": self.relation,
            },
        )


class I5Objective(Objective):
    """Built-in objective minimizing radiation integral I5.

    This is implemented as a regular objective term and can be combined with
    any other custom objective.
    """

    def __init__(
        self,
        reject_unphysical_partition: bool = True,
        **kwargs: Any,
    ):
        super().__init__(**kwargs)
        self.reject_unphysical_partition = bool(reject_unphysical_partition)

    def _i5_data(self, state: MatchState) -> Tuple[float, float, float, float, float]:
        return radiation_integrals(state.lat, state.twiss_start, nsuperperiod=1)

    def residuals(self, state: MatchState) -> np.ndarray:
        try:
            _i1, i2, _i3, i4, i5 = self._i5_data(state)
            if i2 == 0:
                return np.array([1.0e6], dtype=float)
            if self.reject_unphysical_partition:
                je = 2.0 + i4 / i2
                jx = 1.0 - i4 / i2
                jy = 1.0
                if je < 0 or jx < 0 or jy < 0:
                    return np.array([1.0e6], dtype=float)
            return np.array([np.sqrt(max(0.0, float(i5)))], dtype=float)
        except Exception:
            return np.array([1.0e6], dtype=float)

    def report(self, state: MatchState) -> ObjectiveReport:
        residuals = self.residuals(state)
        details: Dict[str, Any] = {"type": "i5"}
        try:
            i1, i2, i3, i4, i5 = self._i5_data(state)
            details.update({"I1": float(i1), "I2": float(i2), "I3": float(i3), "I4": float(i4), "I5": float(i5)})
        except Exception as exc:
            details["error"] = str(exc)
        return ObjectiveReport(
            name=self.name,
            residual_norm=float(np.linalg.norm(residuals)),
            weight=self.weight,
            details=details,
        )


class Target:
    """Base target class."""

    def __init__(
        self,
        name: Optional[str] = None,
        weight: float = 1.0,
        tol: float = 0.0,
        active: bool = True,
        tag: str = "",
    ):
        self.name = name or self.__class__.__name__
        self.weight = float(weight)
        self.tol = float(tol)
        self.active = bool(active)
        self.tag = tag

    def residuals(self, state: MatchState) -> np.ndarray:
        raise NotImplementedError

    def report(self, state: MatchState) -> TargetReport:
        residuals = np.atleast_1d(self.residuals(state)).astype(float)
        rnorm = float(np.linalg.norm(residuals))
        return TargetReport(
            name=self.name,
            residual_norm=rnorm,
            met=bool(np.all(np.abs(residuals) < 1.0e-12)),
            weight=self.weight,
            details={},
        )

    def weighted_residuals(self, state: MatchState) -> np.ndarray:
        return np.sqrt(self.weight) * np.atleast_1d(self.residuals(state)).astype(float)


class TwissTarget(Target):
    """Target one Twiss quantity at one element."""

    def __init__(
        self,
        element: OpticElement,
        quantity: str,
        value: float,
        relation: str = "==",
        **kwargs: Any,
    ):
        super().__init__(**kwargs)
        self.element = element
        self.quantity = quantity
        self.value = _as_float(value, f"{self.quantity}.value")
        self.relation = relation

    def residuals(self, state: MatchState) -> np.ndarray:
        tw = state.twiss_at(self.element)
        actual = float(getattr(tw, self.quantity))
        r = _residual_scalar(actual, self.value, self.relation, self.tol)
        return np.array([r], dtype=float)

    def report(self, state: MatchState) -> TargetReport:
        tw = state.twiss_at(self.element)
        actual = float(getattr(tw, self.quantity))
        r = float(self.residuals(state)[0])
        return TargetReport(
            name=self.name,
            residual_norm=abs(r),
            met=abs(r) < 1.0e-12,
            weight=self.weight,
            details={
                "type": "twiss",
                "element": getattr(self.element, "id", None),
                "quantity": self.quantity,
                "actual": actual,
                "target": self.value,
                "relation": self.relation,
            },
        )


class TwissDifferenceTarget(Target):
    """Target difference of one Twiss quantity between two elements.

    For phase quantities (``mux``, ``muy``), ``wrap_phase=True`` can be used
    to treat differences modulo ``2*pi`` for equality constraints.
    """

    def __init__(
        self,
        start: OpticElement,
        end: OpticElement,
        quantity: str,
        value: float,
        relation: str = "==",
        wrap_phase: bool = False,
        **kwargs: Any,
    ):
        super().__init__(**kwargs)
        self.start = start
        self.end = end
        self.quantity = quantity
        self.value = _as_float(value, f"{self.quantity}.value")
        self.relation = relation
        self.wrap_phase = bool(wrap_phase)

    def _actual(self, state: MatchState) -> float:
        v0 = float(getattr(state.twiss_at(self.start), self.quantity))
        v1 = float(getattr(state.twiss_at(self.end), self.quantity))
        return v1 - v0

    def residuals(self, state: MatchState) -> np.ndarray:
        actual = self._actual(state)
        if self.wrap_phase and self.relation.strip().lower() in {"==", "eq"}:
            raw = actual - self.value
            r = _wrap_to_pi(raw)
            if abs(r) <= self.tol:
                r = 0.0
        else:
            r = _residual_scalar(actual, self.value, self.relation, self.tol)
        return np.array([r], dtype=float)

    def report(self, state: MatchState) -> TargetReport:
        actual = self._actual(state)
        r = float(self.residuals(state)[0])
        return TargetReport(
            name=self.name,
            residual_norm=abs(r),
            met=abs(r) < 1.0e-12,
            weight=self.weight,
            details={
                "type": "twiss_delta",
                "start": getattr(self.start, "id", None),
                "end": getattr(self.end, "id", None),
                "quantity": self.quantity,
                "actual": actual,
                "target": self.value,
                "relation": self.relation,
                "wrap_phase": self.wrap_phase,
            },
        )


class TwissPeriodicityTarget(Target):
    """Target equality of one Twiss quantity between start and end.

    If ``start`` or ``end`` are not provided, the problem's effective initial
    Twiss or final Twiss are used respectively. This allows partial periodic
    constraints without enabling ``MatchProblem(periodic=True)``.
    """

    def __init__(
        self,
        quantity: str,
        start: Optional[OpticElement] = None,
        end: Optional[OpticElement] = None,
        relation: str = "==",
        wrap_phase: bool = False,
        **kwargs: Any,
    ):
        super().__init__(**kwargs)
        self.quantity = quantity
        self.start = start
        self.end = end
        self.relation = relation
        self.wrap_phase = bool(wrap_phase)

    def _start_twiss(self, state: MatchState) -> Twiss:
        if self.start is None:
            return state.twiss_start
        return state.twiss_at(self.start)

    def _end_twiss(self, state: MatchState) -> Twiss:
        if self.end is None:
            return state.twiss_end
        return state.twiss_at(self.end)

    def _actual(self, state: MatchState) -> float:
        v0 = float(getattr(self._start_twiss(state), self.quantity))
        v1 = float(getattr(self._end_twiss(state), self.quantity))
        return v1 - v0

    def residuals(self, state: MatchState) -> np.ndarray:
        actual = self._actual(state)
        if self.wrap_phase and self.relation.strip().lower() in {"==", "eq"}:
            r = _wrap_to_pi(actual)
            if abs(r) <= self.tol:
                r = 0.0
        else:
            r = _residual_scalar(actual, 0.0, self.relation, self.tol)
        return np.array([r], dtype=float)

    def report(self, state: MatchState) -> TargetReport:
        start_tw = self._start_twiss(state)
        end_tw = self._end_twiss(state)
        actual_start = float(getattr(start_tw, self.quantity))
        actual_end = float(getattr(end_tw, self.quantity))
        r = float(self.residuals(state)[0])
        return TargetReport(
            name=self.name,
            residual_norm=abs(r),
            met=abs(r) < 1.0e-12,
            weight=self.weight,
            details={
                "type": "twiss_periodic",
                "start": getattr(self.start, "id", "twiss_start"),
                "end": getattr(self.end, "id", "twiss_end"),
                "quantity": self.quantity,
                "start_value": actual_start,
                "end_value": actual_end,
                "actual": actual_end - actual_start,
                "target": 0.0,
                "relation": self.relation,
                "wrap_phase": self.wrap_phase,
            },
        )


class GlobalTwissTarget(Target):
    """Enforce inequality/equality over all lattice elements."""

    def __init__(
        self,
        quantity: str,
        value: float,
        relation: str = "<=",
        **kwargs: Any,
    ):
        super().__init__(**kwargs)
        self.quantity = quantity
        self.value = _as_float(value, f"{self.quantity}.value")
        self.relation = relation

    def residuals(self, state: MatchState) -> np.ndarray:
        out = []
        for tw in state.twiss_sequence:
            val = float(getattr(tw, self.quantity))
            out.append(_residual_scalar(val, self.value, self.relation, self.tol))
        return np.asarray(out, dtype=float)

    def report(self, state: MatchState) -> TargetReport:
        vals = np.asarray([float(getattr(tw, self.quantity)) for tw in state.twiss_sequence], dtype=float)
        residuals = self.residuals(state)
        return TargetReport(
            name=self.name,
            residual_norm=float(np.linalg.norm(residuals)),
            met=bool(np.all(np.abs(residuals) < 1.0e-12)),
            weight=self.weight,
            details={
                "type": "global_twiss",
                "quantity": self.quantity,
                "target": self.value,
                "relation": self.relation,
                "min": float(np.min(vals)) if len(vals) else np.nan,
                "max": float(np.max(vals)) if len(vals) else np.nan,
            },
        )


class RMatrixElementTarget(Target):
    """Target one R[i,j] entry between two elements."""

    def __init__(
        self,
        start: OpticElement,
        end: OpticElement,
        i: int,
        j: int,
        value: float,
        relation: str = "==",
        **kwargs: Any,
    ):
        super().__init__(**kwargs)
        self.start = start
        self.end = end
        self.i = int(i)
        self.j = int(j)
        self.value = _as_float(value, "r_matrix.value")
        self.relation = relation
        if not (0 <= self.i < 6 and 0 <= self.j < 6):
            raise ValueError(f"R-matrix indices out of range: ({self.i}, {self.j})")

    def residuals(self, state: MatchState) -> np.ndarray:
        r_mat = state.r_matrix(self.start, self.end)
        actual = float(r_mat[self.i, self.j])
        r = _residual_scalar(actual, self.value, self.relation, self.tol)
        return np.array([r], dtype=float)

    def report(self, state: MatchState) -> TargetReport:
        r_mat = state.r_matrix(self.start, self.end)
        actual = float(r_mat[self.i, self.j])
        r = float(self.residuals(state)[0])
        return TargetReport(
            name=self.name,
            residual_norm=abs(r),
            met=abs(r) < 1.0e-12,
            weight=self.weight,
            details={
                "type": "r_matrix",
                "start": getattr(self.start, "id", None),
                "end": getattr(self.end, "id", None),
                "index": (self.i, self.j),
                "actual": actual,
                "target": self.value,
                "relation": self.relation,
            },
        )


class RMatrixBlockTarget(Target):
    """Target a block of R-matrix entries."""

    def __init__(
        self,
        start: OpticElement,
        end: OpticElement,
        target_matrix: np.ndarray,
        rows: Optional[Sequence[int]] = None,
        cols: Optional[Sequence[int]] = None,
        **kwargs: Any,
    ):
        super().__init__(**kwargs)
        self.start = start
        self.end = end
        self.target_matrix = np.asarray(target_matrix, dtype=float)

        if rows is None and cols is None:
            if self.target_matrix.shape == (2, 2):
                rows = [0, 1]
                cols = [0, 1]
            elif self.target_matrix.shape == (4, 4):
                rows = [0, 1, 2, 3]
                cols = [0, 1, 2, 3]
            elif self.target_matrix.shape == (6, 6):
                rows = [0, 1, 2, 3, 4, 5]
                cols = [0, 1, 2, 3, 4, 5]
            else:
                raise ValueError("Implicit rows/cols only supported for (2,2), (4,4), and (6,6) blocks")
        elif rows is None or cols is None:
            raise ValueError("Both rows and cols must be provided for custom R-matrix block")

        self.rows = [int(r) for r in rows]  # type: ignore[arg-type]
        self.cols = [int(c) for c in cols]  # type: ignore[arg-type]
        if self.target_matrix.shape != (len(self.rows), len(self.cols)):
            raise ValueError(
                f"target_matrix shape {self.target_matrix.shape} does not match rows/cols "
                f"{(len(self.rows), len(self.cols))}"
            )

    def residuals(self, state: MatchState) -> np.ndarray:
        r_mat = state.r_matrix(self.start, self.end)
        actual = r_mat[np.ix_(self.rows, self.cols)]
        diff = actual - self.target_matrix
        return diff.ravel().astype(float)

    def report(self, state: MatchState) -> TargetReport:
        residuals = self.residuals(state)
        return TargetReport(
            name=self.name,
            residual_norm=float(np.linalg.norm(residuals)),
            met=bool(np.all(np.abs(residuals) < 1.0e-12)),
            weight=self.weight,
            details={
                "type": "r_matrix_block",
                "start": getattr(self.start, "id", None),
                "end": getattr(self.end, "id", None),
                "shape": self.target_matrix.shape,
            },
        )


class TotalLengthTarget(Target):
    """Target total matched line length (end s)."""

    def __init__(self, value: float, relation: str = "==", **kwargs: Any):
        super().__init__(**kwargs)
        self.value = _as_float(value, "total_length.value")
        self.relation = relation

    def residuals(self, state: MatchState) -> np.ndarray:
        actual = float(state.twiss_end.s)
        r = _residual_scalar(actual, self.value, self.relation, self.tol)
        return np.array([r], dtype=float)

    def report(self, state: MatchState) -> TargetReport:
        actual = float(state.twiss_end.s)
        r = float(self.residuals(state)[0])
        return TargetReport(
            name=self.name,
            residual_norm=abs(r),
            met=abs(r) < 1.0e-12,
            weight=self.weight,
            details={
                "type": "total_length",
                "actual": actual,
                "target": self.value,
                "relation": self.relation,
            },
        )


@dataclass
class MatchResult:
    """Final solve output returned by :meth:`MatchProblem.solve`."""

    success: bool
    message: str
    solver: str
    merit: float
    nfev: int
    nit: Optional[int]
    variables: Dict[str, float]
    target_reports: List[TargetReport]
    objective_reports: List[ObjectiveReport]
    optimize_result: OptimizeResult


class MatchProblem:
    """Main matching problem container.

    A ``MatchProblem`` stores:

    - variables (`Vary`) that the optimizer can change
    - targets (`Target`) that encode constraints
    - objectives (`Objective`) that encode minimization terms

    During evaluation/solve it builds a ``MatchState`` from the current
    variable values, asks all active targets/objectives for residuals, applies
    weights, and minimizes the total squared residual norm.
    """

    def __init__(
        self,
        lat: Any,
        twiss0: Twiss,
        periodic: bool = False,
    ):
        self.lat = lat
        self.twiss0 = Twiss(twiss0)
        self.periodic = bool(periodic)

        self.variables: List[Vary] = []
        self.targets: List[Target] = []
        self.objectives: List[Objective] = []

    # ---- Variables ----
    def add_variable(self, variable: Vary) -> Vary:
        """Register a variable in the optimization problem."""

        self.variables.append(variable)
        return variable

    def vary_element(
        self,
        element: Any,
        quantity: str = "k1",
        name: Optional[str] = None,
        limits: Optional[Tuple[Optional[float], Optional[float]]] = None,
        step: Optional[float] = None,
        weight: float = 1.0,
        max_step: Optional[float] = None,
        active: bool = True,
        tag: str = "",
    ) -> Vary:
        """Add one variable controlling one element quantity.

        Example quantities are ``k1``, ``angle``, ``l``. For ``l`` a default
        lower bound of ``0`` is applied when limits are not provided.
        """

        if limits is None and quantity == "l":
            limits = (0.0, None)

        def _getter() -> float:
            return float(getattr(element, quantity))

        def _setter(x: float) -> None:
            setattr(element, quantity, float(x))

        var = Vary(
            name=name or f"{getattr(element, 'id', element)}.{quantity}",
            getter=_getter,
            setter=_setter,
            limits=limits,
            step=step,
            weight=weight,
            max_step=max_step,
            active=active,
            tag=tag,
        )
        return self.add_variable(var)

    def vary_twiss(
        self,
        quantity: str,
        name: Optional[str] = None,
        limits: Optional[Tuple[Optional[float], Optional[float]]] = None,
        step: Optional[float] = None,
        weight: float = 1.0,
        max_step: Optional[float] = None,
        active: bool = True,
        tag: str = "",
    ) -> Vary:
        """Add one variable controlling an initial Twiss quantity.

        Typical quantities: ``beta_x``, ``alpha_x``, ``beta_y``, ``alpha_y``,
        ``Dx``, ``Dxp``.
        """

        def _getter() -> float:
            return float(getattr(self.twiss0, quantity))

        def _setter(x: float) -> None:
            setattr(self.twiss0, quantity, float(x))

        var = Vary(
            name=name or f"twiss0.{quantity}",
            getter=_getter,
            setter=_setter,
            limits=limits,
            step=step,
            weight=weight,
            max_step=max_step,
            active=active,
            tag=tag,
        )
        return self.add_variable(var)

    def vary_linked_elements(
        self,
        elements: Sequence[Any],
        scales: Optional[Sequence[float]] = None,
        quantity: str = "k1",
        name: str = "linked_elements",
        limits: Optional[Tuple[Optional[float], Optional[float]]] = None,
        step: Optional[float] = None,
        weight: float = 1.0,
        max_step: Optional[float] = None,
        active: bool = True,
        tag: str = "",
    ) -> PowerSupplyVary:
        """Add one shared variable coupled to many elements.

        This models a common hardware channel (for example one power supply).
        Each element gets ``quantity = knob * scale``.
        """

        var = PowerSupplyVary(
            name=name,
            elements=elements,
            quantity=quantity,
            scales=scales,
            limits=limits,
            step=step,
            weight=weight,
            max_step=max_step,
            active=active,
            tag=tag,
        )
        self.variables.append(var)
        return var

    # ---- Targets ----
    def add_target(self, target: Target) -> Target:
        """Register a target (constraint term)."""

        self.targets.append(target)
        return target

    # ---- Objectives ----
    def add_objective(self, objective: Objective) -> Objective:
        """Register an objective (minimization term)."""

        self.objectives.append(objective)
        return objective

    def objective_function(
        self,
        func: Callable[[MatchState], Any],
        mode: str = "minimize",
        target: float = 0.0,
        relation: str = "==",
        tol: float = 0.0,
        name: Optional[str] = None,
        weight: float = 1.0,
        active: bool = True,
        tag: str = "",
    ) -> FunctionObjective:
        """Add a generic objective from a callable ``func(state)``.

        Modes:
        - ``minimize``: minimize ``func(state) - target``
        - ``target``: apply relation/tolerance to ``func(state)``
        - ``residual``: treat returned value directly as residual(s)
        """

        return self.add_objective(
            FunctionObjective(
                func=func,
                mode=mode,
                target=target,
                relation=relation,
                tol=tol,
                name=name or "function_objective",
                weight=weight,
                active=active,
                tag=tag,
            )
        )

    def minimize_function(
        self,
        func: Callable[[MatchState], Any],
        name: Optional[str] = None,
        weight: float = 1.0,
        target: float = 0.0,
        active: bool = True,
        tag: str = "",
    ) -> FunctionObjective:
        """Add a generic function minimization objective.

        ``func`` is called as ``func(state)`` on each evaluation, where
        ``state`` is :class:`MatchState`. It should return a scalar or vector
        quantity that is meaningful to push towards ``target`` (default 0).
        """
        return self.objective_function(
            func=func,
            mode="minimize",
            target=target,
            name=name or "minimize_function",
            weight=weight,
            active=active,
            tag=tag,
        )

    def minimize_i5_integral(
        self,
        name: str = "I5",
        weight: float = 1.0e14,
        reject_unphysical_partition: bool = True,
        active: bool = True,
        tag: str = "",
    ) -> I5Objective:
        """Add built-in radiation integral ``I5`` minimization objective.

        This is a convenience wrapper around :class:`I5Objective`.
        """

        return self.add_objective(
            I5Objective(
                name=name,
                weight=weight,
                active=active,
                tag=tag,
                reject_unphysical_partition=reject_unphysical_partition,
            )
        )

    # Friendly aliases for users coming from "constraint" terminology.
    def add_constraint(self, target: Target) -> Target:
        """Alias of :meth:`add_target`."""

        return self.add_target(target)

    def add_constrain(self, target: Target) -> Target:
        """Legacy spelling alias of :meth:`add_target`."""

        return self.add_target(target)

    def target_twiss(
        self,
        element: OpticElement,
        quantity: str,
        value: float,
        relation: str = "==",
        name: Optional[str] = None,
        weight: float = 1.0,
        tol: float = 0.0,
        active: bool = True,
        tag: str = "",
    ) -> TwissTarget:
        """Constrain one Twiss quantity at one element."""

        return self.add_target(
            TwissTarget(
                element=element,
                quantity=quantity,
                value=value,
                relation=relation,
                name=name or f"{getattr(element, 'id', element)}.{quantity}",
                weight=weight,
                tol=tol,
                active=active,
                tag=tag,
            )
        )

    def target_twiss_delta(
        self,
        start: OpticElement,
        end: OpticElement,
        quantity: str,
        value: float,
        relation: str = "==",
        wrap_phase: bool = False,
        name: Optional[str] = None,
        weight: float = 1.0,
        tol: float = 0.0,
        active: bool = True,
        tag: str = "",
    ) -> TwissDifferenceTarget:
        """Constrain Twiss quantity difference ``end - start``.

        Parameters
        ----------
        wrap_phase:
            Relevant for phase-like quantities (usually ``mux``/``muy``) when
            ``relation == "=="``. If ``True``, residual is wrapped to
            ``[-pi, pi]`` so phase-equivalent targets differing by ``2*pi`` are
            treated as equal.
        """

        return self.add_target(
            TwissDifferenceTarget(
                start=start,
                end=end,
                quantity=quantity,
                value=value,
                relation=relation,
                wrap_phase=wrap_phase,
                name=name or f"{quantity}@{getattr(end, 'id', end)}-{getattr(start, 'id', start)}",
                weight=weight,
                tol=tol,
                active=active,
                tag=tag,
            )
        )

    def target_periodic_twiss(
        self,
        quantity: str,
        start: Optional[OpticElement] = None,
        end: Optional[OpticElement] = None,
        relation: str = "==",
        wrap_phase: bool = False,
        name: Optional[str] = None,
        weight: float = 1.0,
        tol: float = 0.0,
        active: bool = True,
        tag: str = "",
    ) -> TwissPeriodicityTarget:
        """Constrain one Twiss quantity to be equal at start and end.

        This is useful for partial periodicity, for example matching
        ``beta_x`` at the beginning and end while independently targeting
        ``beta_y`` or ``alpha_y``. It does not compute a full periodic Twiss
        solution like ``MatchProblem(periodic=True)``; it adds an ordinary
        residual target ``end - start == 0``.

        If ``start`` is omitted, ``state.twiss_start`` is used. If ``end`` is
        omitted, ``state.twiss_end`` is used.
        """

        start_name = getattr(start, "id", "twiss_start")
        end_name = getattr(end, "id", "twiss_end")
        return self.add_target(
            TwissPeriodicityTarget(
                quantity=quantity,
                start=start,
                end=end,
                relation=relation,
                wrap_phase=wrap_phase,
                name=name or f"periodic.{quantity}@{start_name}->{end_name}",
                weight=weight,
                tol=tol,
                active=active,
                tag=tag,
            )
        )

    def target_global(
        self,
        quantity: str,
        value: float,
        relation: str = "<=",
        name: Optional[str] = None,
        weight: float = 1.0,
        tol: float = 0.0,
        active: bool = True,
        tag: str = "",
    ) -> GlobalTwissTarget:
        """Constrain a Twiss quantity over the full lattice sequence."""

        return self.add_target(
            GlobalTwissTarget(
                quantity=quantity,
                value=value,
                relation=relation,
                name=name or f"global.{quantity}",
                weight=weight,
                tol=tol,
                active=active,
                tag=tag,
            )
        )

    def target_rmatrix(
        self,
        start: OpticElement,
        end: OpticElement,
        i: int,
        j: int,
        value: float,
        relation: str = "==",
        name: Optional[str] = None,
        weight: float = 1.0,
        tol: float = 0.0,
        active: bool = True,
        tag: str = "",
    ) -> RMatrixElementTarget:
        """Constrain one transfer matrix entry ``R[i, j]`` between two elements."""

        return self.add_target(
            RMatrixElementTarget(
                start=start,
                end=end,
                i=i,
                j=j,
                value=value,
                relation=relation,
                name=name or f"R[{i},{j}]@{getattr(start, 'id', start)}->{getattr(end, 'id', end)}",
                weight=weight,
                tol=tol,
                active=active,
                tag=tag,
            )
        )

    def target_rmatrix_block(
        self,
        start: OpticElement,
        end: OpticElement,
        target_matrix: np.ndarray,
        rows: Optional[Sequence[int]] = None,
        cols: Optional[Sequence[int]] = None,
        name: Optional[str] = None,
        weight: float = 1.0,
        tol: float = 0.0,
        active: bool = True,
        tag: str = "",
    ) -> RMatrixBlockTarget:
        """Constrain a selected block of transfer matrix entries."""

        return self.add_target(
            RMatrixBlockTarget(
                start=start,
                end=end,
                target_matrix=target_matrix,
                rows=rows,
                cols=cols,
                name=name or f"Rblock@{getattr(start, 'id', start)}->{getattr(end, 'id', end)}",
                weight=weight,
                tol=tol,
                active=active,
                tag=tag,
            )
        )

    def target_total_length(
        self,
        value: float,
        relation: str = "==",
        name: str = "total_length",
        weight: float = 1.0,
        tol: float = 0.0,
        active: bool = True,
        tag: str = "",
    ) -> TotalLengthTarget:
        """Constrain total matched line length (end ``s``)."""

        return self.add_target(
            TotalLengthTarget(
                value=value,
                relation=relation,
                name=name,
                weight=weight,
                tol=tol,
                active=active,
                tag=tag,
            )
        )

    # ---- Core evaluation ----
    def _active_variables(self) -> List[Vary]:
        return [v for v in self.variables if v.active]

    def _active_targets(self) -> List[Target]:
        return [t for t in self.targets if t.active]

    def _active_objectives(self) -> List[Objective]:
        return [o for o in self.objectives if o.active]

    def _compute_state(self) -> MatchState:
        """Build optics state for current variable values.

        Twiss is propagated element-by-element through each element's first
        order transfer maps and cached in ``MatchState``.
        """

        tw_seed = Twiss(self.twiss0)
        if self.periodic:
            tw_periodic = self.lat.periodic_twiss(tw_seed)
            if tw_periodic is None:
                return MatchState(
                    lat=self.lat,
                    twiss_start=tw_seed,
                    twiss_by_element={},
                    twiss_sequence=[],
                    twiss_end=tw_seed,
                    failed=True,
                    failure_reason="No periodic Twiss solution",
                )
            tw_seed = Twiss(tw_periodic)

        tw = Twiss(tw_seed)
        tw.s = 0.0
        twiss_by_element: Dict[OpticElement, Twiss] = {}
        twiss_sequence: List[Twiss] = []

        try:
            for elem in self.lat.sequence:
                for tm in elem.first_order_tms:
                    tw = tm * tw
                twiss_by_element[elem] = tw
                twiss_sequence.append(tw)
            return MatchState(
                lat=self.lat,
                twiss_start=tw_seed,
                twiss_by_element=twiss_by_element,
                twiss_sequence=twiss_sequence,
                twiss_end=tw,
                failed=False,
                failure_reason="",
            )
        except Exception as exc:
            return MatchState(
                lat=self.lat,
                twiss_start=tw_seed,
                twiss_by_element=twiss_by_element,
                twiss_sequence=twiss_sequence,
                twiss_end=tw,
                failed=True,
                failure_reason=str(exc),
            )

    def _build_residual_vector(self, state: MatchState) -> np.ndarray:
        """Collect weighted residuals from all active targets/objectives."""

        if state.failed:
            return np.array([1.0e9], dtype=float)

        residuals: List[float] = []
        for target in self._active_targets():
            r = target.weighted_residuals(state)
            residuals.extend(r.tolist())

        for objective in self._active_objectives():
            r = objective.weighted_residuals(state)
            residuals.extend(r.tolist())

        if len(residuals) == 0:
            residuals = [0.0]

        arr = np.asarray(residuals, dtype=float)
        if not np.all(np.isfinite(arr)):
            return np.array([1.0e9], dtype=float)
        return arr

    def _build_target_reports(self, state: MatchState) -> List[TargetReport]:
        reports: List[TargetReport] = []
        for target in self._active_targets():
            try:
                reports.append(target.report(state))
            except Exception as exc:
                reports.append(
                    TargetReport(
                        name=target.name,
                        residual_norm=np.inf,
                        met=False,
                        weight=target.weight,
                        details={"error": str(exc)},
                    )
                )
        return reports

    def _build_objective_reports(self, state: MatchState) -> List[ObjectiveReport]:
        reports: List[ObjectiveReport] = []
        for objective in self._active_objectives():
            try:
                reports.append(objective.report(state))
            except Exception as exc:
                reports.append(
                    ObjectiveReport(
                        name=objective.name,
                        residual_norm=np.inf,
                        weight=objective.weight,
                        details={"error": str(exc)},
                    )
                )
        return reports

    def evaluate(self) -> Tuple[float, List[TargetReport], List[ObjectiveReport], MatchState]:
        """Evaluate current problem without running optimization.

        Returns
        -------
        merit:
            Sum of squares of weighted residuals.
        target_reports, objective_reports:
            Per-term diagnostic summaries.
        state:
            Current optics state used for evaluation.
        """

        state = self._compute_state()
        residuals = self._build_residual_vector(state)
        merit = float(np.dot(residuals, residuals))
        target_reports = self._build_target_reports(state)
        objective_reports = self._build_objective_reports(state)
        return merit, target_reports, objective_reports, state

    # ---- Solve ----
    def solve(
        self,
        solver: str = "ls_trf",
        max_iter: int = 300,
        tol: float = 1.0e-8,
        verbose: bool = False,
        restore_if_fail: bool = False,
    ) -> MatchResult:
        """Run numerical optimization and return the best found solution.

        Parameters
        ----------
        solver:
            Supported values:
            ``ls_trf``, ``ls_dogbox``, ``ls_lm``, ``least_squares``,
            ``simplex``, ``nelder-mead``, ``powell``, ``bfgs``, ``cg``,
            ``lbfgsb``, ``l-bfgs-b``, ``slsqp``.
            Bound constraints are not supported by ``ls_lm``, ``bfgs``, ``cg``.
            If any active variable has finite bounds and solver has no bound
            support, :class:`ValueError` is raised.
        max_iter:
            Maximum iterations/function evaluations depending on solver.
        tol:
            Scalar tolerance forwarded to the underlying SciPy solver.
        restore_if_fail:
            If ``True`` and solve fails, restore initial variable values.
        """

        active_vars = self._active_variables()
        if len(active_vars) == 0:
            raise ValueError("No active variables to optimize")
        if len(self._active_targets()) == 0 and len(self._active_objectives()) == 0:
            raise ValueError("No active targets or objectives")

        x0 = np.asarray([v.get() for v in active_vars], dtype=float)
        x_saved = x0.copy()
        lb = np.asarray([v.lower for v in active_vars], dtype=float)
        ub = np.asarray([v.upper for v in active_vars], dtype=float)
        has_bounds = np.any(np.isfinite(lb)) or np.any(np.isfinite(ub))

        # Give users explicit diagnostics before scipy raises a generic error.
        out_of_bounds = []
        eps = 1.0e-12
        for i, var in enumerate(active_vars):
            xi = float(x0[i])
            li = float(lb[i])
            ui = float(ub[i])
            if (np.isfinite(li) and xi < li - eps) or (np.isfinite(ui) and xi > ui + eps):
                out_of_bounds.append((var.name, xi, li, ui))

        if out_of_bounds:
            def _fmt_bound(v: float) -> str:
                if np.isneginf(v):
                    return "-inf"
                if np.isposinf(v):
                    return "+inf"
                return f"{v:.16g}"

            details = "\n".join(
                f"  - {name}: x0={x:.16g}, bounds=({_fmt_bound(lo)}, {_fmt_bound(hi)})"
                for name, x, lo, hi in out_of_bounds
            )
            raise ValueError(
                "Initial guess is outside of provided bounds for variable(s):\n"
                f"{details}"
            )

        diff_step = np.asarray(
            [v.step if v.step is not None else np.nan for v in active_vars],
            dtype=float,
        )
        if np.all(np.isnan(diff_step)):
            diff_step_opt = None
        else:
            diff_step = np.where(np.isnan(diff_step), 0.0, diff_step)
            diff_step_opt = diff_step

        def _apply_x(x: np.ndarray) -> None:
            for i, var in enumerate(active_vars):
                var.set(float(x[i]))

        def _residual_fun(x: np.ndarray) -> np.ndarray:
            _apply_x(x)
            state = self._compute_state()
            return self._build_residual_vector(state)

        solver_key = solver.strip().lower()
        if solver_key in {"ls_trf", "ls_dogbox", "ls_lm", "least_squares"}:
            ls_method = {
                "ls_trf": "trf",
                "ls_dogbox": "dogbox",
                "ls_lm": "lm",
                "least_squares": "trf",
            }[solver_key]

            if ls_method == "lm" and has_bounds:
                raise ValueError(
                    "least_squares method 'lm' does not support bounds. "
                    "Use one of: ls_trf, ls_dogbox, least_squares, "
                    "simplex/nelder-mead, powell, lbfgsb/l-bfgs-b, slsqp."
                )

            opt = least_squares(
                _residual_fun,
                x0,
                method=ls_method,
                bounds=(lb, ub) if ls_method != "lm" else (-np.inf, np.inf),
                max_nfev=int(max_iter),
                xtol=tol,
                ftol=tol,
                gtol=tol,
                diff_step=diff_step_opt,
                verbose=2 if verbose else 0,
            )
            x_best = np.asarray(opt.x, dtype=float)
            merit = float(np.dot(opt.fun, opt.fun))
            success = bool(opt.success)
            message = str(opt.message)
            nfev = int(opt.nfev)
            nit = int(getattr(opt, "nit", -1)) if hasattr(opt, "nit") else None
        else:
            method_map = {
                "simplex": "Nelder-Mead",
                "nelder-mead": "Nelder-Mead",
                "powell": "Powell",
                "bfgs": "BFGS",
                "cg": "CG",
                "lbfgsb": "L-BFGS-B",
                "l-bfgs-b": "L-BFGS-B",
                "slsqp": "SLSQP",
            }
            if solver_key not in method_map:
                raise ValueError(f"Unsupported solver '{solver}'")
            mname = method_map[solver_key]
            if has_bounds and mname not in {"L-BFGS-B", "SLSQP", "Powell", "Nelder-Mead"}:
                raise ValueError(
                    f"Solver '{solver}' does not support bounds. "
                    "Use one of: ls_trf, ls_dogbox, least_squares, "
                    "simplex/nelder-mead, powell, lbfgsb/l-bfgs-b, slsqp."
                )

            def _scalar_fun(x: np.ndarray) -> float:
                r = _residual_fun(x)
                return float(np.dot(r, r))

            options: Dict[str, Any] = {"maxiter": int(max_iter), "disp": bool(verbose)}
            if mname == "Nelder-Mead":
                options["maxfev"] = int(max_iter)

            opt = minimize(
                _scalar_fun,
                x0,
                method=mname,
                bounds=list(zip(lb, ub)) if has_bounds and mname in {"L-BFGS-B", "SLSQP", "Powell", "Nelder-Mead"} else None,
                tol=tol,
                options=options,
            )
            x_best = np.asarray(opt.x, dtype=float)
            merit = float(opt.fun)
            success = bool(opt.success)
            message = str(opt.message)
            nfev = int(getattr(opt, "nfev", -1))
            nit = int(getattr(opt, "nit", -1)) if hasattr(opt, "nit") else None

        _apply_x(x_best)
        merit_eval, target_reports, objective_reports, _state = self.evaluate()
        if not np.isfinite(merit_eval):
            merit_eval = merit

        if restore_if_fail and not success:
            _apply_x(x_saved)
            merit_eval, target_reports, objective_reports, _state = self.evaluate()

        variables_out = {var.name: var.get() for var in active_vars}
        return MatchResult(
            success=success,
            message=message,
            solver=solver,
            merit=float(merit_eval),
            nfev=nfev,
            nit=nit,
            variables=variables_out,
            target_reports=target_reports,
            objective_reports=objective_reports,
            optimize_result=opt,  # type: ignore[arg-type]
        )


__all__ = [
    "Vary",
    "PowerSupplyVary",
    "Target",
    "Objective",
    "TwissTarget",
    "TwissDifferenceTarget",
    "TwissPeriodicityTarget",
    "GlobalTwissTarget",
    "RMatrixElementTarget",
    "RMatrixBlockTarget",
    "TotalLengthTarget",
    "TargetReport",
    "ObjectiveReport",
    "FunctionObjective",
    "I5Objective",
    "MatchState",
    "MatchResult",
    "MatchProblem",
]
