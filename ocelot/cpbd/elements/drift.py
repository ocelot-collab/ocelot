from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.drift_atom import DriftAtom
from ocelot.cpbd.transformations.kick import KickTM
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaGlobalTM, RungeKuttaOcelotTM, RungeKuttaTM
from ocelot.cpbd.transformations.runge_kutta_tr import RungeKuttaTrTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Drift(OpticElement):
    """
    Drift: field-free propagation space.

    Parameters
    ----------
    l : float, default=0
        Length of drift in [m]

    Physics
    -------
    Pure drift (no fields). Simplest element family.
    Always available as reference element.

    Tracking Methods
    ----------------
    Supports multiple tracking methods:
    - TransferMap (default): linear first-order mapping
    - SecondTM: second-order nonlinear mapping
    - KickTM: kick-style tracking
    - RungeKuttaGlobalTM / RungeKuttaTM: fixed-frame RK tracking
    - RungeKuttaOcelotTM: RK tracking converted back to Ocelot coordinates
    - RungeKuttaTrTM: transverse-only fixed-frame RK tracking

    Architecture
    ~~~~~~~~~~~~~
    - Wrapper: Drift (this class)
    - Atom: DriftAtom
    - Edge-aware: No (single MAIN map)

    See Also
    --------
    https://ocelot-collab.github.io/docs/docu/elements/drift/
    DriftAtom : Physics implementation
    """
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM, KickTM, RungeKuttaGlobalTM, RungeKuttaOcelotTM, RungeKuttaTM, RungeKuttaTrTM}

    def __init__(self, l=0., eid=None, tm=None, **kwargs):
        super().__init__(DriftAtom(l=l, eid=eid, **kwargs), tm=tm)
