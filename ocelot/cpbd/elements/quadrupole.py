from __future__ import annotations
import numpy as np

from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.quadrupole_atom import QuadrupoleAtom
from ocelot.cpbd.transformations.kick import KickTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Quadrupole(OpticElement):
    """
    Quadrupole lens: linear focusing/defocusing element.

    Parameters
    ----------
    l : float, default=0
        Length of lens in [m]
    k1 : float, default=0
        Quadrupole focusing strength in [1/m²]
        (positive focuses in x, defocuses in y)
    k2 : float, default=0
        Sextupole coupling strength in [1/m³]
    tilt : float, default=0
        Roll angle around beam axis in [rad]

    Physics
    -------
    Standard magnetic focusing element.
    First-order dynamics are linear (hence R-matrix).
    Second-order effects via sextupole coupling.

    Tracking Methods
    ----------------
    - TransferMap (default): linear R-matrix
    - SecondTM: second-order nonlinear with k2 effects
    - KickTM: kick-based algorithmic tracking

    Architecture
    ~~~~~~~~~~~~~
    - Wrapper: Quadrupole (this class)
    - Atom: QuadrupoleAtom (inherits from Magnet)
    - Edge-aware: No (single MAIN map)

    Convenience Properties
    ~~~~~~~~~~~~~~~~~~~~~~
    - k1l: integrated quadrupole strength (k1 * l)
    - k2l: integrated sextupole strength (k2 * l)

    See Also
    --------
    https://ocelot-collab.github.io/docs/docu/elements/quadrupole/
    QuadrupoleAtom : Physics implementation
    """
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM, KickTM}

    def __init__(self, l: float = 0., k1: float = 0, k2: float = 0., tilt: float = 0., eid: str | None = None, tm=None, **kwargs):
        super().__init__(QuadrupoleAtom(l=l, k1=k1, k2=k2, tilt=tilt, eid=eid, **kwargs), tm=tm)

    @property
    def k1l(self) -> float:
        return self.k1 * self.l

    @k1l.setter
    def k1l(self, value: float) -> None:
        self.k1 = value / self.l

    @property
    def k2l(self) -> float:
        return self.k2 * self.l

    @k2l.setter
    def k2l(self, value: float) -> None:
        self.k2 = value / self.l
