from __future__ import annotations
import numpy as np

from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.quadrupole_atom import QuadrupoleAtom
from ocelot.cpbd.transformations.kick import KickTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Quadrupole(OpticElement):
    """
    quadrupole
    l - length of lens in [m],
    k1 - strength of quadrupole lens in [1/m^2],
    k2 - strength of sextupole lens in [1/m^3],
    tilt - tilt of lens in [rad].
    """
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM, KickTM}

    def __init__(self, l: float = 0., k1: float = 0, k2: float = 0., tilt: float = 0., eid: str | None = None, tm=TransferMap, **kwargs):
        super().__init__(QuadrupoleAtom(l=l, k1=k1, k2=k2, tilt=tilt, eid=eid, **kwargs), tm=tm, default_tm=TransferMap)

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
