import numpy as np

from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.quadrupole_atom import QuadrupoleAtom
from ocelot.cpbd.transformations.transfer_map import TransferMap 

class Quadrupole(OpticElement):
    """
    quadrupole
    l - length of lens in [m],
    k1 - strength of quadrupole lens in [1/m^2],
    k2 - strength of sextupole lens in [1/m^3],
    tilt - tilt of lens in [rad].
    """
    def __init__(self, l=0., k1=0, k2=0., tilt=0., eid=None, tm=TransferMap):
        super().__init__(QuadrupoleAtom(l=l, k1=k1, k2=k2, tilt=tilt, eid=eid), tm=tm, default_tm=TransferMap)
