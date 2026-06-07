import numpy as np

from ocelot.cpbd.elements.xyquadrupole_atom import XYQuadrupoleAtom
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.transformations.transfer_map import TransferMap


class XYQuadrupole(OpticElement):
    """
    Quadrupole with offsets (linear element). The element is to test a transport feature and it is not tested.

    l - length of magnet in [m],
    k1 - strength of quadrupole lens in [1/m^2],
    x_offs - offset in horizontal direction in [m]
    y_offs - offset in vertical direction in [m]
    tilt - tilt of lens in [rad],
    """
    default_tm = TransferMap
    supported_tms = {TransferMap}

    def __init__(self, l=0., x_offs=0.0, y_offs=0.0, k1=0.0, tilt=0.0, eid=None, tm=None, **kwargs):
        super().__init__(XYQuadrupoleAtom(l=l, x_offs=x_offs, y_offs=y_offs, k1=k1, tilt=tilt, eid=eid, **kwargs), tm=tm)
