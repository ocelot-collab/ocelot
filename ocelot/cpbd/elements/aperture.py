import numpy as np

from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.aperture_atom import ApertureAtom
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Aperture(OpticElement):
    """
    Aperture
    xmax - half size in horizontal plane in [m],
    ymax - half size in vertical plane in [m],
    type - "rect" or "ellipt".
    """
    def __init__(self, xmax=np.inf, ymax=np.inf, dx=0, dy=0, type="rect", eid=None, tm=TransferMap):
        super().__init__(ApertureAtom(xmax=xmax, ymax=ymax, dx=dx, dy=dy, type=type, eid=eid), tm=TransferMap, default_tm=TransferMap)