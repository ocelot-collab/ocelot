import numpy as np

from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.hcor_atom import HcorAtom
from ocelot.cpbd.transformations.transfer_map import TransferMap

class Hcor(OpticElement):
    """
    horizontal corrector,
    l - length of magnet in [m],
    angle - angle of bend in [rad],
    """
    def __init__(self, l=0., angle=0., eid=None, tm=TransferMap):
        super().__init__(HcorAtom(l=l, angle=angle, eid=eid), tm=tm, default_tm=TransferMap)
