import numpy as np
import logging

from ocelot.cpbd.transformations.transformation import Transformation
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.hcor_atom import HcorAtom
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap

logger = logging.getLogger(__name__)


class Hcor(OpticElement):
    """
    horizontal corrector,
    l - length of magnet in [m],
    angle - angle of bend in [rad],
    """
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM}

    def __init__(self, l=0., angle=0., eid=None, tm=TransferMap, **kwargs):
        super().__init__(HcorAtom(l=l, angle=angle, eid=eid, **kwargs), tm=tm, default_tm=TransferMap)
