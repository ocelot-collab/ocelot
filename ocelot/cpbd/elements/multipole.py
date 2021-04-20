import logging
from ocelot.cpbd.transformations.transformation import Transformation

import numpy as np

from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.multipole_atom import MultipoleAtom
from ocelot.cpbd.transformations.multipole import MultipoleTM
from ocelot.cpbd.transformations.transfer_map import TransferMap

logger = logging.getLogger(__name__)


class Multipole(OpticElement):
    """
    kn - list of strengths
    """

    def __init__(self, kn=0., eid=None, tm=MultipoleTM):
        if tm != MultipoleTM:
            logger.debug("Multipole Element only support Multipole as its transformation. Set tm to MultipoleTM.")
            tm = MultipoleTM
        super().__init__(MultipoleAtom(kn=kn, eid=eid), tm=tm, default_tm=MultipoleTM)

    def set_tm(self, tm: Transformation):
        logger.debug("Multipole Element only support Multipole as its transformation. Set tm to MultipoleTM.")
