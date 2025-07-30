import logging

from ocelot.cpbd.transformations.transformation import Transformation
from ocelot.cpbd.elements.twcavity_atom import TWCavityAtom
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.transformations.tw_cavity import TWCavityTM

logger = logging.getLogger(__name__)


class TWCavity(OpticElement):
    """
    Traveling wave cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    """
    def __init__(self, l=0., v=0., phi=0., freq=0., eid=None, tm=TWCavityTM):
        if tm != TWCavityTM:
            logger.warning("Cavity Element only support TWCavityTM. Set tm to TWCavityTM.")
            tm = TWCavityTM
        super().__init__(TWCavityAtom(l=l, v=v, phi=phi, freq=freq, eid=eid), tm=tm,
                         default_tm=TWCavityTM)

    def set_tm(self, tm: Transformation, **params):
        logger.debug("Cavity Element only support TWCavityTM. Set tm to TWCavityTM.")

