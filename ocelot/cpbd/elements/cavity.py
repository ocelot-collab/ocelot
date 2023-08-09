import logging

from ocelot.cpbd.transformations.transformation import Transformation
from ocelot.cpbd.elements.cavity_atom import CavityAtom
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.transformations.cavity import CavityTM

logger = logging.getLogger(__name__)


class Cavity(OpticElement):
    """
    Standing wave RF cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    vx_{up/down}, vy_{up/down} - zero order kick of a {up/down}stream coupler
    vxx_{up/down}, vxy_{up/down} - first order kick  a {up/down}stream coupler
    """

    def __init__(self, l=0., v=0., phi=0., freq=0., vx_up=0, vy_up=0, vxx_up=0, vxy_up=0,
                 vx_down=0, vy_down=0, vxx_down=0, vxy_down=0, eid=None, tm=CavityTM):
        if tm != CavityTM:
            logger.warning("Cavity Element only support CavityTM. Set tm to CavityTM.")
            tm = CavityTM
        super().__init__(CavityAtom(l=l, v=v, phi=phi, freq=freq, vx_up=vx_up, vy_up=vy_up, vxx_up=vxx_up, vxy_up=vxy_up,
                                    vx_down=vx_down, vy_down=vy_down, vxx_down=vxx_down, vxy_down=vxy_down, eid=eid), tm=tm, default_tm=CavityTM)

    def set_tm(self, tm: Transformation, **params):
        logger.debug("Cavity Element only support CavityTM. Set tm to CavityTM.")

    def remove_coupler_kick(self):
        self.vx_up = 0
        self.vy_up = 0
        self.vxx_up = 0
        self.vxy_up = 0
        self.vx_down = 0
        self.vy_down = 0
        self.vxx_down = 0
        self.vxy_down = 0
