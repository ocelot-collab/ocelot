import warnings

from ocelot.cpbd.elements.twcavity_atom import TWCavityAtom
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.transformations.tw_cavity import TWCavityTM
from ocelot.cpbd.transformations.transformation import Transformation


class TWCavity(OpticElement):
    """
    Traveling wave cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    """
    default_tm = TWCavityTM
    supported_tms = {TWCavityTM}

    def __init__(self, l=0., v=0., phi=0., freq=0., eid=None, tm=TWCavityTM, **kwargs):
        if tm != TWCavityTM:
            warnings.warn(
                "TWCavity does not declare support for "
                f"{tm.__name__}; falling back to default {TWCavityTM.__name__}.",
                stacklevel=2,
            )
            tm = TWCavityTM
        super().__init__(TWCavityAtom(l=l, v=v, phi=phi, freq=freq, eid=eid, **kwargs), tm=tm,
                         default_tm=TWCavityTM)

    def set_tm(self, tm: Transformation, **params):
        """Pin the wrapper to TWCavityTM until other paths are audited."""
        if tm != TWCavityTM:
            warnings.warn(
                "TWCavity does not declare support for "
                f"{tm.__name__}; falling back to default {TWCavityTM.__name__}.",
                stacklevel=2,
            )
        return super().set_tm(TWCavityTM, **params)
