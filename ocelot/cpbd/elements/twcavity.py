from ocelot.cpbd.elements.twcavity_atom import TWCavityAtom
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.transformations.tw_cavity import TWCavityTM


class TWCavity(OpticElement):
    """
    Traveling wave cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    """
    default_tm = TWCavityTM
    supported_tms = {TWCavityTM}

    def __init__(self, l=0., v=0., phi=0., freq=0., eid=None, tm=None, **kwargs):
        super().__init__(TWCavityAtom(l=l, v=v, phi=phi, freq=freq, eid=eid, **kwargs), tm=tm)
