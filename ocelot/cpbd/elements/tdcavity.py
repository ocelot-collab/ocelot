from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.tdcavity_atom import TDCavityAtom
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


class TDCavity(OpticElement):
    """
    Transverse deflecting cavity - by default kick in horizontal plane
    l - length [m]
    v - total peak voltage for the full nominal length l [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    tilt - tilt of cavity in [rad]
    """
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM}

    def __init__(self, l=0., freq=0.0, phi=0.0, v=0., tilt=0.0, eid=None, tm=None, **kwargs):
        super().__init__(TDCavityAtom(l=l, freq=freq, phi=phi, v=v, tilt=tilt, eid=eid, **kwargs), tm=tm)
