from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.tdcavity_atom import TDCavityAtom
from ocelot.cpbd.transformations.transfer_map import TransferMap


class TDCavity(OpticElement):
    """
    Transverse deflecting cavity - by default kick in horizontal plane
    l - length [m]
    v - voltage [GV/m]
    freq - frequency [Hz]
    phi - phase in [deg]
    tilt - tilt of cavity in [rad]
    """
    def __init__(self, l=0., freq=0.0, phi=0.0, v=0., tilt=0.0, eid=None, tm=TransferMap):
        super().__init__(TDCavityAtom(l=l, freq=freq, phi=phi, v=v, tilt=tilt, eid=eid), tm=tm, default_tm=TransferMap)