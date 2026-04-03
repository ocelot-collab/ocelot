from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.octupole_atom import OctupoleAtom
from ocelot.cpbd.transformations.kick import KickTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Octupole(OpticElement):
    """
    sextupole
    l - length of lens in [m],
    k3 - strength of octupole lens in [1/m^4].
    """
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM, KickTM}

    def __init__(self, l=0., k3=0., tilt=0., eid=None, tm=None, **kwargs):
        super().__init__(OctupoleAtom(l=l, k3=k3, tilt=tilt, eid=eid, **kwargs), tm=tm)
