from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.sextupole_atom import SextupoleAtom
from ocelot.cpbd.transformations.kick import KickTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Sextupole(OpticElement):
    """
    sextupole
    l - length of lens in [m],
    k2 - strength of sextupole lens in [1/m^3].
    """
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM, KickTM}

    def __init__(self, l=0., k2=0., tilt=0., eid=None, tm=None, **kwargs):
        super().__init__(SextupoleAtom(l=l, k2=k2, tilt=tilt, eid=eid, **kwargs), tm=tm)
