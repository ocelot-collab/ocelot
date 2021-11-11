from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.octupole_atom import OctupoleAtom
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Octupole(OpticElement):
    """
    sextupole
    l - length of lens in [m],
    k3 - strength of octupole lens in [1/m^4].
    """

    def __init__(self, l=0., k3=0., tilt=0., eid=None, tm=TransferMap):
        super().__init__(OctupoleAtom(l=l, k3=k3, tilt=tilt, eid=eid), tm=tm, default_tm=TransferMap)