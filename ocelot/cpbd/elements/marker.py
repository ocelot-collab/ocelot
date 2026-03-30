from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.marker_atom import MarkerAtom
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Marker(OpticElement):
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM}

    def __init__(self, eid=None, tm=TransferMap, **kwargs):
        super().__init__(MarkerAtom(eid, **kwargs), tm=tm, default_tm=TransferMap)
