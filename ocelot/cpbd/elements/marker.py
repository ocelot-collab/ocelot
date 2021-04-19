from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.marker_atom import MarkerAtom
from ocelot.cpbd.transformations.transfer_map import TransferMap

class Marker(OpticElement):
    def __init__(self, eid=None, tm=TransferMap):
        super().__init__(MarkerAtom(eid), tm=tm, default_tm=TransferMap)