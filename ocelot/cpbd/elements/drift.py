from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.drift_atom import DriftAtom
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Drift(OpticElement):
    """
    drift - free space
    l - length of drift in [m]
    """
    def __init__(self, l=0., eid=None, tm=TransferMap):
        super().__init__(DriftAtom(l=l, eid=eid), tm=tm, default_tm=TransferMap)
