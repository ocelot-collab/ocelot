from ocelot.cpbd.elements.monitor_atom import MonitorAtom
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.transformations.transfer_map import TransferMap

# to mark locations of bpms and other diagnostics
class Monitor(OpticElement):

    def __init__(self, l=0.0, eid=None, tm=TransferMap):
        super().__init__(MonitorAtom(l=l, eid=eid), tm=tm, default_tm=TransferMap)