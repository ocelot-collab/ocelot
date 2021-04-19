from ocelot.cpbd.elements.matrix_atom import MatrixAtom
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Matrix(OpticElement):
    def __init__(self, l=0., delta_e=0, eid=None, tm=TransferMap, **kwargs):
        super().__init__(MatrixAtom(l=l, delta_e=delta_e, eid=eid, **kwargs), tm=tm, default_tm=TransferMap)
