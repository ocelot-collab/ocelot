from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.undulator_atom import UndulatorAtom
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Undulator(OpticElement):
    def __init__(self, lperiod=0., nperiods=0, Kx=0., Ky=0., field_file=None, eid=None, tm=TransferMap, **params):
        super().__init__(UndulatorAtom(lperiod=lperiod, nperiods=nperiods, Kx=Kx, Ky=Ky, field_file=field_file, eid=eid), tm=tm, default_tm=TransferMap, **params)

    