from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.solenoid_atom import SolenoidAtom
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Solenoid(OpticElement):
    """
    Solenoid
    l - length in m,
    k - strength B0/(2B*rho)
    """

    def __init__(self, l=0., k=0., eid=None, tm=TransferMap):
        super().__init__(SolenoidAtom(l=l, k=k, eid=eid), tm=tm, default_tm=TransferMap)
