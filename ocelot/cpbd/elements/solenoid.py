from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.solenoid_atom import SolenoidAtom
from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.transformations.second_order import SecondTM


class Solenoid(OpticElement):
    """
    Solenoid
    l - length in m,
    k - strength B0/(2B*rho)
    """
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM}
    def __init__(self, l=0., k=0., eid=None, tm=TransferMap, **kwargs):
        super().__init__(SolenoidAtom(l=l, k=k, eid=eid, **kwargs), tm=tm, default_tm=TransferMap)
