from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.drift_atom import DriftAtom
from ocelot.cpbd.transformations.kick import KickTM
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaTM
from ocelot.cpbd.transformations.runge_kutta_tr import RungeKuttaTrTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Drift(OpticElement):
    """
    drift - free space
    l - length of drift in [m]
    """
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM, KickTM, RungeKuttaTM, RungeKuttaTrTM}

    def __init__(self, l=0., eid=None, tm=None, **kwargs):
        super().__init__(DriftAtom(l=l, eid=eid, **kwargs), tm=tm)
