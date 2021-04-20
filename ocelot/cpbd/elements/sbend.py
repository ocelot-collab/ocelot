from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.sbend_atom import SBendAtom
from ocelot.cpbd.transformations.transfer_map import TransferMap

class SBend(OpticElement):
    """
    sector bending magnet,
    l - length of magnet in [m],
    angle - angle of bend in [rad],
    k1 - strength of quadrupole lens in [1/m^2],
    k2 - strength of sextupole lens in [1/m^3],
    tilt - tilt of lens in [rad],
    e1 - the angle of inclination of the entrance face [rad],
    e2 - the angle of inclination of the exit face [rad].
    fint - fringe field integral
    fintx - allows (fintx > 0) to set fint at the element exit different from its entry value.
    gap - the magnet gap [m], NOTE in MAD and ELEGANT: HGAP = gap/2
    h_pole1 - the curvature (1/r) of the entrance face
    h_pole1 - the curvature (1/r) of the exit face
    """

    def __init__(self, l=0., angle=0.0, k1=0.0, k2=0., e1=0.0, e2=0.0, tilt=0.0,
                 gap=0, h_pole1=0., h_pole2=0., fint=0., fintx=None, eid=None, tm=TransferMap):
        super().__init__(SBendAtom(l=l, angle=angle, k1=k1, k2=k2, e1=e1, e2=e2, tilt=tilt,
                      gap=gap, h_pole1=h_pole1, h_pole2=h_pole2, fint=fint, fintx=fintx, eid=eid), tm=tm, default_tm=TransferMap)