from ocelot.cpbd.r_matrix import uni_matrix
from ocelot.cpbd.elements.bend_atom import BendAtom
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams


class RBendAtom(BendAtom):
    """
    rectangular bending magnet,
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

    def __init__(self, l=0., angle=0., k1=0., k2=0., e1=None, e2=None, tilt=0.,
                 gap=0, h_pole1=0., h_pole2=0., fint=0., fintx=None, eid=None):
        if e1 is None:
            e1 = angle / 2.
        else:
            e1 += angle / 2.
        if e2 is None:
            e2 = angle / 2.
        else:
            e2 += angle / 2.

        super().__init__(l=l, angle=angle, e1=e1, e2=e2, k1=k1, k2=k2, tilt=tilt,
                         gap=gap, h_pole1=h_pole1, h_pole2=h_pole2, fint=fint, fintx=fintx, eid=eid)

