from ocelot.cpbd.r_matrix import uni_matrix
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.rbend_atom import RBendAtom
from ocelot.cpbd.transformations.kick import KickTM
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaGlobalTM, RungeKuttaOcelotTM, RungeKuttaTM
from ocelot.cpbd.transformations.runge_kutta_tr import RungeKuttaTrTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


class RBend(OpticElement):
    """
    Rectangular bending magnet: bend with edges perpendicular to entrance/exit.

    Parameters
    ----------
    l : float, default=0
        Length of magnet in [m]
    angle : float, default=0
        Bending angle in [rad] (MAD8 convention)
    k1, k2 : float, default=0
        Quadrupole and sextupole strengths
    e1, e2, tilt : float, default=0
        Entrance/exit face angles and roll in [rad]
    gap, fint, fintx : float, default=0
        Fringe field model parameters
    h_pole1, h_pole2 : float, default=0
        Pole face curvatures (1/r)

    Physics
    -------
    Rectangular (cylindrical) bend.
    Edges perpendicular to reference trajectory at entrance/exit.
    Dispersion element. Edge-aware.

    Tracking Methods
    ----------------
    Same as Bend: TransferMap, SecondTM, KickTM, etc.

    Architecture
    ~~~~~~~~~~~~~
    - Wrapper: RBend (this class)
    - Atom: RBendAtom (inherits from Magnet)
    - Edge-aware: Yes (ENTRANCE → MAIN → EXIT)

    See Also
    --------
    Bend : General bending magnet
    SBend : Sector bend variant
    RBendAtom : Physics implementation
    """
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM, KickTM, RungeKuttaGlobalTM, RungeKuttaOcelotTM, RungeKuttaTM, RungeKuttaTrTM}

    def __init__(self, l=0., angle=0., k1=0., k2=0., e1=None, e2=None, tilt=0.,
                 gap=0, h_pole1=0., h_pole2=0., fint=0., fintx=None, eid=None, tm=None, **kwargs):
        super().__init__(RBendAtom(l=l, angle=angle, e1=e1, e2=e2, k1=k1, k2=k2, tilt=tilt,
                                   gap=gap, h_pole1=h_pole1, h_pole2=h_pole2, fint=fint, fintx=fintx, eid=eid, **kwargs), tm=tm)
