from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.sbend_atom import SBendAtom
from ocelot.cpbd.transformations.kick import KickTM
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaGlobalTM, RungeKuttaOcelotTM, RungeKuttaTM
from ocelot.cpbd.transformations.runge_kutta_tr import RungeKuttaTrTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


class SBend(OpticElement):
    """
    Sector bending magnet: bend with edges perpendicular to reference trajectory.

    Parameters
    ----------
    l : float, default=0
        Length of magnet in [m]
    angle : float, default=0
        Bending angle in [rad]
    k1, k2 : float, default=0
        Quadrupole and sextupole strengths in [1/m²] and [1/m³]
    e1, e2, tilt : float, default=0
        Entrance/exit face angles and roll in [rad]
    gap, fint, fintx : float, default=0
        Fringe field modeling parameters
    h_pole1, h_pole2 : float, default=0
        Pole face curvatures (1/r)

    Physics
    -------
    Sector bend (edges perpendicular to reference trajectory at bend angle).
    Dispersion element with optional focusing. Edge-aware.

    Tracking Methods
    ----------------
    Same as Bend: TransferMap, SecondTM, KickTM, etc.

    Architecture
    ~~~~~~~~~~~~~
    - Wrapper: SBend (this class)
    - Atom: SBendAtom (inherits from Magnet)
    - Edge-aware: Yes (ENTRANCE → MAIN → EXIT)

    See Also
    --------
    Bend : General bending magnet
    RBend : Rectangular variant
    SBendAtom : Physics implementation
    """
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM, KickTM, RungeKuttaGlobalTM, RungeKuttaOcelotTM, RungeKuttaTM, RungeKuttaTrTM}

    def __init__(self, l=0., angle=0.0, k1=0.0, k2=0., e1=0.0, e2=0.0, tilt=0.0,
                 gap=0, h_pole1=0., h_pole2=0., fint=0., fintx=None, eid=None, tm=None, **kwargs):
        super().__init__(SBendAtom(l=l, angle=angle, k1=k1, k2=k2, e1=e1, e2=e2, tilt=tilt,
                      gap=gap, h_pole1=h_pole1, h_pole2=h_pole2, fint=fint, fintx=fintx, eid=eid, **kwargs), tm=tm)
