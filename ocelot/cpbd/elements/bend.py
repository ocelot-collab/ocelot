import numpy as np

from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.bend_atom import BendAtom
from ocelot.cpbd.transformations.kick import KickTM
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaGlobalTM, RungeKuttaOcelotTM, RungeKuttaTM
from ocelot.cpbd.transformations.runge_kutta_tr import RungeKuttaTrTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Bend(OpticElement):
    """
    General bending magnet with optional quadrupole and edge effects.

    Parameters
    ----------
    l : float, default=0
        Length of magnet in [m]
    angle : float, default=0
        Bending angle in [rad] (MAD8 convention: positive = bend right)
    k1 : float, default=0
        Quadrupole (focusing) strength in [1/m²]
    k2 : float, default=0
        Sextupole (coupling) strength in [1/m³]
    e1, e2 : float, default=0
        Entrance and exit face angles in [rad]
    tilt : float, default=0
        Roll angle around beam axis in [rad]
    gap : float, default=0
        Magnet gap in [m] (NOTE: HGAP=gap/2 in MAD/ELEGANT)
    fint, fintx : float, default=0
        Fringe field integrals (entrance and exit)
    h_pole1, h_pole2 : float, default=0
        Curvature (1/r) of entrance and exit face poles

    Physics
    -------
    Bending magnet with dispersion and optional focusing.
    Edge-aware: builds ENTRANCE → MAIN → EXIT map sequence.

    Tracking Methods
    ----------------
    - TransferMap (default): linear maps
    - SecondTM: second-order nonlinear
    - KickTM: kick-style tracking
    - RungeKuttaGlobalTM / RungeKuttaTM: fixed-frame RK tracking
    - RungeKuttaOcelotTM: RK tracking converted back to Ocelot coordinates
    - RungeKuttaTrTM: transverse-only fixed-frame RK tracking

    Architecture
    ~~~~~~~~~~~~~
    - Wrapper: Bend (this class)
    - Atom: BendAtom (inherits from Magnet)
    - Edge-aware: Yes (ENTRANCE → MAIN → EXIT)

    See Also
    --------
    https://ocelot-collab.github.io/docs/docu/elements/bend/
    SBend : Sector bending magnet variant
    RBend : Rectangular bending magnet variant
    BendAtom : Physics implementation with edge modeling
    """
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM, KickTM, RungeKuttaGlobalTM, RungeKuttaOcelotTM, RungeKuttaTM, RungeKuttaTrTM}

    def __init__(self, l=0., angle=0., k1=0., k2=0., e1=0., e2=0., tilt=0.0,
                 gap=0., h_pole1=0., h_pole2=0., fint=0., fintx=None, eid=None, tm=None, **kwargs):
        super().__init__(BendAtom(l=l, angle=angle, k1=k1, k2=k2, e1=e1, e2=e2, tilt=tilt,
                                  gap=gap, h_pole1=h_pole1, h_pole2=h_pole2, fint=fint, fintx=fintx, eid=eid, **kwargs), tm=tm)
