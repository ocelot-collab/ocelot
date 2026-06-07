from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.undulator_atom import UndulatorAtom
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaGlobalTM, RungeKuttaOcelotTM, RungeKuttaTM
from ocelot.cpbd.transformations.runge_kutta_tr import RungeKuttaTrTM
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.transfer_map import TransferMap
from ocelot.cpbd.transformations.undulator_test import UndulatorTestTM


class Undulator(OpticElement):
    """
    Periodic insertion device (undulator/wiggler) for FEL and radiation.

    Parameters
    ----------
    lperiod : float, default=0
        Undulator period length in [m]
    nperiods : int, default=0
        Number of undulator periods
    Kx : float, default=0
        Undulator K parameter in x direction (strength)
    Ky : float, default=0
        Undulator K parameter in y direction
    phase : float, default=0
        Phase angle of the magnetic field [rad]
    end_poles : str, default='1'
        Correction pole configuration at end
    field_file : str, optional
        File path for external field map (can override K parameters)

    Physics
    -------
    Periodic transverse magnetic field for radiation generation.
    Supports field map input or analytic K-parameter model.

    Tracking Methods
    ----------------
    - TransferMap (default): simplified paraxial approximation
    - RungeKuttaGlobalTM / RungeKuttaTM: fixed-frame field integration
    - RungeKuttaOcelotTM: RK field integration converted back to Ocelot coordinates
    - RungeKuttaTrTM: transverse-only fixed-frame RK variant
    - SecondTM: second-order analytic maps
    - UndulatorTestTM: testing/simplified model

    Architecture
    ~~~~~~~~~~~~~
    - Wrapper: Undulator (this class)
    - Atom: UndulatorAtom
    - Edge-aware: No (single MAIN map)

    Notes
    -----
    MagneticLattice handles special length scaling when field map is attached.

    See Also
    --------
    UndulatorAtom : Physics implementation
    """
    default_tm = TransferMap
    supported_tms = {TransferMap, SecondTM, RungeKuttaGlobalTM, RungeKuttaOcelotTM, RungeKuttaTM, RungeKuttaTrTM, UndulatorTestTM}

    def __init__(self, lperiod=0., nperiods=0, Kx=0., Ky=0., phase=0, end_poles='1', field_file=None,
                 mag_field=None, eid=None, tm=None, **params):
        atom = UndulatorAtom(lperiod=lperiod, nperiods=nperiods, Kx=Kx, Ky=Ky, phase=phase,
                             end_poles=end_poles, field_file=field_file, eid=eid)
        atom.mag_field = mag_field
        super().__init__(atom, tm=tm, **params)
