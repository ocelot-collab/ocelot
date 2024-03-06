from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.undulator_atom import UndulatorAtom
from ocelot.cpbd.transformations.transfer_map import TransferMap


class Undulator(OpticElement):
    def __init__(self, lperiod=0., nperiods=0, Kx=0., Ky=0., phase=0, end_poles='1', field_file=None, eid=None, tm=TransferMap, **params):
        super().__init__(UndulatorAtom(lperiod=lperiod, nperiods=nperiods, Kx=Kx, Ky=Ky, phase=phase, end_poles=end_poles, field_file=field_file, eid=eid),
                         tm=tm, default_tm=TransferMap, **params)

        """
        lperiod : 
            undulator period length.
        nperiods : 
            Number of undulator periods. The default is None.
        Kx : 
            undulator K parameter in x direction.
        Ky : 
            undulator K parameter in y direction.
        phase : optional
            Phase from which the magnetic field is stated. The default is 0, which is cos().
        end_poles : optional
            Correction poles to close magnetic field integrals. 
            Might be: 
                '0' - magnetic field starts from 0 value or sin()-like
                '1' - magnetic field starts from maximum value or cos()-like
                '3/4' - magnetic field starts with 0, 1/4, -3/4, +1, -1  (or 0, -1/4, +3/4, -1, +1) poles sequence and finishes with it (fraction is from maximum value of the field)
                '1/2' - magnetic field starts with 1/2 and finishes with -1/2 poles
            The default is '1'.
        """