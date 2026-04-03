import numpy as np

from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.multipole_atom import MultipoleAtom
from ocelot.cpbd.transformations.multipole import MultipoleTM


class Multipole(OpticElement):
    """
    kn - list of strengths
    """
    default_tm = MultipoleTM
    supported_tms = {MultipoleTM}

    def __init__(self, kn=0., eid=None, tm=MultipoleTM, **kwargs):
        super().__init__(MultipoleAtom(kn=kn, eid=eid, **kwargs), tm=tm, default_tm=MultipoleTM)
