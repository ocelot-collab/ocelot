import warnings

import numpy as np

from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.multipole_atom import MultipoleAtom
from ocelot.cpbd.transformations.multipole import MultipoleTM
from ocelot.cpbd.transformations.transformation import Transformation


class Multipole(OpticElement):
    """
    kn - list of strengths
    """
    default_tm = MultipoleTM
    supported_tms = {MultipoleTM}

    def __init__(self, kn=0., eid=None, tm=MultipoleTM, **kwargs):
        if tm != MultipoleTM:
            warnings.warn(
                "Multipole does not declare support for "
                f"{tm.__name__}; falling back to default {MultipoleTM.__name__}.",
                stacklevel=2,
            )
            tm = MultipoleTM
        super().__init__(MultipoleAtom(kn=kn, eid=eid, **kwargs), tm=tm, default_tm=MultipoleTM)

    def set_tm(self, tm: Transformation, **params):
        """Pin the wrapper to MultipoleTM while keeping first-order maps available."""
        if tm != MultipoleTM:
            warnings.warn(
                "Multipole does not declare support for "
                f"{tm.__name__}; falling back to default {MultipoleTM.__name__}.",
                stacklevel=2,
            )
        return super().set_tm(MultipoleTM, **params)
