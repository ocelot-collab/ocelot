import numpy as np

from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.multipole_atom import MultipoleAtom
from ocelot.cpbd.transformations.multipole import MultipoleTM


class Multipole(OpticElement):
    """
    Arbitrary multipole expansion element (dipole through high order).

    Parameters
    ----------
    kn : float or array-like, default=0
        Multipole strength coefficients (list or single value)

    Physics
    -------
    Arbitrary multipole field expansion.
    Tracking uses dedicated MultipoleTM algorithm
    (not standard TransferMap).

    Tracking Method
    ----------------
    Restricted to MultipoleTM only (cannot use TransferMap or SecondTM actively).
    Linear optics path available via first_order_tms for Twiss/optics.

    Architecture
    ~~~~~~~~~~~~~
    - Wrapper: Multipole (this class, restricted active TM)
    - Atom: MultipoleAtom
    - Default TM: MultipoleTM
    - Supported TMs: {MultipoleTM} only (single active method)
    - Edge-aware: No (single MAIN map)

    See Also
    --------
    https://ocelot-collab.github.io/docs/docu/elements/intro/
    MultipoleAtom : Physics implementation
    MultipoleTM : Dedicated multipole tracking algorithm
    """
    default_tm = MultipoleTM
    supported_tms = {MultipoleTM}

    def __init__(self, kn=0., eid=None, tm=None, **kwargs):
        super().__init__(MultipoleAtom(kn=kn, eid=eid, **kwargs), tm=tm)
