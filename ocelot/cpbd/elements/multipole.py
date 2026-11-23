import numpy as np

from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.elements.multipole_atom import MultipoleAtom
from ocelot.cpbd.transformations.multipole import MultipoleTM


class Multipole(OpticElement):
    """Zero-length magnet with integrated multipole strengths.

    Parameters
    ----------
    kn : scalar or array-like, default=0
        A real ``kn[1] = K1L`` [1/m] is the integrated quadrupole strength
        (positive focuses in x). A real ``kn[2] = K2L`` [1/m²] is the
        integrated sextupole strength. Higher orders follow in sequence.
        ``kn[0]`` is the thin-bend angle [rad]. A scalar sets only this angle.
    tilt, dtilt : float, default=0
        Design and mechanical roll in radians about the axis at (dx, dy).
    dx, dy : float, default=0
        Magnetic-axis offsets in nominal beam coordinates, in metres.

    Physics
    -------
    For a normal thin sextupole, use ``Multipole(kn=[0., 0., K2L])``.
    In magnet coordinates the kicks are ``Δpx = -K2L*(x²-y²)/2`` and
    ``Δpy = K2L*x*y``, with momenta in units of the reference momentum.
    A real bend angle supplies angular dispersion and the corresponding
    path-length term. Its transverse kick vanishes on momentum.

    Tracking uses ``MultipoleTM`` automatically. The default ``twiss`` uses
    the first-order design map. To include the focusing of a sextupole at a
    horizontal offset from its axis, pass the entrance particle to
    ``twiss(line, tws0, reference_orbit=particle)``.

    See Also
    --------
    https://ocelot-collab.github.io/docs/docu/elements/intro/
    MultipoleAtom : Physics implementation
    MultipoleTM : Dedicated multipole tracking algorithm
    """
    default_tm = MultipoleTM
    supported_tms = {MultipoleTM}

    def __init__(self, kn=0., eid=None, tm=None, *, tilt=0., **kwargs):
        super().__init__(MultipoleAtom(kn=kn, eid=eid, tilt=tilt, **kwargs), tm=tm)
