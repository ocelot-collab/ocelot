from ocelot.cpbd.elements.twcavity_atom import TWCavityAtom
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.transformations.tw_cavity import TWCavityTM


class TWCavity(OpticElement):
    """
    Travelling-wave accelerating structure (phase velocity = c).

    Parameters
    ----------
    l : float, default=0
        Structure length in [m]
    v : float, default=0
        TOTAL RF voltage in [GV], not a gradient. 0.050 is 50 MV.
    phi : float, default=0
        RF phase in [deg], cosine convention: energy gain is v*cos(phi).
        0 = on crest, 90 = zero crossing, 180 = decelerating.
        Note elegant's RFCA uses sin(PHASE), so its PHASE = 90 - phi.
    freq : float, default=0
        RF frequency in [Hz]

    Physics
    -------
    The body has no ponderomotive focusing: R11 = 1, R21 = 0. The angle only
    damps adiabatically, and all the radial focusing comes from the entrance
    and exit edge lenses -gamma'/(2 gamma0) and +gamma'/(2 gamma1)
    (Rosenzweig & Serafini, PRE 49, 1599 (1994)). The body is written in
    momentum, R12 = p0 l [asinh(p1) - asinh(p0)] / dgamma and R22 = p0/p1,
    which matters below ~10 MeV where the beta -> 1 form is several percent
    off. This reproduces elegant's default RFCA matrix to 1e-12.

    The longitudinal map is identical to the standing-wave Cavity and is
    shared through ocelot.cpbd.rf_maps.

    Not modelled
    ------------
    Phase velocity is fixed at c and one phase is used for the whole
    structure, so there is no phase slippage; the gradient is uniform (no
    constant-impedance taper); there is no transverse RF force and no
    coupler kick. For a sub-MeV injector or a tapered structure these
    matter -- compare elegant's TWLA, which integrates them numerically.

    Tracking Method
    ---------------
    Restricted to TWCavityTM only (cannot be changed by the user).
    Linear optics path still available via first_order_tms for Twiss.

    Architecture
    ~~~~~~~~~~~~
    - Wrapper: TWCavity (this class, restricted active TM)
    - Atom: TWCavityAtom
    - Default TM: TWCavityTM
    - Supported TMs: {TWCavityTM} only
    - Edge-aware: Yes (ENTRANCE -> MAIN -> EXIT)

    See Also
    --------
    Cavity : Standing-wave structure, with coupler kicks
    TWCavityAtom : Physics implementation
    ocelot.cpbd.rf_maps : Longitudinal map shared with Cavity
    """
    default_tm = TWCavityTM
    supported_tms = {TWCavityTM}

    def __init__(self, l=0., v=0., phi=0., freq=0., eid=None, tm=None, **kwargs):
        super().__init__(TWCavityAtom(l=l, v=v, phi=phi, freq=freq, eid=eid, **kwargs), tm=tm)
