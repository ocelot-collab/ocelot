from ocelot.cpbd.elements.cavity_atom import CavityAtom
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.transformations.cavity import CavityTM


class Cavity(OpticElement):
    """
    Standing-wave RF cavity with coupler kicks and energy gain.

    Parameters
    ----------
    l : float, default=0
        Cavity body length in [m]
    v : float, default=0
        RF voltage in [GV]
    phi : float, default=0
        RF phase in [deg]
    freq : float, default=0
        RF frequency in [Hz]
    vx_up, vy_up : complex, default=0
        Zero-order coupler kicks (upstream) in both planes
    vxx_up, vxy_up : complex, default=0
        First-order coupler kicks (upstream)
    vx_down, vy_down, vxx_down, vxy_down : complex, default=0
        Downstream coupler kicks

    Physics
    -------
    RF acceleration with transverse dispersion from coupler kicks.
    Energy gain via voltage and RF phase.
    Edge-aware (ENTRANCE coupler → MAIN cavity → EXIT coupler).

    Tracking Method
    ----------------
    Restricted to CavityTM only (cannot be changed by user).
    Linear optics path still available via first_order_tms for Twiss.

    Architecture
    ~~~~~~~~~~~~~
    - Wrapper: Cavity (this class, restricted active TM)
    - Atom: CavityAtom
    - Default TM: CavityTM
    - Supported TMs: {CavityTM} only (single active method)
    - Edge-aware: Yes (ENTRANCE → MAIN → EXIT)

    Methods
    -------
    remove_coupler_kick() : Remove all coupler kicks (set to zero)

    See Also
    --------
    https://ocelot-collab.github.io/docs/docu/elements/cavity/
    CavityAtom : Physics implementation
    CavityTM : Active tracking method (RF acceleration algorithm)
    """
    default_tm = CavityTM
    supported_tms = {CavityTM}

    def __init__(self, l=0., v=0., phi=0., freq=0., vx_up=0, vy_up=0, vxx_up=0, vxy_up=0,
                 vx_down=0, vy_down=0, vxx_down=0, vxy_down=0, eid=None, tm=None, **kwargs):
        super().__init__(CavityAtom(l=l, v=v, phi=phi, freq=freq, vx_up=vx_up, vy_up=vy_up, vxx_up=vxx_up, vxy_up=vxy_up,
                                    vx_down=vx_down, vy_down=vy_down, vxx_down=vxx_down, vxy_down=vxy_down, eid=eid, **kwargs), tm=tm)

    def remove_coupler_kick(self):
        self.vx_up = 0
        self.vy_up = 0
        self.vxx_up = 0
        self.vxy_up = 0
        self.vx_down = 0
        self.vy_down = 0
        self.vxx_down = 0
        self.vxy_down = 0
