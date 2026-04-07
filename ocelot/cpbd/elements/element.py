import numpy as np

from ocelot.common.math_op import get_tilt_matrix
from ocelot.cpbd.tm_params.second_order_params import SecondOrderParams
from ocelot.cpbd.tm_params.first_order_params import FirstOrderParams
from ocelot.cpbd.tm_params.kick_params import KickParams
from ocelot.cpbd.high_order import t_nnn
from ocelot.cpbd.r_matrix import uni_matrix
from ocelot.cpbd.tm_utils import map_transform_with_offsets


class Element:
    """
    Atom layer: physics state and hook methods for transformations.

    Role in Architecture
    ====================
    Element (and subclasses Magnet, etc.) is the "atom" in the wrapper/atom design:
    - Stores physics parameters (l, k1, k2, offsets, etc.)
    - Implements create_*_params(...) hooks that build TMParams objects
    - Transformations use these TMParams to track particles

    Public wrappers (Quadrupole, Drift, etc.) subclass OpticElement and wrap
    an Element instance, forwarding user-facing API to this physics object.

    Hook Methods
    =============
    Subclasses override hooks to provide physics-specific parameters:

    - create_first_order_main_params(energy) → FirstOrderParams
    - create_second_order_main_params(energy) → SecondOrderParams
    - create_kick_*_params(...) → KickParams (for kick-based tracking)
    - create_delta_e(...) → float (reference energy change)

    Edge-aware families also implement entrance and exit variants.

    Default Behavior
    =================
    Default implementations are generic fallbacks for simple elements.
    Individual families override these with physics-specific code.

    See Also
    --------
    https://ocelot-collab.github.io/docs/docu/elements/element/
    https://ocelot-collab.github.io/docs/docu/how-to/new_element/
    Magnet : Base class for magnetic families (quadrupoles, bends, etc.)
    OpticElement : Public wrapper that uses this atom
    """

    def __init__(self, eid=None, has_edge=False, **kwargs):
        self.has_edge = has_edge
        self.id = eid
        if eid is None:
            self.id = "ID_{0}_".format(np.random.randint(100000000))
        self.l = 0.
        self.angle = 0.
        self.tilt = 0.  # rad, pi/4 to turn positive quad into negative skew
        self.dx = 0.
        self.dy = 0.
        self.width = kwargs.pop('width', 0.05)    # transverse width of the element. it used only in the beamline plotting
        self.height = kwargs.pop('height', 0.05)     # transverse height of the element. it used only in the beamline plotting
        self.color = kwargs.pop('color', None) # color of element. it used only in the beamline plotting
        self.params = {}

    def __hash__(self):
        return hash(id(self))
        # return hash((self.id, self.__class__))

    def __eq__(self, other):
        try:
            # return (self.id, type) == (other.id, type)
            return id(self) == id(other)
        except:
            return False

    def _default_B(self, R):
        return np.dot((np.eye(6) - R), np.array([[self.dx], [0.], [self.dy], [0.], [0.], [0.]]))

    def create_first_order_main_params(self, energy: float, delta_length: float = None) -> FirstOrderParams:
        # Generic straight-element first-order fallback used when a concrete
        # family does not override the main map construction.
        if self.l == 0:
            hx = 0.
        else:
            hx = self.angle / self.l
        if delta_length is not None:
            R = uni_matrix(delta_length, 0., hx=hx, sum_tilts=0, energy=energy)
        else:
            R = uni_matrix(self.l, 0., hx=hx, sum_tilts=0, energy=energy)
        B = self._default_B(R)

        return FirstOrderParams(R, B, self.tilt)

    def create_second_order_main_params(self, energy: float, delta_length: float = 0.0) -> SecondOrderParams:
        # Generic straight-element second-order fallback. Families with custom
        # nonlinear behavior should override this instead of relying on it.
        T = t_nnn(delta_length if delta_length is not None else self.l, 0., 0., 0.,
                  energy)
        first_order_params = self.create_first_order_main_params(energy, delta_length)

        m_off = np.array([[-self.dx], [0], [-self.dy], [0], [0], [0]])
        # alternative way
        # B = first_order_params.B + np.dot(np.dot(m_off.T, T), m_off)[0]
        B, R = map_transform_with_offsets(m_off, first_order_params.R, T)
        return SecondOrderParams(R, B, T, self.tilt, self.dx, self.dy)

    def create_delta_e(self, total_length, delta_length=0.0):
        # Passive default: only active RF-like families are expected to
        # override the reference-energy change.
        return 0.0

    def __repr__(self):
        return f"<{type(self).__name__}: name={self.id} at {hex(id(self))}>"

    def get_transfer_geometry(self):
        """
        Geometry for straight elements.

        Returns:
            R_end: displacement from entrance to exit in the *local* frame [x, y, s]
            S_end: rotation of the *local frame* at exit relative to entrance
            R_mid: displacement from entrance to geometric center in the *local* frame
            S_mid: rotation of the *local frame* at midpoint relative to entrance

        Notes (MAD-style):
            - For straight elements, TILT does NOT change the reference trajectory.
            - Therefore, R_* is purely along +s.
            - S_* is identity (no change of the reference frame).
            - If you want the *magnet body* (box/cylinder) to appear rotated (e.g. skew quad),
              apply TILT in the plotting stage (extra roll around s for the mesh), not here.
        """
        L = getattr(self, "l", 0.0)

        # Straight trajectory along local s
        R_end = np.array([0.0, 0.0, L], dtype=float)
        R_mid = np.array([0.0, 0.0, 0.5 * L], dtype=float)

        # No rotation of the reference frame for a straight segment
        S_end = np.eye(3, dtype=float)
        S_mid = np.eye(3, dtype=float)

        return R_end, S_end, R_mid, S_mid
