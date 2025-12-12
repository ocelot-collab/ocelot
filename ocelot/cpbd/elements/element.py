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
    Element is a basic beamline building element
    Accelerator optics elements are subclasses of Element
    Arbitrary set of additional parameters can be attached if necessary
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
        self.w = kwargs.pop('width', 0.)    # transverse width of the element. it used only in the beamline plotting
        self.h = kwargs.pop('height', 0.)     # transverse width of the element. it used only in the beamline plotting
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
        T = t_nnn(delta_length if delta_length is not None else self.l, 0., 0., 0.,
                  energy)
        first_order_params = self.create_first_order_main_params(energy, delta_length)

        m_off = np.array([[-self.dx], [0], [-self.dy], [0], [0], [0]])
        # alternative way
        # B = first_order_params.B + np.dot(np.dot(m_off.T, T), m_off)[0]
        B, R = map_transform_with_offsets(m_off, first_order_params.R, T)
        return SecondOrderParams(R, B, T, self.tilt, self.dx, self.dy)

    def create_delta_e(self, total_length, delta_length=0.0):
        return 0.0

    def __repr__(self):
        return f"<{type(self).__name__}: name={self.id} at {hex(id(self))}>"

    def get_transfer_geometry(self):
        """
        Returns:
        R_end: Displacement to end [x, y, s]
        S_end: Rotation at end
        R_mid: Displacement to midpoint [x, y, s]
        S_mid: Rotation at midpoint
        """
        # 1. Straight Element Logic
        R_end = np.array([0, 0, self.l])
        R_mid = np.array([0, 0, self.l / 2.0])
        S_end = np.eye(3)
        S_mid = np.eye(3)  # No rotation halfway through a straight line

        # 2. Apply Tilt
        if self.tilt != 0:
            from ocelot.common.math_op import get_tilt_matrix
            T = get_tilt_matrix(self.tilt)

            R_end = T @ R_end
            R_mid = T @ R_mid

            # S = T @ I @ T.T = I, so S doesn't change for straight elements
            # UNLESS you want to track the frame roll explicitly.
            # But usually Tilt is static for the element.
            # We keep S as Identity for straight elements relative to the trajectory.
            pass

        return R_end, S_end, R_mid, S_mid