import numpy as np

from ocelot.cpbd.transformations.first_order import TransferMap
from ocelot.cpbd.transformations.second_order import SecondTM
from ocelot.cpbd.transformations.kick import KickTM
from ocelot.cpbd.transformations.runge_kutta import RungeKuttaTM
from ocelot.cpbd.transformations.runge_kutta_tr import RungeKuttaTrTM
from ocelot.cpbd.high_order import t_nnn
from ocelot.cpbd.r_matrix import rot_mtx, uni_matrix


class Element(object):
    """
    Element is a basic beamline building element
    Accelerator optics elements are subclasses of Element
    Arbitrary set of additional parameters can be attached if necessary
    """

    default_tm = TransferMap
    additional_tms = [SecondTM, KickTM, RungeKuttaTrTM, RungeKuttaTM]

    def __init__(self, eid=None):
        self.id = eid
        if eid is None:
            self.id = "ID_{0}_".format(np.random.randint(100000000))
        self.l = 0.
        self.tilt = 0.  # rad, pi/4 to turn positive quad into negative skew
        self.angle = 0.
        self.k1 = 0.
        self.k2 = 0.
        self.dx = 0.
        self.dy = 0.
        self.params = {}
        self.transfer_map = TransferMap()

    def __hash__(self):
        return hash(id(self))
        # return hash((self.id, self.__class__))

    def __eq__(self, other):
        try:
            # return (self.id, type) == (other.id, type)
            return id(self) == id(other)
        except:
            return False

    def create_r_matrix(self):
        k1 = self.k1
        if self.l == 0:
            hx = 0.
        else:
            hx = self.angle / self.l
        r_z_e = lambda z, energy: uni_matrix(z, k1, hx=hx, sum_tilts=0, energy=energy)
        return r_z_e

    def get_T_z_e_func(self):
        return lambda z, energy: t_nnn(z, 0. if self.l == 0 else self.angle / self.l, self.k1, self.k2,
                                       energy)

    def create_tm(self, method_params=None):
        """
        Creates a transfer map for the element.
        @param method_params: A dictionary with which the transfer map can be set for the element. The user can set a
        global transfer map which will be applied for the element if the element supports this function otherwise the
        default transfer map will be used. Which transfer map is supported can be seen by the attribute "additional_tms"
        and "default_tm".
        @return: None
        """
        params = None
        if method_params:
            if isinstance(method_params, dict):
                params = method_params
                global_tm = params.get('global', None)
                if global_tm is None:
                    global_tm = TransferMap
                if self.__class__.__name__ in params:
                    tm = params[self.__class__.__name__]
                else:
                    tm = global_tm
            else:
                # TODO: Remove this together with MethodTM
                params = method_params.params
                if hasattr(method_params, "global_method"):
                    global_tm = method_params.global_method
                else:
                    global_tm = TransferMap
                if self.__class__ in params:
                    tm = params[self.__class__]
                else:
                    tm = global_tm

        else:
            tm = TransferMap
        # check if the set transfer map is supported otherwise fall back to the default transfer map
        if self._is_tm_supported(tm):
            self.transfer_map = tm.create_from_element(self, params)
        else:
            self.transfer_map = self.default_tm.create_from_element(self, params)
        self._set_general_tm_parameter()

    def _set_general_tm_parameter(self):
        """
        Set general TransferMap parameters from element.
        :return: None
        """
        self.transfer_map.length = self.l
        self.transfer_map.dx = self.dx
        self.transfer_map.dy = self.dy
        tilt = self.tilt
        self.transfer_map.tilt = tilt
        self.transfer_map.R_z = lambda z, energy: np.dot(np.dot(rot_mtx(-tilt), self.create_r_matrix()(z, energy)),
                                                         rot_mtx(tilt))
        self.transfer_map.R = lambda energy: self.transfer_map.R_z(self.l, energy)

    def _is_tm_supported(self, tm):
        if tm == self.default_tm:
            return True
        for add_tm in self.additional_tms:
            if tm == add_tm:
                return True
        return False
