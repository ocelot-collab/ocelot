"""
definition of magnetic lattice
linear dimensions in [m]
"""
from ocelot.cpbd.field_map import FieldMap
from ocelot.cpbd.optics import *

import sys
import numpy as np
import logging

logger = logging.getLogger(__name__)


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
        self.dtilt = 0. # TODO delete
        self.params = {}
        #self.transfer_map = TransferMap()

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

    def _extend_element_def_string(self, params):
        """
        Overload this function to add element specific information to element_def_string.
        Note: This function is a hook for element_def_string.
        @param params:
        @return:
        """
        return params

    def element_def_string(self):
        params = []

        element_type = self.__class__.__name__
        element_ref = getattr(sys.modules[__name__], element_type)()
        params_order = element_ref.__init__.__code__.co_varnames
        argcount = element_ref.__init__.__code__.co_argcount

        for param in params_order[:argcount]:
            if param == 'self':
                continue

            # fix for parameter 'eid'
            if param == 'eid':
                params.append('eid=\'' + self.id + '\'')
                continue

            if isinstance(self.__dict__[param], np.ndarray):

                if not np.array_equal(self.__dict__[param], element_ref.__dict__[param]):
                    params.append(param + '=' + np.array2string(self.__dict__[param], separator=', '))
                continue

            if isinstance(self.__dict__[param], (int, float, complex)):

                # fix for parameters 'e1' and 'e2' in RBend element
                if element_type == 'RBend' and param in ('e1', 'e2'):
                    val = self.__dict__[param] - self.angle / 2.0
                    if val != 0.0:
                        params.append(param + '=' + str(val))
                    continue

                if self.__dict__[param] != element_ref.__dict__[param]:
                    params.append(param + '=' + str(self.__dict__[param]))
                continue

            if isinstance(self.__dict__[param], str):

                if self.__dict__[param] != element_ref.__dict__[param]:
                    params.append(param + '=\'' + self.__dict__[param] + '\'')
                continue

        params = self._extend_element_def_string(params)

        # join all parameters to element definition
        string = self._pprinting(element_type, params)
        return string

    def _pprinting(self, element_type, params):
        string = self.name + ' = ' + element_type + '('
        n0 = len(string)
        n = n0
        for i, param in enumerate(params):
            n += len(params)
            if n > 250:
                string += "\n"
                string += " " * n0 + param + ", "
                n = n0 + len(param) + 2
            else:
                if i == len(params) - 1:
                    string += param
                else:
                    string += param + ", "
        string += ")\n"
        return string

    def create_tm(self, method_params=None):
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
                #TODO: Remove this together with MethodTM
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

        if self._is_tm_supported(tm):
            self.transfer_map = tm.create_from_element(self, params)
        else:
            self.transfer_map = self.default_tm.create_from_element(self, params)
        self._set_general_tm_parameter()

    def _set_general_tm_parameter(self):
        self.transfer_map.length = self.l
        self.transfer_map.dx = self.dx
        self.transfer_map.dy = self.dy
        tilt = self.dtilt + self.tilt
        self.transfer_map.tilt = tilt
        self.transfer_map.R_z = lambda z, energy: np.dot(np.dot(rot_mtx(-tilt), self.create_r_matrix()(z, energy)), rot_mtx(tilt))
        self.transfer_map.R = lambda energy: self.transfer_map.R_z(self.l, energy)

    def _is_tm_supported(self, tm):
        if tm == self.default_tm:
            return True
        for add_tm in self.additional_tms:
            if tm == add_tm:
                return True
        return False

# to mark locations of bpms and other diagnostics
class Monitor(Element):

    def __init__(self, l=0.0, eid=None):
        Element.__init__(self, eid)
        self.l = l
        self.x_ref = 0.
        self.y_ref = 0.
        self.x = 0.
        self.y = 0.

    def __str__(self):
        s = 'Monitor : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l =%8.4f m\n' % self.l
        return s


class Marker(Element):
    def __init__(self, eid=None):
        Element.__init__(self, eid)
        self.l = 0.

    def __str__(self):
        s = 'Marker : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l =%8.4f m\n' % self.l
        return s


class Aperture(Element):
    """
    Aperture
    xmax - half size in horizontal plane in [m],
    ymax - half size in vertical plane in [m],
    type - "rect" or "elliptical".
    """

    def __init__(self, xmax=np.inf, ymax=np.inf, dx=0, dy=0, type="rect", eid=None):
        Element.__init__(self, eid)
        self.l = 0.
        self.xmax = xmax
        self.ymax = ymax
        self.dx = dx
        self.dy = dy
        self.type = type

    def __str__(self):
        s = 'Aperture : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l    =%8.4f m\n' % self.l
        s += 'xmax =%8.5f m\n' % self.xmax
        s += 'ymax =%8.5f m\n' % self.ymax
        s += 'dx   =%8.5f m\n' % self.dx
        s += 'dy   =%8.5f m\n' % self.dy
        s += 'type = %s \n' % self.type
        return s


class Quadrupole(Element):
    """
    quadrupole
    l - length of lens in [m],
    k1 - strength of quadrupole lens in [1/m^2],
    k2 - strength of sextupole lens in [1/m^3],
    tilt - tilt of lens in [rad].
    """

    def __init__(self, l=0., k1=0, k2=0., tilt=0., eid=None):
        # Element.__init__(self, eid)
        super(Quadrupole, self).__init__(eid=eid)
        self.l = l
        self.k1 = k1
        self.k2 = k2
        self.tilt = tilt

    def __str__(self):
        s = 'Quadrupole : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l    =%8.4f m\n' % self.l
        s += 'k1   =%8.3f 1/m^2\n' % self.k1
        s += 'k2   =%8.3f 1/m^3\n' % self.k2
        s += 'tilt =%8.2f deg\n' % (self.tilt * 180.0 / np.pi)
        return s


class Sextupole(Element):
    """
    sextupole
    l - length of lens in [m],
    k2 - strength of sextupole lens in [1/m^3].
    """

    def __init__(self, l=0., k2=0., tilt=0., eid=None):
        Element.__init__(self, eid)
        self.l = l
        self.k2 = k2
        self.tilt = tilt

    def __str__(self):
        s = 'Sextupole : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l    =%8.4f m\n' % self.l
        s += 'k2   =%8.3f 1/m^3\n' % self.k2
        s += 'tilt =%8.2f deg\n' % (self.tilt * 180.0 / np.pi)
        return s


class Octupole(Element):
    """
    octupole
    k3 - strength of octupole lens in [1/m^4],
    l - length of lens in [m].
    """

    def __init__(self, l=0., k3=0., tilt=0., eid=None):
        Element.__init__(self, eid)
        self.l = l
        self.k3 = k3
        self.tilt = tilt

    def __str__(self):
        s = 'Octupole : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l    =%8.4f m\n' % self.l
        s += 'k3   =%8.3f 1/m^4\n' % self.k3
        s += 'tilt =%8.2f deg\n' % (self.tilt * 180.0 / np.pi)
        return s


class Drift(Element):
    """
    drift - free space
    l - length of drift in [m]
    """

    def __init__(self, l=0., eid=None):
        Element.__init__(self, eid)
        self.l = l

    def __str__(self):
        s = 'Drift : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l =%8.4f m\n' % self.l
        return s


class Bend(Element):
    """
    bending magnet
    l - length of magnet in [m],
    angle - angle of bend in [rad],
    k1 - strength of quadrupole lens in [1/m^2],
    k2 - strength of sextupole lens in [1/m^3],
    tilt - tilt of lens in [rad],
    e1 - the angle of inclination of the entrance face [rad],
    e2 - the angle of inclination of the exit face [rad].
    fint - fringe field integral
    fintx - allows (fintx > 0) to set fint at the element exit different from its entry value.
    gap - the magnet gap [m], NOTE in MAD and ELEGANT: HGAP = gap/2
    h_pole1 - the curvature (1/r) of the entrance face
    h_pole1 - the curvature (1/r) of the exit face
    """

    def __init__(self, l=0., angle=0., k1=0., k2=0., e1=0., e2=0., tilt=0.0,
                 gap=0., h_pole1=0., h_pole2=0., fint=0., fintx=None, eid=None):
        Element.__init__(self, eid)
        self.l = l
        self.angle = angle
        self.k1 = k1
        self.k2 = k2
        self.e1 = e1
        self.e2 = e2
        self.gap = gap
        self.h_pole1 = h_pole1
        self.h_pole2 = h_pole2
        self.fint = fint
        self.fintx = fint
        if fintx is not None:
            self.fintx = fintx
        self.tilt = tilt

    def __str__(self):
        s = 'Bend : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l       =%8.4f m\n' % self.l
        s += 'angle   =%8.3f deg\n' % (self.angle * 180.0 / np.pi)
        s += 'e1      =%8.3f deg\n' % (self.e1 * 180.0 / np.pi)
        s += 'e2      =%8.3f deg\n' % (self.e2 * 180.0 / np.pi)
        s += 'tilt    =%8.3f deg\n' % (self.tilt * 180.0 / np.pi)
        s += 'fint    =%8.3f\n' % self.fint
        s += 'fintx   =%8.3f\n' % self.fintx
        s += 'gap     =%8.4f m\n' % self.gap
        s += 'h_pole1 =%8.4f 1/m\n' % self.h_pole1
        s += 'h_pole2 =%8.4f 1/m\n' % self.h_pole2
        return s


class Edge(Bend):
    def __init__(self, l=0., angle=0.0, k1=0., edge=0.,
                 tilt=0.0, dtilt=0.0, dx=0.0, dy=0.0,
                 h_pole=0., gap=0., fint=0., pos=1, eid=None):
        Element.__init__(self, eid)
        if l != 0.:
            self.h = angle / l
        else:
            self.h = 0
        self.l = 0.
        self.k1 = k1
        self.h_pole = h_pole
        self.gap = gap
        self.fint = fint
        self.edge = edge
        self.dx = dx
        self.dy = dy
        self.dtilt = dtilt
        self.tilt = tilt
        self.pos = pos

    def __str__(self):
        s = 'Edge : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'h      =%8.4f 1/m\n' % self.h
        s += 'fint   =%8.3f\n' % self.fint
        s += 'gap    =%8.4f m\n' % self.gap
        s += 'h_pole =%8.4f 1/m\n' % self.h_pole
        s += 'tilt   =%8.3f deg\n' % (self.tilt * 180.0 / np.pi)
        return s

    def create_r_matrix(self):
        sec_e = 1. / np.cos(self.edge)
        phi = self.fint * self.h * self.gap * sec_e * (1. + np.sin(self.edge) ** 2)
        # phi = self.fint * self.h * self.gap * sec_e * (1. + np.sin(2*self.edge) )
        r = np.eye(6)
        r[1, 0] = self.h * np.tan(self.edge)
        r[3, 2] = -self.h * np.tan(self.edge - phi)
        r_z_e = lambda z, energy: r
        return r_z_e

    def get_T_z_e_func(self):
        if self.pos == 1:
            _, T = fringe_ent(h=self.h, k1=self.k1, e=self.edge, h_pole=self.h_pole,
                              gap=self.gap, fint=self.fint)
        else:
            _, T = fringe_ext(h=self.h, k1=self.k1, e=self.edge, h_pole=self.h_pole,
                              gap=self.gap, fint=self.fint)
        return lambda z, energy: T


class SBend(Bend):
    """
    sector bending magnet,
    l - length of magnet in [m],
    angle - angle of bend in [rad],
    k1 - strength of quadrupole lens in [1/m^2],
    k2 - strength of sextupole lens in [1/m^3],
    tilt - tilt of lens in [rad],
    e1 - the angle of inclination of the entrance face [rad],
    e2 - the angle of inclination of the exit face [rad].
    fint - fringe field integral
    fintx - allows (fintx > 0) to set fint at the element exit different from its entry value.
    gap - the magnet gap [m], NOTE in MAD and ELEGANT: HGAP = gap/2
    h_pole1 - the curvature (1/r) of the entrance face
    h_pole1 - the curvature (1/r) of the exit face
    """

    def __init__(self, l=0., angle=0.0, k1=0.0, k2=0., e1=0.0, e2=0.0, tilt=0.0,
                 gap=0, h_pole1=0., h_pole2=0., fint=0., fintx=None, eid=None):
        Bend.__init__(self, l=l, angle=angle, k1=k1, k2=k2, e1=e1, e2=e2, tilt=tilt,
                      gap=gap, h_pole1=h_pole1, h_pole2=h_pole2, fint=fint, fintx=fintx, eid=eid)


class RBend(Bend):
    """
    rectangular bending magnet,
    l - length of magnet in [m],
    angle - angle of bend in [rad],
    k1 - strength of quadrupole lens in [1/m^2],
    k2 - strength of sextupole lens in [1/m^3],
    tilt - tilt of lens in [rad],
    e1 - the angle of inclination of the entrance face [rad],
    e2 - the angle of inclination of the exit face [rad].
    fint - fringe field integral
    fintx - allows (fintx > 0) to set fint at the element exit different from its entry value.
    gap - the magnet gap [m], NOTE in MAD and ELEGANT: HGAP = gap/2
    h_pole1 - the curvature (1/r) of the entrance face
    h_pole1 - the curvature (1/r) of the exit face
    """

    def __init__(self, l=0., angle=0., k1=0., k2=0., e1=None, e2=None, tilt=0.,
                 gap=0, h_pole1=0., h_pole2=0., fint=0., fintx=None, eid=None):
        if e1 is None:
            e1 = angle / 2.
        else:
            e1 += angle / 2.
        if e2 is None:
            e2 = angle / 2.
        else:
            e2 += angle / 2.

        Bend.__init__(self, l=l, angle=angle, e1=e1, e2=e2, k1=k1, k2=k2, tilt=tilt,
                      gap=gap, h_pole1=h_pole1, h_pole2=h_pole2, fint=fint, fintx=fintx, eid=eid)

    def create_r_matrix(self):
        r_z_e = lambda z, energy: uni_matrix(z, 0, hx=0, sum_tilts=0, energy=energy)
        return r_z_e


class XYQuadrupole(SBend):
    """
    Quadrupole with offsets (linear element). The element is to test a transport feature and it is not tested.

    l - length of magnet in [m],
    k1 - strength of quadrupole lens in [1/m^2],
    x_offs - offset in horizontal direction in [m]
    y_offs - offset in vertical direction in [m]
    tilt - tilt of lens in [rad],
    """

    def __init__(self, l=0., x_offs=0.0, y_offs=0.0, k1=0.0, tilt=0.0, eid=None):
        Element.__init__(self, eid)
        self.l = l
        self.k1 = k1
        self.x_offs = x_offs
        self.y_offs = y_offs
        self.tilt = tilt

    def create_r_matrix(self):
        k1 = self.k1

        if self.l == 0:
            hx = 0.
            hy = 0.
        else:
            hx = k1 * self.x_offs
            hy = -k1 * self.y_offs

        def r_mtx(z, k1, hx, hy, sum_tilts=0., energy=0.):
            # r = self.l/self.angle
            #  +K - focusing lens , -K - defoc
            gamma = energy / m_e_GeV

            kx2 = (k1 + hx * hx)
            ky2 = hy * hy - k1
            kx = np.sqrt(kx2 + 0.j)
            ky = np.sqrt(ky2 + 0.j)
            cx = np.cos(z * kx).real
            cy = np.cos(z * ky).real
            sy = (np.sin(ky * z) / ky).real if ky != 0 else z

            igamma2 = 0.

            if gamma != 0:
                igamma2 = 1. / (gamma * gamma)

            beta = np.sqrt(1. - igamma2)

            if kx != 0:
                sx = (np.sin(kx * z) / kx).real
                dx = hx / kx2 * (1. - cx)
                dy = hy / ky2 * (1. - cy)
                r56 = hx * hx * (z - sx) / kx2 / beta ** 2 + hy * hy * (z - sy) / ky2 / beta ** 2
            else:
                sx = z
                dx = z * z * hx / 2.
                dy = z * z * hy / 2.
                r56 = hx * hx * z ** 3 / 6. / beta ** 2 + hy * hy * z ** 3 / 6. / beta ** 2

            r56 -= z / (beta * beta) * igamma2

            u_matrix = np.array([[cx, sx, 0., 0., 0., dx / beta],
                                 [-kx2 * sx, cx, 0., 0., 0., sx * hx / beta],
                                 [0., 0., cy, sy, 0., dy / beta],
                                 [0., 0., -ky2 * sy, cy, 0., sy * hy / beta],
                                 [hx * sx / beta, dx / beta, hy * sy / beta, dy / beta, 1., r56],
                                 [0., 0., 0., 0., 0., 1.]])
            if sum_tilts != 0:
                u_matrix = np.dot(np.dot(rot_mtx(-sum_tilts), u_matrix), rot_mtx(sum_tilts))
            return u_matrix

        r_z_e = lambda z, energy: r_mtx(z, k1, hx=hx, hy=hy, sum_tilts=0, energy=energy)
        return r_z_e

    def get_T_z_e_func(self):
        return lambda z, energy: np.zeros((6, 6, 6))


class Hcor(RBend):
    """
    horizontal corrector,
    l - length of magnet in [m],
    angle - angle of bend in [rad],
    """

    default_tm = HCorrectorTM
    additional_tms = []

    def __init__(self, l=0., angle=0., eid=None):
        RBend.__init__(self, l=l, angle=angle, eid=eid)
        self.l = l
        self.angle = angle
        self.tilt = 0.

    def __str__(self):
        s = 'Hcor : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l     =%8.4f m\n' % self.l
        s += 'angle =%8.3f deg\n' % (self.angle * 180.0 / np.pi)
        s += 'tilt  =%8.2f deg\n' % (self.tilt * 180.0 / np.pi)
        return s


class Vcor(RBend):
    """
    horizontal corrector,
    l - length of magnet in [m],
    angle - angle of bend in [rad],
    """

    default_tm = VCorrectorTM
    additional_tms = []

    def __init__(self, l=0., angle=0., eid=None):
        RBend.__init__(self, l=l, angle=angle, eid=eid)
        self.l = l
        self.angle = angle
        self.tilt = np.pi / 2.

    def __str__(self):
        s = 'Vcor : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l     =%8.4f m\n' % self.l
        s += 'angle =%8.4f deg\n' % (self.angle * 180.0 / np.pi)
        s += 'tilt  =%8.2f deg\n' % (self.tilt * 180.0 / np.pi)
        return s


class Undulator(Element):
    """
    Undulator
    lperiod - undulator period in [m];\n
    nperiod - number of periods;\n
    Kx - undulator paramenter for vertical field; \n
    Ky - undulator parameter for horizantal field;\n
    field_file - absolute path to magnetic field data;\n
    mag_field - None by default, the magnetic field map function - (Bx, By, Bz) = f(x, y, z)
    eid - id of undulator.
    """

    additional_tms = [SecondTM, KickTM, RungeKuttaTrTM, RungeKuttaTM, UndulatorTestTM]

    def __init__(self, lperiod=0., nperiods=0, Kx=0., Ky=0., field_file=None, eid=None):
        Element.__init__(self, eid)
        self.lperiod = lperiod
        self.nperiods = nperiods
        self.l = lperiod * nperiods
        self.Kx = Kx
        self.Ky = Ky
        self.solver = "linear"  # can be "lin" is linear matrix,  "sym" - symplectic method and "rk" is Runge-Kutta
        self.phase = 0.  # phase between Bx and By + pi/4 (spiral undulator)

        self.ax = -1  # width of undulator, when ax is negative undulator width is infinite
        # I need this for analytic description of undulator

        self.field_file = field_file
        self.field_map = FieldMap(self.field_file)
        self.mag_field = None  # the magnetic field map function - (Bx, By, Bz) = f(x, y, z)
        self.v_angle = 0.
        self.h_angle = 0.

    def __str__(self):
        s = 'Undulator : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l        =%8.4f m\n' % self.l
        s += 'nperiods =%8.1f \n' % self.nperiods
        s += 'lperiod  =%8.4f m\n' % self.lperiod
        s += 'Kx       =%8.3f \n' % self.Kx
        s += 'Ky       =%8.3f \n' % self.Ky
        return s

    def create_r_matrix(self):
        """
        in OCELOT coordinates:
        R56 = - Lu/(gamma**2 * beta**2) * (1 + 0.5 * K**2 * beta**2)
        S.Tomin, Varenna, 2017.
        """

        def undulator_r_z(z, lperiod, Kx, Ky, energy):
            gamma = energy / m_e_GeV
            r = np.eye(6)
            r[0, 1] = z
            if gamma != 0 and lperiod != 0 and Kx != 0:
                beta = 1 / np.sqrt(1.0 - 1.0 / (gamma * gamma))

                omega_x = np.sqrt(2.0) * np.pi * Kx / (lperiod * gamma * beta)
                omega_y = np.sqrt(2.0) * np.pi * Ky / (lperiod * gamma * beta)
                r[2, 2] = np.cos(omega_x * z)
                r[2, 3] = np.sin(omega_x * z) / omega_x
                r[3, 2] = -np.sin(omega_x * z) * omega_x
                r[3, 3] = np.cos(omega_x * z)

                r[4, 5] = - z / (gamma * beta) ** 2 * (1 + 0.5 * (Kx * beta) ** 2)

            else:
                r[2, 3] = z
            return r

        r_z_e = lambda z, energy: undulator_r_z(z, lperiod=self.lperiod, Kx=self.Kx, Ky=self.Ky, energy=energy)
        # b_z = lambda z, energy: dot((eye(6) - R_z(z, energy)), array([dx, 0., dy, 0., 0., 0.]))
        return r_z_e


class Cavity(Element):
    """
    Standing wave RF cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    vx_{up/down}, vy_{up/down} - zero order kick of a {up/down}stream coupler
    vxx_{up/down}, vxy_{up/down} - first order kick  a {up/down}stream coupler
    """

    default_tm = CavityTM
    additional_tms = []

    def __init__(self, l=0., v=0., phi=0., freq=0., vx_up=0, vy_up=0, vxx_up=0, vxy_up=0,
                 vx_down=0, vy_down=0, vxx_down=0, vxy_down=0, eid=None):
        Element.__init__(self, eid)
        self.l = l
        self.v = v  # in GV
        self.freq = freq  # Hz
        self.phi = phi  # in grad
        self.E = 0
        self.vx_up = vx_up
        self.vy_up = vy_up
        self.vxx_up = vxx_up
        self.vxy_up = vxy_up
        self.vx_down = vx_down
        self.vy_down = vy_down
        self.vxx_down = vxx_down
        self.vxy_down = vxy_down

    def __str__(self):
        s = 'Cavity : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l    =%8.4f m\n' % self.l
        s += 'v    =%8.5f GV\n' % self.v
        s += 'freq =%8.1e Hz\n' % self.freq
        s += 'phi  =%8.2f deg\n' % self.phi
        s += "\nCoupler kick: \n"
        s += "vx_up    = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vx_up)
        s += "vy_up    = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vy_up)
        s += "vxx_up   = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vxx_up)
        s += "vxy_up   = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vxy_up)
        s += "vx_down  = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vx_down)
        s += "vy_down  = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vy_down)
        s += "vxx_down = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vxx_down)
        s += "vxy_down = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vxy_down)
        return s

    def create_r_matrix(self):

        def cavity_R_z(z, V, E, freq, phi=0.):
            """
            :param z: length
            :param de: delta E
            :param freq: frequency
            :param E: initial energy
            :return: matrix
            """

            phi = phi * np.pi / 180.
            de = V * np.cos(phi)
            # pure pi-standing-wave case
            eta = 1
            # gamma = (E + 0.5 * de) / m_e_GeV
            Ei = E / m_e_GeV
            Ef = (E + de) / m_e_GeV
            Ep = (Ef - Ei) / z  # energy derivative
            if Ei == 0:
                logger.error("CAVITY: Initial energy is 0, check ParticleArray.E or Twiss.E OR cavity.v must be 0")

            cos_phi = np.cos(phi)
            alpha = np.sqrt(eta / 8.) / cos_phi * np.log(Ef / Ei)
            sin_alpha = np.sin(alpha)

            cos_alpha = np.cos(alpha)
            r11 = (cos_alpha - np.sqrt(2. / eta) * cos_phi * sin_alpha)

            if abs(Ep) > 1e-10:
                r12 = np.sqrt(8. / eta) * Ei / Ep * cos_phi * sin_alpha
            else:
                r12 = z
            r21 = -Ep / Ef * (cos_phi / np.sqrt(2. * eta) + np.sqrt(eta / 8.) / cos_phi) * sin_alpha

            r22 = Ei / Ef * (cos_alpha + np.sqrt(2. / eta) * cos_phi * sin_alpha)
            # print(f"z = {z}, V = {V}, E = {E}, phi = {phi}, alpha = {alpha}, r11 = {r11}")

            r56 = 0.
            beta0 = 1
            beta1 = 1

            k = 2. * np.pi * freq / speed_of_light
            r55_cor = 0.
            if V != 0 and E != 0:
                gamma2 = Ei * Ei
                beta0 = np.sqrt(1. - 1 / gamma2)
                gamma2 = Ef * Ef
                beta1 = np.sqrt(1. - 1 / gamma2)

                # r56 = (beta0 / beta1 - 1) * Ei / (Ef - Ei) * z
                r56 = - z / (Ef * Ef * Ei * beta1) * (Ef + Ei) / (beta1 + beta0)
                g0 = Ei
                g1 = Ef
                r55_cor = k * z * beta0 * V / m_e_GeV * np.sin(phi) * (g0 * g1 * (beta0 * beta1 - 1) + 1) / (
                        beta1 * g1 * (g0 - g1) ** 2)

            r66 = Ei / Ef * beta0 / beta1
            r65 = k * np.sin(phi) * V / (Ef * beta1 * m_e_GeV)
            cav_matrix = np.array([[r11, r12, 0., 0., 0., 0.],
                                   [r21, r22, 0., 0., 0., 0.],
                                   [0., 0., r11, r12, 0., 0.],
                                   [0., 0., r21, r22, 0., 0.],
                                   [0., 0., 0., 0., 1. + r55_cor, r56],
                                   [0., 0., 0., 0., r65, r66]]).real

            return cav_matrix

        if self.v == 0.:
            r_z_e = lambda z, energy: uni_matrix(z, 0., hx=0., sum_tilts=self.dtilt + self.tilt, energy=energy)
        else:
            r_z_e = lambda z, energy: cavity_R_z(z, V=self.v * z / self.l, E=energy, freq=self.freq,
                                                 phi=self.phi)
        return r_z_e


class CouplerKick(Element):
    """
    Coupler Kick element for Cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    vx, vy - zero order kick of a stream coupler
    vxx, vxy - first order kick  a stream coupler
    """

    default_tm = CouplerKickTM
    additional_tms = []

    def __init__(self, v=0., phi=0., freq=0., vx=0., vy=0., vxx=0., vxy=0., eid=None):
        Element.__init__(self, eid)
        self.l = 0.
        self.v = v  # in GV
        self.freq = freq  # Hz
        self.phi = phi  # in grad
        self.vx = vx
        self.vy = vy
        self.vxx = vxx
        self.vxy = vxy

    def __str__(self):
        s = 'CouplerKick : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'v    =%8.5f GV\n' % self.v
        s += 'freq =%8.1e Hz\n' % self.freq
        s += 'phi  =%8.2f deg\n' % self.phi
        s += "vx   = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vx)
        s += "vy   = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vy)
        s += "vxx  = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vxx)
        s += "vxy  = {num.real:+9.2e} {num.imag:+9.2e}j\n".format(num=self.vxy)
        return s

    def create_r_matrix(self):

        def ck_matrix(v, phi, vxx, vxy, energy):
            """
            matrix for coupler kick

            :param v: voltage of the cavity in GV
            :param phi: phase [deg] of the cavity
            :param vxx: first order coefficients of the coupler kicks
            :param vxy: first order coefficients of the coupler kicks
            :param energy: beam energy in GeV
            :return:
            """
            phi = phi * np.pi / 180.
            m21 = (vxx * v * np.exp(1j * phi)).real / energy
            m43 = - m21
            m23 = (vxy * v * np.exp(1j * phi)).real / energy

            coupl_kick = np.array([[1, 0., 0., 0., 0., 0.],
                                    [m21, 1, m23, 0., 0., 0.],
                                    [0., 0., 1, 0., 0., 0.],
                                    [m23, 0., m43, 1, 0., 0.],
                                    [0., 0., 0., 0., 1., 0.],
                                    [0., 0., 0., 0., 0., 1]])
            return coupl_kick

        r_z_e = lambda z, energy: ck_matrix(v=self.v, phi=self.phi,
                                            vxx=self.vxx, vxy=self.vxy, energy=energy)
        return r_z_e


class TWCavity(Element):
    """
    Traveling wave cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    """

    default_tm = TWCavityTM
    additional_tms = [SecondTM]

    def __init__(self, l=0., v=0., phi=0., freq=0., eid=None):
        Element.__init__(self, eid)
        self.l = l
        self.v = v  # in GV
        self.freq = freq  # Hz
        self.phi = phi  # in grad
        self.E = 0

    def create_r_matrix(self):

        def tw_cavity_R_z(z, V, E, freq, phi=0.):
            """
            :param z: length
            :param de: delta E
            :param f: frequency
            :param E: initial energy
            :return: matrix
            """
            phi = phi * np.pi / 180.
            de = V * np.cos(phi)
            r12 = z * E / de * np.log(1. + de / E) if de != 0 else z
            r22 = E / (E + de)
            r65 = V * np.sin(phi) / (E + de) * (2 * np.pi / (speed_of_light / freq)) if freq != 0 else 0
            r66 = r22
            cav_matrix = np.array([[1, r12, 0., 0., 0., 0.],
                                   [0, r22, 0., 0., 0., 0.],
                                   [0., 0., 1, r12, 0., 0.],
                                   [0., 0., 0, r22, 0., 0.],
                                   [0., 0., 0., 0., 1., 0],
                                   [0., 0., 0., 0., r65, r66]]).real
            return cav_matrix

        def f_entrance(z, V, E, phi=0.):
            phi = phi * np.pi / 180.
            de = V * np.cos(phi)
            r = np.eye(6)
            r[1, 0] = -de / z / 2. / E
            r[3, 2] = r[1, 0]
            return r

        def f_exit(z, V, E, phi=0.):
            phi = phi * np.pi / 180.
            de = V * np.cos(phi)
            r = np.eye(6)
            r[1, 0] = +de / z / 2. / (E + de)
            r[3, 2] = r[1, 0]
            return r

        def cav(z, V, E, freq, phi):
            R_z = np.dot(tw_cavity_R_z(z, V, E, freq, phi), f_entrance(z, V, E, phi))
            R = np.dot(f_exit(z, V, E, phi), R_z)
            return R

        if self.v == 0.:
            r_z_e = lambda z, energy: uni_matrix(z, 0., hx=0., sum_tilts=self.dtilt + self.tilt, energy=energy)
        else:
            r_z_e = lambda z, energy: cav(z, V=self.v * z / self.l, E=energy, freq=self.freq,
                                          phi=self.phi)
        return r_z_e


class TDCavity(Element):
    """
    Transverse deflecting cavity - by default kick in horizontal plane

    l - length [m]
    v - voltage [GV/m]
    freq - frequency [Hz]
    phi - phase in [deg]
    tilt - tilt of cavity in [rad]
    """
    default_tm = TransferMap
    additional_tms = [SecondTM]

    def __init__(self, l=0., freq=0.0, phi=0.0, v=0., tilt=0.0, eid=None):
        Element.__init__(self, eid)
        self.l = l
        self.v = v  # in GV
        self.freq = freq  # Hz
        self.phi = phi  # in deg
        self.tilt = tilt

    def __str__(self):
        s = 'TDCavity : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l    =%8.4f m\n' % self.l
        s += 'v    =%8.5f GV\n' % self.v
        s += 'freq =%8.1e Hz\n' % self.freq
        s += 'phi  =%8.2f deg\n' % self.phi
        s += 'tilt =%8.2f deg\n' % (self.tilt * 180.0 / np.pi)
        return s

    def create_r_matrix(self):
        """
         R - matrix for TDS - NOT TESTED
         """

        def tds_R_z(z, energy, freq, v, phi):
            """

            :param z:  length [m]
            :param freq: freq [Hz]
            :param v: voltage in [GeV]
            :param phi: phase [deg]
            :param energy: Energy in [GeV]
            :return:
            """
            phi = phi * np.pi / 180.

            gamma = energy / m_e_GeV
            igamma2 = 0.
            k0 = 2 * np.pi * freq / speed_of_light
            if gamma != 0:
                igamma2 = 1. / (gamma * gamma)
            if gamma > 1:
                pref = m_e_GeV * np.sqrt(gamma ** 2 - 1)
                K = v * k0 / pref
            else:
                K = 0.
            cos_phi = np.cos(phi)
            cos2_phi = np.cos(2 * phi)

            rm = np.eye(6)

            rm[0, 1] = z
            rm[0, 4] = -z * K * cos_phi / 2.
            rm[1, 4] = -K * cos_phi
            rm[2, 3] = z
            rm[4, 5] = - z * igamma2 / (1. - igamma2)
            rm[5, 0] = rm[1, 4]
            rm[5, 1] = rm[0, 4]
            rm[5, 4] = -z * K ** 2 * cos2_phi / 6
            return rm

        r_z_e = lambda z, energy: tds_R_z(z, energy, freq=self.freq, v=self.v * z / self.l, phi=self.phi)
        return r_z_e


class Solenoid(Element):
    """
    Solenoid
    l - length in m,
    k - strength B0/(2B*rho)
    """

    def __init__(self, l=0., k=0., eid=None):
        Element.__init__(self, eid)
        self.k = k  # B0/(2B*rho)
        self.l = l

    def __str__(self):
        s = 'Cavity : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l =%8.4f m\n' % self.l
        s += 'k =%8.3f 1/m\n' % self.k
        return s

    def create_r_matrix(self):

        def sol(l, k, energy):
            """
            K.Brown, A.Chao.
            :param l: effective length of solenoid
            :param k: B0/(2*Brho), B0 is field inside the solenoid, Brho is momentum of central trajectory
            :return: matrix
            """
            gamma = energy / m_e_GeV
            c = np.cos(l * k)
            s = np.sin(l * k)
            if k == 0:
                s_k = l
            else:
                s_k = s / k
            r56 = 0.
            if gamma != 0:
                gamma2 = gamma * gamma
                beta = np.sqrt(1. - 1. / gamma2)
                r56 -= l / (beta * beta * gamma2)
            sol_matrix = np.array([[c * c, c * s_k, s * c, s * s_k, 0., 0.],
                                   [-k * s * c, c * c, -k * s * s, s * c, 0., 0.],
                                   [-s * c, -s * s_k, c * c, c * s_k, 0., 0.],
                                   [k * s * s, -s * c, -k * s * c, c * c, 0., 0.],
                                   [0., 0., 0., 0., 1., r56],
                                   [0., 0., 0., 0., 0., 1.]]).real
            return sol_matrix

        r_z_e = lambda z, energy: sol(z, k=self.k, energy=energy)
        return r_z_e


class Multipole(Element):
    """
    kn - list of strengths
    """

    default_tm = MultipoleTM
    additional_tms = []

    def __init__(self, kn=0., eid=None):
        Element.__init__(self, eid)
        kn = np.array([kn]).flatten()
        if len(kn) < 2:
            self.kn = np.append(kn, [0.])
        else:
            self.kn = kn
        self.n = len(self.kn)
        self.l = 0.

    def __str__(self):
        s = 'Multipole : '
        s += 'id = ' + str(self.id) + '\n'
        for i, k in enumerate(self.kn):
            s += 'k%i =%8.4f m\n' % (i, k)
        return s

    def create_r_matrix(self):
        r = np.eye(6)
        r[1, 0] = -self.kn[1]
        r[3, 2] = self.kn[1]
        r[1, 5] = self.kn[0]
        r_z_e = lambda z, energy: r
        return r_z_e


class Matrix(Element):
    """
    Matrix element

    l = 0 - m, length of the matrix element
    r = np.zeros((6, 6)) - R - elements, first order
    t = np.zeros((6, 6, 6)) - T - elements, second order
    delta_e = 0 - GeV, energy gain along the matrix element
    """

    def __init__(self, l=0., delta_e=0, eid=None, **kwargs):
        Element.__init__(self, eid)
        self.l = l

        self.r = np.zeros((6, 6))
        self.t = np.zeros((6, 6, 6))
        # zero order elements - test mode, not implemented yet
        self.b = np.zeros((6, 1))

        for y in kwargs:
            # decode first order arguments in format RXX or rXX where X is number from 1 to 6
            if "r" in y[0].lower() and len(y) > 2:
                if "m" in y[1].lower() and len(y) == 4 and y[2:].isdigit() and (11 <= int(y[2:]) <= 66):
                    self.r[int(y[2]) - 1, int(y[3]) - 1] = float(kwargs[y])
                if len(y) == 3 and y[1:].isdigit() and (11 <= int(y[1:]) <= 66):
                    self.r[int(y[1]) - 1, int(y[2]) - 1] = float(kwargs[y])

            # decode second order arguments in format TXXX or tXXX where X is number from 1 to 6
            if "t" in y[0].lower() and len(y) == 4 and y[1:].isdigit() and (111 <= int(y[1:]) <= 666):
                self.t[int(y[1]) - 1, int(y[2]) - 1, int(y[3]) - 1] = float(kwargs[y])

            # decode zero order arguments in format BX or bX where X is number from 1 to 6
            if "b" in y[0].lower() and len(y) == 2 and y[1:].isdigit() and (1 <= int(y[1:]) <= 6):
                self.b[int(y[1]) - 1, 0] = float(kwargs[y])
        self.delta_e = delta_e

    def __str__(self):
        s = 'Matrix : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l =%8.5f m\n' % self.l
        s += 'R = \n'
        for i in range(6):
            for j in range(6):
                s += '%11.6f' % (self.r[i, j])
            s += "\n"
        return s

    def create_r_matrix(self):
        rm = np.eye(6)
        rm = self.r

        def r_matrix(z, l, rm):
            if z < l:
                r_z = uni_matrix(z, 0, hx=0)
            else:
                r_z = rm
            return r_z

        r_z_e = lambda z, energy: r_matrix(z, self.l, rm)
        return r_z_e

    def get_T_z_e_func(self):
        return lambda z, energy: self.t

    def _extend_element_def_string(self, params):
        for key in self.__dict__:
            if isinstance(self.__dict__[key], np.ndarray):
                # r - elements
                if np.shape(self.__dict__[key]) == (6, 6):
                    for i in range(6):
                        for j in range(6):
                            val = self.__dict__[key][i, j]
                            if np.abs(val) > 1e-9:
                                params.append(key + str(i + 1) + str(j + 1) + '=' + str(val))
                # t - elements
                elif np.shape(self.__dict__[key]) == (6, 6, 6):
                    for i in range(6):
                        for j in range(6):
                            for k in range(6):
                                val = self.__dict__[key][i, j, k]
                                if np.abs(val) > 1e-9:
                                    params.append(key + str(i + 1) + str(j + 1) + str(k + 1) + '=' + str(val))
                # b - elements
                if np.shape(self.__dict__[key]) == (6, 1):
                    for i in range(6):
                        val = self.__dict__[key][i, 0]
                        if np.abs(val) > 1e-9:
                            params.append(key + str(i + 1) + '=' + str(val))
        return params

    def _set_general_tm_parameter(self):
        self.transfer_map.delta_e = self.delta_e
        self.transfer_map.B_z = lambda z, energy: self.b
        self.transfer_map.B = lambda energy: self.b
        super()._set_general_tm_parameter()


class Pulse:
    def __init__(self):
        self.kick_x = lambda tau: 0.0
        self.kick_y = lambda tau: 0.0
        self.kick_z = lambda tau: 0.0


class UnknownElement(Element):
    """
    l - length of lens in [m]
    """

    def __init__(self, l=0, kick=0, xsize=0, ysize=0, volt=0, lag=0, harmon=0, refer=0, vkick=0, hkick=0, eid=None):
        Element.__init__(self, eid)
        self.l = l


class Sequence:
    def __init__(self, l=0, refer=0):
        self.l = l


def survey(lat, ang=0.0, x0=0, z0=0):
    x = []
    z = []
    for e in lat.sequence:
        x.append(x0)
        z.append(z0)
        if e.__class__ in [Bend, SBend, RBend]:
            ang += e.angle * 0.5
            s = 2 * e.l * np.sin(e.angle * 0.5) / e.angle
            x0 += s * np.cos(ang)
            z0 += s * np.sin(ang)
            ang += e.angle * 0.5
        else:
            x0 += e.l * np.cos(ang)
            z0 += e.l * np.sin(ang)
    return x, z, ang


if __name__ == "__main__":
    a = RBend(l=13)
