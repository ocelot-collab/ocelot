"""
definition of magnetic lattice
linear dimensions in [m]
"""

from ocelot.cpbd.field_map import FieldMap
import numpy as np


class Element(object):
    """
    Element is a basic beamline building element
    Accelerator optics elements are subclasses of Element
    Arbitrary set of additional parameters can be attached if necessary
    """

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
        self.dtilt = 0.
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
    type - "rect" (by default) or "ellipt". 
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


class Hcor(RBend):
    """
    horizontal corrector,
    l - length of magnet in [m],
    angle - angle of bend in [rad],
    """

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

    def validate(self):
        pass

        # maybe we will do two functions
        # 1. load data and check magnetic map
        # 2. check all input data (lperiod nperiod ...). something like this we must do for all elements.

        # what do you think about ending poles? We can do several options
        # a) 1/2,-1,1,... -1,1/2
        # b) 1/2,-1,1,... -1,1,-1/2
        # c) 1/4,-3/4,1,-1... -1,3/4,-1/4   I need to check it.

    def __str__(self):
        s = 'Undulator : '
        s += 'id = ' + str(self.id) + '\n'
        s += 'l        =%8.4f m\n' % self.l
        s += 'nperiods =%8.1f \n' % self.nperiods
        s += 'lperiod  =%8.4f m\n' % self.lperiod
        s += 'Kx       =%8.3f \n' % self.Kx
        s += 'Ky       =%8.3f \n' % self.Ky
        return s


class Cavity(Element):
    """
    Standing wave RF cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    vx_{up/down}, vy_{up/down} - zero order kick of a {up/down}stream coupler
    vxx_{up/down}, vxy_{up/down} - first order kick  a {up/down}stream coupler
    """

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


class CouplerKick(Element):
    """
    Coupler Kick element for Cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    vx, vy - zero order kick of a stream coupler
    vxx, vxy - first order kick  a stream coupler
    """

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

class TWCavity(Element):
    """
    Traveling wave cavity
    v - voltage [GV]
    freq - frequency [Hz]
    phi - phase in [deg]
    """

    def __init__(self, l=0., v=0., phi=0., freq=0., eid=None):
        Element.__init__(self, eid)
        self.l = l
        self.v = v  # in GV
        self.freq = freq  # Hz
        self.phi = phi  # in grad
        self.E = 0


class TDCavity(Element):
    """
    Transverse deflecting cavity - by default kick in horizontal plane

    l - length [m]
    v - voltage [GV/m]
    freq - frequency [Hz]
    phi - phase in [deg]
    tilt - tilt of cavity in [rad]
    """

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


class Multipole(Element):
    """
    kn - list of strengths
    """

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
