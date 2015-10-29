__author__ = 'Sergey'

from ocelot.cpbd.beam import Particle, Twiss, ParticleArray
from ocelot.common.globals import *
from copy import copy
from numpy.linalg import inv
from numpy import cosh, sinh
from scipy.misc import factorial
# from scipy import weave
from und_weave import *
from ocelot.cpbd.maps2order import *


def rot_mtx(angle):
    return array([[cos(angle), 0., sin(angle), 0., 0., 0.],
                  [0., cos(angle), 0., sin(angle), 0., 0.],
                  [-sin(angle), 0., cos(angle), 0., 0., 0.],
                  [0., -sin(angle), 0., cos(angle), 0., 0.],
                  [0., 0., 0., 0., 1., 0.],
                  [0., 0., 0., 0., 0., 1.]])


def uni_matrix(z, k1, hx, sum_tilts=0., energy=0.):
    # r = element.l/element.angle
    # - K - focusing lens , +K - defoc
    gamma = energy/m_e_GeV

    kx2 = (k1 + hx*hx)
    ky2 = -k1
    kx = sqrt(kx2 + 0.j)
    ky = sqrt(ky2 + 0.j)
    cx = cos(z*kx).real
    cy = cos(z*ky).real
    sy = (sin(ky*z)/ky).real if ky != 0 else z
    if kx != 0:
        sx = (sin(kx*z)/kx).real
        dx = hx/kx2*(1. - cx)
        r56 = hx*hx*(z - sx)/kx2
    else:
        sx = z
        dx = z*z*hx/2.
        r56 = hx*hx*z**3/6.
    if gamma != 0:
        r56 -= z/(gamma*gamma)

    u_matrix = array([[cx, sx, 0., 0., 0., dx],
                        [-kx2*sx, cx, 0., 0., 0., sx*hx],
                        [0., 0., cy, sy, 0., 0.],
                        [0., 0., -ky2*sy, cy, 0., 0.],
                        [hx*sx, dx, 0., 0., 1., r56],
                        [0., 0., 0., 0., 0., 1.]])
    u_matrix = dot(dot(rot_mtx(-sum_tilts), u_matrix), rot_mtx(sum_tilts))
    return u_matrix


def transform_vec_ent(X, dx, dy, tilt):
    n = len(X)
    for i in range(n/6):
        X0 = X[6*i:6*(i+1)]
        X0 -= array([dx, 0., dy, 0., 0., 0.])
        X[6*i:6*(i+1)] = dot(rot_mtx(tilt), X0)
    return X

def transform_vec_ext(X, dx, dy, tilt):
    n = len(X)
    for i in range(n/6):
        X0 = X[6*i:6*(i+1)]
        X[6*i:6*(i+1)] = dot(rot_mtx(-tilt), X0)
        X0 += array([dx, 0., dy, 0., 0., 0.])
    return X


def t_apply(R, T, X, dx, dy, tilt):

    if dx != 0 or dy != 0 or tilt != 0:
        X = transform_vec_ent(X, dx, dy, tilt)

    n = len(X)
    Xr = transpose(dot(R, transpose(X.reshape(n/6, 6)))).reshape(n)

    Xt = zeros(n)
    x, px, y, py, tau, dp = X[0::6], X[1::6], X[2::6], X[3::6], X[4::6], X[5::6]

    Xt[0::6] = T[0, 0, 0]*x*x + T[0, 0, 1]*x*px + T[0, 0, 5]*x*dp + T[0, 1, 1]*px*px + T[0, 1, 5]*px*dp + \
               T[0, 5, 5]*dp*dp + T[0, 2, 2]*y*y + T[0, 2, 3]*y*py + T[0, 3, 3]*py*py

    Xt[1::6] = T[1, 0, 0]*x*x + T[1, 0, 1]*x*px + T[1, 0, 5]*x*dp + T[1, 1, 1]*px*px + T[1, 1, 5]*px*dp + \
               T[1, 5, 5]*dp*dp + T[1, 2, 2]*y*y + T[1, 2, 3]*y*py + T[1, 3, 3]*py*py

    Xt[2::6] = T[2, 0, 2]*x*y + T[2, 0, 3]*x*py + T[2, 1, 2]*px*y + T[2, 1, 3]*px*py + T[2, 2, 5]*y*dp + T[2, 3, 5]*py*dp

    Xt[3::6] = T[3, 0, 2]*x*y + T[3, 0, 3]*x*py + T[3, 1, 2]*px*y + T[3, 1, 3]*px*py + T[3, 2, 5]*y*dp + T[3, 3, 5]*py*dp

    Xt[4::6] = T[4, 0, 0]*x*x + T[4, 0, 1]*x*px + T[4, 0, 5]*x*dp + T[4, 1, 1]*px*px + T[4, 1, 5]*px*dp + \
               T[4, 5, 5]*dp*dp + T[4, 2, 2]*y*y + T[4, 2, 3]*y*py + T[4, 3, 3]*py*py

    X[:] = Xr[:] + Xt[:]

    if dx != 0 or dy != 0 or tilt != 0:
        X = transform_vec_ext(X, dx, dy, tilt)

    return X


class TransferMap:

    def __init__(self, order=1, identity=False):
        self.identity = identity
        self.type = "none"
        self.order = order
        self.dx = 0.
        self.dy = 0.
        self.tilt = 0.
        self.length = 0
        #self.energy = 0.
        self.k1 = 0.
        self.k2 = 0.
        self.hx = 0.
        # test RF
        self.delta_e = 0.0
        self.delta_e_z = lambda z: 0.0
        self.map = lambda x, energy: x
        # 6x6 linear transfer matrix
        self.R = lambda energy: eye(6)
        self.T = zeros((6, 6, 6))
        self.R_z = lambda z, energy: zeros((6, 6))
        self.B = lambda energy: zeros(6)  # tmp matrix
        # self.map = lambda x: x
        # TODO: implement polynomial transfer maps
        if order > 1:
            pass

    def map_x_twiss(self, tws0):
        E = tws0.E
        M = self.R(E)
        if self.delta_e != 0:
            #M = self.R(E + )
            Ei = tws0.E
            Ef = tws0.E + self.delta_e #* cos(self.phi)
            #print "Ei = ", Ei, "Ef = ", Ef
            k = sqrt(Ef/Ei)
            #k = 1.
            M[0, 0] = M[0, 0]*k
            M[0, 1] = M[0, 1]*k
            M[1, 0] = M[1, 0]*k
            M[1, 1] = M[1, 1]*k
            M[2, 2] = M[2, 2]*k
            M[2, 3] = M[2, 3]*k
            M[3, 2] = M[3, 2]*k
            M[3, 3] = M[3, 3]*k
            E = Ef
        m = tws0
        tws = Twiss(tws0)
        tws.E = E
        tws.p = m.p
        tws.beta_x = M[0, 0]*M[0, 0]*m.beta_x - 2*M[0, 1]*M[0, 0]*m.alpha_x + M[0, 1]*M[0, 1]*m.gamma_x
        # tws.beta_x = ((M[0,0]*tws.beta_x - M[0,1]*m.alpha_x)**2 + M[0,1]*M[0,1])/m.beta_x
        tws.beta_y = M[2, 2]*M[2, 2]*m.beta_y - 2*M[2, 3]*M[2, 2]*m.alpha_y + M[2, 3]*M[2, 3]*m.gamma_y
        # tws.beta_y = ((M[2,2]*tws.beta_y - M[2,3]*m.alpha_y)**2 + M[2,3]*M[2,3])/m.beta_y
        tws.alpha_x = -M[0, 0]*M[1, 0]*m.beta_x + (M[0, 1]*M[1, 0]+M[1, 1]*M[0, 0])*m.alpha_x - M[0, 1]*M[1, 1]*m.gamma_x
        tws.alpha_y = -M[2, 2]*M[3, 2]*m.beta_y + (M[2, 3]*M[3, 2]+M[3, 3]*M[2, 2])*m.alpha_y - M[2, 3]*M[3, 3]*m.gamma_y
    
        tws.gamma_x = (1. + tws.alpha_x*tws.alpha_x)/tws.beta_x
        tws.gamma_y = (1. + tws.alpha_y*tws.alpha_y)/tws.beta_y
    
        tws.Dx = M[0, 0]*m.Dx + M[0, 1]*m.Dxp + M[0, 5]
        tws.Dy = M[2, 2]*m.Dy + M[2, 3]*m.Dyp + M[2, 5]
    
        tws.Dxp = M[1, 0]*m.Dx + M[1, 1]*m.Dxp + M[1, 5]
        tws.Dyp = M[3, 2]*m.Dy + M[3, 3]*m.Dyp + M[3, 5]
        denom_x = M[0, 0]*m.beta_x - M[0, 1]*m.alpha_x
        if denom_x == 0.:
            d_mux = pi/2.*M[0, 1]/np.abs(M[0, 1])
        else:
            d_mux = arctan(M[0, 1]/denom_x)

        if d_mux < 0:
            d_mux += pi
        tws.mux = m.mux + d_mux
        #print M[0, 0]*m.beta_x - M[0, 1]*m.alpha_x, arctan(M[2, 3]/(M[2, 2]*m.beta_y - M[2, 3]*m.alpha_y))
        denom_y = M[2, 2]*m.beta_y - M[2, 3]*m.alpha_y
        if denom_y == 0.:
            d_muy = pi/2.*M[2, 3]/np.abs(M[2, 3])
        else:
            d_muy = arctan(M[2, 3]/denom_y)
        if d_muy < 0:
            d_muy += pi
        tws.muy = m.muy + d_muy
        return tws

    def mul_p_array(self, particles, energy=0., order=1):

        if self.order == 1 and order == 1:
            n = len(particles)
            #print self.R(energy)
            a = np.add(np.transpose(dot(self.R(energy), np.transpose(particles.reshape(n/6, 6)))), self.B(energy)).reshape(n)
            particles[:] = a[:]
        elif order == 3:
            self.sym_map(particles, energy=energy)

        else:
            #print " Brown "
            self.map(particles, energy=energy)
        return particles

    def __mul__(self, m):
        """
        :param m: TransferMap, Particle or Twiss
        :return: TransferMap, Particle or Twiss
        Ma = {Ba, Ra, Ta}
        Mb = {Bb, Rb, Tb}
        X1 = R*(X0 - dX) + dX = R*X0 + B
        B = (E - R)*dX
        """

        if m.__class__ == TransferMap:
            m2 = TransferMap()
            #print m.type
            m2.R = lambda energy: dot(self.R(energy), m.R(energy))
            m2.B = lambda energy: dot(self.R(energy), m.B(energy)) + self.B(energy)  #+dB #check
            m2.length = m.length + self.length
            #m2.delta_e += self.delta_e

            if self.order > 1:
                pass
            return m2

        elif m.__class__ == Particle:
            p = Particle()
            X0 = array([m.x, m.px, m.y, m.py, m.tau, m.p])
            p.x, p.px, p.y, p.py, p.tau, p.p = self.mul_p_array(X0, order=1)
            p.s = m.s + self.length
            #p = self.map_x_particle(m)
            return p

        elif m.__class__ == Twiss:
            tws = self.map_x_twiss(m)
            # trajectory
            X0 = array([m.x, m.xp, m.y, m.yp, m.tau, m.p])
            tws.x, tws.xp, tws.y, tws.yp, tws.tau, tws.dE = self.mul_p_array(X0, energy=tws.E, order=1)
            tws.s = m.s + self.length
            return tws

        else:
            print(m.__class__)
            exit("unknown object in transfer map multiplication (TransferMap.__mul__)")

    def apply(self, prcl_series, order = 1):

        if prcl_series.__class__ == list and prcl_series[0].__class__ == Particle:
            list_e = array([p.E for p in prcl_series])
            if False in (list_e[:] == list_e[0]):
                for p in prcl_series:
                    self.mul_p_array(array([p.x, p.px, p.y, p.py, p.tau, p.p]), energy=p.E, order=order)
                    p.s += self.length
            else:

                pa = ParticleArray()
                pa.list2array(prcl_series)
                pa.E = prcl_series[0].E
                self.mul_p_array(pa.particles, energy=pa.E, order=order)
                pa.s += self.length
                pa.array2ex_list(prcl_series)

        elif prcl_series.__class__ == ParticleArray:
            prcl_series.E += self.delta_e
            self.mul_p_array(prcl_series.particles, energy=prcl_series.E, order=order)
            prcl_series.s += self.length

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = lambda energy: m.R_z(s, energy)
        m.B = lambda energy: m.B_z(s, energy)
        m.T = m.T_z(s)
        m.delta_e = m.delta_e_z(s)
        m.map = lambda u, energy: m.map_z(u, s, energy)
        m.sym_map = lambda u, energy: m.sym_map_z(u, s, energy)
        return m


def create_transfer_map(element, order=1):
    transfer_map = TransferMap()
    transfer_map.type = element.type
    transfer_map.length = element.l
    transfer_map.dx = element.dx
    transfer_map.dy = element.dy
    transfer_map.tilt = element.dtilt + element.tilt
    transfer_map.k1 = element.k1
    transfer_map.k2 = element.k2
    if element.l == 0:
        transfer_map.hx = 0.
    else:
        transfer_map.hx = element.angle/element.l
    #def R_z(z, energy):
    #    return uni_matrix(z, element.k1, hx = transfer_map.hx, sum_tilts=0, energy=energy)
    R_z = lambda z, energy: uni_matrix(z, element.k1, hx = transfer_map.hx, sum_tilts=transfer_map.tilt, energy=energy)
    #def b_z(z, energy):
    #    return dot((eye(6) - R_z(z, energy)), array([element.dx, 0., element.dy, 0., 0., 0.]))
    b_z = lambda z, energy: dot((eye(6) - R_z(z, energy)), array([element.dx, 0., element.dy, 0., 0., 0.]))
    transfer_map.T_z = lambda z: t_nnn(z, transfer_map.hx, element.k1, element.k2)
    transfer_map.T = transfer_map.T_z(element.l)

    transfer_map.map_z = lambda X, z, energy: t_apply(R_z(z, energy), transfer_map.T_z(z), X, element.dx, element.dy, transfer_map.tilt)
    #experiment with symplecticity
    transfer_map.sym_map_z = lambda X, z, energy: sym_map(z, X, transfer_map.hx, element.k1, element.k2, energy)

    if element.type == "quadrupole":
        pass

    elif element.type in ["sbend", "rbend", "bend"]:
        pass

    elif element.type == "drift":
        pass

    elif element.type == "monitor":
        pass

    elif element.type == "marker":
        pass

    elif element.type == "edge":
        tilt = element.tilt + element.dtilt
        if element.pos == 1:
            R, T = fringe_ent(h=element.h, k1=element.k1,  e=element.edge, h_pole=element.h_pole, gap=element.gap, fint=element.fint)
        else:
            R, T = fringe_ext(h=element.h, k1=element.k1,  e=element.edge, h_pole=element.h_pole, gap=element.gap, fint=element.fint)
        R_z = lambda z, energy: dot(dot(rot_mtx(-tilt), R), rot_mtx(tilt))

        transfer_map.T = T
        transfer_map.T_z = lambda z: transfer_map.T

        transfer_map.map_z = lambda X, z, energy: t_apply(R, transfer_map.T_z(z), X, element.dx, element.dy, element.tilt)
        #transfer_map.map_z = lambda X, z, energy: t_apply(R, np.zeros((6, 6, 6)), X, element.dx, element.dy, element.tilt)
        transfer_map.sym_map_z = lambda X, z, energy: t_apply(R, np.zeros((6, 6, 6)), X, element.dx, element.dy, element.tilt)

    elif element.type == "sextupole":

        def map4sextupole(u, z, ms, energy):

            z1 = z/2.
            x = u[0::6] + u[1::6]*z1 - transfer_map.dx
            y = u[2::6] + u[3::6]*z1 - transfer_map.dy

            u[1::6] += -ms/2.*(x*x - y*y)
            u[3::6] += x*y*ms

            u[0::6] = x + u[1::6]*z1 + transfer_map.dx
            u[2::6] = y + u[3::6]*z1 + transfer_map.dy

            return u

        R_z = lambda z, energy: uni_matrix(z, 0., hx=0., energy=energy)
        transfer_map.sym_map_z = lambda X, z, energy: map4sextupole(X, z, element.k2*element.l, energy)
        if element.l == 0:

            transfer_map.ms = element.ms

            transfer_map.T[1, 0, 0] = -element.ms/2.
            transfer_map.T[1, 2, 2] = element.ms/2.
            transfer_map.T[3, 0, 2] = element.ms
            transfer_map.T_z = lambda z: transfer_map.T
            transfer_map.sym_map_z = lambda X, z, energy: map4sextupole(X, z, element.ms, energy)
            #transfer_map.sym_map_z = lambda X, z, energy: t_apply(R_z(z, energy), transfer_map.T_z(z), X, element.dx, element.dy, element.tilt)

        transfer_map.map_z = lambda X, z, energy: t_apply(R_z(z, energy), transfer_map.T_z(z), X, element.dx, element.dy, element.tilt)
        transfer_map.order = 2

    elif element.type == "octupole":

        def map4octupole(u, z, moct):
            #TODO: check expressions
            #v = np.array([transfer_map.dx, transfer_map.dy, z, moct])

            z1 = z/2.
            x = u[0::6] + u[1::6]*z1 - transfer_map.dx
            y = u[2::6] + u[3::6]*z1 - transfer_map.dy

            u[1::6] += -moct/2.*(x*x*x - 3.*y*y*x)
            u[3::6] += moct*(3.*y*x*x-y*y*y)

            u[0::6] = x + u[1::6]*z1 + transfer_map.dx
            u[2::6] = y + u[3::6]*z1 + transfer_map.dy
            return u

        if element.moct == None:
            element.moct = element.k3*element.l
        transfer_map.order = 3

        if element.l == 0:
            transfer_map.map_z = lambda u, z, energy: map4octupole(u, z, element.moct)
        else:
            transfer_map.map_z = lambda u, z, energy: map4octupole(u, z, element.k3*z)

        R_z = lambda z, energy: uni_matrix(z, 0., hx = 0., energy=energy)

        transfer_map.T_z = lambda z: t_nnn(z, h=0., k1=0., k2=0.)
        transfer_map.T = transfer_map.T_z(element.l)

    elif element.type == "undulator":
        def undulator_R_z(z, lperiod, Kx, Ky, energy):
            gamma = energy / m_e_GeV
            R = eye(6)
            R[0,1] = z
            if gamma != 0 and lperiod != 0 and Kx != 0:
                beta = 1 / sqrt(1.0-1.0/(gamma*gamma))
                omega_x = sqrt(2.0) * pi * Kx / (lperiod * gamma*beta)
                omega_y = sqrt(2.0) * pi * Ky / (lperiod * gamma*beta)
                R[2,2] = cos(omega_x * z )
                R[2,3] = sin(omega_x * z ) / omega_x
                R[3,2] = -sin(omega_x * z ) * omega_x
                R[3,3] = cos(omega_x * z )
            else:
                R[2,3] = z
            return R

        R_z = lambda z, energy: undulator_R_z(z, lperiod=element.lperiod, Kx=element.Kx, Ky=element.Ky, energy=energy)

        def map4undulator(u, z, kz, kx, energy, ndiv):
            zi = linspace(0., z, num=ndiv)
            h = zi[1] - zi[0]
            kx2 = kx*kx
            kz2 = kz*kz
            ky2 = kz*kz + kx*kx
            ky = sqrt(ky2)
            gamma = energy/m_e_GeV
            h0 = 0.
            if gamma != 0:
                h0 = 1./(gamma/element.Kx/kz)
            h02 = h0*h0
            h = h/(1.+ u[5::6])
            x = u[::6]
            y = u[2::6]
            for z in range(len(zi)-1):
                chx = cosh(kx*x)
                chy = cosh(ky*y)
                shx = sinh(kx*x)
                shy = sinh(ky*y)
                u[1::6] -= h/2.*chx*shx*(kx*ky2*chy*chy + kx2*kx*shy*shy)/(ky2*kz2)*h02
                u[3::6] -= h/2.*chy*shy*(ky2*chx*chx + kx2*shx*shx)/(ky*kz2)*h02
                u[4::6] -= h/2./(1.+u[5::6]) * ((u[1::6]*u[1::6] + u[3::6]*u[3::6]) + chx*chx*chy*chy/(2.*kz2)*h02 + shx*shx*shy*shy*kx2/(2.*ky2*kz2)*h02)
                u[::6] = x + h*u[1::6]
                u[2::6] = y + h*u[3::6]

            return u

        transfer_map.order = 1
        transfer_map.map_z = lambda X, z, energy: t_apply(R_z(z, energy), transfer_map.T_z(z), X, 0., 0., 0.)
        transfer_map.sym_map_z = lambda X, z, energy: t_apply(R_z(z, energy), transfer_map.T_z(z), X, 0., 0., 0.)

        if element.solver in ["sym", "symplectic"]:
            #print "undulator transfer map is symplectic map! "
            transfer_map.order = 2
            kz = 2.*pi/element.lperiod
            if element.ax == -1:
                kx = 0
            else:
                kx = 2.*pi/element.ax
            transfer_map.map_z = lambda u, z, energy: map4undulator(u, z, kz, kx, energy, ndiv=int(z*10+2))
            transfer_map.sym_map_z = lambda u, z, energy: map4undulator(u, z, kz, kx, energy, ndiv=int(z*10+2))

        transfer_map.map_rk = lambda u, z, energy: rk_field(u, z, N=int(z*10./element.lperiod),
                                                                   energy=energy, mag_field=element.mag_field)

    elif element.type == "hcor" or element.type == "vcor":

        def kick_b(z,l,angle_x, angle_y):
            if l == 0:
                hx = 0.
                hy = 0.
            else:
                hx = angle_x/l
                hy = angle_y/l

            dx = hx*z*z/2.
            dy = hy*z*z/2.
            dx1 = hx*z if l != 0 else angle_x
            dy1 = hy*z if l != 0 else angle_y
            b = array([dx, dx1, dy, dy1, 0., 0.])
            return b

        def map4corr(R, B, X):
            n = len(X)
            X1 = np.add(np.transpose(dot(R, np.transpose(X.reshape(n/6, 6)))), B).reshape(n)
            X[:] = X1[:]
            return X

        if element.type == "hcor":
            b_z = lambda z, energy: kick_b(z, element.l, element.angle, 0)
        else:
            b_z = lambda z, energy: kick_b(z, element.l, 0, element.angle)

        R_z = lambda z, energy: uni_matrix(z, 0., hx = 0., energy=energy)
        transfer_map.T_z = lambda z: t_nnn(z, h=0., k1=0., k2=0.)
        transfer_map.T = transfer_map.T_z(element.l)

        transfer_map.map_z = lambda X, z, energy: map4corr(R_z(z, energy), b_z(z, energy), X)

    elif element.type == "cavity":
                
        def cavity_R_z(z, V, f, E, phi=0.):
            """
            :param z: length
            :param de: delta E
            :param f: frequency
            :param E: initial energy
            :return: matrix
            """
            #phi = 0.
            de = V*cos(phi)
            # pure pi-standing-wave case
            #phi = 0.
            eta = 1.0

            gamma = (E + 0.5*de)/m_e_GeV

            Ep = de/z  # energy derivative

            Ei = E
            Ef = E + de
            
            #print Ei,Ef,Ep, z, de
            if Ei == 0:
                print("Warning! Initial energy is zero and cavity.delta_e != 0! Change Ei or cavity.delta_e must be 0" )
            alpha = sqrt(eta / 8.) / cos(phi) * log(Ef/Ei)
            #alpha = sqrt(1./ 8.) * log(Ef/Ei)
            #print 'cavity map Ei=',E, 'Ef=', Ef, 'alpha=', alpha
            
            r11 = cos(alpha) - sqrt(2./eta) * cos(phi) * sin(alpha)
            r12 = sqrt(8./eta) * Ei / Ep * cos(phi) * sin(alpha)
            r21 = -Ep/Ef * (cos(phi)/ sqrt(2.*eta) + sqrt(eta/8.) / cos(phi) ) * sin(alpha)
            r22 = Ei/Ef * ( cos(alpha) + sqrt(2./eta) * cos(phi) * sin(alpha) )

            #r11 = cos(alpha) - sqrt(2.) * sin(alpha)
            #r12 = sqrt(8.) * Ei / Ep * sin(alpha)
            #r21 = -3.*Ep/sqrt(8.)/Ef * sin(alpha)
            #r22 = Ei/Ef * ( cos(alpha) + sqrt(2.) * sin(alpha) )

            #print -Ep/Ei
            # for twiss parameters go to function -> map_x_twiss()
            #if not track_acceleration:
            #    r11 = r11*sqrt(Ef/Ei)
            #    r12 = r12*sqrt(Ef/Ei)
            #    r21 = r21*sqrt(Ef/Ei)
            #    r22 = r22*sqrt(Ef/Ei)
            #print z, "Ei = ", Ei, "Ef = ", Ef, "de = ", de, "l = ", z
            #if z == 0.346:
            #    print r11, r12,r21 ,r22, "Ei = ", Ei, "Ef = ", Ef, "de = ", de, "l = ", z
            r56 = 0.
            if gamma != 0:
                r56 = z/(gamma*gamma)
            cav_matrix = array([[r11, r12, 0., 0., 0., 0.],
                                [r21, r22, 0., 0., 0., 0.],
                                [0., 0., r11, r12, 0., 0.],
                                [0., 0., r21, r22, 0., 0.],
                                [0., 0., 0., 0., 1., r56],
                                [0., 0., 0., 0., 0., 1.]]).real
    
            return cav_matrix

        def map4cav(R, T, X, dx, dy, tilt, E,  V, freq, phi):
            #print E

            X = t_apply(R, T, X, dx, dy, tilt)
            delta_e = V*cos(phi)
            if E + delta_e > 0:
                k = 2.*pi*freq/speed_of_light
                #phi_rad = phi*np.pi/180.
                #v = delta_e/cos(phi_rad)
                X[5::6] = (X[5::6]*E + V*np.cos(X[4::6]*k + phi) - delta_e)/(E + delta_e)
        transfer_map.phi = element.phi
        transfer_map.order = 2
        #if element.v < 1.e-10 and element.delta_e < 1.e-10:
        if element.delta_e == 0.:
            #transfer_map.order = 1
            R_z = lambda z, energy: uni_matrix(z, 0., hx=0., sum_tilts=element.dtilt + element.tilt, energy=energy)
        else:
            R_z = lambda z, energy: cavity_R_z(z, V=element.v*z/element.l, f=element.f, E=energy, phi=element.phi)
        #def R_z(z, energy):
        #    if energy == 0 or (element.v < 1.e-10 and element.delta_e < 1.e-10 ):
        #        r_z = uni_matrix(z, 0., hx=0., sum_tilts=element.dtilt + element.tilt, energy=energy)
        #    else:
        #        de = element.delta_e*z/element.l
        #        r_z = cavity_R_z(z, de=de, f=element.f, E=energy)
        #    return r_z

        transfer_map.delta_e_z = lambda z: element.v*cos(element.phi) * z / element.l
        transfer_map.delta_e = transfer_map.delta_e_z(element.l)

        #transfer_map.T_z = lambda z: t_nnn(z, h=0., k1=0., k2=0.)
        #transfer_map.T = transfer_map.T_z(element.l)
        transfer_map.map_z = lambda X, z, energy: map4cav(R_z(z, energy), transfer_map.T_z(z), X,
                                                 transfer_map.dx, transfer_map.dy, transfer_map.tilt,
                                                 energy,  element.v*z/element.l, element.f, element.phi)

    elif element.type == "solenoid":
        def sol(l, k, energy):
            """
            K.Brown, A.Chao.
            :param l: efective length of solenoid
            :param k: B0/(2*Brho), B0 is field inside the solenoid, Brho is momentum of central trajectory
            :return: matrix
            """
            gamma = energy/m_e_GeV
            c = cos(l*k)
            s = sin(l*k)
            if k == 0:
                s_k = l
            else:
                s_k = s/k
            r56 = 0.
            if gamma != 0:
                r56 = l/(gamma*gamma)
            sol_matrix = array([[c*c, c*s_k, s*c, s*s_k, 0., 0.],
                                [-k*s*c, c*c, -k*s*s, s*c, 0., 0.],
                                [-s*c, -s*s_k, c*c, c*s_k, 0., 0.],
                                [k*s*s, -s*c, -k*s*c, c*c, 0., 0.],
                                [0., 0., 0., 0., 1., r56],
                                [0., 0., 0., 0., 0., 1.]]).real
            return sol_matrix
        R_z = lambda z, energy: sol(z, k=element.k, energy=energy)
        T = zeros((6, 6, 6))
        transfer_map.map_z = lambda X, z, energy: t_apply(R_z(z, energy), T, X, element.dx, element.dy, element.tilt)

    elif element.type == "matrix":
        Rm = eye(6)
        Rm[0,0] = element.rm11
        Rm[0,1] = element.rm12
        Rm[1,0] = element.rm21
        Rm[1,1] = element.rm22
        Rm[2,2] = element.rm33
        Rm[2,3] = element.rm34
        Rm[3,2] = element.rm43
        Rm[3,3] = element.rm44

        def r_matrix(z, l, Rm):
            if z < l:
                R_z = uni_matrix(z, 0, hx=0)
            else:
                R_z = Rm
            return R_z
        R_z = lambda z, energy: r_matrix(z, element.l, Rm)
        transfer_map.T_z = lambda z: t_nnn(z, h=0., k1=0., k2=0.)
        transfer_map.T = transfer_map.T_z(element.l)
        transfer_map.map_z = lambda X, z, energy: t_apply(R_z(z, energy), transfer_map.T_z(z), X, element.dx, element.dy, element.tilt)
        transfer_map.sym_map_z = lambda X, z, energy: transfer_map.map_z(X, z, energy)

    elif element.type == "multipole":
        def kick(X, kn):
            p = -kn[0]*X[5::6] + 0j
            for n in range(1, len(kn)):
                p += kn[n]*(X[0::6] + 1j*X[2::6])**n/factorial(n)
            X[1::6] = X[1::6] - np.real(p)
            X[3::6] = X[3::6] + np.imag(p)
            X[4::6] = X[4::6] + kn[0::6]*X[0]
            #print X[1], X[3]
            return X

        R = np.eye(6)
        R[1, 0] = -element.kn[1]
        R[3, 2] = element.kn[1]
        R[1, 5] = element.kn[0]
        #transfer_map.order = 2
        R_z = lambda z, energy: R
        if element.n > 2:
            transfer_map.order = 2

        transfer_map.map_z = lambda X, z, energy: kick(X, element.kn)
        transfer_map.sym_map_z = lambda X, z, energy: kick(X, element.kn)
    else:
        print (element.type, " : unknown type of magnetic element. Cannot create transfer map ")

    transfer_map.B_z = lambda z, energy: b_z(z, energy)
    transfer_map.B = lambda energy: transfer_map.B_z(element.l, energy)
    transfer_map.R_z = lambda z, energy: R_z(z, energy)
    transfer_map.R = lambda energy: transfer_map.R_z(element.l, energy)
    transfer_map.map = lambda X, energy: transfer_map.map_z(X, element.l, energy)
    transfer_map.sym_map = lambda X, energy: transfer_map.sym_map_z(X, element.l, energy)
    return transfer_map



def periodic_solution(tws, R):

    """ find periodical twiss  """

    tws = Twiss(tws)
    #print R
    cosmx = (R[0, 0] + R[1, 1])/2.
    cosmy = (R[2, 2] + R[3, 3])/2.

    #print cosmx, cosmy

    if abs(cosmx) >= 1 or abs(cosmy) >= 1:
        print("************ periodic solution does not exist. return None ***********")
        return None
    sinmx = sign(R[0, 1])*sqrt(1.-cosmx*cosmx)
    sinmy = sign(R[2, 3])*sqrt(1.-cosmy*cosmy)
    #sinmx2 = sign(R[0, 1])*sqrt(R[0, 1]*R[1, 0] - 0.25*(R[0, 0] - R[1, 1])**2)
    #print "sinmx = ", sinmx, sinmx2,R[0, 1], R[0, 1]*R[1, 0] - 0.25*(R[0, 0] - R[1, 1])**2

    tws.beta_x = abs(R[0, 1]/sinmx)
    tws.beta_y = abs(R[2, 3]/sinmy)

    tws.alpha_x = (R[0, 0] - R[1, 1])/(2.*sinmx)  # X[0,0]
    #print R[0, 0] - R[1, 1], 2.*sinmx
    tws.gamma_x = (1. + tws.alpha_x*tws.alpha_x)/tws.beta_x  # X[1,0]

    tws.alpha_y = (R[2, 2] - R[3, 3])/(2*sinmy)  # Y[0,0]
    tws.gamma_y = (1. + tws.alpha_y*tws.alpha_y)/tws.beta_y  # Y[1,0]

    Hx = array([[R[0, 0] - 1, R[0, 1]], [R[1, 0], R[1, 1]-1]])
    Hhx = array([[R[0, 5]], [R[1, 5]]])
    hh = dot(inv(-Hx), Hhx)
    tws.Dx = hh[0, 0]
    tws.Dxp = hh[1, 0]
    Hy = array([[R[2, 2] - 1, R[2, 3]], [R[3, 2], R[3, 3]-1]])
    Hhy = array([[R[2, 5]], [R[3, 5]]])
    hhy = dot(inv(-Hy), Hhy)
    tws.Dy = hhy[0, 0]
    tws.Dyp = hhy[1, 0]
    #tws.display()
    return tws


def lattice_transfer_map(lattice, energy):
    """ transfer map for the whole lattice"""
    R = np.eye(6)
    for i, elem in enumerate(lattice.sequence):

        R = dot(elem.transfer_map.R(energy), R)
        #print i, len(lattice.sequence), elem.type, elem.transfer_map.R(6)
    return R


def trace_z(lattice, obj0, z_array):
    """ Z-dependent tracer (twiss(z) and particle(z))
        usage: twiss = trace_z(lattice,twiss_0, [1.23, 2.56, ...]) ,
        to calculate Twiss params at 1.23m, 2.56m etc.
    """
    obj_list = []
    i = 0
    elem = lattice.sequence[i]
    L = elem.l
    obj_elem = obj0
    for z in z_array:
        while z > L:
            obj_elem = lattice.sequence[i].transfer_map*obj_elem
            i += 1
            elem = lattice.sequence[i]
            L += elem.l

        obj_z = elem.transfer_map(z - (L - elem.l))*obj_elem

        obj_list.append(obj_z)
    return obj_list


def trace_obj(lattice, obj, nPoints = None):
    """ track object though lattice
        obj must be Twiss or Particle """

    if nPoints == None:
        obj_list = [obj]
        for e in lattice.sequence:
            obj = e.transfer_map*obj
            obj.id = e.id
            obj_list.append(obj)
    else:
        z_array = linspace(0, lattice.totalLen, nPoints, endpoint=True)
        obj_list = trace_z(lattice, obj, z_array)
    return obj_list


def twiss(lattice, tws0=None, nPoints=None):
    if tws0 == None:
        tws0 = periodic_solution(tws0, lattice_transfer_map(lattice, energy=0.))

    if tws0.__class__ == Twiss:
        if tws0.beta_x == 0  or tws0.beta_y == 0:
            tws0 = periodic_solution(tws0, lattice_transfer_map(lattice, tws0.E))
            if tws0 == None:
                return None
        else:
            tws0.gamma_x = (1. + tws0.alpha_x**2)/tws0.beta_x
            tws0.gamma_y = (1. + tws0.alpha_y**2)/tws0.beta_y

        twiss_list = trace_obj(lattice, tws0, nPoints)
    else:
        exit("unknown object tws0")
    return twiss_list


class Navigator:
    def __init__(self, lattice = None):
        if lattice != None:
            self.lat = lattice
        
    z0 = 0.             # current position if navigator
    n_elem = 0          # current number of the element in lattice
    sum_lengths = 0.    # sum_lengths = Sum[lat.sequence[i].l, {i, 0, n_elem-1}]

    #def check(self, dz):
    #    '''
    #    check if next step exceed the bounds of lattice
    #    '''
    #    if self.z0+dz>self.lat.totalLen:
    #        dz = self.lat.totalLen - self.z0
    #    return dz

def get_map(lattice, dz, navi, order=1):
    #for i, elem in enumerate(lattice.sequence):
    #    print i, elem.type, elem.id
    #order = 2
    TM = []
    tm = TransferMap(identity=True)
    i = navi.n_elem
    z1 = navi.z0 + dz
    elem = lattice.sequence[i]
    L = navi.sum_lengths + elem.l
    rec_count = 0  # counter of recursion in R = lambda energy: dot(R1(energy), R2(energy))
    #print "get_map: order = ", order
    while z1 > L:
        dl = L - navi.z0
        if elem.transfer_map.order > 1 or order > 1:
            if tm.identity == False:
                #print "while, identity == False",  i, elem.type, elem.id,  tm.R
                TM.append(tm)
                rec_count = 0
            #print "while, identity == True", i, elem.type, elem.transfer_map(dl).R(6)
            TM.append(elem.transfer_map(dl))
            tm = TransferMap(identity=True)
        else:
            tm = elem.transfer_map(dl)*tm
            rec_count += 1
            if rec_count > 100:
                TM.append(tm)
                tm = TransferMap(identity=True)
                rec_count = 0
        navi.z0 = L
        dz -= dl
        i += 1
        elem = lattice.sequence[i]
        L += elem.l

    if elem.transfer_map.order > 1 or order > 1:
        if tm.identity == False:
            #print "identity == False", i, elem.type, tm.R(6)
            TM.append(tm)
            rec_count = 0
        #print "identity == True", i, elem.type, tm.R(6)
        TM.append(elem.transfer_map(dz))
        rec_count = 0
        tm = TransferMap(identity=True)
    else:
        tm = elem.transfer_map(dz)*tm

    navi.z0 += dz
    navi.sum_lengths = L - elem.l
    navi.n_elem = i
    if tm.identity == False:
        #print "Final, identity == False", i, elem.type, tm.R(6)
        TM.append(tm)
    return TM


def gaussFromTwiss(emit, beta, alpha):
    phi = 2*pi * random.rand()
    u = random.rand()
    a = sqrt(-2*log( (1-u)) * emit)
    x = a * sqrt(beta) * cos(phi)
    xp = -a / sqrt(beta) * ( sin(phi) + alpha * cos(phi) );
    return (x,xp)


'''
returns two solutions for a periodic fodo, given the mean beta
initial betas are at the center of the focusing quad 
'''
def fodo_parameters(betaXmean = 36.0, L=10.0, verbose = False):
    lquad = 0.001
        
    kap1 = np.sqrt ( 1.0/2.0 * ( (betaXmean/L)*(betaXmean/L) + (betaXmean/L) * np.sqrt(-4.0 + (betaXmean/L)*(betaXmean/L))) )    
    kap2 = np.sqrt ( 1.0/2.0 * ( (betaXmean/L)*(betaXmean/L) - (betaXmean/L) * np.sqrt(-4.0 + (betaXmean/L)*(betaXmean/L))) )
    
    k = 1.0 / (lquad * L * kap2)
    
    f = 1.0 / (k*lquad)
    
    kappa = f / L    
    betaMax = np.array(( L * kap1*(kap1+1)/np.sqrt(kap1*kap1-1), L * kap2*(kap2+1)/np.sqrt(kap2*kap2-1)))
    betaMin = np.array(( L * kap1*(kap1-1)/np.sqrt(kap1*kap1-1), L * kap2*(kap2-1)/np.sqrt(kap2*kap2-1)))
    betaMean = np.array(( L * kap2*kap2 / (np.sqrt(kap2*kap2 - 1.0)),  L * kap1*kap1 / (np.sqrt(kap1*kap1 - 1.0)) ))
    k = np.array((1.0 / (lquad * L * kap1), 1.0 / (lquad * L * kap2) ))
    
    if verbose:
        print ('********* calculating fodo parameters *********')
        print ('fodo parameters:')
        print ('k*l=', k*lquad)
        print ('f=', L * kap1, L*kap2)
        print ('kap1=', kap1)
        print ('kap2=', kap2)
        print ('betaMax=', betaMax)
        print ('betaMin=', betaMin)
        print ('betaMean=', betaMean)
        print ('*********                             *********')
    
    return k*lquad, betaMin, betaMax, betaMean
