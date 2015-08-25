__author__ = 'Sergey'

from ocelot.cpbd.beam import Particle, Twiss, ParticleArray
from ocelot.common.globals import *
from copy import copy
from numpy.linalg import inv
from numpy import cosh, sinh
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


def uni_matrix(z, k1, hx, sum_tilts=0.):
    # r = element.l/element.angle
    # - K - focusing lens , +K - defoc
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
    x, px, y, py, tau, dp = X[0::6], X[1::6],X[2::6], X[3::6], X[4::6], X[5::6]

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
        self.order = order
        self.dx = 0.
        self.dy = 0.
        self.tilt = 0.
        self.length = 0
        self.k1 = 0.
        self.k2 = 0.
        self.hx = 0.
        # 6x6 linear transfer matrix
        self.R = eye(6)
        self.T = zeros((6, 6, 6))
        self.R_z = lambda z: zeros((6, 6))
        self.B = zeros(6)  # tmp matrix
        # self.map = lambda x: x
        # TODO: implement polynomial transfer maps
        if order > 1:
            pass

    def map_x_twiss(self, tws0):
    
        M = self.R
        m = tws0
        tws = Twiss(tws0)
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

        d_mux = arctan(M[0, 1]/(M[0, 0]*m.beta_x - M[0, 1]*m.alpha_x))
        if d_mux < 0:
            d_mux += pi
        tws.mux = m.mux + d_mux

        d_muy = arctan(M[2, 3]/(M[2, 2]*m.beta_y - M[2, 3]*m.alpha_y))
        if d_muy < 0:
            d_muy += pi
        tws.muy = m.muy + d_muy

        return tws

    def mul_p_array(self, particles, order):
        # particles = pa.particles
        if self.order == 1 and order == 1:
            n = len(particles)
            a = np.add(np.transpose(dot(self.R, np.transpose(particles.reshape(n/6, 6)))), self.B).reshape(n)
            particles[:] = a[:]
        else:
            self.map(particles)
            #self.symp_kick2(particles)
        # pa.s += self.length
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
            m2.R = dot(self.R, m.R)
            m2.B = dot(self.R, m.B) + self.B  #+dB #check
            m2.length = m.length + self.length
            if self.order > 1:
                pass
            return m2

        elif m.__class__ == Particle:
            p = Particle()
            X0 = array([m.x, m.px, m.y, m.py, m.tau, m.p])
            p.x, p.px, p.y, p.py, p.tau, p.p = self.mul_p_array(X0, order = 1)
            p.s = m.s + self.length
            #p = self.map_x_particle(m)
            return p

        elif m.__class__ == Twiss:
            tws = self.map_x_twiss(m)
            # trajectory
            X0 = array([m.x, m.xp, m.y, m.yp, m.tau, m.p])
            tws.x, tws.xp, tws.y, tws.yp, tws.tau, tws.dE = self.mul_p_array(X0, order = 1)
            tws.s = m.s + self.length
            """
            p_tmp = Particle(x=m.x, y=m.y, px=m.xp, py=m.yp, s=m.s, p=m.p,  tau=m.tau)
            p = self.map_x_particle(p_tmp)
            tws.x, tws.xp, tws.y, tws.yp, tws.tau, tws.dE = p.x,p.px, p.y, p.py, p.tau, p.p
            tws.s = p.s
            """
            return tws

        else:
            print(m.__class__)
            exit("unknown object in transfer map multiplication (TransferMap.__mul__)")

    def apply(self, prcl_series, order = 1):

        if prcl_series.__class__ == list and prcl_series[0].__class__ == Particle:
            pa = ParticleArray()
            pa.list2array(prcl_series)
            self.mul_p_array(pa.particles, order)
            pa.s += self.length
            pa.array2ex_list(prcl_series)
            """
            n = len(prcl_series)
            if self.order == 1 and order == 1:

                B = self.B #+ dB
                for i in range(n):
                    p = prcl_series[i]
                    V = array([p.x, p.px, p.y, p.py, p.tau, p.p])
                    p.x, p.px, p.y, p.py, p.tau, p.p = dot(self.R, V) + B
                    p.s = p.s + self.length
            else:
                for i in range(n):
                    p = prcl_series[i]
                    V = array([p.x, p.px, p.y, p.py, p.tau, p.p])
                    p.x, p.px, p.y, p.py, p.tau, p.p = self.map(V)
                    p.s += self.length
            """
        elif prcl_series.__class__ == ParticleArray:

            self.mul_p_array(prcl_series.particles, order)
            prcl_series.s += self.length

    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = m.R_z(s)
        m.B = m.B_z(s)
        m.T = m.T_z(s)
        m.map = lambda u: m.map_z(u, s)
        return m


def create_transfer_map(element, order=1, energy=0, track_acceleration=False):
    transfer_map = TransferMap()
    transfer_map.length = element.l
    transfer_map.dx = element.dx
    transfer_map.dy = element.dy
    transfer_map.tilt = element.dtilt + element.tilt
    transfer_map.energy = energy
    transfer_map.k1 = element.k1
    transfer_map.k2 = element.k2
    # transfer_map.map = None

    if element.l == 0:
        transfer_map.hx = 0.
    else:
        transfer_map.hx = element.angle/element.l

    transfer_map.R_z = lambda z: uni_matrix(z, element.k1, hx = transfer_map.hx, sum_tilts=transfer_map.tilt)
    transfer_map.R = transfer_map.R_z(element.l)

    transfer_map.T_z = lambda z: t_nnn(z, transfer_map.hx, element.k1, element.k2)
    transfer_map.T = transfer_map.T_z(element.l)

    transfer_map.B_z = lambda z: dot((eye(6) - transfer_map.R_z(z)), array([element.dx, 0., element.dy, 0., 0., 0.]))
    transfer_map.B = dot((eye(6) - transfer_map.R), array([element.dx, 0., element.dy, 0., 0., 0.]))

    r_z = lambda z: uni_matrix(z, element.k1, hx = transfer_map.hx, sum_tilts=0)

    transfer_map.map_z = lambda X, z: t_apply(r_z(z), transfer_map.T_z(z), X, element.dx, element.dy, transfer_map.tilt)
    transfer_map.map = lambda X: transfer_map.map_z(X, element.l)

    #def symp_kick2(X, transfer_map, ndivs = 1):
    #    h = transfer_map.hx
    #    if h != 0:
    #        ndivs = 100
    #    k1 = transfer_map.k1
    #    if k1 !=0:
    #        ndivs = 30
    #    k2 = transfer_map.k2
    #    beta = 1.
    #    gamma2_inv = 0.
    #    L = transfer_map.length/ndivs
    #    #print h, k1, k2, L
    #    #R_2 = transfer_map.R_z(L/2.)
    #    n = len(X)
    #    #Xt = transpose(dot(R, transpose(X.reshape(n/6, 6)))).reshape(n)
    #    #X[:] = Xt[:]
    #    Xn = zeros(n)
    #    Xn[:] = X[:]
    #    #print R_2
    #    n = len(X)
    #    for i in range(ndivs):
    #        #Xr = transpose(dot(R_2, transpose(Xn.reshape(n/6,6)))).reshape(n)
    #        Xt = zeros(n)
    #        #x, px, y, py, tau, dp = Xr[0::6], Xr[1::6], Xr[2::6], Xr[3::6], Xr[4::6], Xr[5::6]
    #        x, px, y, py, tau, dp = Xn[0::6], Xn[1::6], Xn[2::6], Xn[3::6], Xn[4::6], Xn[5::6]
    #        px2 = px*px
    #        py2 = py*py
    #
    #        Xt[0::6] = x + L*px*(1. + h*x - dp/beta)
    #        Xt[1::6] = px + L*0.5*(-h*(px2 + py2) - 2*(h*h + k1)*x - (2*h*k1 + k2)*x*x + (h*k1+k2)*y*y + (2*h *dp)/beta)
    #
    #        Xt[2::6] = y + L*py*(1. + h*x - dp/beta)
    #        Xt[3::6] = py + L*((k1 + h*k1*x + k2*x)*y)
    #
    #        Xt[4::6] = tau + L*(-(px2 + py2)/(2*beta) - gamma2_inv/(beta*(1+beta)))
    #        Xt[5::6] = dp
    #        Xn[:] = Xt[:]
    #        #Xn = transpose(dot(R_2, transpose(Xt.reshape(n/6, 6)))).reshape(n)
    #    X[:] = Xn[:]
    #    return X
    #transfer_map.symp_kick2 = lambda X: symp_kick2(X, transfer_map.R)
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
        #R, T = fringe_ent(h=element.h, k1=element.k1,  e=element.edge, h_pole=element.h_pole, gap=element.gap, fint=element.fint)
        if element.pos == 1:
            R, T = fringe_ent(h=element.h, k1=element.k1,  e=element.edge, h_pole=element.h_pole, gap=element.gap, fint=element.fint)
        else:
            R, T = fringe_ext(h=element.h, k1=element.k1,  e=element.edge, h_pole=element.h_pole, gap=element.gap, fint=element.fint)
        transfer_map.R = dot(dot(rot_mtx(-tilt), R), rot_mtx(tilt))
        transfer_map.R_z = lambda z: transfer_map.R
        transfer_map.T = T
        transfer_map.T_z = lambda z: transfer_map.T
        transfer_map.map_z = lambda X, z: t_apply(R, transfer_map.T_z(z), X, element.dx, element.dy, element.tilt)
        transfer_map.map = lambda X: transfer_map.map_z(X, element.l)
        #print element.edge, element.pos
        #print T
    elif element.type == "sextupole":

        def nonl_kick_array_z(u, z, ms):

            z1 = z/2.
            x = u[0::6] + u[1::6]*z1 - transfer_map.dx
            y = u[2::6] + u[3::6]*z1 - transfer_map.dy

            u[1::6] += -ms/2.*(x*x - y*y)
            u[3::6] += x*y*ms

            u[0::6] = x + u[1::6]*z1 + transfer_map.dx
            u[2::6] = y + u[3::6]*z1 + transfer_map.dy

            return u

        if element.l == 0:
            transfer_map.nonl_kick_z = lambda u, z: nonl_kick_array_z(u, z, element.ms)
            transfer_map.nonl_kick = lambda u: nonl_kick_array_z(u, element.l, element.ms)

            transfer_map.ms = element.ms
            #transfer_map.symp_kick2 = transfer_map.nonl_kick
            transfer_map.T[1, 0, 0] = -element.ms/2.
            transfer_map.T[1, 2, 2] = element.ms/2.
            transfer_map.T[3, 0, 2] = element.ms
            transfer_map.T_z = lambda z: transfer_map.T

        transfer_map.R_z = lambda z: uni_matrix(z, 0., hx=0.)
        transfer_map.R = transfer_map.R_z(element.l)
        transfer_map.map_z = lambda X, z: t_apply(transfer_map.R_z(z), transfer_map.T_z(z), X, element.dx, element.dy, element.tilt)
        transfer_map.map = lambda X: transfer_map.map_z(X, element.l)
        transfer_map.order = 2
        #print transfer_map.T

    elif element.type == "octupole":

        def nonl_kick_array_z(u, z, moct):
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
            transfer_map.map_z = lambda u, z: nonl_kick_array_z(u, z, element.moct)
            transfer_map.map = lambda u: nonl_kick_array_z(u, element.l, element.moct)

        else:

            transfer_map.map_z = lambda u, z: nonl_kick_array_z(u, z, element.k3*z)
            transfer_map.map = lambda u: nonl_kick_array_z(u, element.l, element.moct)

        transfer_map.R_z = lambda z: uni_matrix(z, 0., hx = 0.)
        transfer_map.R = transfer_map.R_z(element.l)
        transfer_map.T_z = lambda z: t_nnn(z, h=0., k1=0., k2=0.)
        transfer_map.T = transfer_map.T_z(element.l)

    elif element.type == "undulator":
        #transfer_map.T_z = lambda z: t_nnn(z, h=0., k1=0., k2=0.)
        #transfer_map.T = transfer_map.T_z(element.l)

        def undulator_R_z(z, lperiod, Kx, Ky, energy):
            #print energy
            gamma = energy / m_e_GeV

            beta = 1 / sqrt(1.0-1.0/(gamma*gamma))

            omega_x = sqrt(2.0) * pi * Kx / (lperiod * gamma*beta)
            omega_y = sqrt(2.0) * pi * Ky / (lperiod * gamma*beta)

            R = eye(6)
            R[0,1] = z
            R[2,2] = cos(omega_x * z )
            R[2,3] = sin(omega_x * z ) / omega_x
            R[3,2] = -sin(omega_x * z ) * omega_x
            R[3,3] = cos(omega_x * z )
            return R

        def nonl_kick_rk( u):

            z = np.linspace(0, element.l, 10000)
            #x, px,y, py = track_in_undul(u[:4], z, kz, kx, rho)
            x, px,y, py, z, pz = track_und_RK(u, z, kz, kx ,element.Kx, energy)
            u[0] = x[-1]
            u[1] = px[-1]
            u[2] = y[-1]
            u[3] = py[-1]
            return u

        def nonl_kick_array_z(u, z, ndiv):

            if element.solver == "rk":
                # It requires scipy.weave !!!
                track_und_openmp(u, z, 5000, kz, kx ,element.Kx, energy)

            elif element.solver == "rk_func":
                rk_field(u, z, 5000, energy, element.mag_field)

            else:
                zi = linspace(0., z, num=ndiv)
                h = zi[1] - zi[0]
                kx2 = kx*kx
                kz2 = kz*kz
                ky2 = kz*kz + kx*kx
                ky = sqrt(ky2)
                gamma = energy*1957.
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

        #def nonl_kick_array_z(u, z, ndiv):
        #    #N = 10000
        #    track_und_openmp(u, element.l, 10000, kz, kx ,element.Kx, energy)
        #    return u

        if energy == 0 or element.lperiod == 0 or element.Kx == 0:

            transfer_map.order = 1
            transfer_map.R_z = lambda z: uni_matrix(z, 0, hx = 0, sum_tilts = element.dtilt + element.tilt)
            transfer_map.R = transfer_map.R_z(element.l)
            transfer_map.map_z = lambda X, z: t_apply(transfer_map.R_z(z), transfer_map.T_z(z), X, 0, 0, 0)
            transfer_map.map = lambda X: transfer_map.map_z(X, element.l)

        else:
            try:
                element.solver
            except:
                element.solver = "sym"
            transfer_map.order = 2

            transfer_map.R_z = lambda z: undulator_R_z(z, lperiod = element.lperiod, Kx = element.Kx, Ky = element.Ky, energy =energy)
            transfer_map.R = transfer_map.R_z(element.l)

            kz = 2.*pi/element.lperiod
            if element.ax == -1:
                kx = 0
            else:
                kx = 2.*pi/element.ax
            transfer_map.map = lambda u: nonl_kick_array_z(u, element.l, ndiv = 20)
            transfer_map.map_z = lambda u, z: nonl_kick_array_z(u, z, ndiv = 5)

    elif element.type == "hcor" or element.type == "vcor":

        def kick_b(z,l,angle_x, angle_y):
            if l == 0:
                hx = 0.; hy = 0.
            else:
                hx = angle_x/l; hy = angle_y/l

            dx = hx*z*z/2.
            dy = hy*z*z/2.
            dx1 = hx*z if l != 0 else angle_x
            dy1 = hy*z if l != 0 else angle_y
            b = array([dx, dx1, dy, dy1, 0., 0.])
            #ux = 0 if hx != 0 else 1.
            #uy = 0 if hy != 0 else 1.
            # b = array([(1 - cos(z*hx))/(hx+ux), sin(z*hx)+ux*angle_x, (1 - cos(z*hy))/(hy+uy), sin(z*hy)+uy*angle_y, 0, 0])
            return b

        transfer_map.order = 0

        # bend(element,transfer_map)

        if element.type == "hcor":
            transfer_map.B_z = lambda z: kick_b(z, element.l, element.angle, 0)
            transfer_map.B = transfer_map.B_z(element.l)
        else:
            transfer_map.B_z = lambda z: kick_b(z, element.l, 0, element.angle)
            transfer_map.B = transfer_map.B_z(element.l)

        transfer_map.R_z = lambda z: uni_matrix(z, 0., hx = 0.)
        transfer_map.R = transfer_map.R_z(element.l)
        transfer_map.T_z = lambda z: t_nnn(z, h=0., k1=0., k2=0.)
        transfer_map.T = transfer_map.T_z(element.l)

        transfer_map.map_z = lambda X, z: dot(transfer_map.R_z(z), X) + transfer_map.B_z(z)
        transfer_map.map = lambda X: transfer_map.map_z(X, element.l)

    elif element.type == "cavity":
                
        def cavity_R_z(z, de, f, E, track_acceleration):
            """
            :param z: length
            :param de: delta E
            :param f: frequency
            :param E: initial energy
            :return: matrix
            """
            eta = 1.0
            phi = 0
            
            #Ep = de * cos(phi) / (z * 0.000511) # energy derivative
            de = de
            Ep = de / (z) # energy derivative
            #Ep = de  # energy derivative
            #Ef = E + de
            #Ei = E
            Ei = E - de
            Ef = E 
            
            #print Ei,Ef,Ep, z, de
            
            alpha = sqrt(eta / 8.) / cos(phi) * log(Ef/Ei)
            
            #print 'cavity map Ei=',E, 'Ef=', Ef, 'alpha=', alpha
            
            r11 = cos(alpha) - sqrt(2./eta) * cos(phi) * sin(alpha)
            r12 = sqrt(8./eta) * Ei / Ep * cos(phi) * sin(alpha)
            r21 = -Ep/Ei * (cos(phi)/ sqrt(2.*eta) + sqrt(eta/8.) / cos(phi) ) * sin(alpha)
            r22 = Ei/Ef * ( cos(alpha) + sqrt(2./eta) * cos(phi) * sin(alpha) )
            
            # for twiss parameters
            if not track_acceleration:
                r11=r11*sqrt(Ef/Ei)
                r12=r12*sqrt(Ef/Ei)
                r21=r21*sqrt(Ef/Ei)
                r22=r22*sqrt(Ef/Ei)
            
            #print r11, r12, r21, r22
            
            cav_matrix = array([[r11, r12, 0. ,0. ,0., 0.],
                                [r21, r22, 0. ,0. ,0., 0.],
                                [0., 0., r11, r12, 0., 0.],
                                [0., 0., r21, r22, 0., 0.],
                                [0., 0., 0., 0. ,1. ,0.],
                                [0., 0., 0., 0. ,0. ,1.]]).real
    
            return cav_matrix
                
        if element.E == 0 or (element.v < 1.e-10 and element.delta_e < 1.e-10 ):
            #print 'Unit CAVITY MAP:', element.E, 'GeV'
            transfer_map.R_z = lambda z: uni_matrix(z, 0, hx = 0, sum_tilts = element.dtilt + element.tilt)
            transfer_map.R = transfer_map.R_z(element.l)
        else:
            #print 'CAVITY MAP:', element.E, 'GeV'
            #print element.E, element.v, element.delta_e 
            transfer_map.R_z = lambda z: cavity_R_z(z, de = element.delta_e * z / element.l, f=element.f, E=element.E, track_acceleration = track_acceleration)
            transfer_map.R = transfer_map.R_z(element.l)

        transfer_map.T_z = lambda z: t_nnn(z, h=0., k1=0., k2=0.)
        transfer_map.T = transfer_map.T_z(element.l)
        transfer_map.map_z = lambda X, z: t_apply(transfer_map.R_z(z), transfer_map.T_z(z), X, element.dx, element.dy, element.tilt)
        transfer_map.map = lambda X: transfer_map.map_z(X, element.l)

    elif element.type == "solenoid":
        def sol(l, k):
            """
            K.Brown, A.Chao.
            :param l: efective length of solenoid
            :param k: B0/(2*Brho), B0 is field inside the solenoid, Brho is momentum of central trajectory
            :return: matrix
            """
            c = cos(l*k)
            s = sin(l*k)
            if k == 0:
                s_k = l
            else:
                s_k = s/k

            sol_matrix = array([[c*c, c*s_k, s*c, s*s_k, 0., 0.],
                                [-k*s*c, c*c, -k*s*s, s*c, 0., 0.],
                                [-s*c, -s*s_k, c*c, c*s_k, 0., 0.],
                                [k*s*s, -s*c, -k*s*c, c*c, 0., 0.],
                                [0., 0., 0., 0., 1., 0.],
                                [0., 0., 0., 0., 0., 1.]]).real
            return sol_matrix
        T = zeros((6,6,6))
        transfer_map.R_z = lambda z: sol(z, k = element.k)
        transfer_map.R = transfer_map.R_z(element.l)
        transfer_map.map_z = lambda X, z: t_apply(transfer_map.R_z(z), T, X, element.dx, element.dy, element.tilt)
        transfer_map.map = lambda X: transfer_map.map_z(X, element.l)

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

        transfer_map.R = Rm

        def r_matrix(z, l, Rm):
            if z < l:
                R_z = uni_matrix(z, 0, hx = 0)
            else:
                R_z = Rm
            return R_z

        transfer_map.R_z = lambda z: r_matrix(z, element.l, Rm)
        transfer_map.T_z = lambda z: t_nnn(z, h=0., k1=0., k2=0.)
        transfer_map.T = transfer_map.T_z(element.l)
        #transfer_map.map_z = lambda X: t_apply(transfer_map.R, transfer_map.T, X, element.dx, element.dy, element.tilt)
        transfer_map.map_z = lambda X, z: t_apply(transfer_map.R_z(z), transfer_map.T_z(z), X, element.dx, element.dy, element.tilt)
        transfer_map.map = lambda X: transfer_map.map_z(X, element.l)

    else:
        print (element.type, " : unknown type of magnetic element. Cannot create transfer map ")
    #if element.type != "sextupole":
    #    transfer_map.symp_kick2 = lambda X: symp_kick2(X, transfer_map)
    return transfer_map



def periodic_solution(tws, transfer_matrix):

    """ find periodical twiss  """

    tws = Twiss(tws)
    M = transfer_matrix

    cosmx = (M[0, 0] + M[1, 1])/2.
    cosmy = (M[2, 2] + M[3, 3])/2.

    # print cosmx, cosmy

    if abs(cosmx) >= 1 or abs(cosmy) >= 1:
        print("************ periodic solution does not exist. return None ***********")
        return None
    sinmx = sqrt(1.-cosmx*cosmx)
    sinmy = sqrt(1.-cosmy*cosmy)
    tws.beta_x = abs(M[0, 1]/sinmx)
    tws.beta_y = abs(M[2, 3]/sinmy)

    tws.alpha_x = (M[0, 0] - M[1, 1])/(2*sinmx)  # X[0,0]
    tws.gamma_x = (1. + tws.alpha_x*tws.alpha_x)/tws.beta_x  # X[1,0]

    tws.alpha_y = (M[2, 2] - M[3, 3])/(2*sinmy)  # Y[0,0]
    tws.gamma_y = (1. + tws.alpha_y*tws.alpha_y)/tws.beta_y  # Y[1,0]

    Hx = array([[M[0, 0] - 1, M[0, 1]], [M[1, 0], M[1, 1]-1]])
    Hhx = array([[M[0, 5]], [M[1, 5]]])
    hh = dot(inv(-Hx), Hhx)
    tws.Dx = hh[0, 0]
    tws.Dxp = hh[1, 0]

    Hy = array([[M[2, 2] - 1, M[2, 3]], [M[3, 2], M[3, 3]-1]])
    Hhy = array([[M[2, 5]], [M[3, 5]]])
    hhy = dot(inv(-Hy), Hhy)
    tws.Dy = hhy[0, 0]
    tws.Dyp = hhy[1, 0]
    return tws


def lattice_transfer_map(lattice):
    """ transfer map for the whole lattice"""
    lat_transfer_map = TransferMap()
    for elem in lattice.sequence:
        #print elem.type
        lat_transfer_map = elem.transfer_map*lat_transfer_map
    return lat_transfer_map


def trace_z(lattice, obj0, z_array):
    """ Z-dependent tracer (twiss(z) and particle(z))
        usage: twiss = trace_z(lattice,twiss_0, [1.23, 2.56, ...]) ,
        to calculate Twiss params at 1.23m, 2.56m etc.
    """
    obj_list = []
    try:
        E0 = obj0.E
    except:
        E0 = 0
    i = 0
    elem = lattice.sequence[i]
    L = elem.l
    obj_elem = obj0
    for z in z_array:
        while z > L:
            if elem.type == "cavity":
                E0 += elem.delta_e
            obj_elem.E = E0
            elem.E = E0
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
        
        E0 = obj.E
        
        for e in lattice.sequence:
            #print e.type, e.id
            
            obj = e.transfer_map*obj
                        
            if e.type == "cavity":
                E0 += e.delta_e
                
            obj.E = E0
            e.E = E0
            obj.id = e.id


            obj_list.append(obj)
    else:
        z_array = linspace(0, lattice.totalLen, nPoints, endpoint=True)
        obj_list =  trace_z(lattice, obj, z_array)
    return obj_list


def twiss(lattice, tws0, nPoints = None):

    if tws0.__class__ == Twiss:
        if tws0.beta_x == 0  or tws0.beta_y == 0:
            tws0 = periodic_solution(tws0, lattice_transfer_map(lattice).R)
            if tws0 == None:
                return None
        else:
            tws0.gamma_x = (1. + tws0.alpha_x**2)/tws0.beta_x
            tws0.gamma_y = (1. + tws0.alpha_y**2)/tws0.beta_y

        twiss_list = trace_obj(lattice, tws0, nPoints)
    else:
        exit("unknown object tws0")
    return twiss_list


#def trace_obj(lattice, p0, nPoints = None):
#    '''
#    track a particle using generic 'trace_obj' function
#    '''
#    if p0.__class__ == Particle:
#        p_list = trace_obj(lattice, p0, nPoints)
#    else:
#        exit("unknown object p0")
#    return p_list


class Navigator:
    def __init__(self, lattice = None):
        if lattice != None:
            self.lat = lattice
        
    z0 = 0
    n_elem = 0
    sum_lengths = 0

    #def check(self, dz):
    #    '''
    #    check if next step exceed the bounds of lattice
    #    '''
    #    if self.z0+dz>self.lat.totalLen:
    #        dz = self.lat.totalLen - self.z0
    #    return dz

def get_map(lattice, dz, navi, order=1):
    TM = []
    tm = TransferMap(identity=True)
    i = navi.n_elem
    z1 = navi.z0 + dz
    elem = lattice.sequence[i]
    L = navi.sum_lengths + elem.l
    
    delta_e = 0.0
    phi=0.0
    freq=0.0

    while z1 > L:
        dl = L - navi.z0
        #print dl, L, navi.z0, elem.l
        if elem.transfer_map.order > 1 or order > 1:
            if tm.identity == False:
                TM.append(tm)
            TM.append(elem.transfer_map(dl))
            tm = TransferMap(identity=True)
        else:
            tm = elem.transfer_map(dl)*tm
            if elem.type == "cavity":
                delta_e += elem.delta_e * dl / elem.l
                phi=elem.phi
                freq=elem.f

        navi.z0 = L
        dz -= dl
        i += 1
        elem = lattice.sequence[i]
        L += elem.l

    if elem.transfer_map.order > 1 or order > 1:
        TM.append(elem.transfer_map(dz))
    else:
        tm = elem.transfer_map(dz)*tm
        if elem.type == "cavity":
            delta_e += elem.delta_e * dz / elem.l
            phi=elem.phi
            freq=elem.f

    navi.z0 += dz
    navi.sum_lengths = L - elem.l
    navi.n_elem = i
    if tm.identity == False:
        TM.append(tm)
    return TM, delta_e, phi, freq


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
