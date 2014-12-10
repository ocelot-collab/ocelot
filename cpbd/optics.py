__author__ = 'Sergey'

from numpy import sqrt, matrix, cos, sin, log, tan, eye, zeros, pi, array, linspace, dot, abs, random, arctan, sign, cosh, sinh
import numpy as np
from ocelot.cpbd.beam import Beam, Particle, Twiss, ParticleArray
from ocelot.common.globals import *
from numpy.linalg import inv
from scipy.integrate import simps, trapz
from copy import copy
from scipy import weave
from und_weave import track_und_RK, track_und_mag, track_und_sym, track_und_openmp



class TransferMap:

    def __init__(self, order = 1):
        self.order = order
        self.dx = 0
        self.dy = 0
        self.length = 0
        # 6x6 linear transfer matrix
        self.R = eye(6)
        self.R_z = lambda z: zeros((6,6))
        self.B = zeros(6) # tmp matrix

        # TODO: implement polynomial transfer maps
        if order > 1 :
            pass

    def map_x_twiss(self, tws0):
        M = self.R
        m = tws0
        tws = Twiss(tws0)
        tws.p = m.p
        tws.beta_x = M[0,0]*M[0,0]*m.beta_x - 2*M[0,1]*M[0,0]*m.alpha_x + M[0,1]*M[0,1]*m.gamma_x
        #tws.beta_x = ((M[0,0]*tws.beta_x - M[0,1]*m.alpha_x)**2 + M[0,1]*M[0,1])/m.beta_x
        tws.beta_y = M[2,2]*M[2,2]*m.beta_y - 2*M[2,3]*M[2,2]*m.alpha_y + M[2,3]*M[2,3]*m.gamma_y
        #tws.beta_y = ((M[2,2]*tws.beta_y - M[2,3]*m.alpha_y)**2 + M[2,3]*M[2,3])/m.beta_y
        tws.alpha_x = -M[0,0]*M[1,0]*m.beta_x + (M[0,1]*M[1,0]+M[1,1]*M[0,0])*m.alpha_x - M[0,1]*M[1,1]*m.gamma_x
        tws.alpha_y = -M[2,2]*M[3,2]*m.beta_y + (M[2,3]*M[3,2]+M[3,3]*M[2,2])*m.alpha_y - M[2,3]*M[3,3]*m.gamma_y
    
        tws.gamma_x = (1. + tws.alpha_x*tws.alpha_x)/tws.beta_x
        tws.gamma_y = (1. + tws.alpha_y*tws.alpha_y)/tws.beta_y
    
        tws.Dx = M[0,0]*m.Dx + M[0,1]*m.Dxp + M[0,5]
        tws.Dy = M[2,2]*m.Dy + M[2,3]*m.Dyp + M[2,5]
    
        tws.Dxp = M[1,0]*m.Dx + M[1,1]*m.Dxp + M[1,5]
        tws.Dyp = M[3,2]*m.Dy + M[3,3]*m.Dyp + M[3,5]

        d_mux = arctan(M[0,1]/(M[0,0]*m.beta_x - M[0,1]*m.alpha_x))
        if d_mux < 0:
            d_mux += pi
        tws.mux = m.mux + d_mux

        d_muy = arctan(M[2,3]/(M[2,2]*m.beta_y - M[2,3]*m.alpha_y))
        if d_muy < 0:
            d_muy += pi
        tws.muy = m.muy + d_muy

        return tws

    def map_x_particle(self, m, dB):

        p = Particle()
        if self.order <= 1:
            V = dot(self.R, array([m.x, m.px, m.y, m.py, m.tau, m.p])) + self.B + dB
            p.x, p.px, p.y, p.py, p.tau, p.p = V[0], V[1], V[2], V[3], V[4], V[5]
        elif self.order == 2:
            V = array([m.x, m.px, m.y, m.py, m.tau, m.p])
            #print V
            p.x, p.px, p.y, p.py, p.tau, p.p = self.nonl_kick(V)
        p.s = m.s + self.length
        return p

    def __mul__(self, m):

        M = self.R
        dx = self.dx
        dy = self.dy
        dB = array([(M[0,0]-1)*dx + M[0,2]*dy, M[1,0]*dx + M[1,2]*dy, M[2,0]*dx + (M[2,2]-1)*dy, M[3,0]*dx + M[3,2]*dy, M[4,0]*dx + M[4,2]*dy, M[5,0]*dx + M[5,2]*dy])

        if self.order == 0:
            dB = self.b(self.length)

        if m.__class__ == TransferMap:
            m2 = TransferMap()
            m2.R = dot(M , m.R)
            m2.B = dot(M,m.B) + dB + self.B # check!
            m2.length = m.length + self.length
            if self.order > 1:
                pass
            return m2

        elif m.__class__ == Particle:
            p = self.map_x_particle(m, dB)
            return p

        elif m.__class__ == Twiss:
            tws = self.map_x_twiss(m)
            # trajectory
            p_tmp = Particle(x=m.x, y=m.y, px=m.xp, py=m.yp, s=m.s, p=m.p,  tau=m.tau)
            p = self.map_x_particle(p_tmp, dB)
            tws.x, tws.xp, tws.y, tws.yp, tws.tau, tws.dE = p.x,p.px, p.y, p.py, p.tau, p.p
            tws.s = p.s
            return tws

        else:
            print m.__class__
            exit("unknow object in transfer map multiplication (TransferMap.__mul__)")
            
    def apply(self, prcl_series):
        M = self.R
        dx = self.dx
        dy = self.dy
        dB = array([(M[0,0]-1)*dx + M[0,2]*dy, M[1,0]*dx + M[1,2]*dy, M[2,0]*dx + (M[2,2]-1)*dy, M[3,0]*dx + M[3,2]*dy, M[4,0]*dx + M[4,2]*dy, M[5,0]*dx + M[5,2]*dy])

        if prcl_series.__class__ == list and prcl_series[0].__class__== Particle:
            n = len(prcl_series)
            if self.order == 1:
                B = self.B + dB
                for i in xrange(n):
                    p = prcl_series[i]
                    V = array([p.x, p.px, p.y, p.py, p.tau, p.p])
                    p.x, p.px, p.y, p.py, p.tau, p.p = dot(self.R, V) + B
                    p.s = p.s + self.length
            else:
                for i in xrange(n):
                    p = prcl_series[i]
                    V = array([p.x, p.px, p.y, p.py, p.tau, p.p])
                    p.x, p.px, p.y, p.py, p.tau, p.p = self.nonl_kick(V)
                    p.s += self.length

        elif prcl_series.__class__ == ParticleArray:
            particles = prcl_series.particles

            def multy_vect(v, r, b):
                #code = """
                #int n,i,j;
                #for(n = 0; n<Nv[0]/6;n++){
                #    double v2[6] = {V1(n*6+0), V1(n*6+1), V1(n*6+2), V1(n*6+3), V1(n*6+4), V1(n*6+5)};
                #    for(i = 0;i<6;i++){
                #        double tmp = 0.;
                #        for( j =0;j<6;j++){
                #            tmp += R2(i,j)*v2[j];
                #            }
                #        V1(n*6 + i) = tmp + B1(i);
                #        }
                #    }
                #"""
                #weave.inline(code, ["v", "r", "b"])
                n = len(v)
                a = np.add(np.transpose(dot(r, np.transpose(v.reshape(n/6,6)) ) ),b).reshape(n)
                v[:]=a[:]

                return v

            if self.order == 1:
                B = self.B + dB
                particles = multy_vect(v= particles, r = self.R, b = B)

                """
                for i in xrange(len(particles)/6):
                    V = particles[i*6:i*6+6]
                    particles[i*6:i*6+6] = dot(self.R, V)+B
                """
            else:
                self.nonl_kick_array(particles)
                """
                for i in xrange(len(particles)/6):
                    V = particles[i*6:i*6+6]
                    particles[i*6:i*6+6] = self.nonl_kick(V)
                """
            prcl_series.s += self.length

        
            
    def __call__(self, s):
        m = copy(self)
        m.length = s
        m.R = m.R_z(s)
        m.nonl_kick = lambda u: m.nonl_kick_z(u, s)
        m.nonl_kick_array = lambda u: m.nonl_kick_array_z(u, s)
        return m


def create_transfer_map(element, order=1, energy = 0):
    #print 'creating TM', element.id
    transfer_map = TransferMap()
    transfer_map.length = element.l
    transfer_map.dx = element.dx
    transfer_map.dy = element.dy
    transfer_map.energy = energy

    def rot_mtx(angle):
        return array([[cos(angle), 0., sin(angle), 0., 0., 0.],
                            [0., cos(angle), 0., sin(angle), 0., 0.],
                            [-sin(angle), 0., cos(angle), 0., 0., 0.],
                            [0., -sin(angle), 0., cos(angle), 0., 0.],
                            [0., 0., 0., 0., 1., 0.],
                            [0., 0., 0., 0., 0., 1.]])

    def uni_matrix(z, k1, hx, hy = 0, sum_tilts = 0.):
        #r = element.l/element.angle
        # - K - focusing lens , +K - defoc
        dp = 0.
        Kx = (k1 + hx*hx)/(1.+dp)
        Ky = -k1/(1.+dp)

        Qx = sqrt(Kx + 0.j)
        Qy = sqrt(Ky + 0.j)
        ux = 0
        uy = 0
        if Ky == 0:
            uy = 1
        if Kx == 0:
            ux = 1
        cos_x = cos(z*Qx)
        sin_x = sin(z*Qx)
        cos_z = cos(z*Qy)
        sin_z = sin(z*Qy)

        uni_matrix = array([[ cos_x,  sin_x/(Qx+ux) + ux*z, 0. ,0. ,0.,  (1. - cos_x)*hx/(Kx+ux) ],
                            [-Qx*sin_x, cos_x, 0. ,0. ,0., sin_x*hx/(Qx+ux) ],
                            [0., 0., cos_z,  sin_z/(Qy+uy) + uy*z,  0., 0. ],
                            [0., 0., -Qy*sin_z, cos_z, 0., 0.],
                            [hx/(Qx+ux)*sin(z*Qx), hx/(Qx*Qx + ux)*(1. - cos(z*Qx)), 0., 0. ,1. ,hx*hx/(Qx**3+ux)*(Qx*z - sin(Qx*z)) ],
                            [0., 0., 0., 0. ,0. ,1. ]]).real

        uni_matrix = dot(dot(rot_mtx(-sum_tilts),uni_matrix),rot_mtx(sum_tilts))
        return uni_matrix

    def bend(element,transfer_map):
        if element.l == 0:
            hx = 0.
        else:
            hx = element.angle/element.l

        transfer_map.R_z = lambda z: uni_matrix(z, element.k1, hx = hx, sum_tilts = element.dtilt + element.tilt)#lambda z: transfer_map.R1(z)*transfer_map.e_start
        transfer_map.R = transfer_map.R_z(element.l)#transfer_map.e_end*transfer_map.R_z(element.l)


    if element.type == "quadrupole":

        transfer_map.R_z = lambda z: uni_matrix(z, element.k1, hx = 0, sum_tilts = element.dtilt + element.tilt)
        transfer_map.R = transfer_map.R_z(element.l)

    elif element.type == "sbend" or element.type == "rbend" or element.type == "bend":

        bend(element,transfer_map)

    elif element.type == "edge":
        #print 'Edge', element.id, element.edge
        y = tan(element.edge)*element.h
        tilt = element.tilt + element.dtilt
        m1 = y*cos(tilt)*cos(tilt) - y*sin(tilt)*sin(tilt)
        m2 = y*sin(2*tilt)

        transfer_map.R[1,0] = m1
        transfer_map.R[1,2] = m2
        transfer_map.R[3,1] = m2
        transfer_map.R[3,2] = -m1
        transfer_map.R_z = lambda z: transfer_map.R

    elif element.type == "sextupole":

        def nonl_kick_z(u, z, ms):
            #ms = k2*z
            z1 = z/2.
            x = u[0] - transfer_map.dx + u[1]*z1
            y = u[2] - transfer_map.dy + u[3]*z1

            u[1] += -ms/2.*(x*x - y*y)
            u[3] += x*y*ms

            u[0] = x + transfer_map.dx + u[1]*z1
            u[2] = y + transfer_map.dy + u[3]*z1

            # experimental
            #v = np.array([transfer_map.dx, transfer_map.dy, transfer_map.length, transfer_map.ms])
            #code = """
            #double x, y, dx,dy,length, ms;
            #dx = V1(0);
            #dy = V1(1);
            #
            #length = V1(2);
            #ms = V1(3);
            #x = U1(0) - dx + U1(1)*length/2.;
            #y = U1(2) - dy + U1(3)*length/2.;
            #U1(1) = U1(1) - ms/2.*(x*x - y*y);
            #U1(3) = U1(3) + x*y*ms;
            #
            #U1(0) = x + dx + U1(1)*length/2.;
            #U1(2) = y + dy + U1(3)*length/2.;
            #"""
            #weave.inline(code, ["u", "v"])

            return u

        def nonl_kick_array_z(u, z, ms):
            #v = np.array([transfer_map.dx, transfer_map.dy, z, ms])
            #L = transfer_map.length
            #dx = transfer_map.dx
            #dy = transfer_map.dy
            #ms = transfer_map.ms # ms = K2*L
            #code = """
            #double x, y;
            #double dx = V1(0);
            #double dy = V1(1);
            #double L = V1(2);
            #double ms = V1(3);
            #int i;
            #double L1 = L/2.;
            #for(i = 0;i<Nu[0]/6;i++ ){
            #    // symplectic map for sextupol
            #    // The Stoermer-Verlet schemes
            #    x = U1(0+i*6) + U1(1+i*6)*L1 - dx;
            #    y = U1(2+i*6) + U1(3+i*6)*L1 - dy;
            #    U1(1+i*6) -= ms/2.*(x*x - y*y);
            #    U1(3+i*6) += ms*x*y;
            #    U1(0+i*6) = x + U1(1+i*6)*L1 + dx;
            #    U1(2+i*6) = y + U1(3+i*6)*L1 + dy;
            #}
            #"""
            #weave.inline(code, ["u","v"])

            z1 = z/2.

            x = u[0::6] + u[1::6]*z1 - transfer_map.dx
            y = u[2::6] + u[3::6]*z1 - transfer_map.dy

            u[1::6] += -ms/2.*(x*x - y*y)
            u[3::6] += x*y*ms

            u[0::6] = x + u[1::6]*z1 + transfer_map.dx
            u[2::6] = y + u[3::6]*z1 + transfer_map.dy
            #for i in xrange(len(u)/6):
            #    nonl_kick_z(u[i*6:(i+1)*6], z, ms)
            return u

        if element.ms == None:
            element.ms = element.k2*element.l
        transfer_map.order = 2

        #transfer_map.nonl_kick = nonl_kick #lambda V: nonl_kick(transfer_map, V)
        #transfer_map.nonl_kick_array = nonl_kick_array
        if element.l == 0:
            transfer_map.nonl_kick_z = lambda u, z: nonl_kick_z(u, z, element.ms)
            transfer_map.nonl_kick = lambda u: nonl_kick_z(u, element.l, element.ms)
            transfer_map.nonl_kick_array_z = lambda u, z: nonl_kick_array_z(u, z, element.ms)
            transfer_map.nonl_kick_array = lambda u: nonl_kick_array_z(u, element.l, element.ms)
            transfer_map.ms = element.ms
        else:

            transfer_map.nonl_kick_z = lambda u, z: nonl_kick_z(u, z, element.k2*z)
            transfer_map.nonl_kick = lambda u: nonl_kick_z(u, element.l, element.ms)
            transfer_map.nonl_kick_array_z = lambda u, z: nonl_kick_array_z(u, z, element.k2*z)
            transfer_map.nonl_kick_array = lambda u: nonl_kick_array_z(u, element.l, element.ms)
            #transfer_map.ms = element.ms
        #print "transfer_map.ms = ",transfer_map.ms
        transfer_map.R_z = lambda z: uni_matrix(z, 0, hx = 0)
        transfer_map.R = transfer_map.R_z(element.l)

    elif element.type == "octupole":

        def nonl_kick_z( u, z, moct):

            x = u[0] - transfer_map.dx + u[1]*z/2.
            y = u[2] - transfer_map.dy + u[3]*z/2.

            u[1] += -moct*(x*x*x - 3.*y*y*x)
            u[3] += moct*(3.*y*x*x - y*y*y)

            u[0] = x + transfer_map.dx + u[1]*z/2.
            u[2] = y + transfer_map.dy + u[3]*z/2.

            return u

        def nonl_kick_array_z(u, z, moct):
            #TODO: check expressions
            v = np.array([transfer_map.dx, transfer_map.dy, z, moct])

            code = """
            double x, y;
            double dx = V1(0);
            double dy = V1(1);
            double L = V1(2);
            double moct = V1(3);
            int i;

            for(i = 0;i<Nu[0]/6;i++ ){
                // symplectic map for sextupole
                // The Stoermer-Verlet schemes
                x = U1(0+i*6) + U1(1+i*6)*L/2. - dx;
                y = U1(2+i*6) + U1(3+i*6)*L/2. - dy;
                U1(1+i*6) -= moct/2.*(x*x*x - 3.*y*y*x);
                U1(3+i*6) += moct*(3.*y*x*x - y*y*y);
                U1(0+i*6) = x + U1(1+i*6)*L/2. + dx;
                U1(2+i*6) = y + U1(3+i*6)*L/2. + dy;

            }
            """
            weave.inline(code, ["u","v"])
            return u

        if element.moct == None:
            element.moct = element.k3*element.l
        transfer_map.order = 3

        if element.l == 0:
            transfer_map.nonl_kick_z = lambda u, z: nonl_kick_z(u, z, element.moct)
            transfer_map.nonl_kick = lambda u: nonl_kick_z(u, element.l, element.moct)
            transfer_map.nonl_kick_array_z = lambda u, z: nonl_kick_array_z(u, z, element.moct)
            transfer_map.nonl_kick_array = lambda u: nonl_kick_array_z(u, element.l, element.moct)
            transfer_map.moct = element.moct
        else:

            transfer_map.nonl_kick_z = lambda u, z: nonl_kick_z(u, z, element.k3*z)
            transfer_map.nonl_kick = lambda u: nonl_kick_z(u, element.l, element.moct)
            transfer_map.nonl_kick_array_z = lambda u, z: nonl_kick_array_z(u, z, element.k3*z)
            transfer_map.nonl_kick_array = lambda u: nonl_kick_array_z(u, element.l, element.moct)


        transfer_map.R_z = lambda z: uni_matrix(z, 0, hx = 0)
        transfer_map.R = transfer_map.R_z(element.l)

    elif element.type == "drift":

        transfer_map.R_z = lambda z: uni_matrix(z, 0., hx = 0.)
        transfer_map.R = transfer_map.R_z(element.l)

    elif element.type == "undulator":

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

        def nonl_kick_z(u, z, ndiv):
            ky2 = kz*kz + kx*kx
            ky = sqrt(ky2)
            if element.solver == "rk":
                track_und_openmp(u, z, ndiv*50, kz, kx ,element.Kx, energy)
            else:
                zi = np.linspace(0., z, ndiv)
                h = (zi[1] - zi[0])/(1+u[5])
                gamma = energy*1957.
                h0 = 1./(gamma/element.Kx/kz)
                h02 = h0*h0
                kx2 = kx*kx
                kz2 = kz*kz
                #x, px,y, py, tau = track_und_sym(u, zi, kz, kx, element.Kx, energy)
                for z in range(len(zi)-1):
                    chx = cosh(kx*u[0])
                    chy = cosh(ky*u[2])
                    shx = sinh(kx*u[0])
                    shy = sinh(ky*u[2])
                    u[1] -= h/2.*chx*shx*(kx*ky2*chy*chy + kx2*kx*shy*shy)/(ky2*kz2)*h02
                    u[3] -= h/2.*chy*shy*(ky2*chx*chx + kx2*shx*shx)/(ky*kz2)*h02
                    u[4] -= h/2./(1.+u[5]) * ((u[1]*u[1] + u[3]*u[3]) + chx*chx*chy*chy/(2.*kz2)*h02 + shx*shx*shy*shy*kx2/(2.*ky2*kz2)*h02)
                    u[0] += h*u[1]
                    u[2] += h*u[3]
            return u

        def nonl_kick_array_z(u, z, ndiv):
            if element.solver == "rk":

                track_und_openmp(u, z, 5000, kz, kx ,element.Kx, energy)
            else:
                #for i in xrange(len(u)/6):
                #    V = u[i*6:i*6+6]
                #    #u[i*6:i*6+6] = nonl_kick_z(V, z, ndiv)
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
            transfer_map.nonl_kick = lambda u: nonl_kick_z(u, element.l, ndiv = 20)
            #transfer_map.nonl_kick = lambda u: nonl_kick_rk(u)
            transfer_map.nonl_kick_array = lambda u: nonl_kick_array_z(u, element.l, ndiv = 20)
            transfer_map.nonl_kick_array_z = lambda u, z: nonl_kick_array_z(u, z, ndiv = 5)
            transfer_map.nonl_kick_z = lambda u, z: nonl_kick_z(u, z, ndiv = 5)

    elif element.type == "monitor":
        
        transfer_map.R_z = lambda z: uni_matrix(z, 0., hx = 0.)
        transfer_map.R = transfer_map.R_z(element.l)

    elif element.type == "marker":

        transfer_map.R_z = lambda z: eye(6)
        transfer_map.R = transfer_map.R_z(element.l)

    elif element.type == "hcor" or element.type == "vcor":

        def kick_b(z,l,angle_x, angle_y):
            if l == 0:
                hx = 0.; hy = 0.
            else:
                hx = angle_x/l; hy = angle_y/l
            ux = 0
            uy = 0
            if hx == 0:
                ux = 1
            if hy == 0:
                uy = 1
            b = array([(1 - cos(z*hx))/(hx+ux), sin(z*hx)+ux*angle_x, (1 - cos(z*hy))/(hy+uy), sin(z*hy)+uy*angle_y, 0, 0])
            return b

        transfer_map.order = 0
        bend(element,transfer_map)

        if element.type == "hcor":
            transfer_map.b = lambda z: kick_b(z,element.l,element.angle, 0)
        else:
            transfer_map.b = lambda z: kick_b(z,element.l,0, element.angle)

    elif element.type == "cavity":
                
        def cavity_R_z(z, de, f, E):
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
            Ep = de * cos(phi) / (z) # energy derivative
            Ef = E + de 
            
            alpha = sqrt(eta / 8.) / cos(phi) * log(Ef/E)
            
            #print 'cavity map Ei=',E, 'Ef=', Ef, 'alpha=', alpha
            
            r11 = cos(alpha) - sqrt(2./eta) * cos(phi) * sin(alpha)
            r12 = sqrt(8./eta) * Ef / Ep * cos(phi) * sin(alpha)
            r21 = -Ep/E * (cos(phi)/ sqrt(2*eta) + sqrt(eta/8.) / cos(phi) ) * sin(alpha)
            r22 = E/Ef * ( cos(alpha) + sqrt(2./eta) * cos(phi) * sin(alpha) )
            
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
            transfer_map.R_z = lambda z: cavity_R_z(z, de = element.delta_e * z / element.l, f=element.f, E=element.E)
            transfer_map.R = transfer_map.R_z(element.l)


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

            sol_matrix = array([[c*c, c*s_k, s*c ,s*s_k ,0., 0.],
                                [-k*s*c, c*c, -k*s*s ,s*c ,0., 0.],
                                [-s*c, -s*s_k, c*c, c*s_k, 0., 0.],
                                [k*s*s, -s*c, -k*s*c, c*c, 0., 0.],
                                [0., 0., 0., 0. ,1. ,0.],
                                [0., 0., 0., 0. ,0. ,1.]]).real
            return sol_matrix

        transfer_map.R_z = lambda z: sol(z, k = element.k)
        transfer_map.R = transfer_map.R_z(element.l)

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
        def matriz(z,l, Rm):
            if z < l:
                R_z = uni_matrix(z, 0, hx = 0)
            else:
                R_z = Rm
            return R_z
        transfer_map.R_z = lambda z: matriz(z, element.l, Rm)

    else:

        print element.type , " : unknown type of magnetic element. Cannot create transfer map "

    return transfer_map



def periodic_solution(tws, transfer_matrix):

    """ find periodical twiss  """


    tws = Twiss(tws)
    M = transfer_matrix

    cosmx = (M[0,0] + M[1,1])/2.
    cosmy = (M[2,2] + M[3,3])/2.

    #print cosmx, cosmy

    if abs(cosmx) >= 1 or abs(cosmy) >= 1:
        print "************ periodic solution does not exist. return None ***********"
        return None
    sinmx = sqrt(1.-cosmx*cosmx)
    sinmy = sqrt(1.-cosmy*cosmy)
    tws.beta_x = abs(M[0,1]/sinmx)
    tws.beta_y = abs(M[2,3]/sinmy)


    tws.alpha_x = (M[0,0] - M[1,1])/(2*sinmx) #X[0,0]
    tws.gamma_x = (1. + tws.alpha_x*tws.alpha_x)/tws.beta_x#X[1,0]

    tws.alpha_y = (M[2,2] - M[3,3])/(2*sinmy)#Y[0,0]
    tws.gamma_y = (1. + tws.alpha_y*tws.alpha_y)/tws.beta_y# Y[1,0]

    Hx = array([[M[0,0] - 1, M[0,1] ], [M[1,0], M[1,1]-1]])
    Hhx = array([[M[0,5]], [M[1,5]]])
    hh = dot(inv(-Hx),Hhx)
    tws.Dx = hh[0,0]
    tws.Dxp = hh[1,0]

    Hy = array([[M[2,2] - 1, M[2,3] ], [M[3,2], M[3,3]-1]])
    Hhy = array([[M[2,5]], [M[3,5]]])
    hhy = dot(inv(-Hy),Hhy)
    tws.Dy = hhy[0,0]
    tws.Dyp = hhy[1,0]
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
        #print elem.transfer_map(z - (L - elem.l)).dx
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


def trace_particle(lattice, p0, nPoints = None):
    '''
    track a particle using generic 'trace_obj' function
    '''
    if p0.__class__ == Particle:
        p_list = trace_obj(lattice, p0, nPoints)
    else:
        exit("unknown object p0")
    return p_list


class Navigator:
    #def __init__(self, lattice):
    #    self.lat = lattice
        
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

def get_map(lattice, dz, navi):
    TM = []
    tm = TransferMap()
    i = navi.n_elem
    z1 = navi.z0 + dz
    elem = lattice.sequence[i]
    L = navi.sum_lengths + elem.l

    while z1 > L:
        dl = L - navi.z0
        #print dl, L, navi.z0, elem.l
        if elem.transfer_map.order >1:
            TM.append(tm)
            TM.append(elem.transfer_map(dl))
            tm = TransferMap()
        else:
            tm = elem.transfer_map(dl)*tm
        navi.z0 = L
        dz -= dl
        i += 1
        elem = lattice.sequence[i]
        L += elem.l

    if elem.transfer_map.order > 1:
        TM.append(elem.transfer_map(dz))
    else:
        tm = elem.transfer_map(dz)*tm
    navi.z0 += dz
    navi.sum_lengths = L - elem.l
    navi.n_elem = i
    TM.append(tm)
    return TM


def track(lat, particle_list, dz, navi):
    '''
    tracking for a fixed step dz
    '''
    #print navi.z0 + dz , lat.totalLen
    if navi.z0 + dz > lat.totalLen:
        dz = lat.totalLen - navi.z0

    t_maps = get_map(lat, dz, navi)

    for tm in t_maps:
        #print "tm :",  tm.length
        tm.apply(particle_list)
        #tm.apply(plist)

    return

def gaussFromTwiss(emit, beta, alpha):
    phi = 2*pi * random.rand()
    u = random.rand()
    a = sqrt(-log( (1-u)) * emit)
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
        print '********* calculating fodo parameters *********'
        print 'fodo parameters:'
        print 'k*l=', k*lquad
        print 'f=', L * kap1, L*kap2
        print 'kap1=', kap1
        print 'kap2=', kap2
        print 'betaMax=', betaMax
        print 'betaMin=', betaMin
        print 'betaMean=', betaMean
        print '*********                             *********'
    
    return k*lquad, betaMin, betaMax, betaMean