'''
definition of particles, beams and trajectories
'''
import numpy as np
from numpy.core.umath import sqrt, cos, sin
from ocelot.common.globals import *
import pickle
'''
Note:
(A) the reference frame (e.g. co-moving or not with the beam is not fixed) 
(B) xp, yp are in [rad] but the meaning is not specified
'''

from numpy import *

class Twiss:
    """
    class - container for twiss parameters
    """
    def __init__(self, beam = None):
        if beam == None:
            self.emit_x = 0.0 # ???
            self.emit_y = 0.0 # ???
            self.beta_x = 0.0
            self.beta_y = 0.0
            self.alpha_x = 0.0
            self.alpha_y = 0.0
            self.gamma_x = 0.0
            self.gamma_y = 0.0
            self.mux = 0.0
            self.muy = 0.0
            #self.dQ = 0.
            self.Dx = 0.0
            self.Dy = 0.0
            self.Dxp = 0.0
            self.Dyp = 0.0
            self.x = 0.0
            self.y = 0.0
            self.xp = 0.0
            self.yp = 0.0
            self.E = 0.0
            self.p = 0.0
            self.tau = 0.0
            self.s = 0.0 # position along the reference trajectory
            self.id = ""
        else:
            self.emit_x = beam.emit_x
            self.emit_y = beam.emit_y
            self.beta_x = beam.beta_x
            self.beta_y = beam.beta_y
            self.alpha_x = beam.alpha_x
            self.alpha_y = beam.alpha_y
            self.mux = 0.
            self.muy = 0.
            #self.dQ = 0.
            self.Dx = beam.Dx
            self.Dy = beam.Dy
            self.Dxp = beam.Dxp
            self.Dyp = beam.Dyp
            self.x = beam.x
            self.y = beam.y
            self.xp = beam.xp
            self.yp = beam.yp
            if beam.beta_x == 0.0 or beam.beta_y == 0.0:
                self.gamma_x = 0.0
                self.gamma_y = 0.0
            else:
                self.gamma_x = (1 + beam.alpha_x * beam.alpha_x) / beam.beta_x
                self.gamma_y = (1 + beam.alpha_y * beam.alpha_y) / beam.beta_y
            self.E = beam.E
            self.p = 0.0
            self.tau = 0.0
            self.s = 0.0 # position along the reference trajectory
            self.id = ""

    def __str__(self):
        val = ""
        val += "emit_x  = " + str(self.emit_x) + "\n"
        val += "emit_y  = " + str(self.emit_y) + "\n"
        val += "beta_x  = " + str(self.beta_x) + "\n"
        val += "beta_y  = " + str(self.beta_y) + "\n"
        val += "alpha_x = " + str(self.alpha_x) + "\n"
        val += "alpha_y = " + str(self.alpha_y) + "\n"
        val += "gamma_x = " + str(self.gamma_x) + "\n"
        val += "gamma_y = " + str(self.gamma_y) + "\n"
        val += "Dx      = " + str(self.Dx) + "\n"
        val += "Dy      = " + str(self.Dy) + "\n"
        val += "Dxp     = " + str(self.Dxp) + "\n"
        val += "Dyp     = " + str(self.Dyp) + "\n"
        val += "mux     = " + str(self.mux) + "\n"
        val += "muy     = " + str(self.muy) + "\n"
        val += "nu_x    = " + str(self.mux/2./pi) + "\n"
        val += "nu_y    = " + str(self.muy/2./pi) + "\n"
        val += "E       = " + str(self.E) + "\n"
        val += "s        = " + str(self.s) + "\n"
        return val

            
class Particle:
    '''
    particle
    to be used for tracking
    '''
    def __init__(self, x=0.0, y=0.0, px=0.0, py=0.0, s=0.0, p=0.0,  tau=0.0, E=0.0):
        self.x = x
        self.y = y
        self.px = px       # horizontal (generalized) momentum
        self.py = py       # vertical (generalized) momentum 
        self.p = p         # longitudinal momentum
        self.s = s
        self.tau = tau     # time-like coordinate wrt reference particle in the bunch (e.g phase)
        self.E = E        # energy

    def __str__(self):
        val = ""
        val = val + "x = " + str(self.x) + "\n"
        val = val + "px = " + str(self.px) + "\n"
        val = val + "y = " + str(self.y) + "\n"
        val = val + "py = " + str(self.py) + "\n"
        val = val + "tau = " + str(self.tau) + "\n"
        val = val + "p = " + str(self.p) + "\n"
        val = val + "E = " + str(self.E) + "\n"
        val = val + "s = " + str(self.s)
        return val

class Beam:
    def __init__(self,x=0,xp=0,y=0,yp=0):
        # initial conditions
        self.x = x      #[m]
        self.y = y      #[m]
        self.xp = xp    # xprime [rad]
        self.yp = yp    # yprime [rad]

        self.E = 0.0            # electron energy [GeV]
        self.sigma_E = 0.0      # Energy spread [GeV]
        self.I = 0.0            # beam current [A]
        self.emit_x = 0.0       # horizontal emittance [m rad]
        self.emit_y = 0.0       # horizontal emittance [m rad]
        self.tlen = 0.0         # bunch length (rms) in fsec

        # twiss parameters
        self.beta_x = 0.0
        self.beta_y = 0.0
        self.alpha_x = 0.0
        self.alpha_y = 0.0
        self.Dx = 0.0
        self.Dy = 0.0
        self.Dxp = 0.0
        self.Dyp = 0.0

    def sizes(self):
        if self.beta_x != 0:
            self.gamma_x = (1. + self.alpha_x**2)/self.beta_x
        else:
            self.gamma_x = 0.

        if self.beta_y != 0:
            self.gamma_y = (1. + self.alpha_y**2)/self.beta_y
        else:
            self.gamma_y = 0.

        self.sigma_x = sqrt((self.sigma_E/self.E*self.Dx)**2 + self.emit_x*self.beta_x)
        self.sigma_y = sqrt((self.sigma_E/self.E*self.Dy)**2 + self.emit_y*self.beta_y)
        self.sigma_xp = sqrt((self.sigma_E/self.E*self.Dxp)**2 + self.emit_x*self.gamma_x)
        self.sigma_yp = sqrt((self.sigma_E/self.E*self.Dyp)**2 + self.emit_y*self.gamma_y)

    def print_sizes(self):
        self.sizes()
        print("sigma_E/E and Dx/Dy : ", self.sigma_E/self.E, "  and ", self.Dx, "/",self.Dy, " m")
        print("emit_x/emit_y     : ",  self.emit_x*1e9, "/",self.emit_y*1e9, " nm-rad")
        print("sigma_x/y         : ", self.sigma_x*1e6, "/", self.sigma_y*1e6, " um")
        print("sigma_xp/yp       : ", self.sigma_xp*1e6, "/", self.sigma_yp*1e6, " urad")


class Trajectory:
    def __init__(self):
        self.ct = []
        self.E = []
        self.x = []
        self.y = []
        self.xp = []
        self.yp = []
        self.z = []
        self.s = []

    def add(self, ct, x, y, xp, yp, z, s):
        self.ct.append(ct)
        self.x.append(x)
        self.y.append(y)
        self.xp.append(xp)
        self.yp.append(yp)
        self.z.append(z)
        self.s.append(s)

    def last(self):
        p = Particle()
        
        p.ct = self.ct[len(self.ct)-1]
        p.x = self.x[len(self.x)-1]
        p.y = self.y[len(self.y)-1]
        p.xp = self.xp[len(self.xp)-1]
        p.yp = self.yp[len(self.yp)-1]
        try:
            p.E = self.E[len(self.E)-1]
        except IndexError:
            return 0
        p.s = self.s[len(self.s)-1]

        return p


class ParticleArray:
    """
    array of particles of fixed size; for optimized performance
    """
    def __init__(self, n=0):
        self.particles = zeros(n*6)
        self.q_array = np.zeros(n)    # charge
        self.s = 0.0
        self.E = 0.0

    def rm_tails(self, xlim, ylim, px_lim, py_lim):
        """
        comment behaviour and possibly move out of class
        """
        x = abs(self.particles[::6])
        px = abs(self.particles[1::6])
        y = abs(self.particles[2::6])
        py = abs(self.particles[3::6])
        ind_angles = append(argwhere(px > px_lim), argwhere(py > py_lim))
        p_idxs = unique(append(argwhere(x > xlim), append(argwhere(y > ylim), append(argwhere(x != x), append(argwhere(y!= y), ind_angles)) )))
        e_idxs = [append([], x) for x in array([6*p_idxs, 6*p_idxs+1, 6*p_idxs+2, 6*p_idxs+3, 6*p_idxs+4, 6*p_idxs+5])]
        self.particles = delete(self.particles, e_idxs)
        return p_idxs

    def __getitem__(self, idx):
        return Particle(x=self.particles[idx*6], px=self.particles[idx*6 + 1],
                         y=self.particles[idx*6+2], py=self.particles[idx*6+3],
                         tau=self.particles[idx*6+4], p=self.particles[idx*6+5],
                         s=self.s)

    def __setitem__(self, idx, p):
        self.particles[idx*6] = p.x
        self.particles[idx*6 + 1] = p.px
        self.particles[idx*6+2] = p.y
        self.particles[idx*6+3] = p.py
        self.particles[idx*6+4] = p.tau
        self.particles[idx*6+5] = p.p
        self.s = p.s

    def list2array(self, p_list):
        self.particles = zeros(len(p_list)*6)
        for i, p in enumerate(p_list):
            self[i] = p
        self.s = p_list[0].s
        self.E = p_list[0].E

    def array2list(self):
        p_list = []
        for i in range(self.size()):
            p_list.append( self[i])
        return p_list

    def array2ex_list(self, p_list):

        for i, p in enumerate(p_list):
            p.x = self.particles[i*6]
            p.px = self.particles[i*6 + 1]
            p.y = self.particles[i*6+2]
            p.py = self.particles[i*6+3]
            p.tau = self.particles[i*6+4]
            p.p = self.particles[i*6+5]
            p.E = self.E
            p.s = self.s
        return p_list

    def size(self):
        return len(self.particles) / 6

    def x(self): return self.particles[::6]
    def px(self): return self.particles[1::6]
    def y(self): return self.particles[2::6]
    def py(self): return self.particles[3::6]
    def tau(self): return self.particles[4::6]
    def p(self): return self.particles[5::6]


def save_particle_array(p_array, filename):
    output = open(filename, 'wb')
    pickle.dump(p_array, output)
    output.close()


def load_particle_array(filename):
    output = open(filename, 'rb')
    p_array = pickle.load(output)
    output.close()
    return p_array


def get_envelope(p_array, tws_i=Twiss()):
    tws = Twiss()
    p = p_array.p()
    dx = tws_i.Dx*p
    dy = tws_i.Dy*p
    dpx = tws_i.Dxp*p
    dpy = tws_i.Dyp*p
    x = p_array.x() - dx
    px = p_array.px() - dpx

    y = p_array.y() - dy
    py = p_array.py() - dpy
    px = px*(1.-0.5*px*px - 0.5*py*py)
    py = py*(1.-0.5*px*px - 0.5*py*py)
    tws.x = mean(x)
    tws.y = mean(y)
    tws.px =mean(px)
    tws.py =mean(py)
    #print tws.x, tws.y, tws.px,tws.py
    tws.xx = mean((x - tws.x)*(x - tws.x))
    tws.xpx = mean((x-tws.x)*(px-tws.px))
    tws.pxpx = mean((px-tws.px)*(px-tws.px))
    tws.yy = mean((y-tws.y)*(y-tws.y))
    tws.ypy = mean((y-tws.y)*(py-tws.py))
    tws.pypy = mean((py-tws.py)*(py-tws.py))
    tws.p = mean( p_array.p())
    tws.E = p_array.E
    #tws.de = p_array.de

    tws.emit_x = sqrt(tws.xx*tws.pxpx-tws.xpx**2)
    tws.emit_y = sqrt(tws.yy*tws.pypy-tws.ypy**2)
    #print tws.emit_x, sqrt(tws.xx*tws.pxpx-tws.xpx**2), tws.emit_y, sqrt(tws.yy*tws.pypy-tws.ypy**2)
    tws.beta_x = tws.xx/tws.emit_x
    tws.beta_y = tws.yy/tws.emit_y
    tws.alpha_x = -tws.xpx/tws.emit_x
    tws.alpha_y = -tws.ypy/tws.emit_y
    return tws

def get_current(p_array, charge, num_bins = 200):
    """
    return: hist -  current in A
          : bin_edges - points position
    """
    z = p_array.particles[4::6]
    hist, bin_edges = np.histogram(z, bins=num_bins)
    delta_Z = max(z) - min(z)
    delta_z = delta_Z/num_bins
    t_bins = delta_z/speed_of_light
    print( "Imax = ", max(hist)*charge/t_bins)
    hist = np.append(hist, hist[-1])
    return bin_edges, hist*charge/t_bins


def gauss_from_twiss(emit, beta, alpha):
    phi = 2*pi * np.random.rand()
    u = np.random.rand()
    a = sqrt(-2*np.log( (1-u)) * emit)
    x = a * sqrt(beta) * cos(phi)
    xp = -a / sqrt(beta) * ( sin(phi) + alpha * cos(phi) )
    return (x, xp)

def waterbag_from_twiss(emit, beta, alpha):
    phi = 2*pi * np.random.rand()
    a = sqrt(emit) * np.random.rand()
    x = a * sqrt(beta) * cos(phi)
    xp = -a / sqrt(beta) * ( sin(phi) + alpha * cos(phi) )
    return (x, xp)

def ellipse_from_twiss(emit, beta, alpha):
    phi = 2*pi * np.random.rand()
    #u = np.random.rand()
    #a = sqrt(-2*np.log( (1-u)) * emit)
    a = sqrt(emit)
    x = a * sqrt(beta) * cos(phi)
    xp = -a / sqrt(beta) * ( sin(phi) + alpha * cos(phi) )
    return (x, xp)


def moments(x, y, cut=0):
    n = len(x)
    inds = np.arange(n)
    mx = np.mean(x)
    my = np.mean(y)
    x = x - mx
    y = y - my
    x2 = x*x
    mxx = sum(x2)/n
    y2 = y*y
    myy = sum(y2)/n
    xy = x*y
    mxy = sum(xy)/n

    emitt = sqrt(mxx*myy - mxy*mxy)

    if cut>0:
        inds=[]
        beta = mxx/emitt
        gamma = myy/emitt
        alpha = mxy/emitt
        emittp = gamma*x2 + 2.*alpha*xy + beta*y2
        inds0 = np.argsort(emittp)
        n1 = np.round(n*(100-cut)/100)
        inds = inds0[0:n1]
        mx = np.mean(x[inds])
        my = np.mean(y[inds])
        x1 = x[inds] - mx
        y1 = y[inds] - my
        mxx = np.sum(x1*x1)/n1
        myy = np.sum(y1*y1)/n1
        mxy = np.sum(x1*y1)/n1
        emitt = sqrt(mxx*myy - mxy*mxy)
    return mx, my, mxx, mxy, myy, emitt

def m_from_twiss(Tw1, Tw2):
    #% Transport matrix M for two sets of Twiss parameters (alpha,beta,psi)
    b1 = Tw1[1]
    a1 = Tw1[0]
    psi1 = Tw1[2]
    b2 = Tw2[1]
    a2 = Tw2[0]
    psi2 = Tw2[2]

    psi = psi2-psi1
    cosp = cos(psi)
    sinp = sin(psi)
    M = np.zeros((2, 2))
    M[0, 0] = sqrt(b2/b1)*(cosp+a1*sinp)
    M[0, 1] = sqrt(b2*b1)*sinp
    M[1, 0] = ((a1-a2)*cosp-(1+a1*a2)*sinp)/sqrt(b2*b1)
    M[1, 1] = sqrt(b1/b2)*(cosp-a2*sinp)
    return M

def beam_matching(particles, bounds, x_opt, y_opt):
    pd = zeros(( int(len(particles)/6), 6))
    pd[:, 0] = particles[0::6]
    pd[:, 1] = particles[1::6]
    pd[:, 2] = particles[2::6]
    pd[:, 3] = particles[3::6]
    pd[:, 4] = particles[4::6]
    pd[:, 5] = particles[5::6]

    z0 = np.mean(pd[:, 4])
    sig0 = np.std(pd[:, 4])
    #print((z0 + sig0*bounds[0] <= pd[:, 4]) * (pd[:, 4] <= z0 + sig0*bounds[1]))
    inds = np.argwhere((z0 + sig0*bounds[0] <= pd[:, 4]) * (pd[:, 4] <= z0 + sig0*bounds[1]))
    #print(moments(pd[inds, 0], pd[inds, 1]))
    mx, mxs, mxx, mxxs, mxsxs, emitx0 = moments(pd[inds, 0], pd[inds, 1])
    beta = mxx/emitx0
    alpha = -mxxs/emitx0
    print(beta, alpha)
    M = m_from_twiss([alpha, beta, 0], x_opt)
    print(M)
    particles[0::6] = M[0, 0]*pd[:, 0] + M[0, 1]*pd[:, 1]
    particles[1::6] = M[1, 0]*pd[:, 0] + M[1, 1]*pd[:, 1]
    [mx, mxs, mxx, mxxs, mxsxs, emitx0] = moments(pd[inds, 2], pd[inds, 3])
    beta = mxx/emitx0
    alpha = -mxxs/emitx0
    M = m_from_twiss([alpha, beta, 0], y_opt)
    particles[2::6] = M[0, 0]*pd[:, 2] + M[0, 1]*pd[:, 3]
    particles[3::6] = M[1, 0]*pd[:, 2] + M[1, 1]*pd[:, 3]
    return particles


class BeamTransform:
    def __init__(self, x_opt, y_opt):
        """
        Beam matching

        :param x_opt: [alpha, beta, mu (phase advance)]
        :param y_opt: [alpha, beta, mu (phase advance)]
        """
        self.bounds = [-5, 5]  # [start, stop] in sigmas
        self.x_opt = x_opt   # [alpha, beta, mu (phase advance)]
        self.y_opt = y_opt   # [alpha, beta, mu (phase advance)]
        self.step=1

    def prepare(self, lat):
        pass

    def apply(self, p_array, dz):
        print()
        print("BEAM TRANSFORM")
        self.beam_matching(p_array.particles, self.bounds, self.x_opt, self.y_opt)

    def beam_matching(self, particles, bounds, x_opt, y_opt):
        pd = zeros((int(len(particles) / 6), 6))
        pd[:, 0] = particles[0::6]
        pd[:, 1] = particles[1::6]
        pd[:, 2] = particles[2::6]
        pd[:, 3] = particles[3::6]
        pd[:, 4] = particles[4::6]
        pd[:, 5] = particles[5::6]

        z0 = np.mean(pd[:, 4])
        sig0 = np.std(pd[:, 4])
        # print((z0 + sig0*bounds[0] <= pd[:, 4]) * (pd[:, 4] <= z0 + sig0*bounds[1]))
        inds = np.argwhere((z0 + sig0 * bounds[0] <= pd[:, 4]) * (pd[:, 4] <= z0 + sig0 * bounds[1]))
        # print(moments(pd[inds, 0], pd[inds, 1]))
        mx, mxs, mxx, mxxs, mxsxs, emitx0 = self.moments(pd[inds, 0], pd[inds, 1])
        beta = mxx / emitx0
        alpha = -mxxs / emitx0
        #print(beta, alpha)
        M = m_from_twiss([alpha, beta, 0], x_opt)
        #print(M)
        particles[0::6] = M[0, 0] * pd[:, 0] + M[0, 1] * pd[:, 1]
        particles[1::6] = M[1, 0] * pd[:, 0] + M[1, 1] * pd[:, 1]
        [mx, mxs, mxx, mxxs, mxsxs, emitx0] = self.moments(pd[inds, 2], pd[inds, 3])
        beta = mxx / emitx0
        alpha = -mxxs / emitx0
        M = m_from_twiss([alpha, beta, 0], y_opt)
        particles[2::6] = M[0, 0] * pd[:, 2] + M[0, 1] * pd[:, 3]
        particles[3::6] = M[1, 0] * pd[:, 2] + M[1, 1] * pd[:, 3]
        return particles

    def moments(self, x, y, cut=0):
        n = len(x)
        inds = np.arange(n)
        mx = np.mean(x)
        my = np.mean(y)
        x = x - mx
        y = y - my
        x2 = x * x
        mxx = sum(x2) / n
        y2 = y * y
        myy = sum(y2) / n
        xy = x * y
        mxy = sum(xy) / n

        emitt = sqrt(mxx * myy - mxy * mxy)

        if cut > 0:
            inds = []
            beta = mxx / emitt
            gamma = myy / emitt
            alpha = mxy / emitt
            emittp = gamma * x2 + 2. * alpha * xy + beta * y2
            inds0 = np.argsort(emittp)
            n1 = np.round(n * (100 - cut) / 100)
            inds = inds0[0:n1]
            mx = np.mean(x[inds])
            my = np.mean(y[inds])
            x1 = x[inds] - mx
            y1 = y[inds] - my
            mxx = np.sum(x1 * x1) / n1
            myy = np.sum(y1 * y1) / n1
            mxy = np.sum(x1 * y1) / n1
            emitt = sqrt(mxx * myy - mxy * mxy)
        return mx, my, mxx, mxy, myy, emitt

class EmptyProc():
    def __init__(self):
        self.step=1
        self.energy = None
        self.pict_debug = True
        self.traj_step = 0.0002
        #if self.pict_debug:
        #    self.f = plt.figure(figsize=(12, 9))
        #    plt.ion()

    def prepare(self, lat):
        pass

    def apply(self, p_array, dz):
        pass