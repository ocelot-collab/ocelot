'''
definition of particles, beams and trajectories
'''
import numpy as np
from numpy.core.umath import sqrt, cos, sin
from ocelot.common.globals import *
import pickle
from scipy import interpolate
try:
    import numexpr as ne
    ne_flag = True
except:
    print("beam.py: module NUMEXPR is not installed. Install it if you want higher speed calculation.")
    ne_flag = False


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
    (x, x' = px/p0),(y, y' = py/p0),(ds = c*tau, p = dE/(p0*c))
    p0 - momentum
    """
    def __init__(self, n=0):
        #self.particles = zeros(n*6)
        self.rparticles = zeros((6, n))#np.transpose(np.zeros(int(n), 6))
        self.q_array = np.zeros(n)    # charge
        self.s = 0.0
        self.E = 0.0

    def rm_tails(self, xlim, ylim, px_lim, py_lim):
        """
        comment behaviour and possibly move out of class
        """
        x = abs(self.x())
        px = abs(self.px())
        y = abs(self.y())
        py = abs(self.py())
        ind_angles = append(argwhere(px > px_lim), argwhere(py > py_lim))
        p_idxs = unique(append(argwhere(x > xlim), append(argwhere(y > ylim), append(argwhere(x != x), append(argwhere(y!= y), ind_angles)) )))
        #e_idxs = [append([], x) for x in array([6*p_idxs, 6*p_idxs+1, 6*p_idxs+2, 6*p_idxs+3, 6*p_idxs+4, 6*p_idxs+5])]
        self.rparticles = delete(self.rparticles, p_idxs, axis=1)
        return p_idxs

    def __getitem__(self, idx):
        return Particle(x=self.rparticles[0, idx], px=self.rparticles[1, idx],
                         y=self.rparticles[2, idx], py=self.rparticles[3, idx],
                         tau=self.rparticles[4, idx], p=self.rparticles[5, idx],
                         s=self.s)

    def __setitem__(self, idx, p):
        self.rparticles[0, idx] = p.x
        self.rparticles[1, idx] = p.px
        self.rparticles[2, idx] = p.y
        self.rparticles[3, idx] = p.py
        self.rparticles[4, idx] = p.tau
        self.rparticles[5, idx] = p.p
        self.s = p.s

    def list2array(self, p_list):
        self.rparticles = zeros((6, len(p_list)))
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
            p.x =  self.rparticles[0, i]
            p.px = self.rparticles[1, i]
            p.y =  self.rparticles[2, i]
            p.py = self.rparticles[3, i]
            p.tau =self.rparticles[4, i]
            p.p =  self.rparticles[5, i]
            p.E = self.E
            p.s = self.s
        return p_list

    def size(self):
        return self.rparticles.size / 6

    def x(self):  return self.rparticles[0]
    def px(self): return self.rparticles[1]
    def y(self):  return self.rparticles[2]
    def py(self): return self.rparticles[3]
    def tau(self):return self.rparticles[4]
    def p(self):  return self.rparticles[5]


def save_particle_array(filename, p_array):
    np.savez_compressed(filename, rparticles=p_array.rparticles,
                        q_array=p_array.q_array,
                        E=p_array.E, s=p_array.s)

def load_particle_array(filename):
    p_array = ParticleArray()
    with np.load(filename) as data:
        for key in data.keys():
            p_array.__dict__[key] = data[key]
    return p_array


def recalculate_ref_particle(p_array):
    pref = np.sqrt(p_array.E ** 2 / m_e_GeV ** 2 - 1) * m_e_GeV
    Enew = p_array.p()[0]*pref + p_array.E
    s_new = p_array.s - p_array.tau()[0]
    p_array.rparticles[5, :] -= p_array.p()[0]
    p_array.rparticles[4, :] -= p_array.tau()[0]
    p_array.E = Enew
    p_array.s = s_new
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
    if ne_flag:
        px = ne.evaluate('px * (1. - 0.5 * px * px - 0.5 * py * py)')
        py = ne.evaluate('py * (1. - 0.5 * px * px - 0.5 * py * py)')
    else:
        px = px*(1.-0.5*px*px - 0.5*py*py)
        py = py*(1.-0.5*px*px - 0.5*py*py)
    tws.x = mean(x)
    tws.y = mean(y)
    tws.px =mean(px)
    tws.py =mean(py)

    if ne_flag:
        tw_x = tws.x
        tw_y = tws.y
        tw_px = tws.px
        tw_py = tws.py
        tws.xx =  mean(ne.evaluate('(x - tw_x) * (x - tw_x)'))
        tws.xpx = mean(ne.evaluate('(x - tw_x) * (px - tw_px)'))
        tws.pxpx =mean(ne.evaluate('(px - tw_px) * (px - tw_px)'))
        tws.yy =  mean(ne.evaluate('(y - tw_y) * (y - tw_y)'))
        tws.ypy = mean(ne.evaluate('(y - tw_y) * (py - tw_py)'))
        tws.pypy =mean(ne.evaluate('(py - tw_py) * (py - tw_py)'))
    else:
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
    z = p_array.tau()
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
    #inds = np.arange(n)
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
        #inds=[]
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
    pd[:, 0] = particles[0]
    pd[:, 1] = particles[1]
    pd[:, 2] = particles[2]
    pd[:, 3] = particles[3]
    pd[:, 4] = particles[4]
    pd[:, 5] = particles[5]

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
    particles[0] = M[0, 0]*pd[:, 0] + M[0, 1]*pd[:, 1]
    particles[1] = M[1, 0]*pd[:, 0] + M[1, 1]*pd[:, 1]
    [mx, mxs, mxx, mxxs, mxsxs, emitx0] = moments(pd[inds, 2], pd[inds, 3])
    beta = mxx/emitx0
    alpha = -mxxs/emitx0
    M = m_from_twiss([alpha, beta, 0], y_opt)
    particles[2] = M[0, 0]*pd[:, 2] + M[0, 1]*pd[:, 3]
    particles[3] = M[1, 0]*pd[:, 2] + M[1, 1]*pd[:, 3]
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
        self.beam_matching(p_array.rparticles, self.bounds, self.x_opt, self.y_opt)

    def beam_matching(self, particles, bounds, x_opt, y_opt):
        pd = zeros((int(particles.size / 6), 6))
        pd[:, 0] = particles[0]
        pd[:, 1] = particles[1]
        pd[:, 2] = particles[2]
        pd[:, 3] = particles[3]
        pd[:, 4] = particles[4]
        pd[:, 5] = particles[5]

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
        particles[0] = M[0, 0] * pd[:, 0] + M[0, 1] * pd[:, 1]
        particles[1] = M[1, 0] * pd[:, 0] + M[1, 1] * pd[:, 1]
        [mx, mxs, mxx, mxxs, mxsxs, emitx0] = self.moments(pd[inds, 2], pd[inds, 3])
        beta = mxx / emitx0
        alpha = -mxxs / emitx0
        M = m_from_twiss([alpha, beta, 0], y_opt)
        particles[2] = M[0, 0] * pd[:, 2] + M[0, 1] * pd[:, 3]
        particles[3] = M[1, 0] * pd[:, 2] + M[1, 1] * pd[:, 3]
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






def sortrows(x, col):
    return x[:, x[col].argsort()]


def convmode(A, B, mode):
    C=[]
    if mode == 2:
        C = np.convolve(A, B)
    if mode == 1:
        i = np.int_(np.floor(len(B)*0.5))
        n = len(A)
        C1 = np.convolve(A, B)
        C[:n] = C1[i:n+i]
    return C

def s_to_cur(A, sigma, q0, v):
    #  A - s-coordinates of particles
    #  sigma -smoothing parameter
    #  q0 -bunch charge
    #  v mean velocity
    Nsigma = 3
    a = np.min(A) - Nsigma*sigma
    b = np.max(A) + Nsigma*sigma
    s = 0.25*sigma
    #print("s = ", sigma, s)
    N = int(np.ceil((b-a)/s))
    s = (b-a)/N
    #print(a, b, N)
    B = np.zeros((N+1, 2))
    C = np.zeros(N+1)

    #print(len(np.arange(0, (N+0.5)*s, s)))
    B[:, 0] = np.arange(0, (N+0.5)*s, s) + a
    N = np.shape(B)[0]
    #print(N)
    cA = (A - a)/s
    #print(cA[:10])
    I = np.int_(np.floor(cA))
    #print(I)
    xiA = 1 + I - cA
    for k in range(len(A)):
        i = I[k]
        if i > N-1:
            i = N-1
        C[i+0] = C[i+0]+xiA[k]
        C[i+1] = C[i+1]+(1-xiA[k])

    K = np.floor(Nsigma*sigma/s + 0.5)
    G = np.exp(-0.5*(np.arange(-K, K+1)*s/sigma)**2)
    G = G/np.sum(G)
    B[:, 1] = convmode(C, G, 1)
    koef = q0*v/(s*np.sum(B[:, 1]))
    B[:, 1] = koef*B[:, 1]
    return B


def slice_analysis(z, x, xs, M, to_sort):
    z = np.copy(z)
    if to_sort:
        #P=sortrows([z, x, xs])
        indx = z.argsort()
        z = z[indx]
        x = x[indx]
        xs = xs[indx]
        P=[]
    N=len(x)
    mx = np.zeros(N)
    mxs= np.zeros(N)
    mxx = np.zeros(N)
    mxxs = np.zeros(N)
    mxsxs = np.zeros(N)
    emittx = np.zeros(N)
    m = np.max([np.round(M/2), 1])
    xc = np.cumsum(x)
    xsc = np.cumsum(xs)
    for i in range(N):
        n1 = int(max(0, i-m))
        n2 = int(min(N-1, i+m))
        dq = n2 - n1
        mx[i] = (xc[n2] - xc[n1])/dq
        mxs[i] = (xsc[n2] - xsc[n1])/dq

    x = x - mx
    xs = xs - mxs
    x2c = np.cumsum(x*x)
    xs2c = np.cumsum(xs*xs)
    xxsc = np.cumsum(x*xs)
    for i in range(N):
        n1 = int(max(0, i-m))
        n2 = int(min(N-1, i+m))
        dq = n2 - n1
        mxx[i] = (x2c[n2] - x2c[n1])/dq
        mxsxs[i] = (xs2c[n2] - xs2c[n1])/dq
        mxxs[i] = (xxsc[n2] - xxsc[n1])/dq

    emittx = np.sqrt(mxx*mxsxs - mxxs*mxxs)
    return [mx, mxs, mxx, mxxs, mxsxs, emittx]


def simple_filter(x, p, iter):
    n = len(x)
    if iter == 0:
        y = x
        return y
    for k in range(iter):

        y = np.zeros(n)
        for i in range(n):
            i0 = i - p
            if i0 < 0:
                i0 = 0
            i1 = i + p
            if i1 > n-1:
                i1 = n-1
            s = 0
            for j in range(i0, i1+1):
                s = s + x[j]
            y[i] = s / (i1 - i0 + 1)

        x = y
    return y

def interp1(x, y, xnew, k=1):
    if len(xnew) > 0:
        tck = interpolate.splrep(x, y, k=k)
        ynew = interpolate.splev(xnew, tck, der=0)
    else:
        ynew = []
    return ynew

def slice_analysis_transverse(parray, Mslice, Mcur, p, iter):
    q1 = np.sum(parray.q_array)
    print("charge", q1)
    n = np.int_(parray.rparticles.size / 6)
    PD = parray.rparticles
    PD = sortrows(PD, col=4)

    z = np.copy(PD[4])
    mx, mxs, mxx, mxxs, mxsxs, emittx = slice_analysis(z, PD[0], PD[1], Mslice, True)

    my, mys, myy, myys, mysys, emitty = slice_analysis(z, PD[2], PD[3], Mslice, True)

    mm, mm, mm, mm, mm, emitty0 = moments(PD[2], PD[3])
    gamma0 = parray.E / m_e_GeV
    emityn = emitty0*gamma0
    mm, mm, mm, mm, mm, emitt0 = moments(PD[0], PD[1])
    emitxn = emitt0*gamma0

    z, ind = np.unique(z, return_index=True)
    emittx = emittx[ind]
    emitty = emitty[ind]
    smin = min(z)
    smax = max(z)
    n = 1000
    hs = (smax-smin)/(n-1)
    s = np.arange(smin, smax + hs, hs)
    ex = interp1(z, emittx, s)
    ey = interp1(z, emitty, s)

    ex = simple_filter(ex, p, iter)*gamma0*1e6
    ey = simple_filter(ey, p, iter)*gamma0*1e6

    sig0 = np.std(parray.tau())
    B = s_to_cur(z, Mcur*sig0, q1, speed_of_light)
    I = interp1(B[:, 0], B[:, 1], s)
    return [s, I, ex, ey, gamma0, emitxn, emityn]


def global_slice_analysis_extended(parray, Mslice, Mcur, p, iter):
    # %[s, I, ex, ey ,me, se, gamma0, emitxn, emityn]=GlobalSliceAnalysis_Extended(PD,q1,Mslice,Mcur,p,iter)

    q1 = np.sum(parray.q_array)
    print("charge", q1)
    n = np.int_(parray.rparticles.size/6)
    PD = parray.rparticles
    PD = sortrows(PD, col=4)

    z = np.copy(PD[4])
    mx, mxs, mxx, mxxs, mxsxs, emittx = slice_analysis(z, PD[0], PD[1], Mslice, True)
    
    my, mys, myy, myys, mysys, emitty = slice_analysis(z, PD[2], PD[3], Mslice, True)

    pc_0 = np.sqrt(parray.E**2 - m_e_GeV**2)
    E1 = PD[5]*pc_0 + parray.E
    pc_1 = np.sqrt(E1**2 - m_e_GeV**2)
    #print(pc_1[:10])
    mE, mEs, mEE, mEEs, mEsEs, emittE = slice_analysis(z, PD[4], pc_1*1e9, Mslice, True)

    #print(mE, mEs, mEE, mEEs, mEsEs, emittE)
    mE = mEs
    sE = np.sqrt(mEsEs)
    sig0 = np.std(parray.tau())
    B = s_to_cur(z, Mcur*sig0, q1, speed_of_light)
    gamma0 = parray.E/m_e_GeV
    mm, mm, mm, mm, mm, emitty0 = moments(PD[2], PD[3])
    emityn = emitty0*gamma0
    mm, mm, mm, mm, mm, emitt0 = moments(PD[0], PD[1])
    emitxn = emitt0*gamma0

    z, ind = np.unique(z, return_index=True)
    emittx = emittx[ind]
    emitty = emitty[ind]
    sE = sE[ind]
    mE = mE[ind]
    smin = min(z)
    smax = max(z)
    n = 1000
    hs = (smax-smin)/(n-1)
    s = np.arange(smin, smax + hs, hs)
    ex = interp1(z, emittx, s)
    ey = interp1(z, emitty, s)
    se = interp1(z, sE, s)
    me = interp1(z, mE, s)
    ex = simple_filter(ex, p, iter)*gamma0*1e6
    ey = simple_filter(ey, p, iter)*gamma0*1e6
    se = simple_filter(se, p, iter)
    me = simple_filter(me, p, iter)
    I = interp1(B[:, 0], B[:, 1], s)
    return [s, I, ex, ey, me, se, gamma0, emitxn, emityn]