'''
definition of particles, beams and trajectories
'''

'''
Note:
(A) the reference frame (e.g. co-moving or not with the beam is not fixed) 
(B) xp, yp are in [rad] but the meaning is not specified
'''

from numpy import *

class Twiss:
    def __init__(self, beam = None):
        if beam == None:
            self.emit_x = 0 # ???
            self.emit_y = 0 # ???
            self.beta_x = 0
            self.beta_y = 0
            self.alpha_x = 0
            self.alpha_y = 0
            self.gamma_x = 0
            self.gamma_y = 0
            self.mux = 0.
            self.muy = 0.
            #self.dQ = 0.
            self.Dx = 0
            self.Dy = 0
            self.Dxp = 0
            self.Dyp = 0
            self.x = 0
            self.y = 0
            self.xp = 0
            self.yp = 0
            self.E = 0
            self.p = 0
            self.tau = 0
            self.s = 0 # position along the reference trajectory
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
            if beam.beta_x == 0 or beam.beta_y == 0:
                self.gamma_x = 0
                self.gamma_y = 0
            else:
                self.gamma_x = (1 + beam.alpha_x * beam.alpha_x) / beam.beta_x
                self.gamma_y = (1 + beam.alpha_y * beam.alpha_y) / beam.beta_y
            self.E = beam.E
            self.p = 0 
            self.tau = 0
            self.s = 0 # position along the reference trajectory
            self.id = ""

    def display(self):
        print "beta_x  = ", self.beta_x
        print "beta_y  = ", self.beta_y
        print "alpha_x = ", self.alpha_x
        print "alpha_y = ", self.alpha_y
        print "gamma_x = ", self.gamma_x
        print "gamma_y = ", self.gamma_y
        print "Dx      = ", self.Dx
        print "Dy      = ", self.Dy
        print "Dxp     = ", self.Dxp
        print "Dyp     = ", self.Dyp
        print "mux     = ", self.mux
        print "muy     = ", self.muy
        print "nu_x    = ", self.mux/2./pi
        print "nu_y    = ", self.muy/2./pi

            
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
        print ("sigma_E/E and Dx/Dy : ", self.sigma_E/self.E, "  and ", self.Dx, "/",self.Dy, " m")
        print ("emit_x/emit_y     : ",  self.emit_x*1e9, "/",self.emit_y*1e9, " nm-rad")
        print ("sigma_x/y         : ", self.sigma_x*1e6, "/", self.sigma_y*1e6, " um")
        print ("sigma_xp/yp       : ", self.sigma_xp*1e6, "/", self.sigma_yp*1e6, " urad")


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

    def add(self, ct,x,y,xp,yp,z,s):
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
    '''
    array of particles of fixed size; for optimized performance
    '''
    def __init__(self, n = 0):
        self.particles = zeros(n*6)
        self.s = 0
        self.E = 0

    def rm_tails(self, xlim, ylim, px_lim, py_lim):
        '''
        comment behaviour and possibly move out of class
        '''

        x = abs(self.particles[::6])
        px = abs(self.particles[1::6])
        y = abs(self.particles[2::6])
        py = abs(self.particles[3::6])
        ind_angles = append(argwhere(px > px_lim), argwhere(py > py_lim))
        p_idxs = unique(append(argwhere(x > xlim), append(argwhere(y > ylim), append(argwhere(x != x), append(argwhere(y!= y), ind_angles)) )))
        #if len(p_idxs) != 0:
        #e_idxs = map(lambda x: append(array([]), x), array([6*p_idxs, 6*p_idxs+1, 6*p_idxs+2, 6*p_idxs+3, 6*p_idxs+4, 6*p_idxs+5]))
        e_idxs = [append([], x) for x in array([6*p_idxs, 6*p_idxs+1, 6*p_idxs+2, 6*p_idxs+3, 6*p_idxs+4, 6*p_idxs+5])]
        self.particles = delete(self.particles, e_idxs)
        return p_idxs


    def __getitem__(self, idx):
        return Particle(x = self.particles[idx*6], px=self.particles[idx*6 + 1],
                         y=self.particles[idx*6+2], py = self.particles[idx*6+3],
                         tau = self.particles[idx*6+4], p = self.particles[idx*6+5],
                         s = self.s)

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
        for i in xrange(self.size()):
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


def get_envelope(p_array):
    tws = Twiss()
    x = p_array.x()
    px = p_array.px()
    y = p_array.y()
    py = p_array.py()
    tws.x = 0.  #mean(x)
    tws.y = 0.  #mean(y)
    tws.px =0.  #mean(px)
    tws.py =0.  #mean(py)
    #print tws.x, tws.y, tws.px,tws.py
    tws.xx = mean((x-tws.x)*(x-tws.x))
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