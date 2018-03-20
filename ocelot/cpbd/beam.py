'''
definition of particles, beams and trajectories
'''
import numpy as np
from numpy import sqrt, cos, sin, std, mean
from ocelot.common.globals import *
from ocelot.common.math_op import *
# from ocelot.common.math_op import *
from ocelot.common.py_func import filename_from_path
from copy import deepcopy
#import pickle
from scipy import interpolate
from scipy.signal import savgol_filter
import logging

# logger = logging.getLogger(__name__)
_logger = logging.getLogger('ocelot.beam')
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


class Twiss:
    """
    class - container for twiss parameters
    """
    def __init__(self, beam = None):
        if beam == None:
            self.emit_x = 0.0
            self.emit_y = 0.0
            self.emit_xn = 0.0
            self.emit_yn = 0.0
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
            self.emit_xn = beam.emit_xn
            self.emit_yn = beam.emit_yn

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
    def __init__(self, x=0, xp=0, y=0, yp=0):
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
        
        self.shape = 'gaussian' # of 'flattop'

    properties = ['g','dg','emit_xn','emit_yn','p','pz','px','py']

    @property
    def g(self):
        return self.E / m_e_GeV

    @g.setter
    def g(self,value):
        self.E = value * m_e_GeV

    @property
    def dg(self):
        return self.sigma_E / m_e_GeV

    @dg.setter
    def dg(self,value):
        self.sigma_E = value * m_e_GeV

    @property
    def emit_xn(self):
        return self.emit_x * self.g

    @emit_xn.setter
    def emit_xn(self,value):
        self.emit_x = value / self.g

    @property
    def emit_yn(self):
        return self.emit_y * self.g

    @emit_yn.setter
    def emit_yn(self,value):
        self.emit_y = value / self.g

    @property
    def p(self):
        return np.sqrt(self.g**2 - 1)
    @p.setter
    def p(self,value):
        self.g = np.sqrt(value**2 + 1)
    
    @property
    def pz(self):
        return self.p / (self.xp**2 + self.yp**2 + 1)

    @pz.setter
    def pz(self,value):
        self.p = value * (self.xp**2 + self.yp**2 + 1)

    @property
    def px(self):
        return self.p * self.xp

    @px.setter
    def px(self,value):
        self.xp = value / self.p

    @property
    def py(self):
        return self.p * self.yp

    @py.setter
    def py(self,value):
        self.yp = value / self.p

    def len(self):
        return 1
        
    def to_array(self, nslice=100, window_len=None):
        
        if self.tlen in [None, 0] and window_len is None:
            raise ValueError('both self.tlen and window_len are not set')
        if window_len is None:
            if self.shape is 'gaussian':
                window_len = self.tlen * 1e-15 * speed_of_light * 6 #sigmas
            elif self.shape is 'flattop':
                window_len = self.tlen * 1e-15 * speed_of_light * 2 #fwhm
            else:
                raise ValueError('Beam() shape can be either "gaussian" or "flattop"')
        
        beam_arr = BeamArray(nslice)
        for param in beam_arr.params():
            if hasattr(self,param) and len(getattr(beam_arr,param)) == nslice:
                setattr(beam_arr, param, np.ones(nslice) * getattr(self,param))
        beam_arr.s = np.linspace(0,window_len,nslice)
        
        if self.tlen not in [None, 0, np.inf]:
            beam_slen = self.tlen * 1e-15 * speed_of_light
            Ipeak_pos = (np.amax(beam_arr.s) - np.amin(beam_arr.s)) / 2
            if self.shape is 'gaussian':
                beam_arr.I = self.I * np.exp(-(beam_arr.s - Ipeak_pos)**2 / (2 * beam_slen**2))
            elif self.shape is 'flattop':
                beam_arr.I = np.ones_like(beam_arr.s) * self.I
                beam_arr.I[abs(beam_arr.s - Ipeak_pos) > beam_slen/2] = 0
        return beam_arr


    #def sizes(self):
    #    if self.beta_x != 0:
    #        self.gamma_x = (1. + self.alpha_x**2)/self.beta_x
    #    else:
    #        self.gamma_x = 0.
#
    #    if self.beta_y != 0:
    #        self.gamma_y = (1. + self.alpha_y**2)/self.beta_y
    #    else:
    #        self.gamma_y = 0.
#
    #    self.sigma_x = np.sqrt((self.sigma_E/self.E*self.Dx)**2 + self.emit_x*self.beta_x)
    #    self.sigma_y = np.sqrt((self.sigma_E/self.E*self.Dy)**2 + self.emit_y*self.beta_y)
    #    self.sigma_xp = np.sqrt((self.sigma_E/self.E*self.Dxp)**2 + self.emit_x*self.gamma_x)
    #    self.sigma_yp = np.sqrt((self.sigma_E/self.E*self.Dyp)**2 + self.emit_y*self.gamma_y)

    # def print_sizes(self):
        # self.sizes()
        # print("sigma_E/E and Dx/Dy : ", self.sigma_E/self.E, "  and ", self.Dx, "/",self.Dy, " m")
        # print("emit_x/emit_y     : ",  self.emit_x*1e9, "/",self.emit_y*1e9, " nm-rad")
        # print("sigma_x/y         : ", self.sigma_x*1e6, "/", self.sigma_y*1e6, " um")
        # print("sigma_xp/yp       : ", self.sigma_xp*1e6, "/", self.sigma_yp*1e6, " urad")

class BeamArray(Beam):

    def __init__(self, nslice=0):
        super().__init__()
        # initial conditions
        self.s = np.zeros(nslice)
        self.x = np.zeros(nslice)      #[m]
        self.y = np.zeros(nslice)      #[m]
        self.xp = np.zeros(nslice)    # xprime [rad]
        self.yp = np.zeros(nslice)    # yprime [rad]

        self.E = np.zeros(nslice)            # electron energy [GeV]
        self.sigma_E = np.zeros(nslice)      # Energy spread [GeV]
        self.I = np.zeros(nslice)            # beam current [A]
        self.emit_x = np.zeros(nslice)       # horizontal emittance [m rad]
        self.emit_y = np.zeros(nslice)       # horizontal emittance [m rad]
        # self.tlen = 0.0         # bunch length (rms) in fsec

        # twiss parameters
        self.beta_x = np.zeros(nslice)
        self.beta_y = np.zeros(nslice)
        self.alpha_x = np.zeros(nslice)
        self.alpha_y = np.zeros(nslice)
        self.Dx = np.zeros(nslice)
        self.Dy = np.zeros(nslice)
        self.Dxp = np.zeros(nslice)
        self.Dyp = np.zeros(nslice)

        self.eloss = np.zeros(nslice)
        
        del self.shape #inherited, not applicable

    def idx_max(self):
        idx = np.where(self.I == np.nanmax(self.I))[0]
        if len(idx) == 1:
            return idx[0]
        else:
            return idx[int(len(idx)/2)]

    def len(self):
        return np.size(self.s)

    def params(self):
        l = self.len()
        attrs=[]
        for attr in dir(self):
            if attr.startswith('__') or attr in self.properties:
                continue
            # if callable(getattr(self,attr)):
                # print('check')
                # continue
            if np.size(getattr(self,attr)) == l:
                attrs.append(attr)
        return attrs

    def sort(self):
        inds = self.s.argsort()
        for attr in self.params():
            values = getattr(self,attr)
            setattr(self,attr,values[inds])

    def equidist(self):
        dsarr = (self.s - np.roll(self.s,1))[1:]
        dsm = np.mean(dsarr)
        if (np.abs(dsarr-dsm)/dsm > 1/1000).any():
            s_new = np.linspace(np.amin(self.s), np.amax(self.s), self.len())
            for attr in self.params():
                if attr is 's':
                    continue
                #print(attr)
                val = getattr(self,attr)
                val = np.interp(s_new, self.s, val)
                setattr(self,attr,val)
            self.s = s_new
        self.ds = dsm

    def smear(self,sw):
        print('smearing the beam by {:.2e} m'.format(sw))
        self.equidist()
        sn = (sw /self.ds).astype(int)
        if sn<2:
            return
        if not sn%2:
            sn += 1

        for attr in self.params():
            if attr is 's':
                continue
            val = getattr(self,attr)
            val = savgol_filter(val,sn,2,mode='nearest')
            
            if attr in ['E', 'I', 'beta_x', 'beta_y', 'emit_x', 'emit_y', 'sigma_E']:
                # print('attribute {:s} < 0, setting to 0'.format(attr))
                val[val < 0] = 0
        #    val = convolve(val,spike,mode='same')
            setattr(self,attr,val)

    def get_s(self, s):
        idx = find_nearest_idx(self.s, s)
        return self[idx]
    
    def get_E(self):
        return self.E[self.idx_max()]
    
    def center(self, s):
        beam_s = self.get_s(s)
        self.x -= beam_s.x
        self.xp -= beam_s.xp
        self.y -= beam_s.y
        self.yp -= beam_s.yp

    def __getitem__(self,index):
        l = self.len()
        if index.__class__ is not slice:
            if index > l:
                raise IndexError('slice index out of range')
            beam_slice = Beam()
        else:
            beam_slice = deepcopy(self)

        for attr in dir(self):
            if attr.startswith('__') or callable(getattr(self,attr)):
                continue
            value = getattr(self,attr)
            if np.size(value) == l:
                setattr(beam_slice,attr,value[index])
            else:
                setattr(beam_slice,attr,value)
        return beam_slice

    def __delitem__(self,index):
        l = self.len()
        for attr in self.params():
            if attr.startswith('__') or callable(getattr(self, attr)):
                continue
            value = getattr(self, attr)
            if np.size(value) == l:
                setattr(self, attr, np.delete(value, index))

    def pk(self):
        return self[self.idx_max()]
        
    def add_chirp(self, chirp=0):
        '''
        adds linear energy chirp to the beam
        chirp = dE/E/s[um] = dg/g/s[um]
        '''
        if chirp not in [None, 0]:
            center = (np.amax(self.s) - np.amin(self.s)) / 2
            E_center = self.E[find_nearest_idx(self.s, center)]
            self.E += (self.s - center) * chirp * E_center * 1e6
        
    def to_array(self, *args, **kwargs):
        raise NotImplementedError('Method inherited from Beam() class, not applicable for BeamArray objects')


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
        self.rparticles = np.zeros((6, n))#np.transpose(np.zeros(int(n), 6))
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
        ind_angles = np.append(np.argwhere(px > px_lim), np.argwhere(py > py_lim))
        p_idxs = np.unique(np.append(np.argwhere(x > xlim), np.append(np.argwhere(y > ylim), np.append(np.argwhere(x != x), np.append(np.argwhere(y!= y), ind_angles)) )))
        #e_idxs = [append([], x) for x in array([6*p_idxs, 6*p_idxs+1, 6*p_idxs+2, 6*p_idxs+3, 6*p_idxs+4, 6*p_idxs+5])]
        self.rparticles = np.delete(self.rparticles, p_idxs, axis=1)
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
        self.rparticles = np.zeros((6, len(p_list)))
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
    def px(self): return self.rparticles[1] # xp
    def y(self):  return self.rparticles[2]
    def py(self): return self.rparticles[3] # yp
    def tau(self):return self.rparticles[4]
    def p(self):  return self.rparticles[5]
    
    @property
    def t(self):
        return self.rparticles[4]
    
    @t.setter
    def t(self,value):
        self.rparticles[4] = value


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
    tws.x = np.mean(x)
    tws.y = np.mean(y)
    tws.px =np.mean(px)
    tws.py =np.mean(py)

    if ne_flag:
        tw_x = tws.x
        tw_y = tws.y
        tw_px = tws.px
        tw_py = tws.py
        tws.xx =  np.mean(ne.evaluate('(x - tw_x) * (x - tw_x)'))
        tws.xpx = np.mean(ne.evaluate('(x - tw_x) * (px - tw_px)'))
        tws.pxpx =np.mean(ne.evaluate('(px - tw_px) * (px - tw_px)'))
        tws.yy =  np.mean(ne.evaluate('(y - tw_y) * (y - tw_y)'))
        tws.ypy = np.mean(ne.evaluate('(y - tw_y) * (py - tw_py)'))
        tws.pypy =np.mean(ne.evaluate('(py - tw_py) * (py - tw_py)'))
    else:
        tws.xx = np.mean((x - tws.x)*(x - tws.x))
        tws.xpx = np.mean((x-tws.x)*(px-tws.px))
        tws.pxpx = np.mean((px-tws.px)*(px-tws.px))
        tws.yy = np.mean((y-tws.y)*(y-tws.y))
        tws.ypy = np.mean((y-tws.y)*(py-tws.py))
        tws.pypy = np.mean((py-tws.py)*(py-tws.py))
    tws.p = np.mean( p_array.p())
    tws.E = np.copy(p_array.E)
    #tws.de = p_array.de

    tws.emit_x = np.sqrt(tws.xx*tws.pxpx-tws.xpx**2)
    tws.emit_y = np.sqrt(tws.yy*tws.pypy-tws.ypy**2)
    #print tws.emit_x, np.sqrt(tws.xx*tws.pxpx-tws.xpx**2), tws.emit_y, np.sqrt(tws.yy*tws.pypy-tws.ypy**2)
    tws.beta_x = tws.xx/tws.emit_x
    tws.beta_y = tws.yy/tws.emit_y
    tws.alpha_x = -tws.xpx/tws.emit_x
    tws.alpha_y = -tws.ypy/tws.emit_y
    return tws

def get_current(p_array, charge=None, num_bins = 200):
    """
    Function calculates beam current from particleArray
    :param p_array: particleArray
    :param charge: - None, charge of the one macro-particle.
                    If None, charge of the first macro-particle is used
    :param num_bins: number of bins
    :return s, I -  (np.array, np.array) - beam positions [m] and currents in [A]
    """
    if charge == None:
        charge = p_array.q_array[0]
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
    a = np.sqrt(-2*np.log( (1-u)) * emit)
    x = a * np.sqrt(beta) * cos(phi)
    xp = -a / np.sqrt(beta) * ( sin(phi) + alpha * cos(phi) )
    return (x, xp)

def waterbag_from_twiss(emit, beta, alpha):
    phi = 2*pi * np.random.rand()
    a = np.sqrt(emit) * np.random.rand()
    x = a * np.sqrt(beta) * cos(phi)
    xp = -a / np.sqrt(beta) * ( sin(phi) + alpha * cos(phi) )
    return (x, xp)

def ellipse_from_twiss(emit, beta, alpha):
    phi = 2*pi * np.random.rand()
    #u = np.random.rand()
    #a = np.sqrt(-2*np.log( (1-u)) * emit)
    a = np.sqrt(emit)
    x = a * np.sqrt(beta) * cos(phi)
    xp = -a / np.sqrt(beta) * ( sin(phi) + alpha * cos(phi) )
    return (x, xp)


def moments(x, y, cut=0):
    n = len(x)
    #inds = np.arange(n)
    mx = np.mean(x)
    my = np.mean(y)
    x = x - mx
    y = y - my
    x2 = x*x
    mxx = np.sum(x2)/n
    y2 = y*y
    myy = np.sum(y2)/n
    xy = x*y
    mxy = np.sum(xy)/n

    emitt = np.sqrt(mxx*myy - mxy*mxy)

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
        emitt = np.sqrt(mxx*myy - mxy*mxy)
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
    M[0, 0] = np.sqrt(b2/b1)*(cosp+a1*sinp)
    M[0, 1] = np.sqrt(b2*b1)*sinp
    M[1, 0] = ((a1-a2)*cosp-(1+a1*a2)*sinp)/np.sqrt(b2*b1)
    M[1, 1] = np.sqrt(b1/b2)*(cosp-a2*sinp)
    return M

def beam_matching(particles, bounds, x_opt, y_opt):
    pd = np.zeros(( int(len(particles)/6), 6))
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
    #print(beta, alpha)
    M = m_from_twiss([alpha, beta, 0], x_opt)
    #print(M)
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
    """
    Beam matching
    """
    def __init__(self, tws=None, x_opt=None, y_opt=None):
        """
        :param tws : Twiss object
        :param x_opt (obsolete): [alpha, beta, mu (phase advance)]
        :param y_opt (obsolete): [alpha, beta, mu (phase advance)]
        """
        self.bounds = [-5, 5]  # [start, stop] in sigmas
        self.tws = tws       # Twiss
        self.x_opt = x_opt   # [alpha, beta, mu (phase advance)]
        self.y_opt = y_opt   # [alpha, beta, mu (phase advance)]
        self.step = 1

    @property
    def twiss(self):
        if self.tws == None:
            _logger.warning("BeamTransform: x_opt and y_opt are obsolete, use Twiss")
            tws = Twiss()
            tws.alpha_x, tws.beta_x, tws.mux = self.x_opt
            tws.alpha_y, tws.beta_y, tws.muy = self.y_opt
        else:
            tws = self.tws
        return tws


    def prepare(self, lat):
        pass

    def apply(self, p_array, dz):
        _logger.debug("BeamTransform: apply")
        self.x_opt = [self.twiss.alpha_x, self.twiss.beta_x, self.twiss.mux]
        self.y_opt = [self.twiss.alpha_y, self.twiss.beta_y, self.twiss.muy]
        self.beam_matching(p_array.rparticles, self.bounds, self.x_opt, self.y_opt)

    def beam_matching(self, particles, bounds, x_opt, y_opt):
        pd = np.zeros((int(particles.size / 6), 6))
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
        mxx = np.sum(x2) / n
        y2 = y * y
        myy = np.sum(y2) / n
        xy = x * y
        mxy = np.sum(xy) / n

        emitt = np.sqrt(mxx * myy - mxy * mxy)

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
            emitt = np.sqrt(mxx * myy - mxy * mxy)
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
    '''
    returns:
    <x>,<xs>,<x^2>,<x*xs>,<xs^2>,np.sqrt(<x^2> * <xs^2> - <x*xs>^2)
    based on M particles in moving window
    Sergey, check please
    '''
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
        dq = n2 - n1 # window size
        mx[i] = (xc[n2] - xc[n1])/dq # average for over window per particle
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
    _logger.debug("slice_analysis_transverse: charge = " + str(q1))
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
    #print("charge", q1)
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
    mE = mEs #mean energy
    sE = np.sqrt(mEsEs) #energy spread
    sig0 = np.std(parray.tau()) # std pulse duration
    B = s_to_cur(z, Mcur*sig0, q1, speed_of_light)
    gamma0 = parray.E/m_e_GeV
    _, _, _, _, _, emitty0 = moments(PD[2], PD[3])
    emityn = emitty0*gamma0
    _, _, _, _, _, emitt0 = moments(PD[0], PD[1])
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

'''
beam funcions proposed
'''


def parray2beam(parray, step=1e-7):
    '''
    reads ParticleArray()
    returns BeamArray()
    step [m] - long. size ob bin to calculate distribution parameters
    '''
    _logger.info('calculating electron beam distribution from particle array')

    part_c = parray.q_array[0] #fix for general case  # charge per particle
    t_step = step / speed_of_light
    t = parray.tau() / speed_of_light
    t_min = min(t)
    t_max = max(t)
    t_window = t_max - t_min
    npoints = int(t_window / t_step)
    t_step = t_window / npoints
    beam = BeamArray()
    for parm in ['I',
                 's',
                 'emit_x',
                 'emit_y',
                 'beta_x',
                 'beta_y',
                 'alpha_x',
                 'alpha_y',
                 'x',
                 'y',
                 'xp',
                 'yp',
                 'E',
                 'sigma_E',
                 ]:
        setattr(beam, parm, np.zeros((npoints - 1)))

    for i in range(npoints - 1):
        indices = (t > t_min + t_step * i) * (t < t_min + t_step * (i + 1))
        beam.s[i] = (t_min + t_step * (i + 0.5)) * speed_of_light

        if np.sum(indices) > 2:
            e0 = parray.E * 1e9
            p0 = np.sqrt( (e0**2 - m_e_eV**2) / speed_of_light**2 )
            p = parray.rparticles[5][indices] # deltaE / average_impulse / speed_of_light
            dist_e = (p * p0 * speed_of_light + e0)
            dist_x = parray.rparticles[0][indices]
            dist_y = parray.rparticles[2][indices]
            dist_xp = parray.rparticles[1][indices]
            dist_yp = parray.rparticles[3][indices]
            
            beam.I[i] = np.sum(indices) * part_c / t_step
            beam.E[i] = np.mean(dist_e) * 1e-9
            beam.sigma_E[i] = np.std(dist_e) * 1e-9
            
            dist_x_m = np.mean(dist_x)
            dist_y_m = np.mean(dist_y)
            dist_xp_m = np.mean(dist_xp)
            dist_yp_m = np.mean(dist_yp)
            
            beam.x[i] = dist_x_m
            beam.y[i] = dist_y_m
            # g = beam.E[i] / m_e_GeV
            # p = np.sqrt(g**2 - 1)
            beam.xp[i] = dist_xp_m
            beam.yp[i] = dist_yp_m
            
            dist_x -= dist_x_m
            dist_y -= dist_y_m
            dist_xp -= dist_xp_m
            dist_yp -= dist_yp_m
            
            beam.emit_x[i] = np.sqrt(np.mean(dist_x**2) * np.mean(dist_xp**2) - np.mean(dist_x * dist_xp)**2)
            beam.emit_y[i] = np.sqrt(np.mean(dist_y**2) * np.mean(dist_yp**2) - np.mean(dist_y * dist_yp)**2)
            beam.beta_x[i] = np.mean(dist_x**2) / beam.emit_x[i]
            beam.beta_y[i] = np.mean(dist_y**2) / beam.emit_y[i]
            beam.alpha_x[i] = -np.mean(dist_x * dist_xp) / beam.emit_x[i]
            beam.alpha_y[i] = -np.mean(dist_y * dist_yp) / beam.emit_y[i]

    idx = np.where(np.logical_or.reduce((beam.I == 0, beam.E == 0, beam.beta_x > np.mean(beam.beta_x) * 100, beam.beta_y > np.mean(beam.beta_y) * 100)))
    del beam[idx]

    if hasattr(parray,'filePath'):
        beam.filePath = parray.filePath + '.beam'
    return(beam)


def generate_beam(E, I=5000, l_beam=3e-6, **kwargs):
    '''
    Generates BeamArray object
    accepts arguments with the same names as BeamArray().parameters()
    I - current in Amps
    E - beam ebergy in GeV
    emit_x, emit_n(both normalized), emit_xn, etc.
    shape - beam shape ('gaussian' of 'flattop')
    l_beam [m] - beam length in meters
    l_window [m] - window length in um
        by default: l_beam * 2 if flattop,
                    l_beam * 6 if gaussian,
    nslice - number of slices in the beam
    '''
    _logger.info('generating electron beam distribution')
    
    beam = Beam()
    beam.E = E
    beam.tlen = l_beam / speed_of_light * 1e15
    beam.I = I
    nslice=100

    for key, value in kwargs.items():
        if (key in beam.__dict__ or key in beam.properties) and (key not in ['s', 'E', 'tlen', 'I']):
            setattr(beam, key, value)
        if key is 'emit':
            beam.emit_x = value
            beam.emit_y = value
        if key is 'emit_n':
            beam.emit_xn = value
            beam.emit_yn = value
        if key is 'beta':
            beam.beta_x = value
            beam.beta_y = value
        if key is 'nslice':
            nslice = value
            
    if 'l_window' not in kwargs:
        if beam.shape is 'gaussian':
            l_window = l_beam * 6
        elif beam.shape is 'flattop':
            l_window = l_beam * 2
        else:
            raise ValueError('Beam() shape can be either "gaussian" or "flattop"')
    else:
        l_window = kwargs['l_window']
    
    beam_arr = beam.to_array(nslice, l_window)
    
    if 'chirp' in kwargs:
        beam_arr.add_chirp(kwargs['chirp'])
    
    return beam_arr