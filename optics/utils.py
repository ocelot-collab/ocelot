'''
single crystal Bragg reflection
'''

from ocelot.optics.elements import *
from ocelot.optics.wave import *
from ocelot.optics.bragg import *
from ocelot.optics.ray import Ray, trace as trace_ray
from ocelot.gui.optics import *


class Signal(object):
    def __init__(self, n=100):
        self.t = np.linspace(-1,1, n)
        self.f = np.zeros_like(self.t, dtype=np.complex)
        self.n = n

class Signal3D(object):
    def __init__(self, n=100):
        self.t = np.linspace(-1,1, n)
        self.f = np.zeros_like(self.t, dtype=np.complex)
        self.n = n
        
    def field_on_axis(self):
        return self.fs[0]

    def field_sum_abs(self):
        return np.sum(np.abs(self.fs[:])**2)
    
    def free(self):
        pass



def read_signal(file_name, E_ref, npad = 10):
    s = Signal()
    data = np.loadtxt(file_name, dtype = complex)
    s.f = data[:,2]
    s.t = np.real(data[:,0])
    
    ''' spectrum with finer resolution '''
    s.nslice = n = len(s.f)
    s.npad = npad

    s.f_ = np.zeros((2*npad+1)*n, dtype=complex)
    s.f_[npad*n:(npad+1)*n] = s.f  
    s.f = s.f_      
    s.t = (npad+1)*np.linspace(s.t[0], s.t[-1], len(s.f)) 
    
    spec = fft.fft(s.f)
    dt = (s.t[1] - s.t[0]) * 1.e-15
    k0 = E_ref / (hbar * c)  
    s.freq_k = 2*pi*(fftfreq(len(spec), d=dt) / c )
    s.freq_k = -np.roll(s.freq_k, len(spec)/2) + k0 # take into account s/t
    s.freq_ev = s.freq_k * hbar * c
    s.sp = np.roll(spec, len(spec)/2)
    
    return s


def plot_signal(s): 
    plt.plot(s.t, np.abs(s.f))

def plot_signal_spec(s): 
    plt.plot(s.freq_ev, np.abs(s.sp))


def plot_filters(filt, f_test=None, param='tr', ax= None):
    
    if ax == None:
        f = plt.figure()
        ax = f.add_subplot(111)
    
    ax.set_xlabel('Photon Energy [ev]')
    ax2 = ax.twinx()
    plt.grid(True)
    if param == 'tr': 
        ax.set_title('Transmissivity')
        data = filt.tr
        if f_test != None: data_test = f_test.tr
    if param == 'ref': 
        ax.set_title('Reflectivity')
        data = filt.ref
        if f_test != None: data_test = f_test.ref 
    
    l1,=ax.plot(filt.ev, np.abs(data), 'bd')
    #ax2.plot(filt.ev, unfold_angles( np.angle(data)) , 'gd')
    l2,=ax2.plot(filt.ev, np.angle(data) , 'gd')
    
    plt.legend([l1,l2],['abs','phase'])
    
    if f_test != None:
            ax.plot(f_test.ev, np.abs(data_test), 'b--')
            #ax2.plot(f_test.ev, unfold_angles(np.angle(data_test)), 'g--')
            ax2.plot(f_test.ev, np.angle(data_test), 'g--')
            
def plot_spec_filt(s, filt, ax):
    ax.plot(s.freq_ev, np.abs(s.sp), 'b.')    
    tr_r, tr_i = np.real(filt.tr), np.imag(filt.tr)
    tr_mod  = np.real(np.sqrt(tr_r*tr_r + tr_i*tr_i)) #modulus of T
        
    ax.plot(filt.ev, tr_mod / np.max(tr_mod) * np.max(np.abs(s.sp)), 'r.--')
    print(s.freq_k)

    