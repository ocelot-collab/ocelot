'''
single crystal Bragg reflection
'''

# from ocelot.optics.wave import TransferFunction
from ocelot.optics.bragg import c, hbar
from ocelot.optics.ray import Ray, trace as trace_ray
import matplotlib.pyplot as plt
import numpy as np
import numpy.fft as fft
from ocelot.common.globals import pi
from ocelot.optics.wave import calc_ph_sp_dens


class Signal(object):
    def __init__(self, n=100):
        self.t = np.linspace(-1,1, n)
        self.f = np.zeros_like(self.t, dtype=complex)
        self.n = n

class Signal3D(object):
    def __init__(self, n=100):
        self.t = np.linspace(-1,1, n)
        self.f = np.zeros_like(self.t, dtype=complex)
        self.n = n
        
    def field_on_axis(self):
        return self.fs[0]

    def field_sum_abs(self):
        return np.sum(np.abs(self.fs[:])**2)
    
    def free(self):
        pass



def read_signal(file_name, E_ref: float, npad: int = 10) -> Signal:
    """Read time [fs] and complex field from columns 1 and 3 of a text file.

    Keep the time step when zero-padding each end by npad times the input
    length. An envelope proportional to exp(-i omega t) has photon energy
    E_ref + hbar omega [eV].
    """
    data = np.loadtxt(file_name, dtype=complex)
    times = data[:, 0].real
    dt = times[1] - times[0]
    if dt <= 0 or not np.allclose(np.diff(times), dt, rtol=1e-7, atol=0):
        raise ValueError("Signal times must increase with uniform spacing [fs]")
    s = Signal()
    s.nslice = len(times)
    s.npad = npad
    padding = npad * s.nslice
    s.f_ = np.pad(data[:, 2], padding)
    s.f = s.f_
    s.t = times[0] + (np.arange(len(s.f)) - padding) * dt
    frequencies = fft.fftshift(fft.fftfreq(len(s.f), d=dt * 1e-15))
    s.freq_k = E_ref / (hbar * c) - 2 * pi * frequencies / c
    s.freq_ev = s.freq_k * hbar * c
    s.sp = fft.fftshift(fft.fft(s.f))
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
