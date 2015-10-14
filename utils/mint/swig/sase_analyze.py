'''
sase statistic analysis
'''
#import ocelot.utils.mint.mint as mint
#import ocelot.utils.mint.swig.dcs as dcs
import os, sys

from pylab import *
from scipy.optimize import *
from time import sleep
from pickle import dump, load


from time import *

from numpy.fft import *

#print localtime(time()).tm_min, localtime(time()).tm_sec, localtime(time())


def filter_sase(t, gmd):
    sw = fft(gmd)
    w = fftfreq(len(gmd), t[1] - t[0])

    f = zeros_like(sw, dtype=complex)    
    sw_filt = zeros_like(sw, dtype=complex)
    hf_f = zeros_like(sw, dtype=complex)    
    hf_sw = zeros_like(sw, dtype=complex)

    
    
    width = 0.0002
    m = 0.001

    for i in range(len(w)):
        if abs(w[i]) < m: f[i] = 1.0
        elif w[i] >= m: f[i] = exp(-(w[i]-m)**2 / width**2)
        elif w[i] <= m: f[i] = exp(-(w[i]+m)**2 / width**2)
    
        sw_filt[i] = sw[i] * f[i]


    hf_width = 0.001
    hf_m = 0.005
    hf_m0 = 0.02

    for i in range(len(w)):
        if abs(w[i] - hf_m0) < hf_m or abs(w[i] + hf_m0) < hf_m: hf_f[i] = 1.0
        else: hf_f[i] = 0.0
        
        hf_sw[i] = sw[i] * hf_f[i]

    '''
    plt.figure()
    plt.plot(w, abs(hf_f), 'g--')
    plt.plot(w, abs(sw), 'b-')
    plt.plot(w, abs(hf_sw), 'r-')
    plt.show()
    sys.exit(0)
    '''
    
    filt_signal = ifft(sw_filt)
    hw_signal = ifft(hf_sw)
    return filt_signal, hw_signal, sw
    

lines = open('data/jan2015/sase.29jan.log').read().split('\n')

start_time = None

t = []
gmd = []
pyro = []


t_part = []
gmd_part = []

part_start_time = 0


for l in lines:
    tks = l.split()
    if len(tks) < 3: continue
    if start_time == None: start_time = float(tks[0]) 
    #print start_time   
    t.append(float(tks[0]) - start_time)
    gmd.append(float(tks[1]))
    pyro.append(float(tks[2]))

    if float(tks[0]) - start_time > part_start_time:
        t_part.append(float(tks[0]) - start_time)
        gmd_part.append(float(tks[1]))
        
    

fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(t,gmd)
plt.plot(t_part,gmd_part, 'g', lw=2)

gmd_filt, hf_sig, sw_filt = filter_sase(t_part, gmd_part)
plt.plot(t_part,gmd_filt, 'r--', lw=2)
plt.plot(t_part,gmd_filt + hf_sig, 'b--', lw=2)
#plt.plot(t_part,gmd_part - 1.0*hf_sig, color='#000000', lw=2)
ax.set_xlabel('Time [s]')
ax.set_ylabel(r'$[\mu J]$')

sw = fft(gmd_part)
w = fftfreq(len(gmd_part), t[1] - t[0])
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(w, abs(sw), 'gd-', lw=2)
plt.plot(w, abs(sw_filt), 'r--', lw=2)
ax.set_xlabel('Frequency [Hz]')

fig = plt.figure()
ax = fig.add_subplot(111)
sw = fft(pyro)
w = fftfreq(len(pyro), t[1] - t[0])
plt.plot(w, abs(sw), 'gd-', lw=2)
ax.set_xlabel('Frequency [Hz]')


plt.show()
