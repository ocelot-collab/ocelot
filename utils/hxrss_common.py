'''
FEL self-seeding
'''

from ocelot.optics.elements import *
from ocelot.optics.wave import *
from ocelot.optics.bragg import *
from ocelot.optics.ray import Ray, trace as trace_ray
from ocelot.gui.optics import *
from ocelot.common.math_op import peaks
from ocelot.adaptors.genesis import *
#from sim_info import *

from copy import deepcopy
import multiprocessing

#import fftw3
#import sharedmem as shm
import time
import copy
import scipy.fftpack as fft

from ocelot.optics.utils import *
import mpi4py


#import numpy.fft as fft



'''
configurable to e.g. semi-empirical models
'''
class FelSimulator(object):
    
    def __init__(self):
        self.engine = 'genesis'
    
    def run(self):
        if self.engine == 'test_1d':
            w1 = read_signal(file_name=self.input, npad = self.npad , E_ref = self.E_ev)
            return w1, None
        if self.engine == 'test_3d':
            ''' produced  sliced field '''
            w1 = read_signal(file_name=self.input, npad = self.npad , E_ref = self.E_ev)
            s3d = Signal3D()
            s3d.fs = [w1, deepcopy(w1), deepcopy(w1), deepcopy(w1)]
            s3d.mesh_size = (2,2)            
            return s3d, None
        if self.engine == 'test_genesis':
            ''' read test sliced field '''
            g = read_genesis_output(self.input)
            print ('read sliced field ', g('ncar'), g.nSlices)
            slices = readRadiationFile(fileName=self.input + '.dfl', npoints=g('ncar'))            
            s3d = Signal3D()
            s3d.slices = slices
            s3d.mesh_size = (int(g('ncar')),int(g('ncar'))) 
            s3d.g = g           
            return s3d, None



def pulses_from_field(pulse3d, range=None, npad = 2, threaded = False):
    import gc
    gc.collect()
    nx, ny = pulse3d.nx, pulse3d.ny
    n_slices = len(pulse3d.slices) / (nx*ny)
    n_pulses = nx*ny
    
    n = (2*npad+1)*n_slices
    
    tmax = pulse3d.tmax
    E_ref = pulse3d.E_ref
    k0 = E_ref / (hbar * c)
    print ('creating ', n_pulses, ' pulses ', nx, 'x', ny, ' tmax=', tmax)
    
    pulses = Signal()
             
    ''' spectrum with finer resolution '''
    pulses.nslice = n_slices
    pulses.npad = npad
    pulses.n = n
    pulses.n_pulses = n_pulses
    pulses.t = (npad+1)*np.linspace(0, tmax, n)
    
    dt = (pulses.t[1] - pulses.t[0]) * 1.e-15
    pulses.freq_k = 2*pi*(fftfreq(n, d=dt) / c )            
    pulses.freq_k = -np.roll(pulses.freq_k, n/2) + k0
    pulses.freq_ev = pulses.freq_k * hbar * c

    if threaded:
        pulses.f = shm.empty((n_pulses,n,), np.complex)
    else:
        pulses.f = np.zeros([n_pulses,n], dtype=complex)
     
    for i in range(n_pulses): 
        pulses.f[i,npad*n_slices:(npad+1)*n_slices] = pulse3d.slices[i::nx*ny]        

    #del(pulse3d.slices)
    return pulses


def field_from_pulses_test(t, fs):
    s3d = Signal3D()
    s3d.t = t
    s3d.fs = fs
    return s3d

def field_from_pulses(t, fs, mesh_size, slices=None, write_file = None):
    
    nslice = fs.shape[1]
    print (fs.shape)
    print ('creating field', nslice , 'x', mesh_size[0], 'x', mesh_size[1])
    if slices == None: 
        slices = np.zeros([nslice, mesh_size[0], mesh_size[1]], dtype=complex)
        
    for i in range(nslice):
        for j in range(mesh_size[0]):
            for k in range(mesh_size[1]):
                #print i, j, k, j*mesh_size[0] + k, fs[j*mesh_size[0] + k][i]
                slices[i,j,k] = fs[j*mesh_size[0] + k,i]
                
    s3d = Signal3D()
    s3d.t = t
    s3d.slices = slices
    
    return s3d

'''
modulus of wake signal (fel pulse minus core), for seed delay finding
f -- modulus of the original signal (time domain)
'''
def get_wake(f, n_fel_start, smooth_param=120):
    f2 = copy(f)
    fact = 1.0
    #d = np.abs(np.diff(f2))
    for i in range(1,len(f2)):
        if i >= n_fel_start:
            fact = 0.0
        f2[i] *= fact
    f2 = np.convolve(f2, np.ones(smooth_param) / float(smooth_param), mode='same')
    return f2


def update_beam(beam_new, g, beam):
    g0 = np.array(map(lambda x : g.sliceValues[x]['energy'][-1], range(1,g.nSlices+1)) )
    dg = np.array(map(lambda x : g.sliceValues[x]['e-spread'][-1], range(1,g.nSlices+1)) )
    I = np.array(g.I)

    i_diff = len(g0) - len(beam_new.z)
    if i_diff <= 0: i_diff = -len(g0)

    print ('i_diff', i_diff)

    g0 = g0[:-i_diff]
    dg = dg[:-i_diff]
    I = I[:-i_diff]
    beam_new.z = beam_new.z[:-i_diff]

    '''
    plt.figure()
    plt.plot(beam_new.g0)
    plt.plot(g0 + beam.gamma_rel)
    plt.figure()
    plt.plot(beam_new.dg)
    plt.plot(dg)
    plt.figure()
    plt.plot(I)
    plt.plot(beam_new.I)
    plt.show()
    '''

    beam_new.g0 = g0 + beam.gamma_rel
    beam_new.dg = dg


def log_info(sim_info, g, run_id, stage):
    
    r = RunInfo(run_id)
    r.max_power = g.max_power
    r.power = g.power
    r.power_z = g.power_z
    r.stage = stage
    r.z = g.z
    r.t = g.t
    r.spec = g.spec
    r.freq_ev = g.freq_ev

    sim_info.runs[r.id] = r


    print ('saving run ', id)
    f_obj = open(sim_info.log_dir + 'dump.dat', 'wb')
    pickle.dump(sim_info, f_obj)
    f_obj.close()

def log_info_seed(sim_info, g, run_id, stage):
    r = RunInfo(run_id)
    r.max_power = g.max_power
    r.power_z = []
    r.power = g.power
    r.spec = g.spec
    r.freq_ev = g.freq_ev 
    r.power_ref = g.power_ref
    r.stage = stage
    sim_info.runs[r.id] = r
    f_obj = open(sim_info.log_dir + 'dump.dat', 'wb')
    pickle.dump(sim_info, f_obj)
    f_obj.close()


def filter_1d(pulse, transm, i):
    
    n_p = len(pulse.f[i,:])
    if n_p % 100 != 0:
        n_p1 = n_p - n_p % 100 
        pulse.f[i,n_p1: n_p] = 0.0
    else:
        n_p1 = n_p
    
    sp = np.roll(fft.fft(pulse.f[i,0:n_p1]), pulse.n/2)
    sp = sp * transm[0:n_p1]
    sp = np.roll( sp, pulse.n/2)
    pulse.f[i,0:n_p1] = fft.ifft(sp)


def sseed(input_file, E_ev, chicane, run_dir, delay = None, debug=True, 
          output_file = None, wake=None, xt_couple=False, n_peak = 1, 
          npad = 6, threaded = True):

    h = 4.135667516e-15
    c = 299792458.0
    
    t1 = time.time()
    g = read_genesis_output(input_file)
    print ('read sliced field ', g('ncar'), g.nSlices)
    ncar = int(g('ncar'))
    dgrid = float(g('dgrid'))

    idx_max = np.argmax(g.I)
    g.idx_max = idx_max
    print ('ss: idx_max:', g.idx_max)

    slices = readRadiationFile(fileName=input_file + '.dfl', npoints=g('ncar'))  

    print ('field readout :', time.time() - t1, ' sec')
    t1 = time.time()

    s3d = Signal3D()
    s3d.slices = slices
    s3d.mesh_size = (int(g('ncar')),int(g('ncar'))) 
    s3d.g = g           
    pulse3d, bunch = s3d, None

    pulse3d.slices = np.reshape(pulse3d.slices, (1,-1)) # 1d array
    pulse3d_part = Signal3D()
    pulse3d_part.nx = int(pulse3d.g('ncar'))
    pulse3d_part.ny = int(pulse3d.g('ncar'))
    pulse3d_part.tmax = pulse3d.g.nSlices * pulse3d.g('zsep') * pulse3d.g('xlamds') / 2.99792458e8 * 1.e15   
    pulse3d_part.slices = pulse3d.slices[0]
    pulse3d_part.E_ref = 1./ pulse3d.g('xlamds') * 4.135667516e-15 *2.99792458e8

    g.npad = npad
    pulses_1d = pulses_from_field(pulse3d_part, npad=npad, threaded = threaded) 
    g.nslice = pulses_1d.nslice
    print ('*** created', pulses_1d.n_pulses, ' pulses *** ')


    if debug:
        pulse_idx = int(pulse3d_part.nx*pulse3d_part.ny/2.) 
        print ('plotting slice', pulse_idx)
        fig = plt.figure()
        ax = fig.add_subplot(221)
        ax.grid(True)
        ax.plot(pulses_1d.t, np.abs(pulses_1d.f[pulse_idx])+1, '#000000', alpha=0.5)
        ax.set_title('Stage 1 FEL pulse')

    pulse_idx = int(pulse3d_part.nx*pulse3d_part.ny/2.) 
    
    n_fel_start = 0
    for i in range(len(pulses_1d.f[pulse_idx])):
        if np.abs(pulses_1d.f[pulse_idx][i]) > 1.0:
            n_fel_start = i
            break

    n_points = len(pulses_1d.f[pulse_idx]) / (2*npad + 1)
    g.power_ref = np.abs(pulses_1d.f[pulse_idx])[n_points*npad:n_points*(npad+1)]

    r = Ray() 
    r.lamb = 2 * pi * hbar * c / E_ev
    print ('wavelength', r.lamb, '(', E_ev, 'eV), filetring...')


    filt = get_crystal_filter(chicane.cryst, r, ref_idx = chicane.cryst.ref_idx, k = pulses_1d.freq_k)
    print (pulses_1d.freq_k)
    chicane.cryst.filter = filt

    f_av = np.zeros(len(pulses_1d.t)) 

    print ('field/filter preparation: ', time.time() - t1, ' sec' )
    t1 = time.time()
    if threaded:
        print ('filtering (threaded) ... ')
        #ncores = multiprocessing.cpu_count()/8
        pool = multiprocessing.Pool()
        tasks = [pool.apply_async(filter_1d, (pulses_1d, chicane.cryst.filter.tr,i)) for i in range(pulses_1d.n_pulses)]
        for t in tasks:
            t.wait()
        pool.close()
        pool.join()
    else:
        print ('filtering (not threaded) ... ')
        for i in range(pulses_1d.n_pulses):
            #print 'filtering:', i, '/', len(pulses_1d.f[i,:])
            filter_1d(pulses_1d, chicane.cryst.filter.tr, i )

    for i in range(pulses_1d.n_pulses):
        f_av += np.abs(pulses_1d.f[i,:])
    
    print ('filtering time: ', time.time() - t1, ' sec')
    t1 = time.time()

    if wake != None:
        f_av = wake

    n_marg = 20 / (pulses_1d.t[1] - pulses_1d.t[0]) # number of margin slices prior to FEL pulse
    f_wake = get_wake(f_av, n_fel_start - n_marg, smooth_param=2)

    if debug:
        ax2 = fig.add_subplot(222)
        ax2.set_title('Wake')
        ax2.grid(True)
        ax2.set_yscale('log')
        ax2.plot(f_wake, 'r--')
        ax2.plot(f_av, 'g--')


    x,y = peaks(pulses_1d.t, f_wake, n=4)
    print ('peaks', x, y)

    if delay == None:
        delay = pulses_1d.t[pulses_1d.nslice*npad + g.idx_max] - x[n_peak] # based on current maximum

    n_delay = int(delay / (pulses_1d.t[1] - pulses_1d.t[0]))
    n_start = pulses_1d.npad*pulses_1d.nslice - n_delay

    print ('delay', delay, ' fs ', delay * 1.e-15 * 3.e8 / 1.e-6 , 'mu m ', n_delay , ' slices', ' nslice=', pulses_1d.nslice, ' n_start=', n_start )
    print ('max current slice: ', pulses_1d.nslice*npad + g.idx_max)

    i1=0
    i2=pulses_1d.n_pulses

    t = pulses_1d.t[n_start:n_start+pulses_1d.nslice]
    pulse3d.slices = np.reshape(pulse3d.slices, (pulses_1d.nslice, pulse3d_part.nx, pulse3d_part.nx))
    field_from_pulses( t, pulses_1d.f[i1:i2,n_start:n_start+pulses_1d.nslice], pulse3d.mesh_size,pulse3d.slices)

    g.spec = np.fft.fft(pulse3d.slices[:,ncar/2,ncar/2])
    g.max_power = np.max(np.abs(pulse3d.slices[:,ncar/2,ncar/2]))
    g.power =  np.abs(pulse3d.slices[:,ncar/2,ncar/2])
    g.freq_ev = h * fftfreq(len(g.spec), d=(t[1]-t[0])*1.e-15) 

    g.delay = delay

    g.seed_axis = np.abs(pulse3d.slices[:,ncar/2,ncar/2])
    g.wake = f_av
    nslice = len(pulse3d.slices[:,0,0])

    if xt_couple: # spatial-temporal coupling
        dct = (pulses_1d.t[1] - pulses_1d.t[0]) * 1.e-15 * c    

        for i_slice in range(len(pulse3d.slices[:,0,0])):
            #print 'shifting slice', i_slice
            shift_x = (nslice - i_slice)*dct * cos(chicane.cryst.thetaB) / sin(chicane.cryst.thetaB)
            #print 'shift_x' , shift_x
            #print dgrid, ncar
            n_shift = int( shift_x / (dgrid / ncar) )
            #print 'n_shift', n_shift
            pulse3d.slices[i_slice,:,:] = np.roll(pulse3d.slices[i_slice,:,:], n_shift, 0)
    

    if debug:
        ax.set_yscale('log')
        ax.plot(t + 0, np.abs(pulse3d.slices[:,ncar/2,ncar/2]), 'r--')
        ax.plot(t + 0, np.abs(pulse3d.slices[:,ncar/2,ncar/2+5]), 'r--')
        ax.plot(t + delay, np.abs(pulse3d.slices[:,ncar/2,ncar/2]), 'g-')
        ax.plot(pulses_1d.t[pulses_1d.nslice*npad + g.idx_max], f_wake[pulses_1d.nslice*npad + g.idx_max], 'bs')
        print ('seed t', pulses_1d.t[pulses_1d.nslice*npad + g.idx_max], f_wake[pulses_1d.nslice*npad + g.idx_max])

        ax3 = fig.add_subplot(223)
        #ax3.plot(t + delay, np.imag(pulse3d.slices[:,ncar/2,ncar/2]), 'b-', lw=1)
        print ('debug: max slice', g.idx_max,)
        ax3.plot(np.abs(pulse3d.slices[100,:,ncar/2]), 'r--')
        ax3.plot(np.abs(pulse3d.slices[200,:,ncar/2]), 'g--')
        ax3.plot(np.abs(pulse3d.slices[300,:,ncar/2]), 'b--')
        ax3.plot(np.abs(pulse3d.slices[g.idx_max,:,ncar/2]), 'b-')

        ax4 = fig.add_subplot(224)
        ax4.plot(g.freq_ev, np.abs(g.spec), 'b-')
        #plt.title('Spectrum (middle pulse after seed)')

        fig = plt.figure()
        ax5 = fig.add_subplot(111)
        ax5.set_yscale('log')
        ax5.plot(np.abs(pulse3d.slices[:,ncar/2,ncar/2]), 'r-')
        ax5.plot(g.wake[n_start:n_start+pulses_1d.nslice], 'b-')
        ax5.plot(g.I, 'g-')

        plt.show()
    
    print ('creating final field: ', time.time() - t1, ' sec')
    t1 = time.time()

    if output_file != None:
        writeRadiationFile(output_file+'.dfl', pulse3d.slices)    
        print ('written radiation file, slices', len(pulse3d.slices[:,ncar/2,ncar/2]))

    print ('writing final field: ', time.time() - t1, ' sec')

    return g

'''
########################################
#########     Added by G.G.    #########
########################################
'''

def FWHM(X,Y):
    '''
    Function name: FWHM(X,Y)
    
    Description:  
    returns the FWHM of Y(X)


    Arguments:
        -X :   abscissa
        -Y :   ordinate

    Date revised: 2013.8.8
    
    '''
    
    #from pylab import * #removed
    
    half_max = max(Y) / 2.
    print (half_max)
    #find when function crosses line half_max (when sign of diff flips)
    #take the 'derivative' of signum(half_max - Y[])
    d = sign(half_max - array(Y[0:-1])) - sign(half_max - array(Y[1:]))
    #plot(X,d) #if you are interested
    #find the left and right most indexes
    left_idx = find(d > 0)[0]
    right_idx = find(d < 0)[-1]
    print (X[left_idx])
    print (X[right_idx])
    return [left_idx, right_idx, X[right_idx] - X[left_idx]] #return the xpos, left and right and difference (full width)

def readres(namef):
    '''
    Function name: readres(namef)
    
    Description:
                  -Reads files from the Result directory. These are always 2 columns of ascii data.
        -Returns X and Y columns
        -NO comment is written into the log file
        
    Arguments:
        -namef:    file name
        
    Date revised: 2013.8.7
    
    '''
    
  

    dataX = []
    dataY = []
    f = open(namef, 'r')

    for line in f:
        line = line.strip()
        columns = line.split()
        dataX = np.append(dataX,float(columns[0]))
        dataY = np.append(dataY,float(columns[1]))
    
    f.close()
    
    return np.array([dataX,dataY])
     
'''
OLD VERSION!
def update_beam_2(beam_new, g, n_interp):
    beam = deepcopy(beam_new)
    g0 = np.array(map(lambda x : g.sliceValues[x]['energy'][-1], range(1,g.nSlices+1)) )
    dg = np.array(map(lambda x : g.sliceValues[x]['e-spread'][-1], range(1,g.nSlices+1)) )
    
    print len(g0)
    print g.nSlices
    
    print len(beam_new.z)
    
    I = np.array(g.I)
    
    if n_interp == 0: n_interp = g.nSlices
    beam_new.z = np.linspace(beam.z[0], beam.z[-1], n_interp) 
    z2 = np.linspace(beam.z[0], beam.z[-1], g.nSlices)
    beam_new.I = np.interp(beam_new.z, beam.z, beam.I)
            
    zmax, Imax = peaks(beam_new.z, beam_new.I, n=1)
    beam_new.idx_max = np.where(beam_new.z == zmax)[0][0]
            
    beam_new.ex = np.interp(beam_new.z, beam.z, beam.ex) 
    beam_new.ey = np.interp(beam_new.z, beam.z, beam.ey) 
    beam_new.zsep = beam.zsep * len(beam.z) / len(beam_new.z)
    #beam_new.g0 = np.interp(beam_new.z, beam.z, beam.g0) 
    beam_new.g0 = g0 + beam.gamma_rel
    print len(beam_new.z)
    print len(beam_new.g0)
    print len(beam.z)
    beam_new.g0 = np.interp(beam_new.z, z2, beam_new.g0)
    beam_new.dg = dg   
    beam_new.dg = np.interp(beam_new.z, z2, beam_new.dg)
        
    beam_new.eloss = np.interp(beam_new.z, beam.z, beam.eloss)
    
    beam_new.betax = np.interp(beam_new.z, beam.z, beam.betax)
    beam_new.betay = np.interp(beam_new.z, beam.z, beam.betay)
    beam_new.alphax = np.interp(beam_new.z, beam.z, beam.alphax)
    beam_new.alphay = np.interp(beam_new.z, beam.z, beam.alphay)
    
    beam_new.x = np.interp(beam_new.z, beam.z, beam.x) 
    beam_new.px = np.interp(beam_new.z, beam.z, beam.px) 
    beam_new.y = np.interp(beam_new.z, beam.z, beam.y)
    beam_new.py = np.interp(beam_new.z, beam.z, beam.py)
    
'''



def update_beam_2(beam_new, g, n_interp):
    beam = deepcopy(beam_new)
    # g0 = np.array(map(lambda x : g.sliceValues[x]['energy'][-1], range(1,g.nSlices+1)) )
    # dg = np.array(map(lambda x : g.sliceValues[x]['e-spread'][-1], range(1,g.nSlices+1)) )
    g0=g.el_energy[:,-1]# * (0.511e-3)
    dg=g.el_e_spread[:,-1]
    
    print (len(g0))
    print (g.nSlices)
    
    print (len(beam_new.z))
    
    I = np.array(g.I)
    
    if n_interp == 0: n_interp = g.nSlices
    '''
    plt.figure()
    plt.plot(beam_new.g0)
    plt.plot(g0 + beam.gamma_rel)
    plt.figure()
    plt.plot(beam_new.dg)
    plt.plot(dg)
    plt.figure()
    plt.plot(I)
    plt.plot(beam_new.I)
    plt.show()
    '''
    
    
    
    beam_new.z = np.linspace(beam.z[0], beam.z[-1], n_interp) 
    z2 = np.linspace(beam.z[0], beam.z[-1], g.nSlices)
    beam_new.I = np.interp(beam_new.z, beam.z, beam.I)
            
    zmax, Imax = peaks(beam_new.z, beam_new.I, n=1)
    beam_new.idx_max = np.where(beam_new.z == zmax)[0][0]
            
    beam_new.ex = np.interp(beam_new.z, beam.z, beam.ex) 
    beam_new.ey = np.interp(beam_new.z, beam.z, beam.ey) 
    beam_new.zsep = beam.zsep * len(beam.z) / len(beam_new.z)
    #beam_new.g0 = np.interp(beam_new.z, beam.z, beam.g0) 
    # print ("_______________________________")
    # print (g0)
    # print(beam.E)
    # print(beam.E/(0.511e-3))
    # print ("_______________________________")
    beam_new.g0 = g0 + beam.E/(0.511e-3) #potential problem here, no beam.gamma_rel
    print (len(beam_new.z))
    print (len(beam_new.g0))
    print (len(beam.z))
    beam_new.g0 = np.interp(beam_new.z, z2, beam_new.g0)
    beam_new.dg = dg   
    beam_new.dg = np.interp(beam_new.z, z2, beam_new.dg)
        
    beam_new.eloss = np.interp(beam_new.z, beam.z, beam.eloss)
    
    beam_new.betax = np.interp(beam_new.z, beam.z, beam.betax)
    beam_new.betay = np.interp(beam_new.z, beam.z, beam.betay)
    beam_new.alphax = np.interp(beam_new.z, beam.z, beam.alphax)
    beam_new.alphay = np.interp(beam_new.z, beam.z, beam.alphay)
    
    beam_new.x = np.interp(beam_new.z, beam.z, beam.x) 
    beam_new.px = np.interp(beam_new.z, beam.z, beam.px) 
    beam_new.y = np.interp(beam_new.z, beam.z, beam.y)
    beam_new.py = np.interp(beam_new.z, beam.z, beam.py)

def sseed_2(input_file, output_files, E_ev, chicane, run_dir, delay = 0.0, debug=True, output_file = None,  xt_couple=False, filterfilename='',method='hxr_wake_calc',paral=False):

   

    h = 4.135667516e-15
    c = 299792458.0
       
    g = read_genesis_output(input_file, readall=0)
    print ('read sliced field ', g('ncar'), g.nSlices)
    ncar = int(g('ncar'))
    dgrid = float(g('dgrid'))
    xlamds = float(g('xlamds'))
    #nslice = len(g.spec)
    nslice = int(g.nSlices)
    print ('nslice = ',nslice)
    zsep  = float(g('zsep'))
    print (zsep)

    k0 = 2*np.pi/xlamds    
    ds = zsep*xlamds                #interval ds
    srange = nslice*zsep*xlamds     #range in s
    
    
    
    
    
    
    dkold = 2*np.pi/srange    
    krange = dkold*nslice
    #idx_max = np.argmax(g.I)
    #s_imax = -srange/2.0 + idx_max*ds
    
    

    #SHF = np.int((delay-s_imax)/ds)
    
    SHF = int(delay)#np.int(delay/ds) #or divided ds
    
    #print 'imax = ', idx_imax
    #print 'delay = ', delay
    #print 'delay-s_imax', delay-s_imax
    #print 'SHF*ds', SHF*ds
    
    #SHF = np.int((-0.5e-5+s_imax)/ds)
    ####SHF = np.int((delay)/ds)
    #print np.int((delay+s_imax)/ds)*ds
    #print 'due=', SHF*ds
    #exit()

    #g.idx_max = idx_max
    #print 'ss idx_max:', g.idx_max    

    if method=='hxr_wake_calc':
        
        r = Ray() 
        r.lamb = 2 * pi * hbar * c / E_ev
        print ('wavelength', r.lamb, '(', E_ev, 'eV), filetring...')
        
        ref_idx = chicane.cryst.ref_idx
        filt = get_crystal_filter(chicane.cryst, r, nk=1000, ref_idx = ref_idx)
        
        chicane.cryst.filter = filt

        klpos, krpos, cwidth = FWHM(filt.k, 1.0-np.abs(filt.tr))
        cmid_idx   = int((krpos - klpos)/2.0)
        cmid       = filt.k[cmid_idx]        

        H, d, phi = find_bragg(lambd = r.lamb, lattice=chicane.cryst.lattice, ord_max = 15)
        dhkl = d[ref_idx]
        thetaB = phi[ref_idx]* np.pi / 180.0
        
        dk = cwidth/6.0 #The filter transmissivity defines this quantity by taking 5 points on the bottom of T; if not good separation, must tune!!!
        dr = ncar/(dgrid*np.tan(thetaB)) #dr prepares for inclusion of the spatiotemporal coupling; it is transverse pix size/cot(thetaB)
        #dr=1.0
        print ('#################')
        print ('dgrid = ', dgrid)
        print ('ncar = ',ncar)
        print ('dr = ',dr)
        print ('thetaB = ',thetaB)
        print ('ds = ',ds)
        print ('SHF = ',SHF)
        print ('#################')
        
        
        mult = np.int(dkold/dk)
        
         
        if int(mult/2) - mult/2 ==0: mult = mult-1
        if mult == 0: mult =1
        print ('MULT = ',mult)
        dk = dkold/mult #Important! Otherwise the np.int in the line changes the k scale substantially 

        phases = unfold_angles(np.angle(np.conj(filt.tr)))
        
        f1 = open(output_files[1], 'w')
        f2 = open(output_files[2], 'w')

        for i in range(len(filt.k)):
            f1.write('%s ' %(filt.k[i]) + '%s' %np.abs(filt.tr[i]) +'\n')
            f2.write('%s ' %(filt.k[i]) + '%s' %phases[i] +'\n')
        
        f1.close()
        f2.close()

    if method=='sxr_filter_read':
        # f = open(filterfilename, 'r')
        #[abs, dlpl]
        # import numpy as np
        # from math import pi
        #        E_ev=1000.
        #icf_path='d:\Work\!PROJECTS\ocelot_test\ICF_1000.ascii'
        f = open(filterfilename, 'r')
        data = np.genfromtxt(f, delimiter=',')
        #delete(data,0,0) # Erases the first row (i.e. the header)
        #plot(data[:,0],data[:,1],'o')
        f.close()
        
        dlpl=np.flipud(data[:,0])
        Tmod=np.flipud(data[:,1])
        
        
        Tpha=np.zeros(len(data[:,0]))
        lambda_0=1239.8/E_ev*1e-9
        k=2*pi/(lambda_0+lambda_0*dlpl)
        dk_f=k[0]-k[1]
        k=np.concatenate(([k[0]+dk_f],k,[k[-1]-dk_f]))
        Tmod=np.concatenate(([0],Tmod,[0]))
        Tpha=np.concatenate(([0],Tpha,[0]))
        # Tmod=np.insert(Tmod,slice(0),0)
        # Tmod=np.insert(Tmod,slice(-1),0)
        # Tpha=np.insert(Tmod,slice(0),0)
        # Tpha=np.insert(Tmod,slice(-1),0)
        # Tmod=np.insert(Tmod,slice(0,-1),[0,0])
        # Tpha=np.insert(Tpha,slice(0,-1),[0,0])#padding with zeros on edges, so that interpolation was done with zeros as well
        print (np.column_stack((k, Tmod)))
        #f = open(res_dir+'s2.Tmod.dat', 'w')
        #    writepath_Tmod='d:\Work\!PROJECTS\ocelot_test\s2.Tmod.dat'
        #    writepath_Tpha='d:\Work\!PROJECTS\ocelot_test\s2.Tpha.dat'
        np.savetxt(output_files[1],np.column_stack((k, Tmod)))
        np.savetxt(output_files[2],np.column_stack((k, Tpha)))
        
        klpos, krpos, cwidth = FWHM(k, Tmod)
        cwidth=abs(cwidth)
        print ('CWIDTH=',cwidth)
        dk = cwidth/10.0
        mult = np.int(dkold/dk)
        
        #if int(mult/2) - mult/2 ==0: mult = mult-1 %probably it needs to be even
        
        print ('MULT = ',mult)
        dk = dkold/mult #Important! Otherwise the np.int in the line changes the k scale substantially 
        
        # bring to common format !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        dr=0
        
        
        # phases = unfold_angles(np.angle(np.conj(filt.tr)))
        


    if paral == True:
        ARGS   =   ''.join([input_file+'.dfl'+' ', 	        
            output_files[0]+' ',                  
            output_files[1]+' ',                       
            output_files[2]+' ',                       
            output_files[3]+' ',
            output_files[4]+' ',                       
            output_files[5]+' ',                       
            output_files[6]+' ',                       
            output_files[7]+' ',
            output_files[8]+' ', 
            str(xlamds)+' ',                       
            str(ncar)+' ',                       
            str(mult)+' ',                       
            str(ds)+' ',                       
            str(dk)+' ',                       
            str(SHF)+' ',                       
            str(nslice)+' ', 
            str(dr)])
        runpar = '`which mpirun` -x PATH -x MPI_PYTHON_SITEARCH -x PYTHONPATH'
        prog   = ' '+'python /data/netapp/xfel/gianluca/products/ocelot/utils/seed.py '+ARGS+''
        #prog   = ' '+'python /data/netapp/xfel/svitozar/CODE/ocelot/utils/seed.py '+ARGS+''
        
       
        #print 'NEWARGS'
        #print ARGS
        #exit()
        cmd = runpar+prog
        os.system(cmd)
    
    if paral == False:
        import ocelot.utils.seed0 as sd #not found
        sd.filterfield(input_file+'.dfl',         
            output_files[0],                  
            output_files[1],                       
            output_files[2],                       
            output_files[3],
            output_files[4],                       
            output_files[5],                       
            output_files[6],                       
            output_files[7],
            output_files[8],
            #output_files[9],
            #output_files[10],                       
            #output_files[11],                       
            #output_files[12],                       
            #output_files[13],
            #output_files[14],
            #output_files[15], 
            #output_files[16],
            xlamds,                       
            ncar,                       
            mult,                       
            ds,                       
            dk,                       
            SHF,                       
            nslice, 
            dr)
    
    print ('reading filtered file (hxrss_common, lile 728) ', output_files[5])
    
    ssc, Pout = readres(output_files[5])
    lsc, Sout = readres(output_files[7])
    
    return ssc, Pout, lsc, Sout
