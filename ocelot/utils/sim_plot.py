import sys

from ocelot.adaptors.genesis import *

from ocelot.cpbd.elements import Element, Quadrupole, RBend, Drift, Undulator
from ocelot import MagneticLattice
from ocelot.cpbd.beam import Beam
from ocelot.cpbd.optics import *

import numpy.fft as fft
from sim_info import SimInfo, RunInfo

#params = {'backend': 'ps', 'axes.labelsize': 18, 'text.fontsize': 18, 'legend.fontsize': 18, 'xtick.labelsize': 18,  'ytick.labelsize': 18, 'text.usetex': True}
#rcParams.update(params)
#rc('text', usetex=True) # required to have greek fonts on redhat

import argparse

h = 4.135667516e-15
c = 299792458.0

parser = argparse.ArgumentParser(description='FEL simulation postprocessor')
#parser.add_argument('--submit', help='submit to main index file', action='store_true')
parser.add_argument('--path', help='path to the experiment', default='./')
parser.add_argument('--stage', help='undulator/seeding stages 1 through 5', default='1')
parser.add_argument('--range', help='range of runs in the form i1:i2')
parser.add_argument('--field_file', help='read in field file', action='store_true')

args = parser.parse_args()

run_start, run_end =   [int(i) for i in  args.range.split(':') ]


fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.set_xlabel('Time [fs]')
ax1.set_ylabel('Power [W]')

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.set_xlabel('Photon Energy [eV]')
ax2.set_ylabel('Spectrum [arb. units]')
ax2.get_xaxis().get_major_formatter().set_useOffset(False)
ax3 = ax2.twiny()
ax3.set_xlabel('Wavelength [nm]')


power_av = None
spec_av = None


runs = range(run_start, run_end+1)

for run_id in runs:
    run_dir = args.path + '/run_' + str(run_id)
    if args.stage in ['1','3','5']:
        run_file = run_dir + '/run.' + str(run_id) + '.s' + str(args.stage) + '.gout'
        if args.stage == '5' : run_file = run_dir + '/run.' + str(run_id) + '.gout'
        print ('reading', run_file)
        g = readGenesisOutput(run_file)
        field_file = run_file + '.dfl'
        if args.field_file:
            slices = readRadiationFile(fileName=field_file, npoints=g('ncar'))
            P = np.zeros_like(slices[:,0,0])
            for i in xrange(len(P)):            
                P[i] = sum( np.abs(np.multiply(slices[i,:,:], slices[i,:,:].conjugate())) )            
            t = np.linspace(g.t[0], g.t[-1], len(P))
        else:
            P = g.power_int
            t = g.t
        w_l_m  =  g('xlamds')
        w_l_ev = h * c / g('xlamds')
        x = np.roll(g.freq_ev, len(g.freq_ev)/2)+ w_l_ev
        y = np.roll( np.abs(g.spec)**2, len(g.freq_ev)/2)
    else:
        run_file = run_dir + '/run.' + str(run_id) + '.s' + str( int(args.stage) - 1) + '.gout'
        field_file = 'tmp' + str(args.stage) + '.dfl'
        print ('reading', run_file, 'and', field_file)
        g = readGenesisOutput(run_file)
        slices = readRadiationFile(fileName=run_dir + '/' + field_file, npoints=g('ncar'))
        P = np.zeros_like(slices[:,0,0])
        spec = np.zeros_like(slices[:,0,0])
        for i in range(len(P)):
            P[i] = sum( np.abs(np.multiply(slices[i,:,:], slices[i,:,:].conjugate())) )

        t = np.linspace(g.t[0], g.t[-1], len(P))
 
        w_l_m  =  g('xlamds')
        w_l_ev = h * c / g('xlamds')
        
        #x = np.roll(g.freq_ev, len(g.freq_ev)/2)+ w_l_ev
        spec = fft.fft(slices[:,int( g('ncar')/2),int( g('ncar')/2)])

        y = np.abs(spec)**2
        x = h * fftfreq(len(spec), d=g('zsep') * g('xlamds') / c) + w_l_ev
    

    if power_av == None:
        power_av = P / len(runs)
    else:
        power_av += P / len(runs)

    p1, = ax1.plot(t, P, color='black',alpha=0.4)

    if spec_av == None:
        spec_av = y / len(runs)
    else:
        spec_av += y / len(runs)

    p2, = ax2.plot(x, y, color='black', alpha = 0.4)

    ax2.set_xlim(x[0],x[-1])
    ax3.set_xlim(x[0],x[-1])

    x_ticks = ax2.get_xticks()[1:]

    x2 = h*c/(x_ticks) * 1.e9        # coordinates in nm

    ax3.set_xticks(x_ticks)
    ax3.set_xticklabels(["%.4f" % z for z in x2])

ax1.plot(t, power_av, 'b')
ax2.plot(x, spec_av, 'b')

plt.show()

