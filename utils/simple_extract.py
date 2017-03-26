import sys
import ocelot.common.xio as xio
from ocelot.utils.xdb import Xdb

from ocelot.cpbd.elements import Element, Quadrupole, RBend, Drift, Undulator
from ocelot import MagneticLattice
from ocelot.cpbd.beam import Beam
from ocelot.cpbd.optics import *

import numpy.fft as fft

sys.path.append('../utils/')
from xfel_utils import *

def get_field_exit(g):

    xlamds = g('xlamds')
    zsep = g('zsep')

    power = np.zeros(len(g.sliceValues.keys()))
    phase = np.zeros(len(g.sliceValues.keys()))
            
    for i in g.sliceValues.keys():
        power[i-1] = g.sliceValues[i]['power'][-1]
        phase[i-1] = g.sliceValues[i]['phi_mid'][-1]

    t = 1.0e+15 * zsep * xlamds / c * np.arange(0,len(power))

    return power, phase, t

def extract(output_file):

    h = 4.135667516e-15
    c = 299792458.0
    g = readGenesisOutput(output_file)
    npoints = g('ncar')
    zstop = g('zstop')
    delz = g('delz')
    xlamd = g('xlamd')
    xlamds = g('xlamds')
    nslice = int(g('nslice'))
    dgrid = g('dgrid')
    zsep = g('zsep')
    E_gamma = h * c / xlamds
                        
    pulse, phi,  t = get_field_exit(g)
    
    E = np.zeros_like(pulse, dtype=complex)

    for i in xrange(len(E)):
        E[i] = np.sqrt(pulse[i]) * np.exp(1.j*phi[i])
        
    #slices = readRadiationFile(fileName=output_file + '.dfl', npoints=npoints)
    #slices = readRadiationFile_mpi(comm=comm, fileName=output_file+'.dfl', npoints=npoints)
    #slice_files.append(output_file+'.dfl')
    
    return t, pulse, E
    
if __name__ == "__main__":
    print 'simple extraction'

    f_name = '9kev_20fs_seeded.txt'

    runs = xrange(50,120)

    for id in runs:
        path = sys.argv[1] + '/run_'+ str(id) + '/run.' + str(id) + '.gout'
        t, p, E = extract(path)
        f_obj = open(f_name + '.' + str(id), 'w')
        for i in xrange(len(t)):
            f_obj.write("{0}\t{1}\t{2}\n".format(t[i], p[i]/1.e9, E[i]))

    plt.plot(t,p)
    plt.figure()
    plt.plot(t, np.real(E))
    plt.plot(t, np.imag(E))
    plt.figure()
    plt.plot(t, np.abs(E)**2)
    plt.show()

