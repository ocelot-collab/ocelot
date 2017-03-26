'''
functions common to xfel decks
'''


#from ocelot.common.xio import XIO

#from ocelot.adaptors import srwutil as srw


from pylab import *

from ocelot.adaptors.genesis import *
params = {'backend': 'ps', 'axes.labelsize': 18, 'text.fontsize': 16, 'legend.fontsize': 24, 'xtick.labelsize': 32,  'ytick.labelsize': 32, 'text.usetex': True}
rcParams.update(params)
rc('text', usetex=True) # required to have greek fonts on redhat



def show_output(g, show_field = False, output_file = None, show_slice=0):
    h = 4.135667516e-15
    c = 299792458.0
    xrms = np.array(g.sliceValues[g.sliceValues.keys()[show_slice]]['xrms'])
    yrms = np.array(g.sliceValues[g.sliceValues.keys()[show_slice]]['yrms'])
    power = 0*np.array(g.sliceValues[g.sliceValues.keys()[show_slice]]['power'])
    
    nslice = len(g.sliceValues.keys())
    power_s = np.zeros(nslice)
    for i in xrange( nslice ):
        power +=  np.array(g.sliceValues[g.sliceValues.keys()[i]]['power']) / nslice
        pend =  g.sliceValues[g.sliceValues.keys()[i]]['power'][-1]
        power_s[i] = pend
        
    f = plt.figure()
    f.add_subplot(131), plt.plot(g.z, xrms, lw=3), plt.plot(g.z, yrms, lw=3), plt.grid(True)
    f.add_subplot(132), plt.plot(g.z, power, lw=3), plt.grid(True)
    I = np.array(g.I)
    t = 1.0e+15 * float(g('zsep')) * float(g('xlamds')) * np.arange(0,len(I)) / c
    f.add_subplot(133), plt.plot(t,power_s, lw=3), plt.plot(t,np.array(I) * np.max(power_s) / np.max(I), lw=3), plt.grid(True)
    
    npoints = g('ncar')
    zstop = g('zstop')
    delz = g('delz')
    xlamd = g('xlamd')
    xlamds = g('xlamds')
    nslice = g('nslice')
    zsep = g('zsep')
    dgrid = g('dgrid')
    
    smax = nslice * zsep * xlamds   
                 
    print 'wavelength ', xlamds
    
    if show_field:
        #from mpi4py import MPI
        
        #comm = MPI.COMM_WORLD
        #slices = readRadiationFile_mpi(comm=comm, fileName=file+'.dfl', npoints=npoints)
        slices = readRadiationFile(fileName=output_file + '.dfl', npoints=npoints)
        print 'slices:', slices.shape
    
        E = np.zeros_like(slices[0,:,:])
        for i in xrange(slices.shape[0]): E += np.multiply(slices[i,:,:], slices[i,:,:].conjugate())
    
        #Z = E*E.conjugate()
    
        fig = plt.figure()
        fig.add_subplot(111)
        m = plt.imshow(abs(E),cmap='YlOrRd')
        z = abs(slices[100,:,:])
    
        fig = plt.figure()
        P = np.zeros_like(slices[:,0,0])
        for i in xrange(len(P)):
            #s = slices[i,int(npoints/2),int(npoints/2)]
            s = sum( np.abs(np.multiply(slices[i,:,:], slices[i,:,:])) )
            P[i] = abs(s*s.conjugate()) * (dgrid**2 / npoints )**2  
    
        plot(P)
    plt.figure()
    plt.plot(np.linspace(-dgrid/2,dgrid/2,npoints),abs(E[:,int(npoints/2)]), lw=3), plt.grid(True)
            
    plt.show()


def show_output_old(g, show_field = False, output_file = None):

    xrms = np.array(g.sliceValues[g.sliceValues.keys()[0]]['xrms'])
    yrms = np.array(g.sliceValues[g.sliceValues.keys()[0]]['yrms'])
    power = np.array(g.sliceValues[g.sliceValues.keys()[0]]['power'])
    
    nslice = len(g.sliceValues.keys())
    power_s = np.zeros(nslice)
    for i in xrange( nslice ):
        power +=  np.array(g.sliceValues[g.sliceValues.keys()[i]]['power']) / nslice
        pend =  g.sliceValues[g.sliceValues.keys()[i]]['power'][-1]
        power_s[i] = pend
        
    f = plt.figure()
    f.add_subplot(131), plt.plot(g.z, xrms, lw=3), plt.plot(g.z, yrms, lw=3), plt.grid(True)
    f.add_subplot(132), plt.plot(g.z, power, lw=3), plt.grid(True)
    f.add_subplot(133), plt.plot(power_s, lw=3), plt.grid(True)
    
    npoints = g('ncar')
    zstop = g('zstop')
    delz = g('delz')
    xlamd = g('xlamd')
    xlamds = g('xlamds')
    nslice = g('nslice')
    zsep = g('zsep')
    dgrid = g('dgrid')
    
    smax = nslice * zsep * xlamds 
                     
    print 'npoints, nslice', npoints, nslice
    
    
    if show_field:
        #from mpi4py import MPI
        
        #comm = MPI.COMM_WORLD
        #slices = readRadiationFile_mpi(comm=comm, fileName=file+'.dfl', npoints=npoints)
        slices = readRadiationFile(fileName=output_file + '.dfl', npoints=npoints)
        print 'slices:', slices.shape
    
        E = np.zeros_like(slices[0,:,:])
        for i in xrange(slices.shape[0]): E += slices[i]
    
        Z = E*E.conjugate()
    
        fig = plt.figure()
        fig.add_subplot(111)
        m = plt.imshow(abs(E),cmap='YlOrRd')
        z = abs(slices[100,:,:])
    
        fig = plt.figure()
        P = np.zeros_like(slices[:,0,0])
        for i in xrange(len(P)):
            s = slices[i,int(npoints/2),int(npoints/2)]
        #s = sum(slices[i,:,:])
            P[i] = abs(s*s.conjugate()) * (dgrid**2 / npoints )**2  
    
        plot(P)
    
    plt.show()



'''
putting arbitrarily many plots on single figure
'''
def show_plots(displays, fig):
    n1 = (len(displays) -1 )/2 + 1
    n2 = (len(displays) -1) / n1 +1
    #print n1, n2
    fmt = str(n1)+ str(n2)
    print fmt
    
    for i in xrange(len(displays)):
        ax = fig.add_subplot(fmt + str(i+1))
        ax.grid(True)
        for f in displays[i].data:
            x,y = f(x = np.linspace(-10, 10, 100))
            ax.plot(x, y, '.')

    show()

class Display:
    def __init__(self, data=lambda x: (x, 0*x) , xlabel='',ylabel=''):
        self.data = (data,)
        self.xlabel = xlabel
        self.ylabel = ylabel



def plot_beam(fig, beam):
    
    ax = fig.add_subplot(321) 
    plt.grid(True)
    ax.set_xlabel(r'$\mu m$')
    p1,= plt.plot(1.e6 * np.array(beam.z),beam.I,'r',lw=3)
    plt.plot(1.e6 * beam.z[beam.idx_max],beam.I[beam.idx_max],'bs')
    ax.set_xlim([1.e6 * beam.z[0] , 1.e6 * beam.z[-1]])
    ax = ax.twinx()
    
    p2,= plt.plot(1.e6 * np.array(beam.z),1.e-3 * np.array(beam.eloss),'g',lw=3)
    
    ax.legend([p1, p2],['I','Wake [KV/m]'])
    
    ax = fig.add_subplot(322) 
    plt.grid(True)
    ax.set_xlabel(r'$\mu m$')
    #p1,= plt.plot(1.e6 * np.array(beam.z),1.e-3 * np.array(beam.eloss),'r',lw=3)
    p1, = plt.plot(1.e6 * np.array(beam.z),beam.g0,'r',lw=3)
    ax.set_xlim([1.e6 * beam.z[0] , 1.e6 * beam.z[-1]])
    ax = ax.twinx()
    p2, = plt.plot(1.e6 * np.array(beam.z),beam.dg,'g',lw=3)

    ax.legend([p1,p2],[r'$\gamma$',r'$\delta \gamma$'])
    
    ax = fig.add_subplot(323) 
    plt.grid(True)
    ax.set_xlabel(r'$\mu m$')
    p1, = plt.plot(1.e6 * np.array(beam.z),beam.ex, 'r', lw=3)
    p2, = plt.plot(1.e6 * np.array(beam.z),beam.ey, 'g', lw=3)
    plt.plot(1.e6 * beam.z[beam.idx_max],beam.ex[beam.idx_max], 'bs')
    ax.set_xlim([1.e6 * beam.z[0] , 1.e6 * beam.z[-1]])
    ax.legend([p1,p2],[r'$\varepsilon_x$',r'$\varepsilon_y$'])
    #ax3.legend([p3,p4],[r'$\varepsilon_x$',r'$\varepsilon_y$'])
    
    
    ax = fig.add_subplot(324)
    plt.grid(True)
    ax.set_xlabel(r'$\mu m$')
    p1, = plt.plot(1.e6 * np.array(beam.z),beam.betax, 'r', lw=3)
    p2, = plt.plot(1.e6 * np.array(beam.z),beam.betay, 'g', lw=3)
    plt.plot(1.e6 * beam.z[beam.idx_max],beam.betax[beam.idx_max], 'bs')
    ax.set_xlim([1.e6 * beam.z[0] , 1.e6 * beam.z[-1]])
    ax.legend([p1,p2],[r'$\beta_x$',r'$\beta_y$'])


    ax = fig.add_subplot(325)
    plt.grid(True)
    ax.set_xlabel(r'$\mu m$')
    p1, = plt.plot(1.e6 * np.array(beam.z),1.e6 * np.array(beam.x), 'r', lw=3)
    p2, = plt.plot(1.e6 * np.array(beam.z),1.e6 * np.array(beam.y), 'g', lw=3)
    ax.set_xlim([1.e6 * beam.z[0] , 1.e6 * beam.z[-1]])
    ax.legend([p1,p2],[r'$x [\mu m]$',r'$y [\mu m]$'])


    ax = fig.add_subplot(326)
    plt.grid(True)
    ax.set_xlabel(r'$\mu m$')
    p1, = plt.plot(1.e6 * np.array(beam.z),1.e6 * np.array(beam.px), 'r', lw=3)
    p2, = plt.plot(1.e6 * np.array(beam.z),1.e6 * np.array(beam.py), 'g', lw=3)
    ax.set_xlim([1.e6 * beam.z[0] , 1.e6 * beam.z[-1]])
    ax.legend([p1,p2],[r'$p_x [\mu rad]$',r'$p_y [\mu rad]$'])

def plot_beam_2(fig, beam, iplot=0):
    
    ax = fig.add_subplot(111) 
    plt.grid(True)
    
    if iplot == 0:
    
        ax.set_xlabel(r'$\mu m$')
        ax.set_ylabel('A')
        p1,= plt.plot(1.e6 * np.array(beam.z),beam.I,'r',lw=3)
        #plt.plot(1.e6 * beam.z[beam.idx_max],beam.I[beam.idx_max],'bs')
        
        ax = ax.twinx()
        ax.set_ylabel('KV/m')
        
        p2,= plt.plot(1.e6 * np.array(beam.z),1.e-3 * np.array(beam.eloss),'g',lw=3)
        
        ax.legend([p1, p2],['I',r'$E_{wake}$'])

    if iplot == 1:
    
        ax.set_xlabel(r'$\mu m$')
        #p1,= plt.plot(1.e6 * np.array(beam.z),1.e-3 * np.array(beam.eloss),'r',lw=3)
        p1, = plt.plot(1.e6 * np.array(beam.z),beam.g0,'r',lw=3)
        ax = ax.twinx()
        p2, = plt.plot(1.e6 * np.array(beam.z),beam.dg,'g',lw=3)
    
        ax.legend([p1,p2],[r'$\gamma$',r'$\delta \gamma$'])
    
    if iplot == 2:

        ax.set_xlabel(r'$\mu m$')
        ax.set_ylabel(r'$mm \cdot mrad$')
        p1, = plt.plot(1.e6 * np.array(beam.z),1.e6*np.array(beam.ex), 'r', lw=3)
        p2, = plt.plot(1.e6 * np.array(beam.z),1.e6*np.array(beam.ey), 'g', lw=3)
        #plt.plot(1.e6 * beam.z[beam.idx_max],beam.ex[beam.idx_max], 'bs')
        
        ax.legend([p1,p2],[r'$\varepsilon_x$',r'$\varepsilon_y$'])
        #ax3.legend([p3,p4],[r'$\varepsilon_x$',r'$\varepsilon_y$'])
    
    if iplot == 3:
    
        ax.set_xlabel(r'$\mu m$')
        ax.set_ylabel(r'$m$')
        p1, = plt.plot(1.e6 * np.array(beam.z),beam.betax, 'r', lw=3)
        p2, = plt.plot(1.e6 * np.array(beam.z),beam.betay, 'g', lw=3)
        #plt.plot(1.e6 * beam.z[beam.idx_max],beam.betax[beam.idx_max], 'bs')
        
        ax.legend([p1,p2],[r'$\beta_x$',r'$\beta_y$'])
