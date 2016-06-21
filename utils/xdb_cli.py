import sys

import numpy.fft as fft

from ocelot import MagneticLattice
from ocelot.cpbd.optics import *
from ocelot.utils.xdb import Xdb

sys.path.append('../utils/')
from xfel_utils import *

def test_twiss():
    xdb = Xdb(index_file='/home/iagapov/data/xdb/test/index.h5', mode='r')
    f = xdb.read_undulator_config('sase3/1000eV')
    
    exec(f)
    lat = MagneticLattice(sase3)
    
    #rematch(18.0, l_fodo, qdh, lat, extra_fodo, beam, qf, qd) # jeez...
    
    tw0 = Twiss(beam)
    print tw0
    tws=twiss(lat, tw0, nPoints = 1000)
    
    f=plt.figure()
    ax = f.add_subplot(111)
    ax.set_xlim(0, lat.totalLen)
    
    f.canvas.set_window_title('Betas [m]') 
    p1, = plt.plot(map(lambda p: p.s, tws), map(lambda p: p.beta_x, tws), lw=2.0)
    p2, = plt.plot(map(lambda p: p.s, tws), map(lambda p: p.beta_y, tws), lw=2.0)
    plt.grid(True)
    plt.legend([p1,p2], [r'$\beta_x$',r'$\beta_y$', r'$D_x$'])
    plt.show()
    
    
def create_db(index_file):
    #index_file = '/home/iagapov/data/xdb/test/index.h5'
    print 'creating database', index_file 
    xdb = Xdb(index_file=index_file, mode='w')
    xdb.create_index()
        
    xdb.add_undulator('sase1', {'input_file':'/home/iagapov/workspace/xcode/repository/xfel/sase1/sase1.inp'})
        
    xdb.add_undulator_config('sase1', '24000eV/', {'input_file':'/home/iagapov/workspace/xcode/repository/xfel/sase1/sase1.inp'})
        
    xdb.add_undulator('sase3', {'input_file':'/home/iagapov/workspace/xcode/repository/xfel/sase3/sase3.inp'})
    
    xdb.add_undulator_config('sase3', '1000eV/', {'input_file':'/home/iagapov/workspace/xcode/repository/xfel/sase3/sase3.inp'})
    xdb.add_undulator_config('sase3', '1500eV/', {'input_file':'/home/iagapov/workspace/xcode/repository/xfel/sase3/sase3.inp'})
    xdb.add_undulator_config('sase3', '3000eV/', {'input_file':'/home/iagapov/workspace/xcode/repository/xfel/sase3/sase3.inp'})
    
    xdb.add_beam('17.5GeV/nC/80fs/',{})

    xdb.file.close()      

    print 'done'
    
    
def copy_group(in_file, out_file, path, in_group = 'sase', out_group='/FEL'):
    print in_file, out_file, path
    
    f = h5py.File(out_file, 'r+')
    new_group = out_group + '/' + path
    print 'creating group', new_group    
    outfile_root = f.create_group(new_group)

    f.close()

    
    outfile_id  = h5py.h5f.open(out_file,  h5py.h5f.ACC_RDWR)
    infile_id = h5py.h5f.open(in_file, h5py.h5f.ACC_RDONLY)


    infile_root  = h5py.h5g.open(infile_id,  in_group)
    outfile_root = h5py.h5g.open(outfile_id, new_group)

    #outfile_root = h5py.h5g.create(outfile_root, 'iii/zzz/uuu/iii')
    
    #outfile_root = h5py.h5g.open(outfile_id, out_group)

    # Function to copy or link an object in the input file
    def copy_object(name):
        obj_id   = h5py.h5o.open(infile_root, name)
        obj_type = h5py.h5i.get_type(obj_id)
        if obj_type == h5py.h5i.DATASET:
            # If its a dataset, make a link to the original file
            print "Copy object : ",name
            h5py.h5o.copy(infile_root, name, outfile_root, name)
            #outfile_root.links.create_external(name, infilename, name)
        elif obj_type == h5py.h5i.GROUP:
            # If its a group, make corresponding group in the output file
            # (can't use h5ocopy because that would copy the contents too)
            print "Make group  : ",name
            h5py.h5g.create(outfile_root, name)
        else:
            # Anything else gets copied
            print "Copy object : ",name
            h5py.h5o.copy(infile_root, name, outfile_root, name)
        return None
    
    # Execute the function for every object in the input file
    h5py.h5o.visit(infile_root, copy_object)

    outfile_id.close()
    
    del infile_id
    del outfile_id

    
def plot_stats(idx='', base='', root='/FEL/'):

    xdb = Xdb(index_file=idx, mode='r')
    
    try:
        fig = figure(idx)
    except:
        fig = figure()
    
    ax = fig.add_subplot(111)
    ax.grid(True)
    ax.set_xlabel('[fs]')
    ax.set_ylabel('[W]')

    p_med = np.array(xdb.file[root + base + '/pulse_med'])
    p_mean = np.array(xdb.file[root + base + '/pulse_mean'])
    p_std = np.array(xdb.file[root + base + '/pulse_std'])
    t = np.array(xdb.file[root + base + '/t'])
    I = xdb.file[root + base + '/I']
    I = np.concatenate([I,np.zeros(len(t)-len(I))])
    p1, = ax.plot(t, p_med, 'g-', lw=2)
    p2, = ax.plot(t, p_mean, 'b-', lw=5)
    fill_between(t, p_mean - p_std, p_mean + p_std, alpha=0.2)
    #p2, __ , __ , = ax.errorbar(t[::10], p_mean[::10], yerr=p_std[::10], fmt='g--',lw=3, capsize=5)
        
    # gauss fit
    mu, sig = fit_gauss_1d(t, p_mean)
    sig2 = fwhm(t, p_mean)
    print 'pulse duration, fs, (2*rms/ FWHM): ', 2*sig, sig2
    E_pulse = np.sum(p_med) * ( 1.e-15*(t[-1] - t[0]) ) / len(t)
    print 'Pulse energy [mJ]: ', E_pulse * 1.e3
    
    print 'Peak power (mean) [GW]: ', np.max(p_mean) / 1.e9
    print 'Peak power (med) [GW]: ', np.max(p_med) / 1.e9
    
    
    try:
        E_gamma = xdb.file[root + base].attrs["e_gamma"]
        Ng = E_pulse / 1.6e-19 / E_gamma  
        print 'N photons: ', Ng / 1.e12 , ' x 10^12'
    except:
        E_gamma = 1.e24
        print 'N photons: photon energy undefined!!!'
         
    sig = sig2 / 2
    #p3, = ax.plot(t, np.max(p_mean) * np.exp(-(t-mu)**2 / (2.*sig**2)), 'r--', lw=4)
    #plot(t, np.max(p_med) * (I /  np.max(I)), 'b--', lw=3)
    
    '''
    figure()
    plot(t, I, 'b--', lw=3)
    print 'integral', np.sum(I) * (t[1]-t[0])
    '''
    
    #p3, = ax.plot(t, np.sum(p_mean) * (t[-1] - t[0]) / len(t) / (np.sqrt(2*pi) * sig) * np.exp(-(t-mu)**2 / (2*sig**2)), 'r--', lw=3)
    
    
    #ax.legend([p1,p2,p3], ['median','mean/std','gauss fit'])
    
    fig = figure()
    ax = fig.add_subplot(111)
    ax.grid(True)
    ax.set_xlabel('[eV]')
    ax.set_ylabel('[A.U.]')

    s_med = np.array(xdb.file[root + base + '/spec_med'])
    s_mean = np.array(xdb.file[root + base + '/spec_mean'])
    s_std = np.array(xdb.file[root + base + '/spec_std'])
    t = np.array(xdb.file[root + base + '/freq_ev'])
    ax.plot(t, s_med, 'g-', lw=2)
    ax.plot(t, s_mean, 'b-', lw=5)
    mu, sig = fit_gauss_1d(t, p_mean)
    #p2, = ax.plot(t, p_mean, '-', lw=3)
    #p3, = ax.plot(t, np.max(s_mean) * np.exp(-(t-mu)**2 / (2.*sig**2)), 'r--', lw=4)
    #ax.errorbar(t[::10], s_mean[::10], yerr=s_std[::10], fmt='g--',lw=2, capsize=5, alpha=0.9)
    fill_between(t, s_mean-s_std, s_mean +s_std, alpha=0.2)

    # gauss fit
    mu, sig = fit_gauss_1d(t, s_mean)
    sig2 = fwhm(t, s_mean)
    print 'pulse bandwidth, eV, (2*rms/ FWHM): ', 2*sig, sig2, '(', 2 * 100 * sig/E_gamma, 100 * sig2/E_gamma, '%)'
    print 'Max brightness (mean) [AU]: ', np.max(s_mean) / 1.e8
    print 'Max brightness (med) [AU]: ', np.max(s_med) / 1.e8
	
    X = E_gamma * 1.e-3
    print '0.1%bw=', X
    I = np.sum(s_mean) * (t[1] - t[0])
    print 'I [x 10^12] =', I / 1.e12
    Y = Ng / I
    print 'Y' , Ng / (I )
    print 'Max brightness (mean) [N_phot / 0.1% bw]: ', np.max(s_mean) * Y * X
    
    sig = sig2 / 2
    #ax.plot(t, np.max(s_mean) * np.exp(-(t-mu)**2 / (2.*sig**2)), 'r--', lw=3)
    #ax.plot(t, np.max(s_mean) * np.exp(-(t-mu)**2 / (2.*sig2**2)), 'r--', lw=3)
    #ax.plot(t, np.sum(s_mean) * (t[-1] - t[0]) / len(t) / (np.sqrt(2.*pi) * sig) * np.exp(-(t-mu)**2 / (2.*sig**2)), 'r--', lw=3)

    
    #fig = figure()
    #ax = fig.add_subplot(323)
    fig = figure()
    ax = fig.add_subplot(111)

    f_med = xdb.file[root + base + '/field_med']
    try:
        dgrid = xdb.file[root + base].attrs['dgrid']
    except:
        print 'dgrid not defined'
        dgrid = 1
    
    ax.imshow(np.abs(f_med), extent = [-1.e3*dgrid/2,1.e3*dgrid/2,-1.e3*dgrid/2,1.e3*dgrid/2], aspect='auto', cmap='YlOrRd')
    ax.set_ylabel('[mm]')
    ax.set_xlabel('[mm]')
    
    #ax = fig.add_subplot(324)
    
    fig = figure()
    ax = fig.add_subplot(111)

    
    ax.set_ylabel('A.U.')
    ax.set_xlabel('[mm]')
    ax.set_yticks([]) 

    P = np.abs(f_med[:,int(f_med.shape[0]/2)])
    x = np.linspace(-1.e3*dgrid/2.0, 1.e3*dgrid/2.0, len(P))
    ax.plot(x, P, lw=3)

    # gauss fit
    mu, sig = fit_gauss_1d(x, P)
    sig2 = fwhm(x,P)
    #sig /= 1.8
    print 'spot size, mu m, (2*rms / FWHM): ', 2*sig / 1.e-6,'/', sig2 / 1.e-6
    #p3, = ax.plot(x, np.sum(P) * (x[-1] - x[0]) / len(x) / (np.sqrt(2*pi) * sig) * np.exp(-(x-mu)**2 / (2*sig**2)), 'r--', lw=3)

    ax.grid(True)

    #ax = fig.add_subplot(325)
    fig = figure()
    ax = fig.add_subplot(111)

    ax.set_ylabel('A.U.')
    ax.set_xlabel('[mrad]')
    ax.set_yticks([]) 

    P = np.abs(f_med[:,int(f_med.shape[0]/2)])
    x = np.linspace(-1.e3*dgrid/2.0, 1.e3*dgrid/2.0, len(P))
    ax.plot(x, P, lw=3)

    # gauss fit
    mu, sig = fit_gauss_1d(x, P)
    sig2 = fwhm(x,P)
    #sig /= 1.8
    print 'divergence, mu rad, (2*rms / FWHM): ', 2*sig / 1.e-6,'/', sig2 / 1.e-6
    #p3, = ax.plot(x, np.sum(P) * (x[-1] - x[0]) / len(x) / (np.sqrt(2*pi) * sig) * np.exp(-(x-mu)**2 / (2*sig**2)), 'r--', lw=3)
    
    
    #ax = fig.add_subplot(326)
    fig = figure()
    ax = fig.add_subplot(111)
    ax.grid(True)


    ax.set_ylabel('Pulse energy [mJ]')
    ax.set_xlabel('[m]')
    #ax.set_yticks([]) 
    
    z = np.array(xdb.file[root + base + '/z'])

    power_z = np.array(xdb.file[root + base + '/power_z_mean'])
    power_z_std = np.array(xdb.file[root + base + '/power_z_std'])

    fact = E_pulse * 1.e3 / power_z[-1] 
    power_z *= fact
    power_z_std *= fact

    ax.plot(z, power_z, lw=4)
    fill_between(z, power_z - power_z_std, power_z + power_z_std, alpha=0.3)
    #p2, __ , __ , = ax.errorbar(z, power_z, yerr=power_z_std, fmt='g--',lw=3, capsize=5)    
    
    '''
    fig = figure(base + 'HISTO')
    total_power = xdb.file[root + base + '/total_power']
    n, bins = histogram( (total_power/ mean(total_power) - 1) * 10 + 1, 50, normed=False)
    bar( bins[1:] - (bins[1]-bins[0])/2.0 , n / 4.0, width = (bins[1]-bins[0]), alpha=0.5)
    '''
    
    def distr(M,E):
        t = np.linspace(E/20,E*5,1000)
        return t, pow(M,M) / gamma(M) * pow(t/E, M-1) / E * exp(-M*t/E)
    
    #t, E = distr(10.1, 1)
    #plot(t, E, '--', lw=3), grid(True)
    
    show()
    
    
def extract(path, ran, h5_file = None):

    h = 4.135667516e-15
    c = 299792458.0
    verbose = False
    
    powers_z = []
    spectra = []
    pulses = []
    total_power = []
          
    slice_files = []
      
    for run_id in ran:
        
        run_dir =  path + '/run_' + str(run_id)
        t1 = time.time()
        output_file = run_dir + '/run.' + str(run_id) + '.gout'
            
        if rank == 0:
            print 'extracting run ', run_id 
            g = readGenesisOutput(output_file)
        else:
            g = None
            
        g = comm.bcast(g, root=0)
    
        
        t2 = time.time()
        if verbose: print 'file reading:', t2-t1, ' sec'
    
        npoints = g('ncar')
        zstop = g('zstop')
        delz = g('delz')
        xlamd = g('xlamd')
        xlamds = g('xlamds')
        nslice = int(g('nslice'))
        dgrid = g('dgrid')
        zsep = g('zsep')
        E_gamma = h * c / xlamds
            
            
        power_z = get_power_z(g)
        pulse, t = get_power_exit(g)
            
        powers_z.append(power_z)
        
        #slices = readRadiationFile(fileName=output_file + '.dfl', npoints=npoints)
        slices = readRadiationFile_mpi(comm=comm, fileName=output_file+'.dfl', npoints=npoints)
        slice_files.append(output_file+'.dfl')
    
    
        if rank == 0:
            P = np.zeros_like(slices[:,0,0])
            E = np.zeros_like(slices[:,0,0])
            for i in xrange(len(P)):
                P[i] = sum( np.multiply( np.abs(slices[i,:,:]), np.abs(slices[i,:,:]) ) )
                
                #s = sum(  np.multiply(slices[i,:,:]slices[i,:,:])
                #P[i] = abs(s*s.conjugate()) #* (1 / npoints)**2
                E[i] = sum( slices[i,:,:])
        
            
            pulses.append(P)
            
            spec = fft(E)
            spec = np.sqrt( spec * np.conj(spec) )
            spec = np.real( np.roll(spec, len(spec)/2) )
            tmax = nslice * zsep * xlamds / c
            freq_ev = fftfreq(len(spec), d=zsep * xlamds / c)
            freq_ev = h*np.roll(freq_ev, len(spec)/2) 
            
            spectra.append( np.real(spec) )
                    
            total_power.append( sum(P) * zsep * xlamds / c  )
        
        del slices
    
    if rank == 0:
        t = 1.0e+15 * zsep * xlamds * np.arange(0,len(pulses[0])) / c
        
        fig = figure()
        fig.add_subplot(111)
        
        pz_mean, pz_std, pz_med, pz_worst, imed, iworst = stats(powers_z)
        plot(g.z, pz_med, 'b', lw=3)
        errorbar(g.z, pz_mean, yerr=pz_std, fmt='r--',lw=1, capsize=7)
        
        
        fig = figure()
        fig.add_subplot(111)
        
        for i in xrange(len(pulses)):
            plot(t, pulses[i], 'b--', lw=1)
        
        p_mean, p_std, p_med, p_worst, imed, iworst = stats(pulses)
        slices = readRadiationFile(fileName=slice_files[imed], npoints=npoints)
        slices_2d = np.zeros(slices[0,:,:].shape, dtype=complex)
        
        #print slices_2d.shape
        for i in xrange(slices.shape[0]):
            slices_2d[:,:] += slices[i,:,:]
        
        plot(t, p_med, 'g', lw=3)
        errorbar(t, p_mean, yerr=p_std, fmt='r--',lw=1, capsize=7)
    
        n_gamma = np.sum(p_mean)
        
        mu1,mu2, sig1, sig2, sig12 =  fit_gauss_2d(np.linspace(-dgrid/2.0,dgrid/2.0, npoints),np.linspace(-dgrid/2.0,dgrid/2.0, npoints), slices_2d)
        
        field_sig_x, field_sig_y = sig1, sig2
        
    
        fig = figure()
        fig.add_subplot(111)
        
        for i in xrange(len(spectra)):
            plot(freq_ev, spectra[i], 'b--', lw=1)
        
        s_mean, s_std, s_med, s_worst, imed, iworst = stats(spectra)
        plot(freq_ev, s_med, 'g', lw=3)
        #errorbar(t, p_mean, yerr=p_std, fmt='r--',lw=1, capsize=7)
    
    
        figure()
    
        n, bins = histogram(total_power / mean(total_power), 40, normed=False)
        bar( bins[1:] - (bins[1]-bins[0])/2.0 , n, width = (bins[1]-bins[0]), alpha=0.5)
    
        if h5_file != None:
            print 'writing to ', h5_file
            xdb = Xdb(index_file=h5_file, mode='w')
            
            params = {}
            
            params['n_gamma'] = n_gamma
            params['power_z_mean'] = pz_mean
            params['power_z_med'] = pz_med
            params['power_z_std'] = pz_std
            params['pulse_mean'] = np.real(p_mean)
            params['pulse_med'] = np.real(p_med)
            params['pulse_std'] = np.real(p_std)
            params['spec_mean'] = np.real(s_mean)
            params['spec_med'] = np.real(s_med)
            params['spec_std'] = np.real(s_std)
            params['total_power'] = np.real(total_power)
            params['field_med'] = slices_2d
            
            params['field_sig_x'] = field_sig_x
            params['field_sig_y'] = field_sig_y
            
            params['t'] = np.real(t)
            params['z'] = np.array(g.z)
            params['I'] = np.array(g.I)
            params['freq_ev'] = np.real(freq_ev)
            #params['version'] = 'v0'
            params['dgrid'] = dgrid
            params['e_gamma'] = E_gamma
            
            
            xdb.add_fel_calculation('sase/', params, root='/')
                        
        show()
    
def test_print_parameters():
    xdb = Xdb(index_file='/home/iagapov/data/xdb/test/index.h5', mode='r')
    # available undulator configs
     
    global g
    g = xdb.file['Undulators']
    for u in g.keys():
        print u
        for u2 in g[u].keys():
            print '  ', u2

    # available beams
    global s
    s = ''
    def printname(name):
        global s, g
        if len( g[name].keys()) > 1: 
            s = s +  '->' + name + '\n'
        else:
            s = s +  '*' + name + '\n'
    g = xdb.file['Beams']
    g.visit(printname)
    print s
    

import argparse
parser = argparse.ArgumentParser(description='XDB client')
parser.add_argument('--create_db', help='create hdf5 index database in file', action='store_true')
parser.add_argument('--extract', help='process simulation data', action='store_true')
parser.add_argument('--plot', help='display data from hdf5/genesis', action='store_true')
parser.add_argument('--submit', help='submit to main index file', action='store_true')
parser.add_argument('--path', help='various meanings', default='sase')
parser.add_argument('--range', help='various meanings')
parser.add_argument('--format', choices=['hdf5','genesis'], help='file format', default='hdf5')
parser.add_argument('--file', help='file')
parser.add_argument('--dest', help='dest')

args = parser.parse_args()


def get_range(r, dir):
    if r == None:
        return 0
    i = next_run_id(dir) - 1
    i1, i2 = map(int, r.split(':'))
    #print i, i1, i2
    if i < i2:i2 = i
    return i1, i2

if args.create_db:
    create_db(args.file)
    sys.exit(0)
    
if args.extract:
    i1, i2 = get_range(args.range, args.path)
    #print 'extracting from {} to {} range {}:{}'.format(args.path, args.file, i1, i2)
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    extract(args.path,range(i1,i2+1), args.file)


if args.plot:
    print args.format
    if args.format == 'genesis':
        print 'not implemented'
    else:
        plot_stats(idx=args.file, base=args.path, root='/')
    

if args.submit:
    print 'submitting', args.file, args.path, args.dest  
    import h5py, os
    for f in os.listdir(args.path):
        if f.endswith('.h5'):
            #try:
            group_name = args.path + '/' + f.replace('_', '/').replace('.h5','')
            copy_group(args.path + '/' + f, args.file, group_name)
            #except:
            #    print 'error submitting', f
    #xdb = Xdb(index_file=args.dest, mode='r+')
    #xdb = Xdb(index_file='/home/iagapov/data/xdb/test/index.h5', mode='r+')
    #xdb.file['ext link'] = h5py.ExternalLink(args.file, "/test/test/")
    #xdb.file.close()

