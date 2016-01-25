'''
interface to genesis
'''

import struct
from copy import copy

from ocelot.rad.fel import *
from ocelot.cpbd.beam import Beam, gauss_from_twiss
import ocelot.utils.reswake as w
from ocelot.common.math_op import *

h = 4.135667516e-15
c = 299792458.0

inputTemplate = "\
 $newrun \n\
 aw0   =  __AW0__ \n\
 xkx   =  0.000000E+00\n\
 xky   =  1.000000E+00\n\
 wcoefz =  0.000000E+00   0.000000E+00   0.000000E+00\n\
 xlamd =  __XLAMD__\n\
 fbess0 =  __FBESS0__\n\
 delaw =  0.000000E+00\n\
 iertyp =    0\n\
 iwityp =    0\n\
 awd   =  __AW0__ \n\
 awx   =  0.000000E+00\n\
 awy   =  0.000000E+00\n\
 iseed =   __ISEED__\n\
 npart = __NPART__\n\
 gamma0 =  __GAMMA0__\n\
 delgam =  __DELGAM__\n\
 rxbeam =  __RXBEAM__\n\
 rybeam =  __RYBEAM__\n\
 alphax = __ALPHAX__\n\
 alphay = __ALPHAY__\n\
 emitx =  __EMITX__\n\
 emity =  __EMITY__\n\
 xbeam =  __XBEAM__\n\
 ybeam =  __YBEAM__\n\
 pxbeam =  __PXBEAM__\n\
 pybeam =  __PYBEAM__\n\
 conditx =  0.000000E+00\n\
 condity =  0.000000E+00\n\
 bunch =  0.000000E+00\n\
 bunchphase =  0.000000E+00\n\
 emod =  0.000000E+00\n\
 emodphase =  0.000000E+00\n\
 xlamds =  __XLAMDS__\n\
 prad0 =  __PRAD0__\n\
 zrayl =  __ZRAYL__\n\
 zwaist =  0.000000E+00\n\
 ncar  =  __NCAR__\n\
 lbc   =    0\n\
 rmax0 =  __RMAX0__\n\
 dgrid =  __DGRID__\n\
 nscr  =    0\n\
 nscz  =    0\n\
 nptr  =   __NPTR__\n\
 nwig  =   __NWIG__\n\
 zsep  =   __ZSEP__\n\
 delz  =   __DELZ__\n\
 nsec  =   __NSEC__\n\
 iorb  =   __IORB__\n\
 zstop =   __ZSTOP__\n\
 magin =   __MAGIN__\n\
 magout =   __MAGOUT__\n\
 quadf =   __QUADF__\n\
 quadd =   __QUADD__\n\
 fl    =  __FL__\n\
 dl    =  __DL__\n\
 drl   =  __DRL__\n\
 f1st  =  __F1ST__\n\
 qfdx  =  __QFDX__\n\
 qfdy  =  __QFDY__\n\
 solen =  __SOLEN__\n\
 sl    =  __SL__\n\
 ildgam =    5\n\
 ildpsi =    7\n\
 ildx  =    1\n\
 ildy  =    2\n\
 ildpx =    3\n\
 ildpy =    4\n\
 itgaus =    1\n\
 nbins =    4\n\
 igamgaus =    1\n\
 lout  = 1 1 1 1 1 0 1 1 1 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0\n\
 iphsty =    2\n\
 ishsty =    1\n\
 ippart =    0\n\
 ispart =    0\n\
 ipradi =    0\n\
 isradi =    0\n\
 idump =    0\n\
 iotail =    __IOTAIL__\n\
 nharm =    __NHARM__\n\
 iallharm =    1\n\
 iharmsc =    0\n\
 curpeak =  __CURPEAK__\n\
 curlen =  __CURLEN__\n\
 ntail = __NTAIL__\n\
 nslice = __NSLICE__\n\
 __SHOTNOISE__\n\
 isntyp =    0\n\
 iall  =    0\n\
 __ITDP__\n\
 ipseed =   __IPSEED__\n\
 iscan =    0\n\
 nscan =    3\n\
 svar  =  __SVAR__\n\
 isravg =    __SRAVG__\n\
 isrsig =    __SRSIG__\n\
 cuttail = -1.000000E+00\n\
 eloss =  __ELOSS__\n\
 version =  1.0\n\
 ndcut =   -1\n\
 idmpfld =    __DUMP_FIELDS__\n\
 idmppar =    __DUMP_PARTICLES__\n\
 ilog  =    0\n\
 ffspec =    0\n\
 convharm =    1\n\
 ibfield =  0.000000E+00\n\
 imagl =    0.000000E+00\n\
 idril =    0.000000E+00\n\
 alignradf =    1\n\
 offsetradf =    0\n\
 multconv =    0\n\
__BEAMFILE__\n\
__FIELDFILE__\n\
__MAGFILE__\n\
 outputfile ='run.__RUNID__.gout'\n\
 filetype ='ORIGINAL'\n\
 $end\n"

class GenesisInput:
    
    def __init__(self):
        
        # defaults
                
        self.runid = 0
                
        self.iseed = -1  #initial seeding of the random number generator for field errors
        self.ipseed = -1  #initial seeding of the random number generator for shot noise     
        self.emitx = 1.0e-7
        self.emity = 1.0e-7       

        self.alphax = 0.0
        self.alphay = 0.0       
        
        self.gamma0 =  3.424658E+04            

        self.curpeak = 2.500000E+03
        self.curlen = 7E-06
        self.zsep = 20   #separation between slices in terms of radiation length xlamds        
        self.zrayl = 0.00001
        
         
        self.npart = 8192   # number of macroparticles per slice
        self.ncar = 151   # number of grid points for field calculation along one axis
        self.nslice = 1504   # number of slices
        self.delz = 1.0   # time step in terms of undulator periods

        self.ntail = - self.nslice / 2

        self.rmax0 = 9.0   # mesh size in units of radiation+beam sizes
        self.dgrid = 0.0   # exmplicit mesh size
        self.nptr = 40   # space charge mesh points
        
        self.nwig = 98   # 
        self.nsec = 1   #
        self.quadf = 0   #
        self.quadd = 0   # 
        self.fl = 0   #
        self.dl = 0   #
        self.drl = 0   # 
        self.qfdx = 0   #
        self.qfdy = 0   # 
        self.solen = 0   #
        self.sl = 0   #
        self.f1st = 0   # 

        self.eloss = 0
        self.srsig = 1 # energy fluctuations from sr
        self.sravg = 1 # energy loss from sr

        self.iorb = 0   # enforce orbit correction
        
        self.magin = 1   # read in magnetic lattice
        self.magout = 0   # output magnetic lattice
        
        self.zstop = 256.0
        
        self.svar = 0.01  
    
        self.nharm = 1
        self.iotail = 1
    
        self.type = 'steady'
        self.DUMP_FIELDS = 0
        self.DUMP_PARTICLES = 0

        #self.useBeamFile = False
        self.beamfile = None
        self.fieldfile = None

    def input(self):
        input = inputTemplate
        
        if self.type == 'steady':
            input = input.replace("__SHOTNOISE__", "itdp  =    0")
            input = input.replace("__ITDP__", "itdp = 0") 
        else:
            input = input.replace("__SHOTNOISE__", "shotnoise=  1.000000E+00")
            input = input.replace("__ITDP__", "itdp = 1")
            self.prad0 = 0   
            
        if self.beamfile != None:
            input = input.replace("__BEAMFILE__", " beamfile  =  '"+ str(self.beamfile)+ "'")
        else:
            input = input.replace("__BEAMFILE__\n", "")

        if self.fieldfile != None:
            input = input.replace("__FIELDFILE__", " fieldfile  =  '"+ str(self.fieldfile)+ "'")
        else:
            input = input.replace("__FIELDFILE__\n", "")

        
        if self.magin == 0:
            input = input.replace("__MAGFILE__\n", "")
        else:
            input = input.replace("__MAGFILE__\n", "maginfile ='lattice.inp'\n")
        
        for p in self.__dict__.keys():
            input = input.replace("__"  + str(p).upper() + "__", str(self.__dict__[p]))
                                
        return input
    
    def __getattr__(self, name):
        if name not in self.__dict__.keys():
            return 0
        else:
            return self.__dict__[name]


class GenesisOutput:
    
    def __init__(self):
        self.z = []
        self.I = []
        self.n = []
        self.zSlice = []
        self.E = []
        self.aw = []
        self.qfld = []
        
        self.sliceKeys = []
        self.sliceValues = {}
        
        self.parameters = {}
            
    
    def __call__(self, name):
        '''
        if name not in self.__dict__.keys():
            return 0
        else:
            return self.__dict__[name]
        '''
        if name not in self.parameters.keys():
            return 0.0
        else:
            p, = self.parameters[name]
            return float(p.replace('D','E'))
        
        
        
class GenesisBeamDefinition():
    
    def __init__(self):
        self.columns=[]
        self.column_values={}


''' 
   I/O functions
'''

def read_particle_file(file_name, npart, nslice):
    data=open(file_name,'rb').read()
    dataSize = len(data)
    
    print 'sizes ', dataSize, npart, nslice, dataSize / (6.0*npart*nslice) /8.0  
    
    slices = []
    
    nn = npart * 6
    
    for i in range(0, nslice):

        buf=map(lambda x: struct.unpack('d',data[8*x:8*x+8])[0] , range(i*nn,(i+1)*nn))

        E = np.array(buf[0:npart])
        pz = np.array(buf[npart:npart*2])
        x = np.array(buf[npart*2:npart*3])
        y = np.array(buf[npart*3:npart*4])
        px = np.array(buf[npart*4:npart*5])
        py = np.array(buf[npart*5:npart*6])

        slices.append([E,pz,x,px,y,py])

    return slices

def readParticleFile(fileName, npart, nslice):
    
    def read_in_chunks(f, size=1024):
        while True:
            data = f.read(size)
            if not data:
                break
            yield data

    
    f = open(fileName,'rb')
    f.seek(0,2)
    total_data_len = f.tell() / 8 / 6
    f.seek(0,0)
    slice_size = int(6*npart * 8.0)
    
    data_size = total_data_len * 16 

    print 'slice size = ', slice_size
    print 'data size = ', data_size
    print 'total_data_len = ', total_data_len
    
    
    n_slices = total_data_len / (slice_size/ 6 / 8)
    
    print 'n_slices = ', n_slices

    
    
    slices = np.zeros([n_slices,npart*6], dtype=double)
    
    n = 0
    
    for piece in read_in_chunks(f, size  = slice_size):
        
        if True :
            #print 'slice %s' % n
            if ( len(piece) / 8 != npart * 6):
                print 'warning, wrong slice size'
        
            for i in xrange(len(piece) / 8 ):
                slices[n,i] = struct.unpack('d',piece[8*i:8*i+8])[0] 
            n += 1
    
    print 'read slices: ', slices.shape

    return slices


def read_beam_file(fileName):
    beam = GenesisBeamDefinition()
    
    f=open(fileName,'r')
    f.readline()
    for line in f: 
        tokens = line.strip().split()
        
        if len(tokens) < 2:
            continue
        
        #print tokens[0:2]
        
        if tokens[0] == "?" and tokens[1] == "COLUMNS":
            beam.columns = tokens[2:]
            for c in beam.columns:
                beam.column_values[c] = []
                
            print beam.columns
 
        if tokens[0] != "?":
            #print tokens
            for i in range(0,len(tokens)):
                beam.column_values[beam.columns[i]].append( float (tokens[i]) )
            
    #print beam.columns

    beam.z = beam.column_values['ZPOS']
    beam.zsep = beam.z[1] - beam.z[0]
    beam.I = np.array(beam.column_values['CURPEAK'])
    try:
        beam.ex = beam.column_values['EMITX']
        beam.ey = beam.column_values['EMITY']
        beam.betax = beam.column_values['BETAX']
        beam.betay = beam.column_values['BETAY']
        
        beam.alphax = beam.column_values['ALPHAX']
        beam.alphay = beam.column_values['ALPHAY']
        
        beam.x = beam.column_values['XBEAM']
        beam.y = beam.column_values['YBEAM']
        beam.px = beam.column_values['PXBEAM']
        beam.py = beam.column_values['PYBEAM']
        beam.g0 = np.array(beam.column_values['GAMMA0'])
        beam.dg = np.array(beam.column_values['DELGAM'])
    except:
        pass
    
    try:
        beam.eloss = np.array(beam.column_values['ELOSS'])
    except:
        beam.eloss = np.zeros_like(beam.I)
    
    return beam


def readRadiationFile(fileName='simulation.gout.dfl', npoints=51, slice_start=0, slice_end = -1, idx=None):
    
    def read_in_chunks(file, size=1024):
        while True:
            data = file.read(size)
            if not data:
                break
            yield data

    
    f = open(fileName,'rb')
    f.seek(0,2)
    total_data_len = f.tell() / 8 / 2
    f.seek(0,0)
    slice_size = int(2.0*npoints*npoints * 8.0)
    
    if slice_start > 0:
        f.seek(slice_size*slice_start,0)

    if slice_end > 0:
        data_size = (slice_end - slice_start) * slice_size
    else:
        data_size = total_data_len * 16 - (slice_start) * slice_size

    print 'slice size = ', slice_size
    print 'data size = ', data_size
    print 'total_data_len = ', total_data_len
    
    ncar = int(np.sqrt((slice_size/16)))
    print 'ncar=', ncar
    
    n_slices = total_data_len / (slice_size/16)
    if idx == None:
        slices = np.zeros([n_slices,ncar,ncar], dtype=complex)
    else:
        slices = np.zeros([len(idx),ncar,ncar], dtype=complex)
    n = 0
    id = 0
    for piece in read_in_chunks(f, size  = slice_size):
        
        if (idx == None) or (n in idx) :
            #print 'reading', n
            if ( len(piece) / 16 != ncar**2):
                print 'warning, wrong slice size'
        
            for i in xrange(len(piece) / 16 ):
                i2 = i % ncar
                i1 = int(i / ncar)
                #print struct.unpack('d',piece[16*i:16*i+8])
                slices[id,i1,i2] = struct.unpack('d',piece[16*i:16*i+8])[0] + 1j*struct.unpack('d',piece[16*i+8:16*(i+1)])[0] 
            id += 1
        n+= 1
    
    print 'read slices: ', slices.shape

    return slices


def readRadiationFile_mpi(comm=None, fileName='simulation.gout.dfl', npoints=51):
    '''
    not advisable to be used with very small n_proc due to memory overhead ~ file_size / n_proc  
    '''
    from mpi4py import MPI
    
    def read_in_chunks(file, size=1024):
        while True:
            data = file.read(size)
            if not data:
                break
            yield data

    
    rank = comm.Get_rank()
    nproc = comm.Get_size()
        
    f = open(fileName,'rb')
    f.seek(0,2)
    total_data_len = f.tell() / 8 / 2
    f.seek(0,0)
    slice_size = int(2.0*npoints*npoints * 8.0)
    n_slices = int(total_data_len / (npoints**2))

    ncar = int(np.sqrt((slice_size/16)))
    
    local_data_len = int(n_slices / nproc)
    
    n_extra = n_slices - local_data_len * nproc
    
    tmp_buf = np.zeros([local_data_len,ncar,ncar], dtype=complex)

    
    if rank == 0:
        slice_start  = rank * local_data_len
        slices = np.zeros([n_slices,ncar,ncar], dtype=complex)
        slices_to_read = local_data_len + n_extra
    else:
        slice_start = rank * local_data_len + n_extra
        slices = []
        slices_to_read = local_data_len
        
    n = 0
    
    f.seek(slice_start*slice_size, 0)

    print 'rank', rank, ' reading', slice_start, slices_to_read, n_extra

    for piece in read_in_chunks(f, size  = slice_size):
                
        if n >= slices_to_read :
            break
        
        if ( len(piece) / 16 != ncar**2):
            print 'warning, wrong slice size'
    
        for i in xrange(len(piece) / 16 ):
            i2 = i % ncar
            i1 = int(i / ncar)
            if rank == 0:
                #print n, n_extra
                v = struct.unpack('d',piece[16*i:16*i+8])[0] + 1j*struct.unpack('d',piece[16*i+8:16*(i+1)])[0]
                slices[n,i1,i2] = v
                if n - n_extra >= 0:
                    tmp_buf[n-n_extra,i1,i2] = v
            else:
                tmp_buf[n,i1,i2] = struct.unpack('d',piece[16*i:16*i+8])[0] + 1j*struct.unpack('d',piece[16*i+8:16*(i+1)])[0] 
        n += 1
    #print rank, 'tmp_buf=', tmp_buf
    comm.Gather([tmp_buf,  MPI.COMPLEX], [slices[n_extra:], MPI.COMPLEX])
    
    return slices

def writeRadiationFile(filename, slices):
    f=open(filename,'wb')  
    n1, n2 = slices.shape[1], slices.shape[2]
    for i1 in xrange(slices.shape[0]):
        str_bin = ''
        for i2 in xrange(n1):
            for i3 in xrange(n2):
                str_bin += struct.pack('d',slices[i1,i2,i3].real) 
                str_bin += struct.pack('d',slices[i1,i2,i3].imag)
        f.write(str_bin)

    f.close()


def writeRadiationFile_mpi(comm, filename, slices, shape):
    '''
    rank 0 should contain the slices
    '''
    from mpi4py import MPI 
    n_slices, n1, n2 = shape[0], shape[1], shape[2]
       
    rank = comm.Get_rank()
    nproc = comm.Get_size()
    
    if nproc == 1:
        f=open(filename,'wb')  
    else:
        f=open(filename + '.' + str(rank),'wb')
    
    slice_size = n1*n2*8*2
    
    local_data_len = int(n_slices / nproc)
    
    n_extra = n_slices - local_data_len * nproc

    tmp_buf = np.zeros([local_data_len,n1,n2], dtype=complex)    
       
    if rank == 0:
        slice_start  = rank * local_data_len
        slices_to_write = local_data_len + n_extra
    else:
        slice_start = rank * local_data_len + n_extra
        slices = []
        slices_to_write = local_data_len
    
    #print slices
       
    comm.Scatter([slices[n_extra:],  MPI.COMPLEX], [tmp_buf, MPI.COMPLEX])
        
    #print slices
    #print tmp_buf
        
    #f.seek(slice_start*slice_size, 0)
        
    print 'rank=', rank, 'slices_to_write=', slice_start, slices_to_write
        
    
    for i1 in xrange(slices_to_write):
        str_bin = ''
        for i2 in xrange(n1):
            for i3 in xrange(n2):
                
                if rank > 0:
                    #print '>0', tmp_buf[i1,i2,i3]
                    str_bin += struct.pack('d',tmp_buf[i1,i2,i3].real) 
                    str_bin += struct.pack('d',tmp_buf[i1,i2,i3].imag)
                else:
                    #print '0', slices[i1,i2,i3]
                    str_bin += struct.pack('d',slices[i1,i2,i3].real) 
                    str_bin += struct.pack('d',slices[i1,i2,i3].imag)
                    
        #print 'writing', str_bin
        f.write(str_bin)
    
    f.close()
    
    if rank == 0 and nproc>1:
        print 'merging temporary files'
        cmd = 'cat '
        cmd2 = 'rm '
        for i in xrange(nproc): 
            cmd += ' ' + str(filename) + '.' + str(i)
            cmd2 += ' ' + str(filename) + '.' + str(i)
        cmd = cmd + ' > ' + filename
        cmd = cmd + ' ; ' + cmd2
        print cmd
        os.system(cmd)


def readGenesisOutput(fileName):
    out = GenesisOutput()
    out.path = fileName

    chunk = 'header'
    
    nSlice = 0
    
    f=open(fileName,'r')
    f.readline()
    for line in f: 
        tokens = line.strip().split()

        if len(tokens) < 1:
            #chunk = 'none'
            continue
        
        if tokens == ['z[m]', 'aw', 'qfld']:
            chunk = 'optics'
            print 'reading optics '
            continue
        
        if tokens[0] == 'Input':
            chunk = 'input'
            print 'reading input parameters'
            continue
        
        #********** output: slice    10
        if tokens[0] == '**********':
            #print 'slice:', tokens[3]
            chunk = 'slices'
            nSlice = int(tokens[3])
         
        if tokens[0] == 'power':
            chunk = 'slice'
            out.sliceKeys = copy(tokens)
            #out.sliceValues[nSlice] = copy.copy(tokens)
            out.sliceValues[nSlice] = {}
            for i in range(0,len(tokens)):
                out.sliceValues[nSlice][out.sliceKeys[i]] = []
            #print 'reading slices'
            #print out.sliceKeys
            #print out.sliceValues[nSlice]
            continue
            
        if chunk == 'optics':
            z,aw,qfld = map(float,tokens)
            out.z.append(z)
            out.aw.append(aw)
            out.qfld.append(qfld)
            
        if chunk == 'input':
            tokens=line.replace('=','').strip().split()
            out.parameters[tokens[0]] = tokens[1:]
            #print 'input:', tokens


        if chunk == 'slice':
            vals = map(float,tokens)
            #print vals
            for i in range(0,len(vals)):
                out.sliceValues[nSlice][out.sliceKeys[i]].append(vals[i])
            
            #out.zSlice.append(vals[2])
            #out.aw.append(aw)
            #out.qfld.append(qfld)

        if chunk == 'slices':
            if len(tokens) == 2 and tokens[1]=='current':
                #print tokens[1]
                out.I.append(float(tokens[0]))
                out.n.append(nSlice)



    out.nSlices = len(out.sliceValues ) 
    out.nZ = len(out.sliceValues[1][out.sliceKeys[0]])

    print 'nSlice', out.nSlices
    print 'nZ', out.nZ


    out.power = []
    out.phi = []
    out.power_z = 0*np.array(out.sliceValues[out.sliceValues.keys()[1]]['power'])
    out.power_int = []
    out.max_power = 0.0
    for i in xrange(1,out.nSlices):
        pend = out.sliceValues[out.sliceValues.keys()[i]]['power'][-1]
        out.power_int.append(pend)
        out.power.append(out.sliceValues[out.sliceValues.keys()[i]]['p_mid'][-1])
        out.phi.append(out.sliceValues[out.sliceValues.keys()[i]]['phi_mid'][-1])
        out.power_z +=  np.array(out.sliceValues[out.sliceValues.keys()[i]]['power']) / out.nSlices

        if out.max_power < pend: out.max_power = pend

    out.t = 1.0e+15 * out('zsep') * out('xlamds') / c * np.arange(0,len(out.power))
    out.dt = (out.t[1] - out.t[0]) * 1.e-15
    
    out.spec = fft.fft(np.sqrt( np.array(out.power) ) * np.exp( 1.j* np.array(out.phi) ) )
    out.freq_ev = h * fftfreq(len(out.spec), d=out('zsep') * out('xlamds') / c) 

    out.power = np.array(out.power)
    out.phi = np.array(out.phi)
    out.power_int = np.array(out.power_int)

    out.z = np.array(out.z)
    out.I = np.array(out.I)

    return out

def getAverageUndulatorParameter(lattice, unit=1.0, energy = 17.5):
    positions = sorted(lattice.lattice.keys())
          
    prevPos = 0

    ks = []
    ls = []

    for pos in positions:
        if lattice.elements[lattice.lattice[pos].id].type == 'undulator':
            e = lattice.elements[lattice.lattice[pos].id]
            l = float(e.params['nperiods']) * float(e.params['lperiod'])
            
            ks.append( float(e.params['K']) )
            ls.append( l / unit )
            
            #lat += 'AW' +'    '+ e.params['K'] + '   ' + str( l  / unit ) + '  ' + str( (pos - prevPos) / unit ) + '\n'
            
            #if prevPos>0:
            #    drifts.append([str( (pos - prevPos ) / unit ), str(prevLen / unit)])
            
            #prevPos = pos + l
            #prevLen = l 

    return np.mean(ks)
    
def generate_input(up, beam, itdp=False):
    '''
    default input parameters
    '''
    inp = GenesisInput()

    beam.gamma_rel = beam.E / (0.511e-3)
    beam.emit_x = beam.emit_xn / beam.gamma_rel
    beam.emit_y = beam.emit_yn / beam.gamma_rel

    inp.magin = 1
    #inp.zstop = 50
    #inp.nsec = 20
        
    inp.xlamd = up.lw
    inp.aw0 = up.K * np.sqrt(0.5)
    inp.delgam = beam.sigma_E / 0.000510998
    inp.gamma0 = beam.E / 0.000510998
    inp.rxbeam = np.sqrt (beam.emit_x * beam.beta_x )
    inp.rybeam = np.sqrt (beam.emit_y * beam.beta_y )
    
    inp.alphax = beam.alpha_x
    inp.alphay = beam.alpha_y
    
    inp.curpeak = beam.I

    inp.xbeam = beam.x
    inp.ybeam = beam.y
    inp.pxbeam = beam.xp
    inp.pybeam = beam.yp

    inp.emitx = beam.emit_xn
    inp.emity = beam.emit_yn

    felParameters = calculateFelParameters(inp)

    inp.xlamds = felParameters.lambda0
    inp.prad0 = felParameters.power
    inp.fbess0 = felParameters.fc
    inp.zrayl = felParameters.zr
    
    
    if itdp:
        inp.type = "tdp"
        inp.DUMP_FIELDS = 1
        inp.ipseed = 132
        inp.ncar = 151
        inp.nslice = 300
        inp.curlen = beam.tpulse * 3.e-7
        inp.zsep = 8 * int(inp.curlen  / inp.nslice / inp.xlamds )
    
    inp.ntail = - int ( inp.nslice / 2 )
    inp.npart = 2048
    inp.rmax0 = 9.0

    inp.delz = 4.0

    # print out FEL parameter estimates
    #printFelParameters(inp)
    return inp


        

def generate_lattice(lattice, unit=1.0, energy = None, debug = False):
    
    print 'generating lattice file...'
    
    lat = '# header is included\n? VERSION= 1.00  including new format\n? UNITLENGTH= '+str(unit)+' :unit length in header\n'
    undLat = ''
    quadLat = ''
    driftLat = ''

    e0 = lattice.sequence[0]
    prevPos = 0
    prevLen = 0
    prevPosQ = 0
    prevLenQ = 0
    
    pos = 0
    
    drifts = []
    quads = []
    
    for e in lattice.sequence:
        
        l = float(e.l)
        
        #print e.type, pos, prevPos
        
        if e.type == 'undulator':

            l = float(e.nperiods) * float(e.lperiod)
            
            undLat += 'AW' +'    '+ str(e.Kx * np.sqrt(0.5)) + '   ' + str( (l  / unit) ) + '  ' + str( ((pos - prevPos - prevLen) / unit) ) + '\n'

            if debug: print 'added und ', 'pos=',pos, 'prevPos=',prevPos,'prevLen=',prevLen 
                        
            if prevLen>0:
                #drifts.append([str( (pos - prevPos ) / unit ), str(prevLen / unit)])
                if debug: print 'appending drift', str( (prevLen ) / unit )
                driftLat += 'AD' +'    '+ str(e.Kx*np.sqrt(0.5)) + '   ' + str( ((pos - prevPos - prevLen) / unit ) ) + '  ' + str( (prevLen  / unit) ) + '\n'

            prevPos = pos
            prevLen = l 

            
        elif e.type == 'rbend' or e.type == 'sbend' or e.type == 'drift':
            pass
        
        elif e.type == 'quadrupole':
            #k = energy/0.2998 * float(e.k1) *  ( e.l / unit - int(e.l / unit) )
            #k = float(energy) * float(e.k1) / e.l #*  (1 +  e.l / unit - int(e.l / unit) )
            #k = float(energy) * float(e.k1) * 0.2998 / e.l #*  (1 +  e.l / unit - int(e.l / unit) )
            k = float(energy) * float(e.k1) / 0.2998
            if debug: print 'DEBUG', e.k1, k, energy
            quadLat += 'QF' +'    '+ str(k) + '   ' + str( (e.l / unit ) ) + '  ' + str( ( (pos - prevPosQ - prevLenQ)  / unit) ) + '\n'
            prevPosQ = pos
            prevLenQ = l 

            #pass
  
  
        pos = pos + l


    return lat + undLat + driftLat + quadLat

def next_run_id(dir = '.'):
    run_ids = []
    for f in os.listdir(dir):
        if f.startswith('run_'):
            run_ids.append( int( f.replace('run_','')) )

    if len(run_ids) > 0:
        return np.max(run_ids) + 1
    return 0


''' 
   standrard post-processing functions
'''


def get_spectrum(power,phase, smax = 1.0):

    ''' pulse spectrum in eV '''
    
    spec = fft.fft(np.sqrt(power) * np.exp(1.j * phase ))
    spec = np.sqrt( spec * np.conj(spec) )
    spec = np.real( np.roll(spec, len(spec)/2) )
    
    c = 299792458.0
    h = 4.135667516e-15
    
    tmax = smax / c
    freq = h * (np.arange(1,len(spec)+1) / tmax - 0.5 * len(spec) / tmax)
    return freq, spec


def get_power_exit(g):

    xlamds = g('xlamds')
    zsep = g('zsep')


    power = np.zeros(len(g.sliceValues.keys()))
    #phase = np.zeros(len(g.sliceValues.keys()))
            
    for i in g.sliceValues.keys():
        power[i-1] = g.sliceValues[i]['power'][-1]
        #phase[i-1] = g.sliceValues[i]['phi_mid'][iZ]

    t = 1.0e+15 * zsep * xlamds / c * np.arange(0,len(power))

    return power, t

def get_power_z(g):
    #nslice = int(g('nslice'))
    nslice = len(g.sliceValues.keys())
    nz = len(g.sliceValues[g.sliceValues.keys()[0]]['power'])
    power_z = np.zeros(nz)
    for i in xrange(nz):
        for j in xrange(nslice):
            power_z[i] += g.sliceValues[g.sliceValues.keys()[j]]['power'][i]

    return power_z / nslice


def beam_file_str(beam):
    #header = "# \n? VERSION = 1.0\n? SIZE ="+str(len(beam.column_values['ZPOS']))+"\n? COLUMNS ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY XBEAM YBEAM PXBEAM PYBEAM ALPHAX ALPHAY CURPEAK ELOSS\n"
    header = "# \n? VERSION = 1.0\n? SIZE ="+str(len(beam.z))+"\n? COLUMNS"
    #ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY XBEAM YBEAM PXBEAM PYBEAM ALPHAX ALPHAY CURPEAK ELOSS\n"
    for c in beam.columns:
        header = header + " " + c 
    header +="\n"
    
    f_str = header
    
    
    beam.column_values['ZPOS'] = beam.z 
    beam.column_values['CURPEAK'] = beam.I
    
    try:
        beam.column_values['EMITX'] = beam.ex
        beam.column_values['EMITY'] = beam.ey
        
        beam.column_values['BETAX'] = beam.betax
        beam.column_values['BETAY'] = beam.betay
        
        beam.column_values['ALPHAX'] = beam.alphax
        beam.column_values['ALPHAY'] = beam.alphay
	
        beam.column_values['XBEAM'] = beam.x
        beam.column_values['YBEAM'] = beam.y
        
        beam.column_values['PXBEAM'] = beam.px
        beam.column_values['PYBEAM'] = beam.py
        
        beam.column_values['GAMMA0'] = beam.g0
        beam.column_values['DELGAM'] = beam.dg
        
        beam.column_values['ELOSS'] = beam.eloss
    except:
        pass
    
    for i in xrange(len(beam.z)):
        for c in beam.columns:
            buf = str(beam.column_values[c][i])
            f_str = f_str + buf + ' '
        f_str = f_str.rstrip() +  '\n'
    
    return f_str


def add_wake_to_beamf(beamf, new_beamf):
    beam = read_beam_file(beamf)
    s, bunch, wake = w.xfel_pipe_wake(s=array(beam.z), current=array(beam.I))
    print 'read ', len(wake), ' slice values'
    beam.eloss = wake[::-1]

    f=open(new_beamf,'w')
    f.write(beam_file_str(beam))
    f.close()




'''
def plot_beam(fig, beam):
    
    ax = fig.add_subplot(321) 
    plt.grid(True)
    ax.set_xlabel(r'$\mu m$')
    p1,= plt.plot(1.e6 * np.array(beam.z),beam.I,'r',lw=3)
    plt.plot(1.e6 * beam.z[beam.idx_max],beam.I[beam.idx_max],'bs')
    
    ax = ax.twinx()
    
    p2,= plt.plot(1.e6 * np.array(beam.z),1.e-3 * np.array(beam.eloss),'g',lw=3)
    
    ax.legend([p1, p2],['I','Wake [KV/m]'])
    
    ax = fig.add_subplot(322) 
    plt.grid(True)
    ax.set_xlabel(r'$\mu m$')
    #p1,= plt.plot(1.e6 * np.array(beam.z),1.e-3 * np.array(beam.eloss),'r',lw=3)
    p1, = plt.plot(1.e6 * np.array(beam.z),beam.g0,'r',lw=3)
    ax = ax.twinx()
    p2, = plt.plot(1.e6 * np.array(beam.z),beam.dg,'g',lw=3)

    ax.legend([p1,p2],[r'$\gamma$',r'$\delta \gamma$'])
    
    ax = fig.add_subplot(323) 
    plt.grid(True)
    ax.set_xlabel(r'$\mu m$')
    p1, = plt.plot(1.e6 * np.array(beam.z),beam.ex, 'r', lw=3)
    p2, = plt.plot(1.e6 * np.array(beam.z),beam.ey, 'g', lw=3)
    plt.plot(1.e6 * beam.z[beam.idx_max],beam.ex[beam.idx_max], 'bs')
    
    ax.legend([p1,p2],[r'$\varepsilon_x$',r'$\varepsilon_y$'])
    #ax3.legend([p3,p4],[r'$\varepsilon_x$',r'$\varepsilon_y$'])
    
    
    ax = fig.add_subplot(324)
    plt.grid(True)
    ax.set_xlabel(r'$\mu m$')
    p1, = plt.plot(1.e6 * np.array(beam.z),beam.betax, 'r', lw=3)
    p2, = plt.plot(1.e6 * np.array(beam.z),beam.betay, 'g', lw=3)
    plt.plot(1.e6 * beam.z[beam.idx_max],beam.betax[beam.idx_max], 'bs')
    
    ax.legend([p1,p2],[r'$\beta_x$',r'$\beta_y$'])


    ax = fig.add_subplot(325)
    plt.grid(True)
    ax.set_xlabel(r'$\mu m$')
    p1, = plt.plot(1.e6 * np.array(beam.z),1.e6 * np.array(beam.x), 'r', lw=3)
    p2, = plt.plot(1.e6 * np.array(beam.z),1.e6 * np.array(beam.y), 'g', lw=3)
    
    ax.legend([p1,p2],[r'$x [\mu m]$',r'$y [\mu m]$'])


    ax = fig.add_subplot(326)
    plt.grid(True)
    ax.set_xlabel(r'$\mu m$')
    p1, = plt.plot(1.e6 * np.array(beam.z),1.e6 * np.array(beam.px), 'r', lw=3)
    p2, = plt.plot(1.e6 * np.array(beam.z),1.e6 * np.array(beam.py), 'g', lw=3)
    
    ax.legend([p1,p2],[r'$p_x [\mu rad]$',r'$p_y [\mu rad]$'])
'''



def transform_beam_file(beam_file = None, out_file='tmp.beam', transform = [ [25.0,0.1], [21.0, -0.1] ], energy_scale=1, energy_new = None, emit_scale = 1, n_interp = None):
    
    beam = read_beam_file(beam_file)
        
    zmax, Imax = peaks(beam.z, beam.I, n=1)
    idx = beam.z.index(zmax)
    beam.idx_max = idx
    print 'matching to slice', idx
    
    #if plot: plot_beam(plt.figure(), beam)    
    
    
    if transform:
        print 'transforming'
        g1x = np.matrix([[beam.betax[idx], beam.alphax[idx]],
                   [beam.alphax[idx], (1+beam.alphax[idx]**2)/beam.betax[idx]]])

        g1y = np.matrix([[beam.betay[idx], beam.alphay[idx]],
                   [beam.alphay[idx], (1+beam.alphay[idx]**2)/beam.betay[idx]]])


        b2x = transform[0][0]
        a2x = transform[0][1]

        b2y = transform[1][0]
        a2y = transform[1][1]
        
        g2x = np.matrix([[b2x, a2x],
                   [a2x, (1+a2x**2)/b2x]])

        g2y = np.matrix([[b2y, a2y],
                   [a2y, (1+a2y**2)/b2y]])


        Mix, Mx = find_transform(g1x,g2x)
        Miy, My = find_transform(g1y,g2y)
        
        #print Mi
        
        betax_new = []
        alphax_new = []

        betay_new = []
        alphay_new = []

        x_new = []
        px_new = []

        y_new = []
        py_new = []

        
        for i in xrange(len(beam.z)):
            g1x = np.matrix([[beam.betax[i], beam.alphax[i]],
                   [beam.alphax[i], (1+beam.alphax[i]**2)/beam.betax[i]]])
    
            gx = Mix.T * g1x * Mix

            g1y = np.matrix([[beam.betay[i], beam.alphay[i]],
                   [beam.alphay[i], (1+beam.alphay[i]**2)/beam.betay[i]]])
    
            gy = Miy.T * g1y * Miy

            
            #print i, gx[0,1], g1x[0,1]
            
            betax_new.append(gx[0,0])
            alphax_new.append(gx[0,1])

            betay_new.append(gy[0,0])
            alphay_new.append(gy[0,1])
            
            '''
            zx = np.matrix([beam.x[i], beam.px[i]])
            zx = Mix * zx
            x_new.appned(zx[0]) 
            px_new.appned(zx[1])

            zy = np.matrix([beam.y[i], beam.py[i]])
            zy = Miy * zy
            y_new.appned(zy[0]) 
            py_new.appned(zy[1])
            '''

            
        #print betax_new
        beam_new = Beam()
        beam_new.column_values = beam.column_values
        beam_new.columns = beam.columns
        
        if energy_new != None:
            gamma_new=energy_new / (0.511e-3)
            energy_scale=gamma_new/np.mean(np.array(beam.g0))
        
        if n_interp == None:
        
            beam_new.idx_max = idx
            beam_new.ex = beam.ex * emit_scale
            beam_new.ey = beam.ey * emit_scale
            beam_new.zsep = beam.zsep
            beam_new.z = beam.z
            beam_new.I = beam.I
            beam_new.g0 = np.array(beam.g0) * energy_scale
            beam_new.dg = np.array(beam.dg)
            
            beam_new.eloss = beam.eloss
    
            beam_new.betax = betax_new
            beam_new.betay = betay_new
            beam_new.alphax = alphax_new
            beam_new.alphay = alphay_new
    
            beam_new.x = np.array(beam.x) * 0
            beam_new.px = np.array(beam.px) * 0
            beam_new.y = np.array(beam.y) * 0
            beam_new.py = np.array(beam.py) * 0
        else:
            
            
            beam_new.z = np.linspace(beam.z[0], beam.z[-1], n_interp) 
            beam_new.I = np.interp(beam_new.z, beam.z, beam.I)
            
            zmax, Imax = peaks(beam_new.z, beam_new.I, n=1)
            beam_new.idx_max = np.where(beam_new.z == zmax)[0][0]
            
            beam_new.ex = np.interp(beam_new.z, beam.z, beam.ex) * emit_scale
            beam_new.ey = np.interp(beam_new.z, beam.z, beam.ey) * emit_scale
            beam_new.zsep = beam.zsep * len(beam.z) / len(beam_new.z)
            beam_new.g0 = np.interp(beam_new.z, beam.z, beam.g0) * energy_scale
            beam_new.dg = np.interp(beam_new.z, beam.z, beam.dg)
            
            beam_new.eloss = np.interp(beam_new.z, beam.z, beam.eloss)
    
            beam_new.betax = np.interp(beam_new.z, beam.z, betax_new)
            beam_new.betay = np.interp(beam_new.z, beam.z, betay_new)
            beam_new.alphax = np.interp(beam_new.z, beam.z, alphax_new)
            beam_new.alphay = np.interp(beam_new.z, beam.z, alphay_new)
    
            beam_new.x = np.interp(beam_new.z, beam.z, beam.x) 
            beam_new.px = np.interp(beam_new.z, beam.z, beam.px) 
            beam_new.y = np.interp(beam_new.z, beam.z, beam.y)
            beam_new.py = np.interp(beam_new.z, beam.z, beam.py)


        
        #if plot: 
        #    plot_beam(plt.figure(), beam_new)    
        #    plt.show()
    
    if transform != None:
        beam_new.f_str = beam_file_str(beam_new)
    
    return beam_new
    
def find_transform(g1,g2):
    '''
    find transform from twiss matrix g1 to g2: x -> M x, g -> Mi.T g M
    '''
    l1, u1 = np.linalg.eig(g1)
    l2, u2 = np.linalg.eig(g2)

    M1 = np.matrix([[u1[0,0], u1[1,0]],
                    [-u1[1,0], u1[0,0]]])
    
    
    d = sqrt(l1[0]/l2[0]) 
    
    M2 = np.matrix([[d, 0],
                    [0, 1.0 / d]])

    M2d = np.matrix([[1./d, 0],
                    [0, d]])


    M3 = np.matrix([[u2[0,0], -u2[1,0]],
                    [u2[1,0], u2[0,0]]])

    return np.linalg.inv(M3*M2*M1), M3*M2d*M1

    
    
def test_beam_transform(beta1=10.0, alpha1=-0.1, beta2=20, alpha2=2.2):
    
    ex = 1.0
    
    g1 = np.matrix([[beta1, alpha1],
                   [alpha1, (1+alpha1**2)/beta1]])

    g2 = np.matrix([[beta2, alpha2],
                   [alpha2, (1+alpha2**2)/beta2]])

    
    Mi, M = find_transform(g1,g2)
    
    
    g = Mi.T * g1 * Mi
    
    print 'g1=', g1
    print 'g2=', g2
    print 'g=', g

    
    x = []
    xp = []

    x2 = []
    xp2 = []

    x3 = []
    xp3 = []

    
    for i in xrange(5000):
        x_, xp_ = gaussFromTwiss(ex, beta1, alpha1)
        x.append(x_)
        xp.append(xp_)
        
        u = M * np.matrix([[x_,xp_]]).T
        
        x3.append(u[0,0])
        xp3.append(u[1,0])
        
    
        x_, xp_ = gaussFromTwiss(ex, beta2, alpha2)
        x2.append(x_)
        xp2.append(xp_)
        
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    
    plt.plot(x,xp,'r.', alpha =0.2)
    plt.plot(x2,xp2,'g.', alpha=0.2)
    plt.plot(x3,xp3,'b.', alpha=0.2)
    plt.grid()
    
    l1, u1 = np.linalg.eig(g1)
    l2, u2 = np.linalg.eig(g2)


    plt.plot([0,u1[0,0]*np.sqrt(l1[0])], [0,u1[1,0]*np.sqrt(l1[0])], 'b--', lw=3)
    plt.plot([0,u1[0,1]*np.sqrt(l1[1])], [0,u1[1,1]*np.sqrt(l1[1])], 'b--', lw=3)

    plt.plot([0,u2[0,0]*np.sqrt(l2[0])], [0,u2[1,0]*np.sqrt(l2[0])], 'b--', lw=3)
    plt.plot([0,u2[0,1]*np.sqrt(l2[1])], [0,u2[1,1]*np.sqrt(l2[1])], 'b--', lw=3)


    z1 = M*u1[:,0] * sqrt(l1[0])
    z2 = M*u1[:,1] * sqrt(l1[1])

        
    plt.plot([0,z1[0]], [0,z1[1]], color='#000000', lw=5, alpha=0.2)
    plt.plot([0,z2[0]], [0,z2[1]], color='#000000', lw=5, alpha=0.2)

            
    plt.show()

    #49.8131287015 1.12127199531 39.9184728466 -0.897874127701
    
        

#import argparse
#parser = argparse.ArgumentParser(description='Data plotting program.')
#parser.add_argument('--type', choices=["trajectory","intensity", "spectrum"], help='sum the integers (default: find the max)')
#parser.add_argument('input', help='input file')
#args = parser.parse_args()

