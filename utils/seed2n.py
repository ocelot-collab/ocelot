'''
seed.py
w
*Author*
Gianluca Geloni

*Description*
The main python program seeds the output field file from genesis.
Filtering is parallelized.
Program has to be executed with mpirun.
Example: mpirun -n <nproc> python seed_par.py <ARGS>
<ARGS> are defined as follows:
    filename       = ARGS[1]    #Input genesis dfl file
    filenameF      = ARGS[2]    #Output genesis dfl file
    filenameM      = ARGS[3]    #Modulus of filter transmittance - input
    filenameP      = ARGS[4]    #Phase of filter transmittance - input
    filenameMR     = ARGS[5]    #Modulus of filter transmittance - rescaled
    filenamePR     = ARGS[6]    #Phase of filter transmittance - rsecaled
    filenamePout   = ARGS[11]   #Output power after seed (summed over transverse coord)
    filenamePphout = ARGS[12]   #Output t-domain phase (summed over transverse coord)
    filenameSout   = ARGS[13]   #Output spectrum after seed (summed over transverse coord)
    filenameSphout = ARGS[14]   #Output spectral phase after seed (summed over transverse coord)
    
    lam0           = float(ARGS[15])   #resonance frequency from Genesis
    MTR            = int(ARGS[16])     #transverse mesh
    mult           = int(ARGS[17])     #rescaling parameter multiplies the number of slices to reach wanted freq. resolution
    ds             = float(ARGS[18])   #step in the time (s) domain
    dk             = float(ARGS[19])   #step in the frequency (k) domain
    SHF            = int(ARGS[20])     #longitudinal shift of the seed wrt electron beam (chicane tuning)
    nslice         = int(ARGS[21])     #number of slices from Genesis
    dr             = float(ARGS[22])   #Transverse shift parameter to calculate spatio-temporal coupling   
Other functions are used by the program, or are use to filter the field with a single processor (64bit or 32bit), or to extract the field without filtering. In this case, one uses seed.py as a library only     
'''


##########################
#                        #
# Imports used libraries #
#                        #
##########################
import sys, os, time

#sys.path.append('/usr/lib64/python2.6/site-packages/openmpi-intel')
sys.path.append('/data/netapp/xfel/gianluca/products/ocelot') 
    
from mpi4py import MPI
#import ocelot.adaptors.genesis as genesis
    
import numpy as np
import scipy as sp
import scipy.fftpack as fft


import copy
import struct

#import ocelot.adaptors.genesis as genesis
#import ocelot.adaptors.genesis as genesis

###############
#             #
#  FUNCTIONS  #
#             #
###############
   
#BEGIN FUNCTION   

def sscale(Npts=10000, ds = 5.10373e-9):
    '''
    Function name: sscale(Nslice, mult, ds)

    Description  :   Generates an array with s coordinates with number of points Nslice*mult and separation ds, centered at zero

    Arguments:

         -Nslice         : slice number in Genesis simlualtions
         -mult           : number of times the k array should exceed Nslice
         -ds             : separation

    Output:

        -ssc             : array of s

    Date revised: 2012.6.13

    '''

    ssc = []
    sc = np.linspace(-Npts/2, Npts/2, Npts)
    ssc = sc * ds
    return ssc

#END FUNCTION


#BEGIN FUNCTION


def kscale(Npts=10000,k0=4.188790204e10, dk = (5*4.60612e-06 * 21)**(-1)):
    '''
    Function name: kscale(Nslice, mult, k0)

    Description  :   Generates an array with k = 2pi/lambda coordinates with center k0 and number of points Nslice*mult and separation dk*mult

    Arguments:

         -Nslice         : slice number in Genesis simlualtions
         -mult           : number of times the k array should exceed Nslice
         -k0             : central value of k
         -dk             : separation

    Output:

        -k               : array of k 

    Date revised: 2012.6.13

    '''
    
   
    k = []
    sc = np.linspace(-Npts/2, Npts/2, Npts)
    k = k0 + dk*sc #2*np.pi*dk*sc
    #print 'SEED2'
    #print 'k0 = ',k0
    #print 'k-scale from seed_2 is from', k[0], ' to ', k[-1]
    #print 'krange = ',dk*len(sc)
    #print dk*mult
    #exit()
    return k

#END FUNCTION

#BEGIN FUNCTION

def readfilter(filename = 'D:\Python\Scripts\Genesis - scripts\ModT_fig.dat'):
    '''
    Function name: readfilter(filename)

    Description  :   reads a filter in two columns and puts it into an array

    Arguments:

         -filename       : ascii data of the filter
        

    Output:

        -columns         : array of 2 columns: [ 1, 2, 3... ], [ f(1), f(2), f(3)... ]  

    Date revised: 2012.6.13
    '''


    f = file(filename, 'r')  # iterate over the lines in the file
    columns = []
    for line in f:
        # split the line into a list of column values
        #columns1 = line.split('\t')
        columns1 = line.split(' ')
        # clean any whitespace off the items
	#print filename
	#print columns1
        columns.append([float(col.strip()) for col in columns1])
    return columns

#END FUNCTION


#BEGIN FUNCTION

def rescalefilter(filename = 'D:\Python\Scripts\Genesis - scripts\ModT_fig.dat',k=[1]):
    '''
    Function name: rescalefilter(filename,k)

    Description  :   rescales the a filter function in k = 2pi/lambda to the specified range in k

    Arguments:

         -filename       : ascii data of the filter
         -k              : array of values of k= 2pi/lambda to which the filter needs to be rescaled. Interpolation is used
        

    Output:

        -spits           : array with the new data [k, Interdata] with k and Interdata (interpolated data) are 1D arrays

    Date revised: 2012.6.13
    
    '''
    
    
    data = readfilter(filename)
    dataX = [data[i][0] for i in range(len(data))]
    dataY = [data[i][1] for i in range(len(data))]
    
    print 'Filter k range (seed.py, line 197)', dataX[0], dataX[-1]
    
    #print k[2]-k[1]
    #print dataX[2]-dataX[1]

    Interdata = np.interp(k, dataX, dataY)
    spits = [k, Interdata]

    return spits

#END FUNCTION


#BEGIN EXPERIMENTAL FUNCTION

def readRadiationFile_mpi(comm=None, fileName='simulation.gout.dfl', npoints=51):
    #
    #not advisable to be used with very small n_proc due to memory overhead ~ file_size / n_proc  
    #
    
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

    for piece in read_in_chunks(f, size  = slice_size):
                
        if n >= slices_to_read :
            break
        
        if ( len(piece) / 16 != ncar**2):
            print 'warning, wrong slice size'
    
        for i in xrange(len(piece) / 16 ):
            i2 = i % ncar
            i1 = int(i / ncar)
            if rank == 0:
                v = struct.unpack('d',piece[16*i:16*i+8])[0] + 1j*struct.unpack('d',piece[16*i+8:16*(i+1)])[0]
                slices[n,i1,i2] = v
                if n - n_extra > 0:
                    tmp_buf[n-n_extra,i1,i2] = v
            else:
                tmp_buf[n,i1,i2] = struct.unpack('d',piece[16*i:16*i+8])[0] + 1j*struct.unpack('d',piece[16*i+8:16*(i+1)])[0] 
        n += 1
        
    comm.Gather([tmp_buf,  MPI.COMPLEX], [slices[n_extra:], MPI.COMPLEX])
    
    return slices

#END EXPERIMENTAL FUNCTION

#BEGIN EXPERIMENTAL FUNCTION

def writeRadiationFile_mpi(comm, filename, slices, shape):
    '''
    rank 0 should contain the slices
    '''
        
    n_slices, n1, n2 = shape[0], shape[1], shape[2]
       
    rank = comm.Get_rank()
    nproc = comm.Get_size()
    
    confirmed=1
    
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
	    checkex = 0
	    while (checkex ==  0):
	        if os.path.isfile(str(filename) + '.' + str(i)) == 1: 
		    file_i = open(str(filename) + '.' + str(i)).read()
		    if len(file_i) > 100: 
		        checkex = 1
		    else:
		        print 'waiting for '+str(filename) + '.' + str(i)
            cmd += ' ' + str(filename) + '.' + str(i)
            cmd2 += ' ' + str(filename) + '.' + str(i)
        cmd = cmd + ' > ' + filename
        cmd = cmd + ' ; ' + cmd2
        print cmd
	'''
	confirmed = 1
	while (confirmed != 0):
            confirmed = os.system(cmd)
	    if (confirmed != 0): 
	        print 'merging not succeeded... trying again in 5s...'
		time.sleep(5)
	'''
	confirmed = os.system(cmd)    
    confirmed     = comm.bcast(confirmed, root = 0) 
    return confirmed    
	

    
#END EXPERIMENTAL FUNCTION

#BEGIN EXPERIMENTAL FUNCTION


def readRadiationFile_my(fileName, npoints=151):
    import numpy as np
    b=np.fromfile(fileName,dtype=complex)
    slice_num=b.shape[0]/npoints/npoints
    c=b.reshape(slice_num,npoints,npoints)
    return c

#END EXPERIMENTAL FUNCTION
def writeRadiationFile_my(filename,rad):
    print '    writing radiation file' 
    #a new backward compatible version ~10x faster
#    print '        - writing to ', filename
    d=rad.flatten()
    d.tofile(filename,format='complex')


#############################################
#                                           #
#  MAIN PROGRAM TO BE EXECUTED WITH MPIRUN  #
#                                           #
#############################################

if __name__ == "__main__":

    
    #print 'seed.py'
    
    comm = MPI.COMM_WORLD   #Creates one common object "communicator" knwoning about processes
    rank = comm.Get_rank()  #Process number
    size = comm.Get_size()  #Total number of processes
       
    ARGS = sys.argv
    
    #Defines variables from ARGS
    
    filename       = ARGS[1]    #Input genesis dfl file
    filenameF      = ARGS[2]    #Output genesis dfl file
    filenameM      = ARGS[3]    #Modulus of filter transmittance - input
    filenameP      = ARGS[4]    #Phase of filter transmittance - input
    filenameMR     = ARGS[5]    #Modulus of filter transmittance - rescaled
    filenamePR     = ARGS[6]    #Phase of filter transmittance - rsecaled
    filenamePout   = ARGS[7]   #Output power after seed (summed over transverse coord)
    filenamePphout = ARGS[8]   #Output t-domain phase (summed over transverse coord)
    filenameSout   = ARGS[9]   #Output spectrum after seed (summed over transverse coord)
    filenameSphout = ARGS[10]   #Output specftral phase after seed (summed over transverse coord)
    
    lam0           = float(ARGS[11])   #resonance frequency from Genesis
    MTR            = int(ARGS[12])     #transverse mesh
    mult           = int(ARGS[13])     #rescaling parameter multiplies the number of slices to reach wanted freq. resolution
    ds             = float(ARGS[14])   #step in the time (s) domain
    dk             = float(ARGS[15])   #step in the frequency (k) domain
    SHF            = int(ARGS[16])     #longitudinal shift of the seed wrt electron beam (chicane tuning)
    nslice         = int(ARGS[17])     #number of slices from Genesis output
    dr             = float(ARGS[18])   #Transverse shift parameter to calculate spatio-temporal coupling
    
    
    k0 = 2*np.pi/lam0
    
   
    
    #slices = readRadiationFile_mpi(comm=comm, fileName=filename, npoints=MTR)
    slices = readRadiationFile_my(fileName=filename, npoints=MTR)
    #slices = slices[(len(slices)-nslice):len(slices)]
    print len(slices)
    print nslice
    #print len(slices[0])
    
    if rank == 0:
        print 'seed'
        Dslice = len(slices)  - nslice #slices in the rad file are different. we must account for this
	print 'start'
	
	print Dslice*ds
	
	
	print (len(slices) - 1805)*ds
	print 'stop'
    
    
    #Dslice   = comm.bcast(Dslice, root = 0) 
    #comm.barrier()
   
    
    if rank == 0:  #rank 0 is used for input-output
    
        print 'Multiplication factor = ',mult     
	#print filename
	#print filenameF
	#print filenameM
	#print filenameP
	#print filenameMR
	#print filenamePR
	#print filenamePout
	#print filenamePphout
	#print filenameSout
	#print filenameSphout
	         
        print 'Parallelized field filtering.'
	print 'Reading radiation file...'
	
	print 'Read: ',len(slices),' slices vs nslice = ',nslice
	
	
	
	nslice = len(slices)
	print 'NSLICE*DS = ',nslice*ds
	
	
	
				          
	###############################################################################
	# work is now a flat list                                                     #
	# work[2*(n*MTR+m)    *nslice+i] accesses the real      part of the field at  # 
        # work[(2*(n*MTR+m)+1)*nslice+i] accesses the imaginary part of the field at  #
	#                                                                             #
	#  slice: i (0...nslice-1)                                                    # 
	#  transverse position n, m (0..MTR-1)                                        #
	###############################################################################	
	
	print '...Finished reading radiation file.'	
        print 'Prepare for filtering...'
	
	work = np.zeros(MTR*MTR*2*nslice)
	prova = np.real(slices[0:nslice-1,34,33] )
	for n in range(MTR):
	    for m in range(MTR):
	        work[2*(n*MTR+m)    *nslice:2*(n*MTR+m)    *nslice+(nslice-1)] = np.real(slices[0:nslice-1,n,m] )
	        work[(2*(n*MTR+m)+1)*nslice:(2*(n*MTR+m)+1)*nslice+(nslice-1)] = np.imag(slices[0:nslice-1,n,m] )
		
	######################	
	# Defines the filter #
	######################
			
        ntot = nslice*mult
	print ntot
	ntot = pow(2,np.ceil(np.log(1.0*ntot)/np.log(2.0))-1)			
        k  = kscale(ntot, k0, dk)  #the scale in frequency (k) is calculated
        ssc  = sscale(ntot, ds)    #the scale in time (s) is calculated		
        print 'Radiation k range (seed.py, line 480)', k[0], k[-1] 
        Fmod    = rescalefilter(filenameM,k)  #Filter mod is rescaled according to the scale in k
        Fpha    = rescalefilter(filenameP,k)  #Filter pha is rescaled according to the scale in k
        print 'interpolated filter k range (seed.py, line 486)', Fmod[0][0], Fmod[0][-1]
	#Prints the rescaled filter to file
        f1 = open(filenameMR, 'w')
        f3 = open(filenamePR, 'w')
        for i in range(len(k)):
            f1.write('%s ' %(2*np.pi/k[i]) + '%s' %Fmod[1][i] +'\n')
            f3.write('%s ' %(2*np.pi/k[i]) + '%s' %Fpha[1][i] +'\n')
        f1.close()
        f3.close()
	
	##################################################
        # Defines useful quantities for the parallel run #
	##################################################
	
	lenchunck = len(work)/MTR   #work is divited in MTR parts	
        each =  int(MTR/size)+1	    #each processor needs to process <each> parts	
	totw = each*size            #therefore the total number of parts is <totw>
	work = np.append(work,np.zeros((totw-MTR)*lenchunck))  #and one must append zeros to reach the right total length for work to divide work equally between processors
	work = np.reshape(work,(size,-1))  #work is reshaped work[j] gives the part to be processed by rank j	
	
	###################################
	# Initialization of output arrays #
	###################################
	
	
	
        SHX = np.zeros(nslice)
	sumpowafter   = np.zeros(ntot)
        sumphafafter  = np.zeros(ntot)
	sumphatafter  = np.zeros(ntot)
	sumspecafter  = np.zeros(ntot)        
        for i in range(nslice):
            SHX[i] = np.floor(np.abs(ssc[i+ntot/2 -nslice/2 + SHF])*dr)  #Defines the proper shift to account for spatiotemporal coupling
	    #print i, ds, ssc[2]-ssc[1], dr, SHF
        #exit()
        	
	print 'Parallel filtering... '
	print 'MTR = ',MTR
	
    else:
    
        #####################################################
        # Initializes quantities for processes other than 0 #
	#####################################################
    
        flat = []
	Fmod = []
	Fpha = []
	ntot = 1
	each = 0
        work = None
	Mustdo = None

	
    comm.Barrier()  #Makes sure that everything is synchronized now
    
    ##########################################################
    # Broadcasting and scattering of quantities to processes #
    ##########################################################

    Fmod   = comm.bcast(Fmod, root = 0)
    Fpha   = comm.bcast(Fpha, root = 0)
    
    ntot     = comm.bcast(ntot, root = 0)    
    nslice   = comm.bcast(nslice, root = 0)    
    each     = comm.bcast(each, root = 0)    
    Mustdo   = comm.scatter(work, root = 0)  #any <Mustdo> received by any rank, has length = <each>    

    
    # Defines ancyllary quantities for easier access to data
    SH1   = ntot/2-nslice/2
    SHTOT = ntot/2-nslice/2 + SHF       
    imagy  = MTR * MTR
    MTR2   = 2* MTR * MTR
    
   
    # Array initialization
    filament      = np.zeros(ntot,'complex')    
    Dmod2 = np.zeros(ntot)
    Tmod2 = np.zeros(ntot)
    DphaS = np.zeros(ntot)
    TphaS = np.zeros(ntot)
    filmodS = np.zeros(ntot)
    filphaS = np.zeros(ntot)  
    FfilmodS = np.zeros(ntot)
    FfilphaS = np.zeros(ntot)  
    MypartR = np.array([])
    MypartI = np.array([])
    
    comm.Barrier() #Makes sure that everything is synchronized now
    
    
    ##############################TEST NO FILTER
    #if FiltY == 1 :
    #    
    #    for j in range(len(Fmod[1])):
    #        Fmod[1][j] = 1.0
    #	    Fpha[1][j] = 0.0
    ########################################
    
    ##########################
    # Parallelized filtering #
    ##########################
    
    for count in range(each):  #every process runs a number <each> of parts      
     
        for nn in range(MTR):  #every part consists of MTR subparts
                
            filament[SH1:SH1+nslice] = Mustdo[(count)*MTR*nslice*2+nn*nslice*2:(count)*MTR*nslice*2+nn*nslice*2+nslice] + 1j*Mustdo[(count)*MTR*nslice*2+nn*nslice*2+nslice:(count)*MTR*nslice*2+nn*nslice*2+2*nslice]  #filament definition. This is <ntot> long in order to reach resolution; <filament> lives in the time domain	    	    
	    ftfil  = np.roll(fft.fftshift(fft.fft(filament)), -int(mult/2) )                #The final roll is because fftshift treats differently depending on the length of the array. When filtering with increased resolution, an extra-roll is needed                                    #ftfil is the FT of <filament>: we filter in the frequency domain. When mult = 1, roll basically does nothing	    
            Dmod   = np.abs(ftfil)                                                          #Modulus
            Dpha   = np.angle(ftfil)                                                        #Phase 
	    DphaS  = DphaS + np.real(Dpha)
	    Dmod2  = Dmod2 + abs(Dmod)**2
            filmodS = filmodS + abs(filament)**2
            filphaS = filphaS + np.angle(filament) 
            Tmod   = Dmod * Fmod[1]                                                         #Filtered modulus
            Tpha   = Dpha + Fpha[1]                                                         #Filtered phase                     
            Tmod2  = Tmod2 + abs(Tmod)**2		            
            TphaS  = TphaS + np.real(Tpha)			    
            Ffilament = fft.ifft(fft.ifftshift(    np.roll(ftfil * Fmod[1]*np.exp(1j*Fpha[1]), int(mult/2))   ))  #NB: it rolls back too!   #back to the time domain                            
	    FfilmodS = FfilmodS + abs(Ffilament)**2
	    FfilphaS = FfilphaS + np.angle(Ffilament)                
	    MypartR = np.append(MypartR,np.real(Ffilament[SHTOT:SHTOT+nslice]))             #Defines the output field part calculated by a certain process
	    MypartI = np.append(MypartI,np.imag(Ffilament[SHTOT:SHTOT+nslice]))	   
	    if rank==int(size/2.0) and count == int(each/2.0) and nn == int(MTR/2.0):
	        f11 = open(filenamePout+'.axis.dat','w')
	        f12 = open(filenamePphout+'.axix.dat','w')
		ssc  = sscale(ntot, ds)
		for i in range(len(Ffilament)):	
	            f11.write('%s ' %ssc[i] + '%s' %np.abs(Ffilament[i])**2 +'\n')
	            f12.write('%s ' %ssc[i] + '%s' %np.angle(Ffilament[i]) +'\n')
		f11.close()
		f12.close()
    
    #############################################
    # Gets partial quantities back to process 0 #
    #############################################	    
    
    comm.Barrier()
    
    #print MypartR
    
    ResFfilR      = comm.gather(MypartR, root = 0)   
    ResFfilI      = comm.gather(MypartI, root = 0)
    sumpowafter   = comm.reduce(FfilmodS, root = 0)    
    sumphatafter  = comm.reduce(FfilphaS, root = 0)   
    sumspecafter  = comm.reduce(Tmod2, root = 0)    
    sumphafafter  = comm.reduce(TphaS, root = 0) 
    
    #############################################################################
    #                                                                           #
    # Note: ResFfilR[j][(m*MTR+n)*nslice+i] or ResFfilI[j][(e*MTR+n)*nslice+i]  #
    #                                                                           #
    # access the real and imaginary parts of the field filtered                 #
    #                                                                           #
    # by process j (0..size-1)                                                  #
    # at transverse position identified by e (0.. each-1) and n (0.. MTR-1)     #
    # and longitudinal position i (0.. nslice-1)                                #
    #                                                                           #
    #############################################################################
    
        
    if rank == 0: #Process 0 now puts all together; all other processes are idle from now on
        
	print 'Putting data together...'
        
        print 'SHTOT', SHTOT*ds
	print 'SHTOT+nslice',(SHTOT+nslice)*ds
	print 'nslice*ds',nslice*ds
	print 'size ',size
	print 'MTR ',MTR
	print 'each ',each
	print 'nslice ',nslice
	print 'size*MTR*each ', size*MTR*each
	print len(ResFfilR)
	print len(ResFfilR[0])
    	ResFfilR = np.reshape(ResFfilR,(size*MTR*each,nslice))  #Reshapes 
	ResFfilI = np.reshape(ResFfilI,(size*MTR*each,nslice))  #Reshapes 
	
	#############################################################################
        #                                                                           #
        # Note: Now, ResFfilR[m*MTR+n][i] or ResFfilI[m*MTR+n][i]                   #
        #                                                                           #
        # access the real and imaginary parts of the field filtered                 #
        # at transverse position identified by m (0.. MTR-1) and n (0.. MTR-1)      #
        # and longitudinal position i (0.. nslice-1).                               #
        # The other possible values of m are filled up with zeros and are a         # 
	# consequence of the need to divide the work equally to proceess            #
        #                                                                           #
        #############################################################################

	for n in range(MTR):
	    for m in range(MTR):
	        slices[0:nslice,n,m] = ResFfilR[m*MTR+n][0:nslice] + 1j * ResFfilI[m*MTR+n][0:nslice]
	        
	    
        for i in range(nslice):
            slices[i] = np.roll(slices[i],-int(SHX[i]),axis=0)
  
	
    	
	#######################
	# Print data to files #
	#######################
	
	
        print 'Writing radiation file...'
	print nslice
	print len(slices)
	
    #nslice = 1856
    writeRadiationFile_mpi(comm, filenameF, slices, [nslice, MTR, MTR])
    '''
    confirmedok = 1
    while (confirmedok != 0):	    
        confirmedok = writeRadiationFile_mpi(comm, filenameF, slices, [nslice, MTR, MTR])
	comm.Barrier()
	if (rank == 0 and confirmedok != 0): print 'Writing failed; attempting once more...'
    '''	
    if rank == 0:
        print 'Writing output data to file...'
        	
	
	
    	f2 = open(filenameSphout, 'w')
	f8 = open(filenamePout, 'w')
        f9 = open(filenamePphout, 'w')
        f10 = open(filenameSout, 'w')
	
	print 'This is the new version'
	print 'fin s =',ssc[-1]
	print 'in s = ',ssc[0]

	for i in range(len(k)):	              
            f2.write('%s ' %(2*np.pi/k[i]) + '%s' %sumphafafter[i] +'\n')
	    f8.write('%s ' %ssc[i] + '%s' %sumpowafter[i] +'\n')
            f9.write('%s ' %ssc[i] + '%s' %sumphatafter[i] +'\n')
	    f10.write('%s ' %(2*np.pi/k[i]) + '%s' %sumspecafter[i] +'\n')
	    
	f2.close()
    	f8.close()
    	f9.close()
    	f10.close()
	
	print 'End of parallel filtering procedure.'
	
