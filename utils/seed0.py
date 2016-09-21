'''
seed.py

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
    filenamePin    = ARGS[7]    #Input power (summed over transverse coord)
    filenamePphin  = ARGS[8]    #Input t-domain phase (summed over transverse coord)
    filenameSin    = ARGS[9]    #Input spectrum (summed over transverse coord)
    filenameSphin  = ARGS[10]   #Input spectral phase (summed over transverse coord)
    filenamePout   = ARGS[11]   #Output power after seed (summed over transverse coord)
    filenamePphout = ARGS[12]   #Output t-domain phase (summed over transverse coord)
    filenameSout   = ARGS[13]   #Output spectrum after seed (summed over transverse coord)
    filenameSphout = ARGS[14]   #Output specftral phase after seed (summed over transverse coord)
    
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
import sys, os

sys.path.append('/data/netapp/xfel/gianluca/products/ocelot') 
    
    
import numpy as np
import scipy as sp
import scipy.fftpack as fft

import copy
import struct



#BEGIN FUNCTION

def filterfield(filename       = 'D:/Genesis/test_data/simulation.gout.dfl',
                filenameF      = 'D:/Genesis/test_data/simulationFilter.gout.dfl',
                filenameM      = 'D:\Python\Scripts\Genesis - scripts\ModT_fig.dat',
                filenameP      = 'D:\Python\Scripts\Genesis - scripts\PhaseT_fig.dat',
                filenameMR     = 'D:\Python\Scripts\Genesis - scripts\ModT_RES.dat',
                filenamePR     = 'D:\Python\Scripts\Genesis - scripts\PhaseT_RES.dat',
                #filenamePin    = 'D:\Python\Scripts\Genesis - scripts\Pin.dat',
                #filenamePphin  = 'D:\Python\Scripts\Genesis - scripts\Phin.dat',
                #filenameSin    = 'D:\Python\Scripts\Genesis - scripts\Sin.dat',
                #filenameSphin  = 'D:\Python\Scripts\Genesis - scripts\Sphin.dat',
                filenamePout   = 'D:\Python\Scripts\Genesis - scripts\Pout.dat',
                filenamePphout = 'D:\Python\Scripts\Genesis - scripts\Pphout.dat',
                filenameSout   = 'D:\Python\Scripts\Genesis - scripts\Sout.dat',
                filenameSphout = 'D:\Python\Scripts\Genesis - scripts\Sphout.dat',
                #filenamePout_S   = 'D:\Python\Scripts\Genesis - scripts\Pout.dat',
                #filenamePphout_S = 'D:\Python\Scripts\Genesis - scripts\Pphout.dat',
                #filenameSout_S   = 'D:\Python\Scripts\Genesis - scripts\Sout.dat',
                #filenameSphout_S = 'D:\Python\Scripts\Genesis - scripts\Sphout.dat',
                lam0 = 1.5011699558999657e-10,
                MTR = 51,
                mult = 21,
                ds = 1,
                dk = 1,
                SHF = 0,
                nslice=1805,
                dr = 0.0):
    '''
    Function name: filterfieldfilename       = 'D:/Genesis/test_data/simulation.gout.dfl',
                                filenameF      = 'D:/Genesis/test_data/simulationFilter.gout.dfl',
                                filenameM      = 'D:\Python\Scripts\Genesis - scripts\ModT_fig.dat',
                                filenameP      = 'D:\Python\Scripts\Genesis - scripts\PhaseT_fig.dat',
                                filenameMR     = 'D:\Python\Scripts\Genesis - scripts\ModT_RES.dat',
                                filenamePR     = 'D:\Python\Scripts\Genesis - scripts\PhaseT_RES.dat',
                                filenamePin    = 'D:\Python\Scripts\Genesis - scripts\Pin.dat',
                                filenamePphin  = 'D:\Python\Scripts\Genesis - scripts\Phin.dat',
                                filenameSin    = 'D:\Python\Scripts\Genesis - scripts\Sin.dat',
                                filenameSphin  = 'D:\Python\Scripts\Genesis - scripts\Sphin.dat',
                                filenamePout   = 'D:\Python\Scripts\Genesis - scripts\Pout.dat',
                                filenamePphout = 'D:\Python\Scripts\Genesis - scripts\Pphout.dat',
                                filenameSout   = 'D:\Python\Scripts\Genesis - scripts\Sout.dat',
                                filenameSphout = 'D:\Python\Scripts\Genesis - scripts\Sphout.dat',
                                lam0 = 1.5011699558999657e-10,
                                MTR = 51,
                                mult = 21,
                                ds = 1,
                                dk = 1,
                                SHF = 0,
                                nslice=1805,
                                dr = 0.0)
                                
                Description  :  Filters the field file at each transverse position saves it back. The function also obtains spectrum and power before and after filtering, summed over the transverse direction. The filter is rescaled to obtain proper FT.

                Arguments:

                     -filename      : input fld file
                     -filenameF     : output filtered fld file
                     -filenameM     : input modulus of the filter function (two ascii columns [k value])
                     -filenameP     : input phase of the filter function (two ascii columns [k value])
                     -filenameMR    : output modulus of the rescaled filter
                     -filenamePR    : output phase of the rescaled filter
                     -filenamePin   : input power averaged over transverse direction
                     -filenamePphin : input time phase (in the power) averaged over transverse direction
                     -filenameSin   : input spectrum averaged over transverse direction
                     -filenameSphin : input spectral phase averaged over transverse direction
                     -filenamePphout: output time phase (in the power) averaged over transverse direction
                     -filenamePout  : input power averaged over transverse direction
                     -filenameSout  : output spectrum averaged over transverse direction
                     -filenameSphout: output spectral phase averaged over transverse direction
                     -lam0          : central wavelength of the filter
                     -MTR           : transverse number of points
                     -mult          : multiplication factor (see rescalefilter)
                     -ds            : step size in time (abscissa s = ct) domain
                     -dk            : step size in momentum domain
                     -SHF           : field shift in number of slices
                     -nslice        : number of slices
                     -dr            : defines the spatio-temporal cupling; dr must be set as 1/(tmesh*np.tan(thB)),
                                      with tmesh = transverse mesh and thB = Bragg angle. Then the s-dependent transverse shift (in slices) will be
                                      defined as SHX[i] = np.floor(np.abs(ssc[i+nslice*(mult-1)/2 + SHF])*dr), where i is the slice number (see the function body)
                     

                Date revised: 2013.5.24
    '''

    ###################################################################################
    #                                                                                 #
    # FOR FURTHER COMMENTS AND EXPLANATIONS, PLEASE READ COMMENTS IN THE MAIN PROGRAM #
    #                                                                                 #
    ###################################################################################

    print (filename)       
    print (filenameF)   
    print (filenameM)    
    print (filenameP)     
    print (filenameMR)    
    print (filenamePR)    
    #print filenamePin   
    #print filenamePphin  
    #print filenameSin    
    #print filenameSphin  
    print (filenamePout)   
    print (filenamePphout) 
    print (filenameSout)   
    print (filenameSphout) 
    #print filenamePout_S  
    #print filenamePphout_S 
    #print filenameSout_S  
    #print filenameSphout_S
     
    ########################
    # Reads radiation file #
    ########################    

    k0 = 2*np.pi/lam0
    print ('Single process field filtering. Recommended on 64bit processor.')
    print ('Reading radiation file: ',filename,' with ',MTR,'points...')
    slices =  readRadiationFile_my(filename, MTR)   
    print ('...finished reading radiation file.')
    print ('Prepare for filtering...')
    
    nslice = len(slices)   #better use the number in Genesis input file though!
    #flat = np.ravel(slices) #This 64bit version uses a flattened array of data to make the procedure more time-efficient
   
    
    ######################	
    # Defines the filter #
    ######################
    
    ntot = nslice*mult
    ntot_old = ntot
    
    if ntot/2-nslice/2 + SHF < 0: ntot = (-SHF+nslice/2)*2 + 1
    ntot = int(pow(2,np.ceil(np.log(1.0*ntot)/np.log(2.0))))
    
    dk_old = dk
    #print dk, ' ',dk_old
    #print ntot_old,' ',ntot
    
    dk = dk_old * ntot_old/ntot
    #ntot_old -> nslice_old*mult
    #print dk
    
    
    #print ntot 
   
    
 
    k  = kscale(ntot, k0, dk)
    ssc  = sscale(ntot, ds)
    
    Fmod    = rescalefilter(filenameM,k)
    Fpha    = rescalefilter(filenameP,k)
    f1 = open(filenameMR, 'w')
    f3 = open(filenamePR, 'w')
    for i in range(len(k)):
        f1.write('%s ' %(2*np.pi/k[i]) + '%s' %Fmod[1][i] +'\n')
        f3.write('%s ' %(2*np.pi/k[i]) + '%s' %Fpha[1][i] +'\n')
    f1.close()
    f3.close()
    
    ###################################
    # Initialization of output arrays #
    #     and ancillary quantities    #
    ###################################

    filament      = np.zeros(ntot,'complex128')
    sumpowbefore  = np.zeros(ntot)
    sumpowafter   = np.zeros(ntot)
    sumspecbefore = np.zeros(ntot)
    sumspecafter  = np.zeros(ntot)
    sumphatbefore = np.zeros(ntot)
    sumphatafter  = np.zeros(ntot)
    sumphafbefore = np.zeros(ntot)
    sumphafafter  = np.zeros(ntot)
    SHX = np.zeros(nslice)
    for i in range(nslice):
        SHX[i] = np.floor(np.abs(ssc[i-ntot/2-nslice/2 + SHF])*dr)

    #############
    # Filtering #
    #############
        
    print ('Filtering...')
    SH1   = int(ntot/2 - nslice/2)
    SHTOT = int(ntot/2-nslice/2 + SHF)
    print (SHTOT)
    print (SHF)
    
    
    import time
    
    
    
    t2=0
    t=0


    
    import pyfftw
    import multiprocessing
    
    t1 = time.time()
    ncores = multiprocessing.cpu_count()
    pyfftw.interfaces.cache.enable()
    pyfftw.interfaces.cache.set_keepalive_time(30)
    print ('ncores '+str(ncores))
    
    
   
    
    for nn in range(MTR):
        print ('Number '+str (nn)+' of '+str(MTR)+' done.')
        for mm in range(MTR):
           
            #transv = mm * MTR + nn
           
            filament[SH1:SH1+nslice] = slices[0:nslice,nn,mm]
            ftfil = np.fft.fftshift(pyfftw.interfaces.numpy_fft.fft(filament, overwrite_input=False, planner_effort='FFTW_MEASURE', threads=ncores))
            
            #ftfil = np.fft.fftshift(np.fft.fft(filament))
            
           
            Dmod   = np.abs(ftfil)
            Dpha   = np.angle(ftfil)
            Tmod   = Dmod * Fmod[1] 
            Tpha   = Dpha + Fpha[1]	    
            
            finput = fft.ifftshift(ftfil*Fmod[1]*np.exp(1j*Fpha[1]))
           
            Ffilament = pyfftw.interfaces.numpy_fft.ifft(finput, overwrite_input=False, planner_effort='FFTW_MEASURE', threads=ncores)
            #Ffilament = np.fft.ifft(fft.ifftshift(ftfil*Fmod[1]*np.exp(1j*Fpha[1])))


            slices[0:nslice,nn,mm] =  Ffilament[SHTOT:SHTOT+nslice]
            sumpowbefore  = sumpowbefore  + abs(filament)**2
            sumpowafter   = sumpowafter   + abs(Ffilament)**2    
            sumphatbefore = sumphatbefore  + np.angle(filament)
            sumphatafter  = sumphatafter   + np.angle(Ffilament)
            sumspecbefore = sumspecbefore + abs(Dmod)**2
            sumspecafter  = sumspecafter  + abs(Tmod)**2
            sumphafbefore = sumphafbefore + np.real(Dpha)
            sumphafafter  = sumphafafter  + np.real(Tpha)

    t2 = time.time()
    print ('time ')
    print (t2-t1)

    print ('Introducing spatiotemporal coupling...')
    for i in range(nslice):
        slices[i] = np.roll(slices[i],-int(SHX[i]),axis=0)


    #######################
    # Print data to files #
    #######################
    
    print ('Writing radiation file...')
    writeRadiationFile_my(filenameF,slices)
    print ('...Finished writing radiation file.')
    
    print ('Writing output data to file...'  )   
    f2 = open(filenameSphout, 'w')
    #f4 = open(filenamePout, 'w')
    #f5 = open(filenamePphin, 'w')
    #f6 = open(filenameSin, 'w')
    #f7 = open(filenameSphin, 'w')
    f8 = open(filenamePout, 'w')
    f9 = open(filenamePphout, 'w')
    f10 = open(filenameSout, 'w')        
    #f11 = open(filenameSphout_S, 'w')        
    #f12 = open(filenamePout_S, 'w')
    #f13 = open(filenamePphout_S, 'w')
    #f14 = open(filenameSout_S, 'w') 
    
    for i in range(len(k)):
        ##f1.write('%s ' %(2*np.pi/k[i]) + '%s' %Fmod[1][i] +'\n')
        f2.write('%s ' %(2*np.pi/k[i]) + '%s' %sumphafafter[i] +'\n')
        ##f3.write('%s ' %(2*np.pi/k[i]) + '%s' %Fpha[1][i] +'\n')
        #f4.write('%s ' %ssc[i] + '%s' %sumpowbefore[i] +'\n')
        #f5.write('%s ' %ssc[i] + '%s' %sumphatbefore[i] +'\n')
        #f6.write('%s ' %(2*np.pi/k[i]) + '%s' %sumspecbefore[i] +'\n')
        #f7.write('%s ' %(2*np.pi/k[i]) + '%s' %sumphafbefore[i] +'\n')
        f8.write('%s ' %ssc[i] + '%s' %sumpowafter[i] +'\n')
        f9.write('%s ' %ssc[i] + '%s' %sumphatafter[i] +'\n')
        f10.write('%s ' %(2*np.pi/k[i]) + '%s' %sumspecafter[i] +'\n')
        #if i in range(nslice):
        #    f11.write('%s ' %(2*np.pi/k[i]) + '%s' %sumphafafter[i] +'\n')
        #    f12.write('%s ' %ssc[i] + '%s' %sumpowafter[i] +'\n')
        #    f13.write('%s ' %ssc[i] + '%s' %sumphatafter[i] +'\n')
        #    f14.write('%s ' %(2*np.pi/k[i]) + '%s' %sumspecafter[i] +'\n')          
    f2.close()
    #f4.close()
    #f5.close()
    #f6.close()
    #f7.close()
    f8.close()
    f9.close()
    f10.close()
    #f11.close()
    #f12.close()
    #f13.close()
    #f14.close()
     
    print ('End of single process field filtering. Recommended on 64bit processors.')
    
#END FUNCTION



   
#BEGIN FUNCTION   

def sscale(Npts, ds = 5.10373e-9):
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


def kscale(Npts = 10000,k0=4.188790204e10, dk = (5*4.60612e-06 * 21)**(-1)):
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


    # f = file(filename, 'r')  # iterate over the lines in the file
    f = open(filename, 'r')  # iterate over the lines in the file
    columns = []
    for line in f:
        # split the line into a list of column values
        #columns1 = line.split('\t')
        columns1 = line.split(' ')
        # clean any whitespace off the items
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

    Interdata = np.interp(k, dataX, dataY)
    spits = [k, Interdata]

    return spits

#END FUNCTION




#BEGIN FUNCTION


def readRadiationFile_my(fileName, npoints=151):
    import numpy as np
    b=np.fromfile(fileName,dtype=complex)
    slice_num=b.shape[0]/npoints/npoints
    c=b.reshape(slice_num,npoints,npoints)
    return c
#END  FUNCTION

def writeRadiationFile_my(filename,rad):
    print ('    writing radiation file' )
    #a new backward compatible version ~10x faster
    #    print '        - writing to ', filename
    d=rad.flatten()
    d.tofile(filename,format='complex')
