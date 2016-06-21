'''
interface to genesis
'''

import struct
from copy import copy
import time, os
from ocelot.rad.fel import *
from ocelot.cpbd.beam import Beam, gauss_from_twiss
import ocelot.utils.reswake as w
from ocelot.common.math_op import *

h = 4.135667516e-15
c = 299792458.0

inputTemplate = "\
 $newrun \n\
 aw0   =  __AW0__ \n\
 xkx   =  __XKX__\n\
 xky   =  __XKY__\n\
 wcoefz =  __WCOEFZ__\n\
 xlamd =  __XLAMD__\n\
 fbess0 =  __FBESS0__\n\
 delaw =  __DELAW__\n\
 iertyp =  __IERTYP__\n\
 iwityp =  __IWITYP__\n\
 awd   =  __AWD__ \n\
 awx   =  __AWX__\n\
 awy   =  __AWY__\n\
 iseed =  __ISEED__\n\
 npart =  __NPART__\n\
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
 conditx =  __CONDITX__\n\
 condity =  __CONDITY__\n\
 bunch =  __BUNCH__\n\
 bunchphase =  __BUNCHPHASE__\n\
 emod =  __EMOD__\n\
 emodphase =  __EMODPHASE__\n\
 xlamds =  __XLAMDS__\n\
 prad0 =  __PRAD0__\n\
 zrayl =  __ZRAYL__\n\
 zwaist =  __ZWAIST__\n\
 ncar  =  __NCAR__\n\
 lbc   =  __LBC__\n\
 rmax0 =  __RMAX0__\n\
 dgrid =  __DGRID__\n\
 nscr  =  __NSCR__\n\
 nscz  =  __NSCZ__\n\
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
 ildgam =  __ILDGAM__\n\
 ildpsi =  __ILDPSI__\n\
 ildx  =  __ILDX__\n\
 ildy  =  __ILDY__\n\
 ildpx =  __ILDPX__\n\
 ildpy =  __ILDPY__\n\
 itgaus =  __ITGAUS__\n\
 nbins =    __NBINS__\n\
 igamgaus =  __IGAMGAUS__\n\
 lout  = __LOUT__\n\
 iphsty =  __IPHSTY__\n\
 ishsty =  __ISHSTY__\n\
 ippart =  __IPPART__\n\
 ispart =  __ISPART__\n\
 ipradi =  __IPRADI__\n\
 isradi =  __ISRADI__\n\
 idump =  __IDUMP__\n\
 iotail =    __IOTAIL__\n\
 nharm =    __NHARM__\n\
 curpeak =  __CURPEAK__\n\
 curlen =  __CURLEN__\n\
 ntail = __NTAIL__\n\
 nslice = __NSLICE__\n\
 __SHOTNOISE__\n\
 isntyp =  __ISNTYP__\n\
 iall  =  __IALL__\n\
 __ITDP__\n\
 ipseed =   __IPSEED__\n\
 iscan =  __ISCAN__\n\
 nscan =  __NSCAN__\n\
 svar  =  __SVAR__\n\
 isravg =    __ISRAVG__\n\
 isrsig =    __ISRSIG__\n\
 cuttail = __CUTTAIL__\n\
 eloss =  __ELOSS__\n\
 version =  __VERSION__\n\
 ndcut =  __NDCUT__\n\
 idmpfld =    __IDMPFLD__\n\
 idmppar =    __IDMPPAR__\n\
 ilog  =  __ILOG__\n\
 ffspec =  __FFSPEC__\n\
 convharm =    __CONVHARM__\n\
 ibfield =  __IBFIELD__\n\
 imagl =  __IMAGL__\n\
 idril =  __IDRIL__\n\
 alignradf =  __ALIGNRADF__\n\
 offsetradf =  __OFFSETRADF__\n\
 multconv =  __MULTCONV__\n\
__BEAMFILE__\n\
__PARTFILE__\n\
__FIELDFILE__\n\
__RADFILE__\n\
__DISTFILE__\n\
__MAGFILE__\n\
 outputfile ='run.__RUNID__.gout'\n\
 filetype ='ORIGINAL'\n\
 $end\n"

 # iallharm =  __IALLHARM__\n\
 # iharmsc =  __IHARMSC__\n\
 # pradh0 =  __PRADH0__\n\
 
class GenesisInput:
    
    def __init__(self):
        
    # defaults
        
        self.runid = 0
        self.type = 'steady'
        
        #undulator
        self.aw0 = 0.735 #The normalized, dimensionless rms undulator parameter, defined by AW0 = (e/mc)(Bu/ku), where e is the electron charge, m is electron mass, c is speed of light, ku=2pi/lambdau is the undulator wave number, lambdau is the undulator period. Bu is the rms undulator field with Bu = Bp/2 for a planar undulator and Bu = Bp for a helical undulator, where Bp is the on-axis peak field. 
        self.awd = 0.735 #A virtual undulator parameter for the gap between undulator modules.
        self.wcoefz = [0,0,0]   #(1-[m]) Start of undulator tapering.  Note that tapering is applied, even the magnetic lattice is defined by an external file.
                                #(2-[ ]) The relative change of the undulator field over the entire taper length (AW(exit) = (1 -WCOEFZ(2))
                                #(3) The taper model: 1 for linear taper, 2 for quadratic taper,
        self.iertyp =0 # Type of undulator field errors.
        self.iwityp =0 # the undulator type. A value of zero indicates a planar undulator, any other value a helical one. 
        self.xkx   =  0 #Normalized natural focusing of the undulator in x. Common values are XKX = 0.0, XKY = 1.0 for a planar undulator or XKX, XKY = 0.5 for a helical undulator, but might vary if focusing by curved pole faces is simulated. The values should fulfill the constraint XKX + XKY = 1.0. 
        self.xky   =  1 #Normalized natural focusing of the undulator in y 
        self.delaw =  0 #RMS value of the undulator field error distribution. A value of zero disables field errors. 
        self.nwig = 98   # The number of periods within a single undulator module. The product of NWIG and XLAMD defines the length of the undulator module. 
        self.nsec = 1   # The number of sections of the undulator.
        self.awx   =  0 #Maximum offset in x for undulator module misalignment. The error for each individual module follows a uniform distribution
        self.awy   =  0 #Maximum offset in y for undulator module misalignment. The error for each individual module follows a uniform distribution
        self.iseed = -1  #initial seeding of the random number generator for field errors

        #electron beam
        self.curpeak = 2.5E+03 #Peak current of the electron beam. Time-independent simulations enforces a constant current. 
        self.curlen = 7E-06 # Bunch length of the current profile. If CURLEN is positive a Gaussian distribution is generated with an rms length given by CURLEN. A negative or zero value yield a constant profile. 
        self.npart = 8192   # number of macroparticles per slice. NPART must be a multiple of 4*NBINS
        self.gamma0 =  3.424658E+04   # The mean energy of the electron beam in terms of the electron rest mass energy. 
        self.delgam = 5.0E-3  #The RMS value of the energy distribution in terms of electron rest mass energy. 
        self.rxbeam = 100e-6 #The rms value in x of the transverse, spatial distribution. 
        self.rybeam = 100e-6 #The rms value in y of the transverse, spatial distribution.
        self.emitx = 1.0e-7 #The normalized rms emittance in x
        self.emity = 1.0e-7 #The normalized rms emittance in y
        self.alphax = 0.0 #Rotation of the transverse phase space distribution in x according to the standard definition ALPHAX = - < xx' > GAMMA0 / EMITX. 
        self.alphay = 0.0 #Rotation of the transverse phase space distribution in y according to the standard definition ALPHAY = - < yy' > GAMMA0 / EMITY. 
        self.xbeam = 0.0 #Transverse position in x of the electron beam at the undulator entrance with respect to the undulator axis
        self.ybeam = 0.0 #Transverse position in y of the electron beam at the undulator entrance with respect to the undulator axis
        self.pxbeam = 0.0 #Average normalized transverse momentum x of the electron beam at the undulator entrance. The momenta are defined as PXBEAM = betax*gamma where betax = c*v_x is the average transverse velocity and gamma the Lorenz factor of the electron energy. 
        self.pybeam = 0.0 #y
        self.cuttail = -1.0 #Cut in the transverse phase space in measures of the rms size to collimate transverse beam tails/halos. The cut is applied after the loading and beam current is set accordingly to the reduced number of macro particles. It is disabled if the value is negative or the electron beam is imported from an external file. 
        self.conditx =  0.0 #the correlation strength between the amplitude of the eletron's betatron oscillation x and its energy.
        self.condity =  0.0 # y
        self.bunch =  0.0 #Initial value of the bunching factor, when quite loading is used. 
        self.bunchphase =  0.0 # Phase of initial bunching, when quite loading is used. 
        self.emod =  0.0 #Initial energy modulation, when quite loading is used. The value is the modulation amplitude in terms of gamma. 
        self.emodphase =  0.0 #Phase of initial energy modulation, when quite loading is used. 

        #particle loading
        self.ildpsi =    7 #Indices of the Hammersley sequence bases for loading the particle phase. 
        self.ildgam =    5 #Hammersley base for loading the energy distribution. 
        self.ildx  =    1 #Hammersley base for loading the distribution in x.
        self.ildy  =    2 #Hammersley base for loading the distribution in y.
        self.ildpx =    3 #Hammersley base for loading the distribution in px.
        self.ildpy =    4 #Hammersley base for loading the distribution in py.
        self.itgaus =    1 #Defines distribution profile in the transverse variables. The available distributions are: Gaussian (1) Uniform (2) Parabolic (otherwise) 
        self.igamgaus =    1 #Defines distribution profile in energy. A non-zero value generates a Gaussian distribution and a uniform otherwise.
        self.iall  =    0 #A non-zero value of IALL enforces that all slices are starting with the same element of the Hammersley sequences.
        self.ipseed = -1  #Initial seed for the random number generator used for particle phase fluctuation (shot noise). GENESIS 1.3 requires a re-initialization of the random number generator to guarantee the same loading whether magnetic field errors are used or not. 
        self.nbins = 4  # Number of bins in the particle phase. The value has to be at least 4 or larger than (2+2n), depending on whether the bunching at the nth harmonics is needed for space charge calculations or output.
        
        #radiation
        self.xlamds = 1E-9 #The resonant radiation wavelength.
        self.prad0 = 0 #The input power of the radiation field.
        self.pradh0 = 0 #Radiation power for seeding with a harmonics, defined by NHARM.
        self.zrayl = 0.5 #The Rayleigh length of the seeding radiation field.
        self.zwaist = 0 # Position of the waist of the seeding radiation field with respect to the undulator entrance.

        #mesh
        self.ncar = 151   # The number of grid points for the radiation field along a single axis. The total number for the mesh is NCAR^2
        self.lbc   =    0 # Flag to set the boundary condition for the field solver. The options are Direchlet boundary condition (LBC = 0) and Neumann boundary condition (otherwise).
        self.nscr  =    0 #Number of azimuthal modes for space charge calculation. 
        self.nscz  =    0 #Number of longitudinal Fourier components for space charge calculation. NSCZ must be greater than 0 to include space charge but less than (2NBINS+2) for valid results. 
        self.nptr = 40   # Number of radial grid points, on which the space charge field is evaluated.
        self.rmax0 = 9.0   # mesh size in units of radiation+beam sizes
        self.dgrid = 0.0   # explicit trnsverse mesh size overruling the calculation by the RMAX0-parameter. 
        
        
        #focusing
        self.quadf = 0   # The field strength of (F-)quadrupoles, focusing in the x-plane. 
        self.quadd = 0   # The fields strength of (D-)quadrupoles, defocusing in the x-plane.
        self.fl = 0   # Length of the F-quadrupoles in measures of the undulator period.
        self.dl = 0   # Length of the D-quadrupoles in measures of the undulator period. 
        self.drl = 0   # Drift length between F- and D-quadrupoles in measures of undulator period.
        self.qfdx = 0   # Maximum transverse misplacement of the quadrupoles in x-direction. A random offset between [-,] in x is applied to every quadrupole.
        self.qfdy = 0   # Maximum transverse misplacement of the quadrupoles in y-direction, respectively.
        self.solen = 0   # On-axis field strength of a superimposed solenoid field. 
        self.sl = 0   # Length of solenoid field in measures of undulator period. The solenoid is aligned to the beginning of each undulator section. 
        self.f1st = 0   # Position within a FODO cell, where GENESIS 1.3 starts the FODO cell lattice
        
        #simulation
        self.version =  0.1 #Used for backward compatibility of the input decks.
        self.zsep = 20   #Separation of beam slices in measures of the radiation wavelength. ZSEP must be a multiple of DELZ.
        self.nslice = 1000   # Total number of simulated slices. It defines the time window of the simulation with NSLICE * ZSEP * XLAMDS/c
        self.ntail = - self.nslice / 2 #Position of the first simulated slice in measures of ZSEP*XLAMDS. GENESIS 1.3 starts with the tail side of the time window, progressing towards the head. Thus a negative or positive value shifts the slices towards the tail or head region of the beam, respectively.
        self.delz = 1.0   # Integration step size in measure of the undulator period length.
        self.zstop = 256.0 #Defines the total integration length. If the undulator length is shorter than ZSTOP or ZSTOP is zero or negative, the parameter is ignored and the integration is performed over the entire undulator. 
        self.iorb = 0   # enforce orbit correction. For any non-zero value the offset due to the wiggle motion is taken into account for the interaction between electron beam and radiation field.
        self.isravg = 1 # If set to a non-zero value the energy loss due to spontaneous synchrotron radiation is included in the calculation. 
        self.isrsig = 1 # If set to a non-zero value the increase of the energy spread due to the quantum fluctuation of the spontaneous synchrotron radiation is included in the calculation.
        self.eloss = 0 # Externally applied energy loss of the electron beam.
        self.nharm = 1 #Enables the calculation of harmonics up to the one, specified by NHARM. Note that the number of NBINS has to be at least twice as large as NHARM to allow for the correct representation of the harmonics. Note also that this parameter does not enable automatically the output. For that the corresponding bit in LOUT has to be set as well.
        self.iscan =    0 #Selects the parameter for a scan over a certain range of its value 
        #(1.GAMMA0 2.DELGAM 3.CURPEAK 4.XLAMDS 5.AW0 6.ISEED 7.PXBEAM 8.PYBEAM 9.XBEAM 10.YBEAM 11.RXBEAM 12.RYBEAM 13.XLAMD 14.DELAW 15.ALPHAX 16.ALPHAY 17.EMITX 18.EMITY 19.PRAD0 20.ZRAYL 21.ZWAIST 22.AWD 23.BEAMFILE 24.BEAMOPT 25.BEAMGAM)  
        #self.scan = '' #By supplying the parameter name to scan over it overrules the setting of ISCAN
        self.nscan =    3 #Number of steps per scan. 
        self.svar = 0.01  #Defines the scan range of the selected scan parameter. The parameter is varied between (1-SVAR) and (1+SVAR) of its initial value.
        
        #I/O
        self.iphsty =    1 # Generate output in the main output file at each IPHSTYth integration step. To disable output set IPHSTY to zero. 
        self.ishsty =    1 # Generate output in the main output file for each ISHSTYth slice. 
        self.ippart =    0 # Write the particle distribution to file at each IPPARTth integration step. To disable output set IPPART to zero. The filename is the same of the main outputfile + the extension '.par'. 
        self.ispart =    0 # Write the particle distribution to file for every ISPART slice. 
        self.ipradi =    0 # Write the field distribution to file at each IPRADIth integration step. To disable output set IPRADI to zero. The filename is the same of the main outputfile + the extension '.fld'. 
        self.isradi =    0 # Write the field distribution to file for every ISRADI slice. 
        self.iotail = 1 # If set to a non-zero value the output time window is the same as the simulated time window. Otherwise the output for the first slices covered by the slippage length is subpressed.
        self.magin = 1   # read in magnetic lattice (If set to a non-zero value the user is prompted to type in the file name containing a explicit description of the magnetic field.)
        self.magout = 0   # output magnetic lattice
        self.idump =    0 # If set to a non-zero value the complete particle and field distribution is dumped at the undulator exit into two outputfiles.
        self.idmpfld= 0 # Similar to IDUMP but only for the field distribution. 
        self.idmppar = 0 # Similar to IDUMP but only for the particle distribution. 
        self.lout=[1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                #  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
                # 1. radiation power
                # 2. logarithmic derivative of the power growth
                # 3. power density at the undulator axis
                # 4. radiation phase at the undulator axis
                # 5. transverse radiation size
                # 6. rms diffraction angle of the radiation
                # 7. beam energy
                # 8. bunching factor
                # 9. beam size in x
                # 10. beam size in y
                # 11. error in energy conservation
                # 12. beam position in x
                # 13. beam position in y
                # 14. energy spread
                # 15. on-axis field intensity in the far field zone
                # 16. bunching at the 2nd harmonic
                # 17. bunching at the 3rd harmonic
                # 18. bunching at the 4th harmonic
                # 19. bunching at the 5th harmonic
                # 20. bunching phase of the 1st harmonic
                # 21. bunching phase of the 2nd harmonic
                # 22. bunching phase of the 3rd harmonic
                # 23. bunching phase of the 4th harmonic
                # 24. bunching phase of the 5th harmonics

        self.beamfile = None
        self.fieldfile = None
        self.partfile = None
        self.distfile = None
        self.radfile = None
        
        self.ndcut =   -1 # ??? If NDCUT is zero, the time-window is adjusted, so that in average NPART/NBINS particles fall in each slice. 
        self.alignradf =    1 # if zero , Genesis 1.3 aligns the radiation field to the electron beam so that the radiaiton field is one ful slippage behind the electron beam.
        self.offsetradf =    0 # slices to shift the electrron beam with respect to the radiation if ALIGNRADF=1.
        self.convharm = 1 #When the particle distribution is imported from a PARTFILE Genesis 1.3 allows the upconversion to a higher harmonics. The harmonic number is specified with CONVHARM and has a defulat value of 1, corresponding to no up-conversion. The user has to make sure that in the input deck XLAMDS is adjusted, according to the new wavelength.
        self.multconv =    0 # If an imported particle distribution from a PARTFILE is up-converted to a higher harmonics the dault behavior is that the number of slices is preserved. This requires that ZSEPis adjusted together with XLAMDS. However if frequency resolution is a problem then a particle distribution can be converted and used multiple times to keep ZSEP constant. The disadvantage is that the CPU execution time is increased as well. 
        
        self.ibfield =  0.0 #When the PARTFILE features is used the imported particle distribution can be tracked through a generic 4 magnet chicane before running the Genesis simulation. The chicane consists out of 4 bending magnets of the field strength IBFIELD and length IMAGL separated by 5 drifts of length IDRIL.
        self.imagl =    0.0 #The length of each bending magnet of the chicane. If the magnet length is set to zero but IDRIL is not the resulting beam line correspond to a simple drift of the length 5 times IDRIL
        self.idril =    0.0 #The length of the 5 drift lengths of the magnetic chicane (three between the magnets and one before and after the magnets).
        
        self.iallharm =    0 #Setting the value to a non-zero value will also include all harmonics between 1 and NHARM
        self.iharmsc =    0 # setting to a non-zero value includes the coupling of the harmonic radiation back to the electron beam for a self-consistent model of harmonics. Enabling this feature will automatically include all harmonics by setting IALLHARM to one.
        self.isntyp =    0 # Non-zero if the user wants to use the Pennman algorithm for the shot noise (which is not recommended).
        
        self.ilog  =    0 #Create a log file.
        self.ffspec =    0 # amplitude/phase values for spectrum calculation: 0 - on-axis power/phase along the pulse, -1 - the same in far field, 1 - near field total power
        
  
        #self.useBeamFile = False

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
            
        if self.partfile != None:
            input = input.replace("__PARTFILE__", " partfile  =  '"+ str(self.partfile)+ "'")
        else:
            input = input.replace("__PARTFILE__\n", "")
            
        if self.distfile != None:
            input = input.replace("__DISTFILE__", " distfile  =  '"+ str(self.distfile)+ "'")
        else:
            input = input.replace("__DISTFILE__\n", "")

        '''
        Fix for radin file GG - 27.05.2016
        '''
        #print 'self.radfile is equal to ', self.radfile

        if self.radfile != None:
            input = input.replace("__RADFILE__", " radfile  =  '"+ str(self.radfile)+ "'")
        else:
            input = input.replace("__RADFILE__\n", "")



        
        if self.magin == 0:
            input = input.replace("__MAGFILE__\n", "")
        else:
            input = input.replace("__MAGFILE__\n", " maginfile ='lattice.inp'\n")
        
        for p in self.__dict__.keys():
            input = input.replace("__"  + str(p).upper() + "__", str(self.__dict__[p]).replace('[','').replace(']','').replace(',',''))
                                
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
        self.filename = ''
            
    
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
        
class GenesisParticles:
    
    def __init__(self):
        self.e = []
        self.ph = []
        self.x = []
        self.y = []
        self.px = []
        self.py = []
        self.z = []
        self.t = []

        self.filename = ''

        
class GenesisBeamDefinition():
    
    def __init__(self):
        self.columns=[]
        self.column_values={}

# class GenesisOutParm():
#
#     def __init__(self):
#         self.v=[]
#         self.mean_S=[]
#         #self.mean_Z=[]
#         self.max_S=[]
#         #self.max_Z=[]
#         self.end_Z=[]
#
#     def scan(self,value):
#         g=GenesisOutParm()
#         g.v=np.array(value)
#         g.mean_S=np.mean(g.v,axis=0)
#         #g.mean_Z=np.mean(g.v,axis=1)
#         g.end_Z=g.v[:,-1]
#         g.max_S=np.amax(g.v,axis=0)
#         return g



''' 
   I/O functions
'''

#def read_particle_file(filename, npart=[]):
##    print npart
#    #new faster function with different output convenstion    
#    particles=GenesisParticles
##    start_time = time.time()
#    tmp=np.fromfile(filename,dtype=float)
#    
##    if npart!=[] and nslice!=[]    
##    if len(tmp)!=npart*nslice*6:
##    print len(tmp)/npart/6
#    nslice=int(len(tmp)/npart/6)
#    #nslice=600
#    tmp=tmp.reshape(nslice,6,npart)
#    particles.e=tmp[:,0,:] #gamma
#    particles.ph=tmp[:,1,:] 
#    particles.x=tmp[:,2,:]
#    particles.y=tmp[:,3,:]
#    particles.px=tmp[:,4,:]
#    particles.py=tmp[:,5,:]
#    particles.filename=filename
#    
#    return particles
def read_particle_file(file_name, nbins=4, npart=[],debug=0):
    print '    reading praticle file' 
    particles=GenesisParticles    
    
    start_time = time.time()
    b=np.fromfile(file_name,dtype=float)
    print("--- read Particles - %s seconds ---" % (time.time() - start_time))
#    print 'b', b.shape
    nslice=int(len(b)/npart/6)
    if debug:
        print '        nslice',nslice
        print '        npart',npart
        print '        nbins',nbins
#    print 'b=',nslice*npart*6
    b=b.reshape(nslice,6,nbins,npart/nbins)    
    particles.e=b[:,0,:,:] #gamma
    particles.ph=b[:,1,:,:] 
    particles.x=b[:,2,:,:]
    particles.y=b[:,3,:,:]
    particles.px=b[:,4,:,:]
    particles.py=b[:,5,:,:]
    particles.filename=file_name
    return particles
    

def readParticleFile_old(fileName, npart, nslice):
    #old file with different output convenstion
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



def read_rad_file(fileName):
    beam = GenesisBeamDefinition()
    
    f=open(fileName,'r')
    null=f.readline()
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
    beam.prad0 = np.array(beam.column_values['PRAD0'])
        
    return beam
    
    
    

def read_beam_file(fileName):
    
    beam = GenesisBeamDefinition()
    
    f=open(fileName,'r')
    null=f.readline()
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


def readRadiationFile(fileName, npoints=151, slice_start=0, slice_end = -1, vartype=complex128):
    #a new backward compatible version ~100x faster
    print '    reading radiation file'    
    import numpy as np
#    print '        - reading from ', fileName
    b=np.fromfile(fileName,dtype=complex).astype(vartype)
    slice_num=b.shape[0]/npoints/npoints
    b=b.reshape(slice_num,npoints,npoints)
    if slice_end == -1:
        slice_end=None
    print '      done'  
    return b[slice_start:slice_end]

    
def readRadiationFile_old(fileName='simulation.gout.dfl', npoints=51, slice_start=0, slice_end = -1, idx=None):
    # temporarily remains here for backup purposes
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


def writeRadiationFile(filename,rad):
    print '    writing radiation file' 
    #a new backward compatible version ~10x faster
#    print '        - writing to ', filename
    d=rad.flatten()
    d.tofile(filename,format='complex')


def writeRadiationFile_old(filename, slices):
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


def readGenesisOutput_old(fileName , readall=True, debug=None):

    print '    reading output file'    
#    print '        - reading from ', fileName

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
            print '      reading optics '
            continue
        
        if tokens[0] == 'Input':
            chunk = 'input'
            print '      reading input parameters'
            continue
        
        #********** output: slice    10
        if tokens[0] == '**********':
            #print 'slice:', tokens[3]
            chunk = 'slices'
            nSlice = int(tokens[3])
            if debug:
                print '      reading slice # ',nSlice
         
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
            #out.parameters[tokens[0]] = tokens[0:]
            #print 'input:', tokens

        if chunk == 'slice':
            vals = map(float,tokens)
            #print vals
            for i in range(0,len(vals)):
                out.sliceValues[nSlice][out.sliceKeys[i]].append(vals[i])
#            out.sliceValues[nSlice][out.sliceKeys[i]].extend(vals)
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

    print '        nSlice', out.nSlices
    print '        nZ', out.nZ

#    print '      processing 1'
    out.power = []
    out.phi = []
    out.power_z = 0*np.array(out.sliceValues[out.sliceValues.keys()[0]]['power'])
    out.power_int = []
    out.max_power = 0.0
    if readall:
        out.r_size = [] #make an np arrays to save further conversion time?
        out.el_energy = []
        out.bunching = []
        out.xrms = []
        out.yrms = []
        out.error = []
        out.el_e_spread = []
        out.p_mid = []
        out.p_int = []
        out.phi_mid = []
        out.increment = []
        #out.increment=GenesisOutput()
#    print '      processing 2'
    for i in xrange(0,out.nSlices):

        pend = out.sliceValues[out.sliceValues.keys()[i]]['power'][-1]
        out.power_int.append(pend)
        out.power.append(out.sliceValues[out.sliceValues.keys()[i]]['p_mid'][-1])
        out.phi.append(out.sliceValues[out.sliceValues.keys()[i]]['phi_mid'][-1])
        out.power_z +=  np.array(out.sliceValues[out.sliceValues.keys()[i]]['power']) / out.nSlices
        if readall:
            out.p_mid.append(out.sliceValues[out.sliceValues.keys()[i]]['p_mid'])
            out.p_int.append(out.sliceValues[out.sliceValues.keys()[i]]['power'])
            out.phi_mid.append(out.sliceValues[out.sliceValues.keys()[i]]['phi_mid'])
            out.increment.append(out.sliceValues[out.sliceValues.keys()[i]]['increment'])
            out.r_size.append(out.sliceValues[out.sliceValues.keys()[i]]['r_size'])
            out.el_energy.append(out.sliceValues[out.sliceValues.keys()[i]]['energy'])
            out.bunching.append(out.sliceValues[out.sliceValues.keys()[i]]['bunching'])
            out.xrms.append(out.sliceValues[out.sliceValues.keys()[i]]['xrms'])
            out.yrms.append(out.sliceValues[out.sliceValues.keys()[i]]['yrms'])
            out.error.append(out.sliceValues[out.sliceValues.keys()[i]]['error'])
            out.el_e_spread.append(out.sliceValues[out.sliceValues.keys()[i]]['e-spread'])

        if out.max_power < pend: out.max_power = pend
    del out.sliceValues
#    print '      processing 3'
    out.s = out('zsep') * out('xlamds') * np.arange(0,len(out.power))
    out.t = out.s/ c *1.e+15
#    out.t = 1.0e+15 * out('zsep') * out('xlamds') / c * np.arange(0,len(out.power))
    out.dt = (out.t[1] - out.t[0]) * 1.e-15
    
    
    out.spec = fft.fft(np.sqrt(np.array(out.power) ) * np.exp( 1.j* np.array(out.phi) ) )
    out.freq_ev = h * fftfreq(len(out.spec), d=out('zsep') * out('xlamds') / c)

    out.power = np.array(out.power)
    out.phi = np.array(out.phi)
    out.power_int = np.array(out.power_int)
    out.z = np.array(out.z)
    out.I = np.array(out.I)
    out.beam_charge=np.sum(out.I*out('zsep')*out('xlamds')/c)
    if readall:
        out.p_mid = np.array(out.p_mid)
        out.p_int = np.array(out.p_int)
        out.phi_mid = np.array(out.phi_mid)
        out.increment = np.array(out.increment)
        out.r_size = np.array(out.r_size)
        out.el_energy = np.array(out.el_energy)+out('gamma0')
        out.bunching = np.array(out.bunching)
        out.xrms = np.array(out.xrms)
        out.yrms = np.array(out.yrms)
        out.error = np.array(out.error)
        out.el_e_spread = np.array(out.el_e_spread)
        # out.power = np.array(out.power)

        if out('dgrid')==0:
            rbeam=sqrt(out('rxbeam')**2+out('rybeam')**2)
            ray=sqrt(out('zrayl')*out('xlamds')/np.pi*(1+(out('zwaist')/out('zrayl')))**2); #not cross-checked
            out.leng=out('rmax0')*(rbeam+ray)*2
        else:
            out.leng=out('dgrid')*2


        # out.energy_GeV=out.el_energy*0.511e-3
        # out.e_spread_GeV=(out.el_e_spread+out.el_energy)*0.511e-3-out.energy_GeV
        # out.e_spread_GeV=out.el_e_spread*0.511e-3
        # parm_names=out.parameters.keys()
# for val in parm_names:
#     if str.isdigit(val[0]) or  val[0]=='$':
#         pass
#     else:
        
        out.filename = fileName[-fileName[::-1].find('/')::]
        
    return out


def readGenesisOutput(fileName , readall=True, debug=None, precision=float):
    out = GenesisOutput()
    out.path = fileName
    out.filename = fileName[-fileName[::-1].find('/')::]

    print(' ')
    print('    reading output file "'+out.filename+'"')
#    print '        - reading from ', fileName

    chunk = ''
    output_unsorted=[] 
    nSlice = 0
    
    wait_attempt=6
    wait_time=10
    while os.path.isfile(fileName)!=True:
        print('!     waiting for "'+out.filename+'" '+str(wait_time)+'s ['+str(wait_attempt)+']')
        time.sleep(wait_time) #wait for the .out file to be assembled
        wait_attempt-=1
        if wait_attempt==0:
            raise Exception('File '+fileName+' not found')
    
    start_time = time.time()    
    f=open(fileName,'r')
    
    null=f.readline()
    for line in f: 
        tokens = line.strip().split()

        if len(tokens) < 1:
            continue

        if tokens[0] == '**********':
            chunk = 'slices'
            nSlice = int(tokens[3])
            if debug:
                print '      reading slice # ',nSlice

        if tokens[0] == 'power':
            chunk = 'slice'
            if len(out.sliceKeys) == 0: #to record the first instance
                out.sliceKeys = copy(tokens)
                print '      reading slice values '
            continue
            
        if tokens[0] == '$newrun':
            chunk = 'input1'
            print '      reading input parameters'
            continue  

        if tokens[0] == '$end':
            chunk = 'input2'
            continue     
        
        if tokens == ['z[m]', 'aw', 'qfld']:
            chunk = 'magnetic optics'
            print '      reading magnetic optics '
            continue

        if chunk == 'magnetic optics':
            z,aw,qfld = map(precision,tokens)
            out.z.append(z)
            out.aw.append(aw)
            out.qfld.append(qfld)
            
        if chunk == 'input1':
            tokens=line.replace('=','').strip().split()
            out.parameters[tokens[0]] = tokens[1:]
            #out.parameters[tokens[0]] = tokens[0:]
            #print 'input:', tokens
        if chunk == 'input2':
            tokens=line.replace('=','').strip().split()
            out.parameters['_'.join(tokens[1:])] = [tokens[0]]
            #out.parameters[tokens[0]] = tokens[0:]
            #print 'input:', tokens
#
        if chunk == 'slice' and readall:
            vals = map(precision,tokens)
            output_unsorted.append(vals)

        if chunk == 'slices':
            if len(tokens) == 2 and tokens[1]=='current':
                #print tokens[1]
                out.I.append(float(tokens[0]))
                out.n.append(nSlice)  
            elif len(tokens) == 3 and tokens[1]=='scan':
                out.I.append(float(tokens[0]))
                out.n.append(nSlice)  

    for parm in ['z', 'aw', 'qfld', 'I', 'n']:
        exec('out.'+parm+' = np.array(out.'+parm+')')
    
    if readall:
        output_unsorted=np.array(output_unsorted)#.astype(precision)
        # print out.sliceKeys
        for i in range(len(out.sliceKeys)):
            #exec('out.'+out.sliceKeys[i].replace('-','_').replace('<','').replace('>','') + ' = output_unsorted[:,'+str(i)+'].reshape(('+str(out.nSlices)+','+str(len(out.z))+'))')
            exec('out.'+out.sliceKeys[i].replace('-','_').replace('<','').replace('>','') + ' = output_unsorted[:,'+str(i)+'].reshape(('+str(out('history_records'))+','+str(out('entries_per_record'))+'))')
        if hasattr(out,'energy'):
            out.energy+=out('gamma0')
        out.power_z=np.max(out.power,0)

            
    out.nSlices = int(out('history_records'))#number of slices in the output
    out.nZ = int(out('entries_per_record'))#number of records along the undulator
    out.ncar=int(out('ncar')) #number of mesh points

    print '        nSlice', out.nSlices
    print '        nZ', out.nZ

    if out('itdp') == True:

        out.s = out('zsep') * out('xlamds') * np.arange(0,out.nSlices)
        out.t = out.s / c * 1.e+15
        #out.dt = (out.t[1] - out.t[0]) * 1.e-15
        out.dt=out('zsep') * out('xlamds') / c 
        out.beam_charge=np.sum(out.I*out('zsep')*out('xlamds')/c)
        if readall == True:
            out.spec = fft.fft(np.sqrt(np.array(out.power) ) * np.exp( 1.j* np.array(out.phi_mid) ) )
            out.freq_ev = h * fftfreq(len(out.spec), d=out('zsep') * out('xlamds') / c)# d=out.dt


    if out('dgrid')==0:
        rbeam=sqrt(out('rxbeam')**2+out('rybeam')**2)
        ray=sqrt(out('zrayl')*out('xlamds')/np.pi*(1+(out('zwaist')/out('zrayl')))**2);
        out.leng=2*out('rmax0')*(rbeam+ray)
    else:
        out.leng=2*out('dgrid')
        
    if out('iscan')!=0:
        out.scv=out.I #scan value
        out.I=np.linspace(1,1,len(out.scv)) #because used as a weight
    
    #tmp for back_compatibility    
    if readall:    
        out.power_int=out.power[:,-1]
        out.max_power=np.amax(out.power_int)
        for parm in [['power','p_int'],
                     ['energy','el_energy'],
                     ['e_spread','el_e_spread'],
                     ]:
             if hasattr(out,parm[0]):
                 setattr(out,parm[1],getattr(out,parm[0]))
    #             delattr(out,parm[0])
        out.power=out.p_mid[:,-1]
        out.phi=out.phi_mid[:,-1] 
    

    print('      done in %.3f seconds' % (time.time() - start_time))        
    return out


def dpa2dist (gen,file_name_read='',file_name_write='',no_macroparticles=1e5,debug=0):
    
    import random
    import numpy as np
    c = 299792458.0
        
    npart=int(gen('npart'))
    nslice=int(gen('nslice'))
    nbins=int(gen('nbins'))
    xlamds=gen('xlamds')
    zsep=int(gen('zsep'))
    gen_I=gen.I
    gen_t=gen.t    
    if file_name_read=='':
        file_name_read=gen.filename+'.dpa'
    if file_name_write=='':
        file_name_write=gen.filename+'.dist'
    par=read_particle_file(file_name_read, nbins=nbins, npart=npart,debug=debug)
    # print par.e.shape
    #start_time = time.time()
    #for i in range(100):
    m=np.arange(nslice)
    m=np.tile(m,(nbins,npart/nbins,1))
    m=np.rollaxis(m,2,0)
    # m=np.broadcast_to(m,[nbins,npart/nbins,nslice])
    # m=np.rollaxis(m,2,0)
    #print("--- Create matrix - %s seconds ---" % (time.time() - start_time))

    par.z=par.ph*xlamds/2/pi+m*xlamds*zsep+xlamds*zsep*(1-np.random.random((nslice, nbins,npart/nbins)))
    par.t=par.z/c
    
    t_scale=np.linspace(0,nslice*zsep*xlamds/c*1e15,nslice)
    # print 'range_t_scale', np.amin(t_scale), np.amax(t_scale)
    # print 'range_gen_t', np.amin(gen_t), np.amax(gen_t)
    I_scale=np.interp(t_scale,gen_t,gen_I)

    pick_n=I_scale
    pick_n=(pick_n/np.sum(pick_n)*no_macroparticles).astype(int)
    # print sum(pick_n)
    result_filesize=sum(pick_n)
    t_out=[]
    e_out=[]
    x_out=[]
    y_out=[]
    px_out=[]
    py_out=[]
    # print 'max_par.t', np.amax(par.t)
    # print 'par.e', np.amax(par.e),np.amin(par.e)
    par.t=np.reshape(par.t,(nslice,npart))
    par.e=np.reshape(par.e,(nslice,npart))
    par.x=np.reshape(par.x,(nslice,npart))
    par.y=np.reshape(par.y,(nslice,npart))
    par.px=np.reshape(par.px,(nslice,npart))
    par.py=np.reshape(par.py,(nslice,npart))
    # print par.t.shape
    
    # print par.t.shape
    # print 'max_par.t', np.amax(par.t)
    # print 'par.e', np.amax(par.e),np.amin(par.e)
    
    for i in arange(nslice):
        pick_i=random.sample(arange(nslice),pick_n[i])
        t_out=append(t_out,par.t[i,pick_i])
        e_out=append(e_out,par.e[i,pick_i])
        x_out=append(x_out,par.x[i,pick_i])
        y_out=append(y_out,par.y[i,pick_i])
        px_out=append(px_out,par.px[i,pick_i])
        py_out=append(py_out,par.py[i,pick_i])
    t_out=t_out*(-1)+max(t_out)
    # print 'max_y_out', np.amax(t_out)
    # print 'e_out', np.amax(e_out),np.amin(e_out)
    debug=0 #possible problems with pyplot on cluster
    if debug==1:
        import matplotlib.pyplot as plt
        bins=100
        # plt.figure('Time - Enenrgy')
        plt.figure(40001)
    #    plt.clf()
        plt.hist2d(t_out, e_out, bins)
        # plt.figure('Time - X')
        plt.figure(40002)
        plt.hist2d(t_out, x_out, bins)
        # plt.figure('X - Y')
        plt.figure(40003)
        plt.hist2d(x_out, y_out, bins)
        # plt.figure('X - pX')
        plt.figure(40004)
        plt.hist2d(x_out, px_out, bins)
        plt.show()

    #REQUIRES NUMPY 1.7
    # header='? VERSION = 1.0 \n? SIZE = %s \n? CHARGE = %E \n? COLUMNS X XPRIME Y YPRIME T P'%(result_filesize,gen.beam_charge)
    # np.savetxt(file_name_write, np.c_[x_out,px_out/e_out,y_out,py_out/e_out,t_out,e_out],header=header,fmt="%E", newline='\n',comments='')

    
    header='? VERSION = 1.0 \n? SIZE = %s \n? CHARGE = %E \n? COLUMNS X XPRIME Y YPRIME T P\n'%(result_filesize,gen.beam_charge)
    f = file(file_name_write,'w')
    f.write(header)
    f.close()
    f = file(file_name_write,'a')
    np.savetxt(f, np.c_[x_out,px_out/e_out,y_out,py_out/e_out,t_out,e_out],fmt="%E", newline='\n')
    f.close()



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

    #Next line added by GG 27.05.2016: it was in script
    
    #beam.emit_xn, beam.emit_yn = beam.emit[beam.C]
    #beam.gamma_rel = beam.E / (0.511e-3)
    #beam.emit_x = beam.emit_xn / beam.gamma_rel
    #beam.emit_y = beam.emit_yn / beam.gamma_rel
    


    inp.magin = 1
    #inp.zstop = 50
    #inp.nsec = 20
        
    inp.xlamd = up.lw
    inp.aw0 = up.K * np.sqrt(0.5) 
    inp.awd=inp.aw0
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
    printFelParameters(inp)

    inp.xlamds = felParameters.lambda0
    inp.prad0 = felParameters.power
    inp.fbess0 = felParameters.fc
    inp.zrayl = felParameters.zr
    
    if itdp:
        inp.type = "tdp"
        inp.DUMP_FIELDS = 1
        inp.ipseed = 132
        inp.ncar = 151
        #inp.nslice = 300
        inp.curlen = beam.tpulse * c/1e15
        inp.zsep = int(0.25/(4*pi*felParameters.rho)) #0.25 is the additional factor to be "on the safe side"
        # inp.zsep = 8 * int(inp.curlen  / inp.nslice / inp.xlamds )
        inp.nslice = 8 * int(inp.curlen  / inp.zsep / inp.xlamds )
    
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

def rad_file_str(beam):
    #header = "# \n? VERSION = 1.0\n? SIZE ="+str(len(beam.column_values['ZPOS']))+"\n? COLUMNS ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY XBEAM YBEAM PXBEAM PYBEAM ALPHAX ALPHAY CURPEAK ELOSS\n"
    header = "# \n? VERSION = 1.0\n? SIZE ="+str(len(beam.z))+"\n? COLUMNS"
    #ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY XBEAM YBEAM PXBEAM PYBEAM ALPHAX ALPHAY CURPEAK ELOSS\n"
    for c in beam.columns:
        header = header + " " + c 
    header +="\n"
    
    f_str = header
    
    
    beam.column_values['ZPOS'] = beam.z 
    beam.column_values['PRAD0'] = beam.prad0
    
    
    for i in xrange(len(beam.z)):
        for c in beam.columns:
            buf = str(beam.column_values[c][i])
            f_str = f_str + buf + ' '
        f_str = f_str.rstrip() +  '\n'
	
    print beam.z[301]
    	
    
    return f_str
    
    
    
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

def adapt_rad_file(beam = None, rad_file = None, out_file='tmp.rad'):
        
    rad = read_rad_file(rad_file)    
    rad.prad0 = np.interp(beam.z, rad.z, rad.prad0)
    rad.z     = beam.z
    
    #print rad.z[301]    
    #print beam.z[301]
    open(out_file,'w').write(rad_file_str(rad))
    #exit()
    
    
    

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
  

def cut_beam(beam = None, cut_z = [-inf, inf]):
    if np.amin(beam.z)<cut_z[0] or np.amax(beam.z)>cut_z[1]:
        
        condition = (beam.z > cut_z[0]) * (beam.z<cut_z[1])
        beam_new = Beam()
        beam_new.column_values = beam.column_values
        beam_new.columns = beam.columns
        
        for parm in ['x','px','y','py','z','I','ex','ey','g0','dg','eloss','betax','betay','alphax','alphay']:
            if hasattr(beam,parm):
                setattr(beam_new,parm,np.extract(condition,getattr(beam,parm)))
        
        # beam_new.x = np.extract(condition,beam.x)
        # beam_new.px = np.extract(condition,beam.px)
        # beam_new.y = np.extract(condition,beam.y)
        # beam_new.py = np.extract(condition,beam.py)
        # beam_new.z = np.extract(condition,beam.z)
        
        # beam_new.I = np.extract(condition,beam.I)
        # beam_new.ex = np.extract(condition,beam.ex)
        # beam_new.ey = np.extract(condition,beam.ey)
        # beam_new.g0 = np.extract(condition,beam.g0)
        # beam_new.dg = np.extract(condition,beam.dg)
        # beam_new.eloss = np.extract(condition,beam.eloss)
        
        # beam_new.betax = np.extract(condition,beam.betax)
        # beam_new.betay = np.extract(condition,beam.betay)
        # beam_new.alphax = np.extract(condition,beam.alphax)
        # beam_new.alphay = np.extract(condition,beam.alphay)

        beam_new.zsep = beam.zsep
        zmax, Imax = peaks(beam_new.z, beam_new.I, n=1)
        beam_new.idx_max = np.where(beam_new.z == zmax)[0][0]
    else:
        beam_new=beam
    return beam_new

def get_beam_peak(beam = None): #experimental, the code is too inconsistent to introduce suc function yet (e.g. xp <-> px)
    import copy
    #obtains the peak current values
    if len(beam.I)>1:# and np.amax(beam.I)!=np.amin(beam.I):
        pkslice = np.argmax(beam.I)
        
        beam_new=copy.deepcopy(beam)
        
        beam_new.I=beam.I[pkslice]
        beam_new.alpha_x=beam.alphax[pkslice]
        beam_new.alpha_y=beam.alphay[pkslice]
        beam_new.beta_x=beam.betax[pkslice]
        beam_new.beta_y=beam.betay[pkslice]
        beam_new.emit_xn=beam.ex[pkslice]
        beam_new.emit_yn=beam.ey[pkslice]
        beam_new.gamma_rel=beam.g0[pkslice]
        beam_new.sigma_E=beam.dg[pkslice]*(0.000510998)
        
        beam_new.E=beam_new.gamma_rel*(0.511e-3)
        beam_new.emit_x = beam_new.emit_xn / beam_new.gamma_rel
        beam_new.emit_y = beam_new.emit_yn / beam_new.gamma_rel
        
        for parm in ['alphax','alphay','betax','betay','z','ex','ey','g0','dg']:
            if hasattr(beam_new,parm):
                delattr(beam_new,parm)
        
        # beam_new.x = np.extract(condition,beam.x)

    else:
        beam_new=beam
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

