
'''
interface to genesis
'''

import struct
from copy import copy, deepcopy
import time
import os
import socket
import errno
from ocelot.rad.fel import *
from ocelot.cpbd.beam import * # Twiss, Beam, gauss_from_twiss, ParticleArray
from ocelot.cpbd.elements import *
import ocelot.utils.reswake as w
from ocelot.utils.launcher import *
from ocelot.common.math_op import *
from ocelot.common.globals import *  # import of constants like "h_eV_s" and "speed_of_light"

import math
import numpy as np
from numpy import mean, std, inf, shape, append, complex128, complex64

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
 iotail = __IOTAIL__\n\
 nharm = __NHARM__\n\
 iharmsc = __IHARMSC__\n\
 iallharm = __IALLHARM__\n\
 curpeak =  __CURPEAK__\n\
 curlen =  __CURLEN__\n\
 ntail = __NTAIL__\n\
 nslice = __NSLICE__\n\
 shotnoise = __SHOTNOISE__\n\
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
 trama = __TRAMA__\n\
 itram11 = __ITRAM11__\n\
 itram12 = __ITRAM12__\n\
 itram13 = __ITRAM13__\n\
 itram14 = __ITRAM14__\n\
 itram15 = __ITRAM15__\n\
 itram16 = __ITRAM16__\n\
 itram21 = __ITRAM21__\n\
 itram22 = __ITRAM22__\n\
 itram23 = __ITRAM23__\n\
 itram24 = __ITRAM24__\n\
 itram25 = __ITRAM25__\n\
 itram26 = __ITRAM26__\n\
 itram31 = __ITRAM31__\n\
 itram32 = __ITRAM32__\n\
 itram33 = __ITRAM33__\n\
 itram34 = __ITRAM34__\n\
 itram35 = __ITRAM35__\n\
 itram36 = __ITRAM36__\n\
 itram41 = __ITRAM41__\n\
 itram42 = __ITRAM42__\n\
 itram43 = __ITRAM43__\n\
 itram44 = __ITRAM44__\n\
 itram45 = __ITRAM45__\n\
 itram46 = __ITRAM46__\n\
 itram51 = __ITRAM51__\n\
 itram52 = __ITRAM52__\n\
 itram53 = __ITRAM53__\n\
 itram54 = __ITRAM54__\n\
 itram55 = __ITRAM55__\n\
 itram56 = __ITRAM56__\n\
 itram61 = __ITRAM61__\n\
 itram62 = __ITRAM62__\n\
 itram63 = __ITRAM63__\n\
 itram64 = __ITRAM64__\n\
 itram65 = __ITRAM65__\n\
 itram66 = __ITRAM66__\n\
__OUTPUTFILE__\n\
__BEAMFILE__\n\
__PARTFILE__\n\
__FIELDFILE__\n\
__RADFILE__\n\
__DISTFILE__\n\
__MAGFILE__\n\
 filetype ='ORIGINAL'\n\
 $end\n"

# outputfile ='run.__RUNID__.gout'\n\
# iallharm =  __IALLHARM__\n\
# iharmsc =  __IHARMSC__\n\
# pradh0 =  __PRADH0__\n\


class GenesisInput:  
    '''
    Genesis input files storage object
    '''

    def __init__(self):

        # defaults
        self.stageid = None  # optional, handy with multi-stage scripts
        self.runid = 0  # important for statistical runs
        self.type = 'steady'

        # undulator
        self.aw0 = 0.735  # The normalized, dimensionless rms undulator parameter, defined by AW0 = (e/mc)(Bu/ku), where e is the electron charge, m is electron mass, c is speed of light, ku=2pi/lambdau is the undulator wave number, lambdau is the undulator period. Bu is the rms undulator field with Bu = Bp/2 for a planar undulator and Bu = Bp for a helical undulator, where Bp is the on-axis peak field.
        self.awd = 0.735  # A virtual undulator parameter for the gap between undulator modules.
        self.wcoefz = [0, 0, 0]  # (1-[m]) Start of undulator tapering.  Note that tapering is applied, even the magnetic lattice is defined by an external file.
        #(2-[ ]) The relative change of the undulator field over the entire taper length (AW(exit) = (1 -WCOEFZ(2))
        #(3) The taper model: 1 for linear taper, 2 for quadratic taper,
        self.iertyp = 0  # Type of undulator field errors.
        self.iwityp = 0  # the undulator type. A value of zero indicates a planar undulator, any other value a helical one.
        self.xkx = 0  # Normalized natural focusing of the undulator in x. Common values are XKX = 0.0, XKY = 1.0 for a planar undulator or XKX, XKY = 0.5 for a helical undulator, but might vary if focusing by curved pole faces is simulated. The values should fulfill the constraint XKX + XKY = 1.0.
        self.xky = 1  # Normalized natural focusing of the undulator in y
        self.delaw = 0  # RMS value of the undulator field error distribution. A value of zero disables field errors.
        self.nwig = 98   # The number of periods within a single undulator module. The product of NWIG and XLAMD defines the length of the undulator module.
        self.nsec = 1   # The number of sections of the undulator.
        self.awx = 0  # Maximum offset in x for undulator module misalignment. The error for each individual module follows a uniform distribution
        self.awy = 0  # Maximum offset in y for undulator module misalignment. The error for each individual module follows a uniform distribution
        self.iseed = -1  # initial seeding of the random number generator for field errors

        # electron beam
        self.curpeak = 2.5E+03  # Peak current of the electron beam. Time-independent simulations enforces a constant current.
        self.curlen = 7E-06  # Bunch length of the current profile. If CURLEN is positive a Gaussian distribution is generated with an rms length given by CURLEN. A negative or zero value yield a constant profile.
        self.npart = 8192   # number of macroparticles per slice. NPART must be a multiple of 4*NBINS
        self.gamma0 = 3.424658E+04   # The mean energy of the electron beam in terms of the electron rest mass energy.
        self.delgam = 5.0E-3  # The RMS value of the energy distribution in terms of electron rest mass energy.
        self.rxbeam = 100e-6  # The rms value in x of the transverse, spatial distribution.
        self.rybeam = 100e-6  # The rms value in y of the transverse, spatial distribution.
        self.emitx = 1.0e-7  # The normalized rms emittance in x
        self.emity = 1.0e-7  # The normalized rms emittance in y
        self.alphax = 0.0  # Rotation of the transverse phase space distribution in x according to the standard definition ALPHAX = - < xx' > GAMMA0 / EMITX.
        self.alphay = 0.0  # Rotation of the transverse phase space distribution in y according to the standard definition ALPHAY = - < yy' > GAMMA0 / EMITY.
        self.xbeam = 0.0  # Transverse position in x of the electron beam at the undulator entrance with respect to the undulator axis
        self.ybeam = 0.0  # Transverse position in y of the electron beam at the undulator entrance with respect to the undulator axis
        self.pxbeam = 0.0  # Average normalized transverse momentum x of the electron beam at the undulator entrance. The momenta are defined as PXBEAM = betax*gamma where betax = c*v_x is the average transverse velocity and gamma the Lorenz factor of the electron energy.
        self.pybeam = 0.0  # y
        self.cuttail = -1.0  # Cut in the transverse phase space in measures of the rms size to collimate transverse beam tails/halos. The cut is applied after the loading and beam current is set accordingly to the reduced number of macro particles. It is disabled if the value is negative or the electron beam is imported from an external file.
        self.conditx = 0.0  # the correlation strength between the amplitude of the eletron's betatron oscillation x and its energy.
        self.condity = 0.0  # y
        self.bunch = 0.0  # Initial value of the bunching factor, when quite loading is used.
        self.bunchphase = 0.0  # Phase of initial bunching, when quite loading is used.
        self.emod = 0.0  # Initial energy modulation, when quite loading is used. The value is the modulation amplitude in terms of gamma.
        self.emodphase = 0.0  # Phase of initial energy modulation, when quite loading is used.

        # particle loading
        self.ildpsi = 7  # Indices of the Hammersley sequence bases for loading the particle phase.
        self.ildgam = 5  # Hammersley base for loading the energy distribution.
        self.ildx = 1  # Hammersley base for loading the distribution in x.
        self.ildy = 2  # Hammersley base for loading the distribution in y.
        self.ildpx = 3  # Hammersley base for loading the distribution in px.
        self.ildpy = 4  # Hammersley base for loading the distribution in py.
        self.itgaus = 1  # Defines distribution profile in the transverse variables. The available distributions are: Gaussian (1) Uniform (2) Parabolic (otherwise)
        self.igamgaus = 1  # Defines distribution profile in energy. A non-zero value generates a Gaussian distribution and a uniform otherwise.
        self.iall = 0  # A non-zero value of IALL enforces that all slices are starting with the same element of the Hammersley sequences.
        self.ipseed = -1  # Initial seed for the random number generator used for particle phase fluctuation (shot noise). GENESIS 1.3 requires a re-initialization of the random number generator to guarantee the same loading whether magnetic field errors are used or not.
        self.nbins = 4  # Number of bins in the particle phase. The value has to be at least 4 or larger than (2+2n), depending on whether the bunching at the nth harmonics is needed for space charge calculations or output.

        # radiation
        self.xlamds = 1E-9  # The resonant radiation wavelength.
        self.prad0 = 0  # The input power of the radiation field.
        self.pradh0 = 0  # Radiation power for seeding with a harmonics, defined by NHARM.
        self.zrayl = 0.5  # The Rayleigh length of the seeding radiation field.
        self.zwaist = 2  # Position of the waist of the seeding radiation field with respect to the undulator entrance.

        # mesh
        self.ncar = 151   # The number of grid points for the radiation field along a single axis. The total number for the mesh is NCAR^2
        self.lbc = 0  # Flag to set the boundary condition for the field solver. The options are Direchlet boundary condition (LBC = 0) and Neumann boundary condition (otherwise).
        self.nscr = 0  # Number of azimuthal modes for space charge calculation.
        self.nscz = 0  # Number of longitudinal Fourier components for space charge calculation. NSCZ must be greater than 0 to include space charge but less than (2NBINS+2) for valid results.
        self.nptr = 40   # Number of radial grid points, on which the space charge field is evaluated.
        self.rmax0 = 9.0   # mesh size in units of radiation+beam sizes
        self.dgrid = 0.0   # explicit trnsverse mesh size overruling the calculation by the RMAX0-parameter.

        # focusing
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

        # simulation
        self.version = 0.1  # Used for backward compatibility of the input decks.
        self.zsep = 20  # Separation of beam slices in measures of the radiation wavelength. ZSEP must be a multiple of DELZ.
        self.nslice = 0   # Total number of simulated slices. It defines the time window of the simulation with NSLICE * ZSEP * XLAMDS/c
        self.ntail = - self.nslice / 2  # Position of the first simulated slice in measures of ZSEP*XLAMDS. GENESIS 1.3 starts with the tail side of the time window, progressing towards the head. Thus a negative or positive value shifts the slices towards the tail or head region of the beam, respectively.
        self.delz = 1.0   # Integration step size in measure of the undulator period length.
        self.zstop = 256.0  # Defines the total integration length. If the undulator length is shorter than ZSTOP or ZSTOP is zero or negative, the parameter is ignored and the integration is performed over the entire undulator.
        self.iorb = 0   # enforce orbit correction. For any non-zero value the offset due to the wiggle motion is taken into account for the interaction between electron beam and radiation field.
        self.isravg = 1  # If set to a non-zero value the energy loss due to spontaneous synchrotron radiation is included in the calculation.
        self.isrsig = 1  # If set to a non-zero value the increase of the energy spread due to the quantum fluctuation of the spontaneous synchrotron radiation is included in the calculation.
        self.eloss = 0  # Externally applied energy loss of the electron beam.
        self.nharm = 1  # Enables the calculation of harmonics up to the one, specified by NHARM. Note that the number of NBINS has to be at least twice as large as NHARM to allow for the correct representation of the harmonics. Note also that this parameter does not enable automatically the output. For that the corresponding bit in LOUT has to be set as well.
        self.iscan = 0  # Selects the parameter for a scan over a certain range of its value
        #(1.GAMMA0 2.DELGAM 3.CURPEAK 4.XLAMDS 5.AW0 6.ISEED 7.PXBEAM 8.PYBEAM 9.XBEAM 10.YBEAM 11.RXBEAM 12.RYBEAM 13.XLAMD 14.DELAW 15.ALPHAX 16.ALPHAY 17.EMITX 18.EMITY 19.PRAD0 20.ZRAYL 21.ZWAIST 22.AWD 23.BEAMFILE 24.BEAMOPT 25.BEAMGAM)
        # self.scan = '' #By supplying the parameter name to scan over it overrules the setting of ISCAN
        self.nscan = 3  # Number of steps per scan.
        self.svar = 0.01  # Defines the scan range of the selected scan parameter. The parameter is varied between (1-SVAR) and (1+SVAR) of its initial value.

        # I/O
        self.iphsty = 1  # Generate output in the main output file at each IPHSTYth integration step. To disable output set IPHSTY to zero.
        self.ishsty = 1  # Generate output in the main output file for each ISHSTYth slice.
        self.ippart = 0  # Write the particle distribution to file at each IPPARTth integration step. To disable output set IPPART to zero. The filename is the same of the main outputfile + the extension '.par'.
        self.ispart = 0  # Write the particle distribution to file for every ISPART slice.
        self.ipradi = 0  # Write the field distribution to file at each IPRADIth integration step. To disable output set IPRADI to zero. The filename is the same of the main outputfile + the extension '.fld'.
        self.isradi = 0  # Write the field distribution to file for every ISRADI slice.
        self.iotail = 1  # If set to a non-zero value the output time window is the same as the simulated time window. Otherwise the output for the first slices covered by the slippage length is subpressed.
        self.magin = 1   # read in magnetic lattice (If set to a non-zero value the user is prompted to type in the file name containing a explicit description of the magnetic field.)
        self.magout = 0   # output magnetic lattice
        self.idump = 0  # If set to a non-zero value the complete particle and field distribution is dumped at the undulator exit into two outputfiles.
        self.idmpfld = 0  # Similar to IDUMP but only for the field distribution.
        self.idmppar = 0  # Similar to IDUMP but only for the particle distribution.
        self.lout = [1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        #            1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
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

        self.ndcut = -1  # ??? If NDCUT is zero, the time-window is adjusted, so that in average NPART/NBINS particles fall in each slice.
        self.alignradf = 1  # if zero , Genesis 1.3 aligns the radiation field to the electron beam so that the radiaiton field is one ful slippage behind the electron beam.
        self.offsetradf = 0  # slices to shift the electrron beam with respect to the radiation if ALIGNRADF=1.
        self.convharm = 1  # When the particle distribution is imported from a PARTFILE Genesis 1.3 allows the upconversion to a higher harmonics. The harmonic number is specified with CONVHARM and has a defulat value of 1, corresponding to no up-conversion. The user has to make sure that in the input deck XLAMDS is adjusted, according to the new wavelength.
        self.multconv = 0  # If an imported particle distribution from a PARTFILE is up-converted to a higher harmonics the dault behavior is that the number of slices is preserved. This requires that ZSEPis adjusted together with XLAMDS. However if frequency resolution is a problem then a particle distribution can be converted and used multiple times to keep ZSEP constant. The disadvantage is that the CPU execution time is increased as well.

        self.ibfield = 0.0  # When the PARTFILE features is used the imported particle distribution can be tracked through a generic 4 magnet chicane before running the Genesis simulation. The chicane consists out of 4 bending magnets of the field strength IBFIELD and length IMAGL separated by 5 drifts of length IDRIL.
        self.imagl = 0.0  # The length of each bending magnet of the chicane. If the magnet length is set to zero but IDRIL is not the resulting beam line correspond to a simple drift of the length 5 times IDRIL
        self.idril = 0.0  # The length of the 5 drift lengths of the magnetic chicane (three between the magnets and one before and after the magnets).

        self.trama = 0  # Non zero value enables that a transport matrix is applied to the electron distribution when importing it with PARTFILE. The individual matrix is defined by ITRAM$$
        self.itram11 = 1  # The pound signs are place holders for numbers between 1 and 6 (e.g. ITRAM21) and are defining the matrix element for the transport matrix, which is applied when importing a paticle distribution with the PARTFILE option. The matrix is defined in a standard way, acting on the vector (position in X, angle in X, position in Y, angle in Y, position in s, relative energy spread). The default value is the identity matrix.
        self.itram12 = 0
        self.itram13 = 0
        self.itram14 = 0
        self.itram15 = 0
        self.itram16 = 0
        self.itram21 = 0
        self.itram22 = 1
        self.itram23 = 0
        self.itram24 = 0
        self.itram25 = 0
        self.itram26 = 0
        self.itram31 = 0
        self.itram32 = 0
        self.itram33 = 1
        self.itram34 = 0
        self.itram35 = 0
        self.itram36 = 0
        self.itram41 = 0
        self.itram42 = 0
        self.itram43 = 0
        self.itram44 = 1
        self.itram45 = 0
        self.itram46 = 0
        self.itram51 = 0
        self.itram52 = 0
        self.itram53 = 0
        self.itram54 = 0
        self.itram55 = 1
        self.itram56 = 0
        self.itram61 = 0
        self.itram62 = 0
        self.itram63 = 0
        self.itram64 = 0
        self.itram65 = 0
        self.itram66 = 1

        self.iallharm = 0  # Setting the value to a non-zero value will also include all harmonics between 1 and NHARM
        self.iharmsc = 0  # setting to a non-zero value includes the coupling of the harmonic radiation back to the electron beam for a self-consistent model of harmonics. Enabling this feature will automatically include all harmonics by setting IALLHARM to one.
        self.isntyp = 0  # Non-zero if the user wants to use the Pennman algorithm for the shot noise (which is not recommended).

        self.ilog = 0  # Create a log file.
        self.ffspec = 0  # amplitude/phase values for spectrum calculation: 0 - on-axis power/phase along the pulse, -1 - the same in far field, 1 - near field total power
        
        self.shotnoise = 1
        
        # paths to files to import
        self.beamfile = None
        self.fieldfile = None
        self.partfile = None
        self.edistfile = None
        self.outputfile = None
        self.radfile = None
        self.latticefile = None

        # objects to write into e.g. *.dfl or *.dpa with appropriate names
        # and imported to Genesis
        self.beam = None # GenesisBeam()
        self.edist = None # GenesisElectronDist()
        self.lat = None # MagneticLattice()
        self.dfl = None # RadiationField()
        self.dpa = None # GenesisParticlesDump()
        self.rad = None # GenesisRad()

        self.run_dir = None # directory to run simulation in
        self.exp_dir = None # if run_dir==None, it is created based on exp_dir

    def input(self):
        input = inputTemplate

        if self.type == 'steady':
            # input = input.replace("__SHOTNOISE__", "itdp  =    0")
            input = input.replace("__ITDP__", "itdp = 0")
        else:
            # input = input.replace("__SHOTNOISE__", "shotnoise=  1")
            input = input.replace("__ITDP__", "itdp = 1")
            # self.prad0 = 0

        if self.beamfile != None:
            input = input.replace("__BEAMFILE__", " beamfile  =  '" + str(self.beamfile) + "'")
        else:
            input = input.replace("__BEAMFILE__\n", "")

        if self.fieldfile != None:
            input = input.replace("__FIELDFILE__", " fieldfile  =  '" + str(self.fieldfile) + "'")
        else:
            input = input.replace("__FIELDFILE__\n", "")

        if self.partfile != None:
            input = input.replace("__PARTFILE__", " partfile  =  '" + str(self.partfile) + "'")
        else:
            input = input.replace("__PARTFILE__\n", "")

        if self.edistfile != None:
            input = input.replace("__DISTFILE__", " distfile  =  '" + str(self.edistfile) + "'")
        else:
            input = input.replace("__DISTFILE__\n", "")

        if self.outputfile != None:
            input = input.replace("__OUTPUTFILE__", " outputfile  =  '" + str(self.outputfile) + "'")
        else:
            input = input.replace("__OUTPUTFILE__", " outputfile ='run.__RUNID__.gout'")

        # print 'self.radfile is equal to ', self.radfile
        if self.radfile != None:
            input = input.replace("__RADFILE__", " radfile  =  '" + str(self.radfile) + "'")
        else:
            input = input.replace("__RADFILE__\n", "")

        if self.magin == 0:
            input = input.replace("__MAGFILE__\n", "")
        else:
            input = input.replace("__MAGFILE__", " maginfile ='" + str(self.latticefile) + "'")

        # if self.trama == 1:
            # input = input.replace("__TRAMA__\n", "")
        # else:
            # input = input.replace("__TRAMA__\n", "")

        for p in self.__dict__.keys():
            input = input.replace("__" + str(p).upper() + "__", str(self.__dict__[p]).replace('[', '').replace(']', '').replace(',', ''))

        return input

    def __getattr__(self, name):
        if name not in self.__dict__.keys():
            return 0.0
        else:
            return self.__dict__[name]


class GenesisOutput:
    '''
    Genesis output *.out files storage object
    '''

    def __init__(self):
        self.z = []
        self.s = []
        self.I = []
        self.n = []
        self.zSlice = []
        self.E = []
        self.aw = []
        self.qfld = []

        self.sliceKeys = []
        self.sliceValues = {}

        self.parameters = {}
        self.filePath = ''

    def fileName(self):
        return filename_from_path(self.filePath)

    def __call__(self, name):

        if name not in self.parameters.keys():
            return None
        else:
            p, = self.parameters[name]
            return float(p.replace('D', 'E'))

    def calc_spec(self,npad=1):
        if npad <= 1:
            return
        p_int = self.p_int
        phi_mid = self.phi_mid

        zeros = np.zeros((self.nSlices*npad,self.nZ))
        p_int = np.vstack((p_int, zeros))
        phi_mid = np.vstack((phi_mid, zeros))

        spec = abs(np.fft.fft(np.sqrt(np.array(p_int)) * np.exp(1.j * np.array(phi_mid)), axis=0))**2 / (sqrt(self.nSlices) * (2 * self.leng / self('ncar'))**2 * 1e10)
        e_0 = h_eV_s * speed_of_light / self('xlamds')
        freq_ev = h_eV_s * np.fft.fftfreq(len(spec), d=self('zsep') * self('xlamds') * self('ishsty') / speed_of_light) + e_0

        spec = np.fft.fftshift(spec, axes=0)
        freq_ev = np.fft.fftshift(freq_ev, axes=0)
        freq_lamd = h_eV_s * speed_of_light * 1e9 / freq_ev

        self.spec = spec
        self.freq_ev = freq_ev
        self.freq_lamd = freq_lamd

            
    def wig(self,z=inf):
        return wigner_out(self, z=z, method='mp', debug=1)


class GenStatOutput:
    '''
    Genesis statistical output storage object
    '''
    def __init__(self):
        return

    def zi(self, z):
        if z > np.amax(self.z):
            z = np.amax(self.z)
        return np.where(self.z >= z)[0][0]

    def si(self, s):
        if s > np.amax(self.s):
            s = np.amax(self.s)
        return np.where(self.s >= s)[0][0]

    def fi(self, f):
        if f > np.amax(self.f):
            f = np.amax(self.f)
        return np.where(self.f >= f)[0][0]


class GenesisParticlesDump:
    '''
    Genesis particle *.dpa files storage object
    Each particle record in z starts with the energy of all particles 
    followed by the output of the particle phases, 
    positions in x and y and the momenta in x and y. 
    The momenta are normalized to mc
    '''

    def __init__(self):
        self.e = []
        self.ph = []
        self.x = []
        self.y = []
        self.px = []
        self.py = []

        # self.fileName = ''
        self.filePath = ''

    def fileName(self):
        return filename_from_path(self.filePath)


class GenesisElectronDist:
    '''
    Genesis electron beam distribution *.dat files storage object
    GENESIS follows the rule that the longitudinal position 
    is reversed if a time is specified by the T column. 
    In this case smaller numbers correspond to particle 
    in the tail of the distribution.
    '''

    def __init__(self):

        self.x = []  # position in x in meters
        self.y = []  # position in y in meters
        self.xp = []  # divergence in x ### (xprime == angle)
        self.yp = []  # divergence in y
        self.t = []  # longitudinal position in seconds
        self.g = []  # gamma (total energy, normalized mc2) #rename to g?
        self.part_charge = []  # charge per particle
        self.filePath = ''

    def charge(self):  # full charge
        return self.part_charge * self.len()

    def fileName(self):
        return filename_from_path(self.filePath)

    def len(self):
        return len(self.t)

    def center(self):
        edist_out = deepcopy(self)
        edist_out.x -= np.mean(edist_out.x)
        edist_out.y -= np.mean(edist_out.y)
        edist_out.xp -= np.mean(edist_out.xp)
        edist_out.yp -= np.mean(edist_out.yp)
        return edist_out

    def twiss(self):
        tws = Twiss()

        x = self.x
        y = self.y
        xp = self.xp
        yp = self.yp

        mean_x2 = mean(x**2)
        mean_y2 = mean(y**2)
        mean_px2 = mean(xp**2)
        mean_py2 = mean(yp**2)
        mean_xpx = mean(x * xp)
        mean_ypy = mean(y * yp)
        mean_g = mean(self.g)

        tws.emit_x = mean_g * (mean_x2 * mean_px2 - mean_xpx**2)**0.5 / mean_g
        tws.emit_y = mean_g * (mean_y2 * mean_py2 - mean_ypy**2)**0.5 / mean_g
        tws.beta_x = mean_g * mean_x2 / tws.emit_x
        tws.beta_y = mean_g * mean_y2 / tws.emit_y
        tws.alpha_x = -mean_g * mean_xpx / tws.emit_x
        tws.alpha_y = -mean_g * mean_ypy / tws.emit_y
        tws.E = mean_g * m_e_GeV

        return tws


def parray2edist(p_array):
    
    edist = GenesisElectronDist()
    
    e0 = p_array.E * 1e9 #[eV]
    p0 = sqrt( (e0**2 - m_e_eV**2) / speed_of_light**2 )
    
    p_oc = p_array.particles[5::6] # deltaE / average_impulse / speed_of_light
    edist.g = (p_oc * p0 * speed_of_light + e0) / m_e_eV
    edist.x = p_array.particles[::6]  # position in x in meters
    edist.y = p_array.particles[2::6]  # position in y in meters
    edist.xp = p_array.particles[1::6]  # divergence in x
    edist.yp = p_array.particles[3::6]  # divergence in y
    edist.t = p_array.particles[4::6] / speed_of_light  # longitudinal position in seconds

    edist.part_charge = p_array.q_array[0] #fix for general case  # charge per particle
    edist.filePath = ''
    
    return edist
    
def edist2parray(edist):

    p_array = ParticleArray()
    p_array.particles = np.zeros(edist.len() * 6)
    p_array.q_array = np.ones(edist.len()) * edist.part_charge
    
    g0 = np.mean(edist.g) # average gamma
    e0 = g0 * m_e_eV
    p0 = sqrt(g0**2-1) * m_e_eV / speed_of_light # average impulse
#    p0 = sqrt( (e0**2 - m_e_eV**2) / speed_of_light**2 ) # average impulse
    p_array.E = g0 * m_e_GeV # average energy in GeV
    
    p_array.particles[::6] = edist.x # position in x in meters
    p_array.particles[1::6] = edist.xp  # divergence in x
    p_array.particles[2::6] = edist.y # position in x in meters
    p_array.particles[3::6] = edist.yp  # divergence in x
    p_array.particles[4::6] = edist.t * speed_of_light
    p_array.particles[5::6] = (edist.g - g0) * m_e_eV / p0 / speed_of_light
    
    return p_array
    
    
        
        # def twiss(self):#not tested!!!
        # from ocelot.cpbd.beam import Twiss
        # tws=Twiss()
        # tws.x=mean(self.x)
        # tws.y=mean(self.y)
        # tws.px=mean(self.px)
        # tws.py=mean(self.py)
        # self=self.center()

        # self.px = self.px*(1.-0.5*self.px**2 - 0.5*self.py**2)
        # self.py = self.py*(1.-0.5*self.px**2 - 0.5*self.py**2)

        # tws.xx=mean(self.x**2)
        # tws.yy=mean(self.y**2)
        # tws.pxpx=mean(self.px**2)
        # tws.pypy=mean(self.py**2)
        # tws.xpx=mean(self.x*self.px)
        # tws.ypy=mean(self.y*self.py)
        # # mean_g=mean(self.e)
        # tws.E=mean(self.e)*m_e_GeV

        # tws.emit_x= sqrt(tws.x*tws.pxpx-tws.xpx**2)
        # tws.emit_y= sqrt(tws.y*tws.pypy-tws.ypy**2)
        # tws.beta_x=tws.x/tws.emit_x
        # tws.beta_y=tws.y/tws.emit_y
        # tws.alpha_x=-tws.xpx/tws.emit_x
        # tws.alpha_y=-tws.ypy/tws.emit_y

        # return tws


class GenesisBeam():
    '''
    Genesis analytical radiation input files storage object?
    '''

    def __init__(self):
        self.columns = []
        self.column_values = {}
        self.fileName = ''
        # self.filePath=''

    def fileName(self):
        return filename_from_path(self.filePath)

    def len(self):
        return len(self.z)

    def idx_max_refresh(self):
        self.idx_max = np.argmax(self.I)

    def __delitem__(self, indarr):
        self.z = np.delete(self.z, indarr)
        if hasattr(self, 'I'):
            self.I = np.delete(self.I, indarr)
        if hasattr(self, 'g0'):
            self.g0 = np.delete(self.g0, indarr)
        if hasattr(self, 'dg'):
            self.dg = np.delete(self.dg, indarr)
        if hasattr(self, 'x'):
            self.x = np.delete(self.x, indarr)
        if hasattr(self, 'y'):
            self.y = np.delete(self.y, indarr)
        if hasattr(self, 'px'):
            self.px = np.delete(self.px, indarr)
        if hasattr(self, 'py'):
            self.py = np.delete(self.py, indarr)
        if hasattr(self, 'ex'):
            self.ex = np.delete(self.ex, indarr)
        if hasattr(self, 'ey'):
            self.ey = np.delete(self.ey, indarr)
        if hasattr(self, 'betax'):
            self.betax = np.delete(self.betax, indarr)
        if hasattr(self, 'betay'):
            self.betay = np.delete(self.betay, indarr)
        if hasattr(self, 'alphax'):
            self.alphax = np.delete(self.alphax, indarr)
        if hasattr(self, 'alphay'):
            self.alphay = np.delete(self.alphay, indarr)
        return self


class GenesisRad():
    '''
    Genesis analytical radiation input files storage object?
    '''

    def __init__(self):
        self.columns = []
        self.column_values = {}


class RadiationField():
    '''
    3d or 2d coherent radiation distribution, *.fld variable is the same as Genesis dfl structure
    '''

    def __init__(self, shape=(0, 0, 0)):
        # self.fld=np.array([]) #(z,y,x)
        self.fld = np.zeros(shape, dtype=complex128)  # (z,y,x)
        self.dx = []
        self.dy = []
        self.dz = []
        self.xlamds = 0  # wavelength, [nm]
        self.domain_z = 't'  # longitudinal domain (t - time, f - frequency)
        self.domain_xy = 's'  # transverse domain (s - space, k - inverse space)
        self.filePath = ''

    def fileName(self):
        return filename_from_path(self.filePath)

    def copy_param(self, dfl1):
        self.dx = dfl1.dx
        self.dy = dfl1.dy
        self.dz = dfl1.dz
        self.xlamds = dfl1.xlamds
        self.domain_z = dfl1.domain_z
        self.domain_xy = dfl1.domain_xy
        self.filePath = dfl1.filePath

    def __getitem__(self, i):
        return self.fld[i]

    def __setitem__(self, i, fld):
        self.fld[i] = fld

    def shape(self):
        return shape(self.fld)

    def Lz(self):  # full transverse mesh size, 2*dgrid
        return self.dz * self.Nz()

    def Ly(self):  # full transverse mesh size, 2*dgrid
        return self.dy * self.Ny()


    def Lx(self):  # full longitudinal mesh size, nslice*zsep*xlamds
        return self.dx * self.Nx()

    def Nz(self):
        return shape(self.fld)[0]

    def Ny(self):
        return shape(self.fld)[1]

    def Nx(self):
        return shape(self.fld)[2]

    def int(self):  # 3d intensity
        return self.fld.real**2 + self.fld.imag**2

    def int_z(self):  # intensity projection on z (power [W] or spectral density)
        return np.sum(self.int(), axis=(1, 2))

    def ang_z_onaxis(self):
        xn = int((self.Nx() + 1) / 2)
        yn = int((self.Ny() + 1) / 2)
        fld = self[:, yn, xn]
        return np.angle(fld)

    def int_y(self):
        return np.sum(self.int(), axis=(0, 2))

    def int_x(self):
        return np.sum(self.int(), axis=(0, 1))

    def int_xy(self):
        return np.swapaxes(np.sum(self.int(), axis=0), 1, 0)

    def int_zx(self):
        return np.sum(self.int(), axis=1)

    def int_zy(self):
        return np.sum(self.int(), axis=2)

    def E(self):  # energy in the pulse [J]
        if self.Nz() > 1:
            return np.sum(self.int()) * self.Lz() / self.Nz() / speed_of_light
        else:
            return self.int()

    # propper scales in meters or 2 pi / meters
    def scale_kx(self):  # scale in meters or meters**-1
        if self.domain_xy == 's':  # space domain
            return np.linspace(-self.Lx() / 2, self.Lx() / 2, self.Nx())
        elif self.domain_xy == 'k':  # inverse space domain
            k = 2 * np.pi / self.dx
            return np.linspace(-k / 2, k / 2, self.Nx())
        else:
            raise AttributeError('Wrong domain_xy attribute')

    def scale_ky(self):  # scale in meters or meters**-1
        if self.domain_xy == 's':  # space domain
            return np.linspace(-self.Ly() / 2, self.Ly() / 2, self.Ny())
        elif self.domain_xy == 'k':  # inverse space domain
            k = 2 * np.pi / self.dy
            return np.linspace(-k / 2, k / 2, self.Ny())
        else:
            raise AttributeError('Wrong domain_xy attribute')

    def scale_kz(self):  # scale in meters or meters**-1
        if self.domain_z == 't':  # time domain
            return np.linspace(0, self.Lz(), self.Nz())
        elif self.domain_z == 'f':  # frequency domain
            dk = 2 * pi / self.Lz()
            k = 2 * pi / self.xlamds
            return np.linspace(k - dk / 2 * self.Nz(), k + dk / 2 * self.Nz(), self.Nz())
        else:
            raise AttributeError('Wrong domain_z attribute')

    def scale_x(self):  # scale in meters or radians
        if self.domain_xy == 's':  # space domain
            return self.scale_kx()
        elif self.domain_xy == 'k':  # inverse space domain
            return self.scale_kx() * self.xlamds / 2 / np.pi
        else:
            raise AttributeError('Wrong domain_xy attribute')

    def scale_y(self):  # scale in meters or radians
        if self.domain_xy == 's':  # space domain
            return self.scale_ky()
        elif self.domain_xy == 'k':  # inverse space domain
            return self.scale_ky() * self.xlamds / 2 / np.pi
        else:
            raise AttributeError('Wrong domain_xy attribute')

    def scale_z(self):  # scale in meters
        if self.domain_z == 't':  # time domain
            return self.scale_kz()
        elif self.domain_z == 'f':  # frequency domain
            return 2 * pi / self.scale_kz()
        else:
            raise AttributeError('Wrong domain_z attribute')


class WaistScanResults():

    def __init__(self):
        self.filePath = ''
        self.xlamds = None
        self.z_pos = np.array([])
        self.phdens_max = np.array([])
        self.phdens_onaxis = np.array([])
        self.fwhm_x = np.array([])
        self.fwhm_y = np.array([])
        self.std_x = np.array([])
        self.std_y = np.array([])
        self.z_maxint = None

    def fileName(self):
        return filename_from_path(self.filePath)


'''
    Genesis control
'''



def run_genesis(inp, launcher, read_level=2, assembly_ver='pyt', debug=1):
    '''
    Main function for executing Genesis code
    inp               - GenesisInput() object with genesis input parameters
    launcher          - MpiLauncher() object obtained via get_genesis_launcher() function
    read_level           - Parameter to read and calculate values from the output:
                     -1 - do not read
                      0 - read input only (header)
                      1 - read input and current profile
                      2 - read all values
    dfl_slipage_incl  - whether to dedicate time in order to keep the dfl slices, slipped out of the simulation window. if zero, reduces assembly time by ~30%
    assembly_ver      - version of the assembly script: 'sys' - system based, 'pyt' - python based
    '''

    # create experimental directory
    if inp.run_dir == None and inp.exp_dir == None:
        raise ValueError('run_dir and exp_dir are not specified!')

    if inp.run_dir == None:
        if inp.exp_dir[-1]!=os.path.sep:
            inp.exp_dir+=os.path.sep
        inp.run_dir = inp.exp_dir + 'run_' + str(inp.runid)

    try:
        os.makedirs(inp.run_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(inp.run_dir):
            pass
        else:
            raise

    if inp.stageid == None:
        inp_path = inp.run_dir + '/run.' + str(inp.runid) + '.inp'
        out_path = inp.run_dir + '/run.' + str(inp.runid) + '.gout'
    else:
        inp_path = inp.run_dir + '/run.' + str(inp.runid) + '.s' + str(inp.stageid) + '.inp'
        out_path = inp.run_dir + '/run.' + str(inp.runid) + '.s' + str(inp.stageid) + '.gout'

    inp_file = filename_from_path(inp_path)
    out_file = filename_from_path(out_path)

    # cleaning directory
    if debug > 0:
        print ('    removing old files')
    os.system('rm -rf ' + inp.run_dir + '/run.' + str(inp.runid) + '.s' + str(inp.stageid) + '*')  # to make sure all stage files are cleaned
    # os.system('rm -rf ' + out_path+'*') # to make sure out files are cleaned
    # os.system('rm -rf ' + inp_path+'*') # to make sure inp files are cleaned
    os.system('rm -rf ' + inp.run_dir + '/tmp.cmd')

    # create and fill necessary input files
    if inp.latticefile == None:
        if inp.lat != None:
            if debug > 1:
                print ('    writing ' + inp_file + '.lat')
            open(inp_path + '.lat', 'w').write(generate_lattice(inp.lat, unit=inp.xlamd*inp.delz, energy=inp.gamma0 * m_e_GeV))
            inp.latticefile = inp_file + '.lat'

    if inp.beamfile == None:
        if inp.beam != None:
            if debug > 1:
                print ('    writing ' + inp_file + '.beam')
            open(inp_path + '.beam', 'w').write(beam_file_str(inp.beam))
            inp.beamfile = inp_file + '.beam'

    if inp.edistfile == None:
        if inp.edist != None:
            if debug > 1:
                print ('    writing ' + inp_file + '.edist')
            write_edist_file(inp.edist, inp_path + '.edist', debug=1)
            inp.edistfile = inp_file + '.edist'

    if inp.partfile == None:
        if inp.dpa != None:
            if debug > 1:
                print ('    writing ' + inp_file + '.dpa')
            # print ('!!!!!!! no write_particle_file() function')
            write_dpa_file(inp.dpa, inp_path + '.dpa', debug=1)
            inp.partfile = inp_file + '.dpa'

    if inp.fieldfile == None:
        if inp.dfl != None:
            if debug > 1:
                print ('    writing ' + inp_file + '.dfl')
            write_dfl_file(inp.dfl, inp_path + '.dfl', debug=1)
            inp.fieldfile = inp_file + '.dfl'

    if inp.radfile == None:
        if inp.rad != None:
            if debug > 1:
                print ('    writing ' + inp_file + '.rad')
            open(inp_path + '.rad', 'w').write(rad_file_str(inp.rad))
            inp.radfile = inp_file + '.rad'

    if inp.outputfile == None:
        inp.outputfile = out_file
    open(inp_path, 'w').write(inp.input())
    open(inp.run_dir + '/tmp.cmd', 'w').write(inp_file + '\n')

    launcher.dir = inp.run_dir
    launcher.prepare()
    if debug > 1:
        print(inp.input())
    # RUNNING GENESIS ###
    launcher.launch()
    # RUNNING GENESIS ###

    # genesis output slices assembly
    if debug > 1:
        print (' ')
    if debug > 0:
        print ('    assembling slices')
    assembly_time = time.time()

    dfl_slipage_incl = True
    if assembly_ver == 'sys':

        if debug > 0:
            print ('      assembling *.out file')
        start_time = time.time()
        os.system('cat ' + out_path + '.slice* >> ' + out_path)
        os.system('rm ' + out_path + '.slice* 2>/dev/null')
        if debug > 1:
            print ('        done in %.2f seconds' % (time.time() - start_time))

        if debug > 0: # no dfln (n-harmonic) support yet
            print ('      assembling *.dfl file')
        start_time = time.time()
        if dfl_slipage_incl:
            os.system('cat ' + out_path + '.dfl.slice*  >> ' + out_path + '.dfl.tmp')
            #bytes=os.path.getsize(out_path +'.dfl.tmp')
            command = 'dd if=' + out_path + '.dfl.tmp of=' + out_path + '.dfl conv=notrunc conv=notrunc 2>/dev/null'  # obs='+str(bytes)+' skip=1
            os.system(command)
        else:
            os.system('cat ' + out_path + '.dfl.slice*  > ' + out_path + '.dfl')
        os.system('rm ' + out_path + '.dfl.slice* 2>/dev/null')
        os.system('rm ' + out_path + '.dfl.tmp 2>/dev/null')
        if debug > 1:
            print ('        done in %.2f seconds' % (time.time() - start_time))

        if debug > 0:
            print ('      assembling *.dpa file')
        start_time = time.time()
        os.system('cat ' + out_path + '.dpa.slice* >> ' + out_path + '.dpa')
        os.system('rm ' + out_path + '.dpa.slice* 2>/dev/null')
        if debug > 1:
            print ('        done in %.2f seconds' % (time.time() - start_time))
        
        if debug > 0:
            print ('      assembling *.fld file')
        start_time = time.time()
        if dfl_slipage_incl:
            os.system('cat ' + out_path + '.fld.slice*  >> ' + out_path + '.fld.tmp')
            #bytes=os.path.getsize(out_path +'.fld.tmp')
            command = 'dd if=' + out_path + '.fld.tmp of=' + out_path + '.fld conv=notrunc conv=notrunc 2>/dev/null'  # obs='+str(bytes)+' skip=1
            os.system(command)
        else:
            os.system('cat ' + out_path + '.fld.slice*  > ' + out_path + '.fld')
        os.system('rm ' + out_path + '.fld.slice* 2>/dev/null')
        os.system('rm ' + out_path + '.fld.tmp 2>/dev/null')
        if debug > 1:
            print ('        done in %.2f seconds' % (time.time() - start_time))

        if debug > 0:
            print ('      assembling *.par file')
        start_time = time.time()
        os.system('cat ' + out_path + '.par.slice* >> ' + out_path + '.par')
        os.system('rm ' + out_path + '.par.slice* 2>/dev/null')
        if debug > 1:
            print ('        done in %.2f seconds' % (time.time() - start_time))

    elif assembly_ver == 'pyt':
        # there is a bug with dfl assembly
        import glob
        ram = 1

        if debug > 0:
            print ('      assembling *.out file')
        start_time = time.time()
        assemble(out_path, ram=ram, debug=debug)
        os.system('rm ' + out_path + '.slice* 2>/dev/null')
        if debug > 1:
            print ('        done in %.2f seconds' % (time.time() - start_time))
        
        for i in range(10):        #read all possible harmonics (up to 10 now)
            ii=str(i)
            if ii=='0': ii=''
            if os.path.isfile(str(out_path + '.dfl' + ii)):
                if debug > 0:
                    print ('      assembling *.dfl'+ii+' file')
                start_time = time.time()
                assemble(out_path + '.dfl'+ii, overwrite=dfl_slipage_incl, ram=ram, debug=debug)
                os.system('rm ' + out_path + '.dfl'+ii+'.slice* 2>/dev/null')
                os.system('rm ' + out_path + '.dfl'+ii+'.tmp 2>/dev/null')
                if debug > 1:
                    print ('        done in %.2f seconds' % (time.time() - start_time))

        if os.path.isfile(str(out_path + '.dpa')):
            if debug > 0:
                print ('      assembling *.dpa file')
            start_time = time.time()
            assemble(out_path + '.dpa', ram=ram, debug=debug)
            os.system('rm ' + out_path + '.dpa.slice* 2>/dev/null')
            if debug > 1:
                print ('        done in %.2f seconds' % (time.time() - start_time))
                
        if os.path.isfile(str(out_path + '.fld')):
            if debug > 0:
                print ('      assembling *.fld file')
            start_time = time.time()
            assemble(out_path + '.fld', overwrite=dfl_slipage_incl, ram=ram, debug=debug)
            os.system('rm ' + out_path + '.fld.slice* 2>/dev/null')
            os.system('rm ' + out_path + '.fld.tmp 2>/dev/null')
            if debug > 1:
                print ('        done in %.2f seconds' % (time.time() - start_time))

        if os.path.isfile(str(out_path + '.par')):
            if debug > 0:
                print ('      assembling *.par file')
            start_time = time.time()
            assemble(out_path + '.par', ram=ram, debug=debug)
            os.system('rm ' + out_path + '.par.slice* 2>/dev/null')
            if debug > 1:
                print ('        done in %.2f seconds' % (time.time() - start_time))

    else:
        raise ValueError('assembly_ver should be either "sys" or "pyt"')
    # start_time = time.time()

    # print ('        done in %.2f seconds' % (time.time() - start_time))
    if debug > 0:
        print ('      total time %.2f seconds' % (time.time() - assembly_time))

    if read_level >= 0:
        out = read_out_file(out_path, read_level=read_level)
        return out
    else:
        return None


def assemble(fileName, remove=1, overwrite=0, ram=1, debug=1):
    '''
    assembles the fileName.slice* files into fileName
    remove - delete *.slice* files
    overwrite - writes *.slice* files on top of fileName.slice* starting from the beginning. Applicable for genesis dfl file assembly
    ram - store the *.slice* files in ram simultaneously
    '''
    import glob
    import sys
    # try:
        # if overwrite:
        # os.rename(fileName,fileName+'.slice999999')

    fins = glob.glob(fileName + '.slice*')
    fins.sort()

    if overwrite:
        fout = open(fileName, 'r+b')
    else:
        fout = open(fileName, 'ab')
    # else:
    #    fout = file(fileName,'a')
    N = len(fins)
    if ram == 1:
        idata = ''
        data = bytearray()
        if debug > 1:
            print('        reading ' + str(N) + ' slices to RAM...')
        index = 10
        for i, n in enumerate(fins):
            # if i/N>=index:
                # sys.stdout.write(str(index)+'%.')
                # index +=10
            fin = open(n, 'rb')
            while True:
                idata = fin.read(65536)
                if not idata:
                    break
                else:
                    data += idata
        if debug > 1:
            print('        writing...')
        fout.write(data)
        try:
            fin.close()
        except:
            pass

        if remove:
            os.system('rm ' + fileName + '.slice* 2>/dev/null')
            # os.remove(fins)
    else:
        for i, n in enumerate(fins):
            # if i/N>=index:
                # sys.stdout.write(str(index)+'%.')
                # index +=10
            fin = open(n, 'rb')
            while True:
                data = fin.read(65536)
                if not data:
                    break
                fout.write(data)
            fin.close()
            if remove:
                os.remove(fin.name)

    fout.close()
    # except:
    # print('        could not assemble '+fileName)


def create_exp_dir(exp_dir, run_ids):
    '''
    creates the folder structure nested in exp_dir folder.
    run_ids is a list of run numbers.
    resulting folder structure would be
    ..
    ---exp_dir\
    ------run_1\
    ------run_2\
    ------run_3\
    ------...
    '''
    if exp_dir[-1]!=os.path.sep:
        exp_dir+=os.path.sep
    for run_id in run_ids:

        try:
            run_dir = exp_dir + 'run_' + str(run_id)
            os.makedirs(run_dir)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(run_dir):
                pass
            else:
                raise

    try:
        res_dir = exp_dir + 'results'
        os.makedirs(res_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(res_dir):
            pass
        else:
            raise


def generate_input(up, beam, itdp=False):
    '''
    Create Genesis inp object with default input parameters
    '''
    inp = GenesisInput()

    # Next line added by GG 27.05.2016: it was in script


    #beam.emit_xn, beam.emit_yn = beam.emit[beam.C]
    #beam.gamma_rel = beam.E / (0.511e-3)
    #beam.emit_x = beam.emit_xn / beam.gamma_rel
    #beam.emit_y = beam.emit_yn / beam.gamma_rel

    inp.magin = 1
    #inp.zstop = 50
    #inp.nsec = 20

    inp.xlamd = up.lw
    inp.aw0 = up.K * np.sqrt(0.5)
    inp.awd = inp.aw0
    inp.delgam = beam.sigma_E / m_e_GeV
    inp.gamma0 = beam.E / m_e_GeV
    inp.rxbeam = np.sqrt(beam.emit_x * beam.beta_x)
    inp.rybeam = np.sqrt(beam.emit_y * beam.beta_y)

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
    inp.prad0 = 0
    inp.fbess0 = felParameters.fc
    inp.zrayl = felParameters.zr

    if itdp:
        inp.type = "tdp"
        # inp.DUMP_FIELDS = 1
        inp.ipseed = 132
        inp.ncar = 151
        #inp.nslice = 300
        inp.curlen = beam.tpulse * speed_of_light / 1e15
        inp.zsep = int(math.ceil(0.25 / (4 * pi * felParameters.rho)))  # 0.25 is the additional factor to be "on the safe side"

        print(inp.zsep)
        print(inp.xlamds)
        # inp.zsep = 8 * int(inp.curlen  / inp.nslice / inp.xlamds )
        inp.nslice = 8 * int(inp.curlen / inp.zsep / inp.xlamds)

    # inp.ntail = - int ( inp.nslice / 2 )
    inp.ntail = 0
    inp.npart = 2048
    inp.rmax0 = 9

    inp.delz = 1

    # print out FEL parameter estimates
    # printFelParameters(inp)
    return inp


def get_genesis_launcher(launcher_program=None):
    '''
    Returns MpiLauncher() object for given program
    '''
    host = socket.gethostname()

    launcher = MpiLauncher()
    if launcher_program != None:
        launcher.program = launcher_program
    else:

        if host.startswith('kolmogorov'):
            launcher.program = '/home/iagapov/workspace/xcode/codes/genesis/genesis < tmp.cmd | tee log'
        if host.startswith('max'):
            launcher.program = '/data/netapp/xfel/products/genesis/genesis < tmp.cmd | tee log'
        launcher.mpiParameters = '-x PATH -x MPI_PYTHON_SITEARCH -x PYTHONPATH'  # added -n
    #launcher.nproc = nproc
    return launcher

def get_genesis_new_launcher(launcher_program=None, mpi_mode=True):
    '''
    tmp for moga, to be solved in the future
    '''
    host = socket.gethostname()

    launcher = NewLauncher()
    
    if launcher_program != None:
        launcher.program = launcher_program
    else:
        if mpi_mode == True:
            launcher.program = "mpirun" + ' -x PATH -x MPI_PYTHON_SITEARCH -x PYTHONPATH ' + '/data/netapp/xfel/products/genesis/genesis < tmp.cmd | tee log'
        else:
            launcher.program = '/data/netapp/xfel/products/genesis_noparall/genesis_single < tmp.cmd | tee log'

    return launcher

''' 
   I/O functions
'''

''' 
   OUT
'''


def read_out_file(filePath, read_level=2, precision=float, debug=1):
    '''
    reads Genesis output from *.out file.
    returns GenesisOutput() object
    thanks gods Genesis3 out will be in hdf5!

    read_level -    0 = only header is processed. Very fast
                    1 = slice values are not processed. Current information is obtained, ~2x faster
                    2 = all contents are read
    precision - precision of stored values precision
    debug -     0 = no messages printed in console
                1 = basic info and execution time is printed
                2 = most detailed info is printed (real debug)
    '''
    import re
    out = GenesisOutput()
    out.filePath = filePath
    # out.fileName = filename_from_path(filePath)

    if debug > 0:
        print('    reading output file "' + out.fileName() + '"')
#    print '        - reading from ', fileName

    chunk = ''
    output_unsorted = []
    nSlice = 0

    wait_attempt = 6
    wait_time = 10
    while os.path.isfile(out.filePath) != True:
        print('!     waiting for "' + out.fileName() + '" ' + str(wait_time) + 's [' + str(wait_attempt) + ']')
        time.sleep(wait_time)  # wait for the .out file to be assembled
        wait_attempt -= 1
        if wait_attempt == 0:
            raise Exception('File ' + out.fileName() + ' not found')

    start_time = time.time()
    f = open(out.filePath, 'r')

    null = f.readline()
    for line in f:
        tokens = line.strip().split()

        if len(tokens) < 1:
            continue

        if tokens[0] == '**********':
            if read_level == 0:
                break
            chunk = 'slices'
            nSlice = int(tokens[3])
            if debug > 1:
                print ('      reading slice # ' + str(nSlice))

        if tokens[0] == 'power':
            chunk = 'slice'
            if len(out.sliceKeys) == 0:  # to record the first instance
                out.sliceKeys = list(copy(tokens))
                if debug > 0:
                    print ('      reading slice values ')
            continue

        if tokens[0] == '$newrun':
            chunk = 'input1'
            if debug > 0:
                print ('      reading input parameters')
            continue

        if tokens[0] == '$end':
            chunk = 'input2'
            continue

        if tokens == ['z[m]', 'aw', 'qfld']:
            chunk = 'magnetic optics'
            if debug > 0:
                print ('      reading magnetic optics ')
            continue

        if chunk == 'magnetic optics':
            z, aw, qfld = list(map(precision, tokens))
            out.z.append(z)
            out.aw.append(aw)
            out.qfld.append(qfld)

        if chunk == 'input1':
            tokens = line.replace('=', '').strip().split()
            out.parameters[tokens[0]] = tokens[1:]
            #out.parameters[tokens[0]] = tokens[0:]
            # print 'input:', tokens
        if chunk == 'input2':
            tokens = line.replace('=', '').strip().split()
            out.parameters['_'.join(tokens[1:])] = [tokens[0]]
            #out.parameters[tokens[0]] = tokens[0:]
            # print 'input:', tokens
#
        if chunk == 'slice' and read_level == 2:

            # tokens_fixed=re.sub(r'([0-9])\-([0-9])',r'\g<1>E-\g<2>',' '.join(tokens))
            # tokens_fixed=re.sub(r'([0-9])\+([0-9])',r'\g<1>E+\g<2>',tokens_fixed)
            # tokens=tokens_fixed.split()
            try:
                vals = list(map(precision, tokens))
            except ValueError:
                tokens_fixed = re.sub(r'([0-9])\-([0-9])', r'\g<1>E-\g<2>', ' '.join(tokens))
                tokens_fixed = re.sub(r'([0-9])\+([0-9])', r'\g<1>E+\g<2>', tokens_fixed)
                tokens_fixed = tokens_fixed.split()
                vals = list(map(precision, tokens_fixed))

            output_unsorted.append(vals)

        if chunk == 'slices':
            if len(tokens) == 2 and tokens[1] == 'current':
                # print tokens[1]
                out.I.append(float(tokens[0]))
                out.n.append(nSlice)
            elif len(tokens) == 3 and tokens[1] == 'scan':
                out.I.append(float(tokens[0]))
                out.n.append(nSlice)

    for parm in ['z', 'aw', 'qfld', 'I', 'n']:
        exec('out.' + parm + ' = np.array(out.' + parm + ')')

    if out('dgrid') == 0:
        rbeam = sqrt(out('rxbeam')**2 + out('rybeam')**2)
        ray = sqrt(out('zrayl') * out('xlamds') / np.pi * (1 + (out('zwaist') / out('zrayl'))**2))
        out.leng = out('rmax0') * (rbeam + ray)
    else:
        out.leng = 2 * out('dgrid')
    out.ncar = int(out('ncar'))  # number of mesh points
    
    #universal solution?
    out.leng=out('meshsize')*(out.ncar-1)


    if read_level == 0:
        print ('      returning *.out header')
        if debug > 0:
            print('      done in %.2f seconds' % (time.time() - start_time))
        return out

    out.nSlices = len(out.n)
    # int(out('history_records'))#number of slices in the output
    # print(nSlice)
    # print(out.nSlices)
    # if out.nSlices
    assert out('entries_per_record') != None, '.out header is missing!'
    out.nZ = int(out('entries_per_record'))  # number of records along the undulator

    if debug > 1:
        print ('        nSlices ' + str(out.nSlices))
    if debug > 1:
        print ('        nZ ' + str(out.nZ))

    assert nSlice != 0, '.out is empty!'

    assert(out.n[-1] - out.n[0]) == (len(out.n) - 1) * out('ishsty'), '.out is missing at least ' + str((out.n[-1] - out.n[0]) - (len(out.n) - 1) * out('ishsty')) + ' slices!'


    if read_level == 2:
        output_unsorted = np.array(output_unsorted)  # .astype(precision)
        # print out.sliceKeys
        for i in range(len(out.sliceKeys)):
            key = out.sliceKeys[int(i)]
            if key[0].isdigit():
                key='h'+key
            if debug > 1: 
                print ('      assembling',key.replace('-', '_').replace('<', '').replace('>', '')) 
            command = 'out.' + key.replace('-', '_').replace('<', '').replace('>', '') + ' = output_unsorted[:,' + str(i) + '].reshape((' + str(int(out.nSlices)) + ',' + str(int(out.nZ)) + '))'
            # print(command)
            exec(command)
        if hasattr(out, 'energy'):
            out.energy += out('gamma0')
        out.power_z = np.max(out.power, 0)
        out.sliceKeys_used = out.sliceKeys

    if out('itdp') == True:
        out.s = out('zsep') * out('xlamds') * (out.n - out.n[0])  # np.arange(0,out.nSlices)
        out.t = out.s / speed_of_light * 1.e+15
        out.dt = (out.t[1] - out.t[0]) * 1.e-15
        # out.dt=out('zsep') * out('xlamds') / speed_of_light
        out.beam_charge = np.sum(out.I * out.dt)
        out.sn_Imax = np.argmax(out.I)  # slice number with maximum current
        if read_level == 2:
            if debug > 0:
                print ('      calculating spectrum')
            out.spec = abs(np.fft.fft(np.sqrt(np.array(out.power)) * np.exp(1.j * np.array(out.phi_mid)), axis=0))**2 / sqrt(out.nSlices) / (2 * out.leng / out('ncar'))**2 / 1e10
            if debug > 1:
                print ('        done')
            e_0 = h_eV_s * speed_of_light / out('xlamds')
            out.freq_ev = h_eV_s * np.fft.fftfreq(len(out.spec), d=out('zsep') * out('xlamds') * out('ishsty') / speed_of_light) + e_0  # d=out.dt

            out.spec = np.fft.fftshift(out.spec, axes=0)
            out.freq_ev = np.fft.fftshift(out.freq_ev, axes=0)
            out.freq_lamd = h_eV_s * speed_of_light * 1e9 / out.freq_ev
            out.sliceKeys_used.append('spec')
            
            phase_fix = 1  # the way to display the phase, without constant slope caused by different radiation wavelength from xlamds. phase is set to 0 at maximum power slice.
            if phase_fix:
                if debug > 0:
                    print ('      fixing phase display')
                out.phi_mid_disp = deepcopy(out.phi_mid)
                for zi in range(shape(out.phi_mid_disp)[1]):
                    if debug > 1:
                        print ('      fixing phase display: ' + str(zi) + ' of ' + str(range(shape(out.phi_mid_disp)[1])))
                    maxspectrum_index = np.argmax(out.spec[:, zi])
                    maxspower_index = np.argmax(out.power[:, zi])
                    maxspectrum_wavelength = out.freq_lamd[maxspectrum_index] * 1e-9
                    phase = np.unwrap(out.phi_mid[:, zi])
                    phase_cor = np.arange(out.nSlices) * (maxspectrum_wavelength - out('xlamds')) / out('xlamds') * out('zsep') * 2 * pi
                    phase_fixed = phase + phase_cor
                    phase_fixed -= phase_fixed[maxspower_index]
                    n = 1
                    phase_fixed = (phase_fixed + n * pi) % (2 * n * pi) - n * pi
                    out.phi_mid_disp[:, zi] = phase_fixed
                out.sliceKeys_used.append('phi_mid_disp')

            rad_t_size_weighted = 1  # to average the radiation size over slices with radiation power as a weight
            if rad_t_size_weighted and out.nSlices != 1:
                if debug > 1:
                    print ('      averaging the radiation size properly')
                if np.amax(out.power) > 0:
                    weight = out.power + np.amin(out.power[out.power != 0]) / 1e6
                else:
                    weight = np.ones_like(out.power)
                out.rad_t_size_weighted = np.average(out.r_size * 1e6, weights=weight, axis=0)
                out.sliceKeys_used.append('rad_t_size_weighted')

    if out('iscan') != 0:
        out.scv = out.I  # scan value
        out.I = np.linspace(1, 1, len(out.scv))  # because used as a weight

    # tmp for back_compatibility
    if read_level == 2:
        out.power_int = out.power[:, -1]  # remove?
        out.max_power = np.amax(out.power_int)  # remove?
        for parm in [['power', 'p_int'],
                     ['energy', 'el_energy'],
                     ['e_spread', 'el_e_spread'],
                     ]:
            if hasattr(out, parm[0]):
                setattr(out, parm[1], getattr(out, parm[0]))
            for index, parm_key in enumerate(out.sliceKeys_used):
                if parm_key == parm[0]:
                    out.sliceKeys_used[index] = parm[1]
    #             delattr(out,parm[0])
        out.power = out.p_mid[:, -1]
        out.phi = out.phi_mid[:, -1]
        # out.energy=np.mean(out.p_int,axis=0)*out('xlamds')*out('zsep')*out.nSlices/speed_of_light
        if out('itdp'): 
            out.energy = np.sum(out.p_int * out.dt, axis=0)

    if debug > 0:
        print('      done in %.2f seconds' % (time.time() - start_time))
    return out


def read_out_file_stat(proj_dir, stage, run_inp=[], param_inp=[], debug=1):
    '''
    reads statistical info of Genesis simulations,
    returns GenStatOutput() object

    proj_dir - project directory of the following structure:
    proj_dir/run_<run_number>/run.<run_number>.s<stage_number>.gout*
    run_inp - list of genesis runs to be looked for [0:1000] by default
    param_inp - list of genesis output parameters to be processed
    debug - see read_out_file()
    '''
    if debug > 0:
        print ('    reading stat genesis output')
    start_time = time.time()

    if proj_dir[-1] != '/':
        proj_dir += '/'

    outlist = [GenesisOutput() for i in range(1000)]

    if run_inp == []:
        run_range = range(1000)
    else:
        run_range = run_inp

    run_range_good = []

    for irun in run_range:
        out_file = proj_dir + 'run_' + str(irun) + '/run.' + str(irun) + '.s' + str(stage) + '.gout'
        if os.path.isfile(out_file):
            if debug > 0:
                print ('      reading run', irun)
            outlist[irun] = read_out_file(out_file, read_level=2, debug=1)
            run_range_good.append(irun)
            # except:
    run_range = run_range_good
    
    # check if all gout have the same number of slices nSlice and history records nZ
    for irun in run_range[1:]:
        if outlist[irun].nSlices != outlist[run_range[0]].nSlices or outlist[irun].nZ != outlist[run_range[0]].nZ:
            raise ValueError('Non-uniform out objects (run %s)' %(irun))
    
    if debug: print(run_range)

    if param_inp == []:
        if debug > 1:
            print(outlist[run_range[0]].sliceKeys_used)
        param_range = outlist[run_range[0]].sliceKeys_used
    else:
        param_range = param_inp

    out_stat = GenStatOutput()
    for param in param_range:
        param_matrix = []
        for irun in run_range:
            if not hasattr(outlist[irun], param):
                continue
            else:
                param_matrix.append(deepcopy(getattr(outlist[irun], param)))

        param_matrix = np.array(param_matrix)
        if np.ndim(param_matrix) == 3:
            param_matrix = np.swapaxes(param_matrix, 0, 2)
        elif np.ndim(param_matrix) == 2 and shape(param_matrix)[1] == outlist[irun].nZ:
            param_matrix = np.swapaxes(param_matrix, 0, 1)[:, np.newaxis, :]
        else:
            pass
        setattr(out_stat, param, param_matrix)

    out_stat.stage = stage
    out_stat.dir = proj_dir
    out_stat.run = run_range
    out_stat.z = outlist[irun].z
    out_stat.s = outlist[irun].s
    out_stat.f = outlist[irun].freq_lamd
    out_stat.t = outlist[irun].t
    out_stat.dt = outlist[irun].dt
    
    out_stat.xlamds=outlist[irun]('xlamds')
    out_stat.filePath=proj_dir

    if debug > 0:
        print('      done in %.2f seconds' % (time.time() - start_time))
    return out_stat

def read_out_file_stat_u(file_tamplate, run_inp=[], param_inp=[], debug=1):
    '''
    reads statistical info of Genesis simulations,
    universal function for non-standard exp. folder structure
    returns GenStatOutput() object

    file_tamplate = template of the .out file path with # denoting run number
    run_inp - list of genesis runs to be looked for [0:1000] by default
    param_inp - list of genesis output parameters to be processed
    debug - see read_out_file()
    '''
    if debug > 0:
        print ('    reading stat genesis output')
    start_time = time.time()

    # if proj_dir[-1] != '/':
        # proj_dir += '/'

    outlist = [GenesisOutput() for i in range(1000)]

    if run_inp == []:
        run_range = range(1000)
    else:
        run_range = run_inp

    run_range_good = []

    for irun in run_range:
        out_file = file_tamplate.replace('#',str(irun))
        if os.path.isfile(out_file):
            if debug > 0:
                print ('      reading run', irun)
            outlist[irun] = read_out_file(out_file, read_level=2, debug=1)
            run_range_good.append(irun)
            # except:
    run_range = run_range_good
    
    # check if all gout have the same number of slices nSlice and history records nZ
    for irun in run_range[1:]:
        if outlist[irun].nSlices != outlist[run_range[0]].nSlices or outlist[irun].nZ != outlist[run_range[0]].nZ:
            raise ValueError('Non-uniform out objects (run %s)' %(irun))
    
    if debug>0: print(run_range)

    if param_inp == []:
        if debug > 1:
            print('      ',outlist[run_range[0]].sliceKeys_used)
        param_range = outlist[run_range[0]].sliceKeys_used
    else:
        param_range = param_inp

    out_stat = GenStatOutput()
    for param in param_range:
        param_matrix = []
        for irun in run_range:
            if not hasattr(outlist[irun], param):
                continue
            else:
                param_matrix.append(deepcopy(getattr(outlist[irun], param)))

        param_matrix = np.array(param_matrix)
        if np.ndim(param_matrix) == 3:
            param_matrix = np.swapaxes(param_matrix, 0, 2)
        elif np.ndim(param_matrix) == 2 and shape(param_matrix)[1] == outlist[irun].nZ:
            param_matrix = np.swapaxes(param_matrix, 0, 1)[:, np.newaxis, :]
        else:
            pass
        setattr(out_stat, param, param_matrix)

    out_stat.stage = stage
    out_stat.dir = proj_dir
    out_stat.run = run_range
    out_stat.z = outlist[irun].z
    out_stat.s = outlist[irun].s
    out_stat.f = outlist[irun].freq_lamd
    out_stat.t = outlist[irun].t
    out_stat.dt = outlist[irun].dt
    
    out_stat.xlamds=outlist[irun]('xlamds')
    out_stat.filePath=proj_dir

    if debug > 0:
        print('      done in %.2f seconds' % (time.time() - start_time))
    return out_stat
    
'''
    DFL
'''


def read_dfl_file_out(out, filePath=None, debug=1):
    '''
    More compact function than read_dfl_file() to read the file generated with known .out file
    Returns RadiationField object
    No need to pass all parameters (Nxy, Lxy, Lz, zsep, xlamds), they are read from GenesisOutput object
    out     - The relevant GenesisOutput object or path to it
    filePath- Path to *.dfl file.
        if = None, then it is assumed to be *.out.dfl
    '''
    if os.path.isfile(str(out)):
        out = read_out_file(out, read_level=0, debug=0)
    if not isinstance(out, GenesisOutput):
        raise ValueError('out is neither GenesisOutput() nor a valid path')

    if filePath == None:
        filePath = out.filePath + '.dfl'
    dfl = read_dfl_file(filePath, Nxy=out.ncar, Lxy=out.leng, zsep=out('zsep'), xlamds=out('xlamds'), debug=debug)
    return dfl


def read_dfl_file(filePath, Nxy, Lxy=None, zsep=None, xlamds=None, hist_rec=1, vartype=complex, debug=1):
    '''
    Function to read the Genesis output radiation file "dfl".
    Returns RadiationField() object
    filePath - path to *.dfl file
    Nxy - transverse mesh size, e.g. for 151*151 mesh, Nxy=151
    Lxy - transverse mesh size (2*dgrid)
    zsep - separation between slices in terms of wavelengths 
    xlamds  - wavelength of the radiation
    hist_rec - number of dfl records withinb in a single file (*.fld case), not finished!
    '''
    
    if debug > 0:
        print ('    reading radiation file')
    start_time = time.time()
    
    assert (Nxy % 1 == 0), 'Nxy nust be an integer'
    Nxy = int(Nxy)

    if not os.path.isfile(filePath):
        if debug:
            raise IOError('      ! dfl file ' + filePath + ' not found !')
        else:
            print ('      ! dfl file ' + filePath + ' not found !')
    else:
        if debug > 1:
            print ('        - reading from ' + filePath)

        b = np.fromfile(filePath, dtype=complex).astype(vartype)
        Nz = b.shape[0] / Nxy / Nxy / hist_rec
        assert (Nz % 1 == 0), 'Wrong Nxy or corrupted file'
        Nz = int(Nz)
        
        dfl = RadiationField()
        if hist_rec > 1:
            dfl.fld = b.reshape(hist_rec, Nz, Nxy, Nxy)
            dfl.fld = np.rollaxis(dfl.fld, 0, 4)
            print(dfl.shape())
        else:
            dfl.fld = b.reshape(Nz, Nxy, Nxy)
        dfl.dx = Lxy / dfl.Nx()
        dfl.dy = Lxy / dfl.Ny()
        dfl.dz = xlamds * zsep
        dfl.xlamds = xlamds
        dfl.filePath = filePath
        

        
        if debug > 0:
            print('      done in %.2f sec' % (time.time() - start_time))


        return dfl


def write_dfl_file(dfl, filePath=None, debug=1):
    '''
    Function to write the RadiationField object into filePath file
    dfl RadiationField object
    filePath - path top write the file
        if None then filePath = dfl.filePath
    '''
    if debug > 0:
        print ('    writing radiation file')
    start_time = time.time()

    if filePath == None:
        filePath = dfl.filePath

    d = dfl.fld.flatten().astype(complex128)
    d.tofile(filePath, format='complex')

    if debug > 0:
        print('      done in %.2f sec' % (time.time() - start_time))


'''
    DPA
'''


def read_dpa_file_out(out, filePath=None, debug=1):
    '''
    simplifies running the read_dpa_file() function
    reads GenesisOutput() object
    returns GenesisParticlesDump() object
    no need to set nbins and npart parameters as well as file_path (may be overrun)
    all automatically picked up from GenesisOutput() object
    '''

    if os.path.isfile(str(out)):
        out = read_out_file(out, read_level=0, debug=0)
    if not isinstance(out, GenesisOutput):
        raise ValueError('out is neither GenesisOutput() nor a valid path')

    if filePath == None:
        filePath = out.filePath + '.dpa'
    return read_dpa_file(filePath, nbins=out('nbins'), npart=out('npart'), debug=debug)


def read_dpa_file(filePath, nbins=4, npart=None, debug=1):
    '''
    reads genesis particle dump file *.dpa
    returns GenesisParticlesDump() object
    '''
    
    if not os.path.isfile(filePath):
        if debug:
            raise IOError('      ! dpa file ' + filePath + ' not found !')
        else:
            print ('      ! dpa file ' + filePath + ' not found !')
    else:
        if debug > 0:
            print ('    reading particle file')
            # print ('        - reading from ' + filePath)

    
    # if debug > 0:
        # print ('    reading particle file')
    dpa = GenesisParticlesDump()

    start_time = time.time()
    
    if os.path.isfile(filePath):
        b = np.fromfile(filePath, dtype=float)
    else:
        raise IOError('No such file: ' + filePath)
    # if debug: print("     read Particles in %s sec" % (time.time() - start_time))
    assert npart != None, 'number of particles per bin is not defined'
    npart = int(npart)
    nslice = int(len(b) / npart / 6)
    nbins = int(nbins)

    if debug > 1:
        print ('        nslice' + str(nslice))
        print ('        npart' + str(npart))
        print ('        nbins' + str(nbins))
    # print 'b=',nslice*npart*6
    b = b.reshape(nslice, 6, nbins, int(npart / nbins))
    dpa.e = b[:, 0, :, :]  # gamma
    dpa.ph = b[:, 1, :, :]
    dpa.x = b[:, 2, :, :]
    dpa.y = b[:, 3, :, :]
    dpa.px = b[:, 4, :, :]
    dpa.py = b[:, 5, :, :]
    dpa.filePath = filePath
    # dpa.fileName = filename_from_path(filePath)

    if debug > 0:
        print('      done in %.2f sec' % (time.time() - start_time))

    return dpa
    
def write_dpa_file(dpa, filePath=None, debug=1):
    
    if debug > 0:
        print ('    writing particle file')
    start_time = time.time()
    
    if dpa.__class__ != GenesisParticlesDump:
        raise ValueError('wrong particles object: should be GenesisParticlesDump')
    
    if filePath == None:
        filePath = dpa.filePath
        
    nslice,nbins,npart = dpa.e.shape
    b = np.zeros((nslice,6,nbins,npart))
    
    b[:, 0, :, :] = dpa.e
    b[:, 1, :, :] = dpa.ph
    b[:, 2, :, :] = dpa.x
    b[:, 3, :, :] = dpa.y
    b[:, 4, :, :] = dpa.px
    b[:, 5, :, :] = dpa.py
    b = b.flatten()
    b.tofile(filePath)
    
    if debug > 0:
        print('      done in %.2f sec' % (time.time() - start_time))



def max_dpa_dens(out, dpa, slice_pos=None, slice_num=None, repeat=1, bins=(50,50), debug=1):
    y_bins = bins[0]
    z_bins = bins[1]
    if slice_pos == slice_num == None:
        raise ValueError('specify either slice_pos or slice_num')

    if slice_num == None and out.nSlices > 1:
        if type(slice_pos) == str:
            if slice_pos == 'max_I':
                slice_num = np.argmax(out.I)
            elif slice_pos == 'max_P':
                slice_num = np.argmax(out.power)
            elif slice_pos == 'max_B':
                slice_num = np.argmax(out.bunching[:,-1])
            else:
                raise ValueError('slice_pos text should be "max_I" or "max_P"')
        else:
            if slice_pos < np.amin(out.s) or slice_pos > np.amax(out.s):
                raise ValueError('slice_pos outside out.s range')
            else:
                slice_num = np.where(out.s > slice_pos)[0][0]
    else:
        slice_num = 0
    nbins = shape(dpa.ph)[1]
    phase = deepcopy(dpa.ph[slice_num, :, :])
    gamma = deepcopy(dpa.e[slice_num, :, :])
    phase_flat=phase.flatten()
    gamma_flat=gamma.flatten()
    for irep in range(repeat-1):
        phase_flat=np.append(phase_flat,phase.flatten() + 2 * np.pi * (irep+1))
        gamma_flat=np.append(gamma_flat,gamma.flatten())
    # gamma_flat = energy_flat / m_e_eV
    
    hist, phase_scale, gamma_scale = np.histogram2d(phase_flat, gamma_flat, bins=bins)
    hist_max = np.amax(hist)
    gamma_idx = np.where(hist == hist_max)
    gamma_max = gamma_scale[gamma_idx[1]]
    return gamma_max[0], slice_num

def dpa2edist(out, dpa, num_part=1e5, smear=1, debug=1):
    '''
    reads GenesisParticlesDump() object
    returns GenesisElectronDist() object
    smear - whether to shuffle macroparticles smearing microbunching
    '''
    import random
    start_time = time.time()
    if debug > 0:
        print ('    transforming particle to distribution file')

    assert out('itdp') == True, '! steadystate Genesis simulation, dpa2dist() not implemented yet!'

    npart = int(out('npart'))
    # nslice=int(out('nslice'))
    nslice = int(out.nSlices * out('ishsty'))
    nbins = int(out('nbins'))
    xlamds = out('xlamds')
    zsep = int(out('zsep'))
    gen_I = out.I
    gen_t = out.t

    # if dpa==None:
    # dpa=out.filePath+'.dpa'
    # if dpa.__class__==str:
    # try:
    # dpa=read_dpa_file(dpa, nbins=nbins, npart=npart,debug=debug)
    # except IOError:
    # print ('      ERR: no such file "'+dpa+'"')
    # print ('      ERR: reading "'+out.filePath+'.dpa'+'"')
    # dpa=read_dpa_file(out.filePath+'.dpa', nbins=nbins, npart=npart,debug=debug)
    # if dpa.__class__!=GenesisParticles:
    # print('   could not read particle file')

    m = np.arange(nslice)
    m = np.tile(m, (nbins, npart / nbins, 1))
    m = np.rollaxis(m, 2, 0)
    # print('shape_m='+str(shape(m)))
    
    
    if smear:
        z = dpa.ph * xlamds / 2 / pi + m * xlamds * zsep + xlamds * zsep * (1 - np.random.random((nslice, nbins, int(npart / nbins))))
    else:
        z = dpa.ph * xlamds / 2 / pi + m * xlamds * zsep

    t = np.array(z / speed_of_light)

    t_scale = np.linspace(0, nslice * zsep * xlamds / speed_of_light * 1e15, nslice)

    pick_n = np.interp(t_scale, gen_t, gen_I)
    if debug > 1:
        print('sum pick_n=' + str(sum(pick_n)))
    if debug > 1:
        print('npart=' + str(npart))
    if debug > 1:
        print('num_part=' + str(num_part))
    pick_n = pick_n / sum(pick_n) * num_part
    if max(pick_n) > npart:
        pick_n = pick_n / max(pick_n) * npart
    pick_n = pick_n.astype(int)
    # ratio=ceil(num_part/np.sum(pick_n))
    ratio = 1
    
    t = np.reshape(t, (nslice, npart))
    e = np.reshape(dpa.e, (nslice, npart))
    x = np.reshape(dpa.x, (nslice, npart))
    y = np.reshape(dpa.y, (nslice, npart))
    px = np.reshape(dpa.px, (nslice, npart))
    py = np.reshape(dpa.py, (nslice, npart))

    edist = GenesisElectronDist()
    for i in np.arange(nslice):
        for ii in np.arange(int(ratio)):
            pick_i = random.sample(range(npart), pick_n[i])
            edist.t = append(edist.t, t[i, pick_i])
            edist.g = append(edist.g, e[i, pick_i])
            edist.x = append(edist.x, x[i, pick_i])
            edist.y = append(edist.y, y[i, pick_i])
            edist.xp = append(edist.xp, px[i, pick_i])
            edist.yp = append(edist.yp, py[i, pick_i])

    edist.t = edist.t * (-1) + max(edist.t)
    edist.xp /= edist.g
    edist.yp /= edist.g

    edist.x = np.flipud(edist.x)
    edist.y = np.flipud(edist.y)
    edist.xp = np.flipud(edist.xp)
    edist.yp = np.flipud(edist.yp)
    edist.t = np.flipud(edist.t)
    edist.g = np.flipud(edist.g)

    edist.part_charge = out.beam_charge / edist.len()
    # edist.charge=out.beam_charge
    edist.filePath = dpa.filePath + '.edist'
    # print 'max_y_out', np.amax(t_out)
    # print 'e_out', np.amax(e_out),np.amin(e_out)

    if debug > 0:
        print('      done in %.2f sec' % (time.time() - start_time))

    return edist


'''
    EDIST
'''


def read_edist_file_out(out, debug=1):
    return read_dist_file(out.filePath + '.edist', debug=debug)


def read_edist_file(filePath, debug=1):
    '''
    reads particle distribution file (distfile in genesis input)
    returns GenesisElectronDist() 
    '''
    edist = GenesisElectronDist()
    edist.filePath = filePath

    if debug > 0:
        print ('    reading particle distribution file')
    start_time = time.time()

    dist_column_values = {}
    f = open(filePath, 'r')
    null = f.readline()
    for line in f:
        tokens = line.strip().split()

        if len(tokens) < 2:
            continue

        if tokens[0] == "?" and tokens[1] == "CHARGE":
            charge = float(tokens[3])

        if tokens[0] == "?" and tokens[1] == "COLUMNS":
            dist_columns = tokens[2:]
            for col in dist_columns:
                dist_column_values[col] = []
            if debug > 1:
                print(''.join(str(i) + ' ' for i in dist_columns))

        if tokens[0] != "?":
            for i in range(0, len(tokens)):
                dist_column_values[dist_columns[i]].append(float(tokens[i]))

    edist.x = np.array(dist_column_values['X'])
    edist.y = np.array(dist_column_values['Y'])
    edist.xp = np.array(dist_column_values['XPRIME'])
    edist.yp = np.array(dist_column_values['YPRIME'])
    edist.t = np.array(dist_column_values['T'])
    edist.g = np.array(dist_column_values['P'])

    edist.x = np.flipud(edist.x)
    edist.y = np.flipud(edist.y)
    edist.xp = np.flipud(edist.xp)
    edist.yp = np.flipud(edist.yp)
    edist.t = np.flipud(edist.t)
    edist.g = np.flipud(edist.g)

    edist.part_charge = charge / edist.len()

    if debug > 0:
        print('      done in %.2f sec' % (time.time() - start_time))

    return edist


def cut_edist(edist,
              t_lim=(-inf, inf),
              g_lim=(-inf, inf),
              x_lim=(-inf, inf),
              xp_lim=(-inf, inf),
              y_lim=(-inf, inf),
              yp_lim=(-inf, inf), debug=1):
    '''
    cuts GenesisElectronDist() in phase space
    '''

    from numpy import logical_or

    if debug > 0:
        print ('    cutting particle distribution file')
    start_time = time.time()

    index_t = logical_or(edist.t < t_lim[0], edist.t > t_lim[1])
    index_e = logical_or(edist.g < g_lim[0], edist.g > g_lim[1])
    index_x = logical_or(edist.x < x_lim[0], edist.x > x_lim[1])
    index_y = logical_or(edist.y < y_lim[0], edist.y > y_lim[1])
    index_px = logical_or(edist.xp < xp_lim[0], edist.xp > xp_lim[1])
    index_py = logical_or(edist.yp < yp_lim[0], edist.yp > yp_lim[1])

    index = np.logical_or.reduce((index_t, index_e, index_x, index_y, index_px, index_py))
    index = np.where(index)

    edist_f = deepcopy(edist)

    for parm in ['t', 'g', 'xp', 'yp', 'x', 'y']:
        if hasattr(edist_f, parm):
            setattr(edist_f, parm, np.delete(getattr(edist_f, parm), index))

    if debug > 0:
        print('      %.2f percent cut' % ((edist.charge() - edist_f.charge()) / edist.charge() * 100))
    if debug > 0:
        print('      done in %s sec' % (time.time() - start_time))

    return edist_f


def repeat_edist(edist, factor, smear=1):
    '''
    dublicates the GenesisElectronDist() by given factor
    smear - smear new particles by 1e-3 of standard deviation of parameter
    '''

    edist_out = GenesisElectronDist()

    edist_out.x = np.repeat(edist.x, factor)
    edist_out.y = np.repeat(edist.y, factor)
    edist_out.xp = np.repeat(edist.xp, factor)
    edist_out.yp = np.repeat(edist.yp, factor)
    edist_out.t = np.repeat(edist.t, factor)
    edist_out.g = np.repeat(edist.g, factor)
    edist_out.part_charge = edist.part_charge / factor

    if smear:
        n_par = edist_out.len()
        smear_factor = 1e-3  # smear new particles by smear_factor of standard deviation of parameter

        edist_out.x += np.random.normal(scale=np.std(edist_out.x) * smear_factor, size=n_par)
        edist_out.y += np.random.normal(scale=np.std(edist_out.y) * smear_factor, size=n_par)
        edist_out.xp += np.random.normal(scale=np.std(edist_out.xp) * smear_factor, size=n_par)
        edist_out.yp += np.random.normal(scale=np.std(edist_out.yp) * smear_factor, size=n_par)
        edist_out.t += np.random.normal(scale=np.std(edist_out.t) * smear_factor, size=n_par)
        edist_out.g += np.random.normal(scale=np.std(edist_out.g) * smear_factor, size=n_par)

    return edist_out


def write_edist_file(edist, filePath=None, debug=1):
    '''
    writes GenesisElectronDist() into filePath folder
    '''
    # REQUIRES NUMPY 1.7
    # header='? VERSION = 1.0 \n? SIZE = %s \n? CHARGE = %E \n? COLUMNS X XPRIME Y YPRIME T P'%(len(edist.x),charge)
    # np.savetxt(filePath_write, np.c_[edist.x,edist.xp,edist.y,edist.yp,edist.t,edist.g],header=header,fmt="%E", newline='\n',comments='')

    if debug > 0:
        print ('    writing particle distribution file')
    start_time = time.time()

    if filePath == None:
        filePath = edist.filePath

    header = '? VERSION = 1.0 \n? SIZE = %s \n? CHARGE = %E \n? COLUMNS X XPRIME Y YPRIME T P\n' % (edist.len(), edist.charge())
    f = open(filePath, 'w')
    f.write(header)
    f.close()
    f = open(filePath, 'ab')
    np.savetxt(f, np.c_[edist.x, edist.xp, edist.y, edist.yp, edist.t, edist.g], fmt="%e", newline='\n')
    f.close()

    if debug > 0:
        print('      done in %.2f sec' % (time.time() - start_time))


def edist2beam(edist, step=1e-7):
    '''
    reads GenesisElectronDist()
    returns GenesisBeam()
    step [m] - long. size ob bin to calculate distribution parameters
    '''

    from numpy import mean, std

    part_c = edist.part_charge
    t_step = step / speed_of_light
    t_min = min(edist.t)
    t_max = max(edist.t)
    dist_t_window = t_max - t_min
    npoints = int(dist_t_window / t_step)
    t_step = dist_t_window / npoints
    beam = GenesisBeam()
    for parm in ['I',
                 'z',
                 'ex',
                 'ey',
                 'betax',
                 'betay',
                 'alphax',
                 'alphay',
                 'x',
                 'y',
                 'px',
                 'py',
                 'g0',
                 'dg',
                 ]:
        setattr(beam, parm, np.zeros((npoints - 1)))

    for i in range(npoints - 1):
        indices = (edist.t > t_min + t_step * i) * (edist.t < t_min + t_step * (i + 1))
        beam.z[i] = (t_min + t_step * (i + 0.5)) * speed_of_light
        dist_mean_g = beam.g0[i]

        if sum(indices) > 2:
            dist_g = edist.g[indices]
            dist_x = edist.x[indices]
            dist_y = edist.y[indices]
            dist_px = edist.xp[indices]
            dist_py = edist.yp[indices]
            dist_mean_g = mean(dist_g)

            beam.I[i] = sum(indices) * part_c / t_step
            beam.g0[i] = mean(dist_g)
            beam.dg[i] = np.std(dist_g)
            beam.x[i] = mean(dist_x)
            beam.y[i] = mean(dist_y)
            beam.px[i] = mean(dist_px)
            beam.py[i] = mean(dist_py)
            beam.ex[i] = dist_mean_g * (mean(dist_x**2) * mean(dist_px**2) - mean(dist_x * dist_px)**2)**0.5
            # if beam.ex[i]==0: beam.ey[i]=1e-10
            beam.ey[i] = dist_mean_g * (mean(dist_y**2) * mean(dist_py**2) - mean(dist_y * dist_py)**2)**0.5
            # if beam.ey[i]==0: beam.ey[i]=1e-10
            beam.betax[i] = dist_mean_g * mean(dist_x**2) / beam.ex[i]
            beam.betay[i] = dist_mean_g * mean(dist_y**2) / beam.ey[i]
            beam.alphax[i] = -dist_mean_g * mean(dist_x * dist_px) / beam.ex[i]
            beam.alphay[i] = -dist_mean_g * mean(dist_y * dist_py) / beam.ey[i]

    idx = np.where(np.logical_or.reduce((beam.I == 0, beam.g0 == 0, beam.betax > mean(beam.betax) * 10, beam.betay > mean(beam.betay) * 10)))
    del beam[idx]
    # for i in reversed(range(npoints-1)):
    # if beam.I[i]==0:
    # np.delete(beam.I,i)
    # np.delete(beam.g0,i)
    # np.delete(beam.dg,i)
    # np.delete(beam.x,i)
    # np.delete(beam.y,i)
    # np.delete(beam.px,i)
    # np.delete(beam.py,i)
    # np.delete(beam.ex,i)
    # np.delete(beam.ey,i)
    # np.delete(beam.betax,i)
    # np.delete(beam.betay,i)
    # np.delete(beam.alphax,i)
    # np.delete(beam.alphay,i)

    beam.columns = ['ZPOS', 'GAMMA0', 'DELGAM', 'EMITX', 'EMITY', 'BETAX', 'BETAY', 'XBEAM', 'YBEAM', 'PXBEAM', 'PYBEAM', 'ALPHAX', 'ALPHAY', 'CURPEAK', 'ELOSS']

    beam.idx_max = np.argmax(beam.I)
    beam.eloss = np.zeros_like(beam.z)
    # beam.fileName=edist.fileName+'.beam'
    beam.filePath = edist.filePath + '.beam'

    beam.zsep = beam.z[1] - beam.z[0]  # get rid of

    return(beam)


'''
    BEAM
'''


# def read_beam_file_out(out, debug=1):
    # return read_beam_file(out.filePath, debug=debug)


def read_beam_file(filePath, debug=1):
    '''
    reads beam file from filePath folder
    returns GenesisBeam()
    '''
    if debug > 0:
        print ('    reading beam file')
    start_time = time.time()

    beam = GenesisBeam()

    f = open(filePath, 'r')
    null = f.readline()
    for line in f:
        tokens = line.strip().split()

        if len(tokens) < 2:
            continue

        # print tokens[0:2]

        if tokens[0] == "?" and tokens[1] == "COLUMNS":
            beam.columns = tokens[2:]
            for col in beam.columns:
                beam.column_values[col] = []

            print (beam.columns)

        if tokens[0] != "?":
            # print tokens
            for i in range(0, len(tokens)):
                beam.column_values[beam.columns[i]].append(float(tokens[i]))

    # print beam.columns

    beam.z = np.array(beam.column_values['ZPOS'])
    beam.zsep = beam.z[1] - beam.z[0]
    beam.I = np.array(beam.column_values['CURPEAK'])
    beam.idx_max = np.argmax(beam.I)

    dict = {'ex'}

    for parm in [['ex', 'EMITX'],
                 ['ey', 'EMITY'],
                 ['betax', 'BETAX'],
                 ['betay', 'BETAY'],
                 ['alphax', 'ALPHAX'],
                 ['alphay', 'ALPHAY'],
                 ['x', 'XBEAM'],
                 ['y', 'YBEAM'],
                 ['px', 'PXBEAM'],
                 ['py', 'PYBEAM'],
                 ['g0', 'GAMMA0'],
                 ['dg', 'DELGAM'],
                 ]:
        if parm[1] in beam.column_values.keys():
            setattr(beam, parm[0], np.array(beam.column_values[parm[1]]))
        else:
            setattr(beam, parm[0], np.zeros_like(beam.z))

    try:
        beam.eloss = np.array(beam.column_values['ELOSS'])
    except:
        beam.eloss = np.zeros_like(beam.I)

    beam.filePath = filePath

    if debug > 0:
        print('      done in %.2f sec' % (time.time() - start_time))

    return beam


def beam_file_str(beam):
    '''
    reads GenesisBeam()
    returns string of electron beam file, suitable for Genesis
    
    '''
    # header = "# \n? VERSION = 1.0\n? SIZE ="+str(len(beam.column_values['ZPOS']))+"\n? COLUMNS ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY XBEAM YBEAM PXBEAM PYBEAM ALPHAX ALPHAY CURPEAK ELOSS\n"
    header = "# \n? VERSION = 1.0\n? SIZE =" + str(len(beam.z)) + "\n? COLUMNS"
    # ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY XBEAM YBEAM PXBEAM PYBEAM ALPHAX ALPHAY CURPEAK ELOSS\n"
    for col in beam.columns:
        header = header + " " + col
    header += "\n"

    f_str = header

    for parm in [['z', 'ZPOS'],
                 ['I', 'CURPEAK'],
                 ['ex', 'EMITX'],
                 ['ey', 'EMITY'],
                 ['betax', 'BETAX'],
                 ['betay', 'BETAY'],
                 ['alphax', 'ALPHAX'],
                 ['alphay', 'ALPHAY'],
                 ['x', 'XBEAM'],
                 ['y', 'YBEAM'],
                 ['px', 'PXBEAM'],
                 ['py', 'PYBEAM'],
                 ['g0', 'GAMMA0'],
                 ['dg', 'DELGAM'],
                 ['eloss', 'ELOSS'],
                 ]:
        try:
            beam.column_values[parm[1]] = getattr(beam, parm[0])
        except:
            pass
    # beam.column_values['ZPOS'] = beam.z
    # beam.column_values['CURPEAK'] = beam.I

    # try:
        # beam.column_values['EMITX'] = beam.ex
        # beam.column_values['EMITY'] = beam.ey

        # beam.column_values['BETAX'] = beam.betax
        # beam.column_values['BETAY'] = beam.betay

        # beam.column_values['ALPHAX'] = beam.alphax
        # beam.column_values['ALPHAY'] = beam.alphay

        # beam.column_values['XBEAM'] = beam.x
        # beam.column_values['YBEAM'] = beam.y

        # beam.column_values['PXBEAM'] = beam.px
        # beam.column_values['PYBEAM'] = beam.py

        # beam.column_values['GAMMA0'] = beam.g0
        # beam.column_values['DELGAM'] = beam.dg

        # beam.column_values['ELOSS'] = beam.eloss
    # except:
        # pass

    for i in range(len(beam.z)):
        for col in beam.columns:
            buf = str(beam.column_values[col][i])
            f_str = f_str + buf + ' '
        f_str = f_str.rstrip() + '\n'

    return f_str


# def add_wake_to_beamf(beamf, new_beamf):
    # beam = read_beam_file(beamf)
    # s, bunch, wake = w.xfel_pipe_wake(s=array(beam.z), current=array(beam.I))
    # print ('read ' + str(len(wake)) + ' slice values')
    # beam.eloss = wake[::-1]

    # f = open(new_beamf, 'w')
    # f.write(beam_file_str(beam))
    # f.close()


def zero_wake_at_ipk(beamf):
    '''
    reads GenesisBeam()
    shifts the wake pforile so that 
    at maximum current slice wake is zero
    returns GenesisBeam()
    
    allows to account for wake losses without 
    additional linear undulator tapering
    '''
    beamf_new = deepcopy(beamf)
    beamf_new.idx_max_refresh()
    beamf_new.eloss -= beamf_new.eloss[beamf_new.idx_max]
    return beamf_new


def set_beam_energy(beam, E_GeV_new):
    '''
    reads GenesisBeam()
    returns GenesisBeam()
    sets the beam energy with peak current to E_GeV_new
    '''
    beam_peak = get_beam_peak(beam)
    E_GeV_old = beam_peak.E
    g_e_old = E_GeV_old / m_e_GeV
    g_e_new = E_GeV_new / m_e_GeV
    #beam.g0 = beam.g0 / E_GeV_old * E_GeV_new
    beam.g0 = beam.g0 - g_e_old + g_e_new
    return beam


def transform_beam_twiss(beam, s=None, transform=None):
    # transform = [[beta_x,alpha_x],[beta_y, alpha_y]]
    if transform == None:
        return beam
    else:
        beam_peak = get_beam_peak(beam)
        if s == None:
            idx = beam_peak.idx_max
        else:
            idx = np.where(beam.z > s)[0][0]

        g1x = np.matrix([[beam.betax[idx], beam.alphax[idx]],
                         [beam.alphax[idx], (1 + beam.alphax[idx]**2) / beam.betax[idx]]])

        g1y = np.matrix([[beam.betay[idx], beam.alphay[idx]],
                         [beam.alphay[idx], (1 + beam.alphay[idx]**2) / beam.betay[idx]]])

        b2x = transform[0][0]
        a2x = transform[0][1]

        b2y = transform[1][0]
        a2y = transform[1][1]

        g2x = np.matrix([[b2x, a2x],
                         [a2x, (1 + a2x**2) / b2x]])

        g2y = np.matrix([[b2y, a2y],
                         [a2y, (1 + a2y**2) / b2y]])

        Mix, Mx = find_transform(g1x, g2x)
        Miy, My = find_transform(g1y, g2y)

        # print Mi

        betax_new = []
        alphax_new = []
        betay_new = []
        alphay_new = []

        for i in range(len(beam.z)):
            g1x = np.matrix([[beam.betax[i], beam.alphax[i]],
                             [beam.alphax[i], (1 + beam.alphax[i]**2) / beam.betax[i]]])

            gx = Mix.T * g1x * Mix

            g1y = np.matrix([[beam.betay[i], beam.alphay[i]],
                             [beam.alphay[i], (1 + beam.alphay[i]**2) / beam.betay[i]]])

            gy = Miy.T * g1y * Miy

            # print i, gx[0,1], g1x[0,1]

            betax_new.append(gx[0, 0])
            alphax_new.append(gx[0, 1])
            betay_new.append(gy[0, 0])
            alphay_new.append(gy[0, 1])

        beam.betax = betax_new
        beam.betay = betay_new
        beam.alphax = alphax_new
        beam.alphay = alphay_new

        return beam


def add_alpha_beam(beam):

    beam.alphax = beam.column_values['ALPHAX'] = np.zeros_like(beam.g0)
    beam.alphay = beam.column_values['ALPHAY'] = np.zeros_like(beam.g0)

    beam.columns = list(beam.column_values.keys())


def cut_beam(beam=None, cut_z=[-inf, inf]):
    '''
    cuts GenesisBeam() object longitudinally
    cut_z [m] - limits of the cut
    '''
    if np.amin(beam.z) < cut_z[0] or np.amax(beam.z) > cut_z[1]:

        condition = (beam.z > cut_z[0]) * (beam.z < cut_z[1])
        print(sum(condition))
        beam_new = GenesisBeam()
        beam_new.column_values = beam.column_values
        beam_new.columns = beam.columns

        for parm in ['x', 'px', 'y', 'py', 'z', 'I', 'ex', 'ey', 'g0', 'dg', 'eloss', 'betax', 'betay', 'alphax', 'alphay']:
            if hasattr(beam, parm):
                setattr(beam_new, parm, np.extract(condition, getattr(beam, parm)))

        for parm in ['filePath', ]:
            if hasattr(beam, parm):
                setattr(beam_new, parm, getattr(beam, parm))

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

        # beam_new.zsep = beam.zsep
        zmax, Imax = peaks(beam_new.z, beam_new.I, n=1)
        beam_new.idx_max = np.where(beam_new.z == zmax)[0][0]
    else:
        beam_new = beam
    return beam_new


def get_beam_s(beam=None, s=0):
    '''
    obtains values of the beam at s position
    '''
    if len(beam.I) > 1:  # and np.amax(beam.I)!=np.amin(beam.I):
        slice = np.where(beam.z >= s)[0][0]

        # beam_new=deepcopy(beam)
        beam_new = Beam()

        beam_new.idx_max = slice
        beam_new.I = beam.I[slice]
        beam_new.alpha_x = beam.alphax[slice]
        beam_new.alpha_y = beam.alphay[slice]
        beam_new.beta_x = beam.betax[slice]
        beam_new.beta_y = beam.betay[slice]
        beam_new.emit_xn = beam.ex[slice]
        beam_new.emit_yn = beam.ey[slice]
        beam_new.gamma_rel = beam.g0[slice]
        beam_new.sigma_E = beam.dg[slice] * (0.510998e-3)
        beam_new.xp = beam.px[slice]
        beam_new.yp = beam.py[slice]
        beam_new.x = beam.x[slice]
        beam_new.y = beam.y[slice]

        beam_new.E = beam_new.gamma_rel * (0.510998e-3)
        beam_new.emit_x = beam_new.emit_xn / beam_new.gamma_rel
        beam_new.emit_y = beam_new.emit_yn / beam_new.gamma_rel

        beam_new.tpulse = (beam.z[-1] - beam.z[0]) / speed_of_light * 1e15 / 6  # [fs]
        beam_new.C = np.trapz(beam.I, x=np.array(beam.z) / speed_of_light) * 1e9  # bunch charge[nC]

        for parm in ['alphax', 'alphay', 'betax', 'betay', 'z', 'ex', 'ey', 'g0', 'dg', 'px', 'py']:
            if hasattr(beam_new, parm):
                delattr(beam_new, parm)

        # beam_new.x = np.extract(condition,beam.x)

    else:
        beam_new = beam
    return beam_new


def get_beam_peak(beam=None):
    '''
    obtains the peak current values of the beam
    '''
    if len(beam.I) > 1:  # and np.amax(beam.I)!=np.amin(beam.I):
        pkslice = np.argmax(beam.I)

        # beam_new=deepcopy(beam)
        beam_new = Beam()

        beam_new.idx_max = pkslice
        beam_new.I = beam.I[pkslice]
        beam_new.alpha_x = beam.alphax[pkslice]
        beam_new.alpha_y = beam.alphay[pkslice]
        beam_new.beta_x = beam.betax[pkslice]
        beam_new.beta_y = beam.betay[pkslice]
        beam_new.emit_xn = beam.ex[pkslice]
        beam_new.emit_yn = beam.ey[pkslice]
        beam_new.gamma_rel = beam.g0[pkslice]
        beam_new.sigma_E = beam.dg[pkslice] * (0.510998e-3)
        beam_new.xp = beam.px[pkslice]
        beam_new.yp = beam.py[pkslice]
        beam_new.x = beam.x[pkslice]
        beam_new.y = beam.y[pkslice]

        beam_new.E = beam_new.gamma_rel * (0.510998e-3)
        beam_new.emit_x = beam_new.emit_xn / beam_new.gamma_rel
        beam_new.emit_y = beam_new.emit_yn / beam_new.gamma_rel

        beam_new.tpulse = (beam.z[-1] - beam.z[0]) / speed_of_light * 1e15 / 8  # sigma [fs] (8 sigmas=full window)
        beam_new.C = np.trapz(beam.I, x=np.array(beam.z) / speed_of_light) * 1e9  # bunch charge[nC]

        for parm in ['alphax', 'alphay', 'betax', 'betay', 'z', 'ex', 'ey', 'g0', 'dg', 'px', 'py']:
            if hasattr(beam_new, parm):
                delattr(beam_new, parm)

        # beam_new.x = np.extract(condition,beam.x)

    else:
        beam_new = beam
    return beam_new


def find_transform(g1, g2):
    '''
    find transform from twiss matrix g1 to g2: x -> M x, g -> Mi.T g M
    '''
    l1, u1 = np.linalg.eig(g1)
    l2, u2 = np.linalg.eig(g2)

    M1 = np.matrix([[u1[0, 0], u1[1, 0]],
                    [-u1[1, 0], u1[0, 0]]])

    d = sqrt(l1[0] / l2[0])

    M2 = np.matrix([[d, 0],
                    [0, 1.0 / d]])

    M2d = np.matrix([[1. / d, 0],
                     [0, d]])

    M3 = np.matrix([[u2[0, 0], -u2[1, 0]],
                    [u2[1, 0], u2[0, 0]]])

    return np.linalg.inv(M3 * M2 * M1), M3 * M2d * M1


def write_beam_file(filePath, beam, debug=0):
    if debug > 0:
        print ('    writing beam file')
    start_time = time.time()

    fd = open(filePath, 'w')
    fd.write(beam_file_str(beam))
    fd.close()

    if debug > 0:
        print('      done in %.2f sec' % (time.time() - start_time))


'''
    RAD
'''


def create_rad_file(p_duration_s=None, p_intensity=None, beam=None, offset=None, out_file='tmp.rad'):

    def rad_file_str2(beam):
        # header = "# \n? VERSION = 1.0\n? SIZE = "+str(len(beam.z))+"\n? OFFSET = "+str(beam.offset)+"\n? COLUMNS ZPOS PRAD0 \n"
        # header = "# \n? VERSION = 1.0\n? SIZE = "+str(len(beam.z))+"\n? OFFSET = "+str(beam.offset2)+"\n? COLUMNS ZPOS PRAD0 \n"
        header = "? VERSION = 1.0\n? SIZE = " + str(len(beam.z)) + "\n? OFFSET = " + str(beam.offset2) + "\n? COLUMNS ZPOS PRAD0 \n"
        f_str = header

        for i in range(len(beam.z)):
            f_str_tmp = str(beam.z[i]) + ' ' + str(beam.prad0[i]) + '\n'
            f_str += f_str_tmp

        return f_str

    temporal_delay = offset
    extention_twindow = 2
    start = beam.z[0] - extention_twindow * (beam.z[math.floor(len(beam.z) / 2)] - beam.z[0])
    stop = beam.z[len(beam.z) - 1] + extention_twindow * (beam.z[len(beam.z) - 1] - beam.z[math.ceil(len(beam.z) / 2)])
    num = len(beam.z) * (1 + extention_twindow)
    intensity = np.empty(num)
    p_duration_m = p_duration_s * speed_of_light

    class radfileparams:
        offset = temporal_delay
        offset2 = start
        z = np.linspace(start, stop, num)
        for i in range(num):
            intensity[i] = p_intensity * math.exp(-(z[i] - offset)**2 / (2 * p_duration_m**2))
            if intensity[i] < 1e-50:
                intensity[i] = 0
        prad0 = intensity

    rad = radfileparams()
    f = open(out_file, 'w')
    f.write(rad_file_str2(rad))
    f.close()


def read_rad_file(filePath):
    '''
    reads analyticall radiation seed file from Genesis?
    '''
    beam = GenesisBeam()

    f = open(filePath, 'r')
    null = f.readline()
    for line in f:
        tokens = line.strip().split()

        if len(tokens) < 2:
            continue

        # print tokens[0:2]

        if tokens[0] == "?" and tokens[1] == "COLUMNS":
            beam.columns = tokens[2:]
            for col in beam.columns:
                beam.column_values[col] = []

            print (beam.columns)

        if tokens[0] != "?":
            # print tokens
            for i in range(0, len(tokens)):
                beam.column_values[beam.columns[i]].append(float(tokens[i]))

    #print (beam.columns)

    beam.z = beam.column_values['ZPOS']
    beam.prad0 = np.array(beam.column_values['PRAD0'])

    return beam


def adapt_rad_file(beam=None, rad_file=None, out_file='tmp.rad'):

    rad = read_rad_file(rad_file)
    rad.prad0 = np.interp(beam.z, rad.z, rad.prad0)
    rad.z = beam.z

    # print rad.z[301]
    # print beam.z[301]
    open(out_file, 'w').write(rad_file_str(rad))


def rad_file_str(rad):
    # header = "# \n? VERSION = 1.0\n? SIZE ="+str(len(rad.column_values['ZPOS']))+"\n? COLUMNS ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY XBEAM YBEAM PXBEAM PYBEAM ALPHAX ALPHAY CURPEAK ELOSS\n"
    header = "# \n? VERSION = 1.0\n? SIZE =" + str(len(rad.z)) + "\n? COLUMNS"
    # ZPOS GAMMA0 DELGAM EMITX EMITY BETAX BETAY XBEAM YBEAM PXBEAM PYBEAM ALPHAX ALPHAY CURPEAK ELOSS\n"
    for col in rad.columns:
        header = header + " " + col
    header += "\n"

    f_str = header

    rad.column_values['ZPOS'] = rad.z
    rad.column_values['PRAD0'] = rad.prad0

    for i in range(len(rad.z)):
        for col in rad.columns:
            buf = str(rad.column_values[col][i])
            f_str = f_str + buf + ' '
        f_str = f_str.rstrip() + '\n'

    return f_str


'''
    OTHER
'''


def generate_lattice(lattice, unit=1.0, energy=None, debug=False, min_phsh = False):

    print ('generating lattice file...')

    lat = '# header is included\n? VERSION= 1.00  including new format\n? UNITLENGTH= ' + str(unit) + ' :unit length in header\n'
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
    
    gamma = energy / m_e_GeV

    for e in lattice.sequence:

        l = float(e.l)

        # print e.type, pos, prevPos
        if e.__class__ == Undulator:

            l = float(e.nperiods) * float(e.lperiod)

            undLat += 'AW' + '    ' + str(e.Kx * np.sqrt(0.5)) + '   ' + str(round(l / unit, 2)) + '  ' + str(round((pos - prevPos - prevLen) / unit, 2)) + '\n'

            if debug:
                print ('added und ' + 'pos=' + str(pos) + ' prevPos=' + str(prevPos) + ' prevLen=' + str(prevLen))

            if prevLen > 0:
                #drifts.append([str( (pos - prevPos ) / unit ), str(prevLen / unit)])
                if debug:
                    print ('appending drift' + str((prevLen) / unit))
                L = pos - prevPos - prevLen #intersection length [m]
                K_rms = e.Kx * np.sqrt(0.5)
                
                if min_phsh:
                    xlamds = e.lperiod * (1 + K_rms**2) / (2 * gamma**2)
                    slip=(L / gamma**2) / 2 #free space radiation slippage [m]
                    add_slip = xlamds - slip % xlamds #free-space slippage to compensate with undulator K to bring it to integer number of wavelengths
                    K_rms_add = sqrt(2 * add_slip * gamma**2 / L) #compensational K
                    # driftLat += 'AD' + '    ' + str(e.Kx * np.sqrt(0.5)) + '   ' + str(round((pos - prevPos - prevLen) / unit, 2)) + '  ' + str(round(prevLen / unit, 2)) + '\n'
                    driftLat += 'AD' + '    ' + str(K_rms_add) + '   ' + str(round((L) / unit, 2)) + '  ' + str(round(prevLen / unit, 2)) + '\n'
                else:
                    driftLat += 'AD' + '    ' + str(K_rms) + '   ' + str(round((L) / unit, 2)) + '  ' + str(round(prevLen / unit, 2)) + '\n'
            
            prevPos = pos
            prevLen = l

        elif e.__class__ in [RBend, SBend, Drift]:
            pass

        elif e.__class__ == Quadrupole:
            #k = energy/0.2998 * float(e.k1) *  ( e.l / unit - int(e.l / unit) )
            # k = float(energy) * float(e.k1) / e.l #*  (1 +  e.l / unit - int(e.l / unit) )
            # k = float(energy) * float(e.k1) * 0.2998 / e.l #*  (1 +  e.l / unit - int(e.l / unit) )
            k = float(energy) * float(e.k1) / speed_of_light * 1e9
            if debug:
                print ('DEBUG' + str(e.k1) + ' ' + str(k) + ' ' + str(energy))
            quadLat += 'QF' + '    ' + str(k) + '   ' + str(round(e.l / unit, 2)) + '  ' + str(round((pos - prevPosQ - prevLenQ) / unit, 2)) + '\n'
            prevPosQ = pos
            prevLenQ = l
            # pass

        pos = pos + l

    return lat + undLat + driftLat + quadLat


def next_run_id(dir='.'):
    run_ids = []
    for f in os.listdir(dir):
        if f.startswith('run_'):
            run_ids.append(int(f.replace('run_', '')))

    if len(run_ids) > 0:
        return np.max(run_ids) + 1
    return 0


''' 
   standrard post-processing functions

'''

'''
    remove?
'''







def get_spectrum(power, phase, smax=1.0):
    ''' pulse spectrum in eV '''

    spec = np.fft.fft(np.sqrt(power) * np.exp(1.j * phase))
    spec = np.sqrt(spec * np.conj(spec))
    spec = np.real(np.roll(spec, len(spec) / 2))

    tmax = smax / speed_of_light
    freq = h_eV_s * (np.arange(1, len(spec) + 1) / tmax - 0.5 * len(spec) / tmax)
    return freq, spec

# def get_spectrum_n(power,phase, smax = 1.0):
    # spec = abs(np.fft.fft(np.sqrt(np.array(power)) * np.exp( 1.j* np.array(out.phase) ) , axis=0))**2
    # spec = spec / sqrt(out.nSlices)/(2*out.leng/out('ncar'))**2/1e10
    # xlamds=smax /
    # e_0=1239.8/out('xlamds')/1e9
    # out.freq_ev = h_eV_s * np.fft.fftfreq(len(out.spec), d=out('zsep') * out('xlamds') / speed_of_light)+e_0# d=out.dt


def get_power_exit(g):

    xlamds = g('xlamds')
    zsep = g('zsep')

    power = np.zeros(len(g.sliceValues.keys()))
    #phase = np.zeros(len(g.sliceValues.keys()))

    for i in g.sliceValues.keys():
        power[i - 1] = g.sliceValues[i]['power'][-1]
        #phase[i-1] = g.sliceValues[i]['phi_mid'][iZ]

    t = 1.0e+15 * zsep * xlamds / speed_of_light * np.arange(0, len(power))

    return power, t


def get_power_z(g):
    #nslice = int(g('nslice'))
    nslice = len(g.sliceValues.keys())
    nz = len(g.sliceValues[g.sliceValues.keys()[0]]['power'])
    power_z = np.zeros(nz)
    for i in range(nz):
        for j in range(nslice):
            power_z[i] += g.sliceValues[g.sliceValues.keys()[j]]['power'][i]

    return power_z / nslice












def transform_beam_file(beam_file=None, out_file='tmp.beam', s=None, transform=[[25.0, 0.1], [21.0, -0.1]], energy_scale=1, energy_new=None, emit_scale=1, n_interp=None):
    if beam_file.__class__ == str:
        beam = read_beam_file(beam_file)
    elif beam_file.__class__ == GenesisBeam or beam_file.__class__ == Beam:
        beam = beam_file
    else:
        print('Wrong beam input!')



    zmax, Imax = peaks(beam.z, beam.I, n=1)
    if s == None:
        idx = np.where(beam.z == zmax)[0][0]
        beam.idx_max = idx
    else:
        idx = np.where(beam.z > s)[0][0]
        beam.idx_max = idx
    print ('matching to slice ' + str(idx))

    #if plot: plot_beam(plt.figure(), beam)

    if transform:

        # if 'alphax' not in beam.keys() or if 'alphay' not in beam.keys():

        print ('transforming')
        g1x = np.matrix([[beam.betax[idx], beam.alphax[idx]],
                         [beam.alphax[idx], (1 + beam.alphax[idx]**2) / beam.betax[idx]]])

        g1y = np.matrix([[beam.betay[idx], beam.alphay[idx]],
                         [beam.alphay[idx], (1 + beam.alphay[idx]**2) / beam.betay[idx]]])

        b2x = transform[0][0]
        a2x = transform[0][1]

        b2y = transform[1][0]
        a2y = transform[1][1]

        g2x = np.matrix([[b2x, a2x],
                         [a2x, (1 + a2x**2) / b2x]])

        g2y = np.matrix([[b2y, a2y],
                         [a2y, (1 + a2y**2) / b2y]])

        Mix, Mx = find_transform(g1x, g2x)
        Miy, My = find_transform(g1y, g2y)

        # print Mi

        betax_new = []
        alphax_new = []

        betay_new = []
        alphay_new = []

        x_new = []
        px_new = []

        y_new = []
        py_new = []

        for i in range(len(beam.z)):
            g1x = np.matrix([[beam.betax[i], beam.alphax[i]],
                             [beam.alphax[i], (1 + beam.alphax[i]**2) / beam.betax[i]]])

            gx = Mix.T * g1x * Mix

            g1y = np.matrix([[beam.betay[i], beam.alphay[i]],
                             [beam.alphay[i], (1 + beam.alphay[i]**2) / beam.betay[i]]])

            gy = Miy.T * g1y * Miy

            # print i, gx[0,1], g1x[0,1]

            betax_new.append(gx[0, 0])
            alphax_new.append(gx[0, 1])

            betay_new.append(gy[0, 0])
            alphay_new.append(gy[0, 1])

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

        # print betax_new
        beam_new = Beam()
        beam_new.column_values = beam.column_values
        beam_new.columns = beam.columns
        # for parm in ['filePath']:
        # if hasattr(beam,parm):
        # setattr(beam_new,parm,getattr(beam,parm))
        beam_new.filePath = beam.filePath

        if energy_new != None:
            gamma_new = energy_new / m_e_GeV
            energy_scale = gamma_new / np.mean(np.array(beam.g0))

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

            try:
                beam_new.x = np.array(beam.x) * 0
                beam_new.px = np.array(beam.px) * 0
                beam_new.y = np.array(beam.y) * 0
                beam_new.py = np.array(beam.py) * 0
            except:
                pass
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

        # if plot:
        #    plot_beam(plt.figure(), beam_new)
        #    plt.show()

    if transform != None:
        beam_new.f_str = beam_file_str(beam_new)

    return beam_new


def test_beam_transform(beta1=10.0, alpha1=-0.1, beta2=20, alpha2=2.2):

    ex = 1.0

    g1 = np.matrix([[beta1, alpha1],
                    [alpha1, (1 + alpha1**2) / beta1]])

    g2 = np.matrix([[beta2, alpha2],
                    [alpha2, (1 + alpha2**2) / beta2]])

    Mi, M = find_transform(g1, g2)

    g = Mi.T * g1 * Mi

    print ('g1=' + str(g1))
    print ('g2=' + str(g2))
    print ('g=' + str(g))

    x = []
    xp = []

    x2 = []
    xp2 = []

    x3 = []
    xp3 = []

    for i in range(5000):
        x_, xp_ = gaussFromTwiss(ex, beta1, alpha1)
        x.append(x_)
        xp.append(xp_)

        u = M * np.matrix([[x_, xp_]]).T

        x3.append(u[0, 0])
        xp3.append(u[1, 0])

        x_, xp_ = gaussFromTwiss(ex, beta2, alpha2)
        x2.append(x_)
        xp2.append(xp_)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    plt.plot(x, xp, 'r.', alpha=0.2)
    plt.plot(x2, xp2, 'g.', alpha=0.2)
    plt.plot(x3, xp3, 'b.', alpha=0.2)
    plt.grid()

    l1, u1 = np.linalg.eig(g1)
    l2, u2 = np.linalg.eig(g2)

    plt.plot([0, u1[0, 0] * np.sqrt(l1[0])], [0, u1[1, 0] * np.sqrt(l1[0])], 'b--', lw=3)
    plt.plot([0, u1[0, 1] * np.sqrt(l1[1])], [0, u1[1, 1] * np.sqrt(l1[1])], 'b--', lw=3)

    plt.plot([0, u2[0, 0] * np.sqrt(l2[0])], [0, u2[1, 0] * np.sqrt(l2[0])], 'b--', lw=3)
    plt.plot([0, u2[0, 1] * np.sqrt(l2[1])], [0, u2[1, 1] * np.sqrt(l2[1])], 'b--', lw=3)










    z1 = M * u1[:, 0] * sqrt(l1[0])
    z2 = M * u1[:, 1] * sqrt(l1[1])

    plt.plot([0, z1[0]], [0, z1[1]], color='#000000', lw=5, alpha=0.2)
    plt.plot([0, z2[0]], [0, z2[1]], color='#000000', lw=5, alpha=0.2)


    plt.show()

    # 49.8131287015 1.12127199531 39.9184728466 -0.897874127701


#import argparse
#parser = argparse.ArgumentParser(description='Data plotting program.')
#parser.add_argument('--type', choices=["trajectory","intensity", "spectrum"], help='sum the integers (default: find the max)')
#parser.add_argument('input', help='input file')
#args = parser.parse_args()


def read_astra_dist(fileName):
    '''
    reading astra distribution parameters Parameter x y z px py pz clock macro_charge particle_index status_flag
    with units m m m eV/c eV/c eV/c ns nC
    returns numpy array?
    '''
    
    adist = np.loadtxt(fileName)
    adist[1:] = adist[1:] + adist[0]  # normalze to reference particle located at 1-st line
    return adist


def astra2edist(adist, center=1):
    '''
    converts astra particle distribution into GenesisElectronDist() object
    center - centers the distribution transversely
    '''
    edist = GenesisElectronDist()
    edist.x = adist[:, 0]
    edist.y = adist[:, 1]
    edist.t = (adist[:, 2] - np.mean(adist[:, 2])) / speed_of_light  # long position normalized to 0 and converted to time
    edist.xp = adist[:, 3] / adist[:, 5]  # angle of particles in x
    edist.yp = adist[:, 4] / adist[:, 5]  # angle of particles in y
    p_tot = sqrt(adist[:, 3]**2 + adist[:, 4]**2 + adist[:, 5]**2)
    edist.g = p_tot / m_e_eV  # energy to Gamma
    edist.part_charge = abs(adist[0][7]) * 1e-9  # charge of particle from nC

    if center:
        edist.x -= np.mean(edist.x)
        edist.y -= np.mean(edist.y)
        edist.xp -= np.mean(edist.xp)
        edist.yp -= np.mean(edist.yp)
    return edist


def astra2edist_ext(fileName_in, fileName_out='', center=1):
    if fileName_out == '':
        fileName_out = fileName_in + '.edist'
    adist = read_astra_dist(fileName_in)
    edist = astra2edist(adist, center=center)
    write_edist_file(edist, fileName_out, debug=0)


def filename_from_path(path_string):
    # return path_string[-path_string[::-1].find(os.path.sep)::]
    return path_string.split(os.path.sep)[-1]


def rematch_edist(edist, tws):

    from numpy import mean

    betax_n = tws.beta_x
    betay_n = tws.beta_y
    alphax_n = tws.alpha_x
    alphay_n = tws.alpha_y

    edist_out = deepcopy(edist)
    edist_out = edist_out.center()

    x = edist_out.x
    y = edist_out.y
    xp = edist_out.xp
    yp = edist_out.yp

    mean_x2 = mean(x**2)
    mean_y2 = mean(y**2)
    mean_px2 = mean(xp**2)
    mean_py2 = mean(yp**2)
    mean_xpx = mean(x * xp)
    mean_ypy = mean(y * yp)
    mean_g = mean(edist_out.g)

    emitx = mean_g * (mean_x2 * mean_px2 - mean_xpx**2)**0.5
    emity = mean_g * (mean_y2 * mean_py2 - mean_ypy**2)**0.5
    betax = mean_g * mean_x2 / emitx
    betay = mean_g * mean_y2 / emity
    alphax = -mean_g * mean_xpx / emitx
    alphay = -mean_g * mean_ypy / emity

    # remove correlation
    xp = xp + x * alphax / betax
    yp = yp + y * alphay / betay

    # scale beam
    x = x * sqrt(betax_n / betax)
    y = y * sqrt(betay_n / betay)
    xp = xp * sqrt(betax / betax_n)
    yp = yp * sqrt(betay / betay_n)

    # add new correlation
    xp = xp - alphax_n * x / betax_n
    yp = yp - alphay_n * y / betay_n

    edist_out.x = x
    edist_out.y = y
    edist_out.xp = xp
    edist_out.yp = yp

    return edist_out


def cut_lattice(lat, n_cells, elem_in_cell=4):
    '''
    reads MagneticLattice()
    returns MagneticLattice() without first n_cells*elem_in_cell elements
    '''
    n_cells=np.ceil(n_cells).astype(np.int)
    lat_new = deepcopy(lat)
    del lat_new.sequence[0:elem_in_cell * (n_cells)]
    return lat_new
    # returns lattice with #cells elements removed



'''
Scheduled for removal
'''


def getAverageUndulatorParameter(lattice, unit=1.0, energy=17.5):
    positions = sorted(lattice.lattice.keys())

    prevPos = 0

    ks = []
    ls = []

    for pos in positions:
        if lattice.elements[lattice.lattice[pos].id].type == 'undulator':
            e = lattice.elements[lattice.lattice[pos].id]
            l = float(e.params['nperiods']) * float(e.params['lperiod'])

            ks.append(float(e.params['K']))
            ls.append(l / unit)

            #lat += 'AW' +'    '+ e.params['K'] + '   ' + str( l  / unit ) + '  ' + str( (pos - prevPos) / unit ) + '\n'

            # if prevPos>0:
            #    drifts.append([str( (pos - prevPos ) / unit ), str(prevLen / unit)])

            #prevPos = pos + l
            #prevLen = l

    return np.mean(ks)
