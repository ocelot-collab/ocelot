
'''
interface to genesis
'''

import socket
import errno
import math

from ocelot.rad.fel import *
from ocelot.cpbd.beam import * # Twiss, Beam, gauss_from_twiss, ParticleArray
from ocelot.cpbd.elements import *
from ocelot.utils.launcher import *
from ocelot.common.math_op import *
from ocelot.optics.wave import *
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.rad.undulator_params import UndulatorParameters
# from ocelot.optics.utils import calc_ph_sp_dens
from ocelot.common.ocelog import *

_logger = logging.getLogger(__name__) 

_inputTemplate = "\
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
        self.suffix = ''

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
        self.zrayl = 2.5  # The Rayleigh length of the seeding radiation field.
        self.zwaist = 2.5  # Position of the waist of the seeding radiation field with respect to the undulator entrance.

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
        self.lout = [1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
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
        
        self.inp_txt = _inputTemplate
        
        self.int_vals = ('npart', 'nbins','ncar', 'zsep', 'nslice', 'ntail') #continue
        
    def input(self):
        
        inp_txt = deepcopy(self.inp_txt)

        if self.type == 'steady':
            # inp_txt = inp_txt.replace("__SHOTNOISE__", "itdp  =    0")
            inp_txt = inp_txt.replace("__ITDP__", "itdp = 0")
        else:
            # inp_txt = inp_txt.replace("__SHOTNOISE__", "shotnoise=  1")
            inp_txt = inp_txt.replace("__ITDP__", "itdp = 1")
            # self.prad0 = 0

        if self.beamfile != None:
            inp_txt = inp_txt.replace("__BEAMFILE__", " beamfile  =  '" + str(self.beamfile) + "'")
        else:
            inp_txt = inp_txt.replace("__BEAMFILE__\n", "")

        if self.fieldfile != None:
            inp_txt = inp_txt.replace("__FIELDFILE__", " fieldfile  =  '" + str(self.fieldfile) + "'")
        else:
            inp_txt = inp_txt.replace("__FIELDFILE__\n", "")

        if self.partfile != None:
            inp_txt = inp_txt.replace("__PARTFILE__", " partfile  =  '" + str(self.partfile) + "'")
        else:
            inp_txt = inp_txt.replace("__PARTFILE__\n", "")

        if self.edistfile != None:
            inp_txt = inp_txt.replace("__DISTFILE__", " distfile  =  '" + str(self.edistfile) + "'")
        else:
            inp_txt = inp_txt.replace("__DISTFILE__\n", "")

        if self.outputfile != None:
            inp_txt = inp_txt.replace("__OUTPUTFILE__", " outputfile  =  '" + str(self.outputfile) + "'")
        else:
            inp_txt = inp_txt.replace("__OUTPUTFILE__", " outputfile ='run.__RUNID__.gout'")

        if self.radfile != None:
            inp_txt = inp_txt.replace("__RADFILE__", " radfile  =  '" + str(self.radfile) + "'")
        else:
            inp_txt = inp_txt.replace("__RADFILE__\n", "")

        if self.magin == 0:
            inp_txt = inp_txt.replace("__MAGFILE__\n", "")
        else:
            inp_txt = inp_txt.replace("__MAGFILE__", " maginfile ='" + str(self.latticefile) + "'")

        # if self.trama == 1:
            # inp_txt = inp_txt.replace("__TRAMA__\n", "")
        # else:
            # inp_txt = inp_txt.replace("__TRAMA__\n", "")

        for p in self.__dict__.keys():
            if p in self.int_vals:
                val = int(self.__dict__[p])
            else:
                val = self.__dict__[p]
            inp_txt = inp_txt.replace("__" + str(p).upper() + "__", str(val).replace('[', '').replace(']', '').replace(',', ''))

        return inp_txt

    def __getattr__(self, name):
        if name not in self.__dict__.keys():
            return None
        else:
            return self.__dict__[name]

    def copy(self, inp, params):
        '''
        copies list of parameters from another GenesisInput() or GenesisOutput() object
        '''
        _logger.debug('copying input parameters')
        if inp.__class__ is GenesisInput:
            for param in params:
                if np.size(param) == 2:
                    param_r, param_w = param
                else:
                    param_r = param_w = param
                if hasattr(inp, param_r):
                    value = getattr(inp, param_r)
                    setattr(self, param_w, value)
                else:
                    _logger.warning(ind_str + 'could not copy ' + param_r)
        
        if inp.__class__ is GenesisOutput:
            for param in params:
                if np.size(param) == 2:
                    param_r, param_w = param
                else:
                    param_r = param_w = param
                    
                if inp(param_r) is not None:
                    value = inp(param_r)
                    setattr(self, param_w, value)
                else:
                    _logger.warning(ind_str + 'could not copy ' + param)
        
    def copymesh(self, inp, expt=()):
        
        # if inp.__class__ is GenesisInput:
        params = ('npart', 'nbins', 'xlamds', 'ncar', 'dgrid', 'zsep', 'nslice', 'ntail')
        _logger.debug('copying mesh parameters: ' + str(params))
        # elif inp.__class__ is GenesisOutput:
            # params = ('npart', 'nbins', 'xlamds', 'ncar', 'dgrid', 'zsep', 'nslice', 'ntail') #!no nslice in exceptions
        
        
        # if hasattr(inp, 'nslice'):
            # nslice = getattr(inp, 'nslice')
            # if nslice == 0:
                # print('Warning, nslice=0')
        
        # if hasattr(inp, 'ndcut'):
            # ndcut = getattr(inp, 'ndcut')
            # if ndcut == 0:
                # print('Warning, ndcut=0')
        
        
        
        params_exc = list( set(params).difference( set(expt) ) )
        
        
        if inp.__class__ is GenesisInput:
            for param in ['dgrid', 'ndcut', 'nslice']:
                if param in params_exc + ['ndcut']:
                    if hasattr(inp, param):
                        value = getattr(inp, param)
                        if value == 0:
                            _logger.info(ind_str + 'warning, %s=0' %(param))
                
        elif inp.__class__ is GenesisOutput:
            if ('dgrid' in params_exc) and (inp('meshsize') is not None):
                self.dgrid = inp('meshsize')*(inp.ncar-1) / 2
                params_exc.remove('dgrid')
            if ('nslice' in params_exc) and (inp('history_records') is not None):
                self.nslice = int(inp('history_records'))
                params_exc.remove('nslice')
            if ('ndcut' in params_exc) and hasattr(inp, 'ndcut'):
                value = getattr(inp, 'ndcut')
                if value == 0:
                    _logger.warning((ind_str + 'warning, {:}=0'.format('ndcut')))
        
        self.copy(inp, params_exc)



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

    def calc_spec(self, mode='mid', npad=0):
        '''
        calculates the on-axis spectrum at every position along the undulator and writes it into "spec" attirube
        
        if mode = "mid" then on-axis power with on-axis phases is used for calculation
        if mode = "int" then transversely integrated power with on-axis phases is used (strictly speaking inaccurate, but informative)
        npad (integer) if > 0 pads the power with zeros in order to increase resolution of spectrum.
        '''
        if self.nSlices == 1:
            raise AssertionError('Cannot calculate spectrum from steady-state simulation')
        
        if (npad%1 != 0) or npad < 0:
            raise ValueError('npad should be positive integer')
        
        if mode == 'ff':
            try:
                power = self.far_field
            except AttributeError:
                mode = 'mid'
                
        if mode == 'mid':
            power = self.p_mid
        elif mode == 'int':
            power = self.p_int
        elif mode == 'ff':
            pass
        else:
            raise ValueError('mode should be either "mid" or "int"')
            
        _logger.debug('calculating spectrum')
        power = power / (2 * self.leng / self('ncar'))**2
        phi_mid = self.phi_mid
        
        zeros = np.zeros((self.nSlices * npad, self.nZ))
        power = np.vstack((power, zeros))
        phi_mid = np.vstack((phi_mid, zeros))
        
        spec = abs(np.fft.fft(np.sqrt(np.array(power)) * np.exp(1.j * np.array(phi_mid)), axis=0))**2 * self.dt**2 * 1e10
        e_0 = h_eV_s * speed_of_light / self('xlamds')
        freq_ev = h_eV_s * np.fft.fftfreq(len(spec), d=self('zsep') * self('xlamds') * self('ishsty') / speed_of_light) + e_0
        
        spec = np.fft.fftshift(spec, axes=0)
        freq_ev = np.fft.fftshift(freq_ev, axes=0)
        freq_lamd = h_eV_s * speed_of_light * 1e9 / freq_ev
        
        self.spec = spec
        self.freq_ev = freq_ev
        self.freq_lamd = freq_lamd
        self.spec_mode = mode
        self.sliceKeys_used.append('spec')
        
        sum_spec = np.sum(self.spec, axis=0)
        sum_spec[sum_spec == 0] = np.inf
        
        self.freq_ev_mean = np.sum(self.freq_ev[:,np.newaxis] * self.spec, axis=0) / sum_spec
        self.freq_ev_mean[self.freq_ev_mean == 0] = np.inf
        
        self.n_photons = self.pulse_energy / q_e / self.freq_ev_mean
        self.spec_phot_density = calc_ph_sp_dens(self.spec, self.freq_ev, self.n_photons)
        # self.spec_phot_density = self.spec #tmp
        self.sliceKeys_used.append('spec_phot_density')
        # print ('        done')
        
    def phase_fix(self, wav=None, s=None, **kwargs):
        '''
        the way to display the phase, without constant slope caused by different radiation wavelength from xlamds. phase is set to 0 at maximum power slice
        '''
        _logger.debug('rewrapping phase')
        
        if 'spec' not in self.sliceKeys_used:
            raise AssertionError('first spectrum should be calculated')
        
        self.phi_mid_disp = deepcopy(self.phi_mid)
        
        if 'phen' in kwargs:
            wav = (h_eV_s * speed_of_light) / kwargs['phen'] * 1e9
        if wav == None:
            spec_idx = np.argmax(self.spec[:, -1])
        else:
            spec_idx = find_nearest_idx(self.freq_lamd, wav)
        
        _logger.debug(ind_str + 'from {}m to {}m carrier'.format(self('xlamds'), wav))
        
        if s == None:
            pow_idx = np.argmax(self.p_mid[:, -1])
        else:
            pow_idx = find_nearest_idx(self.s,s)
            
        for zi in range(np.shape(self.phi_mid_disp)[1]):
            # if debug > 1:
                # print ('      fixing phase display: ' + str(zi) + ' of ' + str(range(shape(self.phi_mid_disp)[1])))

            maxspectrum_wavelength = self.freq_lamd[spec_idx] * 1e-9
            phase = np.unwrap(self.phi_mid[:, zi])
            phase_cor = np.arange(self.nSlices) * (maxspectrum_wavelength - self('xlamds')) / self('xlamds') * self('zsep') * 2 * pi
            phase_fixed = phase + phase_cor
            phase_fixed -= phase_fixed[pow_idx]
            n = 1
            phase_fixed = (phase_fixed + n * pi) % (2 * n * pi) - n * pi
            self.phi_mid_disp[:, zi] = phase_fixed
        self.sliceKeys_used.append('phi_mid_disp')
        # print ('        done')
        
    def calc_radsize(self, weigh_transv=1):
        '''
        weigh_transv = True to average the transverse radiation size over slices with radiation power as a weight
        '''
        if weigh_transv and self.nSlices != 1:
            _logger.debug('calculating the weighted transverse radiation size')
            if np.amax(self.power) > 0:
                weight = self.power + np.amin(self.power[self.power != 0]) / 1e6
            else:
                weight = np.ones_like(self.power)
            self.rad_t_size_weighted = np.average(self.r_size * 1e6, weights=weight, axis=0)
            self.sliceKeys_used.append('rad_t_size_weighted')
        # print ('        done')
        
    def wig(self,z=np.inf,*args,**kwargs):
        return wigner_out(self, z=z, method='mp', *args, **kwargs)
        
    def re_read(self, read_level=2):
        return read_out_file(self.filePath, read_level=read_level)


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

    def center(self, s='com'):
        _logger.info('centering edist for s = '+str(s))
        if isinstance(s, str):
            if s == 'com': #center of mass
                
                mean_x = np.mean(self.x)
                mean_y = np.mean(self.y)
                mean_xp = np.mean(self.xp)
                mean_yp = np.mean(self.yp)
                
            else:
                _logger.error('unknown s string value')
                return None
        else:
            beam = edist2beam(self, 1e-7)
            beam_s = beam.get_s(s)
            
            mean_x = beam_s.x
            mean_y = beam_s.y
            mean_xp = beam_s.xp
            mean_yp = beam_s.yp
            
        self.x -= mean_x
        self.y -= mean_y
        self.xp -= mean_xp
        self.yp -= mean_yp
        
        _logger.debug(ind_str + 'mean_x  {}'.format(mean_x))
        _logger.debug(ind_str + 'mean_y  {}'.format(mean_y))
        _logger.debug(ind_str + 'mean_xp {}'.format(mean_xp))
        _logger.debug(ind_str + 'mean_yp {}'.format(mean_yp))
            
            
        
        # # edist_out = deepcopy(self)
        # edist_out.x -= np.mean(edist_out.x)
        # edist_out.y -= np.mean(edist_out.y)
        # edist_out.xp -= np.mean(edist_out.xp)
        # edist_out.yp -= np.mean(edist_out.yp)
        # return edist_out
        
    @property
    def s(self):
        return self.t * speed_of_light

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
    
    _logger.info('converting parray to edist')
    
    edist = GenesisElectronDist()
    
    e0 = p_array.E * 1e9 #[eV]
    p0 = np.sqrt( (e0**2 - m_e_eV**2) / speed_of_light**2 )
    
    p_oc = p_array.rparticles[5] # deltaE / average_impulse / speed_of_light
    edist.g = (p_oc * p0 * speed_of_light + e0) / m_e_eV
    edist.x = p_array.rparticles[0]  # position in x in meters
    edist.y = p_array.rparticles[2]  # position in y in meters
    edist.xp = p_array.rparticles[1]  # divergence in x
    edist.yp = p_array.rparticles[3]  # divergence in y
    edist.t = -1 * p_array.rparticles[4] / speed_of_light  # longitudinal position in seconds
    
    edist.part_charge = p_array.q_array[0] #fix for general case  # charge per particle
    edist.filePath = ''
        
    return edist
    
def edist2parray(edist):

    _logger.info('converting edist to parray')
    
    p_array = ParticleArray()
    p_array.rparticles = np.zeros((6,edist.len()))
    p_array.q_array = np.ones(edist.len()) * edist.part_charge
    
    g0 = np.mean(edist.g) # average gamma
    e0 = g0 * m_e_eV
    p0 = np.sqrt(g0**2-1) * m_e_eV / speed_of_light # average impulse
#    p0 = np.sqrt( (e0**2 - m_e_eV**2) / speed_of_light**2 ) # average impulse
    p_array.E = g0 * m_e_GeV # average energy in GeV
    
    p_array.rparticles[0] = edist.x # position in x in meters
    p_array.rparticles[1] = edist.xp  # divergence in x
    p_array.rparticles[2] = edist.y # position in y in meters
    p_array.rparticles[3] = edist.yp  # divergence in y
    p_array.rparticles[4] = -1 * edist.t * speed_of_light
    p_array.rparticles[5] = (edist.g - g0) * m_e_eV / p0 / speed_of_light
    
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

        # tws.emit_x= np.sqrt(tws.x*tws.pxpx-tws.xpx**2)
        # tws.emit_y= np.sqrt(tws.y*tws.pypy-tws.ypy**2)
        # tws.beta_x=tws.x/tws.emit_x
        # tws.beta_y=tws.y/tws.emit_y
        # tws.alpha_x=-tws.xpx/tws.emit_x
        # tws.alpha_y=-tws.ypy/tws.emit_y

        # return tws


# class GenesisBeam():
    # '''
    # Genesis analytical radiation input files storage object?
    # '''

    # def __init__(self):
        # self.columns = []
        # self.column_values = {}
        # self.fileName = ''
        # # self.filePath=''

    # def fileName(self):
        # return filename_from_path(self.filePath)

    # def len(self):
        # return len(self.z)

    # def idx_max_refresh(self):
        # self.idx_max = np.argmax(self.I)

    # def __delitem__(self, indarr):
        # self.z = np.delete(self.z, indarr)
        # if hasattr(self, 'I'):
            # self.I = np.delete(self.I, indarr)
        # if hasattr(self, 'g0'):
            # self.g0 = np.delete(self.g0, indarr)
        # if hasattr(self, 'dg'):
            # self.dg = np.delete(self.dg, indarr)
        # if hasattr(self, 'x'):
            # self.x = np.delete(self.x, indarr)
        # if hasattr(self, 'y'):
            # self.y = np.delete(self.y, indarr)
        # if hasattr(self, 'px'):
            # self.px = np.delete(self.px, indarr)
        # if hasattr(self, 'py'):
            # self.py = np.delete(self.py, indarr)
        # if hasattr(self, 'ex'):
            # self.ex = np.delete(self.ex, indarr)
        # if hasattr(self, 'ey'):
            # self.ey = np.delete(self.ey, indarr)
        # if hasattr(self, 'betax'):
            # self.betax = np.delete(self.betax, indarr)
        # if hasattr(self, 'betay'):
            # self.betay = np.delete(self.betay, indarr)
        # if hasattr(self, 'alphax'):
            # self.alphax = np.delete(self.alphax, indarr)
        # if hasattr(self, 'alphay'):
            # self.alphay = np.delete(self.alphay, indarr)
        # return self
    
    # def getitemold(self,index):
        
        # b_slice = deepcopy(self)
        # if index > b_slice.len():
            # raise IndexError('slice index out of range')
        
        # b_slice.z = b_slice.z[index]
        # if hasattr(b_slice, 'I'):
            # b_slice.I = b_slice.I[index]
        # if hasattr(b_slice, 'g0'):
            # b_slice.g0 = b_slice.g0[index]
        # if hasattr(b_slice, 'dg'):
            # b_slice.dg = b_slice.dg[index]
        # if hasattr(b_slice, 'x'):
            # b_slice.x = b_slice.x[index]
        # if hasattr(b_slice, 'y'):
            # b_slice.y = b_slice.y[index]
        # if hasattr(b_slice, 'px'):
            # b_slice.px = b_slice.px[index]
        # if hasattr(b_slice, 'py'):
            # b_slice.py = b_slice.py[index]
        # if hasattr(b_slice, 'ex'):
            # b_slice.ex = b_slice.ex[index]
        # if hasattr(b_slice, 'ey'):
            # b_slice.ey = b_slice.ey[index]
        # if hasattr(b_slice, 'betax'):
            # b_slice.betax = b_slice.betax[index]
        # if hasattr(b_slice, 'betay'):
            # b_slice.betay = b_slice.betay[index]
        # if hasattr(b_slice, 'alphax'):
            # b_slice.alphax = b_slice.alphax[index]
        # if hasattr(b_slice, 'alphay'):
            # b_slice.alphay = b_slice.alphay[index]
        
        # return b_slice


    # def __getitem__(self,index):
        
        # b_slice = deepcopy(self)
        # if index > b_slice.len():
            # raise IndexError('slice index out of range')
        
        # beam_slice = GenesisBeam()
        # l = self.len()
        # for attr in dir(self):
            # if attr.startswith('__'):
                # continue
            # value = getattr(self,attr)
            # if np.size(value) == l:
                # setattr(beam_slice,attr,value[index])
            # else:
                # setattr(beam_slice,attr,value)
        
        # return beam_slice


class GenesisRad():
    '''
    Genesis analytical radiation input files storage object?
    '''

    def __init__(self):
        self.columns = []
        self.column_values = {}




'''
    Genesis control
'''



def run_genesis(inp, launcher, read_level=2, assembly_ver='pyt', dfl_slipage_incl = True, min_phsh = False, debug=1):
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
    assembly_ver      - version of the assembly script: 'sys' - system based, 'pyt' - python based, None - assembly with Genesis assumed
    '''
    # import traceback
    _logger.info('starting genesis v2 preparation')
    # _logger.warning(len(traceback.extract_stack()))
    
    # create experimental directory
    if inp.run_dir == None and inp.exp_dir == None:
        raise ValueError('run_dir and exp_dir are not specified!')

    if inp.run_dir == None:
        if inp.exp_dir[-1]!=os.path.sep:
            inp.exp_dir+=os.path.sep
        inp.run_dir = inp.exp_dir + 'run_' + str(inp.runid) + '/'

    try:
        os.makedirs(inp.run_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(inp.run_dir):
            pass
        else:
            raise

    if inp.stageid is None:
        inp_path = inp.run_dir + 'run.' + str(inp.runid) + str(inp.suffix) + '.inp'
        out_path = inp.run_dir + 'run.' + str(inp.runid) + str(inp.suffix) + '.gout'
        inp.stageid = ''
        stage_string = ''
    else:
        inp_path = inp.run_dir + 'run.' + str(inp.runid) + '.s' + str(inp.stageid) + str(inp.suffix) + '.inp'
        out_path = inp.run_dir + 'run.' + str(inp.runid) + '.s' + str(inp.stageid) + str(inp.suffix) + '.gout'
        stage_string = '.s' + str(inp.stageid)

    inp_file = filename_from_path(inp_path)
    out_file = filename_from_path(out_path)

    # cleaning directory
    _logger.debug(ind_str + 'removing old files')
    os.system('rm -rf ' + inp.run_dir + 'run.' + str(inp.runid) + stage_string + str(inp.suffix) + '*')  # to make sure all stage files are cleaned
    # os.system('rm -rf ' + out_path+'*') # to make sure out files are cleaned
    # os.system('rm -rf ' + inp_path+'*') # to make sure inp files are cleaned
    os.system('rm -rf ' + inp.run_dir + 'tmp.cmd')

    # create and fill necessary input files
    if inp.latticefile == None:
        if inp.lat != None:
            _logger.debug(ind_str + 'writing ' + inp_file + '.lat')
            if not hasattr(inp, 'lat_unit'):
                lat_unit = inp.xlamd
                _logger.debug(2*ind_str + 'lat_unit_size = xlamds = {} m'.format(lat_unit))
            else:
                if inp.lat_unit is None:
                    lat_unit = inp.xlamd
                    _logger.debug(2*ind_str + 'lat_unit_size = xlamds = {} m'.format(lat_unit))
                else:
                    lat_unit = inp.lat_unit
                    _logger.debug(2*ind_str + 'lat_unit_size = {} m'.format(lat_unit))
            open(inp_path + '.lat', 'w').write(generate_lattice(inp.lat, unit=lat_unit, energy=inp.gamma0 * m_e_GeV, debug = debug, min_phsh = min_phsh))
            inp.latticefile = inp_file + '.lat'

    if inp.beamfile == None:
        if inp.beam != None:
            _logger.debug(ind_str + 'writing ' + inp_file + '.beam')
            open(inp_path + '.beam', 'w').write(beam_file_str(inp.beam))
            inp.beamfile = inp_file + '.beam'

    if inp.edistfile == None:
        if inp.edist != None:
            _logger.debug(ind_str + 'writing ' + inp_file + '.edist')
            write_edist_file(inp.edist, inp_path + '.edist', debug=1)
            inp.edistfile = inp_file + '.edist'

    if inp.partfile == None:
        if inp.dpa != None:
            _logger.debug(ind_str + 'writing ' + inp_file + '.dpa')
            # print ('!!!!!!! no write_particle_file() function')
            write_dpa_file(inp.dpa, inp_path + '.dpa', debug=1)
            inp.partfile = inp_file + '.dpa'

    if inp.fieldfile == None:
        if inp.dfl != None:
            _logger.debug(ind_str + 'writing ' + inp_file + '.dfl')
            write_dfl_file(inp.dfl, inp_path + '.dfl', debug=1)
            inp.fieldfile = inp_file + '.dfl'

    if inp.radfile == None:
        if inp.rad != None:
            _logger.debug(ind_str + 'writing ' + inp_file + '.rad')
            open(inp_path + '.rad', 'w').write(rad_file_str(inp.rad))
            inp.radfile = inp_file + '.rad'

    if inp.outputfile == None:
        inp.outputfile = out_file
    _logger.debug(ind_str + 'writing ' + inp_file)
    open(inp_path, 'w').write(inp.input())
    open(inp.run_dir + 'tmp.cmd', 'w').write(inp_file + '\n')

    launcher.dir = inp.run_dir
    _logger.debug(ind_str + 'preparing launcher')
    launcher.prepare()
    # _logger.debug()
    # RUNNING GENESIS ###
    genesis_time = time.time()
    launcher.launch()
    _logger.info(ind_str + 'genesis simulation time %.2f seconds' % (time.time() - genesis_time))
    # RUNNING GENESIS ###

    if assembly_ver is not None:
        # genesis output slices assembly
        _logger.info(ind_str + 'assembling slices')
        _logger.debug(2 * ind_str + 'assembly_ver = {}'.format(assembly_ver))
        
        assembly_time = time.time()

        
        if assembly_ver == 'sys':

            _logger.info(2 * ind_str + 'assembling *.out file')
            start_time = time.time()
            os.system('cat ' + out_path + '.slice* >> ' + out_path)
            os.system('rm ' + out_path + '.slice* 2>/dev/null')
            _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

            _logger.info(2 * ind_str + 'assembling *.dfl file')
            start_time = time.time()
            if dfl_slipage_incl:
                os.system('cat ' + out_path + '.dfl.slice*  >> ' + out_path + '.dfl.tmp')
                #bytes=os.path.getsize(out_path +'.dfl.tmp')
                command = 'dd if=' + out_path + '.dfl.tmp of=' + out_path + '.dfl conv=notrunc conv=notrunc 2>/dev/null'# obs='+str(bytes)+' skip=1
                os.system(command)
            else:
                os.system('cat ' + out_path + '.dfl.slice*  > ' + out_path + '.dfl')
            os.system('rm ' + out_path + '.dfl.slice* 2>/dev/null')
            os.system('rm ' + out_path + '.dfl.tmp 2>/dev/null')
            _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

            _logger.info(2 * ind_str + 'assembling *.dpa file')
            start_time = time.time()
            os.system('cat ' + out_path + '.dpa.slice* >> ' + out_path + '.dpa')
            os.system('rm ' + out_path + '.dpa.slice* 2>/dev/null')
            _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))
            
            _logger.info(2 * ind_str + 'assembling *.fld file')
            start_time = time.time()
            if dfl_slipage_incl:
                os.system('cat ' + out_path + '.fld.slice*  >> ' + out_path + '.fld.tmp')
                #bytes=os.path.getsize(out_path +'.fld.tmp')
                command = 'dd if=' + out_path + '.fld.tmp of=' + out_path + '.fld conv=notrunc conv=notrunc 2>/dev/null'# obs='+str(bytes)+' skip=1
                os.system(command)
            else:
                os.system('cat ' + out_path + '.fld.slice*  > ' + out_path + '.fld')
            os.system('rm ' + out_path + '.fld.slice* 2>/dev/null')
            os.system('rm ' + out_path + '.fld.tmp 2>/dev/null')
            _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

            _logger.info(2 * ind_str + 'assembling *.par file')
            start_time = time.time()
            os.system('cat ' + out_path + '.par.slice* >> ' + out_path + '.par')
            os.system('rm ' + out_path + '.par.slice* 2>/dev/null')
            _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

        elif assembly_ver == 'pyt':
            # there is a bug with dfl assembly
            import glob
            ram = 1

            _logger.info(2 * ind_str + 'assembling *.out file')
            start_time = time.time()
            assemble(out_path, ram=ram, debug=debug)
            os.system('rm ' + out_path + '.slice* 2>/dev/null')
            _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))
            
            for i in range(10):        #read all possible harmonics (up to 10 now)
                ii=str(i)
                if ii=='0': ii=''
                if os.path.isfile(str(out_path + '.dfl' + ii)):
                    _logger.info(2 * ind_str + 'assembling *.dfl'+ii+' file')
                    start_time = time.time()
                    assemble(out_path + '.dfl'+ii, overwrite=dfl_slipage_incl, ram=ram, debug=debug)
                    os.system('rm ' + out_path + '.dfl'+ii+'.slice* 2>/dev/null')
                    os.system('rm ' + out_path + '.dfl'+ii+'.tmp 2>/dev/null')
                    _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

            if os.path.isfile(str(out_path + '.dpa')):
                _logger.info(2 * ind_str + 'assembling *.dpa file')
                start_time = time.time()
                assemble(out_path + '.dpa', ram=ram, debug=debug)
                os.system('rm ' + out_path + '.dpa.slice* 2>/dev/null')
                _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))
                    
            if os.path.isfile(str(out_path + '.fld')):
                _logger.info(2 * ind_str + 'assembling *.fld file')
                start_time = time.time()
                assemble(out_path + '.fld', overwrite=dfl_slipage_incl, ram=ram, debug=debug)
                os.system('rm ' + out_path + '.fld.slice* 2>/dev/null')
                os.system('rm ' + out_path + '.fld.tmp 2>/dev/null')
                _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

            if os.path.isfile(str(out_path + '.par')):
                _logger.info(2 * ind_str + 'assembling *.par file')
                start_time = time.time()
                assemble(out_path + '.par', ram=ram, debug=debug)
                os.system('rm ' + out_path + '.par.slice* 2>/dev/null')
                _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

        else:
            # raise ValueError('assembly_ver should be either "sys" or "pyt"')
            pass
        _logger.debug(2 * ind_str + 'assembly time %.2f seconds' % (time.time() - assembly_time))
        
    
    
    if read_level >= 0:
        out = read_out_file(out_path, read_level=read_level)
        _logger.debug(ind_str + 'done, time %.2f seconds' % (time.time() - assembly_time))
        return out
    else:
        _logger.debug(ind_str + 'done, time %.2f seconds' % (time.time() - assembly_time))
        return None


def assemble(fileName, remove=1, overwrite=0, ram=1, debug=1):
    '''
    assembles the fileName.slice* files into fileName
    remove - delete *.slice* files
    overwrite - writes *.slice* files on top of fileName.slice* starting from the beginning. Applicable for genesis dfl file assembly
    ram - store the *.slice* files in ram simultaneously
    '''
    import glob
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
        _logger.debug('reading ' + str(N) + ' slices to RAM...')
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
        _logger.debug('writing...')
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
            if debug > 1:
                tot = size(fins)
                _logger.log(5, ind_str + 'slice {} of {}'.format(i,tot))
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
    
    _logger.debug('creating run_n subdirectories')
    
    if exp_dir[-1]!=os.path.sep:
        exp_dir+=os.path.sep
    for run_id in run_ids:

        try:
            run_dir = exp_dir + 'run_' + str(run_id) + '/'
            _logger.log(5, ind_str + run_dir)
            os.makedirs(run_dir)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(run_dir):
                pass
            else:
                raise

    try:
        res_dir = exp_dir + 'results'
        _logger.log(5, ind_str + res_dir)
        os.makedirs(res_dir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(res_dir):
            pass
        else:
            raise


def generate_input(undulator, beam, E_photon = None, itdp=True, *args, **kwargs):
    '''
    Create Genesis inp object with default input parameters
    '''
    inp = GenesisInput()
    
    _logger.info('generating genesis2 input')
    # Next line added by GG 27.05.2016: it was in script
    

    #beam.emit_xn, beam.emit_yn = beam.emit[beam.C]
    #beam.gamma_rel = beam.E / (0.511e-3)
    #beam.emit_x = beam.emit_xn / beam.gamma_rel
    #beam.emit_y = beam.emit_yn / beam.gamma_rel
    
    if undulator.__class__ == UndulatorParameters:
        inp.xlamd = undulator.lw
        inp.aw0 = undulator.K / np.sqrt(2)
    elif undulator.__class__ == MagneticLattice:
        # from ocelot.cpbd.elements import Undulator
        lat = undulator
        indx_u = lat.find_indices(Undulator)
        und = [lat.sequence[i] for i in indx_u if np.any(lat.sequence[i].Kx != 0)][0] #first undulator can be opened
        inp.xlamd = und.lperiod
        inp.aw0 = np.mean(und.Kx) / np.sqrt(2)
        inp.lat = lat
    
    
    if beam.len() > 1:
        inp.beam = beam
        tpulse = np.abs(beam.s[-1] - beam.s[0]) / speed_of_light * 1e15 / 6 #fix to fwhm
        beam = beam.pk()
        beam.tpulse = tpulse
    
    inp.magin = 1
    #inp.zstop = 50
    #inp.nsec = 20


    inp.awd = inp.aw0
    inp.delgam = beam.sigma_E / m_e_GeV
    inp.gamma0 = beam.E / m_e_GeV
    
    inp.betax = beam.beta_x
    inp.betay = beam.beta_y
    
    # inp.rxbeam = np.sqrt(beam.emit_x * beam.beta_x)
    # inp.rybeam = np.sqrt(beam.emit_y * beam.beta_y)

    inp.alphax = beam.alpha_x
    inp.alphay = beam.alpha_y

    inp.curpeak = beam.I

    inp.xbeam = beam.x
    inp.ybeam = beam.y
    inp.pxbeam = beam.xp
    inp.pybeam = beam.yp

    inp.emitx = beam.emit_xn
    inp.emity = beam.emit_yn
    
    # idx_0 = inp.beam.I>0
    rxbeam = np.nanmax(np.sqrt(inp.beam.beta_x * inp.beam.emit_x))
    rybeam = np.nanmax(np.sqrt(inp.beam.beta_y * inp.beam.emit_y))
    inp.dgrid = np.nanmax([rxbeam, rybeam]) * 8 #due to bug in Genesis2 that crashes when electrons leave the mesh

    inp.hn=1 # should be flexible in the future
    
    felParameters = calculateFelParameters(inp)

    if E_photon is None:
        inp.xlamds = felParameters.lambda0
    else:
        inp.xlamds = h_eV_s * speed_of_light / E_photon
    
    if itdp:
        inp.prad0 = 0
    else:
        inp.prad0 = felParameters.P_sn
    
    inp.fbess0 = felParameters.fc
    inp.zrayl = felParameters.zr

    if itdp:
        inp.type = "tdp"
        
        inp.ipseed = 99999
        inp.ncar = 151
        inp.zsep = int(math.ceil(0.25 / (4 * pi * felParameters.rho3)))  # 0.25 is the additional factor to be "on the safe side"
        
        if not hasattr(inp, 'beam'):
            inp.curlen = beam.tpulse * speed_of_light / 1e15
            inp.nslice = 8 * int(inp.curlen / inp.zsep / inp.xlamds)
        else:
            inp.curlen = 0
            inp.nslice = 0
    
    inp.ntail = 0
    inp.npart = 4096
    inp.rmax0 = 9
    inp.delz = 1
    inp.felParameters = felParameters
    # print out FEL parameter estimates
    # printFelParameters(inp)
    return inp


def get_genesis_launcher(launcher_program=None, launcher_argument=''):
    '''
    Returns MpiLauncher() object for given program
    '''
    host = socket.gethostname()

    launcher = MpiLauncher()
    if launcher_program != None:
        launcher.program = launcher_program
    else:
        if host.startswith('max'):
            launcher.program = '/data/netapp/xfel/products/genesis/genesis'
            launcher.argument = ' < tmp.cmd | tee log'
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
                    3 = additional attributed are calculated, like spectrum, transversely averaged radiation size, etc.
    precision - precision of stored values precision
    debug -     0 = no messages printed in console
                1 = basic info and execution time is printed
                2 = most detailed info is printed (real debug)
    '''
    import re
    out = GenesisOutput()
    out.filePath = filePath
    # out.fileName = filename_from_path(filePath)

    _logger.info('reading gen2 {} file'.format(os.path.basename(filePath)))
    _logger.debug(ind_str + 'reading from ' + filePath)
#    print '        - reading from ', fileName

    chunk = ''
    output_unsorted = []
    nSlice = 0

    wait_attempt = 3
    wait_time = 0.5
    while os.path.isfile(out.filePath) != True:
        _logger.warning(ind_str + 'waiting for "' + out.fileName() + '" ' + str(wait_time) + 's [' + str(wait_attempt) + ']')
        time.sleep(wait_time)  # wait for the .out file to be assembled
        wait_attempt -= 1
        if wait_attempt == 0:
            _logger.error(ind_str + 'file "' + out.filePath + '" not found')
            raise IOError('File ' + out.filePath + ' not found')

    if os.path.getsize(out.filePath) == 0:
        _logger.error(ind_str + 'file "' + out.filePath + '" has zero size')
        raise IOError('File ' + out.filePath + ' has zero size')
    
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
            _logger.log(5, ind_str + 'reading slice # ' + str(nSlice))

        if tokens[0] == 'power':
            chunk = 'slice'
            if len(out.sliceKeys) == 0:  # to record the first instance
                out.sliceKeys = list(copy(tokens))
                _logger.debug(ind_str + 'reading slice values ')
            continue

        if tokens[0] == '$newrun':
            chunk = 'input1'
            _logger.debug(ind_str + 'reading input parameters')
            continue

        if tokens[0] == '$end':
            chunk = 'input2'
            continue

        if tokens == ['z[m]', 'aw', 'qfld']:
            chunk = 'magnetic optics'
            _logger.debug(ind_str + 'reading magnetic optics ')
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
        if chunk == 'slice' and read_level >= 2:

            # tokens_fixed=re.sub(r'([0-9])\-([0-9])',r'\g<1>E-\g<2>',' '.join(tokens))
            # tokens_fixed=re.sub(r'([0-9])\+([0-9])',r'\g<1>E+\g<2>',tokens_fixed)
            # tokens=tokens_fixed.split()
            try:
                vals = list(map(precision, tokens))
            except ValueError:
                _logger.log(5, ind_str + 'wrong E value, fixing')
                _logger.log(5, ind_str + str(tokens))
                tokens_fixed = re.sub(r'([0-9])\-([0-9])', r'\g<1>E-\g<2>', ' '.join(tokens))
                tokens_fixed = re.sub(r'([0-9])\+([0-9])', r'\g<1>E+\g<2>', tokens_fixed)
                tokens_fixed = tokens_fixed.split()
                _logger.log(5, ind_str + str(tokens_fixed))
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

    #check for consistency
    if chunk == '':
        _logger.error(ind_str + 'File "' + out.filePath + '" has no genesis output information or is corrupted')
        raise ValueError('File "' + out.filePath + '" has no genesis output information or is corrupted')
    
    for parm in ['z', 'aw', 'qfld', 'I', 'n']:
        exec('out.' + parm + ' = np.array(out.' + parm + ')')

    if out('dgrid') == 0:
        rbeam = np.sqrt(out('rxbeam')**2 + out('rybeam')**2)
        ray = np.sqrt(out('zrayl') * out('xlamds') / np.pi * (1 + (out('zwaist') / out('zrayl'))**2))
        out.leng = out('rmax0') * (rbeam + ray)
    else:
        out.leng = 2 * out('dgrid')
    out.ncar = int(out('ncar'))  # number of mesh points
    
    #universal solution?
    out.leng=out('meshsize')*(out.ncar-1)
    
    
    if read_level == 0:
        _logger.debug(ind_str + 'read_level=0, returning *.out header')
        _logger.debug(ind_str + 'done in %.2f seconds' % (time.time() - start_time))
        return out

    out.nSlices = len(out.n)
    if out('entries_per_record') is None:
        _logger.error(ind_str + 'In file "' + out.filePath + '" file header is missing')
        raise ValueError('In file "' + out.filePath + '" file header is missing')
    
    out.nZ = int(out('entries_per_record'))  # number of records along the undulator
    _logger.debug(ind_str + 'nSlices ' + str(out.nSlices))
    _logger.debug(ind_str + 'nZ ' + str(out.nZ))

    if nSlice == 0:
        _logger.error(ind_str + 'In file "' + out.filePath + '" number of recorded slices is zero')
        raise ValueError('In file "' + out.filePath + '" number of recorded slices is zero')

    n_missing = (out.n[-1] - out.n[0]) - (len(out.n) - 1) * out('ishsty')
    if n_missing != 0:
        _logger.error(ind_str + 'File "' + out.filePath + '" is missing at least ' + str(n_missing) + ' slices')
        raise ValueError('File "' + out.filePath + '" is missing at least ' + str(n_missing) + ' slices')
    
    if read_level >= 2:
        output_unsorted = np.array(output_unsorted)  # .astype(precision)
        _logger.debug('output_unsorted.shape = ' + str(output_unsorted.shape))
        _logger.debug(ind_str + 'out.sliceKeys' + str(out.sliceKeys))
        for i in range(len(out.sliceKeys)):
            key = out.sliceKeys[int(i)]
            _logger.debug(ind_str + 'key = ' + str(key))
            if key[0].isdigit():
                hn = key[0]
                if 'bunch' in key:
                    key = 'bunching'
                elif 'phase' in key:
                    key = 'phi_mid'
                elif 'p-mid' in key:
                    key = 'p_mid'
                elif 'power' in key:
                    key = 'power'
                else:
                    pass
                key='h{:}_'.format(hn)+key
                _logger.debug(2*ind_str + 'key_new = ' + str(key))
                
            _logger.log(5, ind_str + 'assembling') 
            # _logger.debug('output_unsorted[:,' + str(i) + '].shape = ' + str(output_unsorted[:,i].shape))
            command = 'out.' + key.replace('-', '_').replace('<', '').replace('>', '') + ' = output_unsorted[:,' + str(i) + '].reshape((' + str(int(out.nSlices)) + ',' + str(int(out.nZ)) + '))'
            _logger.debug(ind_str + command)
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
        
        if read_level >= 2:
            out.pulse_energy = np.sum(out.power * out.dt, axis=0)
        
        if read_level >= 3:
            out.calc_spec(npad=2)
            # out.phase_fix()
            # out.calc_radsize(weigh_transv=1)
            
    else:
        out.s = [0]
    
    if out('iscan') != 0:
        out.scv = out.I  # scan value
        out.I = np.linspace(1, 1, len(out.scv))  # because used as a weight
    
    # tmp for back_compatibility
    if read_level >= 2:
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
    
    _logger.debug(ind_str + 'done in %.2f seconds' % (time.time() - start_time))
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
    _logger.info('reading stat genesis output')
    _logger.info(ind_str + 'proj_dir = {}'.format(proj_dir))
    _logger.debug(ind_str + 'stage = {}'.format(stage))
    _logger.debug(ind_str + 'run_inp = {}'.format(str(run_inp)))
    _logger.debug(ind_str + 'param_inp = {}'.format(str(param_inp)))
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
            # _logger.debug(ind_str + 'reading run {}'.format(irun))
            _logger.debug(ind_str + 'reading run {}'.format(irun))
            outlist[irun] = read_out_file(out_file, read_level=2, debug=debug)
            outlist[irun].calc_spec()
            run_range_good.append(irun)
            # except:
    run_range = run_range_good
    
    # check if all gout have the same number of slices nSlice and history records nZ
    for irun in run_range[1:]:
        if outlist[irun].nSlices != outlist[run_range[0]].nSlices or outlist[irun].nZ != outlist[run_range[0]].nZ:
            raise ValueError('Non-uniform out objects (run %s)' %(irun))
    
    _logger.debug(ind_str + 'good_run_range = {}'.format(str(run_range)))

    if param_inp == []:
        # if debug > 1:
            # print(outlist[run_range[0]].sliceKeys_used)
        param_range = outlist[run_range[0]].sliceKeys_used
    else:
        param_range = param_inp
    _logger.debug(ind_str + 'param_range = {}'.format(str(param_range)))

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
        elif np.ndim(param_matrix) == 2 and np.shape(param_matrix)[1] == outlist[irun].nZ:
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

    _logger.debug(ind_str + 'done in %.2f seconds' % (time.time() - start_time))
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
    _logger.info('reading stat genesis output')
    _logger.info(ind_str + 'file_tamplate = {}'.format(file_tamplate))
    _logger.debug(ind_str + 'run_inp = {}'.format(str(run_inp)))
    _logger.debug(ind_str + 'param_inp = {}'.format(str(param_inp)))
    start_time = time.time()
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
            outlist[irun].calc_spec(npad=1)
            # print(outlist[irun].freq_lamd[0])
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
            _logger.debug(ind_str + 'run {}'.format(irun))
            if not hasattr(outlist[irun], param):
                continue
            else:
                param_matrix.append(deepcopy(getattr(outlist[irun], param)))

        param_matrix = np.array(param_matrix)
        if np.ndim(param_matrix) == 3:
            param_matrix = np.swapaxes(param_matrix, 0, 2)
        elif np.ndim(param_matrix) == 2 and np.shape(param_matrix)[1] == outlist[irun].nZ:
            param_matrix = np.swapaxes(param_matrix, 0, 1)[:, np.newaxis, :]
        else:
            pass
        setattr(out_stat, param, param_matrix)

    out_stat.stage = None
    out_stat.dir = os.path.dirname(file_tamplate)
    out_stat.run = run_range
    out_stat.z = outlist[irun].z #check if all the same!
    out_stat.s = outlist[irun].s
    out_stat.f = outlist[irun].freq_lamd
    out_stat.t = outlist[irun].t
    out_stat.dt = outlist[irun].dt
    
    out_stat.xlamds=outlist[irun]('xlamds')
    # out_stat.filePath=proj_dir

    _logger.debug(ind_str + 'done in %.2f seconds' % (time.time() - start_time))
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
    _logger.info('reading radiation field file (dfl) from .out.dfl')
    _logger.debug(ind_str + 'opening handle ' + str(out))
    
    if os.path.isfile(str(out)):
        out = read_out_file(out, read_level=0, debug=0)
    if not isinstance(out, GenesisOutput):
        _logger.error('out is neither GenesisOutput() nor a valid path')
        raise ValueError('out is neither GenesisOutput() nor a valid path')

    if filePath is None:
        _logger.debug(ind_str + 'from out.filePath')
        filePath = out.filePath + '.dfl'
    else:
        _logger.debug(ind_str + 'from filepath')
    _logger.debug(2*ind_str + filePath)
    
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
    hist_rec - number of dfl records within a single file (*.fld case), not finished!
    '''
    
    _logger.info('reading radiation field file (dfl)')
    _logger.debug(ind_str + 'reading from ' + filePath)
    start_time = time.time()
    
    assert (Nxy % 1 == 0), 'Nxy nust be an integer'
    Nxy = int(Nxy)
    if not os.path.isfile(filePath):
        _logger.warning(ind_str + 'dfl file ' + filePath + ' not found')
        raise IOError('dfl file ' + filePath + ' not found !')
        
    b = np.fromfile(filePath, dtype=complex).astype(vartype)
    Nz = b.shape[0] / Nxy / Nxy / hist_rec
    assert (Nz % 1 == 0), 'Wrong Nxy or corrupted file'
    Nz = int(Nz)
    
    dfl = RadiationField()
    if hist_rec > 1:
        _logger.debug(ind_str + 'assuming {} history records'.format(hist_rec))
        dfl.fld = b.reshape(hist_rec, Nz, Nxy, Nxy)
        dfl.fld = np.rollaxis(dfl.fld, 0, 4)
    else:
        dfl.fld = b.reshape(Nz, Nxy, Nxy)
    dfl.dx = Lxy / dfl.Nx()
    dfl.dy = Lxy / dfl.Ny()
    dfl.dz = xlamds * zsep
    dfl.xlamds = xlamds
    dfl.filePath = filePath
    
    _logger.debug(ind_str + 'dfl shape = {}'.format(str(dfl.fld.shape)))
    _logger.debug(ind_str + 'done in %.2f sec' % (time.time() - start_time))

    return dfl



def write_dfl_file(dfl, filePath=None, debug=1):
    '''
    Function to write the RadiationField object into filePath file
    dfl RadiationField object
    filePath - path top write the file
        if None then filePath = dfl.filePath
    '''
    _logger.info('writing radiation field file (dfl)')
    _logger.debug(ind_str + 'filePath = {}'.format(filePath))
    start_time = time.time()
    
    if dfl.domains() != ('t', 's'):
        _logger.error("dfl.domains != ('t', 's'), but " + str(dfl.domains()))
        raise ValueError("dfl.domains != ('t', 's')")
    
    if filePath == None:
        filePath = dfl.filePath
        
    if dfl.__class__ != RadiationField:
        raise ValueError('wrong radiation object: should be RadiationField')

    d = dfl.fld.flatten().astype(np.complex128)
    d.tofile(filePath, format='complex')

    _logger.debug(ind_str + 'done in %.2f sec' % (time.time() - start_time))


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
    
    _logger.info('reading particle (dpa) file from .out.dpa')
    
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
    _logger.info('reading particle (dpa) file')
    _logger.debug(ind_str + 'filePath = {}'.format(filePath))
    
    if not os.path.isfile(filePath):
        _logger.warning(ind_str + 'dpa file ' + filePath + ' not found')
        raise IOError('dpa file ' + filePath + ' not found !')

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

    _logger.debug(ind_str + 'nslice = ' + str(nslice))
    _logger.debug(ind_str + 'npart = ' + str(npart))
    _logger.debug(ind_str + 'nbins = ' + str(nbins))
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

    _logger.debug(ind_str + 'done in %.2f sec' % (time.time() - start_time))

    return dpa
    
def write_dpa_file(dpa, filePath=None, debug=1):
    
    _logger.info('writing particle file')
    _logger.debug(ind_str + 'filePath = {}'.format(filePath))
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
    
    _logger.debug(ind_str + 'done in %.2f sec' % (time.time() - start_time))



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
    nbins = np.shape(dpa.ph)[1]
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
    Convert dpa to edist objects
    reads GenesisParticlesDump() object
    returns GenesisElectronDist() object
    num_part - desired approximate number of particles in edist
    smear - whether to shuffle macroparticles smearing microbunching
    '''
    import random
    start_time = time.time()
    _logger.info('transforming particle to distribution file')

    assert out('itdp') == True, '! steadystate Genesis simulation, dpa2dist() not implemented yet!'

    npart = int(out('npart'))
    # nslice=int(out('nslice'))
    nslice = int(out.nSlices * out('ishsty'))
    nbins = int(out('nbins'))
    xlamds = out('xlamds')
    zsep = int(out('zsep'))
    gen_I = out.I
    gen_t = out.t
    num_part = int(num_part)

    if (npart / nbins) % 1 != 0:
        raise ValueError('non-integer number of particles per bin')
    else:
        npart_per_bin = int(npart / nbins)

    m = np.arange(nslice)
    m = np.tile(m, (nbins, npart_per_bin, 1))
    m = np.rollaxis(m, 2, 0)
    # print('shape_m='+str(shape(m)))
    
    if smear:
        z = dpa.ph * xlamds / 2 / pi + m * xlamds * zsep + xlamds * zsep * (1 - np.random.random((nslice, nbins, npart_per_bin)))
    else:
        z = dpa.ph * xlamds / 2 / pi + m * xlamds * zsep

    t = np.array(z / speed_of_light)

    t_scale = np.linspace(0, nslice * zsep * xlamds / speed_of_light * 1e15, nslice)

    pick_n = np.interp(t_scale, gen_t, gen_I)
    _logger.debug(ind_str + 'sum pick_n = ' + str(sum(pick_n)))
    _logger.debug(ind_str + 'npart = ' + str(npart))
    _logger.debug(ind_str + 'num_part = ' + str(num_part))
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
            edist.t = np.append(edist.t, t[i, pick_i])
            edist.g = np.append(edist.g, e[i, pick_i])
            edist.x = np.append(edist.x, x[i, pick_i])
            edist.y = np.append(edist.y, y[i, pick_i])
            edist.xp = np.append(edist.xp, px[i, pick_i])
            edist.yp = np.append(edist.yp, py[i, pick_i])

    # edist.t = edist.t * (-1) + max(edist.t)
    edist.t -= edist.t.min()
    edist.xp /= edist.g
    edist.yp /= edist.g

    edist.x = np.flipud(edist.x)
    edist.y = np.flipud(edist.y)
    edist.xp = np.flipud(edist.xp)
    edist.yp = np.flipud(edist.yp)
    edist.t = np.flipud(edist.t)
    edist.g = np.flipud(edist.g)

    edist.part_charge = out.beam_charge / edist.len()
    _logger.debug(ind_str + 'edist.len() = ' + str(edist.len()))
    # edist.charge=out.beam_charge
    if hasattr(dpa, 'filePath'):
        edist.filePath = dpa.filePath + '.edist'
    _logger.debug(ind_str + 'edist.filePath = ' + edist.filePath)
    # print 'max_y_out', np.amax(t_out)
    # print 'e_out', np.amax(e_out),np.amin(e_out)

    _logger.debug(ind_str + 'done in %.2f sec' % (time.time() - start_time))

    return edist


'''
    EDIST
'''


def read_edist_file_out(out, debug=1):

    return read_dist_file(out.filePath + '.edist', debug=debug)


def read_edist_file(filePath, **kwargs):
    ### MODIFIED, NOT TESTED
    '''
    reads particle distribution file (distfile in genesis input)
    returns GenesisElectronDist() 
    '''
    edist = GenesisElectronDist()
    edist.filePath = filePath

    # if debug > 0:
    _logger.info('reading particle distribution file (edist)')
    _logger.debug(ind_str + 'filePath = {}'.format(filePath))
        # print ('    reading particle distribution file')

    try:
        f = open(filePath, 'r')
    except Exception:
        _logger.error(ind_str + 'edist file ' + filePath + ' not found !')
        raise
    # if not os.path.isfile(filePath):
        # if debug:
            # raise IOError('      ! edist file ' + filePath + ' not found !')
        # else:
            # _logger.error('      ! edist file ' + filePath + ' not found !')
    # else:
        # if debug > 1:
            # print ('        - reading from ' + filePath)

    start_time = time.time()
    dist_column_values = {}
    null = f.readline()
    for line in f:
        tokens = line.strip().split()

        if len(tokens) < 2:
            continue

        if tokens[0] == "?" and tokens[1] == "CHARGE":
            charge = float(tokens[3])
            _logger.debug('  Charge= '+str(charge))

        if tokens[0] == "?" and tokens[1] == "COLUMNS":
            dist_columns = tokens[2:]
            for col in dist_columns:
                dist_column_values[col] = []
            _logger.debug(ind_str + ''.join(str(i) + ' ' for i in dist_columns))

        if tokens[0] != "?":
            for i in range(0, len(tokens)):
                dist_column_values[dist_columns[i]].append(float(tokens[i]))
    
    if 'P' in dist_column_values.keys():
        edist.g = np.sqrt(np.array(dist_column_values.get('P'))**2 + 1)
    elif 'GAMMA' in dist_column_values.keys():
        edist.g = np.array(dist_column_values.get('GAMMA'))
    else:
        _logger.warning(ind_str + 'Neither P nor GAMMA found in edist')
    
    p = np.sqrt(edist.g**2 - 1)
    
    if 'XPRIME' in dist_column_values.keys():
        edist.xp = np.array(dist_column_values.get('XPRIME'))
    elif 'PX' in dist_column_values.keys():
        edist.xp = np.array(dist_column_values.get('PX')) / p
    else:
        _logger.warning(ind_str + 'Neither XPRIME nor PX found in edist')
        
    if 'YPRIME' in dist_column_values.keys():
        edist.yp = np.array(dist_column_values.get('YPRIME'))
    elif 'PY' in dist_column_values.keys():
        edist.yp = np.array(dist_column_values.get('PY')) / p
    else:
        _logger.warning(ind_str + 'Neither YPRIME nor PY found in edist')
        
    edist.x = np.array(dist_column_values.get('X'))
    edist.y = np.array(dist_column_values.get('Y'))
    # edist.y = np.array(dist_column_values['Y'])
    # edist.xp = np.array(dist_column_values['XPRIME'])
    # edist.yp = np.array(dist_column_values['YPRIME'])
    
    if 'T' in dist_column_values.keys():
        edist.t = -np.array(dist_column_values.get('T'))
    elif 'Z' in dist_column_values.keys():
        edist.t = np.array(dist_column_values.get('Z')) / speed_of_light
    else:
        _logger.warning(ind_str + 'Neither T nor Z found in edist')
    
    # edist.t = np.array(dist_column_values['T'])
    
    # edist.x = np.flipud(edist.x)
    # edist.y = np.flipud(edist.y)
    # edist.xp = np.flipud(edist.xp)
    # edist.yp = np.flipud(edist.yp)
    # edist.t = np.flipud(edist.t)
    # edist.g = np.flipud(edist.g)

    edist.part_charge = charge / edist.len()
    
    _logger.debug(ind_str + 'part_n = {:}'.format(edist.len()))
    _logger.debug(ind_str + 'part_charge = {:} C'.format(edist.part_charge))
    _logger.debug(ind_str + 'done in %.2f sec' % (time.time() - start_time))

    return edist


# from inspect import getouterframes, currentframe
# import os

def cut_edist_std(edist, all_std=None, x_std=4, y_std=4, xp_std=4, yp_std=4):
    _logger.info('cutting particle distribution by standard deviation')
    if all_std is not None:
        x_std = y_std = xp_std = yp_std = all_std
    
    x_mean = np.mean(edist.x)
    y_mean = np.mean(edist.y)
    xp_mean = np.mean(edist.xp)
    yp_mean = np.mean(edist.yp)
    _logger.debug(ind_str + 'X  mean {}'.format(x_mean))
    _logger.debug(ind_str + 'Y  mean {}'.format(y_mean))
    _logger.debug(ind_str + 'XP mean {}'.format(xp_mean))
    _logger.debug(ind_str + 'YP mean {}'.format(yp_mean))
    
    x_1std = np.std(edist.x)
    y_1std = np.std(edist.y)
    xp_1std = np.std(edist.xp)
    yp_1std = np.std(edist.yp)
    _logger.debug(ind_str + 'X  std {} * {} sigmas'.format(x_1std, x_std))
    _logger.debug(ind_str + 'Y  std {} * {} sigmas'.format(y_1std, y_std))
    _logger.debug(ind_str + 'XP std {} * {} sigmas'.format(xp_1std, xp_std))
    _logger.debug(ind_str + 'YP std {} * {} sigmas'.format(yp_1std, yp_std))
    
    x_lim = x_std * x_1std
    y_lim = y_std * y_1std
    xp_lim = xp_std * xp_1std
    yp_lim = yp_std * yp_1std
    
    return cut_edist(edist, x_lim=(x_mean-x_lim, x_mean+x_lim), y_lim=(y_mean-y_lim, y_mean+y_lim), xp_lim=(xp_mean-xp_lim, xp_mean+xp_lim), yp_lim=(yp_mean-yp_lim, yp_mean+yp_lim))
    
    

def cut_edist(edist,
              t_lim=(-np.inf, np.inf),
              g_lim=(-np.inf, np.inf),
              x_lim=(-np.inf, np.inf),
              xp_lim=(-np.inf, np.inf),
              y_lim=(-np.inf, np.inf),
              yp_lim=(-np.inf, np.inf), 
              s_lim=None, debug=1):
    '''
    cuts GenesisElectronDist() in phase space
    '''
    from numpy import logical_or

    _logger.info('cutting particle distribution file')
    start_time = time.time()

    _logger.debug(ind_str + 'T  lim {} : {} '.format(*t_lim))
    _logger.debug(ind_str + 'G  lim {} : {} '.format(*g_lim))
    _logger.debug(ind_str + 'X  lim {} : {} '.format(*x_lim))
    _logger.debug(ind_str + 'Y  lim {} : {} '.format(*y_lim))
    _logger.debug(ind_str + 'XP lim {} : {} '.format(*xp_lim))
    _logger.debug(ind_str + 'YP lim {} : {} '.format(*yp_lim))
    
    if s_lim is not None:
        index_t = logical_or(edist.t < s_lim[0]/speed_of_light, edist.t > s_lim[1]/speed_of_light)
    else:
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

    _logger.info(ind_str + '{:.2f} % cut'.format((edist.charge() - edist_f.charge()) / edist.charge() * 100))
    _logger.debug(ind_str + 'done in {:.2f} sec'.format(time.time() - start_time))

    return edist_f


def set_edist_energy(edist, E_GeV, debug=1):
    '''
    adds energy to the electron beam so that new average energy is E_GeV in [GeV]
    '''
    _logger.info('scaling particle distribution file energy to {} GeV'.format(E_GeV))
    if not isinstance(edist, GenesisElectronDist):
        raise ValueError('out is neither GenesisOutput() nor a valid path')
    
    edist_out = deepcopy(edist)
    edist_out.g -= np.mean(edist_out.g)
    edist_out.g += E_GeV * 1e9 / m_e_eV
    
    return edist_out
    
def disperse_edist(edist, R56, debug=1):
    '''
    Introduces dispersion (good for simulating weak chicanes)
    delays or advances time coordinate of the particles depending on ther energy with respect to the averaged energy
    '''
    _logger.info('introducing dispersion to particle distribution file with R56 '+ str(R56) + ' m')
    if not isinstance(edist, GenesisElectronDist):
        raise ValueError('out is neither GenesisOutput() nor a valid path')
    
    edist_out = deepcopy(edist)
    edist_out.t += R56 * (edist_out.g - np.mean(edist_out.g)) / edist_out.g / speed_of_light
    
    return edist_out
    
    
    
def repeat_edist(edist, repeats, smear=1e-3, not_smear=[]):
    '''
    dublicates the GenesisElectronDist() by given factor
    repeats  - the number of repetitions
    smear - smear new particles by x of global standard deviation of parameter
    '''
    
    _logger.info('repeating edist by factor of {}'.format(repeats))
    
    if not isinstance(edist, GenesisElectronDist):
        raise ValueError('out is neither GenesisOutput() nor a valid path')
        
    if repeats < 1:
        raise ValueError('repeats cannot be smaller that 1')
        
    edist_out = GenesisElectronDist()
    edist_out.filePath = edist.filePath

    edist_out.x = np.repeat(edist.x, repeats)
    edist_out.y = np.repeat(edist.y, repeats)
    edist_out.xp = np.repeat(edist.xp, repeats)
    edist_out.yp = np.repeat(edist.yp, repeats)
    edist_out.t = np.repeat(edist.t, repeats)
    edist_out.g = np.repeat(edist.g, repeats)
    edist_out.part_charge = edist.part_charge / repeats

    if smear:
        n_par = edist_out.len()
        smear_factor = 1e-3  # smear new particles by smear_factor of standard deviation of parameter
        
        for attr in ['x','y','xp','yp','t','g']:
            if attr not in not_smear:
                val = getattr(edist_out, attr)
                val += np.random.normal(scale=np.std(val) * smear_factor, size=n_par)
                setattr(edist_out, attr, val)
        # edist_out.x += np.random.normal(scale=np.std(edist_out.x) * smear_factor, size=n_par)
        # edist_out.y += np.random.normal(scale=np.std(edist_out.y) * smear_factor, size=n_par)
        # edist_out.xp += np.random.normal(scale=np.std(edist_out.xp) * smear_factor, size=n_par)
        # edist_out.yp += np.random.normal(scale=np.std(edist_out.yp) * smear_factor, size=n_par)
        # edist_out.t += np.random.normal(scale=np.std(edist_out.t) * smear_factor, size=n_par)
        # edist_out.g += np.random.normal(scale=np.std(edist_out.g) * smear_factor, size=n_par)

    return edist_out


def write_edist_file(edist, filePath=None, debug=1):
    '''
    writes GenesisElectronDist() into filePath folder
    '''
    # REQUIRES NUMPY 1.7
    # header='? VERSION = 1.0 \n? SIZE = %s \n? CHARGE = %E \n? COLUMNS X XPRIME Y YPRIME T P'%(len(edist.x),charge)
    # np.savetxt(filePath_write, np.c_[edist.x,edist.xp,edist.y,edist.yp,edist.t,edist.g],header=header,fmt="%E", newline='\n',comments='')

    _logger.info('writing particle distribution file')
    start_time = time.time()

    if filePath == None:
        filePath = edist.filePath
    
    _logger.debug(ind_str + 'filePath = {}'.format(filePath))

    header = '? VERSION = 1.0 \n? SIZE = %s \n? CHARGE = %E \n? COLUMNS X XPRIME Y YPRIME Z GAMMA\n' % (edist.len(), edist.charge())
    _logger.debug('  '+header)
    f = open(filePath, 'w')
    f.write(header)
    f.close()
    f = open(filePath, 'ab')
    np.savetxt(f, np.c_[edist.x, edist.xp, edist.y, edist.yp, edist.s, edist.g], fmt="%e", newline='\n')
    f.close()

    _logger.debug(ind_str + 'done in %.2f sec' % (time.time() - start_time))


def edist2beam(edist, step=2e-7): #check
    '''
    reads GenesisElectronDist()
    returns BeamArray()
    step [m] - long. size ob bin to calculate distribution parameters
    '''
    
    _logger.info('transforming edist to beamfile')
    
    part_c = edist.part_charge
    t_step = step / speed_of_light
    t_min = min(edist.t)
    t_max = max(edist.t)
    dist_t_window = t_max - t_min
    npoints = int(dist_t_window / t_step)
    t_step = dist_t_window / npoints
    beam = BeamArray(npoints-1)

    for i in range(npoints - 1):
        indices = (edist.t > t_min + t_step * i) * (edist.t < t_min + t_step * (i + 1))
        beam.s[i] = (t_min + t_step * (i + 0.5)) * speed_of_light
        # print(sum(indices))
        if np.sum(indices) > 2:
            dist_g = edist.g[indices]
            dist_E = dist_g * m_e_GeV
            dist_x = edist.x[indices]
            dist_y = edist.y[indices]
            dist_xp = edist.xp[indices]
            dist_yp = edist.yp[indices]
            dist_sigma_E = np.std(dist_g) * m_e_GeV
            dist_p = np.sqrt(dist_g**2 - 1)
            dist_px = dist_xp * dist_p
            dist_py = dist_yp * dist_p
            
            beam.I[i] = np.sum(indices) * part_c / t_step
            beam.E[i] = np.mean(dist_E)
            beam.sigma_E[i] = dist_sigma_E
            
            dist_x_m = np.mean(dist_x)
            dist_y_m = np.mean(dist_y)
            dist_xp_m = np.mean(dist_xp)
            dist_yp_m = np.mean(dist_yp)
            
            beam.x[i] = dist_x_m
            beam.y[i] = dist_y_m
            beam.xp[i] = dist_xp_m
            beam.yp[i] = dist_yp_m
            
            dist_x -= dist_x_m
            dist_y -= dist_y_m
            dist_xp -= dist_xp_m
            dist_yp -= dist_yp_m
            
            beam.emit_x[i] = (np.mean(dist_x**2) * np.mean(dist_xp**2) - np.mean(dist_x * dist_xp)**2)**0.5
            # if beam.ex[i]==0: beam.ey[i]=1e-10
            beam.emit_y[i] = (np.mean(dist_y**2) * np.mean(dist_yp**2) - np.mean(dist_y * dist_yp)**2)**0.5
            # if beam.ey[i]==0: beam.ey[i]=1e-10
            beam.beta_x[i] = np.mean(dist_x**2) / beam.emit_x[i]
            beam.beta_y[i] = np.mean(dist_y**2) / beam.emit_y[i]
            beam.alpha_x[i] = -np.mean(dist_x * dist_xp) / beam.emit_x[i]
            beam.alpha_y[i] = -np.mean(dist_y * dist_yp) / beam.emit_y[i]
    
    idx = np.where(np.logical_or.reduce((beam.I == 0, beam.g == 0)))
    del beam[idx]
    
    if hasattr(edist,'filePath'):
        beam.filePath = edist.filePath + '.beam'
    
    _logger.debug(ind_str + 'done')
    
    return(beam)


'''
    BEAM
'''


# def read_beam_file_out(out, debug=1):
    # return read_beam_file(out.filePath, debug=debug)

def read_beam_file(filePath, *args, **kwargs):
    '''
    reads beam file from filePath folder
    returns BeamArray()
    '''
    # import types
    
    _logger.info('reading beam file')
    start_time = time.time()

    beam = BeamArray()
    
    def fileName(self):
        try:
            str = filename_from_path(self.filePath)
        except:
            str = None
        return str
        
    # setattr(beam, 'fileName', classmethod(fileName))
    # beam.fileName = fileName
    # beam.fileName = types.MethodType( fileName, beam )
    beam.fileName = fileName.__get__(beam)
    
    columns = []
    column_values = {}

    f = open(filePath, 'r')
    null = f.readline()
    for line in f:
        tokens = line.strip().split()

        if len(tokens) < 2:
            continue

        # print tokens[0:2]

        if tokens[0] == "?" and tokens[1] == "COLUMNS":
            columns = tokens[2:]
            for col in columns:
                column_values[col] = []

            _logger.debug(ind_str + str(columns))

        if tokens[0] != "?":
            # print tokens
            for i in range(0, len(tokens)):
                column_values[columns[i]].append(float(tokens[i]))

    # print columns

    beam.s = np.array(column_values['ZPOS'])
    beam.I = np.array(column_values['CURPEAK'])

    for parm in [['g', 'GAMMA0'],
                 ['x', 'XBEAM'],
                 ['y', 'YBEAM'],
                 ['px', 'PXBEAM'],
                 ['py', 'PYBEAM'],
                 ['dg', 'DELGAM'],
                 ['emit_xn', 'EMITX'],
                 ['emit_yn', 'EMITY'],
                 ['beta_x', 'BETAX'],
                 ['beta_y', 'BETAY'],
                 ['alpha_x', 'ALPHAX'],
                 ['alpha_y', 'ALPHAY'],
                 ]:
        if parm[1] in column_values.keys():
            setattr(beam, parm[0], np.array(column_values[parm[1]]))
        else:
            setattr(beam, parm[0], np.zeros_like(beam.s))

    try:
        beam.eloss = np.array(column_values['ELOSS'])
    except:
        beam.eloss = np.zeros_like(beam.I)

    beam.filePath = filePath
    del(column_values, columns)

    _logger.debug(ind_str + 'done in %.2f sec' % (time.time() - start_time))

    return beam


def beam_file_str(beam):
    '''
    reads BeamArray()
    returns string of electron beam file, suitable for Genesis
    '''
    
    dict = {'s' : 'ZPOS',
             'I': 'CURPEAK',
             'emit_xn': 'EMITX',
             'emit_yn': 'EMITY',
             'beta_x': 'BETAX',
             'beta_y': 'BETAY',
             'alpha_x': 'ALPHAX',
             'alpha_y': 'ALPHAY',
             'x': 'XBEAM',
             'y': 'YBEAM',
             'px': 'PXBEAM',
             'py': 'PYBEAM',
             'g': 'GAMMA0',
             'dg': 'DELGAM',
             'eloss': 'ELOSS'}
    
    l = beam.len()
    attrs = beam.params() + list(beam.properties)
    
    f_str = "# \n? VERSION = 1.0\n? SIZE =" + str(l) + "\n? COLUMNS"
    
    for attr in attrs:
        if attr in dict:
            f_str = f_str + " " + dict[attr]
    f_str += "\n"

    # for parm in [['s', 'ZPOS'],
                 # ['I', 'CURPEAK'],
                 # ['emit_xn', 'EMITX'],
                 # ['emit_yn', 'EMITY'],
                 # ['beta_x', 'BETAX'],
                 # ['beta_y', 'BETAY'],
                 # ['alpha_x', 'ALPHAX'],
                 # ['alpha_y', 'ALPHAY'],
                 # ['x', 'XBEAM'],
                 # ['y', 'YBEAM'],
                 # ['px', 'PXBEAM'],
                 # ['py', 'PYBEAM'],
                 # ['g', 'GAMMA0'],
                 # ['dg', 'DELGAM'],
                 # ['eloss', 'ELOSS'],
                 # ]:
        # try:
            # beam.column_values[parm[1]] = getattr(beam, parm[0])
        # except:
            # pass

    for i in range(beam.len()):
        for attr in attrs:
            if attr in dict:
                buf = str(getattr(beam,attr)[i])
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


def zero_wake_at_ipk(beam):
    '''
    reads BeamArray()
    shifts the wake pforile so that 
    at maximum current slice wake is zero
    returns GenesisBeam()
    
    allows to account for wake losses without 
    additional linear undulator tapering
    '''
    if 'eloss' in beam.params():
        beam_new = deepcopy(beam)
        beam_new.eloss -= beam_new.eloss[beam_new.idx_max()]
        return beam_new
    else:
        return beam


def set_beam_energy(beam, E_GeV_new):
    '''
    reads BeamArray()
    returns BeamArray()
    sets the beam energy with peak current to E_GeV_new
    '''
    beam_new = deepcopy(beam)
    idx = beam_new.idx_max()
    beam_new.E += (E_GeV_new - beam_new.E[idx])
    
    return beam_new


def transform_beam_twiss(beam, transform=None, s=None):
    # transform = [[beta_x,alpha_x],[beta_y, alpha_y]] or Twiss()
    if transform == None:
        return beam
    else:
        if s == None:
            idx = beam.idx_max()
        else:
            idx = find_nearest_idx(beam.s, s)

        g1x = np.matrix([[beam.beta_x[idx], beam.alpha_x[idx]],
                         [beam.alpha_x[idx], (1 + beam.alpha_x[idx]**2) / beam.beta_x[idx]]])

        g1y = np.matrix([[beam.beta_y[idx], beam.alpha_y[idx]],
                         [beam.alpha_y[idx], (1 + beam.alpha_y[idx]**2) / beam.beta_y[idx]]])
        
        # temporary legacy support
        if transform.__class__.__name__ == 'Twiss':
            b2x = transform.beta_x
            b2y = transform.beta_y
            a2x = transform.alpha_x
            a2y = transform.alpha_y
        else:
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

        for i in range(beam.len()):
            g1x = np.matrix([[beam.beta_x[i], beam.alpha_x[i]],
                             [beam.alpha_x[i], (1 + beam.alpha_x[i]**2) / beam.beta_x[i]]])

            gx = Mix.T * g1x * Mix

            g1y = np.matrix([[beam.beta_y[i], beam.alpha_y[i]],
                             [beam.alpha_y[i], (1 + beam.alpha_y[i]**2) / beam.beta_y[i]]])

            gy = Miy.T * g1y * Miy

            # print i, gx[0,1], g1x[0,1]

            betax_new.append(gx[0, 0])
            alphax_new.append(gx[0, 1])
            betay_new.append(gy[0, 0])
            alphay_new.append(gy[0, 1])

        beam.beta_x = np.array(betax_new)
        beam.beta_y = np.array(betay_new)
        beam.alpha_x = np.array(alphax_new)
        beam.alpha_y = np.array(alphay_new)

def cut_beam(beam=None, cut_s=[-np.inf, np.inf]):
    '''
    cuts BeamArray() object longitudinally
    cut_z [m] - limits of the cut
    '''
    _logger.info('cutting beam file between {:.2e} and {:.2e}'.format(cut_s[0], cut_s[1]))
    beam_new = deepcopy(beam)
    beam_new.sort()
    idxl, idxr = find_nearest_idx(beam_new.s, cut_s[0]), find_nearest_idx(beam_new.s, cut_s[1])
    _logger.debug(ind_str + 'slice numbers {:.2e} and {:.2e}'.format(idxl, idxr))
    if idxl==idxr:
        _logger.warning(ind_str + 'slice numbers {:.2e} and {:.2e} are the same'.format(idxl, idxr))
    
    return beam_new[idxl:idxr]


def get_beam_s(beam, s=0):
    '''
    obtains values of the beam at s position
    '''
    idx = find_nearest_idx(beam.s,s)
    return(beam[idx])

def get_beam_peak(beam):
    '''
    obtains values of the beam at s position
    '''
    idx = beam.idx_max()
    return beam[idx]
    

def find_transform(g1, g2):
    '''
    find transform from twiss matrix g1 to g2: x -> M x, g -> Mi.T g M
    '''
    l1, u1 = np.linalg.eig(g1)
    l2, u2 = np.linalg.eig(g2)

    M1 = np.matrix([[u1[0, 0], u1[1, 0]],
                    [-u1[1, 0], u1[0, 0]]])

    d = np.sqrt(l1[0] / l2[0])

    M2 = np.matrix([[d, 0],
                    [0, 1.0 / d]])

    M2d = np.matrix([[1. / d, 0],
                     [0, d]])

    M3 = np.matrix([[u2[0, 0], -u2[1, 0]],
                    [u2[1, 0], u2[0, 0]]])

    return np.linalg.inv(M3 * M2 * M1), M3 * M2d * M1


def write_beam_file(filePath, beam, debug=0):
    _logger.info('writing beam file')
    start_time = time.time()
    fd = open(filePath, 'w')
    fd.write(beam_file_str(beam))
    fd.close()
    _logger.debug(ind_str + 'done in %.2f sec' % (time.time() - start_time))

# def write_beam_file(filePath, beam, debug=0):
    # if debug > 0:
        # print ('    writing beam file')
    # start_time = time.time()

    # fd = open(filePath, 'w')
    # fd.write(beam_file_str(beam))
    # fd.close()

    # if debug > 0:
        # print('      done in %.2f sec' % (time.time() - start_time))


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


def generate_lattice(lattice, unit=1.0, energy=None, debug=1, min_phsh = False):
    _logger.info('generating lattice file...')
    _logger.debug(ind_str + 'minimum phase shift = {:}'.format(min_phsh))

    lat = '# header is included\n? VERSION= 1.00  including new format\n? UNITLENGTH= ' + str(unit) + ' :unit length in header\n'
    undLat = ''
    quadLat = ''
    driftLat = ''

    prevPosU = 0
    prevLenU = 0
    prevPosQ = 0
    prevLenQ = 0
    prevPosD = 0
    prevLenD = 0
    
    L_inters = 0
    K_rms = 0
    K_inters = None
    phi_inters = None

    pos = 0

    drifts = []
    quads = []
    
    gamma = energy / m_e_GeV

    for e in lattice.sequence:
        l = e.l
        _logger.log(5, 'L from previous undulator = {}'.format(L_inters))
        _logger.log(5, 'element {} with length {}'.format(e.__class__, l))
        
        if e.__class__ == Undulator:
            
            if L_inters > 0:
                xlamds = l_period * (1 + np.mean(K_rms)**2) / (2 * gamma**2)
                _logger.log(5, 'xlamds = {}'.format(xlamds))
                _logger.log(5, 'gamma = {}'.format(gamma))
                slip_frsp=(L_inters / gamma**2) / 2 # free-space radiation slippage [m]
                phi_frsp = slip_frsp / xlamds * 2 * np.pi # free-space phase shift [m]
                _logger.log(5, 'free-space phase drift is {} rad'.format(phi_frsp))
                
                if K_inters is None:
                    if phi_inters is None:
                        _logger.log(5, 'no phase shifter values in Drift, assuming minimum to bring to n*2pi')
                        slip_add = xlamds - slip_frsp % xlamds #free-space slippage to compensate with undulator K to bring it to integer number of wavelengths
                    elif phi_inters > phi_frsp:
                        _logger.log(5, 'phi intersection {} rad> phi free space {} rad'.format(phi_inters, phi_frsp))
                        phi_add = phi_inters - phi_frsp
                        slip_add = phi_add * xlamds / 2 / np.pi
                        _logger.log(5, '   difference to add {} rad or {} m'.format(phi_add, slip_add))
                    elif phi_inters <= phi_frsp:
                        _logger.log(5, 'phi intersection {} rad < phi free space {} rad'.format(phi_inters, phi_frsp))
                        phi_diff = phi_frsp - phi_inters
                        period_diff = phi_diff // (2*np.pi)
                        phi_phsh1 = phi_inters + (period_diff+1) * 2 * np.pi
                        _logger.log(5, '   raising phi intersection to {} rad'.format(phi_phsh1))
                        phi_add = phi_phsh1 - phi_frsp
                        slip_add = phi_add * xlamds / 2 / np.pi
                        _logger.log(5, '   difference to add {} rad or {} m'.format(phi_add, slip_add))
                    K_rms_add = np.sqrt(2 * slip_add * gamma**2 / L_inters) #compensational K
                else:
                    K_rms_add = np.mean(K_inters)
                    _logger.log(5, '   compensational K {}'.format(K_rms_add))
                
                driftLat += 'AD' + '    ' + str(K_rms_add) + '   ' + str(round((L_inters) / unit, 2)) + '  ' + str(round(prevLenU / unit, 2)) + '\n'
                _logger.log(5, ind_str + 'added DRIFT: pos= {}, len={}, prevPosD={}, prevLenD={}, K_rms={}'.format(pos, L_inters, prevPosD, prevLenD, K_rms_add))
                
                #reset values
                L_inters = 0
                phi_inters = 0
                K_inters = None
                K_rms = 0
                
            L_und = float(e.nperiods) * float(e.lperiod) #remove?
            K_rms = np.sqrt(e.Kx**2 + e.Ky**2) / np.sqrt(2)
            
            # if not hasattr(e, 'K_err'):
                # e.K_err = None
            
            if np.size(K_rms) == 1:
                undLat += 'AW' + '    ' + str(K_rms) + '   ' + str(round(L_und / unit, 2)) + '  ' + str(round((pos - prevPosU - prevLenU) / unit, 2)) + '\n'
                _logger.log(5, ind_str + 'added UND:   pos= {}, len={}, prevPosU={}, prevLenU={}, K_rms={}'.format(pos,l, prevPosU, prevLenU, K_rms))
            else:
                unit_len_mult = e.nperiods / len(K_rms)
                _logger.log(5, ind_str + 'nperiods = {}, lperiod = {}, unit_len_mult = {}'.format(e.nperiods, e.lperiod, unit_len_mult))
                for period, K_rms_i in enumerate(K_rms):
                    if period == 0:
                        # undLat += 'AW' + '    ' + str(K_rms) + '   ' + str(1 * e.lperiod / unit * unit_len_mult) + '  ' + str(round((pos - prevPosU - prevLenU) / unit, 2)) + '\n' 
                        undLat += 'AW' + '    ' + str(K_rms_i) + '   ' + str(1 * e.lperiod / unit * unit_len_mult) + '  ' + str(round((pos - prevPosU - prevLenU) / unit, 2)) + '\n' # 
                        
                        # some bug in Genesis kick the beam if first period has non-resonant K_rms
                    elif period == len(K_rms):
                        undLat += 'AW' + '    ' + str(K_rms_i) + '   ' + str(1 * e.lperiod / unit * unit_len_mult) + '  ' + str(round(0)) + '\n' # some bug in Genesis kick the beam if first period has non-resonant K_rms
                    else:
                        undLat += 'AW' + '    ' + str(K_rms_i) + '   ' + str(1 * e.lperiod / unit * unit_len_mult) + '  ' + str(0) + '\n'
                
                
            
            prevPosU = pos
            prevLenU = L_und
            l_period = e.lperiod
            _logger.log(5, 'prevPosU = {}'.format(prevPosU))
            _logger.log(5, 'prevLenU = {}'.format(prevLenU))
            _logger.log(5, 'l_period = {}'.format(l_period))
        
        elif e.__class__ == UnknownElement and hasattr(e, 'phi'):
            phi_inters = e.phi
            if hasattr(e, 'K'):
                K_inters = e.K
                if K_inters == 'K_und':
                    K_inters = K_rms
            
        elif e.__class__ == Quadrupole:
            #k = energy/0.2998 * float(e.k1) *  ( e.l / unit - int(e.l / unit) )
            # k = float(energy) * float(e.k1) / e.l #*  (1 +  e.l / unit - int(e.l / unit) )
            # k = float(energy) * float(e.k1) * 0.2998 / e.l #*  (1 +  e.l / unit - int(e.l / unit) )
            k_quad = float(energy) * float(e.k1) / speed_of_light * 1e9
            _logger.log(5, ind_str + 'added QUAD:  pos= {}, len={}, prevPosQ={}, prevLenQ={}, K={}'.format(pos,l, prevPosQ, prevLenQ, k_quad))
            quadLat += 'QF' + '    ' + str(k_quad) + '   ' + str(round(e.l / unit, 2)) + '  ' + str(round((pos - prevPosQ - prevLenQ) / unit, 2)) + '\n'
            prevPosQ = pos
            prevLenQ = l
            # pass
        else:
            pass
            
        if e.__class__ != Undulator:
            L_inters = L_inters + e.l
        
        pos = pos + l
        
    full_lat = lat + undLat + driftLat + quadLat
    
    _logger.debug(ind_str + 'full lattice' + str(full_lat))
    return full_lat


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
    # spec = spec / np.sqrt(out.nSlices)/(2*out.leng/out('ncar'))**2/1e10
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
        # beam_new.column_values = beam.column_values
        # beam_new.columns = beam.columns
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










    z1 = M * u1[:, 0] * np.sqrt(l1[0])
    z2 = M * u1[:, 1] * np.sqrt(l1[1])

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
    p_tot = np.sqrt(adist[:, 3]**2 + adist[:, 4]**2 + adist[:, 5]**2)
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

def rematch_edist(edist, tws, s=None):

    betax_n = tws.beta_x
    betay_n = tws.beta_y
    alphax_n = tws.alpha_x
    alphay_n = tws.alpha_y

    edist_out = deepcopy(edist)
    # edist_out.center()

    x = edist_out.x
    y = edist_out.y
    xp = edist_out.xp
    yp = edist_out.yp

    mean_x2 = np.mean(x**2)
    mean_y2 = np.mean(y**2)
    mean_px2 = np.mean(xp**2)
    mean_py2 = np.mean(yp**2)
    mean_xpx = np.mean(x * xp)
    mean_ypy = np.mean(y * yp)
    mean_g = np.mean(edist_out.g)
    
    beam = edist2beam(edist_out)
    
    if s is None:
        tws0 = Twiss(beam.pk())
    else:
        tws0 = Twiss(beam.get_s(s))
        
    
    
    # emitx = mean_g * (mean_x2 * mean_px2 - mean_xpx**2)**0.5
    # emity = mean_g * (mean_y2 * mean_py2 - mean_ypy**2)**0.5
    # betax = mean_g * mean_x2 / emitx
    # betay = mean_g * mean_y2 / emity
    # alphax = -mean_g * mean_xpx / emitx
    # alphay = -mean_g * mean_ypy / emity
    
    betax = tws0.beta_x
    betay = tws0.beta_y
    alphax = tws0.alpha_x
    alphay = tws0.alpha_y

    # remove correlation
    xp = xp + x * alphax / betax
    yp = yp + y * alphay / betay

    # scale beam
    x = x * np.sqrt(betax_n / betax)
    y = y * np.sqrt(betay_n / betay)
    xp = xp * np.sqrt(betax / betax_n)
    yp = yp * np.sqrt(betay / betay_n)

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
