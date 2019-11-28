import time
import os
import socket
import copy

try:
    import h5py
    h5py_avail = True
except ImportError:
    print("wave.py: module h5py is not installed. Install it if you want to use genesis4 adaptor")
    h5py_avail = False

import numpy as np
from ocelot import ParticleArray
from ocelot.optics.wave import calc_ph_sp_dens, RadiationField
from ocelot.common.globals import *
from ocelot.adaptors.genesis import GenesisElectronDist #tmp
from ocelot.common.ocelog import *
from ocelot.utils.launcher import *
import os

_logger = logging.getLogger(__name__)

_inputGen4Template = "\
 $newrun \n\
 aw0   =  __AW0__ \n\
 xkx   =  __XKX__\n\
 \n"


class Namespace:
        def __init__(self):
            pass



class Genesis4Input:
    '''
    Genesis input files storage object
    '''
    def __init__(self):

        self.setup = Namespace()
        self.setup.rootname = ''
        self.setup.lattice = ''
        self.beamline = ''
        self.gamma0 = 11350.3

        self.alter_setup = Namespace()
        self.lattice = Namespace()

        self.exp_dir = None
        self.run_dir = None

        self.template = "\
        rootname=try_genesis_s2_field     \n\
        lattice=lattice.lat     \n\
        beamline=SASE3     \n\
        lambda0=1.77120e-9     \n\
        gamma0=16634.08     \n\
        delz=0.017     \n\
        shotnoise=1     \n\
        &end     \n\
             \n\
        &time     \n\
        slen=2e-6     \n\
        sample=12     \n\
        &end     \n\
             \n\
        #&field     \n\
        #power=0     \n\
        #dgrid=3e-4     \n\
        #ngrid=111     \n\
        #waist_size=30e-6     \n\
        #&end     \n\
             \n\
        &profile_gauss     \n\
        label=prof2     \n\
        c0=3000     \n\
        s0=0.4e-6     \n\
        sig=1e-6     \n\
        &end     \n\
             \n\
        &profile_polynom     \n\
        label=prof_lin     \n\
        c0=0.3e-6     \n\
        c1=0.01     \n\
        &end     \n\
             \n\
        &beam     \n\
        current=2000     \n\
        delgam=3.5     \n\
        ex=0.5e-6     \n\
        ey=0.5e-6     \n\
        alphax = 1.33377002518     \n\
        alphay = -0.787227424221     \n\
        betax = 25     \n\
        betay = 15     \n\
        &end     \n\
             \n\
        &importfield     \n\
        file=try_genesis_s1.fld.h5     \n\
        &end     \n\
             \n\
        #&importbeam     \n\
        #file=try_genesis_s1.par.h5     \n\
        #charge = 34e-12     \n\
        #&end     \n\
             \n\
        #&lattice     \n\
        #zmatch=5     \n\
        #&end     \n\
             \n\
        &track     \n\
        zstop=50     \n\
        &end     \n\
             \n\
        &write     \n\
        field=try_genesis_s2_field     \n\
        beam=try_genesis_s2_field     \n\
        &end     \n\
        \n"

    def input(self):
        return self.template

class Genesis4ParticlesDump:
    '''
    Genesis particle *.dpa files storage object
    Each particle record in z starts with the energy of all particles 
    followed by the output of the particle phases, 
    positions in x and y and the momenta in x and y. 
    The momenta are normalized to mc
    '''

    def __init__(self):
        self.g = []
        self.ph = []
        self.x = []
        self.y = []
        self.px = []
        self.py = []
        self.I = []

        self.one4one = None

        # self.fileName = ''
        self.filePath = ''

    def fileName(self):
        return os.path.basename(self.filePath)
        # return filename_from_path(self.filePath)


class Genesis4Output:
    '''
    Genesis input files storage object
    '''
    def __init__(self):
        self.h5 = None #hdf5 pointer

    def close(self):
        self.h5.close()

    @property
    def filePath(self):
        return self.h5.filename

    def fileName(self):
        return os.path.basename(self.h5.filename)

    @property
    def nZ(self):
        return self.z.size

    @property
    def nSlices(self):
        return self.h5['Beam/current'].size

    @property
    def lambdaref(self):
        return self.h5['Global/lambdaref'][0]

    @property
    def phenref(self):
        return h_eV_s * speed_of_light / self.lambdaref

    @property
    def I(self):
        return self.h5['Beam/current'][0]

    @property
    def beam_charge(self):
        return np.trapz(self.I, self.t)

    @property
    def rad_power(self):
        return self.h5['Field/power']

    @property
    def rad_energy(self):
        return np.trapz(self.rad_power, self.t)

    @property
    def n_photons(self):
        return self.rad_energy / q_e / self.phenref

    @property
    def t(self):
        if not self.tdp:
            return None
        else:
            return self.s / speed_of_light

    def rad_field(self, zi=None, loc='near'):
        if loc == 'far':
            intens = self.h5['Field/intensity-farfield']
            phase = self.h5['Field/phase-farfield']
        elif loc == 'near':
            intens = self.h5['Field/intensity-nearfield']
            phase = self.h5['Field/phase-nearfield']
        else:
            raise ValueError('loc should be either "far" or "near"')

        if zi is not None:
            intens = intens[zi,:]
            phase = phase[zi,:]

        #not scaled properly!!!!
        field = np.sqrt(intens[:]) * np.exp(1j * phase[:])
        return field

    def calc_spec(self, zi=None, loc='near', npad=1, estimate_ph_sp_dens=1):

        field = self.rad_field(zi=zi, loc=loc)
        axis = field.ndim - 1

        spec = np.abs(np.fft.fft(field, axis=axis))**2
        spec = np.fft.fftshift(spec, axes=axis)

        scale_ev = h_eV_s * speed_of_light * (np.fft.fftfreq(self.nSlices, d=self.s[1]-self.s[0]) + 1 / self.lambdaref)
        scale_ev = np.fft.fftshift(scale_ev)

        if estimate_ph_sp_dens:
            tt=np.trapz(spec, scale_ev, axis=axis)
            if axis==1:
                tt[tt==0] = np.inf
                spec *= (self.n_photons / tt)[:, np.newaxis]
            else:
                if tt==0:
                    tt = np.inf
                spec *= (self.n_photons[zi] / tt)

        return scale_ev, spec

    def wig(self,z=np.inf):
        return wigner_out(self, z=z, method='mp', debug=1)

    def close(self):
        _logger.warning('closing h5 file {}'.format(self.filePath))
        self.h5.close()


def get_genesis4_launcher(launcher_program='genesis4', launcher_argument=''):
    '''
    Returns MpiLauncher() object for given program
    '''
    host = socket.gethostname()

    launcher = MpiLauncher()
    launcher.program = launcher_program
    launcher.argument = launcher_argument
    launcher.mpiParameters = '-x PATH -x MPI_PYTHON_SITEARCH -x PYTHONPATH'  # added -n
    # launcher.program = '/data/netapp/xfel/products/genesis/genesis'
    # launcher.argument = ' < tmp.cmd | tee log'

    return launcher


def run_genesis4(inp, launcher, *args, **kwargs):
    '''
    Main function for executing Genesis code
    inp               - GenesisInput() object with genesis input parameters
    launcher          - MpiLauncher() object obtained via get_genesis_launcher() function
    '''
    # import traceback
    _logger.info('Starting genesis v4 preparation')
    # _logger.warning(len(traceback.extract_stack()))

    # create experimental directory
    if inp.run_dir is None and inp.exp_dir is None:
        raise ValueError('run_dir and exp_dir are not specified!')

    if inp.run_dir is None:
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

    inp_file = inp.inp_file
    inp_path = os.path.join(inp.exp_dir, inp.inp_file)
    # if inp.stageid is None:
        # inp_path = inp.run_dir + 'run.' + str(inp.runid) + inp.suffix + '.inp'
        # out_path = inp.run_dir + 'run.' + str(inp.runid) + inp.suffix + '.gout'
        # inp.stageid = ''
        # stage_string = ''
    # else:
        # inp_path = inp.run_dir + 'run.' + str(inp.runid) + '.s' + str(inp.stageid) + inp.suffix + '.inp'
        # out_path = inp.run_dir + 'run.' + str(inp.runid) + '.s' + str(inp.stageid) + inp.suffix + '.gout'
        # stage_string = '.s' + str(inp.stageid)

    # inp_file = filename_from_path(inp_path)
    # out_file = filename_from_path(out_path)

    # # cleaning directory
    # _logger.debug(ind_str + 'removing old files')
    # os.system('rm -rf ' + inp.run_dir + 'run.' + str(inp.runid) + stage_string + inp.suffix + '*')  # to make sure all stage files are cleaned
    # # os.system('rm -rf ' + out_path+'*') # to make sure out files are cleaned
    # # os.system('rm -rf ' + inp_path+'*') # to make sure inp files are cleaned
    # os.system('rm -rf ' + inp.run_dir + 'tmp.cmd')

    # # create and fill necessary input files
    # if inp.latticefile == None:
        # if inp.lat != None:
            # _logger.debug(ind_str + 'writing ' + inp_file + '.lat')
            # if not hasattr(inp, 'lat_unit'):
                # lat_unit = inp.xlamd
                # _logger.debug(2*ind_str + 'lat_unit_size = xlamds = {} m'.format(lat_unit))
            # else:
                # if inp.lat_unit is None:
                    # lat_unit = inp.xlamd
                    # _logger.debug(2*ind_str + 'lat_unit_size = xlamds = {} m'.format(lat_unit))
                # else:
                    # lat_unit = inp.lat_unit
                    # _logger.debug(2*ind_str + 'lat_unit_size = {} m'.format(lat_unit))
            # open(inp_path + '.lat', 'w').write(generate_lattice(inp.lat, unit=lat_unit, energy=inp.gamma0 * m_e_GeV, debug = debug, min_phsh = min_phsh))
            # inp.latticefile = inp_file + '.lat'

    # if inp.beamfile == None:
        # if inp.beam != None:
            # _logger.debug(ind_str + 'writing ' + inp_file + '.beam')
            # open(inp_path + '.beam', 'w').write(beam_file_str(inp.beam))
            # inp.beamfile = inp_file + '.beam'

    # if inp.edistfile == None:
        # if inp.edist != None:
            # _logger.debug(ind_str + 'writing ' + inp_file + '.edist')
            # write_edist_file(inp.edist, inp_path + '.edist', debug=1)
            # inp.edistfile = inp_file + '.edist'

    # if inp.partfile == None:
        # if inp.dpa != None:
            # _logger.debug(ind_str + 'writing ' + inp_file + '.dpa')
            # # print ('!!!!!!! no write_particle_file() function')
            # write_dpa_file(inp.dpa, inp_path + '.dpa', debug=1)
            # inp.partfile = inp_file + '.dpa'

    # if inp.fieldfile == None:
        # if inp.dfl != None:
            # _logger.debug(ind_str + 'writing ' + inp_file + '.dfl')
            # write_dfl_file(inp.dfl, inp_path + '.dfl', debug=1)
            # inp.fieldfile = inp_file + '.dfl'

    # if inp.radfile == None:
        # if inp.rad != None:
            # _logger.debug(ind_str + 'writing ' + inp_file + '.rad')
            # open(inp_path + '.rad', 'w').write(rad_file_str(inp.rad))
            # inp.radfile = inp_file + '.rad'

    # if inp.outputfile == None:
        # inp.outputfile = out_file
    _logger.debug(ind_str + 'writing ' + inp_file)
    open(inp_path, 'w').write(inp.input())
    # open(inp.run_dir + 'tmp.cmd', 'w').write(inp_file + '\n')

    launcher.dir = inp.run_dir
    _logger.debug(ind_str + 'preparing launcher')
    launcher.prepare()
    # _logger.debug()
    # RUNNING GENESIS ###
    # genesis_time = time.time()
    # launcher.launch()
    # _logger.info(ind_str + 'genesis simulation time %.2f seconds' % (time.time() - genesis_time))
    # RUNNING GENESIS ###

    # if assembly_ver is not None:
        # # genesis output slices assembly
        # _logger.info(ind_str + 'assembling slices')
        # _logger.debug(2 * ind_str + 'assembly_ver = {}'.format(assembly_ver))

        # assembly_time = time.time()


        # if assembly_ver == 'sys':

            # _logger.info(2 * ind_str + 'assembling *.out file')
            # start_time = time.time()
            # os.system('cat ' + out_path + '.slice* >> ' + out_path)
            # os.system('rm ' + out_path + '.slice* 2>/dev/null')
            # _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

            # _logger.info(2 * ind_str + 'assembling *.dfl file')
            # start_time = time.time()
            # if dfl_slipage_incl:
                # os.system('cat ' + out_path + '.dfl.slice*  >> ' + out_path + '.dfl.tmp')
                # #bytes=os.path.getsize(out_path +'.dfl.tmp')
                # command = 'dd if=' + out_path + '.dfl.tmp of=' + out_path + '.dfl conv=notrunc conv=notrunc 2>/dev/null'# obs='+str(bytes)+' skip=1
                # os.system(command)
            # else:
                # os.system('cat ' + out_path + '.dfl.slice*  > ' + out_path + '.dfl')
            # os.system('rm ' + out_path + '.dfl.slice* 2>/dev/null')
            # os.system('rm ' + out_path + '.dfl.tmp 2>/dev/null')
            # _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

            # _logger.info(2 * ind_str + 'assembling *.dpa file')
            # start_time = time.time()
            # os.system('cat ' + out_path + '.dpa.slice* >> ' + out_path + '.dpa')
            # os.system('rm ' + out_path + '.dpa.slice* 2>/dev/null')
            # _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

            # _logger.info(2 * ind_str + 'assembling *.fld file')
            # start_time = time.time()
            # if dfl_slipage_incl:
                # os.system('cat ' + out_path + '.fld.slice*  >> ' + out_path + '.fld.tmp')
                # #bytes=os.path.getsize(out_path +'.fld.tmp')
                # command = 'dd if=' + out_path + '.fld.tmp of=' + out_path + '.fld conv=notrunc conv=notrunc 2>/dev/null'# obs='+str(bytes)+' skip=1
                # os.system(command)
            # else:
                # os.system('cat ' + out_path + '.fld.slice*  > ' + out_path + '.fld')
            # os.system('rm ' + out_path + '.fld.slice* 2>/dev/null')
            # os.system('rm ' + out_path + '.fld.tmp 2>/dev/null')
            # _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

            # _logger.info(2 * ind_str + 'assembling *.par file')
            # start_time = time.time()
            # os.system('cat ' + out_path + '.par.slice* >> ' + out_path + '.par')
            # os.system('rm ' + out_path + '.par.slice* 2>/dev/null')
            # _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

        # elif assembly_ver == 'pyt':
            # # there is a bug with dfl assembly
            # import glob
            # ram = 1

            # _logger.info(2 * ind_str + 'assembling *.out file')
            # start_time = time.time()
            # assemble(out_path, ram=ram, debug=debug)
            # os.system('rm ' + out_path + '.slice* 2>/dev/null')
            # _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

            # for i in range(10):        #read all possible harmonics (up to 10 now)
                # ii=str(i)
                # if ii=='0': ii=''
                # if os.path.isfile(str(out_path + '.dfl' + ii)):
                    # _logger.info(2 * ind_str + 'assembling *.dfl'+ii+' file')
                    # start_time = time.time()
                    # assemble(out_path + '.dfl'+ii, overwrite=dfl_slipage_incl, ram=ram, debug=debug)
                    # os.system('rm ' + out_path + '.dfl'+ii+'.slice* 2>/dev/null')
                    # os.system('rm ' + out_path + '.dfl'+ii+'.tmp 2>/dev/null')
                    # _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

            # if os.path.isfile(str(out_path + '.dpa')):
                # _logger.info(2 * ind_str + 'assembling *.dpa file')
                # start_time = time.time()
                # assemble(out_path + '.dpa', ram=ram, debug=debug)
                # os.system('rm ' + out_path + '.dpa.slice* 2>/dev/null')
                # _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

            # if os.path.isfile(str(out_path + '.fld')):
                # _logger.info(2 * ind_str + 'assembling *.fld file')
                # start_time = time.time()
                # assemble(out_path + '.fld', overwrite=dfl_slipage_incl, ram=ram, debug=debug)
                # os.system('rm ' + out_path + '.fld.slice* 2>/dev/null')
                # os.system('rm ' + out_path + '.fld.tmp 2>/dev/null')
                # _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

            # if os.path.isfile(str(out_path + '.par')):
                # _logger.info(2 * ind_str + 'assembling *.par file')
                # start_time = time.time()
                # assemble(out_path + '.par', ram=ram, debug=debug)
                # os.system('rm ' + out_path + '.par.slice* 2>/dev/null')
                # _logger.debug(3 * ind_str + 'done in %.2f seconds' % (time.time() - start_time))

        # else:
            # # raise ValueError('assembly_ver should be either "sys" or "pyt"')
            # pass
        # _logger.debug(2 * ind_str + 'assembly time %.2f seconds' % (time.time() - assembly_time))



    # if read_level >= 0:
        # out = read_out_file(out_path, read_level=read_level)
        # _logger.debug(ind_str + 'done, time %.2f seconds' % (time.time() - assembly_time))
        # return out
    # else:
        # _logger.debug(ind_str + 'done, time %.2f seconds' % (time.time() - assembly_time))
        # return None
    return


def read_gout4(filePath):
    '''
    Reads Genesis1.3 v4 output file with the link to the hdf5 file as out.h5
    to close the file, use out.h5.close()
    
    :param filePath: string, absolute path to .out file
    :returns: Genesis4Output
    '''

    _logger.info('reading gen4 .out file')
    _logger.warning(ind_str + 'in beta')
    _logger.debug(ind_str + 'reading from ' + filePath)

    out = Genesis4Output()
    try:
        out.h5 = h5py.File(filePath, 'r')
    except Exception:
        _logger.error(ind_str + 'no such file ' + filePath)
        raise

    vvv = [int(out.h5['Meta/Version/'+name][0]) for name in ['Major', 'Minor', 'Revision']]
    _logger.debug(ind_str + 'Genesis v{}.{}.{}'.format(*vvv))

    out.z = out.h5['Lattice/zplot'][:]
    out.zlat = out.h5['Lattice/z'][:]

    _logger.debug(ind_str + 'z[0]   , z[-1]   , len(z)    = {}  {}  {}'.format(out.z[0], out.z[-1], len(out.z)))
    _logger.debug(ind_str + 'zlat[0], zlat[-1], len(zlat) = {}  {}  {}'.format(out.zlat[0], out.zlat[-1], len(out.zlat)))


    if 'time' in out.h5['Global'] and out.h5['Global/time'][0] == 1:
        out.tdp = True
        _logger.debug(ind_str + 'tdp = True')
    else:
        out.tdp = False
        _logger.debug(ind_str + 'tdp = False')

    if out.tdp:
        if 's0' in out.h5['Global']:
            s0 = out.h5['Global/s0']
            _logger.debug(ind_str + 's0 = {}'.format(s0))
        else:
            s0 = 0
        sn = out.h5['Beam/current'].size
        out.s = np.linspace(s0, out.h5['Global/slen'][:][0]+s0, sn)
        _logger.debug(ind_str + 's[0], s[-1], len(s) = {}  {}  {}'.format(out.s[0], out.s[-1], sn))

    _logger.debug(ind_str + 'done')

    return out

def read_dfl4(filePath):
    '''
    Reads Genesis1.3 v4 radiation output file
    
    :param filePath: string, absolute path to .fld file
    :returns: RadiationField
    '''
    _logger.info('reading gen4 .dfl file')
    _logger.warning(ind_str + 'in beta')
    _logger.debug(ind_str + 'reading from ' + filePath)

    with h5py.File(filePath, 'r') as h5:

        nslice = h5.get('slicecount')[0]
        lambdaref = h5.get('wavelength')[0]
        sepslice = h5.get('slicespacing')[0]
        gridsize = h5.get('gridsize')[0]
        try:
            ncar = h5.get('gridpoints')[0]
        except:
            # legacy support
            ncar = int(np.sqrt(h5.get('slice000001/field-real').size)) # fix?
        zsep = int(sepslice / lambdaref)
        l_total = sepslice * nslice

        _logger.debug(ind_str + 'nslice = {}'.format(nslice))
        _logger.debug(ind_str + 'sepslice (dz) = {}'.format(sepslice))
        _logger.debug(ind_str + 'lambdaref (xlamds) = {}'.format(lambdaref))
        _logger.debug(ind_str + 'gridsize (dx & dy) = {}'.format(gridsize))
        _logger.debug(ind_str + '')
        _logger.debug(ind_str + 'zsep = {}'.format(zsep))
        _logger.debug(ind_str + 'Nx & Ny = {}'.format(ncar))
        _logger.debug(ind_str + 'Lx & Ly = {}'.format(ncar*gridsize))
        _logger.debug(ind_str + 'Ls_total = {}'.format(l_total))

        field_real = []
        field_imag = []

        for dset in h5:
            if dset.startswith('slice0'):
                _logger.log(5, '{}'.format(dset))
                field_real.append(h5[dset]['field-real'][:].reshape(ncar, ncar))
                field_imag.append(h5[dset]['field-imag'][:].reshape(ncar, ncar))

        dfl = RadiationField()
        dfl.fld = np.array(field_real) + 1j * np.array(field_imag)
        dfl.dx = gridsize
        dfl.dy = gridsize
        dfl.dz = sepslice
        dfl.xlamds = lambdaref
        dfl.domain_z = 't'  # longitudinal domain (t - time, f - frequency)
        dfl.domain_xy = 's'  # transverse domain (s - space, k - inverse space
        dfl.filePath = h5.filename
        # legacy support
        try:
            dfl.refposition = h5.get('refposition')[0]
        except:
            pass

    _logger.debug(ind_str + 'done')

    return dfl


def write_dfl4(dfl: RadiationField, file_path='sample.dfl.h5'):
    """
    Writes ocelot.optics.wave.RadiationField object to Genesis1.3 v4 radiation file

    :param dfl: ocelot.optics.wave.RadiationField object
    :param file_path: path to .dfl file (file will be generate, or data will be rewritten)
    :return:
    """
    _logger.info('writing gen4 file {}'.format(file_path))
    _logger.warning(ind_str + 'in beta')

    _logger.debug(ind_str + 'writing to ' + file_path)
    _logger.debug(ind_str + 'nslice = {}'.format(dfl.Nz()))
    _logger.debug(ind_str + 'sepslice (dz) = {}'.format(dfl.xlamds))
    _logger.debug(ind_str + 'lambdaref (xlamds) = {}'.format(dfl.dz))
    if dfl.dx != dfl.dy:
        _logger.error('dfl.dx is not equal dfl.dy')
        raise ValueError('dfl.dx is not equal dfl.dy')
    else:
        _logger.debug(ind_str + 'gridsize (dx & dy) = {}'.format(dfl.dx))
    _logger.debug(ind_str + '')
    _logger.debug(ind_str + 'zsep = {}'.format(int(dfl.dz / dfl.xlamds)))
    if dfl.Nx() != dfl.Ny():
        _logger.error('dfl.Nx() is not equal dfl.Ny()')
        raise ValueError('dfl.Nx() is not equal dfl.Ny()')
    else:
        _logger.debug(ind_str + 'Nx & Ny = {}'.format(dfl.Nx()))
    _logger.debug(ind_str + 'Lx & Ly = {}'.format(dfl.Lx()))
    _logger.debug(ind_str + 'Ls_total = {}'.format(dfl.Lz()))

    with h5py.File(file_path, 'w') as h5:

        h5.create_dataset('slicecount', data=[dfl.Nz()])
        h5.create_dataset('wavelength', data=[dfl.xlamds])
        h5.create_dataset('slicespacing', data=[dfl.dz])
        h5.create_dataset('gridsize', data=[dfl.dx])
        h5.create_dataset('gridpoints', data=[dfl.Nx()])
        try:
            h5.create_dataset('refposition', data=[dfl.refposition])
        except:
            h5.create_dataset('refposition', data=[0.0])


        for i in range(dfl.Nz()):
            _logger.log(5, 'slice{:06d}'.format(i+1))
            h5.create_dataset('slice{:06d}/field-real'.format(i+1), data=np.real(dfl.fld[i]).flatten())
            h5.create_dataset('slice{:06d}/field-imag'.format(i+1), data=np.imag(dfl.fld[i]).flatten())

    _logger.debug(ind_str + 'done')


def read_dpa4(filePath, start_slice=0, stop_slice=np.inf, estimate_npart=0):
    '''
    Reads Genesis1.3 v4 particle output file
    
    :param filePath: string, absolute path to .par file
    :returns: Genesis4ParticlesDump
    '''

    _logger.info('reading gen4 .dpa file')
    _logger.warning(ind_str + 'in beta')
    _logger.debug(ind_str + 'reading from ' + filePath)

    start_time = time.time()

    _logger.debug(ind_str + 'start_slice : stop_slice = {} : {}'.format(start_slice, stop_slice))

    with h5py.File(filePath, 'r') as h5:

        one4one = bool(h5.get('one4one')[0])
        _logger.info(ind_str + 'one4one = {}'.format(bool(one4one)))

        nslice = int(h5.get('slicecount')[0])
        nbins = int(h5.get('beamletsize')[0])
        lslice = h5.get('slicelength')[0]
        sepslice = h5.get('slicespacing')[0]

        if not one4one:
            npart = int(h5.get('slice000001/gamma').size) # fix?
            _logger.debug(ind_str + 'npart = {}'.format(npart))
        else:
            if estimate_npart:
                I_tmp_full = []
                I_tmp_wind = []
                for dset in h5:
                    if dset.startswith('slice') and type(h5[dset]) == h5py._hl.group.Group:
                        slice = int(dset.replace('slice',''))
                        I_tmp_full.append(h5[dset]['current'][:])
                        if slice >= start_slice and slice <= stop_slice:
                            I_tmp_wind.append(h5[dset]['current'][:])
                npart_full_tmp = np.sum(I_tmp_full) * lslice / speed_of_light / q_e
                npart_wind_tmp = np.sum(I_tmp_wind) * lslice / speed_of_light / q_e
                _logger.info(ind_str + 'estimated npart = {:}M'.format(npart_full_tmp/1e6))
                _logger.info(ind_str + 'estimated npart to be downloaded = {:}.M'.format(npart_wind_tmp/1e6))


        # filePath = h5.filename

        zsep = int(sepslice / lslice)
        l_total = lslice * zsep * nslice


        _logger.debug(ind_str + 'nslice = {}'.format(nslice))
        _logger.debug(ind_str + 'nbins = {}'.format(nbins))
        _logger.debug(ind_str + '')
        _logger.debug(ind_str + 'lslice (aka xlamds) = {} m'.format(lslice))
        _logger.debug(ind_str + 'sepslice = {} m'.format(sepslice))
        _logger.debug(ind_str + 'zsep = {}'.format(zsep))
        _logger.debug(ind_str + 'Ls_total = {}'.format(l_total))

        x = []
        y = []
        px = []
        py = []
        g = []
        ph = []
        # s0 = []
        # s = []
        I = []
        npartpb = []

        _logger.debug(ind_str + 'reading slices between {} and {}'.format(start_slice, stop_slice))
        _logger.debug(2*ind_str + '({} out of {})'.format(stop_slice - start_slice, nslice))
        for dset in sorted(h5):
            # _logger.log(5, ind_str + dset)
            # _logger.log(5, ind_str + str(type(h5[dset])))
            if dset.startswith('slice') and type(h5[dset]) == h5py._hl.group.Group:
                slice = int(dset.replace('slice',''))
                _logger.log(5, 2*ind_str + 'slice number {}'.format(slice))
                if slice >= start_slice and slice <= stop_slice:
                    _logger.log(5, 2*ind_str + 'processing')
                    I.append(h5[dset]['current'][:])
                    ph.append(h5[dset]['theta'][:])
                    _logger.log(5, 2*ind_str + '{} particles'.format(h5[dset]['x'].size))
                    # s.append(s0 + ph0 / 2 / np.pi * lslice)
                    x.append(h5[dset]['x'][:])
                    px.append(h5[dset]['px'][:])
                    y.append(h5[dset]['y'][:])
                    py.append(h5[dset]['py'][:])
                    g.append(h5[dset]['gamma'][:])
                    npartpb.append(h5[dset]['gamma'][:].size)
                    # ph.append(ph0)
                    # s0 += sepslice
        _logger.debug(2*ind_str + 'done')

    dpa = Genesis4ParticlesDump()

    if one4one:
        _logger.debug(ind_str + 'flattening arrays')
        # _logger.log(5, 2*ind_str + 'x.shape {}'.format(np.shape(x)))
        dpa.x = np.hstack(x)
        dpa.px = np.hstack(px)
        dpa.y = np.hstack(y)
        dpa.py = np.hstack(py)
        dpa.g = np.hstack(g)
        dpa.I = np.hstack(I)
        dpa.ph = np.hstack(ph)
        _logger.log(5, 2*ind_str + 'dpa.x.shape {}'.format(dpa.x.shape))

        dpa.npartpb = np.array(npartpb).flatten()

        npart = dpa.x.size
        _logger.info(ind_str + 'npart = {}'.format(npart))
    else:
        npartpb = int(npart/nbins)
        # _logger.debug(ind_str + 'not one4one:')
        _logger.debug(2*ind_str + 'shape of arrays (nslice, npart) {}'.format(np.array(x).shape))

        if np.array(x).shape[0] < nslice:
            nslice = np.array(x).shape[0]

        if nslice == 0:
            _logger.error(2*ind_str + 'nslice == 0')
            return

        _logger.debug(2*ind_str + 'reshaping to (nslice, nbins, npart/bin) ({}, {}, {})'.format(nslice, nbins, npartpb))

        dpa.x = np.array(x).reshape((nslice, nbins, npartpb), order='F')
        dpa.px = np.array(px).reshape((nslice, nbins, npartpb), order='F')
        dpa.y = np.array(y).reshape((nslice, nbins, npartpb), order='F')
        dpa.py = np.array(py).reshape((nslice, nbins, npartpb), order='F')
        dpa.ph = np.array(ph).reshape((nslice, nbins, npartpb), order='F')
        dpa.g = np.array(g).reshape((nslice, nbins, npartpb), order='F')
        dpa.I = np.array(I).flatten()

    _logger.debug(ind_str + 'writing to dpa object')
    dpa.nslice = nslice
    dpa.lslice = lslice
    dpa.npart = npart
    dpa.nbins = nbins
    dpa.zsep = zsep
    dpa.filePath = filePath
    dpa.one4one = one4one

    _logger.debug(ind_str + 'done in %.2f seconds' % (time.time() - start_time))

    return dpa





def dpa42edist(dpa, n_part=None, fill_gaps=False):
    '''
    Convert Genesis1.3 v4 particle output file to ocelot edist object
    
    :param dpa: GenesisParticlesDump
    :param n_part: desired approximate number of particles in edist
    :param fill_gaps: dublicates buckets into gaps
    :returns: GenesisElectronDist
    '''

    _logger.info('converting dpa4 to edist')
    _logger.warning(ind_str + 'in beta')

    # if dpa.one4one is True:
        # _logger.warning(ind_str + 'not tested for one4one yet')
        # return
    if dpa.one4one is None:
        _logger.error(ind_str + 'unknown one4one status')
        return

    # fill_gaps=False #unfinished


    import random
    start_time = time.time()

    npart = dpa.npart
    nslice = dpa.nslice
    nbins = dpa.nbins
    lslice = dpa.lslice
    zsep = dpa.zsep

    _logger.debug(ind_str + 'requested n_part = {}'.format(n_part))
    _logger.debug(ind_str + 'dpa npart = {}'.format(npart))
    _logger.debug(ind_str + 'dpa nslice = {}'.format(nslice))
    _logger.debug(ind_str + 'dpa nbins = {}'.format(nbins))
    _logger.debug(ind_str + 'dpa lslice (aka xlamds) = {} m'.format(lslice))
    _logger.debug(ind_str + 'dpa zsep = {}'.format(zsep))
    _logger.debug(ind_str + '')

    _logger.debug(ind_str + 'fill_gaps = {}'.format(fill_gaps))

    if dpa.one4one:
        t0 = np.hstack([np.ones(n)*i for i, n in enumerate(dpa.npartpb)]) * dpa.lslice / speed_of_light
        t0 += dpa.ph / 2 / np.pi * dpa.lslice / speed_of_light

        edist = GenesisElectronDist()

        C = npart * q_e

        if n_part is not None:
            pick_i = random.sample(range(npart), n_part)
            # particles_kept_ratio = n_part / npart
            # C = n_part * q_e
        else:
            pick_i = range(npart)
            n_part = npart
            # particles_kept_ratio = 1


        _logger.debug(ind_str + 'particles kept = {}%'.format(n_part / npart * 100))
        _logger.debug(ind_str + 'picking {} of {} particles'.format(len(pick_i), npart))

        if n_part == npart:
            edist.g = dpa.g
            edist.xp = dpa.px / edist.g
            edist.yp = dpa.py / edist.g
            edist.x = dpa.x
            edist.y = dpa.y
            edist.t = t0
        else:
            edist.g = dpa.g[pick_i]
            edist.xp = dpa.px[pick_i] / edist.g
            edist.yp = dpa.py[pick_i] / edist.g
            edist.x = dpa.x[pick_i]
            edist.y = dpa.y[pick_i]
            edist.t = t0[pick_i]

        _logger.debug(2*ind_str + 'done')

    else:
        # if fill_gaps:

            # s0 = np.linspace(0, nslice * zsep * lslice, nslice)
            # s = np.linspace(0, nslice * zsep * lslice, nslice * zsep)
            # I = np.interp(s, s0, dpa.I)
            # dt = (s[1] - s[0]) / speed_of_light
        # else:
        s = np.linspace(0, nslice * zsep * lslice, nslice)
        I = dpa.I
        dt = zsep * lslice / speed_of_light
        # dt = (s[1] - s[0]) / speed_of_light
        _logger.debug(ind_str + 'dt zsep = {} s'.format(dt))

        C = np.sum(I) * dt

        _logger.debug(ind_str + 'current max = {} A'.format(I.max()))
        _logger.debug(ind_str + 'bunch charge = {} C'.format(C))

        n_part_max = int(np.floor(np.sum(I / I.max() * dpa.npart)))
        _logger.debug(ind_str + 'maximum possible n_part = {}'.format(n_part_max))
        if n_part is not None:
            _logger.debug(ind_str + 'requested n_part = {}'.format(n_part))
            if n_part > n_part_max:
                _logger.info(ind_str + 'requested n_part > maximum, setting to maximum = {}'.format(n_part_max))
                n_part = n_part_max
        else:
            _logger.info(ind_str + 'requested n_part = None, setting to maximum = {}'.format(n_part_max))
            n_part = n_part_max



        #n_part_bin = (I / np.sum(I) * n_part / nbins).astype(int)
        #n_part_bin[n_part_bin > npart_bin] = npart_bin
        #print(n_part_bin.max())
        n_part_slice = (I / np.sum(I) * n_part).astype(int)
        n_part_slice[n_part_slice > npart] = npart

        _logger.debug(ind_str + 'max particles/slice = {}'.format(n_part_slice.max()))


        # print(n_part_slice.max())


        #pick_i = random.sample(range(n_part), n_part_slice[i])
        # _logger.debug(ind_str + 'reshaping')

        # g = np.reshape(dpa.g, (nslice, npart))
        # x = np.reshape(dpa.x, (nslice, npart))
        # y = np.reshape(dpa.y, (nslice, npart))
        # px = np.reshape(dpa.px, (nslice, npart)) / g
        # py = np.reshape(dpa.py, (nslice, npart)) / g
        # ph = np.reshape(dpa.ph, (nslice, npart))


        g = dpa.g.reshape(nslice, npart)
        x = dpa.x.reshape(nslice, npart)
        y = dpa.y.reshape(nslice, npart)
        px = dpa.px.reshape(nslice, npart) / g
        py = dpa.py.reshape(nslice, npart) / g
        ph = dpa.ph.reshape(nslice, npart)

        t1 = ph / 2 / np.pi * lslice / speed_of_light
        t0 = np.arange(nslice)[:,np.newaxis] * lslice * zsep / speed_of_light

        t = t1+t0
        # _logger.debug(2*ind_str + 'complete')


        edist = GenesisElectronDist()
        #g1 = np.array([])
        #x1 = np.array([])
        #y1 = np.array([])

        gi = []
        xpi = []
        ypi = []
        xi = []
        yi = []
        ti = []

        _logger.debug(ind_str + 'picking random particles')
        for i in np.arange(nslice):
        #    for ii in np.arange(nbins):
        #    pick_i = random.sample(range(npart_bin), n_part_bin[i])
            pick_i = random.sample(range(npart), n_part_slice[i])


            # edist.g = np.append(edist.g, g[i, pick_i])
            # edist.xp = np.append(edist.xp, px[i, pick_i])
            # edist.yp = np.append(edist.yp, py[i, pick_i])
            # edist.x = np.append(edist.x, x[i, pick_i])  
            # edist.y = np.append(edist.y, y[i, pick_i])
            # edist.t = np.append(edist.t, t[i, pick_i])

            gi.append(g[i, pick_i])
            xpi.append(px[i, pick_i])
            ypi.append(py[i, pick_i])
            xi.append(x[i, pick_i])
            yi.append(y[i, pick_i])
            ti.append(t[i, pick_i])

        edist.g = np.array(gi).flatten()
        edist.xp = np.array(xpi).flatten()
        edist.yp = np.array(ypi).flatten()
        edist.x = np.array(xi).flatten()
        edist.y = np.array(yi).flatten()
        edist.t = np.array(ti).flatten()

        if fill_gaps:
            edist.t += np.random.randint(0, dpa.zsep, edist.t.size) * lslice / speed_of_light


        # particles_kept_ratio = 1

    edist.part_charge = C / edist.len()
    _logger.debug('')

    _logger.debug(ind_str + 'edist bunch charge = {} C'.format(C))
    _logger.debug(ind_str + 'edist particle charge = {} C (~{:.2f}*q_e)'.format(edist.part_charge,edist.part_charge/q_e))
    _logger.debug(ind_str + 'edist n_part = {}'.format(edist.len()))

    if hasattr(dpa,'filePath'):
        edist.filePath = dpa.filePath + '.edist'

    _logger.debug(ind_str + 'done')

    return edist

def read_dpa42parray(filePath, N_part=None, fill_gaps=True):

    _logger.info('reading gen4 .dpa file into parray')
    _logger.warning(ind_str + 'in beta')
    _logger.debug(ind_str + 'reading from ' + filePath)


    import random
    N_part = None

    #N_part = 100000
    fill_gaps=True
    _logger.debug('fill_gaps = ' + str(fill_gaps))

    h5 = h5py.File(filePath, 'r')

    nslice = int(h5.get('slicecount')[0])
    lslice = h5.get('slicelength')[0]
    sepslice = h5.get('slicespacing')[0]
    npart = int(h5.get('slice000001/gamma').size)
    nbins = int(h5.get('beamletsize')[0])
    zsep = int(sepslice / lslice)

    one4one = bool(h5.get('one4one')[0])
    if one4one:
        _logger.error('read_dpa42parray does not support one4one, yet')
        raise ValueError('read_dpa42parray does not support one4one, yet')

    _logger.debug('nslice = ' + str(nslice))
    _logger.debug('lslice = ' + str(lslice) + 'm')
    _logger.debug('sepslice = ' + str(sepslice)+'m')
    _logger.debug('zsep = ' + str(zsep))
    _logger.debug('npart = ' + str(npart))
    _logger.debug('nbins = ' + str(nbins))

    I = []
    for dset in h5:
            if dset.startswith('slice0'):
                I.append(h5[dset]['current'][0])
    I = np.array(I)

    dt = zsep * lslice / speed_of_light
    _logger.debug('dt = ' + str(dt) + 'sec')

    N_part_max = np.sum(I / I.max() * npart) # total maximum reasonable number of macroparticles of the same charge that can be extracted

    if N_part is not None:
        if N_part > N_part_max:
            N_part = int(np.floor(N_part_max))
    else:
        N_part = int(np.floor(N_part_max))
    _logger.debug('Number of particles max= ' + str(N_part))

    n_part_slice = (I / np.sum(I) * N_part).astype(int) #array of number of particles per new bin
    n_part_slice[n_part_slice > npart] = npart

    N_part_act = np.sum(n_part_slice) #actual number of particles
    _logger.debug('Number of particles actual= ' + str(N_part_act))

    dt = zsep * lslice / speed_of_light
    C = np.sum(I) * dt #total charge
    c = C / N_part_act # particle charge

    pick_i = [random.sample(range(npart), n_part_slice[i]) for i in range(nslice)]

    x = []
    y = []
    px = []
    py = []
    g = []
    ph = []
    s = []

    for dset in h5:
        if dset.startswith('slice0'):
            i = int(dset.strip('/slice'))-1
            if len(pick_i[i]) > 0:
                ph0 = h5[dset]['theta'].value[pick_i[i]]
                s.append(i*zsep*lslice + ph0 / 2 / np.pi * lslice)
                x.append(np.array(h5[dset]['x'].value[pick_i[i]]))
                px.append(np.array(h5[dset]['px'].value[pick_i[i]]))
                y.append(np.array(h5[dset]['y'].value[pick_i[i]]))
                py.append(np.array(h5[dset]['py'].value[pick_i[i]]))
                g.append(np.array(h5[dset]['gamma'].value[pick_i[i]]))

    p_array = ParticleArray()
    p_array.rparticles = np.empty((6, N_part_act))

    g = np.concatenate(g).ravel()
    g0 = np.mean(g) # average gamma
    p_array.E = g0 * m_e_GeV # average energy in GeV
    p0 = sqrt(g0**2-1) * m_e_eV / speed_of_light

    p_array.rparticles[0] = np.concatenate(x).ravel() # position in x in meters
    p_array.rparticles[1] = np.concatenate(px).ravel() / g0  # divergence in x
    p_array.rparticles[2] = np.concatenate(y).ravel() # position in x in meters
    p_array.rparticles[3] = np.concatenate(py).ravel() / g0  # divergence in x
    p_array.rparticles[4] = -np.concatenate(s).ravel()
    p_array.rparticles[5] = (g - g0) * m_e_eV / p0 / speed_of_light

    if fill_gaps:
        p_array.rparticles[4] -= np.random.randint(0, zsep, N_part_act) * lslice

    p_array.q_array = np.ones(N_part_act) * c

    h5.close()
    return p_array


def write_gen4_lat(lat, filePath, line_name='LINE', l=np.inf):
    from ocelot.cpbd.elements import Undulator, Drift, Quadrupole, UnknownElement
    _logger.info('writing genesis4 lattice')
    _logger.debug(ind_str + 'writing to ' + filePath)

    lat_str = []
    beamline = []
    ll=0

    lat_str.append('# generated with Ocelot\n')

    for element in lat.sequence:

        if ll >= l:
            break

        element_num = str(len(beamline) + 1).zfill(3)

        if hasattr(element,'l'):
            ll += element.l

        if isinstance(element, Undulator):
            element_name = element_num + 'UND'
            s = '{:}: UNDULATOR = {{lambdau = {:}, nwig = {:}, aw = {:.6f}}};'.format(element_name, element.lperiod, element.nperiods, element.Kx/sqrt(2))
#            print(s)

        elif isinstance(element, Drift):
            element_name = element_num + 'DR'
            s = '{:}: DRIFT = {{l={:}}};'.format(element_name, element.l)

        elif isinstance(element, Quadrupole):
            if element.k1>=0:
                element_name = element_num + 'QF'
            else:
                element_name =  element_num + 'QD'
            s = '{:}: QUADRUPOLE = {{l = {:}, k1 = {:.6f} }};'.format(element_name, element.l, element.k1)

        else:
            _logger.debug('Unknown element with length '+ str(element.l))
            continue

        beamline.append(element_name)
        lat_str.append(s)

    lat_str.append('')
    lat_str.append('{:}: LINE = {{{:}}};'.format(line_name, ','.join(beamline)))
    lat_str.append('\n# end of file\n')

    with open(filePath, 'w') as f:
        f.write("\n".join(lat_str))

def write_edist_hdf5(edist, filepath):
    with h5py.File(filepath, 'w') as h5:
        # h5 = h5py.File(filepath, 'w')
        h5.create_dataset('p', data=edist.g)
        h5.create_dataset('t', data=-edist.t)
        h5.create_dataset('x', data=edist.x)
        h5.create_dataset('y', data=edist.y)
        h5.create_dataset('xp', data=edist.xp)
        h5.create_dataset('yp', data=edist.yp)
        h5.create_dataset('charge', data = edist.charge())
        # h5.close()

def read_edist_hdf5(filepath, charge=None):

    edist = GenesisElectronDist()
    with h5py.File(filepath, 'r') as h5:

        edist.g =  h5.get('p')[:]
        edist.t =  -h5.get('t')[:]
        edist.x =  h5.get('x')[:]
        edist.y =  h5.get('y')[:]
        edist.xp =  h5.get('xp')[:]
        edist.yp =  h5.get('yp')[:]

        if charge is not None:
            charge = h5.get('charge')[:]

    edist.part_charge = charge / edist.g.size
    return edist
