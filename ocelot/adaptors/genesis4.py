import time
import os
import sys
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
from ocelot.adaptors.genesis import GenesisElectronDist  # tmp
from ocelot.common.ocelog import *
from ocelot.utils.launcher import *
from ocelot.cpbd.beam import Beam, BeamArray
from ocelot.cpbd.elements import Element, Drift, Quadrupole, Undulator, Marker
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.common.py_func import copy_this_script
from ocelot.rad.undulator_params import lambda2eV, eV2lambda

_logger = logging.getLogger(__name__)

# TODO: move to cpbd (enentually?)
class Chicane(Element):
    """
    chicane (implemented for Genesis4, not used by cpbd)
    l - length of chicane in [m],
    lb - length of an individual dipole [m],
    ld - drift between the outer and inner dipoles [m],
    delay - path length difference between the straight path and the actual trajectory [rad].
    """

    def __init__(self, l=0., lb=0., ld=0., delay=0., eid=None):
        # Element.__init__(self, eid)
        super(Chicane, self).__init__(eid=eid)
        self.l = l
        self.lb = lb
        self.ld = ld
        self.delay = delay


class Phaseshifter(Element):
    """
    phase shifter element (implemented for Genesis4, not used by cpbd)
    l - length of phase shifter in [m]
    phi = introduced phase shift in [rad] 
    """

    def __init__(self, l=0., phi=0., eid=None):
        # Element.__init__(self, eid)
        super(Phaseshifter, self).__init__(eid=eid)
        self.l = l
        self.phi = phi


class Genesis4Simulation:
    """
    Class for storing Genesis4 simulation parameters such as
    input file,
    settings for returning, plotting or cleaning output
    """

    def __init__(self, ginp=None, exp_dir='', **kwargs):

        self.mpiParameters = ''
        if ginp is not None:
            ginp.check_consistency_round2()
            self.ginp = ginp
        else:
            self.ginp = Genesis4Input()

        self.exp_dir = exp_dir  # directory in which simulation will take place, i.e. input and output files will be written

        # self.del_files = ('gout', 'fld', 'par')

        self.plot_output = kwargs.get('plot_output', True)
        self.return_out = kwargs.get('return_out', True)
        self.cleanup_afterwards = kwargs.get('cleanup_afterwards', False)
        self.zstop = kwargs.get('zstop', np.inf)
        self.exec_script_path = kwargs.get('exec_script_path', None)
        self.launcher = kwargs.get('launcher', get_genesis4_launcher())

    @property
    def root_name(self):
        return self.ginp.setup.rootname

    @root_name.setter
    def root_name(self, value):
        self.ginp.setup.rootname = value

    def filepath(self, filename):
        return os.path.join(self.exp_dir, filename)

    def input_filepath(self):
        return self.filepath(self.ginp.filename)

    def root_path(self):
        return self.filepath(self.root_name)

    def create_exp_dir(self):
        if not os.path.isdir(self.exp_dir):
            os.makedirs(self.exp_dir)

    def write_inp_file(self):
        with open(self.input_filepath(), 'w') as file:
            inp_string = str(self.ginp)
            file.write(inp_string)

    def write_lat_file(self):
        write_gen4_lat(self.ginp.attachments.lat, filepath=self.filepath(self.ginp.setup.lattice),
                       zstop=self.zstop)

    def write_dfl_files(self):
        """
        Writes down dfl files for &importfield Genesis1.3-version4 name lists
        :return:
        """
        for element_id in self.ginp.sequence.keys():
            if element_id.startswith('importfield'):
                impfld = self.ginp.sequence[element_id]
                if impfld._name_list_id in self.ginp.attachments.dfl.keys():
                    if impfld.file is not None:
                        filename = self.filepath(impfld.file)
                    else:
                        filename = self.filepath(impfld._name_list_id + 'dfl.h5')
                        impfld.file = filename
                    write_dfl4(self.ginp.attachments.dfl[impfld._name_list_id], filename)

        # for name_list_id, dfl in self.ginp.attachments.dfl: #TODO: iterate over namespaces, not attachments. If there is importing namespace without attachments 
        # if not isinstance(dfl, str):
        # write_dfl4(dfl, self.ginp.sequence[name_list_id].file)

    def write_dpa_files(self):
        """
        NOT IMPLEMENTED YET
        Writes down dpa files for &importbeam Genesis1.3-version4 name lists
        :return:
        """
        for element in self.ginp.sequence:
            if element.startswith('importbeam'):
                imppar = self.ginp.sequence[element]
                if imppar._name_list_id in self.ginp.attachments.dpa:
                    if imppar.file is not None:
                        filename = imppar.file
                    else:
                        filename = imppar._name_list_id + 'par.h5'
                        imppar.file = filename
                    # write_dpa4(self.attachments.dpa[imppar._name_list_id], filename)

    def write_beam_files(self):
        for element_id in self.ginp.sequence.keys():
            if element_id.startswith('beam'):
                impbeam = self.ginp.sequence[element_id]
                if impbeam._name_list_id in self.ginp.attachments.beam.keys():
                    file_path = self.filepath(impbeam._name_list_id + '.beam.h5')
                    write_beamtwiss_hdf5(self.ginp.attachments.beam[impbeam._name_list_id], file_path)

    def write_all_files(self):
        self.create_exp_dir()
        self.write_inp_file()
                             
        self.write_dfl_files()
        self.write_beam_files()
        self.write_lat_file()
        
        if self.exec_script_path is not None:
            copy_this_script(self.exec_script_path, self.exp_dir)
        # TODO: write all attachments

    def prepare_launcher(self, program_path='genesis4'):
        #self.launcher.program = program_path
        self.launcher.dir = self.exp_dir
        self.launcher.argument = ' ' + self.ginp.filename
        #self.mpiParameters = '-np 8'  # TODO: calculate and request minimum reasonable number of cores (to not ovewblow the window) S.S.

        # TODO: implement on Launcher level
        if sys.platform not in ["linux", "linux2"]:
            _logger.error('Only linux platform supported for given launcher (at the moment)')
            pass

    def run(self, launcher=None):
        self.prepare_launcher()
        if launcher is not None:
            self.launcher = launcher
        self.clean_output()
        self.write_all_files()
        
        self.launcher.launch()
        
        if self.return_out:
            out = read_gout4(self.root_path() + '.out.h5')
        else:
            out = None
        
        if self.cleanup_afterwards:
            self.clean_output()
        
        return out

    def clean_output(self):
        os.system('rm -rf {}*'.format(self.root_path()))


class Genesis4Input:
    """
    Ocelot container class to store Genesis4 input parameters
    """

    def __init__(self):
        self.sequence = {}  # sequence of namespaces; each value starts with a supported namespace name, e.g. "field" and if needed is followed by unique suffix after pound symbol, like "field#2"
        self.attachments = Genesis4Attachments()
        # self.add_prefix_attributes = ['setup.lattice', 'setup.rootname', 'filename', 'write.field', 'write.beam']
        self.filename = 'input.in'
        
        # If True, use "GENESIS 1.3" feature 'profile_file_multi', replacing multiple 'profile_file' blocks having shared 's' axis.
        # Feature added to G4 by C. Lechner (available in stable release 4.6.2).
        # By default, continue to emit multiple 'profile_file' blocks to avoid breaking compatibility with older "GENESIS 1.3" releases.
        self.x_use_profilefilemulti=False # True
        
        name_list_subclasses = {}
        for name_list_subclass in Genesis4NameList.__subclasses__():
            name_list_subclasses[name_list_subclass()._name_list_label] = name_list_subclass
        self.name_list_subclasses = name_list_subclasses

    def __str__(self):
        str_out = '# autogenerated with Ocelot\n\n'  # https://github.com/ocelot-collab/ocelot\n\n
        for name_list in self.sequence.values():
            # special procedure needed for verbatim blocks (not writing namelist start/end to have full control over what is written)
            if isinstance(name_list, Genesis4InfileVerbatimNL):
                str_out += '# infileverbatim block\n'
                if name_list.text is not None:
                    str_out += name_list.text
                    # To avoid parsing errors, ensure that there is a new line
                    # before the next block is written (\r for Windows)
                    if not (name_list.text[-1] in ['\n','\r']):
                        print('info: infileverbatim block added \\n line ending')
                        str_out += '\n' # UNIX style line ending
                continue

            # all other types of blocks, these have namelist starts/ends written
            str_out += '&' + name_list._name_list_label + '\n'
            str_out += str(name_list)
            str_out += '&end\n\n'
        return str_out
    
    def append_sequence(self, name_id=None):
        if '#' in name_id:
            temp, _ = name_id.split('#')
        else:
            temp = name_id
        self.sequence[name_id] = self.name_list_subclasses[temp](name_id)
        
    
    def make_sequence(self, name_list_ids=None):
        """
        fills a sequence attribute
        """
        self.sequence = {}
        
        # name_list_subclasses = {}
        # for name_list_subclass in Genesis4NameList.__subclasses__():
            # name_list_subclasses[name_list_subclass()._name_list_label] = name_list_subclass
        
        for name_list_id in name_list_ids:
            self.append_sequence(name_list_id)
            # if '#' in name_list_id:
                # temp, _ = name_list_id.split('#')
            # else:
                # temp = name_list_id
            # self.sequence[name_list_id] = self.name_list_subclasses[temp](name_list_id)
        #self.check_consistency()

    # def make_sequence_simple(self):
    #     attrs = ['setup', 'lattice', 'time', 'importdistribution', 'importbeam', 'beam', 'importfield', 'field',
    #              'track', 'write']
    #     self.make_sequence(attrs)

    @property
    def setup(self):
        for name_list in self.sequence.values():
            if isinstance(name_list, Genesis4SetupNL):
                return name_list

    def check_consistency(self):
        """
        checks if beam of field are not declared few times before tracking
        """
        beam_counter = 0
        field_counter = 0
        for id, name_list in self.sequence.items():
            if isinstance(name_list, (Genesis4BeamNL, Genesis4ImportBeamNL, Genesis4ImportDistributionNL)):
                beam_counter += 1
            if isinstance(name_list, (Genesis4FieldNL, Genesis4ImportFieldNL)):
                field_counter += 1
            if isinstance(name_list, (Genesis4TrackNL,)):
                if beam_counter != 1:
                    _logger.warning(
                        'Genesis4Input.sequence contains more or less than one declaration of beam')  # TODO: specify, 0 or >1
                if field_counter != 1:
                    _logger.warning('Genesis4Input.sequence contains more or less than one declaration of field')
                field_counter = 0
                beam_counter = 0
        # TODO: decide whether error raising is needed if co-exist e.g. importbeam and beam or importfield and field
    def check_consistency_round2(self):
        _logger.info('CL 2023-Aug-02: May need to re-implement this later')
		
    def populate_sequence_beam_array(self, beam_name_list_id, beam=None):
        # Info: profile_file_multi feature was implemented into
        # "GENESIS 1.3" v4 developer version by C. Lechner in Sept-2021.
        # It is included in the stable release 4.6.2
        # Use of this function is controlled by the switch 'x_use_profilefilemulti'.
        _logger.info('Populating beam_array to Genesis4Input.sequence')
        if self.x_use_profilefilemulti:
            _logger.info('    generating G4 input file with &profile_file_multi feature')
        if beam_name_list_id not in self.attachments.beam.keys():
            if beam is None:
                _logger.warning(ind_str + 'no beam_array to populate')
            else:
                self.attachments.beam[beam_name_list_id] = beam
                # Update simulation parameters: gamma0 in '&setup' and slen in '&time'
                self.sequence['setup'].gamma0 = np.mean(beam.g)
                # Special handling required since steady-state simulations typically do not use '&time' block
                got_time_block = 'time' in self.sequence
                if got_time_block:
                    self.sequence['time'].slen = np.amax(beam.s) - np.amin(beam.s)
                else:
                    _logger.warning('Sequence does not contain \'time\' element. This is not an issue for steady-state simulations.')
        ###
        ###
        new_sequence = {}
        for name_list_id, name_list in self.sequence.items():
            # If it is not the 'beam' we are looking for, just propagate to the new sequence
            if name_list_id != beam_name_list_id:
                new_sequence[name_list_id] = name_list
                continue
            ### We have found the beam to be manipulated ###
            # 1) generate profiles, using either &profile_file or &profile_file_multi (results in smaller GENESIS4 infiles)
            if not self.x_use_profilefilemulti:
                # For every attribute to load from HDF5 file:
                # Prepare structure that will eventually result in &profile_file to be written to GENESIS4 infile
                for attr in self.sequence[beam_name_list_id].__dict__.keys():
                    if attr.startswith('_'):
                        continue
                    profile_name_list_id = 'profile_file#' + attr + '_' + beam_name_list_id
                    new_sequence[profile_name_list_id] = Genesis4ProfileFileNL()
                    new_sequence[profile_name_list_id].label = attr + '_' + beam_name_list_id.replace('#', '')
                    new_sequence[profile_name_list_id].xdata = beam_name_list_id + '.beam.h5/s'
                    new_sequence[profile_name_list_id].ydata = beam_name_list_id + '.beam.h5/' + attr
            else:
                # Collect all attributes we need to load from the HDF5 file
                # -- and prepare structure that will result in single &profile_file_multi to be written
                load_attr = []
                for attr in self.sequence[beam_name_list_id].__dict__.keys():
                    if attr.startswith('_'):
                        continue
                    load_attr.append(attr)
                profile_name_list_id = 'profile_file#' + beam_name_list_id
                my_obj = Genesis4ProfileFileMultiNL()
                my_obj.file = beam_name_list_id + '.beam.h5'
                my_obj.label_prefix = beam_name_list_id.replace('#', '')
                my_obj.xdata = 's' # shared by all data sets specified by ydata parameter
                my_obj.ydata = ','.join(load_attr)
                new_sequence[profile_name_list_id] = my_obj

            # (2) add profiles to the beam
            new_sequence[name_list_id] = name_list
            for attr in self.sequence[beam_name_list_id].__dict__.keys():
                if attr.startswith('_'):
                    continue
                if self.x_use_profilefilemulti:
                    # using 'profile_file_multi'
                    profile_label = '@' + beam_name_list_id.replace('#', '') + '.' + attr
                else:
                    # traditional may using 'profile_file' elements
                    profile_label = '@' + attr + '_' + beam_name_list_id.replace('#', '')
                setattr(new_sequence[name_list_id], attr, profile_label)
        # done
        self.sequence = new_sequence


    def populate_sequence_beam(self, beam, name_list_id):
        """
        :param self: Genesis4Input
        :param beam: BeamArray or Beam object
        :param name_list_id:
        :return: None
        """
        if isinstance(beam, BeamArray):
            beam_pk = beam.pk()
            _logger.warning(
                'at the moment method beam_to_sequence_beam parses a single beam slice; peak current value is taken')
        elif isinstance(beam, Beam):
            beam_pk = beam
        else:
            raise TypeError('beam should be an instance of BeamArray or Beam')
        self.sequence['setup'].gamma0 = beam_pk.g
        self.sequence['time'].slen = np.amax(beam.s) - np.amin(beam.s)
        self.sequence[name_list_id].gamma = beam_pk.g  # from [GeV] to [units of the electron rest mass]
        self.sequence[name_list_id].delgam = beam_pk.dg
        self.sequence[name_list_id].current = beam_pk.I
        self.sequence[name_list_id].ex = beam_pk.emit_xn
        self.sequence[name_list_id].ey = beam_pk.emit_yn
        self.sequence[name_list_id].betax = beam_pk.beta_x
        self.sequence[name_list_id].betay = beam_pk.beta_y
        self.sequence[name_list_id].alphax = beam_pk.alpha_x
        self.sequence[name_list_id].alphay = beam_pk.alpha_y
        self.sequence[name_list_id].xcenter = beam_pk.x
        self.sequence[name_list_id].ycenter = beam_pk.y
        self.sequence[name_list_id].pxcenter = beam_pk.px
        self.sequence[name_list_id].pycenter = beam_pk.py
        self.sequence[name_list_id].bunch = 0
        self.sequence[name_list_id].bunchphase = 0
        self.sequence[name_list_id].emod = 0
        self.sequence[name_list_id].emodphase = 0

    # def add_prefixes(self, prefix_str):
    #     def rsetattr(obj, attr, val):
    #         pre, _, post = attr.rpartition('.')
    #         return setattr(rgetattr(obj, pre) if pre else obj, post, val)
    #
    #     def rgetattr(obj, attr, *args):
    #         def _getattr(obj, attr):
    #             return getattr(obj, attr, *args)
    #
    #         return functools.reduce(_getattr, [obj] + attr.split('.'))
    #
    #     for attr in self.add_prefix_attributes:
    #         try:
    #             oldvalue = rgetattr(self, attr)
    #             if oldvalue not in [None, '']:
    #                 # print(rgetattr(self, attr))
    #                 rsetattr(self, attr, prefix_str + oldvalue)
    #                 # print(rgetattr(self, attr))
    #         except:
    #             # _logger.debug...
    #             pass


class Genesis4Attachments:
    # TODO: think about possibility to attach and write several dfls, lattices, etc.
    def __init__(self):
        self.dfl = {}  # e.g. {key1: fld1, key2: fld2} where key is "_name_list_id" of Genesis4ImportFieldNL file should be written as <Genesis4ImportFieldNL.file>.fld
        self.dpa = {}  # e.g. {key1: par1, key2: par2} where key is "_name_list_id" of Genesis4ImportBeamNL file should be written as <Genesis4ImportBeamNL.file>.par
        self.beam = {}
        self.lat = {}  # e.g. {key1: lat1, key2: lat2}, where key is the value of setup.lattice


# def genesis4_input(input_sequence=tuple()):
#    str_out = ''
#    for name_list in input_sequence:
#        str_out += '&' + name_list._name_list_label + '\n'
#        str_out += str(name_list)
#        str_out += '&end\n\n'
#    return str_out

# %%
class Genesis4NameList(object):
    """
    Parent object for Genesis 4 namelists
    str(Genesis4NameList) yields the list text
    attributes not starting with underscore "_" are considered Genesis4NameList variables
    """

    # def attr_list(self):
    #     output_list = dir(self)
    #     for attr in dir(self):
    #         if attr.startswith('_') or callable(getattr(self, attr)):
    #             output_list.remove(attr)
    #     return output_list

    # def all_attrs_is_none(self):
    #    for attr in self.attr_list():
    #        if getattr(self, attr) is not None:
    #            return False
    #    return True

    # def attr_list(self):
    #     output_list = dir(self)
    #     for attr in dir(self):
    #         if attr.startswith('_') or callable(getattr(self, attr)):
    #             output_list.remove(attr)
    #     return output_list

    def name_list_id_set(self, name_list_id=None):
        self._name_list_id = name_list_id if name_list_id is not None else self._name_list_label

    def __str__(self):
        indentstr = '   ' # improve structure of G4 input files by indenting the contents of blocks
        str_out = ''
        for attr, value in self.__dict__.items():
            if attr.startswith('_'):
                continue
            if value is True:
                value = 'true'
            elif value is False:
                value = 'false'
            if value is not None:
                str_out += indentstr + attr + ' = ' + str(value) + '\n'
        return str_out


class Genesis4SetupNL(Genesis4NameList):

    def __init__(self,
                 name_list_id=None):
        # super().__init__()
        self._name_list_label = 'setup'
        self.name_list_id_set(name_list_id)
        self.rootname = 'output'  # name of all files generated by simulation
        self.lattice = None
        self.beamline = None
        self.gamma0 = None
        self.lambda0 = None
        self.delz = None
        self.seed = None
        self.npart = None
        self.nbins = None
        self.one4one = None
        self.shotnoise = None


class Genesis4AlterSetupNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'alter_setup'
        self.name_list_id_set(name_list_id)
        self.rootname = None
        self.beamline = None
        self.delz = None
        self.harmonic = None
        self.subharmonic = None
        self.resample = None


class Genesis4LatticeNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'lattice'
        self.name_list_id_set(name_list_id)
        self.zmatch = None
        self.element = None
        self.field = None
        self.value = None
        self.instance = None
        self.add = None


class Genesis4TimeNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'time'
        self.name_list_id_set(name_list_id)
        self.s0 = None
        self.slen = None
        self.sample = None
        self.time = None


class Genesis4ProfileConstNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'profile_const'
        self.name_list_id_set(name_list_id)
        self.label = None
        self.c0 = None


class Genesis4ProfileGaussNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'profile_gauss'
        self.name_list_id_set(name_list_id)
        self.label = None
        self.c0 = None
        self.s0 = None
        self.sig = None


class Genesis4ProfileStepNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'profile_step'
        self.name_list_id_set(name_list_id)
        self.label = None
        self.c0 = None
        self.s_start = None
        self.s_end = None


class Genesis4ProfilePolynomNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'profile_polynom'
        self.name_list_id_set(name_list_id)
        self.label = None
        self.c0 = None
        self.c1 = None
        self.c2 = None
        self.c3 = None
        self.c4 = None


class Genesis4ProfileFileNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'profile_file'
        self.name_list_id_set(name_list_id)
        self.label = None
        self.xdata = None
        self.ydata = None
        self.isTime = None
        self.reverse = None

# profile_file_multi was added in GENESIS v4 commit id 3fdcd81 (Sept 18, 2021)
class Genesis4ProfileFileMultiNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'profile_file_multi'
        self.name_list_id_set(name_list_id)
        self.file = None
        self.label_prefix = None
        self.xdata = None
        self.ydata = None
        self.isTime = None
        self.reverse = None


class Genesis4BeamNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'beam'
        self.name_list_id_set(name_list_id)
        self.gamma = None
        self.delgam = None
        self.current = None
        self.ex = None
        self.ey = None
        self.betax = None
        self.betay = None
        self.alphax = None
        self.alphay = None
        self.xcenter = None
        self.ycenter = None
        self.pxcenter = None
        self.pycenter = None
        self.bunch = None
        self.bunchphase = None
        self.emod = None
        self.emodphase = None


class Genesis4FieldNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'field'
        self.name_list_id_set(name_list_id)
        # setattr(self, 'lambda', None)
        self.lambda_ = None  # must be lambda, but it is taken by Python for lambda-expressions
        self.power = None
        self.phase = None
        self.waist_pos = None
        self.waist_size = None
        self.xcenter = None
        self.ycenter = None
        self.xangle = None
        self.yangle = None
        self.dgrid = None
        self.ngrid = None
        self.harm = None
        self.nx = None
        self.ny = None
        self.accumulate = None

    def __str__(self):
        return super(Genesis4FieldNL, self).__str__().replace('lambda_', 'lambda')

class Genesis4FieldManipulatorNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'field_manipulator'
        self.name_list_id_set(name_list_id)
        # setattr(self, 'lambda', None)
        self.harm = None
        self.scale_power = None

class Genesis4ImportDistributionNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'importdistribution'
        self.name_list_id_set(name_list_id)
        self.file = None
        self.sdds = None
        self.charge = None
        self.slicewidth = None
        self.output = None
        self.center = None
        self.gamma0 = None
        self.x0 = None
        self.y0 = None
        self.px0 = None
        self.py0 = None
        self.match = None
        self.betax = None
        self.betay = None
        self.alphax = None
        self.alphay = None


class Genesis4ImportFieldNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'importfield'
        self.name_list_id_set(name_list_id)
        self.file = None


class Genesis4ImportBeamNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'importbeam'
        self.name_list_id_set(name_list_id)
        self.file = None


class Genesis4EfieldNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'efield'
        self.name_list_id_set(name_list_id)
        self.rmax = None
        self.nz = None
        self.nphi = None
        self.ngrid = None


class Genesis4SponradNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'sponrad'
        self.name_list_id_set(name_list_id)
        self.seed = None
        self.doLoss = None
        self.doSpread = None


class Genesis4SortNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'sort'
        self.name_list_id_set(name_list_id)


class Genesis4WakeNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'wake'
        self.name_list_id_set(name_list_id)
        self.loss = None
        self.radius = None
        self.roundpipe = None
        self.conductivity = None
        self.relaxation = None
        self.material = None
        self.gap = None
        self.lgap = None
        self.hrough = None
        self.lrough = None
        self.transient = None
        self.ztrans = None


class Genesis4WriteNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'write'
        self.name_list_id_set(name_list_id)
        self.field = None
        self.beam = None


class Genesis4TrackNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'track'
        self.name_list_id_set(name_list_id)
        self.zstop = None
        self.output_step = None
        self.field_dump_step = None
        self.beam_dump_step = None
        self.sort_step = None

# CL, 2023-Jan: Class for text block that is copied to G4 input file
class Genesis4InfileVerbatimNL(Genesis4NameList):
    def __init__(self, name_list_id=None):
        # super().__init__()
        self._name_list_label = 'infileverbatim'
        self.name_list_id_set(name_list_id)
        self.text=None # This is the text that will be written to input file


class Genesis4ParticlesDump:
    """
    Genesis particle *.dpa files storage object
    Each particle record in z starts with the energy of all particles 
    followed by the output of the particle phases, 
    positions in x and y and the momenta in x and y. 
    The momenta are normalized to mc
    """

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


def gen4_lat_str(lat, line_name='LINE', zstop=np.inf):
    """
    Generates a string of lattice 
    in Genesis4 format
    from an ocelot lattice object
    """
    from ocelot.cpbd.elements import Undulator, Drift, Quadrupole, UnknownElement
    lat_str = []
    beamline = []
    location = 0

    for element in lat.sequence:

        if location >= zstop:
            break

        element_num = line_name + '_' + str(len(beamline) + 1).zfill(3)

        if hasattr(element, 'l'):
            location += element.l
        else:
            _logging.warning('[beta] element had no length: {:}'.format(str(element)))

        if isinstance(element, Undulator):
            element_name = element_num + 'UND'

            if element.Kx == 0 or element.Ky == 0:
                is_helical = 'false'
            elif element.Kx == element.Ky:
                is_helical = 'true'
            else:
                _logger.warning(
                    'undulator element with different non-zero Kx and Ky; not implemented; setting to helical')
                is_helical = 'true'
            aw = np.sqrt((element.Kx ** 2 + element.Ky ** 2) / 2)
            # TODO: implement ax,ay,kx,ky,gradx,grady
            s = '{:}: UNDULATOR = {{lambdau = {:}, nwig = {:}, aw = {:.6f}, helical = {:}}};'.format(element_name,
                                                                                                     element.lperiod,
                                                                                                     element.nperiods,
                                                                                                     aw, is_helical)

        elif isinstance(element, Drift):
            element_name = element_num + 'DR'
            s = '{:}: DRIFT = {{l={:}}};'.format(element_name, element.l)

        elif isinstance(element, Quadrupole):
            # TODO: add dx and dy
            if element.k1 >= 0:
                element_name = element_num + 'QF'
            else:
                element_name = element_num + 'QD'
            s = '{:}: QUADRUPOLE = {{l = {:}, k1 = {:.6f} }};'.format(element_name, element.l, element.k1)

        elif isinstance(element, Chicane):
            element_name = element_num + 'CH'
            s = '{:}: CHICANE = {{l = {:}, lb = {:}, ld = {}, delay = {:.5e} }};'.format(element_name, element.l,
                                                                                         element.lb, element.ld,
                                                                                         element.delay)

        elif isinstance(element, Marker):
            element_name = element_num + 'M'
            m_dumpfield = getattr(element, 'dumpfield', 0)
            m_dumpbeam = getattr(element, 'dumpbeam', 0)
            m_sort = getattr(element, 'sort', 0)
            m_stop = getattr(element, 'stop', 0)
            s = '{:}: MARKER = {{dumpfield = {:}, dumpbeam = {:}, sort = {:}, stop = {:} }};'.format(element_name,
                                                                                                     m_dumpfield,
                                                                                                     m_dumpbeam, m_sort,
                                                                                                     m_stop)

        elif isinstance(element, Phaseshifter):
            element_name = element_num + 'PH'
            s = '{:}: PHASESHIFTER = {{l = {:}, phi = {:}}};'.format(element_name, element.l, element.phi)

        else:
            _logger.warning('Unknown element {} with length {}\n replacing with drift'.format(str(element), element.l))
            element_name = element_num + 'UNKNOWN'
            s = '{:}: DRIFT = {{l={:}}};'.format(element_name, element.l)
            continue

        beamline.append(element_name)
        lat_str.append(s)

    lat_str.append('')
    lat_str.append('{:}: LINE = {{{:}}};'.format(line_name, ','.join(beamline)))
    lat_str.append('')  # some white space before the following beamline description (if any)
    lat_str.append('')
    lat_str.append('')
    lat_str = "\n".join(lat_str)
    _logger.debug(ind_str + lat_str)
    return lat_str


def write_gen4_lat(lattices, filepath, zstop=np.inf):
    """
    Writing lattice file for Genesis1.3-version4 simulations
    :param lattices: dictionary: {'line_name': ocelot.cpbd.magnetic_lattice.MagneticLattice(), ...}
    :param filepath: str: path to the file in which lattices information will be writen
    :param zstop: dict with active lengths for each lattice in lattices dictionary (can be a double if len(lattices)==1)
    :return:
    """
    _logger.info('writing genesis4 lattice')
    _logger.debug(ind_str + 'writing to ' + filepath)
    _logger.debug('lattices = {}'.format(str(lattices)))
    f = open(filepath, 'w')  # erasing file content
    f.write('# generated with Ocelot\n')

    if not isinstance(zstop, dict):
        if type(lattices) != dict:
            raise TypeError("type(lattices) != dict: lattices should be a dictionary with keys the same as lattices")
        else:
            lat_name = [key for key in lattices.keys()][0]
            zstop = {lat_name: zstop}


    for line_name, lat in zip(lattices.keys(), lattices.values()):
        _logger.debug(ind_str + "line={}, lat={}, zstop={}".format(line_name, type(lat), zstop.get(line_name, np.inf)))
        lat_str = gen4_lat_str(lat, line_name=line_name, zstop=zstop.get(line_name, np.inf))
        f.write(lat_str)

    f.write('\n# end of file')
    f.close()
    _logger.debug(ind_str + 'done')


class Genesis4Output:
    """
    Genesis input files storage object
    """

    def __init__(self):
        self.h5 = None  # hdf5 pointer

    def close(self):
        self.h5.close()

    @property
    def filePath(self):
        return self.h5.filename
    
    @property
    def filePathRoot(self):
        return self.h5.filename.replace('out.h5','*')

    def fileName(self):
        return os.path.basename(self.h5.filename)

    @property
    def nZ(self):
        return self.z.size

    @property
    def nSlices(self):
        # FIXME 2022-09-16: Newer version of GENESIS v4 can re-compute the current at each integration step (current is matrix then)
        # return self.h5['Beam/current'].size
        return self.h5['Global/s'].size

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
    def rad_power(self): #TODO: refactor from property back to attribute?
        # if harm == 1:
        return self.h5['Field/power']
        # else:
            # return self.h5['Field{:}/power'.format(harm)]

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
            intens = intens[zi, :]
            phase = phase[zi, :]

        # not scaled properly!!!!
        field = np.sqrt(intens[:]) * np.exp(1j * phase[:])
        return field

    def calc_spec(self, zi=None, loc='near', npad=1, estimate_ph_sp_dens=1):
        #TODO: add harm parameter to calculate harmonics
        field = self.rad_field(zi=zi, loc=loc)
        axis = field.ndim - 1
        
        spec = np.abs(np.fft.fft(field, axis=axis)) ** 2
        spec = np.fft.fftshift(spec, axes=axis)
        
        scale_ev = h_eV_s * speed_of_light * (
                np.fft.fftfreq(self.nSlices, d=self.s[1] - self.s[0]) + 1 / self.lambdaref)
        scale_ev = np.fft.fftshift(scale_ev)
        
        if estimate_ph_sp_dens:
            tt = np.trapz(spec, scale_ev, axis=axis)
            if axis == 1:
                tt[tt == 0] = np.inf
                spec *= (self.n_photons / tt)[:, np.newaxis]
            else:
                if tt == 0:
                    tt = np.inf
                spec *= (self.n_photons[zi] / tt)

        return scale_ev, spec

    def wig(self, z=np.inf):
        return wigner_out(self, z=z, method='mp', debug=1)

    def close(self):
        _logger.warning('closing h5 file {}'.format(self.filePath))
        self.h5.close()


def get_genesis4_launcher(launcher_program='genesis4', launcher_argument='', mpiParameters='-x PATH -x MPI_PYTHON_SITEARCH -x PYTHONPATH'):
    """
    Returns MpiLauncher() object for given program
    """
    host = socket.gethostname()

    launcher = MpiLauncher()
    launcher.program = launcher_program
    launcher.argument = launcher_argument
    launcher.mpiParameters = mpiParameters # added -n
    # launcher.program = '/data/netapp/xfel/products/genesis/genesis'
    # launcher.argument = ' < tmp.cmd | tee log'

    return launcher


### USE FOR PILLAGE
def run_genesis4(inp, launcher, *args, **kwargs):
    """
    Main function for executing Genesis code
    inp               - GenesisInput() object with genesis input parameters
    launcher          - MpiLauncher() object obtained via get_genesis_launcher() function
    """
    # import traceback
    _logger.info('Starting genesis v4 preparation')
    # _logger.warning(len(traceback.extract_stack()))

    # create experimental directory
    if inp.run_dir is None and inp.exp_dir is None:
        raise ValueError('run_dir and exp_dir are not specified!')

    if inp.run_dir is None:
        if inp.exp_dir[-1] != os.path.sep:
            inp.exp_dir += os.path.sep
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
    """
    Reads Genesis1.3 v4 output file with the link to the hdf5 file as out.h5
    to close the file, use out.h5.close()
    
    :param filePath: string, absolute path to .out file
    :returns: Genesis4Output
    """

    _logger.info('reading gen4 .out file')
    _logger.warning(ind_str + 'in beta')
    _logger.debug(ind_str + 'reading from ' + filePath)

    out = Genesis4Output()
    try:
        out.h5 = h5py.File(filePath, 'r')
    except Exception:
        _logger.error(ind_str + 'no such file ' + filePath)
        raise

    if not ('Global' and 'Lattice' in out.h5):
        _logger.warning('Not a valid .out file')

    vvv = [int(out.h5['Meta/Version/' + name][0]) for name in ['Major', 'Minor', 'Revision']]
    _logger.debug(ind_str + 'Genesis v{}.{}.{}'.format(*vvv))

    out.z = out.h5['Lattice/zplot'][:]
    out.zlat = out.h5['Lattice/z'][:]

    _logger.debug(ind_str + 'z[0]   , z[-1]   , len(z)    = {}  {}  {}'.format(out.z[0], out.z[-1], len(out.z)))
    _logger.debug(
        ind_str + 'zlat[0], zlat[-1], len(zlat) = {}  {}  {}'.format(out.zlat[0], out.zlat[-1], len(out.zlat)))

    if 'time' in out.h5['Global'] and out.h5['Global/time'][0] == 1:
        out.tdp = True
        _logger.debug(ind_str + 'tdp = True')
    else:
        out.tdp = False
        _logger.debug(ind_str + 'tdp = False')

    if not out.tdp and out.h5['Beam/current'].size > 1:  # same as out.nSlices > 1:
        _logger.error('time independent simulation with {} slices'.format(out.nSlices))

    if out.tdp:
        if 's0' in out.h5['Global']:
            s0 = out.h5['Global/s0']
            _logger.debug(ind_str + 's0 = {}'.format(s0))
        else:
            s0 = 0
        sn = out.h5['Global/s'].size
        out.s = np.linspace(s0, out.h5['Global/slen'][:][0] + s0, sn)
        _logger.debug(ind_str + 's[0], s[-1], len(s) = {}  {}  {}'.format(out.s[0], out.s[-1], sn))

    _logger.debug(ind_str + 'done')

    return out


def read_dfl4(filePath):
    """
    Reads Genesis1.3 v4 radiation output file
    
    :param filePath: string, absolute path to .fld file
    :returns: RadiationField
    """
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
            ncar = int(np.sqrt(h5.get('slice000001/field-real').size))  # fix?
        zsep = int(sepslice / lambdaref)
        l_total = sepslice * nslice

        _logger.debug(ind_str + 'nslice = {}'.format(nslice))
        _logger.debug(ind_str + 'sepslice (dz) = {}'.format(sepslice))
        _logger.debug(ind_str + 'lambdaref (xlamds) = {}'.format(lambdaref))
        _logger.debug(ind_str + 'gridsize (dx & dy) = {}'.format(gridsize))
        _logger.debug(ind_str + '')
        _logger.debug(ind_str + 'zsep = {}'.format(zsep))
        _logger.debug(ind_str + 'Nx & Ny = {}'.format(ncar))
        _logger.debug(ind_str + 'Lx & Ly = {}'.format(ncar * gridsize))
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
            _logger.log(5, 'slice{:06d}'.format(i + 1))
            h5.create_dataset('slice{:06d}/field-real'.format(i + 1), data=np.real(dfl.fld[i]).flatten())
            h5.create_dataset('slice{:06d}/field-imag'.format(i + 1), data=np.imag(dfl.fld[i]).flatten())

    _logger.debug(ind_str + 'done')


def read_dpa4(filePath, start_slice=0, stop_slice=np.inf, estimate_npart=0, partskip=1):
    """
    Reads Genesis1.3 v4 particle output file
    
    :param filePath: string, absolute path to .par file
    :partskip: int, allows to read every nth particle only to save space. (default=1)
    :estimate_npart: reads the header and estimates number of particles and type of simulation. 0 - no estimation. 1 - normal estimation with console output. 2 - estimation only, omitting the actual parsing of particles
    :returns: Genesis4ParticlesDump
    """

    _logger.info('reading gen4 .dpa file')
    _logger.warning(ind_str + 'in beta')
    _logger.debug(ind_str + 'reading from ' + filePath)

    start_time = time.time()

    _logger.debug(ind_str + 'start_slice : stop_slice = {} : {}'.format(start_slice, stop_slice))
    
    sl = slice(None, None, partskip)
    
    if partskip > 1:#TODO: fix partskip for one4one=False
        _logger.debug(ind_str + 'reading every {}st/th particle'.format(partskip))
    
    with h5py.File(filePath, 'r') as h5:

        one4one = bool(h5.get('one4one')[0])
        _logger.info(ind_str + 'one4one = {}'.format(bool(one4one)))

        nslice = int(h5.get('slicecount')[0])
        nbins = int(h5.get('beamletsize')[0])
        lslice = h5.get('slicelength')[0]
        sepslice = h5.get('slicespacing')[0]

        if not one4one:
            npart = int(h5.get('slice000001/gamma').size)  # fix?
            _logger.debug(ind_str + 'npart = {}'.format(npart))
        else:
            if estimate_npart > 0:
                I_tmp_full = []
                I_tmp_wind = []
                for dset in h5:
                    if dset.startswith('slice') and type(h5[dset]) == h5py._hl.group.Group:
                        slicen = int(dset.replace('slice', ''))
                        I_tmp_full.append(h5[dset]['current'][:])
                        if slicen >= start_slice and slicen <= stop_slice:
                            I_tmp_wind.append(h5[dset]['current'][:])
                npart_full_tmp = np.sum(I_tmp_full) * lslice / speed_of_light / q_e
                npart_wind_tmp = np.sum(I_tmp_wind) * lslice / speed_of_light / q_e
                _logger.info(ind_str + 'estimated npart = {:}M'.format(npart_full_tmp / 1e6))
                _logger.info(ind_str + 'estimated npart to be downloaded = {:}.M'.format(npart_wind_tmp / partskip / 1e6))

        # filePath = h5.filename

        zsep = int(sepslice / lslice)
        # if one4one:
            # l_total = lslice * nslice
        # else:
            # l_total = lslice * zsep * nslice
        l_total = sepslice * nslice
        
        
        _logger.debug(ind_str + 'nslice = {}'.format(nslice))
        _logger.debug(ind_str + 'nbins = {}'.format(nbins))
        _logger.debug(ind_str + '')
        _logger.debug(ind_str + 'lslice (aka xlamds) = {} m'.format(lslice))
        _logger.debug(ind_str + 'sepslice = {} m'.format(sepslice))
        _logger.debug(ind_str + 'zsep = {}'.format(zsep))
        _logger.debug(ind_str + 'Ls_total = {}'.format(l_total))
        
        dpa = Genesis4ParticlesDump()
        dpa.nslice = nslice
        dpa.lslice = lslice
        dpa.sepslice = sepslice
        dpa.nbins = nbins
        dpa.zsep = zsep
        dpa.filePath = filePath
        dpa.one4one = one4one
        
        _logger.debug(ind_str + 'nslice = {}'.format(nslice))
        _logger.debug(ind_str + 'nbins = {}'.format(nbins))
        _logger.debug(ind_str + 'lslice (aka xlamds) = {} m'.format(lslice))
        _logger.debug(ind_str + 'sepslice = {} m'.format(sepslice))
        _logger.debug(ind_str + 'zsep = {}'.format(zsep))
        _logger.debug(ind_str + 'Ls_total = {}'.format(l_total))
        
        if estimate_npart == 2:
            _logger.debug(ind_str + 'done in %.2f seconds' % (time.time() - start_time))
            return dpa
            
        
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
        slicenum = [] #slice where the particle belongs
        slicelist = []

        _logger.debug(ind_str + 'reading slices between {} and {}'.format(start_slice, stop_slice))
        _logger.debug(2 * ind_str + '({} out of {})'.format(stop_slice - start_slice, nslice))
        for dset in sorted(h5):
            # _logger.log(5, ind_str + dset)
            # _logger.log(5, ind_str + str(type(h5[dset])))
            if dset.startswith('slice') and type(h5[dset]) == h5py._hl.group.Group:
                slicen = int(dset.replace('slice', ''))
                _logger.log(5, 2 * ind_str + 'slice number {}'.format(slicen))
                if slicen >= start_slice and slicen <= stop_slice:
                    _logger.log(5, 2 * ind_str + 'processing')
                    I.append(h5[dset]['current'][:])
                    ph.append(h5[dset]['theta'][sl])
                    _logger.log(5, 2 * ind_str + '{} particles'.format(h5[dset]['x'].size))
                    # s.append(s0 + ph0 / 2 / np.pi * lslice)
                    x.append(h5[dset]['x'][sl])
                    px.append(h5[dset]['px'][sl])
                    y.append(h5[dset]['y'][sl])
                    py.append(h5[dset]['py'][sl])
                    g.append(h5[dset]['gamma'][sl])
                    npartpbi = h5[dset]['gamma'][sl].size
                    npartpb.append(npartpbi)
                    slicenum.extend(list(np.ones(npartpbi, dtype=int)*slicen))
                    slicelist.append(slicen)
                    # ph.append(ph0)
                    # s0 += sepslice
        _logger.debug(2 * ind_str + 'done')
    
    if one4one:
        _logger.debug(ind_str + 'flattening arrays')
        _logger.debug(ind_str + 'writing to dpa object')
        # _logger.log(5, 2*ind_str + 'x.shape {}'.format(np.shape(x)))
        dpa.x = np.hstack(x)
        dpa.px = np.hstack(px)
        dpa.y = np.hstack(y)
        dpa.py = np.hstack(py)
        dpa.g = np.hstack(g)
        dpa.I = np.hstack(I)
        dpa.ph = np.hstack(ph)
        dpa.slicenum = np.array(slicenum)
        
        _logger.log(5, 2 * ind_str + 'dpa.x.shape {}'.format(dpa.x.shape))

        dpa.npartpb = np.array(npartpb).flatten()

        npart = dpa.x.size
        _logger.info(ind_str + 'npart = {}'.format(npart))
        
        dpa.one4one_qp = q_e * partskip #particle charge, assuming that npart(original) >> partskip
    else:
        npartpb = int(npart / nbins)
        # _logger.debug(ind_str + 'not one4one:')
        _logger.debug(2 * ind_str + 'shape of arrays (nslice, npart) {}'.format(np.array(x).shape))

        if np.array(x).shape[0] < nslice:
            nslice = np.array(x).shape[0]

        if nslice == 0:
            _logger.error(2 * ind_str + 'nslice == 0')
            return

        _logger.debug(
            2 * ind_str + 'reshaping to (nslice, nbins, npart/bin) ({}, {}, {})'.format(nslice, nbins, npartpb))
        _logger.debug(ind_str + 'writing to dpa object')
        dpa.x = np.array(x).reshape((nslice, nbins, npartpb), order='F')
        dpa.px = np.array(px).reshape((nslice, nbins, npartpb), order='F')
        dpa.y = np.array(y).reshape((nslice, nbins, npartpb), order='F')
        dpa.py = np.array(py).reshape((nslice, nbins, npartpb), order='F')
        dpa.ph = np.array(ph).reshape((nslice, nbins, npartpb), order='F')
        dpa.g = np.array(g).reshape((nslice, nbins, npartpb), order='F')
        dpa.I = np.array(I).flatten()
    
    
    dpa.npart = npart
    dpa.slicelist = np.array(slicelist)

    _logger.debug(ind_str + 'done in %.2f seconds' % (time.time() - start_time))
    
    return dpa


def dpa42edist(dpa, n_part=None, fill_gaps=False):
    """
    Convert Genesis1.3 v4 particle output file to ocelot edist object
    
    :param dpa: GenesisParticlesDump
    :param n_part: desired approximate number of particles in edist
    :param fill_gaps: dublicates buckets into gaps
    :returns: GenesisElectronDist
    """

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
        t0 = np.hstack([np.ones(n) * i for i, n in enumerate(dpa.npartpb)]) * dpa.sepslice / speed_of_light
        t0 += dpa.ph / 2 / np.pi * dpa.sepslice / speed_of_light

        edist = GenesisElectronDist()

        #C = npart * q_e
        C = npart * dpa.one4one_qp

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

        _logger.debug(2 * ind_str + 'done')

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

        # n_part_bin = (I / np.sum(I) * n_part / nbins).astype(int)
        # n_part_bin[n_part_bin > npart_bin] = npart_bin
        # print(n_part_bin.max())
        n_part_slice = (I / np.sum(I) * n_part).astype(int)
        n_part_slice[n_part_slice > npart] = npart

        _logger.debug(ind_str + 'max particles/slice = {}'.format(n_part_slice.max()))

        # print(n_part_slice.max())

        # pick_i = random.sample(range(n_part), n_part_slice[i])
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
        t0 = np.arange(nslice)[:, np.newaxis] * lslice * zsep / speed_of_light

        t = t1 + t0
        # _logger.debug(2*ind_str + 'complete')

        edist = GenesisElectronDist()
        # g1 = np.array([])
        # x1 = np.array([])
        # y1 = np.array([])

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
            _logger.log(5, "slice {} pick {} from {}".format(i, n_part_slice[i], npart))
            # _logger.log(1, "  {}".format(str(pick_i)))
            # edist.g = np.append(edist.g, g[i, pick_i])
            # edist.xp = np.append(edist.xp, px[i, pick_i])
            # edist.yp = np.append(edist.yp, py[i, pick_i])
            # edist.x = np.append(edist.x, x[i, pick_i])  
            # edist.y = np.append(edist.y, y[i, pick_i])
            # edist.t = np.append(edist.t, t[i, pick_i])
            _logger.log(5, "  shape before: {}".format(str(g.shape)))
            _logger.log(5, "  picking g[{}, len({})]".format(i, len(pick_i)))
            _logger.log(5, "  shape after: {}".format(str(np.shape(g[i, pick_i]))))

            gi.extend(g[i, pick_i])
            xpi.extend(px[i, pick_i])
            ypi.extend(py[i, pick_i])
            xi.extend(x[i, pick_i])
            yi.extend(y[i, pick_i])
            ti.extend(t[i, pick_i])

            _logger.log(5, "  shape after append: {}".format(str(len(gi))))

        edist.g = np.array(gi)  # .flatten()
        edist.xp = np.array(xpi)  # .flatten()
        edist.yp = np.array(ypi)  # .flatten()
        edist.x = np.array(xi)  # .flatten()
        edist.y = np.array(yi)  # .flatten()
        edist.t = np.array(ti)  # .flatten()
        _logger.log(5, "  shape after flattening: {}".format(str(edist.g.shape)))

        if fill_gaps:
            _logger.info(ind_str + 'randomly distributing particles in time between buckets')
            edist.t += np.random.randint(0, dpa.zsep, edist.t.size) * lslice / speed_of_light

        # particles_kept_ratio = 1

    edist.part_charge = C / edist.len()
    _logger.debug('')

    _logger.debug(ind_str + 'edist bunch charge = {} C'.format(C))
    _logger.debug(
        ind_str + 'edist particle charge = {} C (~{:.2f}*q_e)'.format(edist.part_charge, edist.part_charge / q_e))
    _logger.debug(ind_str + 'edist n_part = {}'.format(edist.len()))

    if hasattr(dpa, 'filePath'):
        edist.filePath = dpa.filePath + '.edist'

    _logger.debug(ind_str + 'done')

    return edist


def read_dpa42parray(filePath, N_part=None, fill_gaps=True):
    _logger.info('reading gen4 .dpa file into parray')
    _logger.warning(ind_str + 'in beta, fix situation with harmonic jump resulting in sepslice > lslice')
    _logger.debug(ind_str + 'reading from ' + filePath)

    import random
    N_part = None

    # N_part = 100000
    fill_gaps = True
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
    _logger.debug('sepslice = ' + str(sepslice) + 'm')
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

    N_part_max = np.sum(
        I / I.max() * npart)  # total maximum reasonable number of macroparticles of the same charge that can be extracted

    if N_part is not None:
        if N_part > N_part_max:
            N_part = int(np.floor(N_part_max))
    else:
        N_part = int(np.floor(N_part_max))
    _logger.debug('Number of particles max= ' + str(N_part))

    n_part_slice = (I / np.sum(I) * N_part).astype(int)  # array of number of particles per new bin
    n_part_slice[n_part_slice > npart] = npart

    N_part_act = np.sum(n_part_slice)  # actual number of particles
    _logger.debug('Number of particles actual= ' + str(N_part_act))

    dt = zsep * lslice / speed_of_light
    C = np.sum(I) * dt  # total charge
    c = C / N_part_act  # particle charge

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
            i = int(dset.strip('/slice')) - 1
            if len(pick_i[i]) > 0:
                ph0 = h5[dset]['theta'].value[pick_i[i]]
                s.append(i * zsep * lslice + ph0 / 2 / np.pi * lslice)
                x.append(np.array(h5[dset]['x'].value[pick_i[i]]))
                px.append(np.array(h5[dset]['px'].value[pick_i[i]]))
                y.append(np.array(h5[dset]['y'].value[pick_i[i]]))
                py.append(np.array(h5[dset]['py'].value[pick_i[i]]))
                g.append(np.array(h5[dset]['gamma'].value[pick_i[i]]))

    p_array = ParticleArray()
    p_array.rparticles = np.empty((6, N_part_act))

    g = np.concatenate(g).ravel()
    g0 = np.mean(g)  # average gamma
    p_array.E = g0 * m_e_GeV  # average energy in GeV
    p0 = np.sqrt(g0 ** 2 - 1) * m_e_eV / speed_of_light

    p_array.rparticles[0] = np.concatenate(x).ravel()  # position in x in meters
    p_array.rparticles[1] = np.concatenate(px).ravel() / g0  # divergence in x
    p_array.rparticles[2] = np.concatenate(y).ravel()  # position in x in meters
    p_array.rparticles[3] = np.concatenate(py).ravel() / g0  # divergence in x
    p_array.rparticles[4] = -np.concatenate(s).ravel()
    p_array.rparticles[5] = (g - g0) * m_e_eV / p0 / speed_of_light

    if fill_gaps:
        p_array.rparticles[4] -= np.random.randint(0, zsep, N_part_act) * lslice

    p_array.q_array = np.ones(N_part_act) * c

    h5.close()
    return p_array


# MOVE TO BEAM
def write_edist_hdf5(edist, filepath, exist_ok=True):
    _logger.info('writing electron distribution to {}'.format(filepath))
    os.makedirs(filepath[:filepath.rindex(os.path.sep)], exist_ok=True)
    with h5py.File(filepath, 'w') as h5:
        h5.create_dataset('p', data=edist.g)
        h5.create_dataset('t', data=-edist.t)
        h5.create_dataset('x', data=edist.x)
        h5.create_dataset('y', data=edist.y)
        h5.create_dataset('xp', data=edist.xp)
        h5.create_dataset('yp', data=edist.yp)
        # h5.create_dataset('charge', data=[edist.charge()])
        # h5.create_dataset('part_charge', data=[edist.part_charge])
    _logger.debug(ind_str + 'done')


def read_edist_hdf5(filepath, charge=None):
    _logger.info('reading electron distribution from {}'.format(filepath))
    edist = GenesisElectronDist()
    with h5py.File(filepath, 'r') as h5:

        edist.g = h5.get('p')[:]
        edist.t = -h5.get('t')[:]
        edist.x = h5.get('x')[:]
        edist.y = h5.get('y')[:]
        edist.xp = h5.get('xp')[:]
        edist.yp = h5.get('yp')[:]
        npart = len(edist.g)
        
        try:
            edist.part_charge = h5.get('part_charge')
            _logger.debug('retrieved edist charge per particle = {}'.format(edist.part_charge))
        except TyprError:
            edist.part_charge = None
            _logger.warn('no part_charge in edist')
        
        if charge != None:
            edist.part_charge = charge / npart
            _logger.debug('particle charge is overriden to {}'.format(edist.part_charge))
        else:
            if edist.part_charge == None:
                _logger.warn('particle charge was not provided neither in the file nor as an argument')
        
    return edist


def write_beamtwiss_hdf5(beam, filepath):
    _logger.info('writing electron beam (twiss) to {}'.format(filepath))
    with h5py.File(filepath, 'w') as h5:
        h5.create_dataset('s', data=beam.s)
        h5.create_dataset('gamma', data=beam.g)
        h5.create_dataset('delgam', data=beam.dg)
        h5.create_dataset('current', data=beam.I)
        h5.create_dataset('ex', data=beam.emit_xn)
        h5.create_dataset('ey', data=beam.emit_yn)
        h5.create_dataset('betax', data=beam.beta_x)
        h5.create_dataset('betay', data=beam.beta_y)
        h5.create_dataset('alphax', data=beam.alpha_x)
        h5.create_dataset('alphay', data=beam.alpha_y)
        h5.create_dataset('xcenter', data=beam.x)
        h5.create_dataset('ycenter', data=beam.y)
        h5.create_dataset('pxcenter', data=beam.px)
        h5.create_dataset('pycenter', data=beam.py)
        if hasattr(beam, 'bunch'):
            h5.create_dataset('bunch', data=beam.bunch)
        else:
            h5.create_dataset('bunch', data=np.zeros_like(beam.s))
        if hasattr(beam, 'bunchphase'):
            h5.create_dataset('bunchphase', data=beam.bunchphase)
        else:
            h5.create_dataset('bunchphase', data=np.zeros_like(beam.s))
        if hasattr(beam, 'emod'):
            h5.create_dataset('emod', data=beam.emod)
        else:
            h5.create_dataset('emod', data=np.zeros_like(beam.s))
        if hasattr(beam, 'emodphase'):
            h5.create_dataset('emodphase', data=beam.emodphase)
        else:
            h5.create_dataset('emodphase', data=np.zeros_like(beam.s))
        if hasattr(beam, 'eloss'):
            h5.create_dataset('eloss', data=beam.eloss)
        else:
            h5.create_dataset('eloss', data=np.zeros_like(beam.s))
    _logger.debug(ind_str + 'done')


def read_beamtwiss_hdf5(filepath):
    _logger.info('reading electron beam (twiss) from {}'.format(filepath))
    beam = BeamArray()
    with h5py.File(filepath, 'r') as h5:
        beam.s = h5.get('s')[:]
        beam.g = h5.get('gamma')[:]
        beam.dg = h5.get('delgam')[:]
        beam.I = h5.get('current')[:]
        beam.emit_xn = h5.get('ex')[:]
        beam.emit_yn = h5.get('ey')[:]
        beam.beta_x = h5.get('betax')[:]
        beam.beta_y = h5.get('betay')[:]
        beam.alpha_x = h5.get('alphax')[:]
        beam.alpha_y = h5.get('alphay')[:]
        beam.x = h5.get('xcenter')[:]
        beam.y = h5.get('ycenter')[:]
        beam.px = h5.get('pxcenter')[:]
        beam.py = h5.get('pycenter')[:]
        #try reading eloss,emodphase,emod etc.
    _logger.debug(ind_str + 'done')
    return beam
