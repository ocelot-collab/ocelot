import os
import functools

import ocelot
from ocelot.gui.beam_plot import *
from ocelot.adaptors.genesis4 import *
from ocelot.gui.genesis4_plot import *
from copy import deepcopy

_logger = logging.getLogger(__name__)
ocelog.setLevel(logging.INFO)

sim_directory = r'/gpfs/exfel/data/group/wp72/sserkez/projects/XFEL/2021_02_Genesis4_test/2023_02_Genesis4_MWE/'
sim_subdirectory = 'with_harmonic3'
# sim_subdirectory = 'with_applex'


E_photon = 1000 #eV

beam = generate_beam(E=14.0, dE=2.5e-3, I=7000, l_beam=5e-6, emit_n=0.5e-6, beta=20, l_window=6e-6, shape='flattop',chirp=0.0001) #shape='gaussian'
# beam = read_beam_file('/data/netapp/xfel/svitozar/projects/XFEL/parameters/beams/beam_0.1nC_sase1_12kev_fresh.txt')
beam.s -= np.amin(beam.s)
# tws_beam = Twiss(beam.pk())
# beam.add_chirp_poly(coeff=[0.0, +5., -1, -0, 0.01], s0=None)

sase_lat_pkg = create_fel_beamline('sase3')

rad_L = 2
inters_L = 1.2
rad_per = 0.09

n_u = rad_L // rad_per
n_is = inters_L // rad_per
# n_u = 332
# n_is = 80
rad_lat_pkg = create_fel_lattice(und_N = 35,
    und_L = rad_per*n_u,
    und_l = rad_per,
    und_K = 0.999,
    inters_L = rad_per*n_is,
    quad_L = rad_per*3,
    quad_K = 0,
    phs_L = rad_per,
    inters_K=0,
    inters_phi=0,
    quad_start = 'd',
    )
    
prepare_el_optics(beam, sase_lat_pkg, E_photon=E_photon, beta_av=20)
prepare_el_optics(deepcopy(beam), rad_lat_pkg, E_photon=E_photon*2, beta_av=20)

# tws_beam = Twiss(beam.pk())
# tws_lat =  twiss(sase_lat_pkg[0], tws_beam, nPoints=1000)
# plot_opt_func(sase_lat_pkg[0], tws_lat, legend=False)
# plt.show()
# plot_beam(beam, showfig=1, savefig=0, fig='Original beam')

sase_lat, _, _ = sase_lat_pkg
rad_lat, _, _ = rad_lat_pkg

# Create input container object
ginp = Genesis4Input()
ginp.filename = 'sase3.in'

# ginp.make_sequence(['setup', 'time', 'lattice', 'field', 'field#h3', 'beam#imported', 'sponrad', 'track', 'sort', 'alter_setup', 'write']) # defines an order of the the namelists

# Edit name lists
ginp.append_sequence('setup') # defines an order of the the namelists
ginp.sequence['setup'].rootname = 'sase3'
ginp.sequence['setup'].lattice = 'sase3.lat'
ginp.sequence['setup'].beamline = 'SASE3'
# ginp.sequence['setup'].gamma0 = 11357.82 # inserted using ginp.populate_sequence_beam()
ginp.sequence['setup'].lambda0 = eV2lambda(E_photon) # reference wavelength
ginp.sequence['setup'].delz = 0.1 # integration period [m]
ginp.sequence['setup'].seed = 123456789 # (ipseed in gen2)
ginp.sequence['setup'].npart = 8192
ginp.sequence['setup'].nbins = 4
ginp.sequence['setup'].one4one = 0
ginp.sequence['setup'].shotnoise = 1

ginp.append_sequence('time')
ginp.sequence['time'].slen = np.amax(beam.s) - np.amin(beam.s) # [m] length of time window
ginp.sequence['time'].sample = 10 # sampling in wavelengths ("zsep" in gen2)
ginp.sequence['time'].time = 1 #time-dependent

ginp.append_sequence('lattice')
#ginp.sequence['lattice'].zmatch = 20 # optional, already matched with OCELOT

ginp.append_sequence('field')
ginp.sequence['field'].power = 0e3 #small value to see both seed and noise
ginp.sequence['field'].dgrid = 2e-4
ginp.sequence['field'].ngrid = 101 # ("ncar" in gen2)
ginp.sequence['field'].waist_size = 30e-6
ginp.sequence['field'].harm = 1

##calculating and plotting field at 3rd harmonic
ginp.append_sequence('field#h3')
ginp.sequence['field#h3'].power = 0e3 #small value to see both seed and noise
ginp.sequence['field#h3'].dgrid = 2e-4
ginp.sequence['field#h3'].ngrid = 101 # ("ncar" in gen2)
ginp.sequence['field#h3'].waist_size = 30e-6
ginp.sequence['field#h3'].harm = 3

# ginp.append_sequence('beam')
# ginp.sequence['beam'].gamma = 27397.321010594387
# ginp.sequence['beam'].delgam = 4.8923787518918544
# ginp.sequence['beam'].current = 4997.704842822041
# ginp.sequence['beam'].ex = 5e-07
# ginp.sequence['beam'].ey = 5e-07
# ginp.sequence['beam'].betax = 24.748339004051143
# ginp.sequence['beam'].betay = 14.54507411359195
# ginp.sequence['beam'].alphax = 1.3348109600614686
# ginp.sequence['beam'].alphay = -0.800438812327674
# ginp.sequence['beam'].xcenter = 0.0
# ginp.sequence['beam'].ycenter = 0.0
# ginp.sequence['beam'].pxcenter = 0.0
# ginp.sequence['beam'].pycenter = 0.0
# ginp.sequence['beam'].bunch = 0
# ginp.sequence['beam'].bunchphase = 0
# ginp.sequence['beam'].emod = 0
# ginp.sequence['beam'].emodphase = 0

ginp.append_sequence('beam#imported')
ginp.populate_sequence_beam_array(beam_name_list_id='beam#imported', beam=beam)
# ginp.populate_sequence_beam(beam, 'beam#imported') # fills in all values from the Beam object (peak value)

ginp.append_sequence('sponrad')
ginp.sequence['sponrad'].seed = 123123 # for quantum fluctuation
ginp.sequence['sponrad'].doLoss = 1 # (isravg in gen2)
ginp.sequence['sponrad'].doSpread = 1 # (isrsig in gen2)

ginp.append_sequence('track')
ginp.sequence['track'].zstop = 40
ginp.sequence['track'].output_step = 3 # populate out file every nth integration step (iphsty in gen2)
ginp.sequence['track'].field_dump_step = 0
ginp.sequence['track'].beam_dump_step = 0
ginp.sequence['track'].sort_step = 0

ginp.append_sequence('sort')


# ginp.append_sequence('alter_setup')
# ginp.sequence['alter_setup'].beamline = 'APPLEX'
# ginp.sequence['alter_setup'].harmonic = 2

# ginp.append_sequence('track#2')
# ginp.sequence['track#2'].zstop = 7
# ginp.sequence['track#2'].output_step = 3 # populate out file every nth integration step (iphsty in gen2)
# ginp.sequence['track#2'].field_dump_step = 0
# ginp.sequence['track#2'].beam_dump_step = 0
# ginp.sequence['track#2'].sort_step = 0


ginp.append_sequence('write')
ginp.sequence['write'].field = ginp.sequence['setup'].rootname
ginp.sequence['write'].beam = ginp.sequence['setup'].rootname

# ginp.attachments.lat = {'SASE3': sase_lat, 'APPLEX': rad_lat}
ginp.attachments.lat = {'SASE3': sase_lat}

sase_sim = Genesis4Simulation(ginp,
                              exp_dir=os.path.join(sim_directory, sim_subdirectory),
                              return_out=1, # after the simulation reads and returns output file as Genesis4Output object
                              cleanup_afterwards=0, # deletes simulation files after the simulation (makes sence when return_out=1)
                              exec_script_path = __file__,
                              #zstop=ginp.sequence['track'].zstop
                              )

# Genesis4 simulation can be started simply using default parameters
sase_gout = sase_sim.run()

##launcher with simulation parameters can be specified:
# launcher=get_genesis4_launcher(launcher_program='genesis4', launcher_argument='', mpiParameters='-x PATH -x MPI_PYTHON_SITEARCH -x PYTHONPATH')
# launcher.dir = sase_sim.exp_dir
# launcher.argument = sase_sim.ginp.filename

##or fully customized via .command attribute, which will be executed via 'os.system(command)'
#launcher = MpiLauncher()
#launcher.command = "mkdir -p /gpfs/exfel/data/group/wp72/sserkez/projects/XFEL/2021_02_Genesis4_test/2023_02_Genesis4_MWE/with_harmonic3_m; cd /gpfs/exfel/data/group/wp72/sserkez/projects/XFEL/2021_02_Genesis4_test/2023_02_Genesis4_MWE/with_harmonic3_m; `which mpirun` -x PATH -x MPI_PYTHON_SITEARCH -x PYTHONPATH genesis4 sase3.in"
# sase_gout = sase_sim.run(launcher)

# If Genesis4 is not installed, one can at least write the files
# sase_sim.clean_output()
# sase_sim.write_all_files()

## if simulation runs successcfully, sase_gout object gets populates with an output from .out file. results can be plotted:
# plot_gen4_out_z(sase_gout, z=np.inf, showfig=0, savefig=1)
# plot_gen4_out_z(sase_gout, z=0, showfig=0, savefig=1)
# plot_gen4_out_e(sase_gout, showfig=0, savefig=1)
# plot_gen4_out_ph(sase_gout, showfig=0, savefig=1)

## Also one can specify just the common path to the files, or their file mask, e.g. missing .out.h5 or .fld.h5, and ocelot will scan through the files and attempt to plot everything.
plot_gen4_out_all(sase_gout.filePath.replace('out.h5','*'), showfig=0, savefig=1, choice='all') #plots reasonable ammount of output from .out, .fld and .par