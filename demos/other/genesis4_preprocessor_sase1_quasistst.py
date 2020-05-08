import os
import functools

import ocelot
from ocelot.gui.beam_plot import *
from ocelot.adaptors.genesis4 import *
from ocelot.gui.genesis4_plot import *

_logger = logging.getLogger(__name__)
ocelog.setLevel(logging.INFO)

sim_directory = r'/gpfs/exfel/data/scratch/svitozar/projects/ocelot_test/genesis4/Ocelot/preprocessor_local/'

E_photon = 7000 #eV

beam = generate_beam(E=14.0, dE=2.5e-3, I=5000, l_beam=1e-6, emit_n=0.5e-6, beta=20, l_window=6e-6, shape='gaussian')
tws_beam = Twiss(beam.pk())
# beam.add_chirp_poly(coeff=[0.0, +5., -1, -0, 0.01], s0=None)

sase_lat_pkg = create_fel_beamline('sase1')
prepare_el_optics(beam, sase_lat_pkg, E_photon=E_photon, beta_av=20)

# tws_beam = Twiss(beam.pk())
# tws_lat =  twiss(sase_lat_pkg[0], tws_beam, nPoints=1000)
# plot_opt_func(sase_lat_pkg[0], tws_lat, legend=False)
# plt.show()
# plot_beam(beam, showfig=1, savefig=0, fig='Original beam')

sase_lat, _, _ = sase_lat_pkg

# Create input container object
ginp = Genesis4Input()
ginp.filename = 'sase1.in'
ginp.make_sequence(['setup', 'time', 'lattice', 'field', 'beam#imported', 'sponrad', 'track', 'write']) # defines an order of the the namelists

# Edit name lists
ginp.sequence['setup'].rootname = 'sase1'
ginp.sequence['setup'].lattice = 'sase1.lat'
ginp.sequence['setup'].beamline = 'SASE1'
# ginp.sequence['setup'].gamma0 = 11357.82 # inserted using ginp.populate_sequence_beam()
ginp.sequence['setup'].lambda0 = eV2lambda(E_photon) # reference wavelength
ginp.sequence['setup'].delz = 0.04 # integration period [m]
ginp.sequence['setup'].seed = 123456789 # (ipseed in gen2)
ginp.sequence['setup'].npart = 8192
ginp.sequence['setup'].nbins = 4
ginp.sequence['setup'].one4one = 0
ginp.sequence['setup'].shotnoise = 1

#ginp.sequence['lattice'].zmatch = 20 # optional, already matched with OCELOT

ginp.sequence['time'].slen = 5e-7 # [m] length of time window
ginp.sequence['time'].sample = 30 # sampling in wavelengths ("zsep" in gen2)
ginp.sequence['time'].time = 1 #time-dependent

ginp.sequence['field'].power = 1e6 #small value to see some seed
ginp.sequence['field'].dgrid = 1e-04
ginp.sequence['field'].ngrid = 151 # ("ncar" in gen2)
ginp.sequence['field'].waist_size = 30e-6

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
ginp.populate_sequence_beam(beam, 'beam#imported') # fills in all values from the Beam object

ginp.sequence['sponrad'].seed = 123123 # for quantum fluctuation
ginp.sequence['sponrad'].doLoss = 1 # (isravg in gen2)
ginp.sequence['sponrad'].doSpread = 1 # (isrsig in gen2)

ginp.sequence['write'].field = ginp.sequence['setup'].rootname
ginp.sequence['write'].beam = ginp.sequence['setup'].rootname

ginp.sequence['track'].zstop = 75
ginp.sequence['track'].output_step = 3 # populate out file every nth integration step (iphsty in gen2)
ginp.sequence['track'].field_dump_step = 0
ginp.sequence['track'].beam_dump_step = 0
ginp.sequence['track'].sort_step = 0

ginp.attachments.lat = {'SASE1': sase_lat}
sase_sim = Genesis4Simulation(ginp,
                              exp_dir=sim_directory + 'Example-SASE1_1',
                              return_out=1, # after the simulation reads and returns output file as Genesis4Output object
                              cleanup_afterwards=0, # deletes simulation files after the simulation (makes sence when return_out=1)
                              zstop=50)

# Genesis4 simulation is usually started simply using
sase_gout = sase_sim.run()

# But as it is not installed here, let's just write the files
# sase_sim.clean_output()
# sase_sim.write_all_files()

# plot_gen4_out_z(sase_gout, z=np.inf, showfig=0, savefig=1)
# plot_gen4_out_z(sase_gout, z=0, showfig=0, savefig=1)
# plot_gen4_out_e(sase_gout, showfig=0, savefig=1)
# plot_gen4_out_ph(sase_gout, showfig=0, savefig=1)
plot_gen4_out_all(sase_gout, showfig=0, savefig=1)