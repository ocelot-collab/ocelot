#!/usr/bin/env python3

# Christoph Lechner, European XFEL, 2024-03-01
#
# Demo of upgraded create_fel_lattice function.
# . Matches a FODO cell containing two undulators. This is
#   possible thanks to new functionality in create_fel_lattice.
# . Finally, a GENESIS4 lattice file is written.

# import from Ocelot main modules and functions
import ocelot
from ocelot.adaptors.genesis4 import *
from ocelot.gui.beam_plot import *
from ocelot.utils.xfel_utils import create_fel_beamline

from copy import deepcopy

### Test that we always get the requested number of undulators ###
for Nundreq in range(2,38):
    test_lat_pkg = create_fel_beamline('sase1', override_und_N=Nundreq)
    test_lat, _, _ = test_lat_pkg
    ucntr=0
    for ele in test_lat.sequence:
        if isinstance(ele,Undulator):
            ucntr+=1
    print(f'{Nundreq} {ucntr}')
    if Nundreq!=ucntr:
	    raise ValueError('error: undulator counts do not match')

print('ok, we always get the requested number of undulators')

### Match demo ###
# Shorter lattice with SASE1-type FODO (added by CL to OCELOT in Feb 2024)
# New feature added in March-2024: Request complete final half-cell, allowing matching two-undulator-cell lattice
sase_lat_pkg = create_fel_beamline('sase1', override_und_N=2, final_halfcell_complete=True)
# sase_lat_pkg = create_fel_beamline('sase1')
sase_lat, _, _ = sase_lat_pkg
print(str(sase_lat))

betaavg = 20 # m
Ephot = 9000 # eV
beam = generate_beam(E=14, dE=0.003, I=5000, l_beam=1e-6, l_window=1.5e-6, shape='flattop', emit_n=0.4e-6, beta=betaavg, chirp=0)
prepare_el_optics(beam, sase_lat_pkg, E_photon=Ephot, beta_av=betaavg)

# drop lattice file
write_gen4_lat({'LM':sase_lat}, filepath='./lat_matched.lat')



### Draw plot ###
tws0 = Twiss()
tws0.beta_x = beam.beta_x[0]
tws0.beta_y = beam.beta_y[0]
tws0.alpha_x = beam.alpha_x[0]
tws0.alpha_y = beam.alpha_y[0]
print(tws0)

# import from Ocelot graphical modules
from ocelot.gui.accelerator import *

tws = twiss(sase_lat, tws0, nPoints=500) # nPoints gives beta functions more natual appearance (for instance: parabolic in drifts)
plot_opt_func(sase_lat, tws, legend=False)
plt.show()


