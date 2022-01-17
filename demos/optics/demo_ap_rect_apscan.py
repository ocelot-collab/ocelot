#!/usr/bin/env python3.6

# C. Lechner, 2021-10-01
# demo for scan of energy vs. rectangular aperture

import math
import matplotlib.pyplot as plt
from ocelot.optics.wave import *
from ocelot import ocelog
plt.ion()
ocelog.setLevel(logging.DEBUG)

# generate a beam
lambda0=1e-9
field_at_waist = generate_gaussian_dfl(xlamds=lambda0, shape=(101,101,100), dgrid=(2e-3,2e-3,50e-6))

# run aperture scan
en,cutoffs=dfl_ap_rect_apscan(field_at_waist)

# compare with OCELOT function...
dfl=field_at_waist
stuff=deepcopy(dfl)
dfl_ap_rect(stuff, ap_x=2*cutoffs[2], ap_y=2*cutoffs[2]) # function modifies field!
print('Comparing relative energy reduction by cut (values should be in agreement):')
print('       OCELOT function dfl_ap_rect: {}'.format(stuff.E()/dfl.E()))
print('OCELOT function dfl_ap_rect_apscan: {}'.format(en[2]/dfl.E()))


### plot the result
fig,axs = plt.subplots(2,1)
axs[0].plot(cutoffs,np.array(en)/dfl.E())
axs[0].set_ylabel('rel. energy transmission E/E0')
axs[1].semilogy(cutoffs,1-np.array(en)/dfl.E())
axs[1].set_ylabel('rel. energy reduction, 1-E/E0')
axs[1].set_xlabel('half-size of quadratic aperture [m]')
axs[1].set_ylim([1e-6,1])
plt.show()
