from ocelot import *
from ocelot.gui.accelerator import *
from ocelot.test.xfel_lat.I1B1 import *
from ocelot.test.xfel_lat.B1D import *


# undu_49_i1.lperiod = 0.074
# undu_49_i1.nperiods = 10
# undu_49_i1.Kx = 1.935248
# undu_49_i1.l = 0.74
# Undulator(lperiod=0.074, nperiods=10, Kx=1.935248, Ky=0.0, eid='UNDU.49.I1')
lat = MagneticLattice(cell)

tws0 = Twiss()
tws0.E = 0.005000000
tws0.beta_x = 29.171000000
tws0.alpha_x = 10.955000000
tws0.beta_y = 29.171000
tws0.alpha_y = 10.955000

tws = twiss(lat, tws0)

plot_opt_func(lat, tws, top_plot=["Dx", "Dy"])
plt.show()