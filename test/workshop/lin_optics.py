__author__ = 'Sergey Tomin'
"""
unofficial linear optics of the XFEL injector
"""
from ocelot import *
from ocelot.gui.accelerator import *
from ocelot.test.workshop.injector_sim import *

tws0 = Twiss()
tws0.beta_x = 29.171
tws0.beta_y = 29.171
tws0.alpha_x = 10.955
tws0.alpha_y = 10.955
tws0.E = 0.005


method = MethodTM()
method.global_method = TransferMap

lat = MagneticLattice(cell, method=method)

tws = twiss(lat, tws0, nPoints=None)

plot_opt_func(lat, tws, top_plot=["Dx"], fig_name="i1", legend=False)
plt.show()