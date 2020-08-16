from ocelot.cpbd.beam import BeamArray, generate_beam
from ocelot.adaptors.genesis import edist2beam, beam2edist
from ocelot.gui.beam_plot import plot_beam
from ocelot.gui.genesis_plot import plot_edist

"""
Created on Mon Jan  6 20:29:09 2020

@author: andrei
"""

# Create electron beam file (slice parameters) with polynomial energy chirp, quadratic beta, emittance, energy spread and non-zero alpha. Plot it.
beam = BeamArray()
beam = generate_beam(E=20.0, dE=2.5e-3, I=5000, l_beam=1e-6, emit_n=0.5e-6, beta=20, l_window=6e-6, shape='gaussian')
s0 = (np.max(beam.s) - np.min(beam.s))/ 2
beta_x0 = 15
beta_x0 = 15
emit_x0 = beam.emit_x[1]
emit_y0 = beam.emit_y[1]
dg0 = 6
beam.beta_x = beta_x0*(1 + 0.008*(beam.s - s0)**2/2 / (speed_of_light * 1e-15)**2 )
beam.beta_y = beta_x0*(1-0.4 + 0.006*(beam.s - s0)**2/2 / (speed_of_light * 1e-15)**2 )
beam.emit_x = emit_x0*(1+0.2 + 0.05*(beam.s - s0)**2/2 / (speed_of_light * 1e-15)**2 )
beam.emit_y = emit_y0*(1 + 0.01*(beam.s - s0)**2/2 / (speed_of_light * 1e-15)**2 )

beam.dg = dg0*(1 + 0.01*(beam.s - s0)**2/2 / (speed_of_light * 1e-15)**2 )

alpha_x = np.empty(np.shape(beam.s)[0])
alpha_x.fill(1)
alpha_y = np.empty(np.shape(beam.s)[0])
alpha_y.fill(-1)
beam.alpha_x = alpha_x
beam.alpha_y = alpha_y
beam.add_chirp_poly(coeff=[0.0, +5., -1, -0, 0.01], s0=None)
plot_beam(beam, showfig=1, savefig=0, fig='Original beam')

# Convert it to randomly distributed particle distribution and plot
edist = beam2edist(beam, npart=200000)
plot_edist(edist, fig_name='Generated distribution')

# recalculate slice parameters and generate another beam file
beam1 = edist2beam(edist)
plot_beam(beam1, showfig=1, savefig=0, fig='Re-calculated beam')

#%% Plot XXp and YYp distribution to check correlation introduced by non-zero alpha
bins=(150, 150, 150, 150)
cmin=0

plt.figure('X_Xp')
plt.hist2d(edist.x * 1e6, edist.xp * 1e6, [bins[0], bins[1]])
plt.xlabel('x [$\mu$m]')
plt.ylabel('xp [$\mu$rad]')
plt.show()

plt.figure('Y_Yp')
plt.hist2d(edist.y * 1e6, edist.yp * 1e6, [bins[0], bins[1]])
plt.xlabel('y [$\mu$m]')
plt.ylabel('yp [$\mu$rad]')
plt.show()