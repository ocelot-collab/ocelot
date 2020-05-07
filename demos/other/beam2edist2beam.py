from ocelot.cpbd.beam import BeamArray, generate_beam
from ocelot.adaptors.genesis import edist2beam, beam2edist
from ocelot.gui.beam_plot import plot_beam
from ocelot.gui.genesis_plot import plot_edist

"""
Created on Mon Jan  6 20:29:09 2020

@author: andrei
"""

# Create electron beam file (slice parameters) with chirp and plot it
beam = BeamArray()
beam = generate_beam(E=20.0, dE=2.5e-3, I=5000, l_beam=1e-6, emit_n=0.5e-6, beta=20, l_window=6e-6, shape='gaussian')
beam.add_chirp_poly(coeff=[0.0, +5., -1, -0, 0.01], s0=None)
plot_beam(beam, showfig=1, savefig=0, fig='Original beam')

# Convert it to randomly distributed particle distribution and plot
edist = beam2edist(beam, npart=10000)
plot_edist(edist, fig_name='Generated distribution')

# recalculate slice parameters and generate another beam file
beam1 = edist2beam(edist)
plot_beam(beam1, showfig=1, savefig=0, fig='Re-calculated beam')