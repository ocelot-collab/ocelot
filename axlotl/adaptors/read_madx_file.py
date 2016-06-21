from ocelot import MagneticLattice

__author__ = 'Sergey Tomin'

from ocelot.adaptors.madx import lattice_str_from_madx, madx_seq2ocelot_seq
from ocelot.cpbd import io
from ocelot.gui.accelerator import *



# the first method - recommended
# MAD-X can generate sequence of elements and the sequence can be used for Ocelot lattice construction


def RFcavity(l, volt, lag, harmon):
    rf = Cavity(l = l, eid= id)
    rf.volt = volt
    rf.lag = lag
    rf.harmon = harmon
    return rf

lines_seq = lattice_str_from_madx("data/p3x_v16.seq")
exec("\n".join(lines_seq))
seq = madx_seq2ocelot_seq(lattice, tot_length = ring.l, exclude_elems = [])

lat = MagneticLattice(seq)
print "ring circumstance = ", lat.totalLen
tws = twiss(lat, Twiss())
plot_opt_func(lat, tws)
plt.show()

#save lattice to file
io.write_lattice(lat, "data/petra.inp")


# second method more difficult and probability of errors is higher

lines_data = lattice_str_from_madx("data/quadsex_p3x_v16.dat")
lines_geo = lattice_str_from_madx("data/petra3_upv16.geo")

#save_lattice_str(lines_data, "quadsex_p3x_v16.inp")
#save_lattice_str(lines_geo, "petra3_upv16.inp")
exec("\n".join(lines_data))
exec("\n".join(lines_geo))

#lines_seq = lattice_str_from_madx("data/p3x_v16.seq")

seq = madx_seq2ocelot_seq(lattice, tot_length = ring.l, exclude_elems = [])

lat = MagneticLattice(seq)
print "ring circumstance = ", lat.totalLen
tws = twiss(lat, Twiss())
plot_opt_func(lat, tws)
plt.show()


