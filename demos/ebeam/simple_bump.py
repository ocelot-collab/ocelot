from ocelot import *
from ocelot.utils.bump_utils import *
from ocelot.gui import *

d = Drift(l=0.35)
d1 = Drift(l=0.6)
qf = Quadrupole(l=0.2, k1=4)
qd = Quadrupole(l=0.2, k1=-4)

c1 = Vcor(l=0.1)
c2 = Vcor(l=0.1)
c3 = Vcor(l=0.1)
c4 = Vcor(l=0.1)

m = Marker()

cell = (d, qf, c1, d1, qd, c2, d1, m, qf, c3, d1, qd, c4, d1, qf, d)

lat = MagneticLattice(cell)
tws0 = Twiss()
tws0.beta_x = 10
tws0.beta_y = 10

tws = twiss(lat, tws0)
plot_opt_func(lat, tws)
plt.show()
cor_list = [c1, c2, c3, c4]
a = bump_4cors(lat, cor_list, marker=m, x=0.001, xp=-0.00, energy=1)
print("corrector, strength: ", a * 1000, " mrad")

plist = lattice_track(lat, Particle(E=14))

s = np.array([t.s for t in plist])
y = np.array([t.y for t in plist])
x = np.array([t.x for t in plist])
fig, ax = plot_API(lat, legend=False, font_size=15, fig_name="Bump")
ax.plot(s, y * 1000, lw=3, label="Trajectory: Y")
ax.plot(s, x * 1000, lw=3, label="Trajectory: X")
ax.set_ylabel("Y [mm]")
ax.legend()
plt.show()

lat_new = convert_cors2dipoles(lat, cor_list, energy=14)

tws = twiss(lat_new, tws0=tws0)
plot_opt_func(lat_new, tws, top_plot=["Dx", "Dy"], fig_name="1")
plt.show()