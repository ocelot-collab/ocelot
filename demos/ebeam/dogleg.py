"""
DogLeg. Second order achromat

S.Tomin. 10.2018
"""

from ocelot import *
from ocelot.gui.accelerator import *
import dogleg_lattice as dl
from ocelot.cpbd.beam import *


# create and plot dogleg lattice

method = MethodTM()
method.global_method = SecondTM

lat = MagneticLattice(dl.cell,  method=method)

tws = twiss(lat, tws0=dl.tws0)
plot_opt_func(lat, tws, top_plot=["Dy"])
plt.show()



def track_though_dl(emitt, energy_shift, tws0, energy=0.13, nparticles=100000):
    """
    Function generates the particleArray, matches the beam at the entrance of the dogleg and tracks it though lattice

    :param emitt: emittance for X and Y planes
    :param energy_shift: energy shift in dE/pc
    :param tws0: initial twiss parameters for matching
    :param energy: beam energy in GeV
    :param nparticles: number of particles
    :return: Twiss - Twiss parameters at the end of the lattice which was gotten from the Particle distribution
    """

    # generation of the ParticleArray
    p_array = generate_parray(sigma_x=np.sqrt(emitt*dl.tws0.beta_x), sigma_px=np.sqrt(emitt*dl.tws0.gamma_x),
                               energy=energy, nparticles=nparticles)

    # introduce energy shift
    p_array.p()[:] += energy_shift

    # Beam transform for the matching the particleArray at the entrance of the DL
    bt = BeamTransform(tws0)
    navi = Navigator(lat)
    navi.add_physics_proc(bt, lat.sequence[0], lat.sequence[0])

    # tracking
    tws_track, p_array = track(lat, p_array, navi=navi, calc_tws=True, print_progress=True)

    return tws_track[-1]


def phase_ellipse(alpha, beta, emitt):
    t = np.arange(0, 1, 0.01)*2*pi
    x = np.cos(t)
    y = np.sin(t)
    M = m_from_twiss([0, 1, 0], [alpha, beta, 0])
    x1 = (M[0, 0]*x + M[0, 1] * y)*np.sqrt(emitt)
    y1 = (M[1, 0]*x + M[1, 1] * y)*np.sqrt(emitt)
    return x1, y1



colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

fig, axs = plt.subplots(ncols=3, nrows=2)
fig.suptitle('Phase space at the end of Dogleg')
axs[0, 0].set(title="Sext ON. Bend in Y plane. X plane", xlabel=r"$x$, $\mu m$", ylabel=r"$p_x$, $\mu rad$")
axs[1, 0].set(title="Sext ON. Bend in Y plane. Y plane", xlabel=r"$y$, $\mu m$", ylabel=r"$p_y$, $\mu rad$")

for i in range(3):

    energy_shift = (-1 + i) * 0.05
    tws = track_though_dl(emitt=3.9308e-09, energy_shift=energy_shift, tws0=dl.tws0, energy=0.13, nparticles=100000)

    x, px = phase_ellipse(tws.alpha_x, tws.beta_x, tws.emit_x)
    y, py = phase_ellipse(tws.alpha_y, tws.beta_y, tws.emit_y)

    x += tws.x
    y += tws.y
    px += tws.px
    py += tws.py

    axs[0, 0].plot(x*1e6, px*1e6, linestyle="-", color=colors[i], label=r"$\frac{\Delta E}{pc}=$"+ str(energy_shift))
    axs[0, 0].plot(tws.x*1e6, tws.px*1e6, marker='o',color=colors[i])

    axs[1, 0].plot(y*1e6, py*1e6, linestyle="-", color=colors[i], label=r"$\frac{\Delta E}{pc}=$"+ str(energy_shift))
    axs[1, 0].plot(tws.y*1e6, tws.py*1e6, marker='o',color=colors[i])


for elem in lat.sequence:
    if elem.__class__ in [Sextupole, SBend]:
        elem.tilt = 0
    if elem.__class__ == Quadrupole:
        elem.k1 *= -1
lat.update_transfer_maps()

axs[0, 1].set(title="Sext ON. Bend in X plane. X plane", xlabel=r"$x$, $\mu m$", ylabel=r"$p_x$, $\mu rad$")
axs[1, 1].set(title="Sext ON. Bend in X plane. Y plane", xlabel=r"$y$, $\mu m$", ylabel=r"$p_y$, $\mu rad$")

for i in range(3):

    energy_shift = (-1 + i) * 0.05
    tws0 = Twiss()
    tws0.beta_x = dl.tws0.beta_y
    tws0.beta_y = dl.tws0.beta_x
    tws0.alpha_x =dl.tws0.alpha_y
    tws0.alpha_y =dl.tws0.alpha_x
    tws0.gamma_x =dl.tws0.gamma_y
    tws0.gamma_y =dl.tws0.gamma_x


    tws = track_though_dl(emitt=3.9308e-09, energy_shift=energy_shift, tws0=tws0, energy=0.13, nparticles=100000)

    x, px = phase_ellipse(tws.alpha_x, tws.beta_x, tws.emit_x)
    y, py = phase_ellipse(tws.alpha_y, tws.beta_y, tws.emit_y)

    x += tws.x
    y += tws.y
    px += tws.px
    py += tws.py

    axs[0, 1].plot(x*1e6, px*1e6, linestyle="-", color=colors[i], label=r"$\frac{\Delta E}{pc}=$"+ str(energy_shift))
    axs[0, 1].plot(tws.x*1e6, tws.px*1e6, marker='o',color=colors[i])

    axs[1, 1].plot(y*1e6, py*1e6, linestyle="-", color=colors[i], label=r"$\frac{\Delta E}{pc}=$"+ str(energy_shift))
    axs[1, 1].plot(tws.y*1e6, tws.py*1e6, marker='o',color=colors[i])

for elem in lat.sequence:
    if elem.__class__ == Sextupole:
        elem.k2 = 0
lat.update_transfer_maps()

axs[0, 2].set(title="Sext OFF. Bend in X plane. X plane", xlabel=r"$x$, $\mu m$", ylabel=r"$p_x$, $\mu rad$")
axs[1, 2].set(title="Sext OFF. Bend in X plane. Y plane", xlabel=r"$y$, $\mu m$", ylabel=r"$p_y$, $\mu rad$")

for i in range(3):

    energy_shift = (-1 + i) * 0.05
    tws0 = Twiss()
    tws0.beta_x = dl.tws0.beta_y
    tws0.beta_y = dl.tws0.beta_x
    tws0.alpha_x =dl.tws0.alpha_y
    tws0.alpha_y =dl.tws0.alpha_x
    tws0.gamma_x =dl.tws0.gamma_y
    tws0.gamma_y =dl.tws0.gamma_x


    tws = track_though_dl(emitt=3.9308e-09, energy_shift=energy_shift, tws0=tws0, energy=0.13, nparticles=100000)

    x, px = phase_ellipse(tws.alpha_x, tws.beta_x, tws.emit_x)
    y, py = phase_ellipse(tws.alpha_y, tws.beta_y, tws.emit_y)

    x += tws.x
    y += tws.y
    px += tws.px
    py += tws.py

    axs[0, 2].plot(x*1e6, px*1e6, linestyle="-", color=colors[i], label=r"$\frac{\Delta E}{pc}=$"+ str(energy_shift))
    axs[0, 2].plot(tws.x*1e6, tws.px*1e6, marker='o',color=colors[i])

    axs[1, 2].plot(y*1e6, py*1e6, linestyle="-", color=colors[i], label=r"$\frac{\Delta E}{pc}=$"+ str(energy_shift))
    axs[1, 2].plot(tws.y*1e6, tws.py*1e6, marker='o',color=colors[i])

plt.show()