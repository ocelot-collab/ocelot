__author__ = 'Sergey Tomin'

import matplotlib.pyplot as plt
from ocelot import *
from ocelot.cpbd.optics import *
from ocelot.cpbd.transformations.second_order import SecondTM

Q1 = Quadrupole(l= 0.4, k1=-1.3, eid= "Q1")
Q2 = Quadrupole(l= 0.8, k1=1.4, eid= "Q2")
Q3 = Quadrupole(l= 0.4, k1=-1.7, eid= "Q3")
Q4 = Quadrupole(l= 0.5, k1=1.19250444829 , eid= "Q4")

B  = Bend(l=2.7, k1=-.06, angle=2*pi/16., e1=pi/16., e2=pi/16., eid= "B")

SF = Sextupole(l=0.01, k2 = 150, eid= "SF") #random value
SD = Sextupole(l=0.01, k2 = -150, eid= "SD") #random value

D1 = Drift(l=2., eid= "D1")
D2 = Drift(l=0.6, eid= "D2")
D3 = Drift(l=0.3, eid= "D3")
D4 = Drift(l=0.7, eid= "D4")
D5 = Drift(l=0.9, eid= "D5")
D6 = Drift(l=0.2, eid= "D6")

cell = (D1, Q1, D2, Q2, D3, Q3, D4, B, D5, SD, D5, SF, D6, Q4, D6, SF, D5, SD,D5, B, D4, Q3, D3, Q2, D2, Q1, D1)


method = {'global': SecondTM}
lat = MagneticLattice(cell, method=method)


def track_ellipse(lat, dz):
    t = np.linspace(0, 2*pi, num=100)
    x, xp = 0.1*np.cos(t), 0.1*np.sin(t)
    plist = [Particle(x=xi, px=xpi) for xi, xpi in zip(x, xp)]

    navi = Navigator(lat)
    tracking_step(lat, plist, dz=dz, navi=navi)

    x2 = [p.x for p in plist]
    xp2 = [p.px for p in plist]
    return x, xp, x2, xp2


def plot_ellipse(title, x, xp, x2, xp2, dz):
    fig, (ax0, ax1) = plt.subplots(1, 2, num=title, figsize=(10, 4), clear=True)
    fig.suptitle(title)

    ax0.set_title("S = 0 m")
    ax0.plot(x, xp, "r.-", label="X")
    ax0.set_xlabel("X, m")
    ax0.set_ylabel("Xp, rad")
    ax0.grid(True)

    ax1.set_title(f"S = {dz:g} m")
    ax1.plot(x2, xp2, "r.-", label="X")
    ax1.set_xlabel("X, m")
    ax1.set_ylabel("Xp, rad")
    ax1.grid(True)

    fig.tight_layout()


dz = 10.
x, xp, x2, xp2 = track_ellipse(lat, dz)
plot_ellipse("Tracking with sextupoles", x, xp, x2, xp2, dz)
plt.show()


for element in lat.sequence:
    if element.__class__ == Sextupole:
        element.k2 = 0.

x, xp, x2, xp2 = track_ellipse(lat, dz)
plot_ellipse("Tracking without sextupoles", x, xp, x2, xp2, dz)
plt.show()
