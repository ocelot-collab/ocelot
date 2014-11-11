__author__ = 'Sergey Tomin'

from ocelot.cpbd.match import *
#from ocelot.gui.accelerator import *
#import matplotlib.pyplot as plt
from pylab import *

Q1 = Quadrupole(l= 0.4, k1=-1.3, id = "Q1")
Q2 = Quadrupole(l= 0.8, k1=1.4, id = "Q2")
Q3 = Quadrupole(l= 0.4, k1=-1.7, id = "Q3")
Q4 = Quadrupole(l= 0.5, k1=1.19250444829 , id = "Q4")

B  = Bend(l=2.7, k1=-.06, angle=2*pi/16., e1=pi/16., e2=pi/16., id = "B")

SF = Sextupole(l=0., ms = 1.5, id = "SF") #random value
SD = Sextupole(l=0., ms = -1.5, id = "SD") #random value

D1 = Drift(l=2., id = "D1")
D2 = Drift(l=0.6, id = "D2")
D3 = Drift(l=0.3, id = "D3")
D4 = Drift(l=0.7, id = "D4")
D5 = Drift(l=0.9, id = "D5")
D6 = Drift(l=0.2, id = "D6")

cell = (D1, Q1, D2, Q2, D3, Q3, D4, B, D5, SD, D5, SF, D6, Q4, D6, SF, D5, SD,D5, B, D4, Q3, D3, Q2, D2, Q1, D1)



lat = MagneticLattice(cell)


t = linspace(0, 2*pi, num = 100)
x, xp = 0.1*cos(t), 0.1*sin(t)
plist = []
for xi, xpi in zip(x,xp):
    plist.append(Particle(x = xi, px= xpi))

plot(x, xp)
navi = Navigator(lat)
dz = 10.
track(lat, plist, dz = dz, navi = navi)


x2 = map(lambda f: f.x, plist)
xp2 = map(lambda f: f.px, plist)
suptitle("Tracking with sextupoles")
subplot(121)
plt.title("S = 0 m")
plt.plot(x, xp, "r.-", label = "X")
plt.xlabel("X, m")
plt.ylabel("Xp, rad")
plt.grid(True)

subplot(122)
plt.title("S = 10 m")
plt.plot(x2, xp2, "r.-", label = "X")
#plt.legend()
plt.xlabel("X, m")
plt.ylabel("Xp, rad")
plt.grid(True)
plt.show()


for element in lat.sequence:
    if element.type == "sextupole":
        element.ms = 0.
lat.update_transfer_maps()

t = linspace(0, 2*pi, num = 100)
x, xp = 0.1*cos(t), 0.1*sin(t)
plist = []
for xi, xpi in zip(x,xp):
    plist.append(Particle(x = xi, px= xpi))

plot(x, xp)
navi = Navigator(lat)
dz = 10.
track(lat, plist, dz = dz, navi = navi)


x2 = map(lambda f: f.x, plist)
xp2 = map(lambda f: f.px, plist)
plt.plot(x, xp, "r.-", label = "X")


suptitle("Tracking without sextupoles")
subplot(121)
plt.title("S = 0 m")
plt.plot(x, xp, "r.-", label = "X")
plt.xlabel("X, m")
plt.ylabel("Xp, rad")
plt.grid(True)

subplot(122)
plt.title("S = 10 m")
plt.plot(x2, xp2, "r.-", label = "X")
#plt.legend()
plt.xlabel("X, m")
plt.ylabel("Xp, rad")
plt.grid(True)
plt.show()

