__author__ = 'Sergey Tomin'


from pylab import *
from ocelot import *
from ocelot.cpbd.optics import *

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


method = MethodTM(params={"global":SecondTM})
lat = MagneticLattice(cell, method=method)


t = linspace(0, 2*pi, num = 100)
x, xp = 0.1*cos(t), 0.1*sin(t)
plist = []
for xi, xpi in zip(x, xp):
    plist.append(Particle(x=xi, px=xpi))


plot(x, xp)
navi = Navigator()
dz = 10.
tracking_step(lat, plist, dz=dz, navi=navi)


x2 = [f.x for f in plist] #map(lambda f: f.x, plist)
xp2 = [f.px for f in plist] #map(lambda f: f.px, plist)
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
    if element.__class__ == Sextupole:
        element.k2 = 0.
lat.update_transfer_maps()

t = linspace(0, 2*pi, num = 100)
x, xp = 0.1*cos(t), 0.1*sin(t)
plist = []
for xi, xpi in zip(x,xp):
    plist.append(Particle(x = xi, px= xpi))

plot(x, xp)
navi = Navigator()
dz = 10.
tracking_step(lat, plist, dz = dz, navi = navi)

x2 = [f.x for f in plist] #map(lambda f: f.x, plist)
xp2 = [f.px for f in plist] #map(lambda f: f.px, plist)

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

