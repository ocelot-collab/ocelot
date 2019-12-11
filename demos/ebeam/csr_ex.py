__author__ = 'Sergey Tomin'

from ocelot import *
from ocelot.gui.accelerator import *
from ocelot.adaptors import *
import time


# generate beam file
sigma_x = 0.000121407185261
sigma_px = 1.80989470506e-05
sigma_y = 0.000165584800564
sigma_py = 4.00994225888e-05

x = np.random.randn(200000)*sigma_x
px = np.random.randn(200000)*sigma_px
y = np.random.randn(200000)*sigma_y
py = np.random.randn(200000)*sigma_py

# covariance matrix for [tau, p] for beam compression in BC
cov_t_p = [[1.30190131e-06, 2.00819771e-05],
       [2.00819771e-05, 3.09815718e-04]]
long_dist = np.random.multivariate_normal((0, 0), cov_t_p, 200000)
tau = long_dist[:, 0]
dp = long_dist[:, 1]


p_array = ParticleArray(n=200000)
p_array.E = 0.130 # GeV
p_array.rparticles[0] = x
p_array.rparticles[1] = px
p_array.rparticles[2] = y
p_array.rparticles[3] = py
p_array.rparticles[4] = tau
p_array.rparticles[5] = dp

Q = 5e-9

p_array.q_array = np.ones(200000)*Q/200000

sI1, I1 = get_current(p_array, num_bins=200)

plt.figure(1)
plt.title("energy distribution: start")
plt.plot(-p_array.tau()*1000, p_array.p(), 'r.')
plt.xlabel("s, mm")
plt.ylabel("dp/p")
plt.grid(True)

plt.figure(2)
plt.title("current: start")
plt.plot(sI1*1000, I1)
plt.xlabel("s, mm")
plt.ylabel("I, A")
plt.grid(True)


# lattice

b1 = Bend(l = 0.501471, angle = 0.132729704703,  e1 = 0.0, e2 = 0.132729704703,   tilt = 0., eid= "b")
b2 = Bend(l = 0.501471, angle = -0.132729704703, e1 = -0.132729704703, e2 = 0.0,  tilt = 0., eid= "b")
b3 = Bend(l = 0.501471, angle = -0.132729704703, e1 = 0.0, e2 = -0.132729704703,  tilt = 0., eid= "b")
b4 = Bend(l = 0.501471, angle = 0.132729704703,  e1 = 0.132729704703, e2 = 0.0,   tilt = 0., eid= "b")
d1 = Drift(l=1.)
d2 = Drift(l=1.5)

cell = [Marker(eid="m1"), Drift(l=0.1), b1, d1, b2, d2, b3, d1, b4, d1, Marker(eid="m2")]

method = MethodTM()
method.global_method = SecondTM

lat = MagneticLattice(cell, method=method)


csr = CSR()
csr.traj_step = 0.0002
csr.apply_step = 0.0005



navi = Navigator(lat)
navi.add_physics_proc(csr, lat.sequence[0], lat.sequence[-1])

navi.unit_step = 0.05

start = time.time()
tws_track, p_array = track(lat, p_array, navi)
print("time execution = ",  time.time() - start)


plot_opt_func(lat, tws_track, top_plot=["E"], fig_name=0, legend=False)

sI1, I1 = get_current(p_array, num_bins=200)

plt.figure(3)
plt.title("energy distribution: end")
plt.plot(-p_array.tau()*1000, p_array.p(), 'r.')
plt.xlabel("s, mm")
plt.ylabel("dp/p")
plt.grid(True)

plt.figure(4)
plt.title("current: end")

plt.plot(sI1*1000, I1)
plt.xlabel("s, mm")
plt.ylabel("I, A")
plt.grid(True)
plt.show()

start = time.time()
show_e_beam(p_array)
print(time.time() - start, " sec")
plt.show()
show_e_beam(p_array, inverse_tau=True, headtail=False)
plt.show()