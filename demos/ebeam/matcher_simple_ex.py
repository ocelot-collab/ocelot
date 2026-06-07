from ocelot import *
from ocelot.gui import *
from ocelot.cpbd.matcher import MatchProblem


start = Marker(eid="start")
d = Drift(l=0.5, eid="d")
end = Marker(eid="end")
q1 = Quadrupole(l=0.2, k1=0.5, eid="q1")
q2 = Quadrupole(l=0.2, k1=-0.5, eid="q2")
cell = [start, d, q1, d, q2, d, end]

lat = MagneticLattice(cell)

tw0 = Twiss()
tw0.beta_x = 9.0
tw0.beta_y = 10.0
tw0.E = 1.0

problem = MatchProblem(lat, tw0, periodic=False)

# Variables
problem.vary_element(q1, quantity="k1", limits=(-5, 5))
problem.vary_element(q2, quantity="k1", limits=(-5, 5))
problem.vary_twiss("alpha_x", limits=(-5, 5))
problem.vary_twiss("alpha_y", limits=(-5, 5))

# Twiss targets
problem.target_twiss(end, "alpha_x", 0.0, weight=1e6)
problem.target_twiss(end, "beta_y", 9.0, weight=1e6)

# Partial periodic targets
problem.target_periodic_twiss("beta_x", start, end, weight=1e6)
problem.target_periodic_twiss("alpha_y", start, end, weight=1e6)

result = problem.solve(solver="ls_trf", max_iter=300)
print(result.success, result.merit)

print("solved twiss0:")
print("  beta_x =", problem.twiss0.beta_x)
print("  alpha_x =", problem.twiss0.alpha_x)
print("  beta_y =", problem.twiss0.beta_y)
print("  alpha_y =", problem.twiss0.alpha_y)

merit, target_reports, objective_reports, state = problem.evaluate()
print("beta_x start/end =", state.twiss_at(start).beta_x, state.twiss_at(end).beta_x)
print("alpha_y start/end =", state.twiss_at(start).alpha_y, state.twiss_at(end).alpha_y)

for report in result.target_reports:
    print(report.name, report.details)

print(result.variables)
tws = twiss(lat, problem.twiss0)
plot_opt_func(lat,tws)
plt.show()
