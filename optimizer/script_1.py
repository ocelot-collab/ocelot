from ocelot.optimizer.mint import mint
from ocelot.optimizer.mint import opt_objects as obj

d1 = obj.TestDevice(eid="d1")
d2 = obj.TestDevice(eid="d2")
d3 = obj.TestDevice(eid="d3")
target = obj.TestTarget()

minimizer = mint.Simplex()
minimizer.max_iter = 300


opt = mint.Optimizer()
opt.minimizer = minimizer
seq = [mint.Action(func=opt.max_target_func, args=[target, [d1, d2, d3]])]
opt.seq = seq

#opt.eval(seq)
opt.start()