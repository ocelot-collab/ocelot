'''
main tuning script, LCLS
'''
import ocelot
from pylab import *
import numpy as np

from ocelot.utils.mint.mint import Optimizer, Action, TestInterface
from ocelot.utils.mint.flash1_interface_pydoocs3 import FLASH1MachineInterface, FLASH1DeviceProperties
from ocelot.utils.mint.machine_setup import *


class Optimization:
    def __init__(self, lattice, mi, dp):
        self.high_level_mint = HighLevelInterface(lattice, mi, dp)
        self.mi = mi
        self.dp = dp
        self.opt = Optimizer(mi, dp)

    def run(self, seq_dict, opt_params):
        self.set_limits(seq_dict)
        self.set_params(opt_params)
        seq = self.create_seq(seq_dict)
        self.opt.eval(seq)

    def set_params(self, opt_params):
        self.opt.debug = opt_params["debug"]             #True
        self.opt.logging = opt_params["logging"]         #True
        self.opt.log_file = opt_params['log file']       #'test.log'
        self.opt.timeout = opt_params['timeout']         #1.2

    def create_seq(self, seq_dict):

        sequence = []
        for act in seq_dict:
            func = self.opt.max_sase
            if act["func"] == "min_orbit":
                func = self.opt.min_orbit

            args = [act["devices"], act["method"], {'maxiter': act["maxiter"]}]
            print(args)
            action = Action(func=func, args=args)
            sequence.append(action)
        return sequence


    def read_dev_values(self, seq_dict):
        for act in seq_dict:
            for i, devname in enumerate(act["devices"]):
                self.high_level_mint.get_value(devname)


    def set_limits(self, seq_dict):
        limits = [0, 0]
        for act in seq_dict:
            for i, devname in enumerate(act["devices"]):
                tol = act["tol"][i]
                #print(tol,act['values'][i])
                lim = [float(s) for s in tol.split(',')]
                #print (lim,lim[0]/100.*act['values'][i])
                if len(lim) == 1:
                    limits[0] = float(act['values'][i])*(1. - lim[0]/100.)
                    limits[1] = float(act['values'][i])*(1. + lim[0]/100.)
                else:
                    limits[0] = lim[0]
                    limits[1] = lim[1]
                print(devname, limits)
                self.dp.set_limits( dev_name=devname, limits=limits)


"""

#import json
def get_dict(lat, bpms):
    dict_bpms = {}
    for elem in lat.sequence:
        if elem.type == "monitor" and elem.mi_id in bpms:
            dict_bpms[elem.mi_id] = {}
            dict_bpms[elem.mi_id]["x"] = elem.x
            dict_bpms[elem.mi_id]["y"] = elem.y
    return dict_bpms

#dp = FLASH1DeviceProperties()

mi = FLASH1MachineInterface()
dp = FLASH1DeviceProperties()
#opt = Optimizer(mi, dp)
opt = Optimizer(TestInterface(), dp)
opt.debug = True
opt.logging = True
opt.log_file = 'test.log'
opt.timeout = 1.2

seq1 = [Action(func=opt.max_sase, args=[ ['H10SMATCH','H12SMATCH'], 'simplex'] ) ]

seq2 = [Action(func=opt.max_sase, args=[ ['V14SMATCH','V7SMATCH'], 'simplex' ] )]
seq3 = [Action(func=opt.max_sase, args=[ ['V14SMATCH','V7SMATCH','H10SMATCH','H12SMATCH'], 'simplex' ] )]

seq4 = [Action(func=opt.max_sase, args=[ ['Q13SMATCH','Q15SMATCH'], 'simplex' ] )]
seq5 = [Action(func=opt.max_sase, args=[ ['H3DBC3','V3DBC3'], 'simplex' ] )]

seq6 = [Action(func=opt.max_sase, args=[ ['H3DBC3','V3DBC3','H10ACC7','V10ACC7'], 'simplex' ] )]

seq7 = [Action(func=opt.max_sase, args=[ ['Q5UND1.3.5','Q5UND2.4'], 'simplex' ] )]

seq8 = [Action(func=opt.max_sase, args=[ ['H3UND1','H3UND3','H3UND4','H3UND5'], 'simplex' ] )]

seq9 = [Action(func=opt.max_sase, args=[ ['H8TCOL','V8TCOL'], 'simplex' ] )]
seq10 = [Action(func=opt.max_sase, args=[ ['H3DBC3'], 'simplex' ] )]


seq0 = [Action(func=opt.max_sase, args=[ ['H10SMATCH','H12SMATCH'], 'cg', {'maxiter':15}] ), 
        Action(func=opt.max_sase, args=[ ['H10SMATCH','H12SMATCH'], 'simplex', {'maxiter':25}] )]


def apply_bump(names, currents, dIs, alpha):
        mi.set_value(names, currents+dIs*alpha)

cors = ['H3DBC3', 'H10ACC4','H9ACC5', 'H10ACC5', 'H9ACC6', 'H10ACC6', 'H10ACC7']
dI =  np.array([-0.0114768844711, -0.183727960466, 0.325959042831, 0.318743893708, 0.15280311903, 0.130996600233, -0.831909116508])
currents = np.array([ -0.0229914523661, 0.0250000003725, 0.985000014305, 0.0, -1.17299997807,  0.0, 0.148000001907])

bump = {"correctors":cors, "dI": dI, "currents":currents}
alpha = 0.1
seq_bump = [Action(func=opt.max_sase_bump, args=[ bump, alpha, 'simplex' ] )]


orbit = {}
orbit["correctors"] = ['H3SFELC', 'H4SFELC', 'H10SMATCH', 'D11SMATCH', 'H12SMATCH']

setup = log.MachineSetup()
#setup.save_lattice(lat, "init.txt")
lat_all = MagneticLattice(lattice)
setup.load_lattice("init.txt", lat_all)

orbit["bpms"] = get_dict(lat, bpms)




seq_min_orb = [Action(func=opt.min_orbit, args=[orbit, 'simplex' ] )]

opt.eval(seq_bump)
apply_bump(cors, currents, dI, alpha=0.1)
"""