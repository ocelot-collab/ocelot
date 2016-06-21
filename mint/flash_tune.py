'''
main tuning script, LCLS
'''
import numpy as np

from ocelot.mint.mint import Optimizer, Action
from ocelot.mint.flash1_interface import FLASH1MachineInterface, FLASH1DeviceProperties, TestInterface




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


opt.eval(seq1)

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
