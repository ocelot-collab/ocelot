'''
main tuning script, LCLS
'''

from pylab import *
import numpy as np

from ocelot.utils.mint.mint import Optimizer, Action, TestInterface
from flash1_interface import FLASH1MachineInterface, FLASH1DeviceProperties
#from flash1_interface import FLASH1DeviceProperties


#import json


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

opt.eval(seq0)

#opt.eval(seq5 + seq3 + seq6 + seq8 + seq9)
 

