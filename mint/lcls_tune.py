'''
main tuning script, LCLS
'''

from ocelot.mint.mint import Optimizer, Action, TestInterface
from lcls_interface import LCLSMachineInterface, LCLSDeviceProperties

mi = LCLSMachineInterface()
dp = LCLSDeviceProperties()
opt = Optimizer(TestInterface(), dp)
opt.debug = True

seq = [Action(func=opt.max_sase, args=[ ['YCOR:IN20:952','XCOR:IN20:951'] ] )]
opt.eval(seq)
#json.dumps(seq) # TODO: add json or other serialization to store mutated sequences
 
