__author__ = 'Sergey Tomin'




from ocelot.mint.mint import Optimizer
from ocelot.mint.xfel_interface import XFELMachineInterface, XFELDeviceProperties
#from flash1_interface import FLASH1DeviceProperties


##import json
def bpms_dict(bpms):
    dict_bpms = {}
    for bpm in bpms:
        dict_bpms[bpm] = {}
        dict_bpms[bpm]["x"] = 0.
        dict_bpms[bpm]["y"] = 0.
    return dict_bpms

def checking():
    print("@@@@@@@@@@@@   checking correctors ..... ")
    for name in horizantal + vertical:
        print(name)
        try:
            val = mi.get_value(name)
            print(val)
        except:
            print("not available")

    print("@@@@@@@@@@@@   checking BPM ..... ")
    for bpm in bpms:
        print(bpm)
        try:
            x, y = mi.get_bpms_xy([bpm])
            print(x, y)
        except:
            print("not available")

    vals = mi.init_corrector_vals(horizantal + vertical)
    print(vals)

    print("@@@@@@@@@@@@   checking BLM ..... ")
    alarms = mi.get_alarms()
    print("alarms", alarms)


def set_limits(corrs, limit):
    for name in corrs:
        val = mi.get_value(name)
        dp.set_limits(name, [val - limit, val + limit])
        print(name," I = ", val, "A; ",  "limits:", dp.get_limits(name))

#dp = FLASH1DeviceProperties()

mi = XFELMachineInterface()
dp = XFELDeviceProperties()
#opt = Optimizer(mi, dp)
#mi = TestInterface()
opt = Optimizer(mi, dp)
opt.debug = True
opt.logging = True
opt.log_file = 'test.log'
opt.timeout = 1.2


orbit = {}

horizantal = [
                'CKX.23.I1',
                'CKX.24.I1',
                'CKX.25.I1',
                'CX.37.I1',
                'CX.39.I1',
                'CIX.51.I1',
                'CIX.57.I1',
                'CIX.65.I1',
                #'CIX.73I.I1',
                #'CIX.73II.I1',
                #'CIX.76.I1',
                #'CIX.78.I1',
                #'CIX.81.I1',
                #'CIX.83.I1',
               ]

vertical = ['CKY.23.I1',
            'CKY.24.I1',
            'CKY.25.I1',
            'CY.37.I1',
            'CY.39.I1',
            'CIY.51.I1',
            'CIY.55.I1',
            'CIY.58.I1',
            'CIY.63.I1',
            #'CIY.72.I1',
            ]

bpms = [
        'BPMG.24.I1',
        'BPMG.25I.I1',
        'BPMC.38I.I1',
        'BPMR.38II.I1',
        'BPMF.47.I1',
        'BPMF.48.I1',
        'BPMF.52.I1',
        'BPMA.55.I1',
        'BPMA.57.I1',
        'BPMA.59.I1',
        #'BPMATEST.60.I1',
        #'BPMATEST.61.I1',
        'BPMA.63.I1',
        #'BPMA.72.I1',
        #'BPMA.75.I1',
        #'BPMA.77.I1',
        #'BPMA.80.I1',
        #'BPMA.82.I1',
        #'BPMA.85.I1',
        ]



orbit["correctors"]  =  horizantal + vertical #['V3DBC3', 'V10ACC4', 'H10ACC5', 'H10ACC6', 'H10ACC7']

orbit["bpms"] = bpms_dict(bpms)#bpms

checking()

mi.set_value('CKX.24.I1', 0.228)
#set_limits(horizantal, limit= 0.01)
#set_limits(vertical, limit= 0.01)

#seq_min_orb = [Action(func=opt.min_orbit, args=[orbit, 'simplex' ] )]

#opt.eval(seq_min_orb)


#opt.eval(seq5 + seq3 + seq6 + seq8 + seq9)
