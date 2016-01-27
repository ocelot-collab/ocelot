'''
machine interface
includes online optimizer, response measurement and other stuff
'''
try:
    # in server "doocsdev12" set environment
    #  $ export PYTHONPATH=/home/ttflinac/user/python-2.7/Debian/
    import doocs
except:
    print 'error importing doocs library'
    
import re
from pylab import *



class FLASH1MachineInterface():
    def __init__(self):
        
        self.debug = False
        """
        self.blm_names = ['14L.SMATCH','14R.SMATCH',
                          '1L.UND1', '1R.UND1',
                          '1L.UND2', '1R.UND2', 
                          '1L.UND3', '1R.UND3', 
                          '1L.UND4', '1R.UND4',
                          '1L.UND5', '1R.UND5',
                          '1L.UND6', '1R.UND6',
                          '10SMATCH','3SDUMP']
        """
        self.blm_names = ['1L.UND1', '1R.UND1',
                          '1L.UND2', '1R.UND2',
                          '1L.UND3', '1R.UND3',
                          '1L.UND4', '1R.UND4',
                          '1L.UND5', '1R.UND5',
                          '1L.UND6', '1R.UND6',
                          '10SMATCH','3SDUMP']



    def init_corrector_vals(self, correctors):
        vals = np.zeros(len(correctors))#np.zeros(len(correctors))
        for i in range(len(correctors)):
            #print correctors[i],
            mag_channel = 'TTF2.MAGNETS/STEERER/' + correctors[i] + '/PS'
            vals[i] = doocs.read(mag_channel)
            #print vals[i], doocs.read(mag_channel), mag_channel
        return vals

    def get_cavity_info(self, cavs):
        ampls = [0.0]*len(cavs)#np.zeros(len(correctors))
        phases = [0.0]*len(cavs)#np.zeros(len(correctors))
        for i in range(len(cavs)):
            #ampl_channel = 'FLASH.RF/LLRF.CONTROLLER/CTRL.' + cavs[i] + '/SP.AMPL'
            #phase_channel = 'FLASH.RF/LLRF.CONTROLLER/CTRL.' + cavs[i] + '/SP.PHASE'
            ampl_channel = "FLASH.RF/LLRF.CONTROLLER/PVS." + cavs[i] + "/AMPL.SAMPLE"
            phase_channel = "FLASH.RF/LLRF.CONTROLLER/PVS." + cavs[i]+ "/PHASE.SAMPLE"
            #print(ampl_channel)
            #print(phase_channel)
            ampls[i] = doocs.read(ampl_channel)
            phases[i] = doocs.read(phase_channel)
            #print cavs[i], ampls[i], phases[i]
        return ampls, phases

    def get_gun_energy(self):
        gun_energy = doocs.read("FLASH.RF/LLRF.ENERGYGAIN.ML/GUN/ENERGYGAIN.FLASH1")
        gun_energy = gun_energy*0.001 # MeV -> GeV
        return gun_energy

    def get_bpms_xy(self, bpms):
        X = [0.0]*len(bpms)#np.zeros(len(correctors))
        Y = [0.0]*len(bpms)
        for i in range(len(bpms)):
            mag_channel = 'TTF2.DIAG/ORBIT/' + bpms[i]# + '/PS'
            #print mag_channel
            X[i] = doocs.read(mag_channel + "/X.FLASH1")*0.001 # mm -> m
            Y[i] = doocs.read(mag_channel + "/Y.FLASH1")*0.001 # mm -> m
            #print X, Y
        return X, Y

    def get_quads_current(self, quads):
        vals = [0.0]*len(quads)#np.zeros(len(correctors))
        for i in range(len(quads)):
            mag_channel = 'TTF2.MAGNETS/QUAD/' + quads[i]# + '/PS'
            vals[i] = doocs.read(mag_channel + "/PS")
        return vals

    def get_bends_current(self, bends):
        vals = [0.0]*len(bends)#np.zeros(len(correctors))
        for i in range(len(bends)):
            mag_channel = 'TTF2.MAGNETS/DIPOLE/' + bends[i]# + '/PS'
            vals[i] = doocs.read(mag_channel + "/PS")
        return vals

    def get_sext_current(self, sext):
        vals = [0.0]*len(sext)#np.zeros(len(correctors))
        for i in range(len(sext)):
            mag_channel = "TTF2.MAGNETS/SEXT/" + sext[i]
            vals[i] = doocs.read(mag_channel + "/PS")
        return vals

    def get_alarms(self):
        alarm_vals = np.zeros(len(self.blm_names))
        for i in xrange(len(self.blm_names)):
            blm_channel = 'TTF2.DIAG/BLM/'+self.blm_names[i]+'/CH00.TD'
            blm_alarm_ch  = ('TTF2.DIAG/BLM/'+self.blm_names[i]).replace('BLM', 'BLM.ALARM') + '/THRFHI'
            if (self.debug): print 'reading alarm channel', blm_alarm_ch
            alarm_val = doocs.read(blm_alarm_ch) * 1.25e-3 # alarm thr. in Volts
            if (self.debug): print 'alarm:', alarm_val
            sample = doocs.read(blm_channel)
            h = np.array([x[1] for x in sample])

            alarm_vals[i] = np.max( np.abs(h) ) / alarm_val 
            
        return alarm_vals

    def get_sase(self, detector='gmd_default'):
        
        if detector == 'mcp':
            # incorrect
            return doocs.read('TTF2.DIAG/MCP.HV/MCP.HV1/HV_CURRENT')
            #return np.abs( np.mean(h) )
        if detector == 'gmd_fl1_slow':
            return doocs.read('TTF2.FEL/BKR.FLASH.STATE/BKR.FLASH.STATE/SLOW.INTENSITY' )

        # default 'BKR' gmd
        h = np.array(doocs.read('TTF2.FEL/BKR.FLASH.STATE/BKR.FLASH.STATE/ENERGY.CLIP.SPECT'))
        val = np.mean(np.array([x[1] for x in h]))
        return val



    def get_sase_pos(self):

        x1 = doocs.read('TTF2.FEL/GMDPOSMON/TUNNEL/IX.POS')
        y1 = doocs.read('TTF2.FEL/GMDPOSMON/TUNNEL/IY.POS')

        x2 = doocs.read('TTF2.FEL/GMDPOSMON/BDA/IX.POS')
        y2 = doocs.read('TTF2.FEL/GMDPOSMON/BDA/IY.POS')
    
        return [ (x1,y1), (x2,y2) ] 

    def get_spectrum(self, f=None, detector='tunnel_default'):

        f_min = 13.0 # spectrum window (nm). TODO: replace with readout
        f_max = 14.0
        
        spec = np.array(doocs.read('TTF2.EXP/PBD.PHOTONWL.ML/WAVE_LENGTH/VAL.TD'))
    
        if f == None:
            f = np.linspace(f_min, f_max, len(spec))
    
        return f, spec
 
    def get_value(self, device_name):
        ch = 'TTF2.MAGNETS/STEERER/' + device_name + '/PS.RBV'
        return doocs.read(ch)
    
    def set_value(self, device_name, val):
        ch = 'TTF2.MAGNETS/STEERER/' + device_name + '/PS'
        return doocs.write(ch, str(val))
 
 
class FLASH1DeviceProperties:
    def __init__(self):
        self.patterns = {}
        self.limits = {}
        self.patterns['launch_steerer'] = re.compile('[HV][0-9]+SMATCH')
        self.limits['launch_steerer'] = [-4,4]
        
        self.patterns['intra_steerer'] = re.compile('H3UND[0-9]')
        self.limits['intra_steerer'] = [-5.0,-2.0]
        
        self.patterns['QF'] = re.compile('Q5UND1.3.5')
        self.limits['QF'] = [1,7]
        
        self.patterns['QD'] = re.compile('Q5UND2.4')
        self.limits['QD'] = [-5,-1]
        
        self.patterns['Q13MATCH'] = re.compile('Q13SMATCH')
        self.limits['Q13MATCH'] = [17.0,39.0]

        self.patterns['Q15MATCH'] = re.compile('Q15SMATCH')
        self.limits['Q15MATCH'] = [-16.0,-2.0]

        self.patterns['H3DBC3'] = re.compile('H3DBC3')
        self.limits['H3DBC3'] = [-0.035, -0.018]

        self.patterns['V3DBC3'] = re.compile('V3DBC3')
        self.limits['V3DBC3'] = [-0.1, 0.10]

        self.patterns['H10ACC7'] = re.compile('H10ACC7')
        self.limits['H10ACC7'] = [0.12, 0.17]

        self.patterns['H10ACC6'] = re.compile('H10ACC6')
        self.limits['H10ACC6'] = [-0.9, -0.4]

        self.patterns['H10ACC5'] = re.compile('H10ACC5')
        self.limits['H10ACC5'] = [0.8, 1.2]

        self.patterns['H10ACC4'] = re.compile('H10ACC4')
        self.limits['H10ACC4'] = [-0.2, 0.1]

        self.patterns['V10ACC7'] = re.compile('V10ACC7')
        self.limits['V10ACC7'] = [-2.6,-1.8]

        self.patterns['V10ACC4'] = re.compile('V10ACC4')
        self.limits['V10ACC4'] = [0.9,1.2]

        self.patterns['V10ACC5'] = re.compile('V10ACC5')
        self.limits['V10ACC5'] = [-0.7,-0.5]


        self.patterns['H8TCOL'] = re.compile('H8TCOL')
        self.limits['H8TCOL'] = [0.02,0.06]

        self.patterns['V8TCOL'] = re.compile('V8TCOL')
        self.limits['V8TCOL'] = [0.09,0.15]
        
        
    def get_limits(self, device):
        for k in self.patterns.keys():
            #print 'testing', k
            if self.patterns[k].match(device) != None:
                return self.limits[k]
        return [-2, 2]

    def get_polarity(self, quads):
        vals = [0.0]*len(quads)#np.zeros(len(correctors))
        for i in range(len(quads)):
            mag_channel = 'TTF2.MAGNETS/QUAD/' + quads[i]# + '/PS'
            vals[i] = doocs.read(mag_channel + "/PS.Polarity")
        return vals

    def get_type_magnet(self, quads):
        vals = [0.0]*len(quads)#np.zeros(len(correctors))
        for i in range(len(quads)):
            mag_channel = 'TTF2.MAGNETS/QUAD/' + quads[i]# + '/PS'
            vals[i] = doocs.get(mag_channel + "/DEVTYPE")
        return vals
