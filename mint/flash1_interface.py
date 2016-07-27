'''
machine interface
includes online optimizer, response measurement and other stuff
'''
try:
    # in server "doocsdev12" set environment
    #  $ export PYTHONPATH=/home/ttflinac/user/python-2.7/Debian/
    import pydoocs
except:
    print ('error importing doocs library')
"""
import sys
sys.path.append("/home/ttflinac/user/tomins/repository/pyqtgraph-0.9.10/build/lib")
sys.path.append("/home/ttflinac/user/tomins/repository")
"""
import re
from pylab import *
import time
import pickle
from threading import Thread, Lock
from ocelot.utils.db import *
from ocelot.common.globals import *

class SaveOptParams:
    def __init__(self, mi, dp, lat=None, dbname='flash.db'):
        self.mi = mi
        self.dp = dp
        self.dbname = dbname
        self.db = PerfDB(dbname = self.dbname)
        self.data = []
        self.wasSaved = False
        if lat != None:
            self.lat = lat

    def read_magnets(self, types):
        id2I_dict = {}
        for elem in self.lat.sequence:
            if elem.type in types:
                try:
                    id2I_dict[elem.id] = self.mi.get_value(elem.id)
                except:
                    id2I_dict[elem.id] = None
        return id2I_dict

    def read_bpms(self):
        orbit = []
        L = 0.
        for elem in self.lat.sequence:
            L += elem.l
            if elem.type == "monitor":
                try:
                    x, y = self.mi.get_bpms_xy([elem.id])
                except:
                    x, y = [None], [None]
                orbit.append([elem.id, L - elem.l/2., x[0], y[0]])
        return orbit

    def read_cavs(self):
        dict_cavity = {}
        for elem in self.lat.sequence:
            if elem.type == "cavity":
                try:
                    ampl = self.mi.get_cav_ampl(elem.id)
                    phi = self.mi.get_cav_phase(elem.id)
                except:
                    ampl, phi = None, None
                dict_cavity[elem.id + ".ampl"] = ampl
                dict_cavity[elem.id + ".phi"] = phi
        return dict_cavity

    def send_to_db(self):
        self.db = PerfDB(dbname = self.dbname)
        tune_id = self.db.current_tuning_id()
        print ('new action for tune_id', tune_id)
        self.db.new_action(tune_id, start_sase = self.data[0]["sase_slow"], end_sase = self.data[1]["sase_slow"])
        #print ('current actions in tuning', [(t.id, t.tuning_id, t.sase_start, t.sase_end) for t in self.db.get_actions()])
        action_id = self.db.current_action_id()
        print ('updating', tune_id, action_id)
        names = []
        start_data = []
        stop_data = []
        for name in self.data[0].keys():
            names.append(name)
            start_data.append(self.data[0][name])
            stop_data.append(self.data[1][name])
            print('{0}, {1}, {2}'.format(name, self.data[0][name], self.data[1][name]) )
        self.db.add_action_parameters(tune_id, action_id, param_names=names, start_vals=start_data, end_vals=stop_data)
        """
        names = np.array([self.data[0]["devices"]])
        names = np.append(names, ["tol." + name for name in self.data[0]["devices"]])

        sase_pos_name = ["sase_pos.tun.x", "sase_pos.tun.y", "sase_pos.bda.x","sase_pos.bda.y"]
        names = np.append(names, sase_pos_name)

        names = np.append(names, "timestamp")
        names = np.append(names, "method")
        names = np.append(names, "maxiter")
        names = np.append(names, "sase")
        names = np.append(names, "sase_slow")
        names = np.append(names, "flag")
        names = np.append(names, "niter")
        names = np.append(names, "pBPM.tun")
        names = np.append(names, "pBPM.bda")
        names = np.append(names, "charge")
        names = np.append(names, "wavelength")
        print("NAMES ", names, self.data[0].keys(), self.data[0]['devices'])

        names2 = []
        for name in self.data[0].keys():
            if self.data[0][name].__class__ == dict:
                for name2 in self.data[0][name].keys():
                    print("inside", name2)
            else:
                print(name)

        vals = [np.array([]), np.array([])]
        for i, val in enumerate(vals):
            #print(val,  self.data[i]["currents"])
            vals[i] = np.append(vals[i], self.data[i]["currents"])
            vals[i] = np.append(vals[i], [(lim[1]-lim[0])/2. for lim in self.data[i]["limits"]])
            pos = self.data[i]["sase_pos"]
            vals[i] = np.append(vals[i], [pos[0][0],  pos[0][1],  pos[1][0], pos[1][1]])

            vals[i] = np.append(vals[i], self.data[i]["timestamp"])
            vals[i] = np.append(vals[i], self.data[i]["method"])
            vals[i] = np.append(vals[i], self.data[i]["maxiter"])
            vals[i] = np.append(vals[i], self.data[i]["sase"])
            vals[i] = np.append(vals[i], self.data[i]["sase_slow"])
            vals[i] = np.append(vals[i], self.data[i]["flag"])
            vals[i] = np.append(vals[i], self.data[i]["niter"])
            vals[i] = np.append(vals[i], self.data[i]["pBPM"]["Tunnel"])
            vals[i] = np.append(vals[i], self.data[i]["pBPM"]["BDA"])
            vals[i] = np.append(vals[i], self.data[i]["charge"])
            vals[i] = np.append(vals[i], self.data[i]["wavelength"])
        for name, val in zip(names, vals[0]):
            print("db ", name, val)
        #orbit !!!
        """
        """
        orbit1 = self.data[0]["orbit"]
        orbit2 = self.data[1]["orbit"]
        for i, bpm in enumerate(orbit1):
            names = np.append(names, ["bpm.x." + bpm[0]])
            names = np.append(names, ["bpm.y." + bpm[0]])
            vals[0] = np.append(vals[0], bpm[2])
            vals[0] = np.append(vals[0], bpm[3])
            bpm2 = orbit2[i]
            vals[1] = np.append(vals[1], bpm2[2])
            vals[1] = np.append(vals[1], bpm2[3])
        """
        #print("shape datda", len(self.data) )
        #print(self.data[0]["sase_pos"])
        #print(names, vals[0], vals[1])
        #print(len(names), len(vals[0]), len(vals[1]))

        #self.db.add_action_parameters(tune_id, action_id, param_names = names, start_vals = vals[0], end_vals=vals[1])
        return True


    def save(self, args, time, niter, flag="start"):
        data_base = {}
        data_base["flag"] = flag
        data_base["timestamp"] = time
        #data_base["devices"] = args[0]
        data_base["method"] = args[1]
        data_base["maxiter"] = args[2]["maxiter"]
        #print('data_base["pBPM"] = args[2]["pBPM"]', args[2]["pBPM"])
        #data_base["pBPM"] = args[2]["pBPM"]

        for dev in args[0]:
            data_base[dev] = self.mi.get_value(dev)
            lim = self.dp.get_limits(dev)
            data_base["lim."+dev] = (lim[1] - lim[0])/2.

        pos_names = ["sase_pos.tun.x", "sase_pos.tun.y", "sase_pos.bda.x", "sase_pos.bda.y"]
        pos_sase = np.ravel(self.mi.get_sase_pos())
        for name, sase_pos in zip(pos_names, pos_sase):
            data_base[name] = sase_pos

        #data_base["sase_pos"] = self.mi.get_sase_pos()
        data_base["niter"] = niter
        data_base["sase"] = self.mi.get_sase()
        data_base["sase_slow"] = self.mi.get_sase(detector='gmd_fl1_slow')
        data_base['charge'] = self.mi.get_charge()
        data_base['wavelength'] = self.mi.get_wavelangth()
        # orbit !!!
        """
        orbit = []
        if self.lat != None:
            orbit = self.read_bpms()
        data_base["orbit"] = orbit
        """

        if flag == "start":
            print(flag)
            self.data = []
            self.data.append(data_base)
            print("#### ***** ", flag, len(self.data) #, data_base
                  )
            return True
        else:
            self.data.append(data_base)
            print("#### ***** ", flag,  len(self.data)  #, data_base
                  )
            self.wasSaved = self.send_to_db()
            self.data = []
            return self.wasSaved


    def new_tuning(self):
        self.db = PerfDB(dbname=self.dbname)
        sexts =  self.read_magnets(["sextupole"])
        quands = self.read_magnets(["quadrupole"])
        cors =   self.read_magnets(["hcor", "vcor"])
        bends =  self.read_magnets(["bend", "rbend", "sbend"])
        cavs = self.read_cavs()

        charge = self.mi.get_charge()
        wl = self.mi.get_wavelangth()

        self.db.new_tuning({'wl': wl, 'charge': charge, 'comment': 'test tuning'}) # creates new tuning record (e.g. for each shift);
        tunings = self.db.get_tunings()
        print ('current tunings', [(t.id, t.time, t.charge, t.wl) for t in tunings])
        tune_id = self.db.current_tuning_id()
        print ('current id', tune_id)

        mach_par = {}
        mach_par["bc2_pyros"] = self.mi.get_bc2_pyros()[0]
        mach_par["bc2_pyros.fine"] = self.mi.get_bc2_pyros()[1]
        mach_par["bc3_pyros"] = self.mi.get_bc3_pyros()[0]
        mach_par["bc3_pyros.fine"] = self.mi.get_bc3_pyros()[1]
        mach_par["final_energy"] = self.mi.get_final_energy()
        mach_par["solenoid"] = self.mi.get_sol_value()
        mach_par["nbunches"] = self.mi.get_nbunches()
        mach_par["gun_energy"] = self.mi.get_gun_energy()
        mach_par["time"] = time.time()
        mach_par.update(sexts)
        mach_par.update(quands)
        mach_par.update(cors)
        mach_par.update(bends)
        mach_par.update(cavs)
        #print(mach_par)
        self.db.add_machine_parameters(tune_id, params = mach_par)
        print ('current machine parameters', self.db.get_machine_parameters(tune_id))




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
                          '10SMATCH','3SDUMP',
                          "5R.UND6"
                          ]
        self.mutex = Lock()


    def init_corrector_vals(self, correctors):
        vals = np.zeros(len(correctors))
        for i in range(len(correctors)):
            #self.mutex.acquire()
            vals[i] = self.get_value(correctors[i])
            #self.mutex.release()
        return vals

    def get_cav_ampl(self, cav):
        return pydoocs.read("FLASH.RF/LLRF.CONTROLLER/PVS." + cav + "/AMPL.SAMPLE")["data"]

    def get_cav_phase(self, cav):
        return pydoocs.read("FLASH.RF/LLRF.CONTROLLER/PVS." + cav + "/PHASE.SAMPLE")["data"]

    def get_cavity_info(self, cavs):
        ampls = np.zeros(len(cavs))
        phases = np.zeros(len(cavs))
        for i in range(len(cavs)):
            #ampl_channel = 'FLASH.RF/LLRF.CONTROLLER/CTRL.' + cavs[i] + '/SP.AMPL'
            #phase_channel = 'FLASH.RF/LLRF.CONTROLLER/CTRL.' + cavs[i] + '/SP.PHASE'
            ampls[i] = self.get_cav_ampl(cavs[i])
            phases[i] = self.get_cav_phase(cavs[i])
        return ampls, phases

    def get_gun_energy(self):
        gun_energy = pydoocs.read("FLASH.RF/LLRF.ENERGYGAIN.ML/GUN/ENERGYGAIN.FLASH1")['data']
        gun_energy = gun_energy*0.001 # MeV -> GeV
        return gun_energy

    def get_bpms_xy(self, bpms):
        X = [0.0]*len(bpms)#np.zeros(len(correctors))
        Y = [0.0]*len(bpms)
        for i in range(len(bpms)):
            mag_channel = 'TTF2.DIAG/ORBIT/' + bpms[i]# + '/PS'
            X[i] = pydoocs.read(mag_channel + "/X.FLASH1.PULSE.MEAN")['data']*0.001 # mm -> m
            Y[i] = pydoocs.read(mag_channel + "/Y.FLASH1.PULSE.MEAN")['data']*0.001 # mm -> m
        return X, Y

    def get_quads_current(self, quads):
        vals = np.zeros(len(quads))
        for i in range(len(quads)):
            mag_channel = 'TTF2.MAGNETS/QUAD/' + quads[i]# + '/PS'
            vals[i] = pydoocs.read(mag_channel + "/PS")['data']
        return vals

    def get_bends_current(self, bends):
        vals = [0.0]*len(bends)#np.zeros(len(correctors))
        for i in range(len(bends)):
            mag_channel = 'TTF2.MAGNETS/DIPOLE/' + bends[i]# + '/PS'
            vals[i] = pydoocs.read(mag_channel + "/PS")['data']
        return vals

    def get_sext_current(self, sext):
        vals = [0.0]*len(sext)#np.zeros(len(correctors))
        for i in range(len(sext)):
            mag_channel = "TTF2.MAGNETS/SEXT/" + sext[i]
            vals[i] = pydoocs.read(mag_channel + "/PS")['data']
        return vals

    def get_alarm_ACC39(self):
        val = self.get_cav_ampl("M1.ACC39")
        if val >= 20.2:
            return 1.
        elif val >= 20:
            return 0.8
        else:
            return 0.

    def get_alarms(self):
        alarm_vals = np.zeros(len(self.blm_names))
        for i in range(len(self.blm_names)):
            blm_channel = 'TTF2.DIAG/BLM/'+self.blm_names[i]+'/CH00.TD'
            blm_alarm_ch  = ('TTF2.DIAG/BLM/'+self.blm_names[i]).replace('BLM', 'BLM.ALARM') + '/THRFHI'
            if (self.debug): print('reading alarm channel', blm_alarm_ch)
            alarm_val = pydoocs.read(blm_alarm_ch)['data'] * 1.25e-3 # alarm thr. in Volts
            if (self.debug): print ('alarm:', alarm_val)
            sample = pydoocs.read(blm_channel)['data']
            h = np.array([x[1] for x in sample])

            alarm_vals[i] = np.max( np.abs(h) ) / alarm_val

        #alarm_vals = np.append(alarm_vals, self.get_alarm_ACC39() )
        #print(alarm_vals)
        return alarm_vals

    def get_sase(self, detector='gmd_default'):
        #print(detector)
        if detector == 'mcp':
            # incorrect
            return pydoocs.read('TTF2.DIAG/MCP.HV/MCP.HV1/HV_CURRENT')['data']
            #return np.abs( np.mean(h) )
        if detector == 'gmd_fl1_slow':
            return pydoocs.read('TTF2.FEL/BKR.FLASH.STATE/BKR.FLASH.STATE/SLOW.INTENSITY' )['data']
        if detector == "bkr":
            # default 'BKR' gmd
            h = np.array(pydoocs.read('TTF2.FEL/BKR.FLASH.STATE/BKR.FLASH.STATE/ENERGY.CLIP.SPECT')['data'])
            val = np.mean(np.array([x[1] for x in h]))
            return val
        return pydoocs.read("TTF2.FEL/PHFLUX/TUNNEL/ENPULSEIC")['data']



    def get_sase_pos(self):

        x1 = pydoocs.read('TTF2.FEL/GMDPOSMON/TUNNEL/IX.POS')['data']
        y1 = pydoocs.read('TTF2.FEL/GMDPOSMON/TUNNEL/IY.POS')['data']

        x2 = pydoocs.read('TTF2.FEL/GMDPOSMON/BDA/IX.POS')['data']
        y2 = pydoocs.read('TTF2.FEL/GMDPOSMON/BDA/IY.POS')['data']
    
        return [ (x1,y1), (x2,y2) ] 

    def get_spectrum(self, f=None, detector='tunnel_default'):

        f_min = 13.0 # spectrum window (nm). TODO: replace with readout
        f_max = 14.0
        
        spec = np.array(pydoocs.read('TTF2.EXP/PBD.PHOTONWL.ML/WAVE_LENGTH/VAL.TD')['data'])
    
        if f == None:
            f = np.linspace(f_min, f_max, len(spec))
    
        return f, spec

    def get_charge(self):
        charge = pydoocs.read('TTF2.FEEDBACK/LONGITUDINAL/MONITOR2/MEAN_AVG')
        #print("charge = ", charge["data"], " nQ")
        return charge["data"]

    def get_wavelangth(self):
        wavelength = pydoocs.read('TTF2.DAQ/ENERGY.DOGLEG/LAMBDA_MEAN/VAL')
        #print("wavelength = ", wavelength["data"], "nm")
        return wavelength["data"]

    def get_bc2_pyros(self):
        bc2_pyro = pydoocs.read( "FLASH.DIAG/BCM/9DBC2.1/CH00.FLASH1")
        bc2_pyro_fine = pydoocs.read("FLASH.DIAG/BCM/9DBC2.2/CH00.FLASH1")
        #print("BC2 comp_fine = ", bc2_pyro["data"], bc2_pyro_fine["data"])
        return (bc2_pyro["data"], bc2_pyro_fine["data"])

    def get_bc3_pyros(self):
        bc3_pyro = pydoocs.read("FLASH.DIAG/BCM/4DBC3.1/CH00.FLASH1")
        bc3_pyro_fine = pydoocs.read("FLASH.DIAG/BCM/4DBC3.2/CH00.FLASH1")
        #print("BC3 comp_fine = ", bc3_pyro["data"], bc3_pyro_fine["data"])
        return (bc3_pyro["data"], bc3_pyro_fine["data"])


    def get_final_energy(self):
        final_energy = pydoocs.read("TTF2.FEEDBACK/LONGITUDINAL/MONITOR11/MEAN_AVG")
        #print("final_energy = ", final_energy["data"], "MeV")
        return final_energy["data"]

    def get_sol_value(self):
        sol = pydoocs.read("TTF2.MAGNETS/SOL/1GUN/PS")
        #print("sol = ", sol["data"], "A")
        return sol["data"]

    def get_nbunches(self):
        nbunches = pydoocs.read("FLASH.DIAG/TIMER/FLASHCPUTIME1.0/NUMBER_BUNCHES.1")
        #print("nbunches = ", nbunches["data"])
        return nbunches["data"]

    def sase_fast(self):
        nb = self.get_nbunches()
        ch = 'TTF2.FEL/PHFLUX/TUNNEL/PHOTONPULSE.FF.SPECT'
        data = pydoocs.read(ch)['data'][700:700 + nb - 1]
        v = np.mean(np.array([x[1] for x in data]))
        return v


    def get_value_ps(self, device_name):
        ch = 'TTF2.MAGNETS/STEERER/' + device_name + '/PS'
        self.mutex.acquire()
        val = pydoocs.read(ch)['data']
        self.mutex.release()
        return val

    def get_value(self, device_name):
        #print(device_name, device_name.find("ACC"))
        if device_name.find("ACC")==0:
            ch = 'FLASH.RF/LLRF.CONTROLLER/CTRL.' + device_name
        elif device_name.find("knob")>=0:
            pass
            #return get_knob_value(device_name)
        elif device_name.find("SUMVOLTAGE")>=0:
            ch = 'FLASH.RF/LLRF.SUMVOLTAGE_CTRL/ACC139/' + device_name + '.SP.FLASH1'
        else:
            ch = 'TTF2.MAGNETS/STEERER/' + device_name + '/PS.RBV'
        #print("getting value = ", ch)
        #self.mutex.acquire()
        try:
            val = pydoocs.read(ch)['data']
            return val
        except:
            print("Error in pydoocs.read(): ch = ", ch)
        #self.mutex.release()

    
    def set_value(self, device_name, val):
        if device_name.find("ACC")>=0:
            ch = 'FLASH.RF/LLRF.CONTROLLER/CTRL.' + device_name
        elif device_name.find("SUMVOLTAGE")>=0:
            ch = 'FLASH.RF/LLRF.SUMVOLTAGE_CTRL/ACC139/' + device_name + '.SP.FLASH1'
        else:
            ch = 'TTF2.MAGNETS/STEERER/' + device_name + '/PS'
        print (ch, val)
        self.mutex.acquire()
        try:
            pydoocs.write(ch, str(val))
            #pass
        except:
            print("Error in pydoocs.write(): ch = ", ch, " val = ", val)
        self.mutex.release()
        return 0


class FLASH1DeviceProperties:
    def __init__(self):
        self.stop_exec = False
        self.save_machine = False
        self.patterns = {}
        self.limits = {}
        """
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
        self.limits['H8TCOL'] = [0.04,0.05]

        self.patterns['V8TCOL'] = re.compile('V8TCOL')
        self.limits['V8TCOL'] = [0.033,0.043]
        
        self.patterns['V1ORS'] = re.compile('V1ORS')
        self.limits['V1ORS'] = [0.1, 0.15]

        self.patterns['H5ORS'] = re.compile('H5ORS')
        self.limits['H5ORS'] = [-0.06, -0.01]

        self.patterns['H10ORS'] = re.compile('H10ORS')
        self.limits['H10ORS'] = [-0.2, -0.1]

        self.patterns['V12ORS'] = re.compile('V12ORS')
        self.limits['V12ORS'] = [0.04, 0.15]
        """

    def set_limits(self, dev_name, limits):
        self.patterns[dev_name] = re.compile(dev_name)
        #print(self.patterns[dev_name])
        self.limits[dev_name] = limits
        #print("inside dp set = ", self.patterns[dev_name], self.limits)

    def get_limits(self, device):
        #print(self.limits)
        for k in self.patterns.keys():
            #print('testing', k)
            if self.patterns[k].match(device) != None:
                #print("inside dp get = ", device, self.limits[k])
                return self.limits[k]
        return [-2, 2]

    def get_polarity(self, quads):
        vals = [0.0]*len(quads)#np.zeros(len(correctors))
        for i in range(len(quads)):
            mag_channel = 'TTF2.MAGNETS/QUAD/' + quads[i]# + '/PS'
            vals[i] = pydoocs.read(mag_channel + "/PS.Polarity")['data']
        return vals

    def get_type_magnet(self, quads):
        vals = [0.0]*len(quads)#np.zeros(len(correctors))
        for i in range(len(quads)):
            mag_channel = 'TTF2.MAGNETS/QUAD/' + quads[i]# + '/PS'
            vals[i] = pydoocs.get(mag_channel + "/DEVTYPE")['data']
        return vals


'''
test interface
'''
class TestInterface:
    def __init__(self):
        self.dev = {}

    def get_alarms(self):
        return np.random.rand(4)*0.2 #0.0, 0.0, 0.0, 0.0]

    def get_sase(self, detector=None):
        #print("Detector", detector)
        S = 0.
        for name in self.dev:
            S += (self.dev[name] - 0.3)**2 + np.random.rand(1)[0]*0.1
        return S

    def init_corrector_vals(self, correctors):
        vals = []
        for name in correctors:
            vals.append(self.get_value( name))
        return vals

    def get_cor_value(self, devname):
        return np.random.rand(1)[0]

    def get_value(self, device_name):
        #print("get value", device_name)
        try:
            val = self.dev[device_name]
        except:
            val = 0.
        return val#np.random.rand(1)[0]

    def set_value(self, device_name, val):
        self.dev[device_name] = val
        return 0.0

    def get_quads_current(self, device_names):
        return np.random.rand(len(device_names))

    def get_bpms_xy(self, bpms):
        X = np.zeros(len(bpms))
        Y = np.zeros(len(bpms))
        return X, Y

    def get_sase_pos(self):
        return [(0,0), (0, 0)]

    def get_gun_energy(self):
        return 0.

    def get_cavity_info(self, devnames):
        X = np.zeros(len(devnames))
        Y = np.zeros(len(devnames))
        return X, Y

    def get_charge(self):
        return 0.

    def get_wavelangth(self):
        return 0.

    def get_bc2_pyros(self):
        return (0, 0)

    def get_bc3_pyros(self):
        return (0, 0)

    def get_final_energy(self):
        return 0

    def get_sol_value(self):
        time.sleep(0.1)
        return 0

    def get_nbunches(self):
        return 0
    def get_value_ps(self, device_name):
        return 0.
