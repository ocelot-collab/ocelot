"""
Sergey Tomin, XFEL/DESY, 2017
"""
from ocelot.optimizer.mint.opt_objects import Device
import numpy as np
import time
from ocelot import *

class MIMagnets(Device):
    def __init__(self, eid=None):
        super(MIMagnets, self).__init__(eid=eid)
        self.time_delay = 0.1         # sec
        self.data = []

    def read(self):
        ch = "XFEL.MAGNETS/MAGNET.ML/*/KICK_MRAD.SP"
        #ch = "XFEL.MAGNETS/MAGNET.ML/*/KICK.RBV"
        data = self.mi.get_value(ch)
        return data

    def data2dict(self, data):
        mag_dict = {}
        for line in data:
            mag_dict[line["str"]] = line["float1"]
        return mag_dict

    def get_magnets(self, magnets):
        """
        All bpm works with [m] but doocs gives position in [mm]

        :param bpms: list of BPM objects
        :param charge_threshold:
        :return:
        """
        data = self.read()
        if len(data) == 0:
            return False
        mag_dict = self.data2dict(data)
        mag_names = mag_dict.keys()

        for i, mag in enumerate(magnets):
            if mag.__class__ == Quadrupole:
                if mag.id not in mag_names:
                    #print("Magnet is not in DOOCS: ", mag.id)
                    mag.ui.set_alarm(True)
                else:
                    val = mag_dict[mag.id]/1000.
                    mag.ui.set_value(val/mag.l)
                    mag.ui.set_alarm(False)
            if mag.__class__ in [SBend, Bend, RBend]:
                if mag.id not in mag_names:
                    #print("Magnet is not in DOOCS: ", mag.id)
                    mag.ui.set_alarm(True)
                else:
                    val = mag_dict[mag.id]
                    mag.ui.set_value(val)
                    mag.ui.set_alarm(False)
        return True

class MIBendBC(MIMagnets):
    def __init__(self, eid=None):
        super(MIBendBC, self).__init__(eid=eid)
        self.time_delay = 0.1         # sec
        self.data = []

    def get_magnets(self, magnets):
        """
        All bpm works with [m] but doocs gives position in [mm]

        :param bpms: list of BPM objects
        :param charge_threshold:
        :return:
        """
        data = self.read()
        if len(data) == 0:
            return False
        mag_dict = self.data2dict(data)
        mag_names = mag_dict.keys()

        for i, mag in enumerate(magnets):
            if mag.id not in mag_names:
                #print("Magnet is not in DOOCS: ", mag.id)
                mag.ui.set_alarm(True)
            else:
                val = mag_dict[mag.id]/1000.
                # mrad -> rad
                mag.ui.set_value(val)
                mag.ui.set_alarm(False)
        return True


class MICavity(MIMagnets):
    def __init__(self, eid=None):
        super(MIMagnets, self).__init__(eid=eid)
        self.time_delay = 0.1         # sec
        self.data = []

    def read_v(self):
        ch = "XFEL.RF/LLRF.CONTROLLER/*/SP.AMPL"
        data = self.mi.get_value(ch)
        return data
    
    def read_phi(self):
        ch = "XFEL.RF/LLRF.CONTROLLER/*/SP.PHASE"
        data = self.mi.get_value(ch)
        return data

    def data2dict(self, data):
        mag_dict = {}
        for line in data:
            #  use only name of the cav e.g. A1, AH1, A2, A5, ... 
            mag_dict[line["str"].split(".")[1]] = line["float1"]
        return mag_dict

    def get_cavities(self, cavs):
        """
        All bpm works with [m] but doocs gives position in [mm]

        :param bpms: list of BPM objects
        :param charge_threshold:
        :return:
        """
        data_v = self.read_v()
        data_phi = self.read_phi()
        if len(data_v) == 0 or len(data_phi) == 0:
            return False
        volt_dict = self.data2dict(data_v)
        phi_dict = self.data2dict(data_phi)
        cav_names = volt_dict.keys()
        print(cav_names)
        for i, cav in enumerate(cavs):
            if cav.id not in cav_names:
                print("Magnet is not in DOOCS: ", cav.id)
                cav.ui.set_volt(0)
                cav.ui.set_phi(0)
                cav.ui.set_alarm(True)
            else:
                v = volt_dict[cav.id]
                phi = phi_dict[cav.id]
                cav.ui.set_volt(v/1000.)
                cav.ui.set_phi(phi)
                cav.ui.phi_was_changed(np.abs(np.abs(cav.ui.get_init_phi() - phi)) > 0.01)
                cav.ui.volt_was_changed(np.abs(np.abs(cav.ui.get_init_volt() - v)) > 0.01)
                cav.ui.set_alarm(False)
        return True


class MIOrbit(Device):
    def __init__(self, eid=None):
        super(MIOrbit, self).__init__(eid=eid)
        self.bpm_server = "BPM" # "ORBIT"     # or "BPM"
        self.time_delay = 0.1         # sec
        self.charge_threshold = 0.005 # nC
        self.bpm_names = []
        self.x = []
        self.y = []
        self.mean_x = []
        self.mena_y = []
        self.mean_charge = []
        #self.charge = []

    def read_positions(self):
        #try:
        orbit_x = self.mi.get_value("XFEL.DIAG/" + self.bpm_server + "/*/X.SA1")
        orbit_y = self.mi.get_value("XFEL.DIAG/" + self.bpm_server + "/*/Y.SA1")
        #except:
        #    print("ERROR: reading from DOOCS")
        #    return False
        #print(orbit_x)
        names_x = [data["str"] for data in orbit_x]
        names_y = [data["str"] for data in orbit_y]
        if not np.array_equal(names_x, names_y):
            print("X and Y orbits are not equal")
        self.x = np.array([data["float"] for data in orbit_x])
        self.y = np.array([data["float"] for data in orbit_y])
        return [names_x, self.x, self.y]

    def read_charge(self):
        #try:
        charge = self.mi.get_value("XFEL.DIAG/BPM/*/CHARGE.SA1")
        #except:
        #    print("ERROR: reading from DOOCS")
        #    return False
        names = [data["str"] for data in charge]
        values = np.array([data["float"] for data in charge])
        return names, values

    def read_orbit(self):
        names_xy, x, y = self.read_positions()
        names_charge, charge = self.read_charge()
        #print(names_xy)
        #print(names_charge)
        if not np.array_equal(names_xy, names_charge):
            print("CHARGE reading and POSITIONS are not equal")
            #return False
        return names_xy, x, y, charge


    def read_and_average(self, nreadings, take_last_n):
        orbits_x = []
        orbits_y = []
        orbits_charge = []
        saved_names = []
        for i in range(nreadings):
            names, x, y, charge = self.read_orbit()
            orbits_x.append(x)
            orbits_y.append(y)
            orbits_charge.append(charge)
            if i > 0:
                if not np.array_equal(saved_names, names):
                    print("error: different ")
            saved_names = names
            time.sleep(self.time_delay)
        self.bpm_names = saved_names
        self.mean_x = np.mean(orbits_x[-take_last_n:], axis=0)
        self.mean_y = np.mean(orbits_y[-take_last_n:], axis=0)
        self.mean_charge = np.mean(orbits_charge[-take_last_n:], axis=0)
        return self.bpm_names, self.mean_x, self.mean_y, self.mean_charge

    def get_bpms(self, bpms):
        """
        All bpm works with [m] but doocs gives position in [mm]

        :param bpms: list of BPM objects
        :param charge_threshold:
        :return:
        """
        if len(self.bpm_names) == 0:
            return False
        #bpm_names = [bpm.id for bpm in bpms]
        print("bpm", len(bpms))
        indxs = []
        valid_bpm_inx = []
        for i, bpm in enumerate(bpms):
            if bpm.id not in self.bpm_names:
                bpm.ui.uncheck()
            else:
                valid_bpm_inx.append(i)
                indxs.append(self.bpm_names.index(bpm.id))
        print("bpm", len(bpms), len(indxs))
        bpms = [bpms[indx] for indx in valid_bpm_inx]
        for i, bpm in enumerate(bpms):
            inx = indxs[i]
            bpm.x = self.mean_x[inx]/1000      # [mm] -> [m]
            bpm.y = self.mean_y[inx]/1000      # [mm] -> [m]
            bpm.charge = self.mean_charge[inx] # nC
        return True
