__author__ = 'Sergey Tomin'
import pickle
from ocelot.utils.mint.flash1_converter import *
import numpy as np
import time

def show_currents( elems, alpha):
    print("******* displaying currents - START ********")
    for elem in elems:
        if elem.dI == 0:
            continue
        n = len(elem.id)
        n2 = len(str(elem.I + elem.dI))
        n3 = len(str(elem.I))
        print(elem.id, " "*(10-n) + "<-- ", elem.I + elem.dI,  " "*(18-n2)+ " was = ", elem.I, " "*(18-n3) + " dI = ", elem.dI, "x", alpha)
    print("******* displaying currents - END ********")

def set_currents(mi, elems, alpha):
    for elem in elems:
        if elem.dI == 0:
            continue
        n = len(elem.id)
        new_I = elem.I + elem.dI*alpha
        n2 = len(str(new_I))
        print(elem.id,  " "*(10-n) + "<-- ", new_I,  " "*(18-n2)+ " was = ", elem.I)
        mi.set_value(elem.mi_id, new_I)

def restore_current(mi, elems):
    for elem in elems:
        if elem.dI == 0:
            continue
        n = len(elem.id)
        print(elem.id, " "*(10-n) +"<-- ", elem.I)
        mi.set_value(elem.mi_id, elem.I)

def currents2angles(orb):
    for elem in np.append(orb.hcors, orb.vcors):

        angle = tpi2k(elem.dev_type, elem.E, elem.dI)*0.001

        #print elem.id, elem.dI, angle, elem.E
        elem.angle = angle
        #dI = tpk2i(elem.dev_type, elem.E, elem.angle*1000.)
        if abs(angle) > 1e-10:
            elem.angle = angle
            if elem.id in ['H10ACC5', 'H10ACC6', 'V10ACC5', 'V10ACC6', 'V2SFELC']:
                elem.angle = angle*2.
            #print elem.id, "angle=", elem.angle, " dI = ", elem.dI, " I = ", elem.I
        else:
            elem.dI = 0.
            elem.angle = 0.
        if abs(elem.angle) > 0.005:
            print(elem.id, " @@@@@@@@@@@@@@@@ HIGH CURRENT @@@@@@@@@@@@@@@ = ", elem.angle)

def angles2currents(orb):
    for elem in np.append(orb.hcors, orb.vcors):

        #print elem.id
        dI = tpk2i(elem.dev_type, elem.E, elem.angle*1000.)
        #print elem.id, dI, elem.angle, elem.E
        if abs(dI) > 0.00005:
            elem.dI = dI
            if elem.id in ['H10ACC5', 'H10ACC6', 'V10ACC5', 'V10ACC6', 'V2SFELC']:
                elem.dI = dI/2.
            #print elem.id, "angle=", elem.angle, " dI = ", elem.dI, " I = ", elem.I
        else:
            elem.dI = 0.
            elem.angle = 0.
        if abs(elem.dI) > 0.5:
            print(elem.id, " @@@@@@@@@@@@@@@@ HIGH CURRENT @@@@@@@@@@@@@@@ = ", elem.dI)


class HighLevelInterface:

    def __init__(self, lattice, mi, dp):
        self.lat = lattice
        self.mi = mi
        self.dp = dp

    def read_all(self):
        self.lat.gun_energy = self.mi.get_gun_energy()
        self.lat.sase = self.mi.get_sase(detector="gmd_fl1_slow")
        self.timestamp = time.time()
        self.read_cavs()
        self.read_quads()
        self.read_bends()
        self.read_cors()
        self.read_sexts()
        self.read_bpms()

    def read_quads(self):
        id2I_dict = {}
        for elem in self.lat.sequence:
            if elem.type == "quadrupole":
                name = elem.id
                name = name.replace("_U", "")
                name = name.replace("_D", "")
                name = name.replace("_", ".")
                try:
                    elem.mi_id
                except:
                    elem.mi_id = name
                elem.I = 0
                elem.polarity = 1
                if elem.mi_id in id2I_dict.keys():
                    elem.I = id2I_dict[elem.mi_id]["I"]
                    elem.polarity = id2I_dict[elem.mi_id]["polarity"]
                else:
                    try:
                        elem.I = self.mi.get_quads_current([elem.mi_id])[0]
                        #elem.polarity = dp.get_polarity([elem.mi_id])[0]
                        id2I_dict[elem.mi_id] = {}
                        id2I_dict[elem.mi_id]["I"] = elem.I
                        id2I_dict[elem.mi_id]["polarity"] = elem.polarity
                    except:
                        print("* ", name, "  CAN MOT FIND")

    def read_bends(self):
        id2I_dict = {}
        for elem in self.lat.sequence:
            if elem.type in ["bend", "sbend", "rbend"]:
                #name = elem.id
                if "BC2" in elem.id:
                    elem.mi_id = "D1BC2"
                elif "BC3" in elem.id:
                    elem.mi_id = "D1BC3"
                elif "ECOL" in elem.id:
                    elem.mi_id = "D1ECOL"
                elif "SMATCH" in elem.id:
                    elem.mi_id = "D9SMATCH"
                else:
                    continue

                elem.I = 0
                elem.polarity = 1
                if elem.mi_id in id2I_dict.keys():
                    elem.I = id2I_dict[elem.mi_id]["I"]
                else:
                    try:
                        elem.I = self.mi.get_bends_current([elem.mi_id])[0]
                        id2I_dict[elem.mi_id] = {}
                        id2I_dict[elem.mi_id]["I"] = elem.I
                    except:
                        print("* ", elem.id, "  CAN MOT FIND")

    def read_cavs(self):
        dict_cavity = {}
        for elem in self.lat.sequence:
            if elem.type == "cavity":
                name = elem.id.split("_")
                elem.mi_id = name[-2] + "." + name[-1]
                try:
                    if elem.mi_id in dict_cavity.keys():
                        ampls = dict_cavity[elem.mi_id]["ampl"]
                        phases = dict_cavity[elem.mi_id]["phi"]
                    else:
                        ampls, phases = self.mi.get_cavity_info([elem.mi_id])
                        dict_cavity[elem.mi_id] = {}
                        dict_cavity[elem.mi_id]["ampl"] = ampls[0]
                        dict_cavity[elem.mi_id]["phi"] = phases[0]
                except:
                    print("* UNKNOWN cav", elem.mi_id, elem.id)
                    continue
                ampls = dict_cavity[elem.mi_id]["ampl"]
                phases = dict_cavity[elem.mi_id]["phi"]

                elem.v = ampls*0.001 # MeV -> GeV
                elem.phi = phases
                if elem.mi_id == "M1.ACC1":
                    elem.v = elem.v/8.
                    if abs(elem.phi) > 10:
                        print("* too large phase on ", elem.mi_id, elem.phi)
                elif elem.mi_id == "M1.ACC39":
                    # deaccelerator
                    elem.v = elem.v/4.
                    elem.phi = phases + 180.
                elif "ACC23" in elem.mi_id:
                    elem.v = elem.v/8.
                    if abs(elem.phi) > 10:
                        print("* too large phase on ", elem.mi_id, elem.phi, "  # set zero phase: ", elem.mi_id,".phi <-- 0")
                        elem.phi = 0.
                elif "ACC45" in elem.mi_id :
                    elem.v = elem.v/8.
                    if abs(elem.phi) > 10:
                        print("* too large phase on ", elem.mi_id, elem.phi)
                elif "ACC67" in elem.mi_id:
                    elem.v = elem.v/8.
                    if abs(elem.phi) > 10:
                        print("* too large phase on ", elem.mi_id, elem.phi)
        for cav in dict_cavity.keys():
            print (cav, " E = ", dict_cavity[cav]["ampl"], "MeV, phi =  ", dict_cavity[cav]["phi"], " grad")
        self.lat.update_transfer_maps()
        return self.lat

    def read_cors(self):
        id2I_dict = {}
        for elem in self.lat.sequence:
            if elem.type in ["hcor", "vcor"]:
                name = elem.id
                name = name.replace("_", ".")
                try:
                    elem.mi_id
                except:
                    elem.mi_id = name
                if elem.mi_id in id2I_dict.keys():
                    elem.I = id2I_dict[elem.mi_id]
                else:
                    try:
                        #print elem.mi_id, elem.id
                        vals = self.mi.init_corrector_vals([elem.mi_id])
                        elem.I = vals[0]
                        id2I_dict[elem.mi_id] = elem.I
                        #print elem.I
                    except:
                        print("* ", elem.mi_id, "UNKNOW")
                        elem.type = "drift"

    def read_sexts(self):
        id2I_dict = {}
        for elem in self.lat.sequence:
            if elem.type =="sextupole":
                if elem.id == "S2ECOL":
                    elem.mi_id = "S2.6ECOL"
                    vals = self.mi.get_sext_current([elem.mi_id])
                    elem.I = vals[0]
                    id2I_dict[elem.mi_id] = elem.I
                elif elem.id == "S6ECOL":
                    elem.mi_id = "S2.6ECOL"
                    if elem.mi_id in id2I_dict.keys():
                        elem.I = -id2I_dict[elem.mi_id]
                    else:
                        vals = self.mi.get_sext_current([elem.mi_id])
                        elem.I = -vals[0]

    def read_bpms(self):
        X = []
        Y = []
        for elem in self.lat.sequence:
            if elem.type == "monitor":
                name = elem.id.replace("BPM", "")
                elem.mi_id = name
                try:

                    x, y = self.mi.get_bpms_xy([elem.mi_id])
                    elem.x = x[0]
                    elem.y = y[0]
                    X.append(elem.x)
                    Y.append(elem.y)
                except:
                    print("* ", name, "  CAN MOT FIND")
        return np.array(X), np.array(Y)


class MachineSetup:
    def __init__(self, lattice, mi=None, dp=None):
        self.lat = lattice
        self.mi = mi
        self.dp = dp
        self.hli = HighLevelInterface(lattice, mi, dp)
        pass

    def read_save_lattice(self, filename):
        print("reading currents from DOOCS system ... ",)
        self.hli.read_all()
        print("OK")

        print("getting currents of quad ... ")
        self.dict_quad = self.get_elem_type_currents(self.lat, ["quadrupole"])
        print("OK")

        print("getting currents of bend ... ")
        self.dict_bend = self.get_elem_type_currents(self.lat, ["sbend", "rbend", "bend"])
        print("OK")

        print("getting currents of correctors ... ")
        self.dict_cor = self.get_elem_type_currents(self.lat, ["hcor", "vcor"])
        print("OK")

        print("getting params of cavities ... ")
        self.dict_cav = self.get_cav_params(self.lat)
        print("OK")

        print("getting beam positions ...  ")
        self.dict_orbit = self.get_orbit(self.lat)
        print("OK")

        print("getting currents of sext  ...  ")
        self.dict_sext = self.get_elem_type_currents(self.lat, ["sextupole"])
        print("OK")

        data = {}
        data["quad"] = self.dict_quad
        data["bend"] = self.dict_bend
        data["cor"] = self.dict_cor
        data["cav"] = self.dict_cav
        data["orbit"] = self.dict_orbit
        data["sext"] = self.dict_sext
        data["sase"] = self.lat.sase
        data["gun_energy"] = self.lat.gun_energy
        data["timestamp"] = self.hli.timestamp
        pickle.dump(data, open(filename, "wb"))

    def load_lattice(self, filename, lattice):
        data = pickle.load(open(filename, "rb"))
        self.dict_quad = data["quad"]
        self.dict_bend = data["bend"]
        self.dict_cor = data["cor"]
        self.dict_cav = data["cav"]
        self.dict_orbit = data["orbit"]
        self.dict_sext = data["sext"]
        try:
            lattice.sase = data["sase"]
        except:
            print("SASE value was not logged")
        try:
            lattice.gun_energy = data["gun_energy"]
        except:
            print("gun energy was not logged")

        for elem in lattice.sequence:

            if elem.type == "quadrupole" and elem.id in self.dict_quad.keys():
                elem.mi_id = self.dict_quad[elem.id]["mi_id"]
                elem.dev_type = self.dict_quad[elem.id]["dev_type"]
                elem.I = self.dict_quad[elem.id]["I"]

            if elem.type in ["bend", "rbend", "sbend"] and elem.id in self.dict_bend.keys():
                elem.mi_id = self.dict_bend[elem.id]["mi_id"]
                elem.dev_type = self.dict_bend[elem.id]["dev_type"]
                elem.I = self.dict_bend[elem.id]["I"]

            if elem.type in ["hcor", "vcor"] and elem.id in self.dict_cor.keys():
                elem.mi_id = self.dict_cor[elem.id]["mi_id"]
                elem.dev_type = self.dict_cor[elem.id]["dev_type"]
                elem.I = self.dict_cor[elem.id]["I"]

            if elem.type == "cavity" and elem.id in self.dict_cav.keys():
                elem.mi_id = self.dict_cav[elem.id]["mi_id"]
                elem.v = self.dict_cav[elem.id]["v"]
                elem.phi = self.dict_cav[elem.id]["phi"]

            if elem.type == "monitor" and elem.id in self.dict_orbit.keys():
                elem.mi_id = self.dict_orbit[elem.id]["mi_id"]
                elem.x = self.dict_orbit[elem.id]["x"]
                elem.y = self.dict_orbit[elem.id]["y"]

            if elem.type == "sextupole" and elem.id in self.dict_sext.keys():
                elem.mi_id = self.dict_sext[elem.id]["mi_id"]
                elem.dev_type = self.dict_sext[elem.id]["dev_type"]
                elem.I = self.dict_sext[elem.id]["I"]
        return lattice.update_transfer_maps()

    def set_elem_energy(self, lat, init_energy):
        E = init_energy
        for elem in lat.sequence:
            E += elem.transfer_map.delta_e
            elem.E = E


    def get_elem_type_currents(self, lat, type):
        data = {}
        for elem in lat.sequence:
            if elem.type in type:
                try:
                    I = elem.I
                except:
                    print (elem.type, elem.id, " No current!")
                    continue
                data[elem.id] = {}
                data[elem.id]["type"] = elem.type
                data[elem.id]["mi_id"] = elem.mi_id
                data[elem.id]["dev_type"] = elem.dev_type
                data[elem.id]["I"] = I
        #pickle.dump(data, open(filename, "wb"))
        return data

    def get_cav_params(self, lat):
        data = {}
        for elem in lat.sequence:
            if elem.type == "cavity":
                try:
                    v = elem.v
                    phi = elem.phi
                except:
                    print (elem.type, elem.id, " No elem.v or elem.phi!")
                    continue
                data[elem.id] = {}
                data[elem.id]["type"] = elem.type
                data[elem.id]["mi_id"] = elem.mi_id
                data[elem.id]["phi"] = phi
                data[elem.id]["v"] = v
        #pickle.dump(data, open(filename, "wb"))
        return data

    def get_orbit(self, lat):
        data = {}
        for bpm in lat.sequence:
            if bpm.type == "monitor":
                try:
                    mi_id = bpm.mi_id
                except:
                    print ("bpm: ", bpm.id, " No mi_id")
                    continue
                data[bpm.id] = {}
                data[bpm.id]["type"] = "monitor"
                data[bpm.id]["mi_id"] = bpm.mi_id
                data[bpm.id]["x"] = bpm.x
                data[bpm.id]["y"] = bpm.y
        #pickle.dump(data, open(filename, "wb"))
        return data

    def save_orbit(self, lat, filename):
        print("getting beam positions ... ")
        self.dict_orbit = self.get_orbit(lat)
        print("OK")
        pickle.dump(self.dict_orbit, open(filename, "wb"))

    def load_orbit(self, filename, lat):
        data = pickle.load(open(filename, "rb"))
        data = data["orbit"]
        #print data
        for elem in lat.sequence:
            if elem.type == "monitor":
                if elem.id in data.keys():
                    elem.mi_id = data[elem.id]["mi_id"]
                    elem.x = data[elem.id]["x"]
                    elem.y = data[elem.id]["y"]
        return lat

    def set_orbit(self, lat):
        data = self.dict_orbit
        for elem in lat.sequence:
            if elem.type == "monitor":
                if elem.id in data.keys():
                    elem.mi_id = data[elem.id]["mi_id"]
                    elem.x = data[elem.id]["x"]
                    elem.y = data[elem.id]["y"]
        return lat


    def convert_currents(self, lat, init_energy):
        E = init_energy
        for elem in lat.sequence:
            E += elem.transfer_map.delta_e
            elem.E = E
            if elem.type == "quadrupole":

                k1 = tpi2k(elem.dev_type, elem.E, elem.I)
                k1 = abs(k1)*np.sign(elem.k1)
                #if elem.mi_id in ["Q4DBC2","Q9ACC2", 'Q3.5ECOL', 'Q5UND1.3.5', "Q5UND2.4", 'Q6UND1']:
                #    k1 = abs(k1)*sign(elem.k1)
                #K1 = k1
                #print elem.id,  "i.k1=", elem.k1, " r.k1=", k1, "I=", elem.I, "E=", E
                #print(elem.id,  "ideal: k1 = ", elem.k1, " real k1 = ", K1, " dk/k = ", (K1-elem.k1)/elem.k1*100.)
                elem.k1 = k1
            elif elem.type == "sextupole":
                #pass

                k2 = tpi2k(elem.dev_type, elem.E, elem.I)
                #print elem.id,  "i.k1=", elem.k2, " r.k1=", k2, "I=", elem.I, "E=", E
                elem.k2 = k2
                #print elem.id, elem.k2
            elif elem.type in ["bend", "sbend", "rbend"]:
                try:
                    elem.dev_type
                except:
                    continue
                angle = tpi2k(elem.dev_type, elem.E, elem.I)
                angle = abs(angle)*np.sign(elem.angle)*np.pi/180.
                #print elem.id,  "i.a=", elem.angle, " r.a=", angle, "I=", elem.I, "E=", E
                elem.angle = angle

            #elif elem.type in ["hcor", "vcor"]:
            #    angle = tpi2k(elem.dev_type, E, elem.I)
            #    if angle == None:
            #        print(elem.id,  elem.I, E, angle, elem.dev_type)
            #    else:
            #        elem.angle = angle*0.001
        lat.update_transfer_maps()