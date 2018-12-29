"""
Section class for s2e tracking.
S.Tomin. XFEL/DESY. 2017
"""
from ocelot import *
from ocelot.adaptors.astra2ocelot import *
import numpy as np
import copy

#data_dir = "N:/4all/xxl/zagor/mpy_xxl"
data_dir = "/Users/zagor"
data_dir = "/Users/tomins/ownCloud/DESY/repository/forSergey"
#data_dir = "/Volumes/Promise RAID/UserFolders/zagor_xxl"

class SectionLattice:
    """
    High level class to work with SectionTrack()
    """
    def __init__(self, sequence, tws0=None, data_dir="."):
        """

        :param sequence: list of SectionTrack()
        """
        self.sec_seq = sequence
        self.elem_seq = None
        self.tws = None
        self.data_dir = data_dir
        self.initialize(tws0=tws0)

    def initialize(self, tws0=None):
        self.init_sections()
        self.tws = self.calculate_twiss(tws0)


    def init_sections(self):
        """
        Method initiates section and return dictionary with initialized sections

        :return: self.dict_sections - dictionary
        """
        self.dict_sections = {}
        self.elem_seq = []
        for sec in self.sec_seq:
            s = sec(self.data_dir)
            #s.data_dir = self.data_dir
            #s.update()
            self.dict_sections[sec] = s
            self.elem_seq.append(s.lattice.sequence)
        return self.dict_sections

    def calculate_twiss(self, tws0=None):
        """
        Method calculates twiss parameters for whole lattice [sequence of lattices]
        and assigns for each section initial Twiss() - section.tws0

        :param tws0: Twiss() - initial twiss parameters
        :return: return list of Twiss()
        """
        tws_whole = []
        for sec in self.sec_seq:
            s = self.dict_sections[sec]
            s.tws0 = tws0
            if s.tws0 != None:
                tws = twiss(s.lattice, tws0)
                tws0 = tws[-1]
                tws_whole = np.append(tws_whole, tws)
        return tws_whole

    def update_sections(self, sections, config=None):
        #blank_sections = sections.sections
        #sections = []
        new_sections = []
        for sec in sections:
            sec = self.dict_sections[sec]
            #sec = blank_sec()
            #sections.append(sec)
            if config != None and sec.__class__ in config.keys():
                conf = config[sec.__class__]
                if "rho" in conf.keys():
                    sec.update_bunch_compressor(rho=conf["rho"])
                if "phi" in conf.keys() and "v" in conf.keys():
                    sec.update_cavity(phi=conf["phi"], v=conf["v"])
                if "match" in conf.keys() and conf["match"] == True:
                    sec.apply_matching()
                if "SC" in conf.keys():
                    sec.sc_flag = conf["SC"]
                if "CSR" in conf.keys():
                    sec.csr_flag = conf["CSR"]
                if "smooth" in conf.keys():
                    sec.smooth_flag = conf["smooth"]
                if "wake" in conf.keys():
                    sec.wake_flag = conf["wake"]
            new_sections.append(sec)
        return new_sections

    def track_sections(self, sections, p_array, config=None, force_ext_p_array=False):
        self.update_sections(sections, config=config)
        #twis_track = []
        for i, sec in enumerate(sections):
            sec = self.dict_sections[sec]
            if i == 0 and sec.__class__ != self.sec_seq[0] and not force_ext_p_array:
                p_array = None
            p_array = sec.tracking(particles=p_array)
            #twis_track.append(sec.load_twiss_file())
        return p_array

    def load_twiss_track(self, sections):
        s = []
        bx = []
        by = []
        s0 = 0
        #tws_list = []
        for i, sec in enumerate(sections):
            sec = self.dict_sections[sec]
            tws = sec.load_twiss_file()
            s_tws = np.array(tws["s"])
            s_tws = s_tws - s_tws[0] + s0
            s0 = s_tws[-1]
            s = np.append(s, s_tws)
            bx = np.append(bx, tws["beta_x"])
            by = np.append(by, tws["beta_y"])
        return s, bx, by


class SectionTrack:
    def __init__(self, data_dir):

        self.lattice_name = ""
        self.lattice = None
        self.dipoles = None
        self.tws0 = None
        self.dipole_len = None
        self.bc_gap = None
        self.cav_name_pref = None

        self.unit_step = 1.0

        self.input_beam_file = None
        self.output_beam_file = None
        self.tws_file = None
        #self.data_dir = "."
        self.particle_dir = data_dir + "/particles/"
        self.tws_dir = data_dir + "/tws/"

        self.physics_processes_array = []  # list of physics process
        self.method = MethodTM()
        self.method.global_method = SecondTM
        self.translator = {SpaceCharge: "sc_flag", CSR: "csr_flag",
                           Wake: "wake_flag", BeamTransform:"bt_flag", SmoothBeam: "smooth_flag"}
        self.sc_flag = True
        self.csr_flag = True
        self.wake_flag = True
        self.bt_flag = True
        self.smooth_flag = True

        self.print_progress = True
        self.calc_tws = True
        self.kill_track = False
        #self.update()

    #def update(self):
    #    self.particle_dir = self.data_dir + "/particles/"
    #    self.tws_dir = self.data_dir + "/tws/"

    def apply_matching(self):
        if self.tws0 == None:
            print("TWISS is not defined")
            return
        self.tws0.mux = 0
        self.tws0.muy = 0
        bt = BeamTransform(tws=self.tws0)
        bt.bounds = [-0.5, 0.5]
        self.add_physics_process(bt, self.lattice.sequence[0], self.lattice.sequence[0])


    def bc_analysis(self):
        #print("BC analysis")
        # find positions
        L = 0.
        for elem in self.lattice.sequence:

            if elem in self.dipoles:
                elem.s_pos = L
            L += elem.l
        self.dipoles = np.array(self.dipoles)
        # sorting - put dipoles in right order
        sort_index = np.argsort(np.array([d.s_pos for d in self.dipoles]))
        #print(sort_index)
        self.dipoles = self.dipoles[sort_index]

        #self.bc_gap = self.dipoles[2].s_pos - self.dipoles[1].s_pos - self.dipoles[1].l
        #print("bc_gap", self.bc_gap)
        self.get_bc_shoulders()

    def get_bc_shoulders(self):
        self.left_shoulder = []
        self.right_shoulder = []
        left_flag = False
        right_flag = False
        for i, elem in enumerate(self.lattice.sequence):
            if elem == self.dipoles[0]:
                left_flag = True
            elif elem == self.dipoles[1]:
                left_flag = False

            if elem == self.dipoles[2]:
                right_flag = True
            elif elem == self.dipoles[3]:
                right_flag = False

            if left_flag and elem != self.dipoles[0] and (elem.__class__ not in [Edge]):
                self.lattice.sequence[i] = copy.deepcopy(elem)
                self.left_shoulder.append(self.lattice.sequence[i])
            if right_flag and elem != self.dipoles[2] and (elem.__class__ not in [Edge]):
                self.lattice.sequence[i] = copy.deepcopy(elem)
                self.right_shoulder.append(self.lattice.sequence[i])

        self.bc_gap_left = np.sum([d.l for d in self.left_shoulder])
        right_len = np.sum([d.l for d in self.right_shoulder])

        for d in self.left_shoulder:
            d.len_coef = d.l / self.bc_gap_left
        for d in self.right_shoulder:
            d.len_coef = d.l / right_len

    def change_bc_shoulders(self, drift):
        for d in self.left_shoulder:
            d.l = drift * d.len_coef

        for d in self.right_shoulder:
            d.l = drift * d.len_coef

    def update_bunch_compressor(self, rho):
        if self.dipoles is None:
            print(self.__class__.__name__ + " No BC")
            return
        self.bc_analysis()
        if self.dipole_len == None:
            self.dipole_len = copy.copy(self.dipoles[0].l)
        angle = np.arcsin(self.dipole_len / rho)
        ds = angle * rho
        if self.bc_gap == None:
            self.bc_gap = self.bc_gap_left*np.cos(self.dipoles[0].angle)
        drift = self.bc_gap / np.cos(angle)
        self.change_bc_shoulders(drift)
        # d.l=drift
        for i, dip in enumerate(self.dipoles):
            dip.angle = angle * np.sign(dip.angle)
            dip.l = ds
            if i in [0, 2]:
                dip.e2 = angle * np.sign(dip.angle)
            else:
                dip.e1 = angle * np.sign(dip.angle)
        self.lattice.update_transfer_maps()
        #for elem in self.lattice.sequence:
        #    print(elem.id, elem.l)


    def update_cavity(self, phi, v):
        for elem in self.lattice.sequence:
            if elem.__class__ == Cavity:
                if self.cav_name_pref == None:
                    elem.v = v
                    elem.phi = phi
                else:
                    if self.cav_name_pref in elem.id:
                        elem.v = v
                        elem.phi = phi

        self.lattice.update_transfer_maps()


    def init_navigator(self):

        # init navigator
        self.navigator = Navigator(self.lattice)
        self.navigator.unit_step = self.unit_step

        # init physics processes
        for physics_process in self.physics_processes_array:
            if (physics_process[0].__class__ == SpaceCharge or physics_process[0].__class__ == LSC) and self.sc_flag:
                self.navigator.add_physics_proc(physics_process[0], physics_process[1], physics_process[2])

            if physics_process[0].__class__ == CSR and self.csr_flag:
                self.navigator.add_physics_proc(physics_process[0], physics_process[1], physics_process[2])

            if (physics_process[0].__class__ == Wake or physics_process[0].__class__ == WakeKick) and self.wake_flag:
                self.navigator.add_physics_proc(physics_process[0], physics_process[1], physics_process[2])

            if physics_process[0].__class__ == BeamTransform and self.bt_flag:
                self.navigator.add_physics_proc(physics_process[0], physics_process[1], physics_process[2])

            if physics_process[0].__class__ == SmoothBeam and self.smooth_flag:
                self.navigator.add_physics_proc(physics_process[0], physics_process[1], physics_process[2])

            if physics_process[0].__class__ not in [SpaceCharge, CSR, Wake, WakeKick, BeamTransform, SmoothBeam, LSC]:
                #print(physics_process, " ADDED")
                self.navigator.add_physics_proc(physics_process[0], physics_process[1], physics_process[2])

    def add_physics_process(self, physics_process, start, stop):

        self.physics_processes_array.append([physics_process, start, stop])

    def read_beam_file(self):

        particles = None
        #print(self.input_beam_file)
        extension = self.input_beam_file.split(".")[-1]
        #print(extension)
        if extension == "ast":
        
            try:
                particles = astraBeam2particleArray(filename=self.input_beam_file)
            except:
                print(self.lattice_name + ' - #### ERROR #### - NO START PARTICLES FILE.')
        
        else:
            
            try:
                particles = load_particle_array(self.input_beam_file)
        
            except:
                print(self.lattice_name + ' - #### ERROR #### - NO START PARTICLES FILE.')

        return particles

    def save_beam_file(self, particles):

        save_particle_array(self.output_beam_file, particles)

    def save_twiss_file(self, twiss_list):
        if self.tws_file == None:
            tws_file_name = self.output_beam_file.replace("particles", "tws")
        else:
            tws_file_name = self.tws_file

        bx = np.array([tw.beta_x for tw in twiss_list])
        by = np.array([tw.beta_y for tw in twiss_list])
        ax = np.array([tw.alpha_x for tw in twiss_list])
        ay = np.array([tw.alpha_x for tw in twiss_list])
        s = np.array([tw.s for tw in twiss_list])
        E = np.array([tw.E for tw in twiss_list])

        emit_x = np.array([tw.emit_x for tw in twiss_list])
        emit_y = np.array([tw.emit_y for tw in twiss_list])

        np.savez_compressed(tws_file_name, beta_x=bx, beta_y=by, alpha_x=ax, alpha_y=ay, E=E, s=s,
                            emit_x=emit_x, emit_y=emit_y)

    def load_twiss_file(self):

        #with np.load(self.tws_file) as data:
        #    #for key in data.keys():
        #    #    p_array.__dict__[key] = data[key]
        return np.load(self.tws_file)

    def tracking(self, particles=None):

        # read beam file
        if particles == None:
            particles = self.read_beam_file()
            #print("Particle.s", particles.s, particles.q_array)
            if particles == None:
                return None

        # init navigator
        self.init_navigator()

        # tracking
        print()
        print(self.lattice_name + ' TRACKING')
        print("std1 = ", np.std(particles.tau()))
        tws_track, particles = track(self.lattice, particles, self.navigator,
                                     print_progress=self.print_progress, calc_tws=self.calc_tws)
        self.tws_track = tws_track
        # save tracking results
        if self.output_beam_file != None and not self.kill_track:
            self.save_beam_file(particles)
            self.save_twiss_file(tws_track)

        return particles

 

