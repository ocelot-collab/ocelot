from ocelot import *
import numpy as np

class SectionTrack:
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):

        self.lattice_name = ""
        self.lattice = None
        self.unit_step = 1.0

        self.input_beam_file = None
        self.output_beam_file = None
        self.tws_file = None

        self.physics_processes_array = []  # list of physics process

        self.method = MethodTM()
        self.method.global_method = SecondTM
        self.translator = {SpaceCharge: "sc_flag", CSR:"csr_flag",
                           WakeKick:"wake_flag", BeamTransform:"bt_flag"}
        self.sc_flag = sc_flag
        self.csr_flag = csr_flag
        self.wake_flag = wake_flag
        self.bt_flag = True

        self.print_progress = True
        self.calc_tws = True
        self.kill_track = False

    def init_navigator(self):

        # init navigator
        self.navigator = Navigator(self.lattice)
        self.navigator.unit_step = self.unit_step

        # init physics processes
        for physics_process in self.physics_processes_array:

            if physics_process[0].__class__ == SpaceCharge and self.sc_flag:
                self.navigator.add_physics_proc(physics_process[0], physics_process[1], physics_process[2])

            if physics_process[0].__class__ == CSR and self.csr_flag:
                self.navigator.add_physics_proc(physics_process[0], physics_process[1], physics_process[2])

            if (physics_process[0].__class__ == Wake or physics_process[0].__class__ == WakeKick) and self.wake_flag:
                self.navigator.add_physics_proc(physics_process[0], physics_process[1], physics_process[2])

            if physics_process[0].__class__ == BeamTransform and self.bt_flag:
                self.navigator.add_physics_proc(physics_process[0], physics_process[1], physics_process[2])

            if physics_process[0].__class__ not in [SpaceCharge, CSR, Wake, WakeKick, BeamTransform]:
                self.navigator.add_physics_proc(physics_process[0], physics_process[1], physics_process[2])

    def add_physics_process(self, physics_process, start, stop):

        self.physics_processes_array.append([physics_process, start, stop])

    def read_beam_file(self):

        particles = None
        #print(self.input_beam_file)
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

        np.savez_compressed(tws_file_name, beta_x=bx, beta_y=by, alpha_x=ax, alpha_y=ay, E=E, s=s)

    def load_twiss_file(self):

        #with np.load(self.tws_file) as data:
        #    #for key in data.keys():
        #    #    p_array.__dict__[key] = data[key]
        return np.load(self.tws_file)

    def tracking(self, particles=None):

        # read beam file
        if particles == None:
            particles = self.read_beam_file()

            if particles == None:
                return None

        # init navigator
        self.init_navigator()

        # tracking
        print(self.lattice_name + ' TRACKING')
        tws_track, particles = track(self.lattice, particles, self.navigator,
                                     print_progress=self.print_progress, calc_tws=self.calc_tws)

        # save tracking results
        if self.output_beam_file != None and not self.kill_track:
            self.save_beam_file(particles)
            self.save_twiss_file(tws_track)

        return particles


