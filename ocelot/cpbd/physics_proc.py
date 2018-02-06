from ocelot.cpbd.beam import save_particle_array, load_particle_array

class PhysProc:
    def __init__(self):
        self.step=1
        self.energy = None

    def prepare(self, lat):
        pass

    def apply(self, p_array, dz):
        pass


class EmptyProc(PhysProc):
    def __init__(self):
        PhysProc.__init__(self)
        self.step=1
        self.energy = None
        self.pict_debug = True
        self.traj_step = 0.0002

    def prepare(self, lat):
        pass

    def apply(self, p_array, dz):
        pass


class SaveBeam(PhysProc):
    def __init__(self, filename):
        PhysProc.__init__(self)
        self.step=1
        self.energy = None
        self.filename = filename

    def prepare(self, lat):
        pass

    def apply(self, p_array, dz):
        save_particle_array(filename=self.filename, p_array=p_array)