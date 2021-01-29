from ocelot.utils.section_track import *
import accelerator.lattice.i1 as i1
import accelerator.lattice.l1 as l1
from ocelot.cpbd.physics_proc import *
import os


Sig_Z = (0.0012761713706630296, 0.00043356707510439394, 6.229583056423405e-05, 5.132214278750506e-06)

SmoothPar=1000
LHE=1.4*10000e-9 # GeV


class A1(SectionTrack):
    def __init__(self, data_dir, *args, **kwargs):
        super().__init__(data_dir)
        # setting parameters
        self.lattice_name = 'A1'
        self.unit_step = 0.2
        self.input_beam_file = None # self.particle_dir + 'out.ast'
        self.output_beam_file = self.particle_dir + 'section_A1.npz'
        self.tws_file = self.tws_dir + "tws_section_A1.npz"
        # init tracking lattice
        start_sim = i1.start_sim
        acc1_stop = i1.a1_sim_stop
        if "coupler_kick" in kwargs:
            self.coupler_kick = kwargs["coupler_kick"]
        else:
            self.coupler_kick = True

        if "suffix" in kwargs:
            filename, file_extension = os.path.splitext(self.output_beam_file)
            self.output_beam_file = filename + str(kwargs["suffix"]) + file_extension
            filename, file_extension = os.path.splitext(self.tws_file)
            self.tws_file = filename + str(kwargs["suffix"]) + file_extension

        self.lattice = MagneticLattice(i1.cell, start=start_sim, stop=acc1_stop, method=self.method)
        # init physics processes
        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [63, 63, 63]
        sc2 = SpaceCharge()
        sc2.step = 1
        sc2.nmesh_xyz = [63, 63, 63]
        wake = Wake()
        wake.wake_table = WakeTable('./unit_tests/ebeam_test/section_track/accelerator/wakes/RF/mod_TESLA_MODULE_WAKE_TAYLOR.dat')
        wake.factor = 1
        wake.step = 50
        smooth = SmoothBeam()
        smooth.mslice = SmoothPar
        # adding physics processes
        acc1_1_stop = i1.a1_1_stop
        acc1_wake_kick = acc1_stop
        self.add_physics_process(smooth, start=start_sim, stop=start_sim)
        self.add_physics_process(sc, start=start_sim, stop=acc1_1_stop)
        self.add_physics_process(sc2, start=acc1_1_stop, stop=acc1_wake_kick)
        self.add_physics_process(wake, start=i1.c_a1_1_1_i1, stop=acc1_wake_kick)


class AH1(SectionTrack):
    def __init__(self, data_dir, *args, **kwargs):
        super().__init__(data_dir)
        # setting parameters
        self.lattice_name = 'Injector AH1'
        self.unit_step = 0.2
        self.input_beam_file = self.particle_dir + 'section_A1.npz'
        self.output_beam_file =  self.particle_dir + 'section_AH1.npz'
        self.tws_file = self.tws_dir + "tws_section_AH1.npz"
        # init tracking lattice
        acc1_stop = i1.a1_sim_stop
        acc39_stop = i1.stlat_47_i1

        self.lattice = MagneticLattice(i1.cell, start=acc1_stop, stop=acc39_stop, method=self.method)
        # init physics processes
        sc = SpaceCharge()
        sc.step = 5
        sc.nmesh_xyz = [63, 63, 63]
        wake = Wake()
        wake.wake_table = WakeTable('./unit_tests/ebeam_test/section_track/accelerator/wakes/RF/mod_THIRD_HARMONIC_SECTION_WAKE_TAYLOR.dat')
        wake.factor = 2
        wake.step = 50

        # adding physics processes
        match_acc39 = acc1_stop
        acc39_wake_kick = i1.stlat_47_i1
        self.add_physics_process(sc, start=match_acc39, stop=acc39_wake_kick)
        self.add_physics_process(wake, start=i1.c3_ah1_1_1_i1, stop=acc39_wake_kick)


class LH(SectionTrack):
    def __init__(self, data_dir, *args, **kwargs):
        super().__init__(data_dir)
        # setting parameters
        self.lattice_name = 'LASER HEATER MAGNETS'
        self.unit_step = 0.2
        self.input_beam_file = self.particle_dir + 'section_AH1.npz'
        self.output_beam_file = self.particle_dir + 'section_LH.npz'
        self.tws_file = self.tws_dir + "tws_section_LH.npz"

        if "suffix" in kwargs:
            filename, file_extension = os.path.splitext(self.input_beam_file)
            self.input_beam_file = filename + str(kwargs["suffix"]) + file_extension
            filename, file_extension = os.path.splitext(self.output_beam_file)
            self.output_beam_file = filename + str(kwargs["suffix"]) + file_extension
            filename, file_extension = os.path.splitext(self.tws_file)
            self.tws_file = filename + str(kwargs["suffix"]) + file_extension


        # init tracking lattice
        acc39_stop = i1.stlat_47_i1
        lhm_stop = l1.dl_start
        self.lattice = MagneticLattice(i1.cell + l1.cell, start=acc39_stop, stop=lhm_stop, method=self.method)
        # init physics processes
        sigma=Sig_Z[0]
        csr = CSR()
        csr.sigma_min = sigma * 0.1
        csr.traj_step = 0.0005
        csr.apply_step = 0.005
        sc = SpaceCharge()
        sc.step = 50
        sc.nmesh_xyz = [63, 63, 63]

        lh = LaserModulator()
        lh.Lu = 0.74
        lh.dE = LHE
        lh.sigma_l=300
        lh.sigma_x = 300e-6
        lh.sigma_y = 300e-6
        lh.z_waist = None
        self.add_physics_process(sc, start=acc39_stop, stop=lhm_stop)
        self.add_physics_process(csr, start=acc39_stop, stop=lhm_stop)
        self.add_physics_process(lh, start=i1.lh_start, stop=i1.lh_stop)


class DL(SectionTrack):
    def __init__(self, data_dir, *args, **kwargs):
        super().__init__(data_dir)
        # setting parameters
        self.lattice_name = 'DOGLEG'
        self.unit_step = 0.2
        self.input_beam_file = self.particle_dir + 'section_LH.npz'
        self.output_beam_file = self.particle_dir + 'section_DL.npz'
        self.tws_file = self.tws_dir + "tws_section_DL.npz"

        if "suffix" in kwargs:
            filename, file_extension = os.path.splitext(self.input_beam_file)
            self.input_beam_file = filename + str(kwargs["suffix"]) + file_extension
            filename, file_extension = os.path.splitext(self.output_beam_file)
            self.output_beam_file = filename + str(kwargs["suffix"]) + file_extension
            filename, file_extension = os.path.splitext(self.tws_file)
            self.tws_file = filename + str(kwargs["suffix"]) + file_extension

        # init tracking lattice
        st2_stop = l1.dl_start
        dogleg_stop = l1.stlat_96_i1
        self.lattice = MagneticLattice(l1.cell, start=st2_stop, stop=dogleg_stop, method=self.method)
        # init physics processes
        sigma=Sig_Z[0]
        csr = CSR()
        csr.sigma_min = sigma*0.1
        csr.traj_step = 0.0005
        csr.apply_step = 0.005

        sc = SpaceCharge()
        sc.step = 25
        sc.nmesh_xyz = [63, 63, 63]
        self.add_physics_process(csr, start=st2_stop, stop=dogleg_stop)
        self.add_physics_process(sc, start=st2_stop, stop=dogleg_stop)


class BC0(SectionTrack):

    def __init__(self, data_dir, *args, **kwargs):
        super().__init__(data_dir)

        # setting parameters
        self.lattice_name = 'BC0'
        self.unit_step = 0.5


        self.input_beam_file = self.particle_dir + 'section_DL.npz'
        self.output_beam_file = self.particle_dir + 'section_BC0.npz'
        self.tws_file = self.tws_dir + "tws_section_BC0.npz"


        if "suffix" in kwargs:
            filename, file_extension = os.path.splitext(self.input_beam_file)
            self.input_beam_file = filename + str(kwargs["suffix"]) + file_extension
            filename, file_extension = os.path.splitext(self.output_beam_file)
            self.output_beam_file = filename + str(kwargs["suffix"]) + file_extension
            filename, file_extension = os.path.splitext(self.tws_file)
            self.tws_file = filename + str(kwargs["suffix"]) + file_extension

        # init tracking lattice
        st4_stop = l1.stlat_96_i1
        bc0_stop = l1.enlat_101_i1
        self.lattice = MagneticLattice(l1.cell, start=st4_stop, stop=bc0_stop, method=self.method)

        # init physics processes

        sigma=Sig_Z[0]
        csr = CSR()
        csr.sigma_min = sigma*0.1
        csr.traj_step = 0.0005
        csr.apply_step = 0.005

        sc = SpaceCharge()
        sc.step = 10
        sc.nmesh_xyz = [63, 63, 63]
        sc.low_order_kick = False
        match_bc0 = st4_stop
        self.add_physics_process(sc, start=match_bc0, stop=bc0_stop)
        self.add_physics_process(csr, start=match_bc0, stop=bc0_stop)
        self.dipoles = [l1.bb_96_i1, l1.bb_98_i1, l1.bb_100_i1, l1.bb_101_i1]
        self.dipole_len = 0.5
        self.bc_gap=1.0


class L1(SectionTrack):
    
    def __init__(self, data_dir, *args, **kwargs):
        super().__init__(data_dir)

        # setting parameters
        self.lattice_name = 'L1'
        self.unit_step = 0.2


        self.input_beam_file = self.particle_dir + 'section_BC0.npz'
        self.output_beam_file = self.particle_dir + 'section_L1.npz'
        self.tws_file = self.tws_dir + "tws_section_L1.npz"

        if "suffix" in kwargs:
            filename, file_extension = os.path.splitext(self.input_beam_file)
            suff = str(kwargs["suffix"])
            indx = suff.find("_chirpL1_")
            input_suff = suff[:indx]
            self.input_beam_file = filename + input_suff + file_extension
            print("SECTION L1: ", self.input_beam_file)
            filename, file_extension = os.path.splitext(self.output_beam_file)
            self.output_beam_file = filename + str(kwargs["suffix"]) + file_extension
            filename, file_extension = os.path.splitext(self.tws_file)
            self.tws_file = filename + str(kwargs["suffix"]) + file_extension


        bc0_stop = l1.enlat_101_i1
        acc2_stop = l1.stlat_182_b1

        if "coupler_kick" in kwargs:
            self.coupler_kick = kwargs["coupler_kick"]
        else:
            self.coupler_kick = True

        # init tracking lattice
        self.lattice = MagneticLattice(l1.cell, start=bc0_stop, stop=acc2_stop, method=self.method)

        # init physics processes
        smooth = SmoothBeam()
        smooth.mslice = SmoothPar

        sc = SpaceCharge()
        sc.step = 50
        sc.nmesh_xyz = [31, 31, 31]
        wake = Wake()
        wake.wake_table = WakeTable('./unit_tests/ebeam_test/section_track/accelerator/wakes/RF/mod_TESLA_MODULE_WAKE_TAYLOR.dat')
        wake.factor = 4
        wake.step = 100

        match_acc2 = bc0_stop
        L1_wake_kick = acc2_stop
        self.add_physics_process(smooth, start=match_acc2, stop=match_acc2)
        self.add_physics_process(sc, start=match_acc2, stop=L1_wake_kick)
        self.add_physics_process(wake, start=l1.c_a2_1_1_l1, stop=l1.c_a2_4_8_l1)


class BC1(SectionTrack):

    def __init__(self, data_dir, *args, **kwargs):
        super().__init__(data_dir)

        # setting parameters
        self.lattice_name = 'BC1'
        self.unit_step = 0.2

        self.input_beam_file = self.particle_dir + 'section_L1.npz'
        self.output_beam_file = self.particle_dir + 'section_BC1.npz'
        self.tws_file = self.tws_dir + "tws_section_BC1.npz"

        if "suffix" in kwargs:
            filename, file_extension = os.path.splitext(self.input_beam_file)
            self.input_beam_file = filename + str(kwargs["suffix"]) + file_extension
            print("SECTION B1: ", self.input_beam_file)
            filename, file_extension = os.path.splitext(self.output_beam_file)
            self.output_beam_file = filename + str(kwargs["suffix"]) + file_extension
            filename, file_extension = os.path.splitext(self.tws_file)
            self.tws_file = filename + str(kwargs["suffix"]) + file_extension

        acc2_stop = l1.stlat_182_b1
        bc1_stop = l1.tora_203_b1
        # init tracking lattice
        self.lattice = MagneticLattice(l1.cell, start=acc2_stop, stop=bc1_stop, method=self.method)

        # init physics processes

        sigma = Sig_Z[1]
        csr = CSR()
        csr.sigma_min = sigma*0.1
        csr.traj_step = 0.0005
        csr.apply_step = 0.005

        sc = SpaceCharge()
        sc.step = 20
        sc.nmesh_xyz = [31, 31, 31]
        match_bc1 = acc2_stop
        self.add_physics_process(csr, start=match_bc1, stop=bc1_stop)
        self.add_physics_process(sc, start=match_bc1, stop=bc1_stop)
        self.dipoles = [l1.bb_182_b1, l1.bb_191_b1, l1.bb_193_b1, l1.bb_202_b1]
        self.dipole_len = 0.5
        self.bc_gap=8.5
