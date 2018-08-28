from ocelot import *
import accelerator.lattice.muy_5pi_mode_A.i1_full as i1
import accelerator.lattice.muy_5pi_mode_A.l1_full as l1
import accelerator.lattice.muy_5pi_mode_A.l2_full as l2
import accelerator.lattice.muy_5pi_mode_A.l3_full as l3
import accelerator.lattice.muy_5pi_mode_A.cl_full as cl
import accelerator.lattice.muy_5pi_mode_A.tl34 as tl34
import accelerator.lattice.muy_5pi_mode_A.sase1 as sase1
import accelerator.lattice.muy_5pi_mode_A.t4 as t4
import accelerator.lattice.muy_5pi_mode_A.sase3 as sase3
import accelerator.lattice.muy_5pi_mode_A.t4d as t4d
import numpy as np
from accelerator.s2e_sections.section_track import *
from ocelot.cpbd.physics_proc import *



class A1(SectionTrack):

    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'A1'
        self.unit_step = 0.02

        #self.input_beam_file = self.particle_dir + 'Exfel.0320_200k.ast'
        self.input_beam_file = self.particle_dir + 'Exfel.0320_250pC.ast'
        #self.input_beam_file = self.particle_dir + 'XFEL_500pC_setup_D_03p2_1M.ast'
        self.output_beam_file = self.particle_dir + 'section_A1.npz'
        self.tws_file = self.tws_dir + "tws_section_A1.npz"

        # init tracking lattice
        start_sim = i1.id_22433449_
        acc1_stop = i1.id_68749308_

        self.lattice = MagneticLattice(i1.cell, start=start_sim, stop=acc1_stop, method=self.method)

        # init physics processes
        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [63, 63, 63]

        sc2 = SpaceCharge()
        sc2.step = 5
        sc2.nmesh_xyz = [63, 63, 63]

        wake = WakeKick()
        wake.wake_table = WakeTable('accelerator/wakes/mod_TESLA_MODULE_WAKE_TAYLOR.dat')
        wake.factor = 1

        smooth = SmoothBeam()
        smooth.mslice = 10000
        # adding physics processes
        acc1_1_stop = i1.id_75115473_
        acc1_wake_kick = acc1_stop


        self.add_physics_process(smooth, start=start_sim, stop=start_sim)
        self.add_physics_process(sc, start=start_sim, stop=acc1_1_stop)
        self.add_physics_process(sc2, start=acc1_1_stop, stop=acc1_wake_kick)
        self.add_physics_process(wake, start=acc1_wake_kick, stop=acc1_wake_kick)


class AH1(SectionTrack):

    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'AH1'
        self.unit_step = 0.02


        self.input_beam_file = self.particle_dir + 'section_A1.npz'
        self.output_beam_file = self.particle_dir + 'section_AH1.npz'
        self.tws_file = self.tws_dir + "tws_section_AH1.npz"
        # init tracking lattice
        acc1_stop = i1.id_68749308_
        acc39_stop = i1.stlat_47_i1

        self.lattice = MagneticLattice(i1.cell, start=acc1_stop, stop=acc39_stop, method=self.method)

        # init physics processes
        sc = SpaceCharge()
        sc.step = 5
        sc.nmesh_xyz = [63, 63, 63]

        wake = WakeKick()
        wake.wake_table = WakeTable('accelerator/wakes/mod_THIRD_HARMONIC_SECTION_WAKE_TAYLOR.dat')
        wake.factor = 2

        # adding physics processes
        match_acc39 = acc1_stop
        acc39_wake_kick = i1.stlat_47_i1
        smooth = SmoothBeam()
        smooth.mslice = 10000
        #self.add_physics_process(smooth, start=acc39_wake_kick, stop=acc39_wake_kick)
        self.add_physics_process(sc, start=match_acc39, stop=acc39_wake_kick)
        self.add_physics_process(wake, start=acc39_wake_kick, stop=acc39_wake_kick)


class LH(SectionTrack):

    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'LH'
        self.unit_step = 0.02

        self.input_beam_file = self.particle_dir + 'section_AH1.npz'
        self.output_beam_file = self.particle_dir + 'section_LH.npz'
        self.tws_file = self.tws_dir + "tws_section_LH.npz"
        # init tracking lattice
        acc39_stop = i1.stlat_47_i1
        lhm_stop = l1.id_90904668_

        self.lattice = MagneticLattice(i1.cell + l1.cell, start=acc39_stop, stop=lhm_stop, method=self.method)

        # init physics processes
        sigma=0.001852509071807481

        csr = CSR()
        csr.sigma_min = sigma * 0.1
        csr.traj_step = 0.0005
        csr.apply_step = 0.005

        sc = SpaceCharge()
        sc.step = 50
        sc.nmesh_xyz = [63, 63, 63]
        tws = Twiss()
        tws.beta_x = 2.36238404123
        tws.beta_y = 2.90712039319
        tws.alpha_x = 1.23079453323
        tws.alpha_y = -1.45354874634
        tws.E = 0.13

        bt = BeamTransform(tws=tws)
        self.add_physics_process(bt, start=i1.match_55_i1, stop=i1.match_55_i1)
        self.add_physics_process(sc, start=acc39_stop, stop=lhm_stop)
        self.add_physics_process(csr, start=acc39_stop, stop=lhm_stop)




class DL(SectionTrack):
    
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)
        
        # setting parameters
        self.lattice_name = 'DL'
        self.unit_step = 0.02

        self.input_beam_file = self.particle_dir + 'section_LH.npz'
        self.output_beam_file = self.particle_dir + 'section_DL.npz'
        self.tws_file = self.tws_dir + "tws_section_DL.npz"

        # init tracking lattice
        st2_stop = l1.id_90904668_
        dogleg_stop = l1.stlat_96_i1

        self.lattice = MagneticLattice(l1.cell, start=st2_stop, stop=dogleg_stop, method=self.method)

        # init physics processes
        sigma=0.001852509071807481
        csr = CSR()
        #csr.step=10
        #csr.n_bin = 100
        csr.sigma_min = sigma*0.1
        csr.traj_step = 0.0005
        csr.apply_step = 0.005

        sc = SpaceCharge()
        sc.step = 25
        sc.nmesh_xyz = [31, 31, 31]
        self.add_physics_process(csr, start=st2_stop, stop=dogleg_stop)
        self.add_physics_process(sc, start=st2_stop, stop=dogleg_stop)


class BC0(SectionTrack):

    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'BC0'
        self.unit_step = 0.1


        self.input_beam_file = self.particle_dir + 'section_DL.npz'
        self.output_beam_file = self.particle_dir + 'section_BC0.npz'
        self.tws_file = self.tws_dir + "tws_section_BC0.npz"
        # init tracking lattice
        st4_stop = l1.stlat_96_i1
        bc0_stop = l1.enlat_101_i1
        self.lattice = MagneticLattice(l1.cell, start=st4_stop, stop=bc0_stop, method=self.method)

        # init physics processes

        sigma=0.0015288820707131323
        csr = CSR()
        #csr.step=10
        #csr.n_bin = 100
        csr.sigma_min = sigma*0.1
        csr.traj_step = 0.0005
        csr.apply_step = 0.005

        smooth = SmoothBeam()
        smooth.mslice = 10000

        #sc = SpaceCharge()
        #sc.step = 10
        #sc.nmesh_xyz = [63, 63, 63]
        #sc.low_order_kick = False
        match_bc0 = st4_stop
        self.add_physics_process(smooth, start=bc0_stop, stop=bc0_stop)
        self.add_physics_process(csr, start=match_bc0, stop=bc0_stop)

        self.dipoles = [l1.bb_96_i1, l1.bb_98_i1, l1.bb_100_i1, l1.bb_101_i1]

class L1(SectionTrack):
    
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'L1'
        self.unit_step = 0.02


        self.input_beam_file = self.particle_dir + 'section_BC0.npz'
        self.output_beam_file = self.particle_dir + 'section_L1.npz'
        self.tws_file = self.tws_dir + "tws_section_L1.npz"
        bc0_stop = l1.enlat_101_i1
        acc2_stop = l1.stlat_182_b1
        # init tracking lattice
        self.lattice = MagneticLattice(l1.cell, start=bc0_stop, stop=acc2_stop, method=self.method)

        # init physics processes

        sc = SpaceCharge()
        sc.step = 50
        sc.nmesh_xyz = [31, 31, 31]
        wake = WakeKick()
        wake.wake_table = WakeTable('accelerator/wakes/mod_TESLA_MODULE_WAKE_TAYLOR.dat')
        wake.factor = 4
        match_acc2 = bc0_stop
        L1_wake_kick = acc2_stop
        self.add_physics_process(sc, start=match_acc2, stop=L1_wake_kick)
        self.add_physics_process(wake, start=L1_wake_kick, stop=L1_wake_kick)


class BC1(SectionTrack):

    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):

        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'BC1'
        self.unit_step = 0.05

        self.input_beam_file = self.particle_dir + 'section_L1.npz'
        self.output_beam_file = self.particle_dir + 'section_BC1.npz'
        self.tws_file = self.tws_dir + "tws_section_BC1.npz"

        acc2_stop = l1.stlat_182_b1
        bc1_stop = l1.tora_203_b1
        # init tracking lattice
        self.lattice = MagneticLattice(l1.cell, start=acc2_stop, stop=bc1_stop, method=self.method)

        # init physics processes

        smooth = SmoothBeam()
        smooth.mslice = 10000

        sigma = 0.0005748235512972583
        csr = CSR()
        #csr.step = 10
        #csr.n_bin = 100
        csr.sigma_min = sigma*0.1
        csr.traj_step = 0.0005
        csr.apply_step = 0.005

        sc = SpaceCharge()
        sc.step = 20
        sc.nmesh_xyz = [31, 31, 31]
        match_bc1 = acc2_stop
        self.add_physics_process(smooth, start=bc1_stop, stop=bc1_stop)
        self.add_physics_process(csr, start=match_bc1, stop=bc1_stop)
        self.add_physics_process(sc, start=match_bc1, stop=bc1_stop)
        self.dipoles = [l1.bb_182_b1, l1.bb_191_b1, l1.bb_193_b1, l1.bb_202_b1]



class L2(SectionTrack):

    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'L2'
        self.unit_step = 0.02

        self.input_beam_file = self.particle_dir + 'section_BC1.npz'
        self.output_beam_file = self.particle_dir + 'section_L2.npz'
        self.tws_file = self.tws_dir + "tws_section_L2.npz"

        bc1_stop = l1.tora_203_b1
        acc3t5_stop = l2.stlat_393_b2
        # init tracking lattice
        self.lattice = MagneticLattice(l1.cell + l2.cell, start=bc1_stop, stop=acc3t5_stop, method=self.method)

        # init physics processes

        sc = SpaceCharge()
        sc.step = 100
        sc.nmesh_xyz = [31, 31, 31]

        wake = WakeKick()
        wake.wake_table = WakeTable('accelerator/wakes/mod_TESLA_MODULE_WAKE_TAYLOR.dat')
        wake.factor = 4 * 3

        self.add_physics_process(sc, start=bc1_stop, stop=acc3t5_stop)
        self.add_physics_process(wake, start=acc3t5_stop, stop=acc3t5_stop)


class BC2(SectionTrack):
    
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'BC2'
        self.dipoles = [l2.bb_393_b2, l2.bb_402_b2, l2.bb_404_b2, l2.bb_413_b2]

        self.unit_step = 0.02

        self.input_beam_file = self.particle_dir + 'section_L2.npz'
        self.output_beam_file = self.particle_dir + 'section_BC2.npz'
        self.tws_file = self.tws_dir + "tws_section_BC2.npz"

        acc3t5_stop = l2.stlat_393_b2
        bc2_stop = l2.tora_415_b2
        # init tracking lattice
        self.lattice = MagneticLattice(l2.cell, start=acc3t5_stop, stop=bc2_stop, method=self.method)

        # init physics processes

        csr = CSR()
        csr.step=2
        csr.n_bin = 100
        csr.sigma_min = 7.e-6
        csr.traj_step = 0.0005
        csr.apply_step = 0.005

        smooth = SmoothBeam()
        smooth.mslice = 10000

        sc = SpaceCharge()
        sc.step = 50
        sc.nmesh_xyz = [31, 31, 31]

        self.add_physics_process(smooth, start=bc2_stop, stop=bc2_stop)
        self.add_physics_process(csr, start=acc3t5_stop, stop=bc2_stop)
        self.add_physics_process(sc, start=acc3t5_stop, stop=bc2_stop)


class L3(SectionTrack):
    
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'L3'
        self.unit_step = 5.0


        self.input_beam_file = self.particle_dir + 'section_BC2.npz'
        self.output_beam_file = self.particle_dir + 'section_L3.npz'
        self.tws_file = self.tws_dir + "tws_section_L3.npz"

        bc2_stop = l2.tora_415_b2
        acc6t26_stop = cl.match_1673_cl
        # init tracking lattice
        self.lattice = MagneticLattice(l2.cell + l3.cell + cl.cell, start=bc2_stop, stop=acc6t26_stop, method=self.method)
        #self.cav_list = [".A23.", ".A22.", ".A21.", ".A20.", ".A19.", ".A18.", ".A17.", ".A16."]
        #for elem in self.lattice.sequence:
        #    if elem.__class__ == Cavity:
        #        for cav_num in self.cav_list:
        #            if cav_num in elem.id:
        #                print(elem.id)
        #                elem.v = 0.
        #self.lattice.update_transfer_maps()

        # init physics processes

        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [31, 31, 31]

        wake = WakeKick()
        wake.wake_table = WakeTable('accelerator/wakes/mod_TESLA_MODULE_WAKE_TAYLOR.dat')
        wake.factor = 4 * 21

        self.add_physics_process(wake, start=acc6t26_stop, stop=acc6t26_stop)
        self.add_physics_process(sc, start=bc2_stop, stop=acc6t26_stop)


class CL1(SectionTrack):
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'CL1'
        self.unit_step = 0.2


        self.input_beam_file = self.particle_dir + 'section_L3.npz'
        self.output_beam_file = self.particle_dir + 'section_CL1.npz'
        self.tws_file = self.tws_dir + "tws_section_CL1.npz"

        acc6t26_stop = cl.match_1673_cl
        collimator1_stop = cl.bpma_1746_cl
        # init tracking lattice
        self.lattice = MagneticLattice(cl.cell, start=acc6t26_stop, stop=collimator1_stop, method=self.method)

        # init physics processes

        sc = SpaceCharge()
        sc.step = 10
        sc.nmesh_xyz = [31, 31, 31]
        sc.low_order_kick = False

        csr = CSR()
        csr.traj_step = 0.0005
        csr.apply_step = 0.005
        csr.n_bin = 300
        csr.sigma_min = 0.2e-6

        self.add_physics_process(csr, start=acc6t26_stop, stop=collimator1_stop)
        self.add_physics_process(sc, start=acc6t26_stop, stop=collimator1_stop)


class CL2(SectionTrack):
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'CL2'
        self.unit_step = 1


        self.input_beam_file = self.particle_dir + 'section_CL1.npz'
        self.output_beam_file = self.particle_dir + 'section_CL2.npz'
        self.tws_file = self.tws_dir + "tws_section_CL2.npz"

        collimator1_stop = cl.bpma_1746_cl
        collimator2_stop = cl.bpma_1783_cl
        # init tracking lattice
        self.lattice = MagneticLattice(cl.cell, start=collimator1_stop, stop=collimator2_stop, method=self.method)

        # init physics processes

        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [31, 31, 31]
        sc.low_order_kick = False
        self.add_physics_process(sc, start=collimator1_stop, stop=collimator2_stop)


class CL3(SectionTrack):
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'CL3'
        self.unit_step = 0.2


        self.input_beam_file = self.particle_dir + 'section_CL2.npz'
        self.output_beam_file = self.particle_dir + 'section_CL3.npz'
        self.tws_file = self.tws_dir + "tws_section_CL3.npz"

        collimator2_stop = cl.bpma_1783_cl
        collimator3_stop = cl.ensec_1854_cl
        # init tracking lattice
        self.lattice = MagneticLattice(cl.cell, start=collimator2_stop, stop=collimator3_stop, method=self.method)

        # init physics processes

        sc = SpaceCharge()
        sc.step = 10
        sc.nmesh_xyz = [31, 31, 31]
        sc.low_order_kick = False

        csr = CSR()
        csr.traj_step = 0.0005
        csr.apply_step = 0.005
        csr.n_bin = 300
        csr.sigma_min = 0.2e-6

        self.add_physics_process(csr, start=collimator2_stop, stop=collimator3_stop)
        self.add_physics_process(sc, start=collimator2_stop, stop=collimator3_stop)


class STN10(SectionTrack):
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'ST10'
        self.unit_step = 10

        self.input_beam_file = self.particle_dir + 'section_CL3.npz'
        self.output_beam_file = self.particle_dir + 'section_STN10.npz'
        self.tws_file = self.tws_dir + "tws_section_STN10.npz"

        collimator3_stop = cl.ensec_1854_cl
        stN10_stop = sase1.ensec_2235_t2
        # init tracking lattice
        self.lattice = MagneticLattice(cl.cell + tl34.cell + sase1.cell, start=collimator3_stop, stop=stN10_stop, method=self.method)

        # init physics processes

        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [31, 31, 31]
        sc.low_order_kick = False

        self.add_physics_process(sc, start=collimator3_stop, stop=stN10_stop)


class SASE1(SectionTrack):
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)
        # setting parameters
        self.lattice_name = 'SASE1'
        self.unit_step = 5

        self.input_beam_file = self.particle_dir + 'section_STN10.npz'
        self.output_beam_file = self.particle_dir + 'section_SASE1.npz'
        self.tws_file = self.tws_dir + "tws_section_SASE1.npz"
        # last element sase1 - stsec_2461_t4
        stN10_stop = sase1.ensec_2235_t2
        sase1_stop = sase1.stsec_2461_t4
        # init tracking lattice
        self.lattice = MagneticLattice(sase1.cell, start=stN10_stop, stop=sase1_stop, method=self.method)

        # init physics processes
        wake = Wake()
        wake.wake_table = WakeTable('accelerator/wakes/wake_undulator_OCELOT.txt')
        wake.step = 15
        wake.w_sampling = 500
        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [31, 31, 31]
        sc.low_order_kick = False
        self.add_physics_process(wake, start=sase1.match_2247_sa1, stop=sase1_stop)
        self.add_physics_process(sc, start=stN10_stop, stop=sase1_stop)


class T4(SectionTrack):
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)
        # setting parameters
        self.lattice_name = 'T4'
        self.unit_step = 0.2

        self.input_beam_file = self.particle_dir + 'section_SASE1.npz'
        self.output_beam_file = self.particle_dir + 'section_T4.npz'
        self.tws_file = self.tws_dir + "tws_section_T4.npz"
        # last element sase1 - stsec_2461_t4
        sase1_stop = sase1.stsec_2461_t4
        t4_stop = t4.ensub_2800_t4
        csr_start = t4.t4_start_csr
        csr_stop = t4.bpma_2606_t4
        # init tracking lattice
        self.lattice = MagneticLattice(sase1.cell + t4.cell, start=sase1_stop, stop=t4_stop, method=self.method)

        # init physics processes

        sc = SpaceCharge()
        sc.step = 25
        sc.nmesh_xyz = [31, 31, 31]

        csr = CSR()
        csr.traj_step = 0.0005
        csr.apply_step = 0.005
        csr.n_bin = 300
        csr.sigma_min = 0.2e-6

        sc2 = SpaceCharge()
        sc2.step = 25
        sc2.nmesh_xyz = [31, 31, 31]


        # creation of wake object with parameters
        wake = Wake()
        wake.wake_table = WakeTable('accelerator/wakes/wake_hor_1m_500mkm.txt')
        #wake.wake_table = WakeTable('accelerator/wakes/wake_hor_1m.txt')
        # w_sampling - defines the number of the equidistant sampling points for the one-dimensional
        # wake coefficients in the Taylor expansion of the 3D wake function.
        wake.w_sampling = 500
        wake.factor = 1
        wake.step = 1  # step in Navigator.unit_step, dz = Navigator.unit_step * wake.step [m]


        # creation of wake object with parameters
        wake_vert = Wake()
        wake_vert.factor = 1
        wake_vert.wake_table = WakeTable('accelerator/wakes/wake_vert_1m_500mkm.txt')
        #wake_vert.wake_table = WakeTable('accelerator/wakes/wake_vert_1m.txt')
        wake_vert.w_sampling = 500
        wake_vert.step = 1  # step in Navigator.unit_step, dz = Navigator.unit_step * wake.step [m]




        self.add_physics_process(wake, start=t4.wake_start, stop=t4.m_tds)

        svb4 = SaveBeam(filename=self.particle_dir + "before_structure.npz")
        self.add_physics_process(svb4, start=t4.wake_start, stop=t4.wake_start)

        self.add_physics_process(wake_vert, start=t4.m_tds, stop=t4.wake_stop)


        svb3 = SaveBeam(filename=self.particle_dir + "after_structure.npz")
        self.add_physics_process(svb3, start=t4.wake_stop, stop=t4.wake_stop)

        svb1 = SaveBeam(filename=self.particle_dir + "screen1.npz")
        self.add_physics_process(svb1, start=t4.m_img1, stop=t4.m_img1)

        svb2 = SaveBeam(filename=self.particle_dir + "screen2.npz")
        self.add_physics_process(svb2, start=t4.m_img2, stop=t4.m_img2)

        #self.add_physics_process(sc, start=sase1_stop, stop=csr_start)
        self.add_physics_process(sc, start=sase1_stop, stop=csr_start)
        self.add_physics_process(csr, start=csr_start, stop=csr_stop)
        self.add_physics_process(sc2, start=csr_stop, stop=t4.ensub_2800_t4)

        sc_in_bend = SpaceCharge()
        sc_in_bend.step = 25
        sc_in_bend.nmesh_xyz = [31, 31, 31]
        #self.add_physics_process(sc_in_bend, start=csr_start, stop=csr_stop)


class SASE3(SectionTrack):
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)
        # setting parameters
        self.lattice_name = 'SASE3'
        self.unit_step = 1

        self.input_beam_file = self.particle_dir + 'section_T4.npz'
        self.output_beam_file = self.particle_dir + 'section_SASE3.npz'
        self.tws_file = self.tws_dir + "tws_section_SASE3.npz"

        start = sase3.ensec_2800_t4
        stop = sase3.ensec_2940_sa3
        # init tracking lattice
        self.lattice = MagneticLattice(sase3.cell, start=start, stop=stop, method=self.method)

        # init physics processes

        sc = SpaceCharge()
        sc.step = 10
        sc.nmesh_xyz = [31, 31, 31]
        self.add_physics_process(sc, start=start, stop=stop)


class T4D(SectionTrack):
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)
        # setting parameters
        self.lattice_name = 'SASE3'
        self.unit_step = 1

        self.input_beam_file = self.particle_dir + 'section_SASE3.npz'
        self.output_beam_file = self.particle_dir + 'section_T4D.npz'
        self.tws_file = self.tws_dir + "tws_section_tT4D.npz"

        start = t4d.stsec_2940_t4d
        stop = t4d.ensec_3106_t4d
        # init tracking lattice
        self.lattice = MagneticLattice(t4d.cell, start=start, stop=stop, method=self.method)

        # init physics processes

        sc = SpaceCharge()
        sc.step = 10
        sc.nmesh_xyz = [31, 31, 31]
        self.add_physics_process(sc, start=start, stop=stop)

        csr = CSR()
        csr.traj_step = 0.0005
        csr.apply_step = 0.005
        csr.n_bin = 300
        csr.sigma_min = 0.2e-6

        self.add_physics_process(csr, start=t4d.tora_3065_t4d, stop=t4d.qk_3090_t4d)


class T4_short(SectionTrack):
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)
        # setting parameters
        self.lattice_name = 'T4'
        self.unit_step = 1
        self.calc_tws = False

        self.input_beam_file = self.particle_dir + 'before_structure.npz'
        self.output_beam_file = self.particle_dir + 'section_T4.npz'
        self.tws_file = self.tws_dir + "tws_section_T4.npz"
        # last element sase1 - stsec_2461_t4
        sase1_stop = t4.wake_start# sase1.stsec_2461_t4
        t4_stop = t4.m_img1 #t4.ensub_2800_t4
        csr_start = t4.t4_start_csr
        csr_stop = t4.bpma_2606_t4
        # init tracking lattice
        self.lattice = MagneticLattice(sase1.cell + t4.cell, start=sase1_stop, stop=t4_stop, method=self.method)

        # init physics processes

        sc = SpaceCharge()
        sc.step = 25
        sc.nmesh_xyz = [31, 31, 31]

        csr = CSR()
        csr.traj_step = 0.0005
        csr.apply_step = 0.005
        csr.n_bin = 300
        csr.sigma_min = 0.2e-6

        sc2 = SpaceCharge()
        sc2.step = 25
        sc2.nmesh_xyz = [31, 31, 31]


        # creation of wake object with parameters
        wake = Wake()
        wake.wake_table = WakeTable('accelerator/wakes/wake_hor_1m_500mkm.txt')

        # w_sampling - defines the number of the equidistant sampling points for the one-dimensional
        # wake coefficients in the Taylor expansion of the 3D wake function.
        wake.w_sampling = 500
        wake.factor = 1
        wake.step = 1  # step in Navigator.unit_step, dz = Navigator.unit_step * wake.step [m]


        # creation of wake object with parameters
        wake_vert = Wake()
        wake_vert.factor = 1
        wake_vert.wake_table = WakeTable('accelerator/wakes/wake_vert_1m_500mkm.txt')
        wake_vert.w_sampling = 500
        wake_vert.step = 1  # step in Navigator.unit_step, dz = Navigator.unit_step * wake.step [m]


        #svb1 = SaveBeam(filename=self.particle_dir + "screen1.npz")

        #self.add_physics_process(svb1, start=t4.m_img1, stop=t4.m_img1)

        #svb2 = SaveBeam(filename=self.particle_dir + "screen2.npz")
        svb3 = SaveBeam(filename=self.particle_dir + "after_structure.npz")
        #svb4 = SaveBeam(filename=self.particle_dir + "before_structure.npz")
        #self.add_physics_process(svb2, start=t4.m_img2, stop=t4.m_img2)

        self.add_physics_process(wake_vert, start=t4.wake_start, stop=t4.m_tds)
        #self.add_physics_process(wake, start=t4.m_tds, stop=t4.wake_stop)

        #self.add_physics_process(svb3, start=t4.wake_stop, stop=t4.wake_stop)
        #self.add_physics_process(svb4, start=t4.wake_start, stop=t4.wake_start)

        #self.add_physics_process(sc, start=sase1_stop, stop=csr_start)
        #self.add_physics_process(sc, start=sase1_stop, stop=csr_start)
        #self.add_physics_process(csr, start=csr_start, stop=csr_stop)
        #self.add_physics_process(sc2, start=csr_stop, stop=t4.ensub_2800_t4)

        sc_in_bend = SpaceCharge()
        sc_in_bend.step = 25
        sc_in_bend.nmesh_xyz = [31, 31, 31]
        #self.add_physics_process(sc_in_bend, start=csr_start, stop=csr_stop)



sections = [A1, AH1, LH, DL,  BC0, L1, BC1, L2, BC2, L3, CL1,
            CL2, CL3, STN10]