from ocelot import *
from accelerator.lattice.xfel_i1_mad import *
from accelerator.lattice.xfel_l1_mad import *
from accelerator.lattice.xfel_l2_mad import *
from accelerator.lattice.xfel_l3_no_cl_mode_B import *
from accelerator.lattice.xfel_cl_mode_B import *
from accelerator.lattice.xfel_sase1_mode_B import *
import numpy as np
from accelerator.s2e_sections.section_track import *


class ACC1(SectionTrack):

    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'Injector ACC1'
        self.unit_step = 0.01

        self.input_beam_file = 'track_data/particles/particles0_Injector_Gun_end_200kp.npz'
        self.output_beam_file = 'track_data/particles/particles1_Injector_ACC1_end.npz'
        self.tws_file = "track_data/tws/tws1_Injector_ACC1_end.npz"

        # init tracking lattice
        start_sim = id_22433449_
        acc1_stop = id_68749308_

        self.lattice = MagneticLattice(cell_i1, start=start_sim, stop=acc1_stop, method=self.method)

        # init physics processes
        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [63, 63, 63]

        sc2 = SpaceCharge()
        sc2.step = 10
        sc2.nmesh_xyz = [63, 63, 63]

        wake = WakeKick()
        wake.wake_table = WakeTable('accelerator/wakes/mod_TESLA_MODULE_WAKE_TAYLOR.dat')
        wake.factor = 1

        # adding physics processes
        acc1_1_stop = id_75115473_
        acc1_wake_kick = id_68749308_
        self.add_physics_process(sc, start=start_sim, stop=acc1_1_stop)
        self.add_physics_process(sc2, start=acc1_1_stop, stop=acc1_wake_kick)
        self.add_physics_process(wake, start=acc1_wake_kick, stop=acc1_wake_kick)


class ACC39(SectionTrack):

    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'Injector ACC39'
        self.unit_step = 0.01 * 50.0

        self.input_beam_file = 'track_data/particles/particles1_Injector_ACC1_end.npz'
        self.output_beam_file = 'track_data/particles/particles2_Injector_ACC39_end.npz'
        self.tws_file = "track_data/tws/tws2_Injector_ACC39_end.npz"
        # init tracking lattice
        acc1_stop = id_68749308_
        acc39_stop = stlat_47_i1
        self.lattice = MagneticLattice(cell_i1, start=acc1_stop, stop=acc39_stop, method=self.method)

        # init physics processes
        bt = BeamTransform(x_opt=[-2.2407528, 18.45501, 0], y_opt=[-2.2407528, 18.45501, 0])
        bt.bounds = [-5.0, 5.0]

        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [63, 63, 63]

        wake = WakeKick()
        wake.wake_table = WakeTable('accelerator/wakes/mod_THIRD_HARMONIC_SECTION_WAKE_TAYLOR.dat')
        wake.factor = 2

        # adding physics processes
        match_acc39 = acc1_stop
        acc39_wake_kick = stlat_47_i1
        self.add_physics_process(bt, start=match_acc39, stop=match_acc39)
        self.add_physics_process(sc, start=match_acc39, stop=acc39_wake_kick)
        self.add_physics_process(wake, start=acc39_wake_kick, stop=acc39_wake_kick)


class LH(SectionTrack):

    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'LASER HEATER MAGNETS'
        self.unit_step = 0.1

        self.input_beam_file = 'track_data/particles/particles2_Injector_ACC39_end.npz'
        self.output_beam_file = 'track_data/particles/particles3_LaserHeaterMagnets1_end.npz'
        self.tws_file = "track_data/tws/tws3_LaserHeaterMagnets1_end.npz"
        # init tracking lattice
        acc39_stop = stlat_47_i1
        lhm_stop = enlat_50_i1
        self.lattice = MagneticLattice(cell_i1, start=acc39_stop, stop=lhm_stop, method=self.method)

        # init physics processes
        bt = BeamTransform(x_opt=[2.713046129, 8.8913627, 0], y_opt=[-2.736921396, 10.95398564, 0])
        bt.bounds = [-3.0, 3.0]

        csr = CSR()
        csr.traj_step = 0.0005
        csr.apply_step = 0.005

        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [63, 63, 63]

        match_lhm = stlat_47_i1
        self.add_physics_process(bt, start=match_lhm, stop=match_lhm)
        self.add_physics_process(sc, start=match_lhm, stop=match_lhm)
        self.add_physics_process(csr, start=match_lhm, stop=lhm_stop)


class ST2(SectionTrack):
    
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)
        
        # setting parameters
        self.lattice_name = 'STRAIGHT SECTION 2'
        self.unit_step = 1.0

        self.input_beam_file = 'track_data/particles/particles3_LaserHeaterMagnets1_end.npz'
        self.output_beam_file = 'track_data/particles/particles4_StraighSection2_end.npz'
        self.tws_file = "track_data/tws/tws4_StraighSection2_end.npz"

        # init tracking lattice
        lhm_stop = enlat_50_i1
        st2_stop = id_90904668_
        self.lattice = MagneticLattice(cell_i1 + cell_l1, start=lhm_stop, stop=st2_stop, method=self.method)

        # init physics processes
        bt = BeamTransform(x_opt=[0.010627616729, 1.063601468, 0], y_opt=[2.157218065, 15.9609948, 0])
        bt.bounds = [-4.0, 4.0]

        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [63, 63, 63]
        match_st2 = lhm_stop
        self.add_physics_process(bt, start=match_st2, stop=match_st2)
        self.add_physics_process(sc, start=match_st2, stop=st2_stop)


class DL(SectionTrack):
    
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)
        
        # setting parameters
        self.lattice_name = 'DOGLEG'
        self.unit_step = 0.1

        self.input_beam_file = 'track_data/particles/particles4_StraighSection2_end.npz'
        self.output_beam_file = 'track_data/particles/particles5_dogleg_end.npz'
        self.tws_file = "track_data/tws/tws5_dogleg_end.npz"
        # init tracking lattice
        st2_stop = id_90904668_
        dogleg_stop = enlat_93_i1
        self.lattice = MagneticLattice(cell_l1, start=st2_stop, stop=dogleg_stop, method=self.method)

        # init physics processes
        bt = BeamTransform(x_opt=[2.51620, 4.65599, 0], y_opt=[-1.1327, 1.18986, 0])
        bt.bounds = [-4.0, 4.0]

        csr = CSR()
        csr.traj_step = 0.0005
        csr.apply_step = 0.005

        #sc = SpaceCharge()
        #sc.step = 1
        #sc.nmesh_xyz = [63, 63, 63]
        match_dogleg = st2_stop
        self.add_physics_process(bt, start=match_dogleg, stop=match_dogleg)
        self.add_physics_process(csr, start=st2_stop, stop=dogleg_stop)


class ST4(SectionTrack):
    
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'STRAIGHT SECTION 4'
        self.unit_step = 1.0

        self.input_beam_file = 'track_data/particles/particles5_dogleg_end.npz'
        self.output_beam_file = 'track_data/particles/particles6_StraighSection4_end.npz'
        self.tws_file = "track_data/tws/tws6_StraighSection4_end.npz"

        # init tracking lattice
        dogleg_stop = enlat_93_i1
        st4_stop = stlat_96_i1
        self.lattice = MagneticLattice(cell_l1, start=dogleg_stop, stop=st4_stop, method=self.method)
        
        # init physics processes
        bt = BeamTransform(x_opt=[-2.5868502, 3.9779798, 0], y_opt=[0.5395814718, 1.004727395, 0])
        bt.bounds = [-4.0, 4.0]

        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [63, 63, 63]
        match_st4 = enlat_93_i1
        self.add_physics_process(bt, start=match_st4, stop=match_st4)
        self.add_physics_process(sc, start=match_st4, stop=st4_stop)


class BC0(SectionTrack):

    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'BUNCH COMPRESSOR 0'
        self.unit_step = 0.1

        self.input_beam_file = 'track_data/particles/particles6_StraighSection4_end.npz'
        self.output_beam_file = 'track_data/particles/particles7_BunchCompressor0_end.npz'
        self.tws_file = "track_data/tws/tws7_BunchCompressor0_end.npz"
        # init tracking lattice
        st4_stop = stlat_96_i1
        bc0_stop = enlat_101_i1
        self.lattice = MagneticLattice(cell_l1, start=st4_stop, stop=bc0_stop, method=self.method)

        # init physics processes
        bt = BeamTransform(x_opt=[-0.270371393, 8.104624, 0], y_opt=[1.92407299, 15.803186796, 0])
        bt.bounds = [-4.0, 4.0]

        csr = CSR()
        csr.traj_step = 0.0005
        csr.apply_step = 0.005

        #sc = SpaceCharge()
        #sc.step = 10
        #sc.nmesh_xyz = [63, 63, 63]
        #sc.low_order_kick = False
        match_bc0 = st4_stop
        self.add_physics_process(bt, start=match_bc0, stop=match_bc0)
        self.add_physics_process(csr, start=match_bc0, stop=bc0_stop)


class ACC2(SectionTrack):
    
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'L1 ACC2'
        self.unit_step = 1.0

        self.input_beam_file = 'track_data/particles/particles7_BunchCompressor0_end.npz'
        self.output_beam_file = 'track_data/particles/particles8_ACC2_end.npz'
        self.tws_file = "track_data/tws/tws8_ACC2_end.npz"

        bc0_stop = enlat_101_i1
        acc2_stop = stlat_182_b1
        # init tracking lattice
        self.lattice = MagneticLattice(cell_l1, start=bc0_stop, stop=acc2_stop, method=self.method)

        # init physics processes
        bt = BeamTransform(x_opt=[0.228156038, 8.58270789, 0], y_opt=[0.23745129, 3.55040786, 0])
        bt.bounds = [-4.0, 4.0]

        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [15, 15, 15]

        wake = WakeKick()
        wake.wake_table = WakeTable('accelerator/wakes/mod_TESLA_MODULE_WAKE_TAYLOR.dat')
        wake.factor = 4
        match_acc2 = bc0_stop
        L1_wake_kick = acc2_stop
        self.add_physics_process(bt, start=match_acc2, stop=match_acc2)
        self.add_physics_process(sc, start=match_acc2, stop=L1_wake_kick)
        self.add_physics_process(wake, start=L1_wake_kick, stop=L1_wake_kick)


class BC1(SectionTrack):

    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):

        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'BUNCH COMPRESSOR 1'
        self.unit_step = 0.1

        self.input_beam_file = 'track_data/particles/particles8_ACC2_end.npz'
        self.output_beam_file = 'track_data/particles/particles9_BunchCompressor1_end.npz'
        self.tws_file = "track_data/tws/tws9_BunchCompressor1_end.npz"

        acc2_stop = stlat_182_b1
        bc1_stop = enlat_202_b1
        # init tracking lattice
        self.lattice = MagneticLattice(cell_l1, start=acc2_stop, stop=bc1_stop, method=self.method)

        # init physics processes
        bt = BeamTransform(x_opt=[1.47614849, 42.4595599, 0], y_opt=[2.48688943, 50.76637119, 0])
        bt.bounds = [-4.0, 4.0]

        csr = CSR()
        csr.traj_step = 0.0005
        csr.apply_step = 0.005
        csr.n_bin = 300
        csr.sigma_min = 0.2e-6

        #sc = SpaceCharge()
        #sc.step = 10
        #sc.nmesh_xyz = [63, 63, 63]
        #sc.low_order_kick = False
        match_bc1 = acc2_stop
        self.add_physics_process(bt, start=match_bc1, stop=match_bc1)
        self.add_physics_process(csr, start=match_bc1, stop=bc1_stop)


class ACC3t5(SectionTrack):

    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'L2 ACC3t5'
        self.unit_step = 1.0

        self.input_beam_file = 'track_data/particles/particles9_BunchCompressor1_end.npz'
        self.output_beam_file = 'track_data/particles/particles10_ACC3t5_end.npz'
        self.tws_file = "track_data/tws/tws10_ACC3t5_end.npz"


        bc1_stop = enlat_202_b1
        acc3t5_stop = stlat_393_b2
        # init tracking lattice
        self.lattice = MagneticLattice(cell_l1 + cell_l2, start=bc1_stop, stop=acc3t5_stop, method=self.method)

        # init physics processes
        bt = BeamTransform(x_opt=[-0.14292329, 9.248336532, 0], y_opt=[-0.43871120, 8.425951098, 0])
        bt.bounds = [-3.0, 3.0]

        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [15, 15, 15]

        wake = WakeKick()
        wake.wake_table = WakeTable('accelerator/wakes/mod_TESLA_MODULE_WAKE_TAYLOR.dat')
        wake.factor = 4 * 3
        L2_wake_kick = acc3t5_stop
        match_acc3t5 = bc1_stop
        self.add_physics_process(bt, start=match_acc3t5, stop=match_acc3t5)
        self.add_physics_process(sc, start=match_acc3t5, stop=L2_wake_kick)
        self.add_physics_process(wake, start=L2_wake_kick, stop=L2_wake_kick)


class BC2(SectionTrack):
    
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'BUNCH COMPRESSOR 2'
        self.unit_step = 0.1

        self.input_beam_file = 'track_data/particles/particles10_ACC3t5_end.npz'
        self.output_beam_file = 'track_data/particles/particles11_BunchCompressor2_end.npz'
        self.tws_file = "track_data/tws/tws11_BunchCompressor2_end.npz"


        acc3t5_stop = stlat_393_b2
        bc2_stop = enlat_414_b2
        # init tracking lattice
        self.lattice = MagneticLattice(cell_l2, start=acc3t5_stop, stop=bc2_stop, method=self.method)

        # init physics processes
        bt = BeamTransform(x_opt=[1.5668977, 37.9362938, 0], y_opt=[3.499375186, 79.65474747, 0])
        bt.bounds = [-4.0, 4.0]

        csr = CSR()
        csr.traj_step = 0.0005
        csr.apply_step = 0.005
        csr.n_bin = 300
        csr.sigma_min = 0.2e-6

        sc = SpaceCharge()
        sc.step = 10
        sc.nmesh_xyz = [63, 63, 63]

        match_bc2 = acc3t5_stop
        self.add_physics_process(bt, start=match_bc2, stop=match_bc2)
        self.add_physics_process(csr, start=match_bc2, stop=bc2_stop)
        self.add_physics_process(sc, start=match_bc2, stop=bc2_stop)


class ACC6t26(SectionTrack):
    
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'L3 ACC6t26'
        self.unit_step = 5.0

        self.input_beam_file = 'track_data/particles/particles11_BunchCompressor2_end.npz'
        self.output_beam_file = 'track_data/particles/particles12_ACC6t26_end.npz'
        self.tws_file = "track_data/tws/tws12_ACC6t26_end.npz"


        bc2_stop = enlat_414_b2
        acc6t26_stop = match_1673_cl
        #L3_wake_kick, acc6t26_stop, match_collimator1
        # init tracking lattice
        self.lattice = MagneticLattice(cell_l2 + cell_l3_no_cl + cell_cl, start=bc2_stop, stop=acc6t26_stop, method=self.method)

        # init physics processes
        bt = BeamTransform(x_opt=[-0.34317185, 10.8758710, 0], y_opt=[0.067963199, 6.041441, 0])
        bt.bounds = [-3.0, 3.0]

        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [31, 31, 31]

        wake = WakeKick()
        wake.wake_table = WakeTable('accelerator/wakes/mod_TESLA_MODULE_WAKE_TAYLOR.dat')
        wake.factor = 4 * 21
        match_acc6t26 = bc2_stop
        self.add_physics_process(bt, start=match_acc6t26, stop=match_acc6t26)
        self.add_physics_process(sc, start=match_acc6t26, stop=acc6t26_stop)


class CL1(SectionTrack):
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'CL1'
        self.unit_step = 0.2

        self.input_beam_file = 'track_data/particles/particles12_ACC6t26_end.npz'
        self.output_beam_file = 'track_data/particles/particles13_Collimator1_end.npz'
        self.tws_file = "track_data/tws/tws13_Collimator1_end.npz"


        acc6t26_stop = match_1673_cl
        collimator1_stop = bpma_1746_cl
        # init tracking lattice
        self.lattice = MagneticLattice(cell_cl, start=acc6t26_stop, stop=collimator1_stop, method=self.method)

        # init physics processes
        bt = BeamTransform(x_opt=[-3.46547756, 1.924270, 0], y_opt=[-0.0353452348356, 3.8597009238, 0])
        bt.bounds = [-4.0, 4.0]

        sc = SpaceCharge()
        sc.step = 10
        sc.nmesh_xyz = [63, 63, 63]
        sc.low_order_kick = False

        csr = CSR()
        csr.traj_step = 0.0005
        csr.apply_step = 0.005
        csr.n_bin = 300
        csr.sigma_min = 0.2e-6
        match_collimator1 = acc6t26_stop
        self.add_physics_process(bt, start=acc6t26_stop, stop=match_collimator1)
        self.add_physics_process(csr, start=match_collimator1, stop=collimator1_stop)
        self.add_physics_process(sc, start=match_collimator1, stop=collimator1_stop)


class CL2(SectionTrack):
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'CL2'
        self.unit_step = 1

        self.input_beam_file = 'track_data/particles/particles13_Collimator1_end.npz'
        self.output_beam_file = 'track_data/particles/particles14_Collimator2_end.npz'
        self.tws_file = "track_data/tws/tws14_Collimator2_end.npz"


        collimator1_stop = bpma_1746_cl
        collimator2_stop = bpma_1783_cl
        # init tracking lattice
        self.lattice = MagneticLattice(cell_cl, start=collimator1_stop, stop=collimator2_stop, method=self.method)

        # init physics processes
        bt = BeamTransform(x_opt=[-1.6362969895, 0.543989864, 0], y_opt=[0.03577769164, 3.856943, 0])
        bt.bounds = [-3.0, 3.0]

        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [15, 15, 15]
        sc.low_order_kick = False
        match_collimator2 = collimator1_stop
        self.add_physics_process(bt, start=collimator1_stop, stop=match_collimator2)
        self.add_physics_process(sc, start=match_collimator2, stop=collimator2_stop)


class CL3(SectionTrack):
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'CL3'
        self.unit_step = 0.2

        self.input_beam_file = 'track_data/particles/particles14_Collimator2_end.npz'
        self.output_beam_file = 'track_data/particles/particles15_Collimator3_end.npz'
        self.tws_file = "track_data/tws/tws15_Collimator3_end.npz"


        collimator2_stop = bpma_1783_cl
        collimator3_stop = ensec_1854_cl
        # init tracking lattice
        self.lattice = MagneticLattice(cell_cl, start=collimator2_stop, stop=collimator3_stop, method=self.method)

        # init physics processes
        bt = BeamTransform(x_opt=[-2.079721869, 1.37738079, 0], y_opt=[0.1627387392799, 2.21406649259, 0])
        bt.bounds = [-3.0, 3.0]

        sc = SpaceCharge()
        sc.step = 10
        sc.nmesh_xyz = [63, 63, 63]
        sc.low_order_kick = False

        csr = CSR()
        csr.traj_step = 0.0005
        csr.apply_step = 0.005
        csr.n_bin = 300
        csr.sigma_min = 0.2e-6

        match_collimator3 = collimator2_stop

        self.add_physics_process(bt, start=collimator2_stop, stop=match_collimator3)
        self.add_physics_process(csr, start=match_collimator3, stop=collimator3_stop)
        self.add_physics_process(sc, start=match_collimator3, stop=collimator3_stop)


class STN10(SectionTrack):
    def __init__(self, sc_flag=True, csr_flag=True, wake_flag=True):
        SectionTrack.__init__(self, sc_flag, csr_flag, wake_flag)

        # setting parameters
        self.lattice_name = 'ST10'
        self.unit_step = 10

        self.input_beam_file = 'track_data/particles/particles15_Collimator3_end.npz'
        self.output_beam_file = 'track_data/particles/particles16_STN10_end.npz'
        self.tws_file = "track_data/tws/tws16_STN10_end.npz"

        collimator3_stop = ensec_1854_cl
        stN10_stop = ensec_2235_t2
        # init tracking lattice
        self.lattice = MagneticLattice(cell_cl + cell_sase1, start=collimator3_stop, stop=stN10_stop, method=self.method)

        # init physics processes
        bt = BeamTransform(x_opt=[-1.03338305, 0.534912526, 0], y_opt=[0.290016573, 2.336701478, 0])
        bt.bounds = [-3.0, 3.0]

        sc = SpaceCharge()
        sc.step = 1
        sc.nmesh_xyz = [15, 15, 15]
        sc.low_order_kick = False

        match_stN10 = collimator3_stop
        self.add_physics_process(bt, start=collimator3_stop, stop=match_stN10)
        self.add_physics_process(sc, start=match_stN10, stop=stN10_stop)

sections = [ACC1, ACC39, LH, ST2, DL, ST4, BC0, ACC2, BC1, ACC3t5, BC2, ACC6t26, CL1,
            CL2, CL3, STN10]