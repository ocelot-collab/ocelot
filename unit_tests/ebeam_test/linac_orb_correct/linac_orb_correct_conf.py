"""Test parameters description"""

import pytest

from ocelot import *

"""lattice elements description"""



# drifts
d_13 = Drift(l=0.182428, eid='D_13')
d_14 = Drift(l=5.1e-05, eid='D_14')
cbl_73_i1 = Drift(l=0.0, eid='CBL.73.I1')
d_15 = Drift(l=0.135, eid='D_15')
d_16 = Drift(l=0.239432, eid='D_16')
d_17 = Drift(l=0.0751, eid='D_17')
d_19 = Drift(l=0.13494, eid='D_19')
d_20 = Drift(l=0.158616, eid='D_20')
d_21 = Drift(l=0.1325, eid='D_21')
d_22 = Drift(l=0.09865, eid='D_22')
d_23 = Drift(l=0.11615, eid='D_23')
d_24 = Drift(l=0.273615, eid='D_24')
d_28 = Drift(l=0.324432, eid='D_28')
d_29 = Drift(l=0.150051, eid='D_29')
d_31 = Drift(l=0.182429, eid='D_31')
d_37 = Drift(l=0.134939, eid='D_37')
d_42 = Drift(l=0.273616, eid='D_42')
d_46 = Drift(l=0.324431, eid='D_46')
d_47 = Drift(l=0.150052, eid='D_47')
d_50 = Drift(l=5.2e-05, eid='D_50')
d_52 = Drift(l=0.239431, eid='D_52')
d_74 = Drift(l=8e-06, eid='D_74')
d_75 = Drift(l=0.158608, eid='D_75')
d_78 = Drift(l=0.489766, eid='D_78')
d_82 = Drift(l=0.091932, eid='D_82')

# quadrupoles
qi_74_i1 = Quadrupole(l=0.2377, k1=-4.99137522, tilt=0.0, eid='QI.74.I1')
qi_74_i1.dx = 0.001
qi_74_i1.dy = 0.001
qi_75_i1 = Quadrupole(l=0.2377, k1=5.027338701, tilt=0.0, eid='QI.75.I1')
qi_77_i1 = Quadrupole(l=0.2377, k1=-4.99137522, tilt=0.0, eid='QI.77.I1')
qi_78_i1 = Quadrupole(l=0.2377, k1=4.632555992, tilt=0.0, eid='QI.78.I1')
qi_79_i1 = Quadrupole(l=0.2377, k1=-4.99137522, tilt=0.0, eid='QI.79.I1')
qi_80_i1 = Quadrupole(l=0.2377, k1=5.027338701, tilt=0.0, eid='QI.80.I1')
qi_82_i1 = Quadrupole(l=0.2377, k1=-4.99137522, tilt=0.0, eid='QI.82.I1')
qi_83_i1 = Quadrupole(l=0.2377, k1=4.632555992, tilt=0.0, eid='QI.83.I1')
qi_84_i1 = Quadrupole(l=0.2377, k1=-4.99137522, tilt=0.0, eid='QI.84.I1')
qi_85_i1 = Quadrupole(l=0.2377, k1=5.027338701, tilt=0.0, eid='QI.85.I1')
qi_86_i1 = Quadrupole(l=0.2377, k1=-4.99137522, tilt=0.0, eid='QI.86.I1')
qi_88_i1 = Quadrupole(l=0.2377, k1=4.632555992, tilt=0.0, eid='QI.88.I1')
qi_89_i1 = Quadrupole(l=0.2377, k1=-4.99137522, tilt=0.0, eid='QI.89.I1')
qi_90_i1 = Quadrupole(l=0.2377, k1=5.027338701, tilt=0.0, eid='QI.90.I1')
qi_92_i1 = Quadrupole(l=0.2377, k1=-4.99137522, tilt=0.0, eid='QI.92.I1')

# bending magnets
bl_73_i1 = SBend(l=0.200102, angle=-0.1109740393, e1=-0.05548702, e2=-0.05548702, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.73.I1')
bl_75_i1 = SBend(l=0.200102, angle=0.0426524581, e1=0.021326229, e2=0.021326229, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.75.I1')
bl_76_i1 = SBend(l=0.200102, angle=0.0426524581, e1=0.021326229, e2=0.021326229, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.76.I1')
bl_77_i1 = SBend(l=0.200102, angle=-0.1109740393, e1=-0.05548702, e2=-0.05548702, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.77.I1')
bl_78_i1 = SBend(l=0.200102, angle=-0.1109740393, e1=-0.05548702, e2=-0.05548702, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.78.I1')
bl_80_i1 = SBend(l=0.200102, angle=0.0426524581, e1=0.021326229, e2=0.021326229, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.80.I1')
bl_81_i1 = SBend(l=0.200102, angle=0.0426524581, e1=0.021326229, e2=0.021326229, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.81.I1')
bl_82_i1 = SBend(l=0.200102, angle=-0.1109740393, e1=-0.05548702, e2=-0.05548702, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.82.I1')
bl_83_i1 = SBend(l=0.200102, angle=0.1109740393, e1=0.05548702, e2=0.05548702, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.83.I1')
bl_85_i1 = SBend(l=0.200102, angle=-0.0426524581, e1=-0.021326229, e2=-0.021326229, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.85.I1')
bl_86_i1 = SBend(l=0.200102, angle=-0.0426524581, e1=-0.021326229, e2=-0.021326229, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.86.I1')
bl_87_i1 = SBend(l=0.200102, angle=0.1109740393, e1=0.05548702, e2=0.05548702, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.87.I1')
bl_88_i1 = SBend(l=0.200102, angle=0.1109740393, e1=0.05548702, e2=0.05548702, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.88.I1')
bl_90_i1 = SBend(l=0.200102, angle=-0.0426524581, e1=-0.021326229, e2=-0.021326229, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.90.I1')
bl_91_i1 = SBend(l=0.200102, angle=-0.0426524581, e1=-0.021326229, e2=-0.021326229, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.91.I1')
bl_92_i1 = SBend(l=0.200102, angle=0.1109740393, e1=0.05548702, e2=0.05548702, gap=0, tilt=1.570796327, fint=0.0, fintx=0.0, eid='BL.92.I1')

# correctors
cix_73ii_i1 = Hcor(l=0.1, angle=0.0, eid='CIX.73II.I1')
ciy_75_i1 = Vcor(l=0.1, angle=0.0, eid='CIY.75.I1')
cix_76_i1 = Hcor(l=0.1, angle=0.0, eid='CIX.76.I1')
cix_78_i1 = Hcor(l=0.1, angle=0.0, eid='CIX.78.I1')
ciy_80_i1 = Vcor(l=0.1, angle=0.0, eid='CIY.80.I1')
cix_81_i1 = Hcor(l=0.1, angle=0.0, eid='CIX.81.I1')
cix_83_i1 = Hcor(l=0.1, angle=0.0, eid='CIX.83.I1')
ciy_85_i1 = Vcor(l=0.1, angle=0.0, eid='CIY.85.I1')
cix_86_i1 = Hcor(l=0.1, angle=0.0, eid='CIX.86.I1')
cix_88_i1 = Hcor(l=0.1, angle=0.0, eid='CIX.88.I1')
cix_90_i1 = Hcor(l=0.1, angle=0.0, eid='CIX.90.I1')
ciy_92_i1 = Vcor(l=0.1, angle=0.0, eid='CIY.92.I1')

# markers
id_90904668_ = Marker(eid='ID_90904668_')
id_60455493_ = Marker(eid='ID_60455493_')
id_74558584_ = Marker(eid='ID_74558584_')

# monitor
bpma_75_i1 = Monitor(eid='BPMA.75.I1')
bpma_77_i1 = Monitor(eid='BPMA.77.I1')
bpma_80_i1 = Monitor(eid='BPMA.80.I1')
bpma_82_i1 = Monitor(eid='BPMA.82.I1')
bpma_85_i1 = Monitor(eid='BPMA.85.I1')
bpma_87_i1 = Monitor(eid='BPMA.87.I1')
bpma_90_i1 = Monitor(eid='BPMA.90.I1')
bpma_92_i1 = Monitor(eid='BPMA.92.I1')

# sextupoles
sc_74i_i1 = Sextupole(l=0.1121, k2=-87.57825836, tilt=1.570796327, eid='SC.74I.I1')
sc_74ii_i1 = Sextupole(l=0.1121, k2=-53.061653289999995, tilt=1.570796327, eid='SC.74II.I1')
sc_76_i1 = Sextupole(l=0.1121, k2=-53.061653289999995, tilt=1.570796327, eid='SC.76.I1')
sc_77_i1 = Sextupole(l=0.1121, k2=-87.57825836, tilt=1.570796327, eid='SC.77.I1')
sc_79i_i1 = Sextupole(l=0.1121, k2=-87.57825836, tilt=1.570796327, eid='SC.79I.I1')
sc_79ii_i1 = Sextupole(l=0.1121, k2=-53.061653289999995, tilt=1.570796327, eid='SC.79II.I1')
sc_81_i1 = Sextupole(l=0.1121, k2=-53.061653289999995, tilt=1.570796327, eid='SC.81.I1')
sc_82_i1 = Sextupole(l=0.1121, k2=-87.57825836, tilt=1.570796327, eid='SC.82.I1')
sc_84i_i1 = Sextupole(l=0.1121, k2=87.57825836, tilt=1.570796327, eid='SC.84I.I1')
sc_84ii_i1 = Sextupole(l=0.1121, k2=53.061653289999995, tilt=1.570796327, eid='SC.84II.I1')
sc_86_i1 = Sextupole(l=0.1121, k2=53.061653289999995, tilt=1.570796327, eid='SC.86.I1')
sc_87_i1 = Sextupole(l=0.1121, k2=87.57825836, tilt=1.570796327, eid='SC.87.I1')
sc_89i_i1 = Sextupole(l=0.1121, k2=87.57825836, tilt=1.570796327, eid='SC.89I.I1')
sc_89ii_i1 = Sextupole(l=0.1121, k2=53.061653289999995, tilt=1.570796327, eid='SC.89II.I1')
sc_91_i1 = Sextupole(l=0.1121, k2=53.061653289999995, tilt=1.570796327, eid='SC.91.I1')
sc_92_i1 = Sextupole(l=0.1121, k2=87.57825836, tilt=1.570796327, eid='SC.92.I1')



"""pytest fixtures description"""

@pytest.fixture(scope='module')
def cell():
    dogleg = (id_90904668_, d_13, id_60455493_, bl_73_i1, d_14, cbl_73_i1, d_15,
            cix_73ii_i1, d_16, sc_74i_i1, d_17, qi_74_i1, d_17, sc_74ii_i1, d_19,
            bl_75_i1, d_20, bpma_75_i1, d_21, ciy_75_i1, d_22, qi_75_i1, d_23,
            cix_76_i1, d_24, bl_76_i1, d_19, sc_76_i1, d_17, qi_77_i1, d_17,
            sc_77_i1, d_28, bpma_77_i1, d_29, bl_77_i1, d_13, qi_78_i1, d_31,
            bl_78_i1, d_14, cbl_73_i1, d_15, cix_78_i1, d_16, sc_79i_i1, d_17,
            qi_79_i1, d_17, sc_79ii_i1, d_37, bl_80_i1, d_20, bpma_80_i1, d_21,
            ciy_80_i1, d_22, qi_80_i1, d_23, cix_81_i1, d_42, bl_81_i1, d_19,
            sc_81_i1, d_17, qi_82_i1, d_17, sc_82_i1, d_46, bpma_82_i1, d_47,
            bl_82_i1, d_13, qi_83_i1, d_13, bl_83_i1, d_50, cbl_73_i1, d_15,
            cix_83_i1, d_52, sc_84i_i1, d_17, qi_84_i1, d_17, sc_84ii_i1, d_19,
            bl_85_i1, d_20, bpma_85_i1, d_21, ciy_85_i1, d_22, qi_85_i1, d_23,
            cix_86_i1, d_42, bl_86_i1, d_37, sc_86_i1, d_17, qi_86_i1, d_17,
            sc_87_i1, d_28, bpma_87_i1, d_47, bl_87_i1, d_13, qi_88_i1, d_13,
            bl_88_i1, d_14, cbl_73_i1, d_15, cix_88_i1, d_16, sc_89i_i1, d_17,
            qi_89_i1, d_17, sc_89ii_i1, d_37, bl_90_i1, d_74, cbl_73_i1, d_75,
            bpma_90_i1, d_21, cix_90_i1, d_22, qi_90_i1, d_78, bl_91_i1, d_19,
            sc_91_i1, d_17, qi_92_i1, d_17, sc_92_i1, d_82, ciy_92_i1, d_21,
            bpma_92_i1, d_29, bl_92_i1, d_13, id_74558584_)
    return dogleg


@pytest.fixture(scope='module')
def method():
    return MethodTM()


@pytest.fixture(scope='module')
def lattice(cell, method):
    return MagneticLattice(cell, method=method)
