"""Test parameters description file"""

import pytest

from ocelot import *

"""Lattice defenition"""
# drifts 
D_1 = Drift(l=0.276, eid='D_1')
D_2 = Drift(l=0.474, eid='D_2')
D_3 = Drift(l=0.25, eid='D_3')
D_4 = Drift(l=0.15, eid='D_4')
D_5 = Drift(l=0.198782, eid='D_5')
D_6 = Drift(l=0.329218, eid='D_6')
D_7 = Drift(l=0.06, eid='D_7')
D_8 = Drift(l=0.225, eid='D_8')
D_9 = Drift(l=0.088, eid='D_9')
D_10 = Drift(l=0.111, eid='D_10')
D_11 = Drift(l=0.05, eid='D_11')
D_12 = Drift(l=0.3864, eid='D_12')

start_sim = Marker()
D_13 = (Drift(l=0.5512-0.0996, eid='D_13a'), start_sim, Drift(l=0.0996, eid='D_13b') )

D_14 = Drift(l=0.2216, eid='D_14')
D_15 = Drift(l=0.3459, eid='D_15')
D_16 = Drift(l=0.3459, eid='D_16')
D_17 = Drift(l=0.3459, eid='D_17')
D_18 = Drift(l=0.3459, eid='D_18')
D_19 = Drift(l=0.3459, eid='D_19')
D_20 = Drift(l=0.3459, eid='D_20')
D_21 = Drift(l=0.3459, eid='D_21')
D_22 = Drift(l=0.2043, eid='D_22')
D_23 = Drift(l=0.085, eid='D_23')
D_24 = Drift(l=0.4579, eid='D_24')
D_25 = Drift(l=0.2211, eid='D_25')
D_26 = Drift(l=0.085, eid='D_26')
D_27 = Drift(l=0.202, eid='D_27')
D_28 = Drift(l=0.262, eid='D_28')
D_29 = Drift(l=0.262, eid='D_29')
D_30 = Drift(l=0.262, eid='D_30')
D_31 = Drift(l=0.262, eid='D_31')
D_32 = Drift(l=0.262, eid='D_32')
D_33 = Drift(l=0.262, eid='D_33')
D_34 = Drift(l=0.262, eid='D_34')
D_35 = Drift(l=0.155, eid='D_35')
D_36 = Drift(l=1.04, eid='D_36')
D_37 = Drift(l=1.2964, eid='D_37')
D_38 = Drift(l=0.3519, eid='D_38')
D_39 = Drift(l=0.105, eid='D_39')
D_40 = Drift(l=0.2175, eid='D_40')
D_41 = Drift(l=0.096, eid='D_41')
D_42 = Drift(l=0.1177, eid='D_42')
D_43 = Drift(l=0.1423, eid='D_43')
D_44 = Drift(l=0.1, eid='D_44')
D_45 = Drift(l=0.0509998582, eid='D_45')
D_46 = Drift(l=0.1004967165, eid='D_46')
D_47 = Drift(l=0.0922498582, eid='D_47')
D_48 = Drift(l=0.096, eid='D_48')
D_49 = Drift(l=0.27175, eid='D_49')
D_50 = Drift(l=0.194, eid='D_50')
D_51 = Drift(l=0.194, eid='D_51')
D_52 = Drift(l=0.1199998582, eid='D_52')
D_53 = Drift(l=0.1004977165, eid='D_53')
D_54 = Drift(l=0.0509998582, eid='D_54')
D_55 = Drift(l=0.1, eid='D_55')
D_56 = Drift(l=0.20849, eid='D_56')
D_57 = Drift(l=0.2, eid='D_57')
D_58 = Drift(l=0.303, eid='D_58')
D_59 = Drift(l=0.1975, eid='D_59')
D_60 = Drift(l=0.096, eid='D_60')
D_61 = Drift(l=0.1065, eid='D_61')
D_62 = Drift(l=0.07, eid='D_62')
D_63 = Drift(l=0.07, eid='D_63')
D_64 = Drift(l=0.752, eid='D_64')
D_65 = Drift(l=0.74, eid='D_65')
D_66 = Drift(l=0.035, eid='D_66')
D_67 = Drift(l=0.175, eid='D_67')
D_68 = Drift(l=0.15, eid='D_68')
D_69 = Drift(l=0.15, eid='D_69')
D_70 = Drift(l=0.775, eid='D_70')
D_71 = Drift(l=0.175, eid='D_71')
D_72 = Drift(l=0.15, eid='D_72')
D_73 = Drift(l=0.15, eid='D_73')
D_74 = Drift(l=0.775, eid='D_74')
D_75 = Drift(l=0.275, eid='D_75')
D_76 = Drift(l=0.2, eid='D_76')
D_77 = Drift(l=0.22, eid='D_77')
D_78 = Drift(l=0.38, eid='D_78')
D_79 = Drift(l=0.3, eid='D_79')
D_80 = Drift(l=0.2817, eid='D_80')
D_81 = Drift(l=0.3643, eid='D_81')
D_82 = Drift(l=0.15, eid='D_82')
D_83 = Drift(l=0.254, eid='D_83')
D_84 = Drift(l=0.55, eid='D_84')
D_85 = Drift(l=1.24, eid='D_85')
D_86 = Drift(l=0.1, eid='D_86')
D_87 = Drift(l=0.15, eid='D_87')
D_88 = Drift(l=1.45, eid='D_88')
D_89 = Drift(l=0.65, eid='D_89')
D_90 = Drift(l=2.9, eid='D_90')
D_91 = Drift(l=1.8, eid='D_91')
D_92 = Drift(l=0.7, eid='D_92')
D_93 = Drift(l=0.2, eid='D_93')
D_94 = Drift(l=0.2, eid='D_94')
D_95 = Drift(l=0.2, eid='D_95')
D_96 = Drift(l=0.1, eid='D_96')
D_97 = Drift(l=0.2012266681, eid='D_97')
D_98 = Drift(l=0.1349996681, eid='D_98')
D_99 = Drift(l=0.245482, eid='D_99')
D_100 = Drift(l=0.1, eid='D_100')
D_101 = Drift(l=0.1, eid='D_101')
D_102 = Drift(l=0.1409824195, eid='D_102')
D_103 = Drift(l=0.1586084195, eid='D_103')
D_104 = Drift(l=0.1325, eid='D_104')
D_105 = Drift(l=0.1175, eid='D_105')
D_106 = Drift(l=0.135, eid='D_106')
D_107 = Drift(l=0.2736074195, eid='D_107')
D_108 = Drift(l=0.1409824195, eid='D_108')
D_109 = Drift(l=0.1, eid='D_109')
D_110 = Drift(l=0.1, eid='D_110')
D_111 = Drift(l=0.330482, eid='D_111')
D_112 = Drift(l=0.1499996681, eid='D_112')
D_113 = Drift(l=0.2012266681, eid='D_113')
D_114 = Drift(l=0.2012276681, eid='D_114')
D_115 = Drift(l=0.1349996681, eid='D_115')
D_116 = Drift(l=0.245482, eid='D_116')
D_117 = Drift(l=0.1, eid='D_117')
D_118 = Drift(l=0.1, eid='D_118')
D_119 = Drift(l=0.1409814195, eid='D_119')
D_120 = Drift(l=0.1586084195, eid='D_120')
D_121 = Drift(l=0.1325, eid='D_121')
D_122 = Drift(l=0.1175, eid='D_122')
D_123 = Drift(l=0.135, eid='D_123')
D_124 = Drift(l=0.2736084195, eid='D_124')
D_125 = Drift(l=0.1409824195, eid='D_125')
D_126 = Drift(l=0.1, eid='D_126')
D_127 = Drift(l=0.1, eid='D_127')
D_128 = Drift(l=0.330481, eid='D_128')
D_129 = Drift(l=0.1500006681, eid='D_129')
D_130 = Drift(l=0.2012266681, eid='D_130')
D_131 = Drift(l=0.2012266681, eid='D_131')
D_132 = Drift(l=6.681e-07, eid='D_132')
D_133 = Drift(l=0.135, eid='D_133')
D_134 = Drift(l=0.245481, eid='D_134')
D_135 = Drift(l=0.1, eid='D_135')
D_136 = Drift(l=0.1, eid='D_136')
D_137 = Drift(l=0.1409824195, eid='D_137')
D_138 = Drift(l=0.1586084195, eid='D_138')
D_139 = Drift(l=0.1325, eid='D_139')
D_140 = Drift(l=0.1175, eid='D_140')
D_141 = Drift(l=0.135, eid='D_141')
D_142 = Drift(l=0.2736084195, eid='D_142')
D_143 = Drift(l=0.1409814195, eid='D_143')
D_144 = Drift(l=0.1, eid='D_144')
D_145 = Drift(l=0.1, eid='D_145')
D_146 = Drift(l=0.330482, eid='D_146')
D_147 = Drift(l=0.1500006681, eid='D_147')
D_148 = Drift(l=0.2012266681, eid='D_148')
D_149 = Drift(l=0.2012266681, eid='D_149')
D_150 = Drift(l=0.1349996681, eid='D_150')
D_151 = Drift(l=0.245482, eid='D_151')
D_152 = Drift(l=0.1, eid='D_152')
D_153 = Drift(l=0.1, eid='D_153')
D_154 = Drift(l=0.1409814195, eid='D_154')
D_155 = Drift(l=4.195e-07, eid='D_155')
D_156 = Drift(l=0.158608, eid='D_156')
D_157 = Drift(l=0.1325, eid='D_157')
D_158 = Drift(l=0.1175, eid='D_158')
D_159 = Drift(l=0.5086084195, eid='D_159')
D_160 = Drift(l=0.1409824195, eid='D_160')
D_161 = Drift(l=0.1, eid='D_161')
D_162 = Drift(l=0.1, eid='D_162')
D_163 = Drift(l=0.097982, eid='D_163')
D_164 = Drift(l=0.1325, eid='D_164')
D_165 = Drift(l=0.1499996681, eid='D_165')
D_166 = Drift(l=0.2012266681, eid='D_166')
D_167 = Drift(l=0.90013, eid='D_167')
D_168 = Drift(l=0.2, eid='D_168')
D_169 = Drift(l=0.11, eid='D_169')
D_170 = Drift(l=0.3, eid='D_170')
D_171 = Drift(l=0.096, eid='D_171')
D_172 = Drift(l=0.204, eid='D_172')
D_173 = Drift(l=0.11, eid='D_173')
D_174 = Drift(l=0.4, eid='D_174')
D_175 = Drift(l=0.1000004395, eid='D_175')
D_176 = Drift(l=1.0088728791, eid='D_176')
D_177 = Drift(l=4.395e-07, eid='D_177')
D_178 = Drift(l=0.865, eid='D_178')
D_179 = Drift(l=0.31, eid='D_179')
D_180 = Drift(l=0.3250004395, eid='D_180')
D_181 = Drift(l=1.0088728791, eid='D_181')
D_182 = Drift(l=4.395e-07, eid='D_182')
D_183 = Drift(l=0.1, eid='D_183')
D_184 = Drift(l=0.41, eid='D_184')
D_185 = Drift(l=0.11, eid='D_185')
D_186 = Drift(l=0.25, eid='D_186')
D_187 = Drift(l=0.096, eid='D_187')
D_188 = Drift(l=0.144, eid='D_188')
D_189 = Drift(l=0.11, eid='D_189')
D_190 = Drift(l=0.2, eid='D_190')
D_191 = Drift(l=0.56, eid='D_191')
D_192 = Drift(l=0.11, eid='D_192')
D_193 = Drift(l=0.2, eid='D_193')
D_194 = Drift(l=1.89, eid='D_194')
D_195 = Drift(l=0.11, eid='D_195')
D_196 = Drift(l=0.2, eid='D_196')
D_197 = Drift(l=1.89, eid='D_197')
D_198 = Drift(l=0.11, eid='D_198')
D_199 = Drift(l=0.2, eid='D_199')
D_200 = Drift(l=1.89, eid='D_200')
D_201 = Drift(l=0.11, eid='D_201')
D_202 = Drift(l=0.2, eid='D_202')
D_203 = Drift(l=2.04, eid='D_203')
D_204 = Drift(l=0.11, eid='D_204')
D_205 = Drift(l=0.2, eid='D_205')
D_206 = Drift(l=1.1242, eid='D_206')
D_207 = Drift(l=0.1408, eid='D_207')
D_208 = Drift(l=0.11, eid='D_208')
D_209 = Drift(l=0.1585, eid='D_209')
D_210 = Drift(l=1.09132, eid='D_210')
D_211 = Drift(l=0.34903, eid='D_211')
D_212 = Drift(l=0.09115, eid='D_212')
D_213 = Drift(l=0.11, eid='D_213')
D_214 = Drift(l=0.1628, eid='D_214')
D_215 = Drift(l=0.3522, eid='D_215')

# quadrupoles 
Q_37_I1 = Quadrupole(l=0.3, k1=-1.537886, tilt=0.0, eid='Q.37.I1')
Q_38_I1 = Quadrupole(l=0.3, k1=1.435078, tilt=0.0, eid='Q.38.I1')
QI_46_I1 = Quadrupole(l=0.2, k1=-0.0932136, tilt=0.0, eid='QI.46.I1')
QI_47_I1 = Quadrupole(l=0.2, k1=0.8875443, tilt=0.0, eid='QI.47.I1')
QI_50_I1 = Quadrupole(l=0.2, k1=-0.7654263, tilt=0.0, eid='QI.50.I1')
QI_52_I1 = Quadrupole(l=0.2, k1=-0.07117866, tilt=0.0, eid='QI.52.I1')
QI_53_I1 = Quadrupole(l=0.2, k1=2.494601, tilt=0.0, eid='QI.53.I1')
QI_54_I1 = Quadrupole(l=0.2, k1=0.9480791, tilt=0.0, eid='QI.54.I1')
QI_55_I1 = Quadrupole(l=0.2, k1=-4.15374, tilt=0.0, eid='QI.55.I1')
QI_57_I1 = Quadrupole(l=0.2, k1=4.15355, tilt=0.0, eid='QI.57.I1')
QI_59_I1 = Quadrupole(l=0.2, k1=-4.15374, tilt=0.0, eid='QI.59.I1')
QI_60_I1 = Quadrupole(l=0.2, k1=2.128, tilt=0.0, eid='QI.60.I1')
QI_61_I1 = Quadrupole(l=0.2, k1=0.933, tilt=0.0, eid='QI.61.I1')
QI_63_I1 = Quadrupole(l=0.2, k1=-2.2531364, tilt=0.0, eid='QI.63.I1')
QI_66_I1 = Quadrupole(l=0.2, k1=2.441008, tilt=0.0, eid='QI.66.I1')
QI_69_I1 = Quadrupole(l=0.2, k1=-2.559564, tilt=0.0, eid='QI.69.I1')
QI_71_I1 = Quadrupole(l=0.2, k1=3.653292, tilt=0.0, eid='QI.71.I1')
QI_72_I1 = Quadrupole(l=0.2, k1=-4.341087, tilt=0.0, eid='QI.72.I1')
QI_73_I1 = Quadrupole(l=0.2, k1=5.444766974, tilt=0.0, eid='QI.73.I1')
QI_74_I1 = Quadrupole(l=0.2, k1=-5.865493569, tilt=0.0, eid='QI.74.I1')
QI_75_I1 = Quadrupole(l=0.2, k1=5.906955582, tilt=0.0, eid='QI.75.I1')
QI_77_I1 = Quadrupole(l=0.2, k1=-5.865493569, tilt=0.0, eid='QI.77.I1')
QI_78_I1 = Quadrupole(l=0.2, k1=5.444766974, tilt=0.0, eid='QI.78.I1')
QI_79_I1 = Quadrupole(l=0.2, k1=-5.865493569, tilt=0.0, eid='QI.79.I1')
QI_80_I1 = Quadrupole(l=0.2, k1=5.906955582, tilt=0.0, eid='QI.80.I1')
QI_82_I1 = Quadrupole(l=0.2, k1=-5.865493569, tilt=0.0, eid='QI.82.I1')
QI_83_I1 = Quadrupole(l=0.2, k1=5.444766974, tilt=0.0, eid='QI.83.I1')
QI_84_I1 = Quadrupole(l=0.2, k1=-5.865493569, tilt=0.0, eid='QI.84.I1')
QI_85_I1 = Quadrupole(l=0.2, k1=5.906955582, tilt=0.0, eid='QI.85.I1')
QI_86_I1 = Quadrupole(l=0.2, k1=-5.865493569, tilt=0.0, eid='QI.86.I1')
QI_88_I1 = Quadrupole(l=0.2, k1=5.444766974, tilt=0.0, eid='QI.88.I1')
QI_89_I1 = Quadrupole(l=0.2, k1=-5.865493569, tilt=0.0, eid='QI.89.I1')
QI_90_I1 = Quadrupole(l=0.2, k1=5.906955582, tilt=0.0, eid='QI.90.I1')
QI_92_I1 = Quadrupole(l=0.2, k1=-5.865493569, tilt=0.0, eid='QI.92.I1')
QI_93_I1 = Quadrupole(l=0.2, k1=-0.525034487, tilt=0.0, eid='QI.93.I1')
QI_94_I1 = Quadrupole(l=0.2, k1=3.854329117, tilt=0.0, eid='QI.94.I1')
QI_95_I1 = Quadrupole(l=0.2, k1=-3.5693961, tilt=0.0, eid='QI.95.I1')
QI_102_I1 = Quadrupole(l=0.2, k1=0.7805790045, tilt=0.0, eid='QI.102.I1')
QI_103_I1 = Quadrupole(l=0.2, k1=-1.631481454, tilt=0.0, eid='QI.103.I1')
QI_104_I1 = Quadrupole(l=0.2, k1=1.761852655, tilt=0.0, eid='QI.104.I1')
QI_107_I1 = Quadrupole(l=0.2, k1=-1.8, tilt=0.0, eid='QI.107.I1')
QI_109_I1 = Quadrupole(l=0.2, k1=1.8, tilt=0.0, eid='QI.109.I1')
QI_112_I1 = Quadrupole(l=0.2, k1=-1.8, tilt=0.0, eid='QI.112.I1')
QI_114_I1 = Quadrupole(l=0.2, k1=1.184851277, tilt=0.0, eid='QI.114.I1')
QI_116_I1 = Quadrupole(l=0.2, k1=0.6386913113, tilt=0.0, eid='QI.116.I1')
QI_118_I1 = Quadrupole(l=0.2, k1=-1.118664445, tilt=0.0, eid='QI.118.I1')

# bending magnets 
BL_48I_I1 = SBend(l=0.200330283531, angle=-0.099484, e1=0.0, e2=-0.099484, tilt=0.0, fint=0.0, eid='BL.48I.I1')
BL_48II_I1 = SBend(l=0.200330283531, angle=0.099484, e1=0.099484, e2=0.0, tilt=0.0, fint=0.0, eid='BL.48II.I1')
BL_50I_I1 = SBend(l=0.200330283531, angle=0.099484, e1=0.0, e2=0.099484, tilt=0.0, fint=0.0, eid='BL.50I.I1')
BL_50II_I1 = SBend(l=0.200330283531, angle=-0.099484, e1=-0.099484, e2=0.0, tilt=0.0, fint=0.0, eid='BL.50II.I1')
BL_73_I1 = SBend(l=0.200102663853, angle=-0.1109740393, e1=-0.05548702, e2=-0.05548702, tilt=1.570796327, fint=0.0, eid='BL.73.I1')
BL_75_I1 = SBend(l=0.200015161073, angle=0.0426524581, e1=0.021326229, e2=0.021326229, tilt=1.570796327, fint=0.0, eid='BL.75.I1')
BL_76_I1 = SBend(l=0.200015161073, angle=0.0426524581, e1=0.021326229, e2=0.021326229, tilt=1.570796327, fint=0.0, eid='BL.76.I1')
BL_77_I1 = SBend(l=0.200102663853, angle=-0.1109740393, e1=-0.05548702, e2=-0.05548702, tilt=1.570796327, fint=0.0, eid='BL.77.I1')
BL_78_I1 = SBend(l=0.200102663853, angle=-0.1109740393, e1=-0.05548702, e2=-0.05548702, tilt=1.570796327, fint=0.0, eid='BL.78.I1')
BL_80_I1 = SBend(l=0.200015161073, angle=0.0426524581, e1=0.021326229, e2=0.021326229, tilt=1.570796327, fint=0.0, eid='BL.80.I1')
BL_81_I1 = SBend(l=0.200015161073, angle=0.0426524581, e1=0.021326229, e2=0.021326229, tilt=1.570796327, fint=0.0, eid='BL.81.I1')
BL_82_I1 = SBend(l=0.200102663853, angle=-0.1109740393, e1=-0.05548702, e2=-0.05548702, tilt=1.570796327, fint=0.0, eid='BL.82.I1')
BL_83_I1 = SBend(l=0.200102663853, angle=0.1109740393, e1=0.05548702, e2=0.05548702, tilt=1.570796327, fint=0.0, eid='BL.83.I1')
BL_85_I1 = SBend(l=0.200015161073, angle=-0.0426524581, e1=-0.021326229, e2=-0.021326229, tilt=1.570796327, fint=0.0, eid='BL.85.I1')
BL_86_I1 = SBend(l=0.200015161073, angle=-0.0426524581, e1=-0.021326229, e2=-0.021326229, tilt=1.570796327, fint=0.0, eid='BL.86.I1')
BL_87_I1 = SBend(l=0.200102663853, angle=0.1109740393, e1=0.05548702, e2=0.05548702, tilt=1.570796327, fint=0.0, eid='BL.87.I1')
BL_88_I1 = SBend(l=0.200102663853, angle=0.1109740393, e1=0.05548702, e2=0.05548702, tilt=1.570796327, fint=0.0, eid='BL.88.I1')
BL_90_I1 = SBend(l=0.200015161073, angle=-0.0426524581, e1=-0.021326229, e2=-0.021326229, tilt=1.570796327, fint=0.0, eid='BL.90.I1')
BL_91_I1 = SBend(l=0.200015161073, angle=-0.0426524581, e1=-0.021326229, e2=-0.021326229, tilt=1.570796327, fint=0.0, eid='BL.91.I1')
BL_92_I1 = SBend(l=0.200102663853, angle=0.1109740393, e1=0.05548702, e2=0.05548702, tilt=1.570796327, fint=0.0, eid='BL.92.I1')
BB_96_I1 = SBend(l=0.501471120927, angle=0.1327297047, e1=0.0, e2=0.132729705, tilt=1.570796327, fint=0.0, eid='BB.96.I1')
BB_98_I1 = SBend(l=0.501471120927, angle=-0.1327297047, e1=-0.132729705, e2=0.0, tilt=1.570796327, fint=0.0, eid='BB.98.I1')
BB_100_I1 = SBend(l=0.501471120927, angle=-0.1327297047, e1=0.0, e2=-0.132729705, tilt=1.570796327, fint=0.0, eid='BB.100.I1')
BB_101_I1 = SBend(l=0.501471120927, angle=0.1327297047, e1=0.132729705, e2=0.0, tilt=1.570796327, fint=0.0, eid='BB.101.I1')

# correctors 
CKX_23_I1 = Hcor(l=0.025, angle=0.0, eid='CKX.23.I1')
CKY_23_I1 = Vcor(l=0.025, angle=0.0, eid='CKY.23.I1')
CKX_24_I1 = Hcor(l=0.025, angle=0.0, eid='CKX.24.I1')
CKY_24_I1 = Vcor(l=0.025, angle=0.0, eid='CKY.24.I1')
CKX_25_I1 = Hcor(l=0.025, angle=0.0, eid='CKX.25.I1')
CKY_25_I1 = Vcor(l=0.025, angle=0.0, eid='CKY.25.I1')
CX_37_I1 = Hcor(l=0.0, angle=0.0, eid='CX.37.I1')
CY_37_I1 = Vcor(l=0.0, angle=0.0, eid='CY.37.I1')
CX_39_I1 = Hcor(l=0.0, angle=0.0, eid='CX.39.I1')
CY_39_I1 = Vcor(l=0.0, angle=0.0, eid='CY.39.I1')
CIY_51_I1 = Vcor(l=0.1, angle=0.0, eid='CIY.51.I1')
CIX_51_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.51.I1')
CIY_55_I1 = Vcor(l=0.1, angle=0.0, eid='CIY.55.I1')
CIX_57_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.57.I1')
CIY_58_I1 = Vcor(l=0.1, angle=0.0, eid='CIY.58.I1')
CIY_63_I1 = Vcor(l=0.1, angle=0.0, eid='CIY.63.I1')
CIX_65_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.65.I1')
CIY_72_I1 = Vcor(l=0.1, angle=0.0, eid='CIY.72.I1')
CIX_73I_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.73I.I1')
CBL_73_I1 = Vcor(l=0.0, angle=0.0, eid='CBL.73.I1')
CIX_73II_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.73II.I1')
CIY_75_I1 = Vcor(l=0.1, angle=0.0, eid='CIY.75.I1')
CIX_76_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.76.I1')
CBL_78_I1 = Vcor(l=0.0, angle=0.0, eid='CBL.78.I1')
CIX_78_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.78.I1')
CIY_80_I1 = Vcor(l=0.1, angle=0.0, eid='CIY.80.I1')
CIX_81_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.81.I1')
CBL_83_I1 = Vcor(l=0.0, angle=0.0, eid='CBL.83.I1')
CIX_83_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.83.I1')
CIY_85_I1 = Vcor(l=0.1, angle=0.0, eid='CIY.85.I1')
CIX_86_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.86.I1')
CBL_88_I1 = Vcor(l=0.0, angle=0.0, eid='CBL.88.I1')
CIX_88_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.88.I1')
CBL_90_I1 = Vcor(l=0.0, angle=0.0, eid='CBL.90.I1')
CIX_90_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.90.I1')
CIY_92_I1 = Vcor(l=0.1, angle=0.0, eid='CIY.92.I1')
CIY_94_I1 = Vcor(l=0.1, angle=0.0, eid='CIY.94.I1')
CIX_95_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.95.I1')
CBB_98_I1 = Vcor(l=0.0, angle=0.0, eid='CBB.98.I1')
CBB_100_I1 = Vcor(l=0.0, angle=0.0, eid='CBB.100.I1')
CBB_101_I1 = Vcor(l=0.0, angle=0.0, eid='CBB.101.I1')
CIX_102_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.102.I1')
CIY_103_I1 = Vcor(l=0.1, angle=0.0, eid='CIY.103.I1')
CIX_104_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.104.I1')
CIY_107_I1 = Vcor(l=0.1, angle=0.0, eid='CIY.107.I1')
CIX_109_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.109.I1')
CIY_112_I1 = Vcor(l=0.1, angle=0.0, eid='CIY.112.I1')
CIX_114_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.114.I1')
CIY_116_I1 = Vcor(l=0.1, angle=0.0, eid='CIY.116.I1')
CIX_118_I1 = Hcor(l=0.1, angle=0.0, eid='CIX.118.I1')

# markers 
STSEC_23_I1 = Marker(eid='STSEC.23.I1')
STSUB_23_I1 = Marker(eid='STSUB.23.I1')
GUN_23_I1 = Marker(eid='GUN.23.I1')
ENSUB_24_I1 = Marker(eid='ENSUB.24.I1')
STSUB_24_I1 = Marker(eid='STSUB.24.I1')
START_25_I1 = Marker(eid='START.25.I1')
STAC_26_I1 = Marker(eid='STAC.26.I1')
ENAC_38_I1 = Marker(eid='ENAC.38.I1')
STAC_38_I1 = Marker(eid='STAC.38.I1')
ENAC_44_I1 = Marker(eid='ENAC.44.I1')
END_45_I1 = Marker(eid='END.45.I1')
MPBPMF_47_I1 = Marker(eid='MPBPMF.47.I1')
STLAT_47_I1 = Marker(eid='STLAT.47.I1')
MPBPMF_48_I1 = Marker(eid='MPBPMF.48.I1')
ENLAT_50_I1 = Marker(eid='ENLAT.50.I1')
MPBPMF_52_I1 = Marker(eid='MPBPMF.52.I1')
START_55_I1 = Marker(eid='START.55.I1')
ENSUB_62_I1 = Marker(eid='ENSUB.62.I1')
START_73_I1 = Marker(eid='START.73.I1')
END_93_I1 = Marker(eid='END.93.I1')
MPBPMF_95_I1 = Marker(eid='MPBPMF.95.I1')
START_96_I1 = Marker(eid='START.96.I1')
END_101_I1 = Marker(eid='END.101.I1')
MPBPMF_103_I1 = Marker(eid='MPBPMF.103.I1')
START_104_I1 = Marker(eid='START.104.I1')
END_119_I1 = Marker(eid='END.119.I1')

# monitor 
BPMG_24_I1 = Monitor(eid='BPMG.24.I1')
SCRN_24_I1 = Monitor(eid='SCRN.24.I1')
FCUP_24_I1 = Monitor(eid='FCUP.24.I1')
TORA_25_I1 = Monitor(eid='TORA.25.I1')
SCRN_25I_I1 = Monitor(eid='SCRN.25I.I1')
FCUP_25I_I1 = Monitor(eid='FCUP.25I.I1')
BPMG_25I_I1 = Monitor(eid='BPMG.25I.I1')
DCM_25_I1 = Monitor(eid='DCM.25.I1')
BPMC_38I_I1 = Monitor(eid='BPMC.38I.I1')
BPMR_38II_I1 = Monitor(eid='BPMR.38II.I1')
TORA_46_I1 = Monitor(eid='TORA.46.I1')
BAM_47_I1 = Monitor(eid='BAM.47.I1')
BPMF_47_I1 = Monitor(eid='BPMF.47.I1')
DCM_47_I1 = Monitor(eid='DCM.47.I1')
BPMF_48_I1 = Monitor(eid='BPMF.48.I1')
OTRL_48_I1 = Monitor(eid='OTRL.48.I1')
OTRL_50_I1 = Monitor(eid='OTRL.50.I1')
EOD_51_I1 = Monitor(eid='EOD.51.I1')
BPMF_52_I1 = Monitor(eid='BPMF.52.I1')
OTRC_55_I1 = Monitor(eid='OTRC.55.I1')
BPMA_55_I1 = Monitor(eid='BPMA.55.I1')
OTRC_56_I1 = Monitor(eid='OTRC.56.I1')
BPMA_57_I1 = Monitor(eid='BPMA.57.I1')
OTRC_58_I1 = Monitor(eid='OTRC.58.I1')
BPMA_59_I1 = Monitor(eid='BPMA.59.I1')
OTRC_59_I1 = Monitor(eid='OTRC.59.I1')
TORA_60_I1 = Monitor(eid='TORA.60.I1')
BPMATEST_60_I1 = Monitor(eid='BPMATEST.60.I1')
BPMATEST_61_I1 = Monitor(eid='BPMATEST.61.I1')
BPMA_63_I1 = Monitor(eid='BPMA.63.I1')
BPMA_72_I1 = Monitor(eid='BPMA.72.I1')
BPMA_75_I1 = Monitor(eid='BPMA.75.I1')
BPMA_77_I1 = Monitor(eid='BPMA.77.I1')
BPMA_80_I1 = Monitor(eid='BPMA.80.I1')
BPMA_82_I1 = Monitor(eid='BPMA.82.I1')
BPMA_85_I1 = Monitor(eid='BPMA.85.I1')
BPMA_87_I1 = Monitor(eid='BPMA.87.I1')
BPMA_90_I1 = Monitor(eid='BPMA.90.I1')
BPMA_92_I1 = Monitor(eid='BPMA.92.I1')
TORA_94_I1 = Monitor(eid='TORA.94.I1')
BPMF_95_I1 = Monitor(eid='BPMF.95.I1')
BPMS_99_I1 = Monitor(eid='BPMS.99.I1')
OTRS_99_I1 = Monitor(eid='OTRS.99.I1')
BPMF_103_I1 = Monitor(eid='BPMF.103.I1')
BPMA_103_I1 = Monitor(eid='BPMA.103.I1')
BPMA_105_I1 = Monitor(eid='BPMA.105.I1')
BPMA_107_I1 = Monitor(eid='BPMA.107.I1')
BPMA_110_I1 = Monitor(eid='BPMA.110.I1')
BPMA_112_I1 = Monitor(eid='BPMA.112.I1')
BPMA_115_I1 = Monitor(eid='BPMA.115.I1')
TORA_116_I1 = Monitor(eid='TORA.116.I1')
BPMA_117_I1 = Monitor(eid='BPMA.117.I1')
OTRA_118_I1 = Monitor(eid='OTRA.118.I1')
DCM_118_I1 = Monitor(eid='DCM.118.I1')
BPMA_119_I1 = Monitor(eid='BPMA.119.I1')

# sextupoles 
SC_74I_I1 = Sextupole(l=0.1, k2=-9.817522762, tilt=1.570796327, eid='SC.74I.I1')
SC_74II_I1 = Sextupole(l=0.1, k2=-5.948211334, tilt=1.570796327, eid='SC.74II.I1')
SC_76_I1 = Sextupole(l=0.1, k2=-5.948211334, tilt=1.570796327, eid='SC.76.I1')
SC_77_I1 = Sextupole(l=0.1, k2=-9.817522762, tilt=1.570796327, eid='SC.77.I1')
SC_79I_I1 = Sextupole(l=0.1, k2=-9.817522762, tilt=1.570796327, eid='SC.79I.I1')
SC_79II_I1 = Sextupole(l=0.1, k2=-5.948211334, tilt=1.570796327, eid='SC.79II.I1')
SC_81_I1 = Sextupole(l=0.1, k2=-5.948211334, tilt=1.570796327, eid='SC.81.I1')
SC_82_I1 = Sextupole(l=0.1, k2=-9.817522762, tilt=1.570796327, eid='SC.82.I1')
SC_84I_I1 = Sextupole(l=0.1, k2=9.817522762, tilt=1.570796327, eid='SC.84I.I1')
SC_84II_I1 = Sextupole(l=0.1, k2=5.948211334, tilt=1.570796327, eid='SC.84II.I1')
SC_86_I1 = Sextupole(l=0.1, k2=5.948211334, tilt=1.570796327, eid='SC.86.I1')
SC_87_I1 = Sextupole(l=0.1, k2=9.817522762, tilt=1.570796327, eid='SC.87.I1')
SC_89I_I1 = Sextupole(l=0.1, k2=9.817522762, tilt=1.570796327, eid='SC.89I.I1')
SC_89II_I1 = Sextupole(l=0.1, k2=5.948211334, tilt=1.570796327, eid='SC.89II.I1')
SC_91_I1 = Sextupole(l=0.1, k2=5.948211334, tilt=1.570796327, eid='SC.91.I1')
SC_92_I1 = Sextupole(l=0.1, k2=9.817522762, tilt=1.570796327, eid='SC.92.I1')

# octupole 

# undulator 
UNDU_49_I1 = Undulator(lperiod=0.074, nperiods=10, Kx=1.36*1.414213, Ky=0.0, eid='UNDU.49.I1')

# cavity 
C_A1_1_1_I1 = Cavity(l=1.0377, v=0.01815975, freq=1.3e9, phi=0.0, eid='C.A1.1.1.I1')
C_A1_1_2_I1 = Cavity(l=1.0377, v=0.01815975, freq=1.3e9, phi=0.0, eid='C.A1.1.2.I1')
C_A1_1_3_I1 = Cavity(l=1.0377, v=0.01815975, freq=1.3e9, phi=0.0, eid='C.A1.1.3.I1')
C_A1_1_4_I1 = Cavity(l=1.0377, v=0.01815975, freq=1.3e9, phi=0.0, eid='C.A1.1.4.I1')
C_A1_1_5_I1 = Cavity(l=1.0377, v=0.01815975, freq=1.3e9, phi=0.0, eid='C.A1.1.5.I1')
C_A1_1_6_I1 = Cavity(l=1.0377, v=0.01815975, freq=1.3e9, phi=0.0, eid='C.A1.1.6.I1')
C_A1_1_7_I1 = Cavity(l=1.0377, v=0.01815975, freq=1.3e9, phi=0.0, eid='C.A1.1.7.I1')
C_A1_1_8_I1 = Cavity(l=1.0377, v=0.01815975, freq=1.3e9, phi=0.0, eid='C.A1.1.8.I1')
C3_AH1_1_1_I1 = Cavity(l=0.346, v=0.0024999884, freq=3.9e9, phi=180.0, eid='C3.AH1.1.1.I1')
C3_AH1_1_2_I1 = Cavity(l=0.346, v=0.0024999884, freq=3.9e9, phi=180.0, eid='C3.AH1.1.2.I1')
C3_AH1_1_3_I1 = Cavity(l=0.346, v=0.0024999884, freq=3.9e9, phi=180.0, eid='C3.AH1.1.3.I1')
C3_AH1_1_4_I1 = Cavity(l=0.346, v=0.0024999884, freq=3.9e9, phi=180.0, eid='C3.AH1.1.4.I1')
C3_AH1_1_5_I1 = Cavity(l=0.346, v=0.0024999884, freq=3.9e9, phi=180.0, eid='C3.AH1.1.5.I1')
C3_AH1_1_6_I1 = Cavity(l=0.346, v=0.0024999884, freq=3.9e9, phi=180.0, eid='C3.AH1.1.6.I1')
C3_AH1_1_7_I1 = Cavity(l=0.346, v=0.0024999884, freq=3.9e9, phi=180.0, eid='C3.AH1.1.7.I1')
C3_AH1_1_8_I1 = Cavity(l=0.346, v=0.0024999884, freq=3.9e9, phi=180.0, eid='C3.AH1.1.8.I1')


TDSA_52_I1 = Cavity(l=0.7, v=0.0, freq=0.0, phi=0.0, eid='TDSA.52.I1')

# Matrices 

# Solenoids 
SOLB_23_I1 = Solenoid(l=0.0, k=0.0, eid='SOLB.23.I1')

und_start = Marker()
und_stop = Marker()
stop_A1 = Marker()


"""pytest fixtures defenition"""

@pytest.fixture(scope='module')
def cell():
    cell = (STSEC_23_I1, STSUB_23_I1, GUN_23_I1, D_1, SOLB_23_I1, D_2, CKX_23_I1, 
    CKY_23_I1, D_3, BPMG_24_I1, D_4, SCRN_24_I1, FCUP_24_I1, D_5, ENSUB_24_I1, 
    STSUB_24_I1, D_6, CKX_24_I1, CKY_24_I1, D_7, TORA_25_I1, D_8, SCRN_25I_I1, 
    FCUP_25I_I1, D_9, BPMG_25I_I1, D_10, DCM_25_I1, D_11, CKX_25_I1, CKY_25_I1, 
    D_12, START_25_I1, D_13, STAC_26_I1, D_14, C_A1_1_1_I1, D_15, C_A1_1_2_I1, 
    D_16, C_A1_1_3_I1, D_17, C_A1_1_4_I1, D_18, C_A1_1_5_I1, D_19, C_A1_1_6_I1, 
    D_20, C_A1_1_7_I1, D_21, C_A1_1_8_I1, D_22, stop_A1, Q_37_I1, CX_37_I1, CY_37_I1,
    D_23, BPMC_38I_I1, D_24, ENAC_38_I1, STAC_38_I1, D_25, BPMR_38II_I1, D_26, 
    Q_38_I1, CX_39_I1, CY_39_I1, D_27, C3_AH1_1_1_I1, D_28, C3_AH1_1_2_I1, D_29, 
    C3_AH1_1_3_I1, D_30, C3_AH1_1_4_I1, D_31, C3_AH1_1_5_I1, D_32, C3_AH1_1_6_I1, D_33, 
    C3_AH1_1_7_I1, D_34, C3_AH1_1_8_I1, D_35, ENAC_44_I1, D_36, END_45_I1, D_37, 
    TORA_46_I1, D_38, QI_46_I1, D_39, BAM_47_I1, D_40, BPMF_47_I1, D_41, 
    MPBPMF_47_I1, D_42, DCM_47_I1, D_43, QI_47_I1, D_44, STLAT_47_I1, D_45, 
    BL_48I_I1, D_46, BL_48II_I1, D_47, MPBPMF_48_I1, D_48, BPMF_48_I1, D_49, 
    OTRL_48_I1, D_50, und_start, UNDU_49_I1, und_stop, D_51, OTRL_50_I1, D_52, BL_50I_I1, D_53,
    BL_50II_I1, D_54, ENLAT_50_I1, D_55, QI_50_I1, D_56, EOD_51_I1, D_57, 
    CIY_51_I1, D_58, CIX_51_I1, D_59, BPMF_52_I1, D_60, MPBPMF_52_I1, D_61, 
    QI_52_I1, D_62, TDSA_52_I1, D_63, QI_53_I1, D_64, QI_54_I1, D_65, 
    START_55_I1, D_66, OTRC_55_I1, D_67, CIY_55_I1, D_68, BPMA_55_I1, D_69, 
    QI_55_I1, D_70, OTRC_56_I1, D_71, CIX_57_I1, D_72, BPMA_57_I1, D_73, 
    QI_57_I1, D_74, OTRC_58_I1, D_75, CIY_58_I1, D_76, QI_59_I1, D_77, 
    BPMA_59_I1, D_78, OTRC_59_I1, D_79, QI_60_I1, D_80, TORA_60_I1, D_81, 
    BPMATEST_60_I1, D_82, BPMATEST_61_I1, D_83, QI_61_I1, D_84, ENSUB_62_I1, D_85, 
    CIY_63_I1, D_86, QI_63_I1, D_87, BPMA_63_I1, D_88, CIX_65_I1, D_89, 
    QI_66_I1, D_90, QI_69_I1, D_91, QI_71_I1, D_92, BPMA_72_I1, D_93, 
    QI_72_I1, D_94, CIY_72_I1, D_95, CIX_73I_I1, D_96, START_73_I1, QI_73_I1, 
    D_97, BL_73_I1, CBL_73_I1, D_98, CIX_73II_I1, D_99, SC_74I_I1, D_100, 
    QI_74_I1, D_101, SC_74II_I1, D_102, BL_75_I1, D_103, BPMA_75_I1, D_104, 
    CIY_75_I1, D_105, QI_75_I1, D_106, CIX_76_I1, D_107, BL_76_I1, D_108, 
    SC_76_I1, D_109, QI_77_I1, D_110, SC_77_I1, D_111, BPMA_77_I1, D_112, 
    BL_77_I1, D_113, QI_78_I1, D_114, BL_78_I1, CBL_78_I1, D_115, CIX_78_I1, 
    D_116, SC_79I_I1, D_117, QI_79_I1, D_118, SC_79II_I1, D_119, BL_80_I1, 
    D_120, BPMA_80_I1, D_121, CIY_80_I1, D_122, QI_80_I1, D_123, CIX_81_I1, 
    D_124, BL_81_I1, D_125, SC_81_I1, D_126, QI_82_I1, D_127, SC_82_I1, 
    D_128, BPMA_82_I1, D_129, BL_82_I1, D_130, QI_83_I1, D_131, BL_83_I1, 
    D_132, CBL_83_I1, D_133, CIX_83_I1, D_134, SC_84I_I1, D_135, QI_84_I1, 
    D_136, SC_84II_I1, D_137, BL_85_I1, D_138, BPMA_85_I1, D_139, CIY_85_I1, 
    D_140, QI_85_I1, D_141, CIX_86_I1, D_142, BL_86_I1, D_143, SC_86_I1, 
    D_144, QI_86_I1, D_145, SC_87_I1, D_146, BPMA_87_I1, D_147, BL_87_I1, 
    D_148, QI_88_I1, D_149, BL_88_I1, CBL_88_I1, D_150, CIX_88_I1, D_151, 
    SC_89I_I1, D_152, QI_89_I1, D_153, SC_89II_I1, D_154, BL_90_I1, D_155, 
    CBL_90_I1, D_156, BPMA_90_I1, D_157, CIX_90_I1, D_158, QI_90_I1, D_159, 
    BL_91_I1, D_160, SC_91_I1, D_161, QI_92_I1, D_162, SC_92_I1, D_163, 
    CIY_92_I1, D_164, BPMA_92_I1, D_165, BL_92_I1, D_166, QI_93_I1, END_93_I1, 
    D_167, TORA_94_I1, D_168, CIY_94_I1, D_169, QI_94_I1, D_170, BPMF_95_I1, 
    D_171, MPBPMF_95_I1, D_172, CIX_95_I1, D_173, QI_95_I1, D_174, START_96_I1, 
    D_175, BB_96_I1, D_176, BB_98_I1, D_177, CBB_98_I1, D_178, BPMS_99_I1, 
    D_179, OTRS_99_I1, D_180, BB_100_I1, CBB_100_I1, D_181, BB_101_I1, D_182, 
    CBB_101_I1, D_183, END_101_I1, D_184, QI_102_I1, D_185, CIX_102_I1, D_186, 
    BPMF_103_I1, D_187, MPBPMF_103_I1, D_188, CIY_103_I1, D_189, QI_103_I1, D_190, 
    BPMA_103_I1, D_191, CIX_104_I1, D_192, QI_104_I1, START_104_I1, D_193, BPMA_105_I1, 
    D_194, CIY_107_I1, D_195, QI_107_I1, D_196, BPMA_107_I1, D_197, CIX_109_I1, 
    D_198, QI_109_I1, D_199, BPMA_110_I1, D_200, CIY_112_I1, D_201, QI_112_I1, 
    D_202, BPMA_112_I1, D_203, CIX_114_I1, D_204, QI_114_I1, D_205, BPMA_115_I1, 
    D_206, TORA_116_I1, D_207, CIY_116_I1, D_208, QI_116_I1, D_209, BPMA_117_I1, 
    D_210, OTRA_118_I1, D_211, DCM_118_I1, D_212, CIX_118_I1, D_213, QI_118_I1, 
    D_214, BPMA_119_I1, D_215, END_119_I1)
    
    return cell
    
    
@pytest.fixture(scope='module')
def method():

    m = MethodTM()
    m.global_method = SecondTM
    
    return m
    
    
@pytest.fixture(scope='module')
def lattice(cell, method):
    return MagneticLattice(cell, method=method)

    
@pytest.fixture(scope='function')
def tws0():
    
    t = Twiss()
    t.beta_x = 29.171
    t.beta_y = 29.171
    t.alpha_x = 10.955
    t.alpha_y = 10.955
    t.E = 0.005
    
    return t
