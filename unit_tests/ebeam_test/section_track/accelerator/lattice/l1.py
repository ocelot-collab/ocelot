from ocelot import * 
tws = Twiss()
tws.beta_x  = 3.203506642234328
tws.beta_y  = 4.580189307026297
tws.alpha_x = 0.3111465345215999
tws.alpha_y = -1.7973623618150731

tws.beta_x  = 2.6096907242276925
tws.beta_y  = 7.150678422205259
tws.alpha_x = 0.22820424990918614
tws.alpha_y = -2.165836718254254
tws.s = 62.089004999999936


tws.E = 0.13
# Drifts
d_1 = Drift(l=0.5, eid='D_1')
d_2 = Drift(l=0.74, eid='D_2')
d_3 = Drift(l=0.08115, eid='D_3')
d_4 = Drift(l=0.13115, eid='D_4')
d_5 = Drift(l=1.45, eid='D_5')
d_6 = Drift(l=0.63115, eid='D_6')
d_7 = Drift(l=2.8623, eid='D_7')
d_8 = Drift(l=1.7623, eid='D_8')
d_9 = Drift(l=0.68115, eid='D_9')
d_10 = Drift(l=0.18115, eid='D_10')
d_12 = Drift(l=0.2, eid='D_12')
d_14 = Drift(l=0.182428, eid='D_14')
d_15 = Drift(l=5.1e-05, eid='D_15')
d_16 = Drift(l=0.135, eid='D_16')
d_17 = Drift(l=0.239432, eid='D_17')
d_18 = Drift(l=0.0751, eid='D_18')
d_20 = Drift(l=0.13494, eid='D_20')
d_21 = Drift(l=0.158616, eid='D_21')
d_22 = Drift(l=0.1325, eid='D_22')
d_23 = Drift(l=0.09865, eid='D_23')
d_24 = Drift(l=0.11615, eid='D_24')
d_25 = Drift(l=0.273615, eid='D_25')
d_29 = Drift(l=0.324432, eid='D_29')
d_30 = Drift(l=0.150051, eid='D_30')
d_32 = Drift(l=0.182429, eid='D_32')
d_38 = Drift(l=0.134939, eid='D_38')
d_43 = Drift(l=0.273616, eid='D_43')
d_47 = Drift(l=0.324431, eid='D_47')
d_48 = Drift(l=0.150052, eid='D_48')
d_51 = Drift(l=5.2e-05, eid='D_51')
d_53 = Drift(l=0.239431, eid='D_53')
d_75 = Drift(l=8e-06, eid='D_75')
d_76 = Drift(l=0.158608, eid='D_76')
d_79 = Drift(l=0.489766, eid='D_79')
d_83 = Drift(l=0.091932, eid='D_83')
d_87 = Drift(l=0.88128, eid='D_87')
d_89 = Drift(l=0.09115, eid='D_89')
d_90 = Drift(l=0.37715, eid='D_90')
d_91 = Drift(l=0.204, eid='D_91')
d_93 = Drift(l=0.481747, eid='D_93')
d_94 = Drift(l=1.008384, eid='D_94')
d_95 = Drift(l=0.000597, eid='D_95')
d_96 = Drift(l=0.865, eid='D_96')
d_97 = Drift(l=0.31, eid='D_97')
d_98 = Drift(l=0.325597, eid='D_98')
d_100 = Drift(l=1.007788, eid='D_100')
d_101 = Drift(l=0.000596, eid='D_101')
d_102 = Drift(l=0.49115, eid='D_102')
d_104 = Drift(l=0.346, eid='D_104')
d_105 = Drift(l=0.144, eid='D_105')
d_108 = Drift(l=0.56, eid='D_108')
d_111 = Drift(l=1.89, eid='D_111')
d_120 = Drift(l=2.04, eid='D_120')
d_123 = Drift(l=1.1242, eid='D_123')
d_124 = Drift(l=0.1408, eid='D_124')
d_126 = Drift(l=0.13965, eid='D_126')
d_127 = Drift(l=1.09132, eid='D_127')
d_128 = Drift(l=0.44018, eid='D_128')
d_130 = Drift(l=0.14395, eid='D_130')
d_131 = Drift(l=3.923801, eid='D_131')
d_132 = Drift(l=0.345899, eid='D_132')
d_133 = Drift(l=0.3459, eid='D_133')
d_138 = Drift(l=0.345901, eid='D_138')
d_139 = Drift(l=0.247499, eid='D_139')
d_140 = Drift(l=0.043201, eid='D_140')
d_141 = Drift(l=0.084999, eid='D_141')
d_142 = Drift(l=0.679501, eid='D_142')
d_150 = Drift(l=0.2475, eid='D_150')
d_151 = Drift(l=0.0432, eid='D_151')
d_152 = Drift(l=0.085001, eid='D_152')
d_153 = Drift(l=0.679499, eid='D_153')
d_162 = Drift(l=0.043199, eid='D_162')
d_163 = Drift(l=0.085, eid='D_163')
d_175 = Drift(l=3.7079, eid='D_175')
d_176 = Drift(l=1.545999, eid='D_176')
d_177 = Drift(l=0.248001, eid='D_177')
d_178 = Drift(l=0.15215, eid='D_178')
d_179 = Drift(l=0.082149, eid='D_179')
d_180 = Drift(l=1.215001, eid='D_180')
d_181 = Drift(l=1.613, eid='D_181')
d_182 = Drift(l=0.15265, eid='D_182')
d_183 = Drift(l=0.131649, eid='D_183')
d_184 = Drift(l=0.650001, eid='D_184')
d_185 = Drift(l=0.231649, eid='D_185')
d_186 = Drift(l=0.10165, eid='D_186')
d_187 = Drift(l=0.4515, eid='D_187')
d_188 = Drift(l=0.428605, eid='D_188')
d_189 = Drift(l=8.510829, eid='D_189')
d_190 = Drift(l=0.000104, eid='D_190')
d_193 = Drift(l=0.325104, eid='D_193')
d_195 = Drift(l=8.510725, eid='D_195')
d_197 = Drift(l=0.1, eid='D_197')
d_198 = Drift(l=0.579, eid='D_198')
d_199 = Drift(l=0.4508, eid='D_199')
d_200 = Drift(l=0.1542, eid='D_200')
d_201 = Drift(l=0.09715, eid='D_201')
d_202 = Drift(l=0.66515, eid='D_202')
d_204 = Drift(l=0.62315, eid='D_204')
d_205 = Drift(l=0.536, eid='D_205')
d_206 = Drift(l=0.15315, eid='D_206')
d_210 = Drift(l=0.53115, eid='D_210')
d_211 = Drift(l=0.379, eid='D_211')
d_213 = Drift(l=0.13165, eid='D_213')
d_214 = Drift(l=1.03115, eid='D_214')
d_215 = Drift(l=1.23015, eid='D_215')
d_216 = Drift(l=0.18, eid='D_216')
d_219 = Drift(l=1.329, eid='D_219')
d_222 = Drift(l=1.066, eid='D_222')
d_223 = Drift(l=0.147, eid='D_223')
d_225 = Drift(l=0.83115, eid='D_225')
d_226 = Drift(l=0.699, eid='D_226')
d_227 = Drift(l=0.13265, eid='D_227')
d_228 = Drift(l=0.83165, eid='D_228')
d_229 = Drift(l=0.679, eid='D_229')
d_231 = Drift(l=0.10765, eid='D_231')
d_232 = Drift(l=0.624, eid='D_232')
d_235 = Drift(l=0.18165, eid='D_235')
d_236 = Drift(l=0.55, eid='D_236')
d_237 = Drift(l=0.34615, eid='D_237')
d_238 = Drift(l=1.22815, eid='D_238')
d_239 = Drift(l=0.152151, eid='D_239')
d_240 = Drift(l=0.081151, eid='D_240')
d_241 = Drift(l=0.15, eid='D_241')
d_244 = Drift(l=1.494401, eid='D_244')

# Quadrupoles
qi_63_i1 = Quadrupole(l=0.2377, k1=-1.824087238956668, eid='QI.63.I1')
qi_66_i1 = Quadrupole(l=0.2377, k1=2.1425757568363486, eid='QI.66.I1')
qi_69_i1 = Quadrupole(l=0.2377, k1=-2.1265116390408076, eid='QI.69.I1')
qi_71_i1 = Quadrupole(l=0.2377, k1=3.4033787980647876, eid='QI.71.I1')
qi_72_i1 = Quadrupole(l=0.2377, k1=-4.434528949936896, eid='QI.72.I1')
qi_73_i1 = Quadrupole(l=0.2377, k1=4.632555992006731, eid='QI.73.I1')
qi_74_i1 = Quadrupole(l=0.2377, k1=-4.991375220025242, eid='QI.74.I1')
qi_75_i1 = Quadrupole(l=0.2377, k1=5.027338700883466, eid='QI.75.I1')
qi_77_i1 = Quadrupole(l=0.2377, k1=-4.991375220025242, eid='QI.77.I1')
qi_78_i1 = Quadrupole(l=0.2377, k1=4.632555992006731, eid='QI.78.I1')
qi_79_i1 = Quadrupole(l=0.2377, k1=-4.991375220025242, eid='QI.79.I1')
qi_80_i1 = Quadrupole(l=0.2377, k1=5.027338700883466, eid='QI.80.I1')
qi_82_i1 = Quadrupole(l=0.2377, k1=-4.991375220025242, eid='QI.82.I1')
qi_83_i1 = Quadrupole(l=0.2377, k1=4.632555992006731, eid='QI.83.I1')
qi_84_i1 = Quadrupole(l=0.2377, k1=-4.991375220025242, eid='QI.84.I1')
qi_85_i1 = Quadrupole(l=0.2377, k1=5.027338700883466, eid='QI.85.I1')
qi_86_i1 = Quadrupole(l=0.2377, k1=-4.991375220025242, eid='QI.86.I1')
qi_88_i1 = Quadrupole(l=0.2377, k1=4.632555992006731, eid='QI.88.I1')
qi_89_i1 = Quadrupole(l=0.2377, k1=-4.991375220025242, eid='QI.89.I1')
qi_90_i1 = Quadrupole(l=0.2377, k1=5.027338700883466, eid='QI.90.I1')
qi_92_i1 = Quadrupole(l=0.2377, k1=-4.991375220025242, eid='QI.92.I1')
qi_93_i1 = Quadrupole(l=0.2377, k1=-0.7125210100967607, eid='QI.93.I1')
qi_94_i1 = Quadrupole(l=0.2377, k1=3.345427061842659, eid='QI.94.I1')
qi_95_i1 = Quadrupole(l=0.2377, k1=-3.013328963819941, eid='QI.95.I1')
qi_102_i1 = Quadrupole(l=0.2377, k1=0.38443867522086667, eid='QI.102.I1')
qi_103_i1 = Quadrupole(l=0.2377, k1=-1.0684325170382836, eid='QI.103.I1')
qi_104_i1 = Quadrupole(l=0.2377, k1=1.4759984968447624, eid='QI.104.I1')
qi_107_i1 = Quadrupole(l=0.2377, k1=-1.522641254101809, eid='QI.107.I1')
qi_109_i1 = Quadrupole(l=0.2377, k1=1.5226410681531344, eid='QI.109.I1')
qi_112_i1 = Quadrupole(l=0.2377, k1=-1.522641254101809, eid='QI.112.I1')
qi_114_i1 = Quadrupole(l=0.2377, k1=0.9967996613378207, eid='QI.114.I1')
qi_116_i1 = Quadrupole(l=0.2377, k1=0.5375414030290282, eid='QI.116.I1')
qi_118_i1 = Quadrupole(l=0.2377, k1=-0.9412199200673117, eid='QI.118.I1')
q_134_l1 = Quadrupole(l=0.2136, k1=0.2815559414794007, eid='Q.134.L1')
q_146_l1 = Quadrupole(l=0.2136, k1=-0.299064888576779, eid='Q.146.L1')
q_158_l1 = Quadrupole(l=0.2136, k1=0.21347748782771533, eid='Q.158.L1')
q_170_l1 = Quadrupole(l=0.2136, k1=-0.22026272752808987, eid='Q.170.L1')
qi_176_b1 = Quadrupole(l=0.2377, eid='QI.176.B1')
qd_179_b1 = Quadrupole(l=0.2367, k1=0.7946818453738911, eid='QD.179.B1')
qd_181_b1 = Quadrupole(l=0.2367, k1=-0.7289505606252641, eid='QD.181.B1')
qi_204_b1 = Quadrupole(l=0.2377, k1=-0.9082538216238958, eid='QI.204.B1')
qi_205_b1 = Quadrupole(l=0.2377, k1=0.015226619267984857, eid='QI.205.B1')
qi_206_b1 = Quadrupole(l=0.2377, k1=0.691130037442154, eid='QI.206.B1')
qi_209_b1 = Quadrupole(l=0.2377, k1=1.0443862490534286, eid='QI.209.B1')
qd_210_b1 = Quadrupole(l=0.2367, k1=-2.1831505251373047, eid='QD.210.B1')
qi_211_b1 = Quadrupole(l=0.2377, k1=0.6386118847286495, eid='QI.211.B1')
qi_213_b1 = Quadrupole(l=0.2377, k1=1.1869696491375683, eid='QI.213.B1')
qi_215_b1 = Quadrupole(l=0.2377, k1=-1.1237338750525874, eid='QI.215.B1')
qi_217_b1 = Quadrupole(l=0.2377, k1=-1.4350758199411022, eid='QI.217.B1')
qd_219_b1 = Quadrupole(l=0.2367, k1=2.8595907038445287, eid='QD.219.B1')
qd_221_b1 = Quadrupole(l=0.2367, k1=-2.8595909290240815, eid='QD.221.B1')
qd_223_b1 = Quadrupole(l=0.2367, k1=2.8595907038445287, eid='QD.223.B1')
qi_224_b1 = Quadrupole(l=0.2377, k1=-0.8565833407656711, eid='QI.224.B1')
qi_226_b1 = Quadrupole(l=0.2377, k1=-1.5618567639882204, eid='QI.226.B1')
qi_227_b1 = Quadrupole(l=0.2377, k1=1.523532902818679, eid='QI.227.B1')

# SBends
bl_73_i1 = SBend(l=0.2, angle=-0.1109740393, e1=-0.05548702, e2=-0.05548702, tilt=1.570796327, eid='BL.73.I1')
bl_75_i1 = SBend(l=0.2, angle=0.0426524581, e1=0.021326229, e2=0.021326229, tilt=1.570796327, eid='BL.75.I1')
bl_76_i1 = SBend(l=0.2, angle=0.0426524581, e1=0.021326229, e2=0.021326229, tilt=1.570796327, eid='BL.76.I1')
bl_77_i1 = SBend(l=0.2, angle=-0.1109740393, e1=-0.05548702, e2=-0.05548702, tilt=1.570796327, eid='BL.77.I1')
bl_78_i1 = SBend(l=0.2, angle=-0.1109740393, e1=-0.05548702, e2=-0.05548702, tilt=1.570796327, eid='BL.78.I1')
bl_80_i1 = SBend(l=0.2, angle=0.0426524581, e1=0.021326229, e2=0.021326229, tilt=1.570796327, eid='BL.80.I1')
bl_81_i1 = SBend(l=0.2, angle=0.0426524581, e1=0.021326229, e2=0.021326229, tilt=1.570796327, eid='BL.81.I1')
bl_82_i1 = SBend(l=0.2, angle=-0.1109740393, e1=-0.05548702, e2=-0.05548702, tilt=1.570796327, eid='BL.82.I1')
bl_83_i1 = SBend(l=0.2, angle=0.1109740393, e1=0.05548702, e2=0.05548702, tilt=1.570796327, eid='BL.83.I1')
bl_85_i1 = SBend(l=0.2, angle=-0.0426524581, e1=-0.021326229, e2=-0.021326229, tilt=1.570796327, eid='BL.85.I1')
bl_86_i1 = SBend(l=0.2, angle=-0.0426524581, e1=-0.021326229, e2=-0.021326229, tilt=1.570796327, eid='BL.86.I1')
bl_87_i1 = SBend(l=0.2, angle=0.1109740393, e1=0.05548702, e2=0.05548702, tilt=1.570796327, eid='BL.87.I1')
bl_88_i1 = SBend(l=0.2, angle=0.1109740393, e1=0.05548702, e2=0.05548702, tilt=1.570796327, eid='BL.88.I1')
bl_90_i1 = SBend(l=0.2, angle=-0.0426524581, e1=-0.021326229, e2=-0.021326229, tilt=1.570796327, eid='BL.90.I1')
bl_91_i1 = SBend(l=0.2, angle=-0.0426524581, e1=-0.021326229, e2=-0.021326229, tilt=1.570796327, eid='BL.91.I1')
bl_92_i1 = SBend(l=0.2, angle=0.1109740393, e1=0.05548702, e2=0.05548702, tilt=1.570796327, eid='BL.92.I1')
bb_96_i1 = SBend(l=0.5, angle=0.1366592804, e2=0.1366592804, tilt=1.570796327, eid='BB.96.I1')
bb_98_i1 = SBend(l=0.5, angle=-0.1366592804, e1=-0.1366592804, tilt=1.570796327, eid='BB.98.I1')
bb_100_i1 = SBend(l=0.5, angle=-0.1366592804, e2=-0.1366592804, tilt=1.570796327, eid='BB.100.I1')
bb_101_i1 = SBend(l=0.5, angle=0.1366592804, e1=0.1366592804, tilt=1.570796327, eid='BB.101.I1')
bb_182_b1 = SBend(l=0.5, angle=0.0532325422, e2=0.0532325422, tilt=1.570796327, eid='BB.182.B1')
bb_191_b1 = SBend(l=0.5, angle=-0.0532325422, e1=-0.0532325422, tilt=1.570796327, eid='BB.191.B1')
bb_193_b1 = SBend(l=0.5, angle=-0.0532325422, e2=-0.0532325422, tilt=1.570796327, eid='BB.193.B1')
bb_202_b1 = SBend(l=0.5, angle=0.0532325422, e1=0.0532325422, tilt=1.570796327, eid='BB.202.B1')

# Sextupoles
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

# Hcors
cbb_62_i1d = Hcor(eid='CBB.62.I1D')
cix_65_i1 = Hcor(l=0.1, eid='CIX.65.I1')
cix_73i_i1 = Hcor(l=0.1, eid='CIX.73I.I1')
cix_73ii_i1 = Hcor(l=0.1, eid='CIX.73II.I1')
cix_76_i1 = Hcor(l=0.1, eid='CIX.76.I1')
cix_78_i1 = Hcor(l=0.1, eid='CIX.78.I1')
cix_81_i1 = Hcor(l=0.1, eid='CIX.81.I1')
cix_83_i1 = Hcor(l=0.1, eid='CIX.83.I1')
cix_86_i1 = Hcor(l=0.1, eid='CIX.86.I1')
cix_88_i1 = Hcor(l=0.1, eid='CIX.88.I1')
cix_90_i1 = Hcor(l=0.1, eid='CIX.90.I1')
cix_95_i1 = Hcor(l=0.1, eid='CIX.95.I1')
cix_102_i1 = Hcor(l=0.1, eid='CIX.102.I1')
cix_104_i1 = Hcor(l=0.1, eid='CIX.104.I1')
cix_109_i1 = Hcor(l=0.1, eid='CIX.109.I1')
cix_114_i1 = Hcor(l=0.1, eid='CIX.114.I1')
cix_118_i1 = Hcor(l=0.1, eid='CIX.118.I1')
cx_134_l1 = Hcor(eid='CX.134.L1')
cx_146_l1 = Hcor(eid='CX.146.L1')
cx_158_l1 = Hcor(eid='CX.158.L1')
cx_170_l1 = Hcor(eid='CX.170.L1')
cix_177_b1 = Hcor(l=0.1, eid='CIX.177.B1')
ccx_179_b1 = Hcor(l=0.1, eid='CCX.179.B1')
cix_205_b1 = Hcor(l=0.1, eid='CIX.205.B1')
cix_209_b1 = Hcor(l=0.1, eid='CIX.209.B1')
cix_213_b1 = Hcor(l=0.1, eid='CIX.213.B1')
cix_216_b1 = Hcor(l=0.1, eid='CIX.216.B1')
cfx_223_b1 = Hcor(l=0.1, eid='CFX.223.B1')
cfx_226_b1 = Hcor(l=0.1, eid='CFX.226.B1')

# Vcors
ciy_63_i1 = Vcor(l=0.1, eid='CIY.63.I1')
ciy_72_i1 = Vcor(l=0.1, eid='CIY.72.I1')
cbl_73_i1 = Vcor(eid='CBL.73.I1')
ciy_75_i1 = Vcor(l=0.1, eid='CIY.75.I1')
cbl_78_i1 = Vcor(eid='CBL.78.I1')
ciy_80_i1 = Vcor(l=0.1, eid='CIY.80.I1')
cbl_83_i1 = Vcor(eid='CBL.83.I1')
ciy_85_i1 = Vcor(l=0.1, eid='CIY.85.I1')
cbl_88_i1 = Vcor(eid='CBL.88.I1')
cbl_90_i1 = Vcor(eid='CBL.90.I1')
ciy_92_i1 = Vcor(l=0.1, eid='CIY.92.I1')
ciy_94_i1 = Vcor(l=0.1, eid='CIY.94.I1')
cbb_98_i1 = Vcor(eid='CBB.98.I1')
cbb_100_i1 = Vcor(eid='CBB.100.I1')
cbb_101_i1 = Vcor(eid='CBB.101.I1')
ciy_103_i1 = Vcor(l=0.1, eid='CIY.103.I1')
ciy_107_i1 = Vcor(l=0.1, eid='CIY.107.I1')
ciy_112_i1 = Vcor(l=0.1, eid='CIY.112.I1')
ciy_116_i1 = Vcor(l=0.1, eid='CIY.116.I1')
cy_134_l1 = Vcor(eid='CY.134.L1')
cy_146_l1 = Vcor(eid='CY.146.L1')
cy_158_l1 = Vcor(eid='CY.158.L1')
cy_170_l1 = Vcor(eid='CY.170.L1')
ciy_176_b1 = Vcor(l=0.1, eid='CIY.176.B1')
ccy_181_b1 = Vcor(l=0.1, eid='CCY.181.B1')
cbb_191_b1 = Vcor(eid='CBB.191.B1')
cbb_193_b1 = Vcor(eid='CBB.193.B1')
cbb_202_b1 = Vcor(eid='CBB.202.B1')
ciy_204_b1 = Vcor(l=0.1, eid='CIY.204.B1')
ccy_210_b1 = Vcor(l=0.1, eid='CCY.210.B1')
ciy_214_b1 = Vcor(l=0.1, eid='CIY.214.B1')
ccy_217_b1 = Vcor(l=0.1, eid='CCY.217.B1')
ccy_221_b1 = Vcor(l=0.1, eid='CCY.221.B1')
ciy_226_b1 = Vcor(l=0.1, eid='CIY.226.B1')

# Cavitys
dip_kick = 0
c_a2_1_1_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.1.1.L1')
c_a2_1_2_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.1.2.L1')
c_a2_1_3_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.1.3.L1')
c_a2_1_4_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.1.4.L1')
c_a2_1_5_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.1.5.L1')
c_a2_1_6_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.1.6.L1')
c_a2_1_7_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.1.7.L1')
c_a2_1_8_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.1.8.L1')
c_a2_2_1_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.2.1.L1')
c_a2_2_2_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.2.2.L1')
c_a2_2_3_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.2.3.L1')
c_a2_2_4_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.2.4.L1')
c_a2_2_5_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.2.5.L1')
c_a2_2_6_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.2.6.L1')
c_a2_2_7_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.2.7.L1')
c_a2_2_8_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.2.8.L1')
c_a2_3_1_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.3.1.L1')
c_a2_3_2_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.3.2.L1')
c_a2_3_3_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.3.3.L1')
c_a2_3_4_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.3.4.L1')
c_a2_3_5_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.3.5.L1')
c_a2_3_6_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.3.6.L1')
c_a2_3_7_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.3.7.L1')
c_a2_3_8_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.3.8.L1')
c_a2_4_1_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.4.1.L1')
c_a2_4_2_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.4.2.L1')
c_a2_4_3_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.4.3.L1')
c_a2_4_4_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.4.4.L1')
c_a2_4_5_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.4.5.L1')
c_a2_4_6_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.4.6.L1')
c_a2_4_7_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.4.7.L1')
c_a2_4_8_l1 = Cavity(l=1.0377, v=0.0196539, freq=1300000000.0, phi=25, vx_up=dip_kick * (-5.6813e-05 + 1.0751e-05j), vy_up=dip_kick * (-4.1091e-05 + 5.739e-07j), vxx_up=(0.00099943 - 0.00081401j), vxy_up=(0.0034065 - 0.0004146j), vx_down=dip_kick * (-2.4014e-05 + 1.2492e-05j), vy_down=dip_kick * (3.6481e-05 + 7.9888e-06j), vxx_down=(-0.004057 - 0.0001369j), vxy_down=(0.0029243 - 1.2891e-05j), eid='C.A2.4.8.L1')
tdsb_208_b1 = TDCavity(l=1.5, freq=2800000000.0, eid='TDSB.208.B1')

# Monitors
bpma_63_i1 = Monitor(eid='BPMA.63.I1')
bpma_72_i1 = Monitor(eid='BPMA.72.I1')
bpma_75_i1 = Monitor(eid='BPMA.75.I1')
bpma_77_i1 = Monitor(eid='BPMA.77.I1')
bpma_80_i1 = Monitor(eid='BPMA.80.I1')
bpma_82_i1 = Monitor(eid='BPMA.82.I1')
bpma_85_i1 = Monitor(eid='BPMA.85.I1')
bpma_87_i1 = Monitor(eid='BPMA.87.I1')
bpma_90_i1 = Monitor(eid='BPMA.90.I1')
bpma_92_i1 = Monitor(eid='BPMA.92.I1')
bpmf_95_i1 = Monitor(eid='BPMF.95.I1')
bpms_99_i1 = Monitor(eid='BPMS.99.I1')
bpmf_103_i1 = Monitor(eid='BPMF.103.I1')
bpma_103_i1 = Monitor(eid='BPMA.103.I1')
bpma_105_i1 = Monitor(eid='BPMA.105.I1')
bpma_107_i1 = Monitor(eid='BPMA.107.I1')
bpma_110_i1 = Monitor(eid='BPMA.110.I1')
bpma_112_i1 = Monitor(eid='BPMA.112.I1')
bpma_115_i1 = Monitor(eid='BPMA.115.I1')
bpma_117_i1 = Monitor(eid='BPMA.117.I1')
bpma_119_i1 = Monitor(eid='BPMA.119.I1')
bpmc_134_l1 = Monitor(eid='BPMC.134.L1')
bpmr_146_l1 = Monitor(eid='BPMR.146.L1')
bpmc_158_l1 = Monitor(eid='BPMC.158.L1')
bpmr_170_l1 = Monitor(eid='BPMR.170.L1')
bpma_175_b1 = Monitor(eid='BPMA.175.B1')
bpma_179_b1 = Monitor(eid='BPMA.179.B1')
bpmf_181_b1 = Monitor(eid='BPMF.181.B1')
bpms_192_b1 = Monitor(eid='BPMS.192.B1')
bpmf_203_b1 = Monitor(eid='BPMF.203.B1')
bpma_206_b1 = Monitor(eid='BPMA.206.B1')
bpma_210_b1 = Monitor(eid='BPMA.210.B1')
bpma_213_b1 = Monitor(eid='BPMA.213.B1')
bpma_215_b1 = Monitor(eid='BPMA.215.B1')
bpma_217_b1 = Monitor(eid='BPMA.217.B1')
bpma_219_b1 = Monitor(eid='BPMA.219.B1')
bpma_221_b1 = Monitor(eid='BPMA.221.B1')
bpma_223_b1 = Monitor(eid='BPMA.223.B1')
bpma_226_b1 = Monitor(eid='BPMA.226.B1')
bpma_227_b1 = Monitor(eid='BPMA.227.B1')

# Markers
match_73_i1 = Marker(eid='MATCH.73.I1')
tora_94_i1 = Marker(eid='TORA.94.I1')
otrs_99_i1 = Marker(eid='OTRS.99.I1')
match_104_i1 = Marker(eid='MATCH.104.I1')
tora_116_i1 = Marker(eid='TORA.116.I1')
otra_118_i1 = Marker(eid='OTRA.118.I1')
match_174_l1 = Marker(eid='MATCH.174.L1')
tora_175_b1 = Marker(eid='TORA.175.B1')
otra_180_b1 = Marker(eid='OTRA.180.B1')
otrs_192_b1 = Marker(eid='OTRS.192.B1')
match_202_b1 = Marker(eid='MATCH.202.B1')
tora_203_b1 = Marker(eid='TORA.203.B1')
otra_206_b1 = Marker(eid='OTRA.206.B1')
match_207_b1 = Marker(eid='MATCH.207.B1')
match_218_b1 = Marker(eid='MATCH.218.B1')
otrb_218_b1 = Marker(eid='OTRB.218.B1')
otrb_220_b1 = Marker(eid='OTRB.220.B1')
otrb_222_b1 = Marker(eid='OTRB.222.B1')
otrb_224_b1 = Marker(eid='OTRB.224.B1')
ensub_229_b1 = Marker(eid='ENSUB.229.B1')

dl_start = Marker()
stlat_96_i1 = Marker()
enlat_101_i1 = Marker()

d_188_2 = Drift(l=0.100105, eid='D_208')
stlat_182_b1 = Marker()
d_188_1 = Drift(l=d_188.l - d_188_2.l, eid='D_188')
d_188_n = (d_188_1, stlat_182_b1, d_188_2)



# Lattice 
cell = (d_1, cbb_62_i1d, d_2, ciy_63_i1, d_3, qi_63_i1, d_4, bpma_63_i1, d_5, 
cix_65_i1, d_6, qi_66_i1, d_7, qi_69_i1, d_8, qi_71_i1, d_9, bpma_72_i1, d_10, 
qi_72_i1, d_10, ciy_72_i1, d_12, cix_73i_i1, d_3, match_73_i1, qi_73_i1,dl_start, d_14, bl_73_i1,
d_15, cbl_73_i1, d_16, cix_73ii_i1, d_17, sc_74i_i1, d_18, qi_74_i1, d_18, sc_74ii_i1, 
d_20, bl_75_i1, d_21, bpma_75_i1, d_22, ciy_75_i1, d_23, qi_75_i1, d_24, cix_76_i1, 
d_25, bl_76_i1, d_20, sc_76_i1, d_18, qi_77_i1, d_18, sc_77_i1, d_29, bpma_77_i1, 
d_30, bl_77_i1, d_14, qi_78_i1, d_32, bl_78_i1, d_15, cbl_78_i1, d_16, cix_78_i1, 
d_17, sc_79i_i1, d_18, qi_79_i1, d_18, sc_79ii_i1, d_38, bl_80_i1, d_21, bpma_80_i1, 
d_22, ciy_80_i1, d_23, qi_80_i1, d_24, cix_81_i1, d_43, bl_81_i1, d_20, sc_81_i1, 
d_18, qi_82_i1, d_18, sc_82_i1, d_47, bpma_82_i1, d_48, bl_82_i1, d_14, qi_83_i1, 
d_14, bl_83_i1, d_51, cbl_83_i1, d_16, cix_83_i1, d_53, sc_84i_i1, d_18, qi_84_i1, 
d_18, sc_84ii_i1, d_20, bl_85_i1, d_21, bpma_85_i1, d_22, ciy_85_i1, d_23, qi_85_i1, 
d_24, cix_86_i1, d_43, bl_86_i1, d_38, sc_86_i1, d_18, qi_86_i1, d_18, sc_87_i1, 
d_29, bpma_87_i1, d_48, bl_87_i1, d_14, qi_88_i1, d_14, bl_88_i1, d_15, cbl_88_i1, 
d_16, cix_88_i1, d_17, sc_89i_i1, d_18, qi_89_i1, d_18, sc_89ii_i1, d_38, bl_90_i1, 
d_75, cbl_90_i1, d_76, bpma_90_i1, d_22, cix_90_i1, d_23, qi_90_i1, d_79, bl_91_i1, 
d_20, sc_91_i1, d_18, qi_92_i1, d_18, sc_92_i1, d_83, ciy_92_i1, d_22, bpma_92_i1, 
d_30, bl_92_i1, d_14, qi_93_i1, d_87, tora_94_i1, d_12, ciy_94_i1, d_89, qi_94_i1, 
d_90, bpmf_95_i1, d_91, cix_95_i1, d_89, qi_95_i1, d_93,stlat_96_i1, bb_96_i1, d_94, bb_98_i1,
d_95, cbb_98_i1, d_96, bpms_99_i1, d_97, otrs_99_i1, d_98, bb_100_i1, d_95, cbb_100_i1, 
d_100, bb_101_i1, d_101, cbb_101_i1, d_102, qi_102_i1, d_89, cix_102_i1, d_104, bpmf_103_i1, 
d_105, ciy_103_i1, d_89, enlat_101_i1, qi_103_i1, d_10, bpma_103_i1, d_108, cix_104_i1, d_89, qi_104_i1,
match_104_i1, d_10, bpma_105_i1, d_111, ciy_107_i1, d_89, qi_107_i1, d_10, bpma_107_i1, d_111, 
cix_109_i1, d_89, qi_109_i1, d_10, bpma_110_i1, d_111, ciy_112_i1, d_89, qi_112_i1, d_10, 
bpma_112_i1, d_120, cix_114_i1, d_89, qi_114_i1, d_10, bpma_115_i1, d_123, tora_116_i1, d_124, 
ciy_116_i1, d_89, qi_116_i1, d_126, bpma_117_i1, d_127, otra_118_i1, d_128, cix_118_i1, d_89, 
qi_118_i1, d_130, bpma_119_i1, d_131, c_a2_1_1_l1, d_132, c_a2_1_2_l1, d_133, c_a2_1_3_l1, d_133, 
c_a2_1_4_l1, d_133, c_a2_1_5_l1, d_133, c_a2_1_6_l1, d_133, c_a2_1_7_l1, d_138, c_a2_1_8_l1, d_139, 
q_134_l1, d_140, cx_134_l1, cy_134_l1, d_141, bpmc_134_l1, d_142, c_a2_2_1_l1, d_133, c_a2_2_2_l1, 
d_133, c_a2_2_3_l1, d_133, c_a2_2_4_l1, d_133, c_a2_2_5_l1, d_132, c_a2_2_6_l1, d_133, c_a2_2_7_l1, 
d_133, c_a2_2_8_l1, d_150, q_146_l1, d_151, cx_146_l1, cy_146_l1, d_152, bpmr_146_l1, d_153, 
c_a2_3_1_l1, d_133, c_a2_3_2_l1, d_133, c_a2_3_3_l1, d_138, c_a2_3_4_l1, d_133, c_a2_3_5_l1, d_133, 
c_a2_3_6_l1, d_133, c_a2_3_7_l1, d_133, c_a2_3_8_l1, d_150, q_158_l1, d_162, cx_158_l1, cy_158_l1, 
d_163, bpmc_158_l1, d_142, c_a2_4_1_l1, d_133, c_a2_4_2_l1, d_133, c_a2_4_3_l1, d_133, c_a2_4_4_l1, 
d_133, c_a2_4_5_l1, d_133, c_a2_4_6_l1, d_133, c_a2_4_7_l1, d_133, c_a2_4_8_l1, d_139, q_170_l1, 
d_151, cx_170_l1, cy_170_l1, d_152, bpmr_170_l1, d_175, match_174_l1, d_176, tora_175_b1, d_177, 
bpma_175_b1, d_178, qi_176_b1, d_179, ciy_176_b1, d_180, cix_177_b1, d_181, bpma_179_b1, d_182, 
qd_179_b1, d_183, ccx_179_b1, d_184, otra_180_b1, d_185, qd_181_b1, d_186, ccy_181_b1, d_187, 
bpmf_181_b1, d_188_n, bb_182_b1, d_189, bb_191_b1, d_190, cbb_191_b1, d_96, bpms_192_b1, d_97,
otrs_192_b1, d_193, bb_193_b1, d_190, cbb_193_b1, d_195, bb_202_b1, d_190, cbb_202_b1, d_197, 
match_202_b1, d_198, bpmf_203_b1, d_199, tora_203_b1, d_200, ciy_204_b1, d_201, qi_204_b1, d_202, 
cix_205_b1, d_201, qi_205_b1, d_204, otra_206_b1, d_205, bpma_206_b1, d_206, qi_206_b1, d_3, 
match_207_b1, d_197, tdsb_208_b1, d_10, qi_209_b1, d_210, cix_209_b1, d_211, bpma_210_b1, d_182, 
qd_210_b1, d_213, ccy_210_b1, d_214, qi_211_b1, d_215, cix_213_b1, d_216, bpma_213_b1, d_178, 
qi_213_b1, d_3, ciy_214_b1, d_219, bpma_215_b1, d_178, qi_215_b1, d_201, cix_216_b1, d_222, 
ccy_217_b1, d_223, bpma_217_b1, d_178, qi_217_b1, d_225, match_218_b1, otrb_218_b1, d_226, bpma_219_b1, 
d_227, qd_219_b1, d_228, otrb_220_b1, d_229, bpma_221_b1, d_182, qd_221_b1, d_231, ccy_221_b1, 
d_232, otrb_222_b1, d_229, bpma_223_b1, d_182, qd_223_b1, d_235, cfx_223_b1, d_236, otrb_224_b1, 
d_237, qi_224_b1, d_238, bpma_226_b1, d_239, qi_226_b1, d_240, ciy_226_b1, d_241, cfx_226_b1, 
d_216, bpma_227_b1, d_239, qi_227_b1, d_244, ensub_229_b1)

# power supplies 

#  
qi_63_i1.ps_id = 'QI.13.I1'
qi_66_i1.ps_id = 'QI.14.I1'
qi_69_i1.ps_id = 'QI.15.I1'
qi_71_i1.ps_id = 'QI.16.I1'
qi_72_i1.ps_id = 'QI.17.I1'
qi_73_i1.ps_id = 'QI.18.I1'
qi_74_i1.ps_id = 'QI.19.I1'
qi_75_i1.ps_id = 'QI.20.I1'
qi_77_i1.ps_id = 'QI.19.I1'
qi_78_i1.ps_id = 'QI.18.I1'
qi_79_i1.ps_id = 'QI.19.I1'
qi_80_i1.ps_id = 'QI.20.I1'
qi_82_i1.ps_id = 'QI.21.I1'
qi_83_i1.ps_id = 'QI.22.I1'
qi_84_i1.ps_id = 'QI.21.I1'
qi_85_i1.ps_id = 'QI.24.I1'
qi_86_i1.ps_id = 'QI.23.I1'
qi_88_i1.ps_id = 'QI.22.I1'
qi_89_i1.ps_id = 'QI.23.I1'
qi_90_i1.ps_id = 'QI.24.I1'
qi_92_i1.ps_id = 'QI.23.I1'
qi_93_i1.ps_id = 'QI.25.I1'
qi_94_i1.ps_id = 'QI.26.I1'
qi_95_i1.ps_id = 'QI.27.I1'
qi_102_i1.ps_id = 'QI.28.I1'
qi_103_i1.ps_id = 'QI.29.I1'
qi_104_i1.ps_id = 'QI.30.I1'
qi_107_i1.ps_id = 'QI.31.I1'
qi_109_i1.ps_id = 'QI.32.I1'
qi_112_i1.ps_id = 'QI.31.I1'
qi_114_i1.ps_id = 'QI.33.I1'
qi_116_i1.ps_id = 'QI.34.I1'
qi_118_i1.ps_id = 'QI.35.I1'
q_134_l1.ps_id = 'Q.A2.1.L1'
q_146_l1.ps_id = 'Q.A2.2.L1'
q_158_l1.ps_id = 'Q.A2.3.L1'
q_170_l1.ps_id = 'Q.A2.4.L1'
qi_176_b1.ps_id = 'QI.1.B1'
qd_179_b1.ps_id = 'QD.3.B1'
qd_181_b1.ps_id = 'QD.4.B1'
qi_204_b1.ps_id = 'QI.5.B1'
qi_205_b1.ps_id = 'QI.6.B1'
qi_206_b1.ps_id = 'QI.7.B1'
qi_209_b1.ps_id = 'QI.8.B1'
qd_210_b1.ps_id = 'QD.9.B1'
qi_211_b1.ps_id = 'QI.10.B1'
qi_213_b1.ps_id = 'QI.11.B1'
qi_215_b1.ps_id = 'QI.12.B1'
qi_217_b1.ps_id = 'QI.13.B1'
qd_219_b1.ps_id = 'QD.14.B1'
qd_221_b1.ps_id = 'QD.15.B1'
qd_223_b1.ps_id = 'QD.16.B1'
qi_224_b1.ps_id = 'QI.17.B1'
qi_226_b1.ps_id = 'QI.18.B1'
qi_227_b1.ps_id = 'QI.19.B1'

#  
sc_74i_i1.ps_id = 'SC.1.I1'
sc_74ii_i1.ps_id = 'SC.2.I1'
sc_76_i1.ps_id = 'SC.2.I1'
sc_77_i1.ps_id = 'SC.1.I1'
sc_79i_i1.ps_id = 'SC.1.I1'
sc_79ii_i1.ps_id = 'SC.2.I1'
sc_81_i1.ps_id = 'SC.2.I1'
sc_82_i1.ps_id = 'SC.1.I1'
sc_84i_i1.ps_id = 'SC.1.I1'
sc_84ii_i1.ps_id = 'SC.2.I1'
sc_86_i1.ps_id = 'SC.2.I1'
sc_87_i1.ps_id = 'SC.1.I1'
sc_89i_i1.ps_id = 'SC.1.I1'
sc_89ii_i1.ps_id = 'SC.2.I1'
sc_91_i1.ps_id = 'SC.2.I1'
sc_92_i1.ps_id = 'SC.1.I1'

#  

#  
c_a2_1_1_l1.ps_id = 'C.A2.L1'
c_a2_1_2_l1.ps_id = 'C.A2.L1'
c_a2_1_3_l1.ps_id = 'C.A2.L1'
c_a2_1_4_l1.ps_id = 'C.A2.L1'
c_a2_1_5_l1.ps_id = 'C.A2.L1'
c_a2_1_6_l1.ps_id = 'C.A2.L1'
c_a2_1_7_l1.ps_id = 'C.A2.L1'
c_a2_1_8_l1.ps_id = 'C.A2.L1'
c_a2_2_1_l1.ps_id = 'C.A2.L1'
c_a2_2_2_l1.ps_id = 'C.A2.L1'
c_a2_2_3_l1.ps_id = 'C.A2.L1'
c_a2_2_4_l1.ps_id = 'C.A2.L1'
c_a2_2_5_l1.ps_id = 'C.A2.L1'
c_a2_2_6_l1.ps_id = 'C.A2.L1'
c_a2_2_7_l1.ps_id = 'C.A2.L1'
c_a2_2_8_l1.ps_id = 'C.A2.L1'
c_a2_3_1_l1.ps_id = 'C.A2.L1'
c_a2_3_2_l1.ps_id = 'C.A2.L1'
c_a2_3_3_l1.ps_id = 'C.A2.L1'
c_a2_3_4_l1.ps_id = 'C.A2.L1'
c_a2_3_5_l1.ps_id = 'C.A2.L1'
c_a2_3_6_l1.ps_id = 'C.A2.L1'
c_a2_3_7_l1.ps_id = 'C.A2.L1'
c_a2_3_8_l1.ps_id = 'C.A2.L1'
c_a2_4_1_l1.ps_id = 'C.A2.L1'
c_a2_4_2_l1.ps_id = 'C.A2.L1'
c_a2_4_3_l1.ps_id = 'C.A2.L1'
c_a2_4_4_l1.ps_id = 'C.A2.L1'
c_a2_4_5_l1.ps_id = 'C.A2.L1'
c_a2_4_6_l1.ps_id = 'C.A2.L1'
c_a2_4_7_l1.ps_id = 'C.A2.L1'
c_a2_4_8_l1.ps_id = 'C.A2.L1'
tdsb_208_b1.ps_id = 'TDSB.B1'

#  
bl_73_i1.ps_id = 'BL.6.I1'
bl_75_i1.ps_id = 'BL.7.I1'
bl_76_i1.ps_id = 'BL.7.I1'
bl_77_i1.ps_id = 'BL.6.I1'
bl_78_i1.ps_id = 'BL.6.I1'
bl_80_i1.ps_id = 'BL.7.I1'
bl_81_i1.ps_id = 'BL.7.I1'
bl_82_i1.ps_id = 'BL.6.I1'
bl_83_i1.ps_id = 'BL.8.I1'
bl_85_i1.ps_id = 'BL.7.I1'
bl_86_i1.ps_id = 'BL.7.I1'
bl_87_i1.ps_id = 'BL.8.I1'
bl_88_i1.ps_id = 'BL.8.I1'
bl_90_i1.ps_id = 'BL.7.I1'
bl_91_i1.ps_id = 'BL.7.I1'
bl_92_i1.ps_id = 'BL.8.I1'
bb_96_i1.ps_id = 'BB.1.I1'
bb_98_i1.ps_id = 'BB.1.I1'
bb_100_i1.ps_id = 'BB.1.I1'
bb_101_i1.ps_id = 'BB.1.I1'
bb_182_b1.ps_id = 'BB.1.B1'
bb_191_b1.ps_id = 'BB.1.B1'
bb_193_b1.ps_id = 'BB.1.B1'
bb_202_b1.ps_id = 'BB.1.B1'
