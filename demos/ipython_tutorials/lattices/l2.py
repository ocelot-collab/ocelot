# Converted from component_list_2026.01.21.xls

from ocelot import (
    Twiss,
    Cavity,
    Drift,
    Hcor,
    Marker,
    Monitor,
    Quadrupole,
    RBend,
    SBend,
    TDCavity,
    Vcor,
)

twiss0 = Twiss()
twiss0.E = 0.7000000007706801
twiss0.alpha_x = -0.6943989709726064
twiss0.alpha_y = -1.2901254625463643
twiss0.beta_x = 7.488039009254501
twiss0.beta_y = 8.698683623541857
twiss0.s = 206.100754


# Drifts:
d_0 = Drift(l=0.7527499999999918, eid="D_0")
d_1 = Drift(l=1.9499960000000272, eid="D_1")
d_2 = Drift(l=0.03164999999998486, eid="D_2")
d_3 = Drift(l=0.11164999999998218, eid="D_3")
d_4 = Drift(l=0.15000000000002842, eid="D_4")
d_5 = Drift(l=0.1516499999999894, eid="D_5")
d_6 = Drift(l=0.1216500000000015, eid="D_6")
d_7 = Drift(l=0.15, eid="D_7")
d_8 = Drift(l=0.1416499999999985, eid="D_8")
d_9 = Drift(l=0.1316499999999924, eid="D_9")
d_10 = Drift(l=0.07423000000002844, eid="D_10")
vcdst_234_b1 = Drift(l=1.215, eid="VCDST.234.B1")
d_11 = Drift(l=0.12576999999996352, eid="D_11")
cfb_236_l2 = Drift(l=3.35, eid="CFB.236.L2")
d_12 = Drift(l=0.22159999999999513, eid="D_12")
d_13 = Drift(l=0.3459000000000012, eid="D_13")
d_14 = Drift(l=0.3459000000000012, eid="D_14")
d_15 = Drift(l=0.3459000000000012, eid="D_15")
d_16 = Drift(l=0.3459000000000012, eid="D_16")
d_17 = Drift(l=0.3459000000000012, eid="D_17")
d_18 = Drift(l=0.3459000000000012, eid="D_18")
d_19 = Drift(l=0.3459000000000012, eid="D_19")
d_20 = Drift(l=0.24750000000000316, eid="D_20")
d_21 = Drift(l=0.04319999999999835, eid="D_21")
d_22 = Drift(l=0.08499999999997954, eid="D_22")
d_23 = Drift(l=0.4578999999999951, eid="D_23")
d_24 = Drift(l=0.22160000000002356, eid="D_24")
d_25 = Drift(l=0.3459000000000012, eid="D_25")
d_26 = Drift(l=0.3459000000000012, eid="D_26")
d_27 = Drift(l=0.3459000000000012, eid="D_27")
d_28 = Drift(l=0.3459000000000012, eid="D_28")
d_29 = Drift(l=0.3459000000000012, eid="D_29")
d_30 = Drift(l=0.3459000000000012, eid="D_30")
d_31 = Drift(l=0.3459000000000012, eid="D_31")
d_32 = Drift(l=0.24749999999997474, eid="D_32")
d_33 = Drift(l=0.04319999999999835, eid="D_33")
d_34 = Drift(l=0.08500000000000796, eid="D_34")
d_35 = Drift(l=0.4578999999999951, eid="D_35")
d_36 = Drift(l=0.22160000000002356, eid="D_36")
d_37 = Drift(l=0.3459000000000012, eid="D_37")
d_38 = Drift(l=0.3459000000000012, eid="D_38")
d_39 = Drift(l=0.3459000000000012, eid="D_39")
d_40 = Drift(l=0.3459000000000012, eid="D_40")
d_41 = Drift(l=0.3459000000000012, eid="D_41")
d_42 = Drift(l=0.3459000000000012, eid="D_42")
d_43 = Drift(l=0.3459000000000012, eid="D_43")
d_44 = Drift(l=0.24749999999997474, eid="D_44")
d_45 = Drift(l=0.04319999999999835, eid="D_45")
d_46 = Drift(l=0.08500000000000796, eid="D_46")
d_47 = Drift(l=0.4578999999999951, eid="D_47")
d_48 = Drift(l=0.22160000000002356, eid="D_48")
d_49 = Drift(l=0.3458999999999728, eid="D_49")
d_50 = Drift(l=0.3459000000000012, eid="D_50")
d_51 = Drift(l=0.34590000000002963, eid="D_51")
d_52 = Drift(l=0.3459000000000012, eid="D_52")
d_53 = Drift(l=0.3459000000000012, eid="D_53")
d_54 = Drift(l=0.3459000000000012, eid="D_54")
d_55 = Drift(l=0.3459000000000012, eid="D_55")
d_56 = Drift(l=0.2474999999999179, eid="D_56")
d_57 = Drift(l=0.043200000000055194, eid="D_57")
d_58 = Drift(l=0.08499999999997954, eid="D_58")
d_59 = Drift(l=0.4578999999999951, eid="D_59")
d_60 = Drift(l=0.22160000000002356, eid="D_60")
d_61 = Drift(l=0.3459000000000012, eid="D_61")
d_62 = Drift(l=0.3459000000000012, eid="D_62")
d_63 = Drift(l=0.3459000000000012, eid="D_63")
d_64 = Drift(l=0.3459000000000012, eid="D_64")
d_65 = Drift(l=0.3459000000000012, eid="D_65")
d_66 = Drift(l=0.3459000000000012, eid="D_66")
d_67 = Drift(l=0.3459000000000012, eid="D_67")
d_68 = Drift(l=0.24749999999997474, eid="D_68")
d_69 = Drift(l=0.04319999999999835, eid="D_69")
d_70 = Drift(l=0.08500000000003638, eid="D_70")
d_71 = Drift(l=0.4578999999999951, eid="D_71")
d_72 = Drift(l=0.2215999999999667, eid="D_72")
d_73 = Drift(l=0.3459000000000012, eid="D_73")
d_74 = Drift(l=0.3459000000000012, eid="D_74")
d_75 = Drift(l=0.3459000000000012, eid="D_75")
d_76 = Drift(l=0.3459000000000012, eid="D_76")
d_77 = Drift(l=0.3459000000000012, eid="D_77")
d_78 = Drift(l=0.3459000000000012, eid="D_78")
d_79 = Drift(l=0.3459000000000012, eid="D_79")
d_80 = Drift(l=0.24749999999997474, eid="D_80")
d_81 = Drift(l=0.043200000000055194, eid="D_81")
d_82 = Drift(l=0.08499999999997954, eid="D_82")
d_83 = Drift(l=0.4578999999999951, eid="D_83")
d_84 = Drift(l=0.22160000000002356, eid="D_84")
d_85 = Drift(l=0.3459000000000012, eid="D_85")
d_86 = Drift(l=0.3459000000000012, eid="D_86")
d_87 = Drift(l=0.3459000000000012, eid="D_87")
d_88 = Drift(l=0.3459000000000012, eid="D_88")
d_89 = Drift(l=0.3459000000000012, eid="D_89")
d_90 = Drift(l=0.3459000000000012, eid="D_90")
d_91 = Drift(l=0.3459000000000012, eid="D_91")
d_92 = Drift(l=0.2474999999999179, eid="D_92")
d_93 = Drift(l=0.043200000000055194, eid="D_93")
d_94 = Drift(l=0.08499999999997954, eid="D_94")
d_95 = Drift(l=0.4578999999999951, eid="D_95")
d_96 = Drift(l=0.22160000000002356, eid="D_96")
d_97 = Drift(l=0.3459000000000012, eid="D_97")
d_98 = Drift(l=0.3459000000000012, eid="D_98")
d_99 = Drift(l=0.3459000000000012, eid="D_99")
d_100 = Drift(l=0.3459000000000012, eid="D_100")
d_101 = Drift(l=0.3459000000000012, eid="D_101")
d_102 = Drift(l=0.3459000000000012, eid="D_102")
d_103 = Drift(l=0.3459000000000012, eid="D_103")
d_104 = Drift(l=0.24749999999997474, eid="D_104")
d_105 = Drift(l=0.04319999999999835, eid="D_105")
d_106 = Drift(l=0.08500000000003638, eid="D_106")
d_107 = Drift(l=0.4578999999999951, eid="D_107")
d_108 = Drift(l=0.22160000000002356, eid="D_108")
d_109 = Drift(l=0.34589999999994436, eid="D_109")
d_110 = Drift(l=0.3459000000000012, eid="D_110")
d_111 = Drift(l=0.3459000000000012, eid="D_111")
d_112 = Drift(l=0.3459000000000012, eid="D_112")
d_113 = Drift(l=0.3459000000000012, eid="D_113")
d_114 = Drift(l=0.3459000000000012, eid="D_114")
d_115 = Drift(l=0.3459000000000012, eid="D_115")
d_116 = Drift(l=0.24749999999997474, eid="D_116")
d_117 = Drift(l=0.043200000000055194, eid="D_117")
d_118 = Drift(l=0.08499999999997954, eid="D_118")
d_119 = Drift(l=0.4578999999999951, eid="D_119")
d_120 = Drift(l=0.22160000000002356, eid="D_120")
d_121 = Drift(l=0.3459000000000012, eid="D_121")
d_122 = Drift(l=0.3459000000000012, eid="D_122")
d_123 = Drift(l=0.3459000000000012, eid="D_123")
d_124 = Drift(l=0.3459000000000012, eid="D_124")
d_125 = Drift(l=0.3459000000000012, eid="D_125")
d_126 = Drift(l=0.3459000000000012, eid="D_126")
d_127 = Drift(l=0.3459000000000012, eid="D_127")
d_128 = Drift(l=0.24749999999997474, eid="D_128")
d_129 = Drift(l=0.04319999999999835, eid="D_129")
d_130 = Drift(l=0.08499999999997954, eid="D_130")
d_131 = Drift(l=0.4578999999999951, eid="D_131")
d_132 = Drift(l=0.22160000000002356, eid="D_132")
d_133 = Drift(l=0.3459000000000012, eid="D_133")
d_134 = Drift(l=0.3459000000000012, eid="D_134")
d_135 = Drift(l=0.3459000000000012, eid="D_135")
d_136 = Drift(l=0.3459000000000012, eid="D_136")
d_137 = Drift(l=0.3459000000000012, eid="D_137")
d_138 = Drift(l=0.3459000000000012, eid="D_138")
d_139 = Drift(l=0.3459000000000012, eid="D_139")
d_140 = Drift(l=0.24749999999997474, eid="D_140")
d_141 = Drift(l=0.04319999999999835, eid="D_141")
d_142 = Drift(l=0.08500000000003638, eid="D_142")
d_143 = Drift(l=0.4578999999999951, eid="D_143")
d_144 = Drift(l=0.22160000000002356, eid="D_144")
d_145 = Drift(l=0.3459000000000012, eid="D_145")
d_146 = Drift(l=0.34589999999994436, eid="D_146")
d_147 = Drift(l=0.3459000000000012, eid="D_147")
d_148 = Drift(l=0.3459000000000012, eid="D_148")
d_149 = Drift(l=0.3459000000000012, eid="D_149")
d_150 = Drift(l=0.3459000000000012, eid="D_150")
d_151 = Drift(l=0.3459000000000012, eid="D_151")
d_152 = Drift(l=0.24749999999997474, eid="D_152")
d_153 = Drift(l=0.043200000000055194, eid="D_153")
d_154 = Drift(l=0.08499999999997954, eid="D_154")
d_155 = Drift(l=0.4578999999999951, eid="D_155")
ctb_383_l2 = Drift(l=3.25, eid="CTB.383.L2")
d_156 = Drift(l=0.12574999999998226, eid="D_156")
vcdst_386_b2 = Drift(l=1.215, eid="VCDST.386.B2")
d_157 = Drift(l=0.07425000000003812, eid="D_157")
d_158 = Drift(l=0.13087999999999056, eid="D_158")
d_159 = Drift(l=0.24799999999999045, eid="D_159")
d_160 = Drift(l=0.15276999999997543, eid="D_160")
d_161 = Drift(l=0.13165000000002083, eid="D_161")
d_162 = Drift(l=0.3000800000000027, eid="D_162")
d_163 = Drift(l=0.3275499999999738, eid="D_163")
d_164 = Drift(l=0.13567000000004628, eid="D_164")
d_165 = Drift(l=1.0499999999999772, eid="D_165")
d_166 = Drift(l=0.10000000000002274, eid="D_166")
d_167 = Drift(l=0.288880000000006, eid="D_167")
d_168 = Drift(l=0.339999999999975, eid="D_168")
d_169 = Drift(l=0.15276999999997543, eid="D_169")
d_170 = Drift(l=0.10065000000001492, eid="D_170")
d_171 = Drift(l=0.2918800000000147, eid="D_171")
d_172 = Drift(l=0.2707699999999704, eid="D_172")
d_173 = Drift(l=0.08163000000000312, eid="D_173")
d_174 = Drift(l=0.34390000000004195, eid="D_174")
d_175 = Drift(l=0.19399999999995998, eid="D_175")
d_176 = Drift(l=0.1564999999999941, eid="D_176")
d_177 = Drift(l=0.09600000000000364, eid="D_177")
d_178 = Drift(l=0.15962000000001808, eid="D_178")
d_179 = Drift(l=0.30000000000001137, eid="D_179")
d_180 = Drift(l=0.10000000000002274, eid="D_180")
d_181 = Drift(l=7.299999998622297e-05, eid="D_181")
d_182 = Drift(l=7.199999998874773e-05, eid="D_182")
d_183 = Drift(l=8.50737399999997, eid="D_183")
d_184 = Drift(l=7.200000004559115e-05, eid="D_184")
d_185 = Drift(l=7.199999998874773e-05, eid="D_185")
d_186 = Drift(l=0.39749999999997954, eid="D_186")
d_187 = Drift(l=0.08500000000003638, eid="D_187")
d_188 = Drift(l=0.3824999999999932, eid="D_188")
d_189 = Drift(l=0.3100000000000023, eid="D_189")
d_190 = Drift(l=0.32499999999998863, eid="D_190")
d_191 = Drift(l=7.299999998622297e-05, eid="D_191")
d_192 = Drift(l=7.199999998874773e-05, eid="D_192")
d_193 = Drift(l=1.0531890000000317, eid="D_193")
d_194 = Drift(l=7.454184999999995, eid="D_194")
d_195 = Drift(l=7.199999998874773e-05, eid="D_195")
d_196 = Drift(l=7.199999998874773e-05, eid="D_196")
d_197 = Drift(l=0.10000000000002274, eid="D_197")
d_198 = Drift(l=0.17000000000001592, eid="D_198")
d_199 = Drift(l=0.15197999999998046, eid="D_199")
d_200 = Drift(l=0.1564999999999941, eid="D_200")
d_201 = Drift(l=0.09600000000000364, eid="D_201")
d_202 = Drift(l=0.30680000000000973, eid="D_202")
d_203 = Drift(l=0.15734999999995125, eid="D_203")
d_204 = Drift(l=0.14464999999999764, eid="D_204")
d_205 = Drift(l=0.15002000000000634, eid="D_205")
d_206 = Drift(l=0.34888000000003105, eid="D_206")
d_207 = Drift(l=0.5827699999999822, eid="D_207")
d_208 = Drift(l=0.7816499999999981, eid="D_208")
d_209 = Drift(l=0.15, eid="D_209")
d_210 = Drift(l=0.1788800000000151, eid="D_210")
d_211 = Drift(l=0.15276999999997543, eid="D_211")
d_212 = Drift(l=0.13165000000002083, eid="D_212")
d_213 = Drift(l=4.0, eid="D_213")
d_214 = Drift(l=2.399979999999971, eid="D_214")
d_215 = Drift(l=0.17165000000002237, eid="D_215")
d_216 = Drift(l=0.1916500000000231, eid="D_216")
d_217 = Drift(l=0.3898999999999774, eid="D_217")
d_218 = Drift(l=0.2470000000000141, eid="D_218")
d_219 = Drift(l=0.4347699999999577, eid="D_219")
d_220 = Drift(l=0.09165000000005721, eid="D_220")
d_221 = Drift(l=0.19999999999998863, eid="D_221")
d_222 = Drift(l=0.09427999999996928, eid="D_222")
d_223 = Drift(l=0.4000000000000341, eid="D_223")
d_224 = Drift(l=0.17834999999996626, eid="D_224")
d_225 = Drift(l=0.09065000000002402, eid="D_225")
d_226 = Drift(l=0.322300000000007, eid="D_226")
d_227 = Drift(l=0.17772000000002208, eid="D_227")
d_228 = Drift(l=0.19999999999998863, eid="D_228")
d_229 = Drift(l=1.5316499999999564, eid="D_229")
d_230 = Drift(l=0.13165000000002083, eid="D_230")
d_231 = Drift(l=2.2316499999999677, eid="D_231")
d_232 = Drift(l=0.13165000000002083, eid="D_232")
d_233 = Drift(l=0.20000000000004547, eid="D_233")
d_234 = Drift(l=1.4498799999999505, eid="D_234")
d_235 = Drift(l=1.629000000000019, eid="D_235")
d_236 = Drift(l=0.15276999999997543, eid="D_236")
d_237 = Drift(l=0.13165000000002083, eid="D_237")
d_238 = Drift(l=0.9, eid="D_238")
d_239 = Drift(l=0.19999999999998863, eid="D_239")
d_240 = Drift(l=2.0788800000000265, eid="D_240")
d_241 = Drift(l=0.15276999999997543, eid="D_241")
d_242 = Drift(l=0.22553000000003368, eid="D_242")
d_243 = Drift(l=0.0699999999999591, eid="D_243")
d_244 = Drift(l=0.7860000000000241, eid="D_244")
d_245 = Drift(l=0.350120000000004, eid="D_245")
d_246 = Drift(l=0.19999999999998863, eid="D_246")
d_247 = Drift(l=0.35000000000002274, eid="D_247")
d_248 = Drift(l=0.6288800000000038, eid="D_248")
d_249 = Drift(l=0.15276999999997543, eid="D_249")
d_250 = Drift(l=0.13165000000002083, eid="D_250")
d_251 = Drift(l=1.549879999999996, eid="D_251")
d_252 = Drift(l=1.350120000000004, eid="D_252")
d_253 = Drift(l=0.19999999999998863, eid="D_253")
d_254 = Drift(l=0.1048799999999801, eid="D_254")
d_255 = Drift(l=0.15285000000000082, eid="D_255")
d_256 = Drift(l=0.19944999999999782, eid="D_256")
d_257 = Drift(l=0.07000000000001594, eid="D_257")
d_258 = Drift(l=0.7860000000000241, eid="D_258")
d_259 = Drift(l=0.3500999999999408, eid="D_259")
d_260 = Drift(l=1.1789000000000214, eid="D_260")
d_261 = Drift(l=0.15276999999997543, eid="D_261")
d_262 = Drift(l=0.13165000000002083, eid="D_262")
d_263 = Drift(l=0.10000000000004547, eid="D_263")
d_264 = Drift(l=0.19999999999998863, eid="D_264")
d_265 = Drift(l=1.249879999999962, eid="D_265")
d_266 = Drift(l=1.629000000000019, eid="D_266")
d_267 = Drift(l=0.15276999999997543, eid="D_267")
d_268 = Drift(l=0.13165000000002083, eid="D_268")
d_269 = Drift(l=1.0999999999999885, eid="D_269")
d_270 = Drift(l=0.20000000000004547, eid="D_270")
d_271 = Drift(l=0.2498799999999619, eid="D_271")
d_272 = Drift(l=1.0790000000000077, eid="D_272")
d_273 = Drift(l=0.15274999999996908, eid="D_273")
d_274 = Drift(l=0.3705500000000218, eid="D_274")
d_275 = Drift(l=0.461119999999994, eid="D_275")
d_276 = Drift(l=0.33164999999999056, eid="D_276")
d_277 = Drift(l=0.13165000000002083, eid="D_277")
d_278 = Drift(l=0.34999999999998865, eid="D_278")
d_279 = Drift(l=0.1788800000000151, eid="D_279")
d_280 = Drift(l=0.15276999999997543, eid="D_280")
d_281 = Drift(l=0.13165000000002083, eid="D_281")
d_282 = Drift(l=0.4000000000000341, eid="D_282")

# Quadrupoles:
qd_231_b1 = Quadrupole(l=0.2367, k1=1.6557307729995776, eid="QD.231.B1")
qd_232_b1 = Quadrupole(l=0.2367, k1=-1.2501718110012676, eid="QD.232.B1")
qd_233_b1 = Quadrupole(l=0.2367, k1=-0.4397248283016477, eid="QD.233.B1")
q_249_l2 = Quadrupole(l=0.2136, k1=0.424912906900749, eid="Q.249.L2")
q_261_l2 = Quadrupole(l=0.2136, k1=-0.39384852710205986, eid="Q.261.L2")
q_273_l2 = Quadrupole(l=0.2136, k1=0.387775731900749, eid="Q.273.L2")
q_285_l2 = Quadrupole(l=0.2136, k1=-0.4247537369990636, eid="Q.285.L2")
q_297_l2 = Quadrupole(l=0.2136, k1=0.4881898113998127, eid="Q.297.L2")
q_309_l2 = Quadrupole(l=0.2136, k1=-0.6058097569990636, eid="Q.309.L2")
q_321_l2 = Quadrupole(l=0.2136, k1=0.4881898113998127, eid="Q.321.L2")
q_333_l2 = Quadrupole(l=0.2136, k1=-0.6058097569990636, eid="Q.333.L2")
q_345_l2 = Quadrupole(l=0.2136, k1=0.6525829807022472, eid="Q.345.L2")
q_357_l2 = Quadrupole(l=0.2136, k1=-0.38908541679775277, eid="Q.357.L2")
q_369_l2 = Quadrupole(l=0.2136, k1=0.43128398199906365, eid="Q.369.L2")
q_381_l2 = Quadrupole(l=0.2136, k1=-0.38391488749999997, eid="Q.381.L2")
qd_387_b2 = Quadrupole(l=0.2367, k1=0.335173384799324, eid="QD.387.B2")
qd_388_b2 = Quadrupole(l=0.2367, k1=0.35599643220109844, eid="QD.388.B2")
qd_391_b2 = Quadrupole(l=0.2367, k1=-0.7255245525982257, eid="QD.391.B2")
qd_392_b2 = Quadrupole(l=0.2367, k1=0.19699609059991552, eid="QD.392.B2")
qd_415_b2 = Quadrupole(l=0.2367, k1=0.1806864174017744, eid="QD.415.B2")
qd_417_b2 = Quadrupole(l=0.2367, k1=-0.7502347655006337, eid="QD.417.B2")
qd_418_b2 = Quadrupole(l=0.2367, k1=0.6491929171989861, eid="QD.418.B2")
qd_425_b2 = Quadrupole(l=0.2367, k1=-1.3008034, eid="QD.425.B2")
qd_427_b2 = Quadrupole(l=0.2367, k1=0.9414835846007604, eid="QD.427.B2")
qd_431_b2 = Quadrupole(l=0.2367, k1=0.43518302749894383, eid="QD.431.B2")
qd_434_b2 = Quadrupole(l=0.2367, k1=-0.5278581910012674, eid="QD.434.B2")
qd_437_b2 = Quadrupole(l=0.2367, k1=0.4055492834980989, eid="QD.437.B2")
qd_440_b2 = Quadrupole(l=0.2367, k1=-0.6685246719983101, eid="QD.440.B2")
qd_444_b2 = Quadrupole(l=0.2367, k1=-0.4582186614997888, eid="QD.444.B2")
qd_448_b2 = Quadrupole(l=0.2367, k1=0.8960955489987326, eid="QD.448.B2")
qd_452_b2 = Quadrupole(l=0.2367, k1=-1.263284384000845, eid="QD.452.B2")
qd_456_b2 = Quadrupole(l=0.2367, k1=0.8960955489987326, eid="QD.456.B2")
qd_459_b2 = Quadrupole(l=0.2367, k1=-1.263284384000845, eid="QD.459.B2")
qd_463_b2 = Quadrupole(l=0.2367, k1=-0.569607097600338, eid="QD.463.B2")
qd_464_b2 = Quadrupole(l=0.2367, k1=1.2982678500000002, eid="QD.464.B2")
qd_465_b2 = Quadrupole(l=0.2367, k1=-0.24686100550063372, eid="QD.465.B2")

# SBends:
bb_393_b2 = SBend(
    l=0.5, angle=0.04163879892, e2=0.041638799, tilt=1.570796327, eid="BB.393.B2"
)
bb_402_b2 = SBend(
    l=0.5, angle=-0.04163879892, e1=-0.041638799, tilt=1.570796327, eid="BB.402.B2"
)
bb_404_b2 = SBend(
    l=0.5, angle=-0.04163879892, e2=-0.041638799, tilt=1.570796327, eid="BB.404.B2"
)
bb_413_b2 = SBend(
    l=0.5, angle=0.04163879892, e1=0.041638799, tilt=1.570796327, eid="BB.413.B2"
)

# RBends:
kdy_445_b2 = RBend(l=0.35, e1=0.0, e2=0.0, eid="KDY.445.B2")
kdy_446_b2 = RBend(l=0.35, e1=0.0, e2=0.0, eid="KDY.446.B2")
kdy_452_b2 = RBend(l=0.35, e1=0.0, e2=0.0, eid="KDY.452.B2")
kdy_453_b2 = RBend(l=0.35, e1=0.0, e2=0.0, eid="KDY.453.B2")

# Hcors:
ccx_232_b1 = Hcor(l=0.1, eid="CCX.232.B1")
cx_249_l2 = Hcor(eid="CX.249.L2")
cx_273_l2 = Hcor(eid="CX.273.L2")
cx_297_l2 = Hcor(eid="CX.297.L2")
cx_321_l2 = Hcor(eid="CX.321.L2")
cx_345_l2 = Hcor(eid="CX.345.L2")
cx_369_l2 = Hcor(eid="CX.369.L2")
ccx_388_b2 = Hcor(l=0.1, eid="CCX.388.B2")
ccx_392_b2 = Hcor(l=0.1, eid="CCX.392.B2")
ccx_415_b2 = Hcor(l=0.1, eid="CCX.415.B2")
ccx_418_b2 = Hcor(l=0.1, eid="CCX.418.B2")
ccx_425_b2 = Hcor(l=0.1, eid="CCX.425.B2")
ccx_431_b2 = Hcor(l=0.1, eid="CCX.431.B2")
ccx_441_b2 = Hcor(l=0.1, eid="CCX.441.B2")
ccx_447_b2 = Hcor(l=0.1, eid="CCX.447.B2")
ccx_454_b2 = Hcor(l=0.1, eid="CCX.454.B2")
ccx_456_b2 = Hcor(l=0.1, eid="CCX.456.B2")
ccx_464_b2 = Hcor(l=0.1, eid="CCX.464.B2")
ccx_465_b2 = Hcor(l=0.1, eid="CCX.465.B2")

# Vcors:
cbb_229_b1d = Vcor(eid="CBB.229.B1D")
ccy_232_b1 = Vcor(l=0.1, eid="CCY.232.B1")
cy_261_l2 = Vcor(eid="CY.261.L2")
cy_285_l2 = Vcor(eid="CY.285.L2")
cy_309_l2 = Vcor(eid="CY.309.L2")
cy_333_l2 = Vcor(eid="CY.333.L2")
cy_357_l2 = Vcor(eid="CY.357.L2")
cy_381_l2 = Vcor(eid="CY.381.L2")
ccy_387_b2 = Vcor(l=0.1, eid="CCY.387.B2")
ccy_391_b2 = Vcor(l=0.1, eid="CCY.391.B2")
cbb_403_b2 = Vcor(eid="CBB.403.B2")
cbb_405_b2 = Vcor(eid="CBB.405.B2")
cbb_414_b2 = Vcor(eid="CBB.414.B2")
ccy_416_b2 = Vcor(l=0.1, eid="CCY.416.B2")
ccy_418_b2 = Vcor(l=0.1, eid="CCY.418.B2")
ccy_426_b2 = Vcor(l=0.1, eid="CCY.426.B2")
ccy_434_b2 = Vcor(l=0.1, eid="CCY.434.B2")
ccy_448_b2 = Vcor(l=0.1, eid="CCY.448.B2")
ccy_460_b2 = Vcor(l=0.1, eid="CCY.460.B2")
ccy_464_b2 = Vcor(l=0.1, eid="CCY.464.B2")

# Cavitys:
c_a3_1_1_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.1.1.L2")
c_a3_1_2_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.1.2.L2")
c_a3_1_3_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.1.3.L2")
c_a3_1_4_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.1.4.L2")
c_a3_1_5_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.1.5.L2")
c_a3_1_6_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.1.6.L2")
c_a3_1_7_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.1.7.L2")
c_a3_1_8_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.1.8.L2")
c_a3_2_1_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.2.1.L2")
c_a3_2_2_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.2.2.L2")
c_a3_2_3_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.2.3.L2")
c_a3_2_4_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.2.4.L2")
c_a3_2_5_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.2.5.L2")
c_a3_2_6_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.2.6.L2")
c_a3_2_7_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.2.7.L2")
c_a3_2_8_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.2.8.L2")
c_a3_3_1_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.3.1.L2")
c_a3_3_2_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.3.2.L2")
c_a3_3_3_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.3.3.L2")
c_a3_3_4_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.3.4.L2")
c_a3_3_5_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.3.5.L2")
c_a3_3_6_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.3.6.L2")
c_a3_3_7_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.3.7.L2")
c_a3_3_8_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.3.8.L2")
c_a3_4_1_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.4.1.L2")
c_a3_4_2_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.4.2.L2")
c_a3_4_3_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.4.3.L2")
c_a3_4_4_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.4.4.L2")
c_a3_4_5_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.4.5.L2")
c_a3_4_6_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.4.6.L2")
c_a3_4_7_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.4.7.L2")
c_a3_4_8_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A3.4.8.L2")
c_a4_1_1_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.1.1.L2")
c_a4_1_2_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.1.2.L2")
c_a4_1_3_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.1.3.L2")
c_a4_1_4_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.1.4.L2")
c_a4_1_5_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.1.5.L2")
c_a4_1_6_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.1.6.L2")
c_a4_1_7_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.1.7.L2")
c_a4_1_8_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.1.8.L2")
c_a4_2_1_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.2.1.L2")
c_a4_2_2_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.2.2.L2")
c_a4_2_3_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.2.3.L2")
c_a4_2_4_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.2.4.L2")
c_a4_2_5_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.2.5.L2")
c_a4_2_6_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.2.6.L2")
c_a4_2_7_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.2.7.L2")
c_a4_2_8_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.2.8.L2")
c_a4_3_1_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.3.1.L2")
c_a4_3_2_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.3.2.L2")
c_a4_3_3_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.3.3.L2")
c_a4_3_4_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.3.4.L2")
c_a4_3_5_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.3.5.L2")
c_a4_3_6_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.3.6.L2")
c_a4_3_7_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.3.7.L2")
c_a4_3_8_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.3.8.L2")
c_a4_4_1_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.4.1.L2")
c_a4_4_2_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.4.2.L2")
c_a4_4_3_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.4.3.L2")
c_a4_4_4_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.4.4.L2")
c_a4_4_5_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.4.5.L2")
c_a4_4_6_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.4.6.L2")
c_a4_4_7_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.4.7.L2")
c_a4_4_8_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A4.4.8.L2")
c_a5_1_1_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.1.1.L2")
c_a5_1_2_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.1.2.L2")
c_a5_1_3_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.1.3.L2")
c_a5_1_4_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.1.4.L2")
c_a5_1_5_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.1.5.L2")
c_a5_1_6_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.1.6.L2")
c_a5_1_7_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.1.7.L2")
c_a5_1_8_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.1.8.L2")
c_a5_2_1_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.2.1.L2")
c_a5_2_2_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.2.2.L2")
c_a5_2_3_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.2.3.L2")
c_a5_2_4_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.2.4.L2")
c_a5_2_5_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.2.5.L2")
c_a5_2_6_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.2.6.L2")
c_a5_2_7_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.2.7.L2")
c_a5_2_8_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.2.8.L2")
c_a5_3_1_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.3.1.L2")
c_a5_3_2_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.3.2.L2")
c_a5_3_3_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.3.3.L2")
c_a5_3_4_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.3.4.L2")
c_a5_3_5_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.3.5.L2")
c_a5_3_6_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.3.6.L2")
c_a5_3_7_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.3.7.L2")
c_a5_3_8_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.3.8.L2")
c_a5_4_1_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.4.1.L2")
c_a5_4_2_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.4.2.L2")
c_a5_4_3_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.4.3.L2")
c_a5_4_4_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.4.4.L2")
c_a5_4_5_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.4.5.L2")
c_a5_4_6_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.4.6.L2")
c_a5_4_7_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.4.7.L2")
c_a5_4_8_l2 = Cavity(l=1.0377, v=0.01770833333, freq=1300000000.0, eid="C.A5.4.8.L2")

# TDCavitys:
tdsb_428_b2 = TDCavity(l=1.5, freq=2800000000.0, eid="TDSB.428.B2")
tdsb_430_b2 = TDCavity(l=1.5, freq=2800000000.0, eid="TDSB.430.B2")

# Monitors:
bpma_233_b1 = Monitor(eid="BPMA.233.B1")
bpmc_249_l2 = Monitor(eid="BPMC.249.L2")
bpmc_261_l2 = Monitor(eid="BPMC.261.L2")
bpmr_273_l2 = Monitor(eid="BPMR.273.L2")
bpmr_285_l2 = Monitor(eid="BPMR.285.L2")
bpmc_297_l2 = Monitor(eid="BPMC.297.L2")
bpmr_309_l2 = Monitor(eid="BPMR.309.L2")
bpmc_321_l2 = Monitor(eid="BPMC.321.L2")
bpmr_333_l2 = Monitor(eid="BPMR.333.L2")
bpmc_345_l2 = Monitor(eid="BPMC.345.L2")
bpmc_357_l2 = Monitor(eid="BPMC.357.L2")
bpmr_369_l2 = Monitor(eid="BPMR.369.L2")
bpmr_381_l2 = Monitor(eid="BPMR.381.L2")
bpma_387_b2 = Monitor(eid="BPMA.387.B2")
bpma_390_b2 = Monitor(eid="BPMA.390.B2")
bpmf_393_b2 = Monitor(eid="BPMF.393.B2")
bpms_404_b2 = Monitor(eid="BPMS.404.B2")
bpmf_414_b2 = Monitor(eid="BPMF.414.B2")
bpma_418_b2 = Monitor(eid="BPMA.418.B2")
bpma_426_b2 = Monitor(eid="BPMA.426.B2")
bpma_432_b2 = Monitor(eid="BPMA.432.B2")
bpma_440_b2 = Monitor(eid="BPMA.440.B2")
bpma_444_b2 = Monitor(eid="BPMA.444.B2")
bpma_448_b2 = Monitor(eid="BPMA.448.B2")
bpma_452_b2 = Monitor(eid="BPMA.452.B2")
bpma_455_b2 = Monitor(eid="BPMA.455.B2")
bpma_459_b2 = Monitor(eid="BPMA.459.B2")
bpma_462_b2 = Monitor(eid="BPMA.462.B2")
bpma_465_b2 = Monitor(eid="BPMA.465.B2")

# Markers:
stsub_229_b1 = Marker(eid="STSUB.229.B1")
stgrd_231_b1 = Marker(eid="STGRD.231.B1")
tora_232_b1 = Marker(eid="TORA.232.B1")
enblock_233_b1 = Marker(eid="ENBLOCK.233.B1")
engrd_233_b1 = Marker(eid="ENGRD.233.B1")
stac_238_l2 = Marker(eid="STAC.238.L2")
enac_250_l2 = Marker(eid="ENAC.250.L2")
stac_250_l2 = Marker(eid="STAC.250.L2")
enac_262_l2 = Marker(eid="ENAC.262.L2")
stac_262_l2 = Marker(eid="STAC.262.L2")
enac_274_l2 = Marker(eid="ENAC.274.L2")
stac_274_l2 = Marker(eid="STAC.274.L2")
enac_286_l2 = Marker(eid="ENAC.286.L2")
stac_286_l2 = Marker(eid="STAC.286.L2")
enac_298_l2 = Marker(eid="ENAC.298.L2")
stac_298_l2 = Marker(eid="STAC.298.L2")
enac_310_l2 = Marker(eid="ENAC.310.L2")
mid_310_l2 = Marker(eid="MID.310.L2")
stac_310_l2 = Marker(eid="STAC.310.L2")
enac_322_l2 = Marker(eid="ENAC.322.L2")
stac_322_l2 = Marker(eid="STAC.322.L2")
enac_334_l2 = Marker(eid="ENAC.334.L2")
stac_334_l2 = Marker(eid="STAC.334.L2")
enac_346_l2 = Marker(eid="ENAC.346.L2")
stac_346_l2 = Marker(eid="STAC.346.L2")
enac_358_l2 = Marker(eid="ENAC.358.L2")
stac_358_l2 = Marker(eid="STAC.358.L2")
enac_370_l2 = Marker(eid="ENAC.370.L2")
stac_370_l2 = Marker(eid="STAC.370.L2")
enac_382_l2 = Marker(eid="ENAC.382.L2")
vcst78t40_385_l2 = Marker(eid="VCST78T40.385.L2")
ensec_385_l2 = Marker(eid="ENSEC.385.L2")
stsec_385_b2 = Marker(eid="STSEC.385.B2")
stsub_385_b2 = Marker(eid="STSUB.385.B2")
match_385_b2 = Marker(eid="MATCH.385.B2")
stgrd_386_b2 = Marker(eid="STGRD.386.B2")
tora_387_b2 = Marker(eid="TORA.387.B2")
dcm_388_b2 = Marker(eid="DCM.388.B2")
engrd_390_b2 = Marker(eid="ENGRD.390.B2")
stgrd_390_b2 = Marker(eid="STGRD.390.B2")
eod_390_b2 = Marker(eid="EOD.390.B2")
bcm_391_b2 = Marker(eid="BCM.391.B2")
otra_392_b2 = Marker(eid="OTRA.392.B2")
bam_392_b2 = Marker(eid="BAM.392.B2")
midbpmf_393_b2 = Marker(eid="MIDBPMF.393.B2")
engrd_393_b2 = Marker(eid="ENGRD.393.B2")
l2_stop_bc2_start = Marker(eid="l2_stop_bc2_start")
stlat_393_b2 = Marker(eid="STLAT.393.B2")
mbb_393a_b2 = Marker(eid="MBB.393a.B2")
mbb_393d_b2 = Marker(eid="MBB.393d.B2")
vcst40t400y_394_b2 = Marker(eid="VCST40T400Y.394.B2")
vcst400y_402_b2 = Marker(eid="VCST400Y.402.B2")
mbb_402a_b2 = Marker(eid="MBB.402a.B2")
mbb_402d_b2 = Marker(eid="MBB.402d.B2")
colo_403_b2 = Marker(eid="COLO.403.B2")
colu_403_b2 = Marker(eid="COLU.403.B2")
otrs_404_b2 = Marker(eid="OTRS.404.B2")
mbb_404a_b2 = Marker(eid="MBB.404a.B2")
mbb_404d_b2 = Marker(eid="MBB.404d.B2")
vcst400yt40_405_b2 = Marker(eid="VCST400YT40.405.B2")
srm_406_b2 = Marker(eid="SRM.406.B2")
vcst40y_413_b2 = Marker(eid="VCST40Y.413.B2")
mbb_413a_b2 = Marker(eid="MBB.413a.B2")
mbb_413d_b2 = Marker(eid="MBB.413d.B2")
match_414_b2 = Marker(eid="MATCH.414.B2")
enlat_414_b2 = Marker(eid="ENLAT.414.B2")
stlat_414_b2 = Marker(eid="STLAT.414.B2")
stgrd_414_b2 = Marker(eid="STGRD.414.B2")
bam_414_b2 = Marker(eid="BAM.414.B2")
midbpmf_414_b2 = Marker(eid="MIDBPMF.414.B2")
tora_415_b2 = Marker(eid="TORA.415.B2")
bc2_stop_b2d_start = Marker(eid="bc2_stop_b2d_start")
bcm_416_b2 = Marker(eid="BCM.416.B2")
engrd_419_b2 = Marker(eid="ENGRD.419.B2")
stgrd_423_b2 = Marker(eid="STGRD.423.B2")
otra_426_b2 = Marker(eid="OTRA.426.B2")
engrd_427_b2 = Marker(eid="ENGRD.427.B2")
stgrd_427_b2 = Marker(eid="STGRD.427.B2")
stblock_428_b2 = Marker(eid="STBLOCK.428.B2")
match_428_b2 = Marker(eid="MATCH.428.B2")
engrd_432_b2 = Marker(eid="ENGRD.432.B2")
stgrd_432_b2 = Marker(eid="STGRD.432.B2")
engrd_437_b2 = Marker(eid="ENGRD.437.B2")
stgrd_437_b2 = Marker(eid="STGRD.437.B2")
otra_438_b2 = Marker(eid="OTRA.438.B2")
engrd_442_b2 = Marker(eid="ENGRD.442.B2")
stgrd_442_b2 = Marker(eid="STGRD.442.B2")
match_446_b2 = Marker(eid="MATCH.446.B2")
chr_446_b2 = Marker(eid="CHR.446.B2")
engrd_446_b2 = Marker(eid="ENGRD.446.B2")
stgrd_447_b2 = Marker(eid="STGRD.447.B2")
otrb_450_b2 = Marker(eid="OTRB.450.B2")
engrd_451_b2 = Marker(eid="ENGRD.451.B2")
stgrd_451_b2 = Marker(eid="STGRD.451.B2")
otrb_454_b2 = Marker(eid="OTRB.454.B2")
engrd_456_b2 = Marker(eid="ENGRD.456.B2")
stgrd_456_b2 = Marker(eid="STGRD.456.B2")
otrb_457_b2 = Marker(eid="OTRB.457.B2")
engrd_461_b2 = Marker(eid="ENGRD.461.B2")
stgrd_461_b2 = Marker(eid="STGRD.461.B2")
otrb_461_b2 = Marker(eid="OTRB.461.B2")
crd_463_b2 = Marker(eid="CRD.463.B2")
engrd_466_b2 = Marker(eid="ENGRD.466.B2")
enlat_466_b2 = Marker(eid="ENLAT.466.B2")
ensub_466_b2 = Marker(eid="ENSUB.466.B2")

# Sequence:
cell = (
    stsub_229_b1,
    d_0,
    cbb_229_b1d,
    d_1,
    stgrd_231_b1,
    d_2,
    qd_231_b1,
    d_3,
    ccx_232_b1,
    d_4,
    tora_232_b1,
    d_5,
    qd_232_b1,
    d_6,
    ccy_232_b1,
    d_7,
    bpma_233_b1,
    d_8,
    qd_233_b1,
    d_9,
    enblock_233_b1,
    engrd_233_b1,
    d_10,
    vcdst_234_b1,
    d_11,
    cfb_236_l2,
    stac_238_l2,
    d_12,
    c_a3_1_1_l2,
    d_13,
    c_a3_1_2_l2,
    d_14,
    c_a3_1_3_l2,
    d_15,
    c_a3_1_4_l2,
    d_16,
    c_a3_1_5_l2,
    d_17,
    c_a3_1_6_l2,
    d_18,
    c_a3_1_7_l2,
    d_19,
    c_a3_1_8_l2,
    d_20,
    q_249_l2,
    d_21,
    cx_249_l2,
    d_22,
    bpmc_249_l2,
    d_23,
    enac_250_l2,
    stac_250_l2,
    d_24,
    c_a3_2_1_l2,
    d_25,
    c_a3_2_2_l2,
    d_26,
    c_a3_2_3_l2,
    d_27,
    c_a3_2_4_l2,
    d_28,
    c_a3_2_5_l2,
    d_29,
    c_a3_2_6_l2,
    d_30,
    c_a3_2_7_l2,
    d_31,
    c_a3_2_8_l2,
    d_32,
    q_261_l2,
    d_33,
    cy_261_l2,
    d_34,
    bpmc_261_l2,
    d_35,
    enac_262_l2,
    stac_262_l2,
    d_36,
    c_a3_3_1_l2,
    d_37,
    c_a3_3_2_l2,
    d_38,
    c_a3_3_3_l2,
    d_39,
    c_a3_3_4_l2,
    d_40,
    c_a3_3_5_l2,
    d_41,
    c_a3_3_6_l2,
    d_42,
    c_a3_3_7_l2,
    d_43,
    c_a3_3_8_l2,
    d_44,
    q_273_l2,
    d_45,
    cx_273_l2,
    d_46,
    bpmr_273_l2,
    d_47,
    enac_274_l2,
    stac_274_l2,
    d_48,
    c_a3_4_1_l2,
    d_49,
    c_a3_4_2_l2,
    d_50,
    c_a3_4_3_l2,
    d_51,
    c_a3_4_4_l2,
    d_52,
    c_a3_4_5_l2,
    d_53,
    c_a3_4_6_l2,
    d_54,
    c_a3_4_7_l2,
    d_55,
    c_a3_4_8_l2,
    d_56,
    q_285_l2,
    d_57,
    cy_285_l2,
    d_58,
    bpmr_285_l2,
    d_59,
    enac_286_l2,
    stac_286_l2,
    d_60,
    c_a4_1_1_l2,
    d_61,
    c_a4_1_2_l2,
    d_62,
    c_a4_1_3_l2,
    d_63,
    c_a4_1_4_l2,
    d_64,
    c_a4_1_5_l2,
    d_65,
    c_a4_1_6_l2,
    d_66,
    c_a4_1_7_l2,
    d_67,
    c_a4_1_8_l2,
    d_68,
    q_297_l2,
    d_69,
    cx_297_l2,
    d_70,
    bpmc_297_l2,
    d_71,
    enac_298_l2,
    stac_298_l2,
    d_72,
    c_a4_2_1_l2,
    d_73,
    c_a4_2_2_l2,
    d_74,
    c_a4_2_3_l2,
    d_75,
    c_a4_2_4_l2,
    d_76,
    c_a4_2_5_l2,
    d_77,
    c_a4_2_6_l2,
    d_78,
    c_a4_2_7_l2,
    d_79,
    c_a4_2_8_l2,
    d_80,
    q_309_l2,
    d_81,
    cy_309_l2,
    d_82,
    bpmr_309_l2,
    d_83,
    enac_310_l2,
    mid_310_l2,
    stac_310_l2,
    d_84,
    c_a4_3_1_l2,
    d_85,
    c_a4_3_2_l2,
    d_86,
    c_a4_3_3_l2,
    d_87,
    c_a4_3_4_l2,
    d_88,
    c_a4_3_5_l2,
    d_89,
    c_a4_3_6_l2,
    d_90,
    c_a4_3_7_l2,
    d_91,
    c_a4_3_8_l2,
    d_92,
    q_321_l2,
    d_93,
    cx_321_l2,
    d_94,
    bpmc_321_l2,
    d_95,
    enac_322_l2,
    stac_322_l2,
    d_96,
    c_a4_4_1_l2,
    d_97,
    c_a4_4_2_l2,
    d_98,
    c_a4_4_3_l2,
    d_99,
    c_a4_4_4_l2,
    d_100,
    c_a4_4_5_l2,
    d_101,
    c_a4_4_6_l2,
    d_102,
    c_a4_4_7_l2,
    d_103,
    c_a4_4_8_l2,
    d_104,
    q_333_l2,
    d_105,
    cy_333_l2,
    d_106,
    bpmr_333_l2,
    d_107,
    enac_334_l2,
    stac_334_l2,
    d_108,
    c_a5_1_1_l2,
    d_109,
    c_a5_1_2_l2,
    d_110,
    c_a5_1_3_l2,
    d_111,
    c_a5_1_4_l2,
    d_112,
    c_a5_1_5_l2,
    d_113,
    c_a5_1_6_l2,
    d_114,
    c_a5_1_7_l2,
    d_115,
    c_a5_1_8_l2,
    d_116,
    q_345_l2,
    d_117,
    cx_345_l2,
    d_118,
    bpmc_345_l2,
    d_119,
    enac_346_l2,
    stac_346_l2,
    d_120,
    c_a5_2_1_l2,
    d_121,
    c_a5_2_2_l2,
    d_122,
    c_a5_2_3_l2,
    d_123,
    c_a5_2_4_l2,
    d_124,
    c_a5_2_5_l2,
    d_125,
    c_a5_2_6_l2,
    d_126,
    c_a5_2_7_l2,
    d_127,
    c_a5_2_8_l2,
    d_128,
    q_357_l2,
    d_129,
    cy_357_l2,
    d_130,
    bpmc_357_l2,
    d_131,
    enac_358_l2,
    stac_358_l2,
    d_132,
    c_a5_3_1_l2,
    d_133,
    c_a5_3_2_l2,
    d_134,
    c_a5_3_3_l2,
    d_135,
    c_a5_3_4_l2,
    d_136,
    c_a5_3_5_l2,
    d_137,
    c_a5_3_6_l2,
    d_138,
    c_a5_3_7_l2,
    d_139,
    c_a5_3_8_l2,
    d_140,
    q_369_l2,
    d_141,
    cx_369_l2,
    d_142,
    bpmr_369_l2,
    d_143,
    enac_370_l2,
    stac_370_l2,
    d_144,
    c_a5_4_1_l2,
    d_145,
    c_a5_4_2_l2,
    d_146,
    c_a5_4_3_l2,
    d_147,
    c_a5_4_4_l2,
    d_148,
    c_a5_4_5_l2,
    d_149,
    c_a5_4_6_l2,
    d_150,
    c_a5_4_7_l2,
    d_151,
    c_a5_4_8_l2,
    d_152,
    q_381_l2,
    d_153,
    cy_381_l2,
    d_154,
    bpmr_381_l2,
    d_155,
    enac_382_l2,
    ctb_383_l2,
    vcst78t40_385_l2,
    ensec_385_l2,
    stsec_385_b2,
    stsub_385_b2,
    match_385_b2,
    d_156,
    vcdst_386_b2,
    d_157,
    stgrd_386_b2,
    d_158,
    tora_387_b2,
    d_159,
    bpma_387_b2,
    d_160,
    qd_387_b2,
    d_161,
    ccy_387_b2,
    d_162,
    dcm_388_b2,
    d_163,
    qd_388_b2,
    d_164,
    ccx_388_b2,
    d_165,
    engrd_390_b2,
    d_166,
    stgrd_390_b2,
    d_167,
    eod_390_b2,
    d_168,
    bpma_390_b2,
    d_169,
    qd_391_b2,
    d_170,
    ccy_391_b2,
    d_171,
    bcm_391_b2,
    d_172,
    qd_392_b2,
    d_173,
    ccx_392_b2,
    d_174,
    otra_392_b2,
    d_175,
    bam_392_b2,
    d_176,
    bpmf_393_b2,
    d_177,
    midbpmf_393_b2,
    d_178,
    engrd_393_b2,
    d_179,
    l2_stop_bc2_start,
    stlat_393_b2,
    d_180,
    mbb_393a_b2,
    d_181,
    bb_393_b2,
    d_182,
    mbb_393d_b2,
    vcst40t400y_394_b2,
    d_183,
    vcst400y_402_b2,
    mbb_402a_b2,
    d_184,
    bb_402_b2,
    d_185,
    mbb_402d_b2,
    cbb_403_b2,
    d_186,
    colo_403_b2,
    d_187,
    colu_403_b2,
    d_188,
    bpms_404_b2,
    d_189,
    otrs_404_b2,
    d_190,
    mbb_404a_b2,
    d_191,
    bb_404_b2,
    d_192,
    mbb_404d_b2,
    cbb_405_b2,
    vcst400yt40_405_b2,
    d_193,
    srm_406_b2,
    d_194,
    vcst40y_413_b2,
    mbb_413a_b2,
    d_195,
    bb_413_b2,
    d_196,
    mbb_413d_b2,
    cbb_414_b2,
    d_197,
    match_414_b2,
    enlat_414_b2,
    stlat_414_b2,
    d_198,
    stgrd_414_b2,
    d_199,
    bam_414_b2,
    d_200,
    bpmf_414_b2,
    d_201,
    midbpmf_414_b2,
    d_202,
    tora_415_b2,
    bc2_stop_b2d_start,
    d_203,
    qd_415_b2,
    d_204,
    ccx_415_b2,
    d_205,
    ccy_416_b2,
    d_206,
    bcm_416_b2,
    d_207,
    qd_417_b2,
    d_208,
    ccx_418_b2,
    d_209,
    ccy_418_b2,
    d_210,
    bpma_418_b2,
    d_211,
    qd_418_b2,
    d_212,
    engrd_419_b2,
    d_213,
    stgrd_423_b2,
    d_214,
    ccx_425_b2,
    d_215,
    qd_425_b2,
    d_216,
    ccy_426_b2,
    d_217,
    otra_426_b2,
    d_218,
    bpma_426_b2,
    d_219,
    qd_427_b2,
    d_220,
    engrd_427_b2,
    d_221,
    stgrd_427_b2,
    d_222,
    stblock_428_b2,
    match_428_b2,
    tdsb_428_b2,
    d_223,
    tdsb_430_b2,
    d_224,
    qd_431_b2,
    d_225,
    ccx_431_b2,
    d_226,
    bpma_432_b2,
    d_227,
    engrd_432_b2,
    d_228,
    stgrd_432_b2,
    d_229,
    qd_434_b2,
    d_230,
    ccy_434_b2,
    d_231,
    qd_437_b2,
    d_232,
    engrd_437_b2,
    d_233,
    stgrd_437_b2,
    d_234,
    otra_438_b2,
    d_235,
    bpma_440_b2,
    d_236,
    qd_440_b2,
    d_237,
    ccx_441_b2,
    d_238,
    engrd_442_b2,
    d_239,
    stgrd_442_b2,
    d_240,
    bpma_444_b2,
    d_241,
    qd_444_b2,
    d_242,
    kdy_445_b2,
    d_243,
    kdy_446_b2,
    d_244,
    match_446_b2,
    chr_446_b2,
    d_245,
    engrd_446_b2,
    d_246,
    stgrd_447_b2,
    d_247,
    ccx_447_b2,
    d_248,
    bpma_448_b2,
    d_249,
    qd_448_b2,
    d_250,
    ccy_448_b2,
    d_251,
    otrb_450_b2,
    d_252,
    engrd_451_b2,
    d_253,
    stgrd_451_b2,
    d_254,
    bpma_452_b2,
    d_255,
    qd_452_b2,
    d_256,
    kdy_452_b2,
    d_257,
    kdy_453_b2,
    d_258,
    otrb_454_b2,
    d_259,
    ccx_454_b2,
    d_260,
    bpma_455_b2,
    d_261,
    qd_456_b2,
    d_262,
    ccx_456_b2,
    d_263,
    engrd_456_b2,
    d_264,
    stgrd_456_b2,
    d_265,
    otrb_457_b2,
    d_266,
    bpma_459_b2,
    d_267,
    qd_459_b2,
    d_268,
    ccy_460_b2,
    d_269,
    engrd_461_b2,
    d_270,
    stgrd_461_b2,
    d_271,
    otrb_461_b2,
    d_272,
    bpma_462_b2,
    d_273,
    qd_463_b2,
    d_274,
    crd_463_b2,
    d_275,
    ccx_464_b2,
    d_276,
    qd_464_b2,
    d_277,
    ccy_464_b2,
    d_278,
    ccx_465_b2,
    d_279,
    bpma_465_b2,
    d_280,
    qd_465_b2,
    d_281,
    engrd_466_b2,
    d_282,
    enlat_466_b2,
    ensub_466_b2,
)

# Power Supply IDs:
# Quadrupole power supplies:
qd_231_b1.ps_id = "QD.20.B1"
qd_232_b1.ps_id = "QD.21.B1"
qd_233_b1.ps_id = "QD.22.B1"
q_249_l2.ps_id = "Q.A3.1.L2"
q_261_l2.ps_id = "Q.A3.2.L2"
q_273_l2.ps_id = "Q.A3.3.L2"
q_285_l2.ps_id = "Q.A3.4.L2"
q_297_l2.ps_id = "Q.A4.1.L2"
q_309_l2.ps_id = "Q.A4.2.L2"
q_321_l2.ps_id = "Q.A4.3.L2"
q_333_l2.ps_id = "Q.A4.4.L2"
q_345_l2.ps_id = "Q.A5.1.L2"
q_357_l2.ps_id = "Q.A5.2.L2"
q_369_l2.ps_id = "Q.A5.3.L2"
q_381_l2.ps_id = "Q.A5.4.L2"
qd_387_b2.ps_id = "QD.1.B2"
qd_388_b2.ps_id = "QD.2.B2"
qd_391_b2.ps_id = "QD.3.B2"
qd_392_b2.ps_id = "QD.4.B2"
qd_415_b2.ps_id = "QD.6.B2"
qd_417_b2.ps_id = "QD.7.B2"
qd_418_b2.ps_id = "QD.8.B2"
qd_425_b2.ps_id = "QD.9.B2"
qd_427_b2.ps_id = "QD.10.B2"
qd_431_b2.ps_id = "QD.11.B2"
qd_434_b2.ps_id = "QD.12.B2"
qd_437_b2.ps_id = "QD.13.B2"
qd_440_b2.ps_id = "QD.14.B2"
qd_444_b2.ps_id = "QD.15.B2"
qd_448_b2.ps_id = "QD.16.B2"
qd_452_b2.ps_id = "QD.17.B2"
qd_456_b2.ps_id = "QD.18.B2"
qd_459_b2.ps_id = "QD.19.B2"
qd_463_b2.ps_id = "QD.21.B2"
qd_464_b2.ps_id = "QD.22.B2"
qd_465_b2.ps_id = "QD.23.B2"

# SBend power supplies:
bb_393_b2.ps_id = "BB.1.B2"
bb_402_b2.ps_id = "BB.1.B2"
bb_404_b2.ps_id = "BB.1.B2"
bb_413_b2.ps_id = "BB.1.B2"

# RBend power supplies:
kdy_445_b2.ps_id = "KDY.445.B2"
kdy_446_b2.ps_id = "KDY.446.B2"
kdy_452_b2.ps_id = "KDY.452.B2"
kdy_453_b2.ps_id = "KDY.453.B2"

# Cavity power supplies:
c_a3_1_1_l2.ps_id = "C.A3.L2"
c_a3_1_2_l2.ps_id = "C.A3.L2"
c_a3_1_3_l2.ps_id = "C.A3.L2"
c_a3_1_4_l2.ps_id = "C.A3.L2"
c_a3_1_5_l2.ps_id = "C.A3.L2"
c_a3_1_6_l2.ps_id = "C.A3.L2"
c_a3_1_7_l2.ps_id = "C.A3.L2"
c_a3_1_8_l2.ps_id = "C.A3.L2"
c_a3_2_1_l2.ps_id = "C.A3.L2"
c_a3_2_2_l2.ps_id = "C.A3.L2"
c_a3_2_3_l2.ps_id = "C.A3.L2"
c_a3_2_4_l2.ps_id = "C.A3.L2"
c_a3_2_5_l2.ps_id = "C.A3.L2"
c_a3_2_6_l2.ps_id = "C.A3.L2"
c_a3_2_7_l2.ps_id = "C.A3.L2"
c_a3_2_8_l2.ps_id = "C.A3.L2"
c_a3_3_1_l2.ps_id = "C.A3.L2"
c_a3_3_2_l2.ps_id = "C.A3.L2"
c_a3_3_3_l2.ps_id = "C.A3.L2"
c_a3_3_4_l2.ps_id = "C.A3.L2"
c_a3_3_5_l2.ps_id = "C.A3.L2"
c_a3_3_6_l2.ps_id = "C.A3.L2"
c_a3_3_7_l2.ps_id = "C.A3.L2"
c_a3_3_8_l2.ps_id = "C.A3.L2"
c_a3_4_1_l2.ps_id = "C.A3.L2"
c_a3_4_2_l2.ps_id = "C.A3.L2"
c_a3_4_3_l2.ps_id = "C.A3.L2"
c_a3_4_4_l2.ps_id = "C.A3.L2"
c_a3_4_5_l2.ps_id = "C.A3.L2"
c_a3_4_6_l2.ps_id = "C.A3.L2"
c_a3_4_7_l2.ps_id = "C.A3.L2"
c_a3_4_8_l2.ps_id = "C.A3.L2"
c_a4_1_1_l2.ps_id = "C.A4.L2"
c_a4_1_2_l2.ps_id = "C.A4.L2"
c_a4_1_3_l2.ps_id = "C.A4.L2"
c_a4_1_4_l2.ps_id = "C.A4.L2"
c_a4_1_5_l2.ps_id = "C.A4.L2"
c_a4_1_6_l2.ps_id = "C.A4.L2"
c_a4_1_7_l2.ps_id = "C.A4.L2"
c_a4_1_8_l2.ps_id = "C.A4.L2"
c_a4_2_1_l2.ps_id = "C.A4.L2"
c_a4_2_2_l2.ps_id = "C.A4.L2"
c_a4_2_3_l2.ps_id = "C.A4.L2"
c_a4_2_4_l2.ps_id = "C.A4.L2"
c_a4_2_5_l2.ps_id = "C.A4.L2"
c_a4_2_6_l2.ps_id = "C.A4.L2"
c_a4_2_7_l2.ps_id = "C.A4.L2"
c_a4_2_8_l2.ps_id = "C.A4.L2"
c_a4_3_1_l2.ps_id = "C.A4.L2"
c_a4_3_2_l2.ps_id = "C.A4.L2"
c_a4_3_3_l2.ps_id = "C.A4.L2"
c_a4_3_4_l2.ps_id = "C.A4.L2"
c_a4_3_5_l2.ps_id = "C.A4.L2"
c_a4_3_6_l2.ps_id = "C.A4.L2"
c_a4_3_7_l2.ps_id = "C.A4.L2"
c_a4_3_8_l2.ps_id = "C.A4.L2"
c_a4_4_1_l2.ps_id = "C.A4.L2"
c_a4_4_2_l2.ps_id = "C.A4.L2"
c_a4_4_3_l2.ps_id = "C.A4.L2"
c_a4_4_4_l2.ps_id = "C.A4.L2"
c_a4_4_5_l2.ps_id = "C.A4.L2"
c_a4_4_6_l2.ps_id = "C.A4.L2"
c_a4_4_7_l2.ps_id = "C.A4.L2"
c_a4_4_8_l2.ps_id = "C.A4.L2"
c_a5_1_1_l2.ps_id = "C.A5.L2"
c_a5_1_2_l2.ps_id = "C.A5.L2"
c_a5_1_3_l2.ps_id = "C.A5.L2"
c_a5_1_4_l2.ps_id = "C.A5.L2"
c_a5_1_5_l2.ps_id = "C.A5.L2"
c_a5_1_6_l2.ps_id = "C.A5.L2"
c_a5_1_7_l2.ps_id = "C.A5.L2"
c_a5_1_8_l2.ps_id = "C.A5.L2"
c_a5_2_1_l2.ps_id = "C.A5.L2"
c_a5_2_2_l2.ps_id = "C.A5.L2"
c_a5_2_3_l2.ps_id = "C.A5.L2"
c_a5_2_4_l2.ps_id = "C.A5.L2"
c_a5_2_5_l2.ps_id = "C.A5.L2"
c_a5_2_6_l2.ps_id = "C.A5.L2"
c_a5_2_7_l2.ps_id = "C.A5.L2"
c_a5_2_8_l2.ps_id = "C.A5.L2"
c_a5_3_1_l2.ps_id = "C.A5.L2"
c_a5_3_2_l2.ps_id = "C.A5.L2"
c_a5_3_3_l2.ps_id = "C.A5.L2"
c_a5_3_4_l2.ps_id = "C.A5.L2"
c_a5_3_5_l2.ps_id = "C.A5.L2"
c_a5_3_6_l2.ps_id = "C.A5.L2"
c_a5_3_7_l2.ps_id = "C.A5.L2"
c_a5_3_8_l2.ps_id = "C.A5.L2"
c_a5_4_1_l2.ps_id = "C.A5.L2"
c_a5_4_2_l2.ps_id = "C.A5.L2"
c_a5_4_3_l2.ps_id = "C.A5.L2"
c_a5_4_4_l2.ps_id = "C.A5.L2"
c_a5_4_5_l2.ps_id = "C.A5.L2"
c_a5_4_6_l2.ps_id = "C.A5.L2"
c_a5_4_7_l2.ps_id = "C.A5.L2"
c_a5_4_8_l2.ps_id = "C.A5.L2"
