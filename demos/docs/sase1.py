from ocelot import * 

#Initial Twiss parameters
tws0 = Twiss()
tws0._E = 14
tws0._beta_x = 10.776023018690426
tws0._beta_y = 39.48730762141566
tws0._alpha_x = 0.6698245190372976
tws0._alpha_y = -2.0320682898456615
tws0.s = 2025.3865970000204

# Drifts
id_35371550_ = Drift(l=1.472401, eid='ID_35371550_')
d_3 = Drift(l=13.997421, eid='D_3')
d_4 = Drift(l=0.208951, eid='D_4')
d_5 = Drift(l=0.153951, eid='D_5')
d_6 = Drift(l=11.37, eid='D_6')
d_7 = Drift(l=0.17, eid='D_7')
d_8 = Drift(l=0.27, eid='D_8')
d_9 = Drift(l=1.742401, eid='D_9')
id_23741906_ = Drift(l=13.997441, eid='ID_23741906_')
d_12 = Drift(l=0.20895, eid='D_12')
d_13 = Drift(l=0.15395, eid='D_13')
d_14 = Drift(l=6.405, eid='D_14')
d_15 = Drift(l=6.6, eid='D_15')
d_18 = Drift(l=14.005, eid='D_18')
d_21 = Drift(l=14.0051, eid='D_21')
id_6624594_ = Drift(l=14.705, eid='ID_6624594_')
d_28 = Drift(l=5.53795, eid='D_28')
d_29 = Drift(l=0.153944, eid='D_29')
d_30 = Drift(l=2.297556, eid='D_30')
d_31 = Drift(l=1.85645, eid='D_31')
d_34 = Drift(l=15.205, eid='D_34')
d_37 = Drift(l=0.153946, eid='D_37')
id_52103571_ = Drift(l=13.551054, eid='ID_52103571_')
d_42 = Drift(l=0.4075, eid='D_42')
d_43 = Drift(l=6.029, eid='D_43')
d_44 = Drift(l=6.73645, eid='D_44')
d_47 = Drift(l=0.153948, eid='D_47')
d_48 = Drift(l=0.25395, eid='D_48')
id_38309706_ = Drift(l=2.765102, eid='ID_38309706_')
d_52 = Drift(l=0.21, eid='D_52')
id_32688183_ = Drift(l=5.2620000000000005, eid='ID_32688183_')
d_63 = Drift(l=0.253949, eid='D_63')
id_71170601_ = Drift(l=9.621051, eid='ID_71170601_')
id_62721345_ = Drift(l=3.866817, eid='ID_62721345_')
d_70 = Drift(l=5.772483, eid='D_70')
d_71 = Drift(l=0.19975, eid='D_71')
d_72 = Drift(l=0.21815, eid='D_72')
d_73 = Drift(l=0.17815, eid='D_73')
d_74 = Drift(l=0.125, eid='D_74')
d_75 = Drift(l=2.1675, eid='D_75')
d_76 = Drift(l=1.8267, eid='D_76')
d_77 = Drift(l=0.15, eid='D_77')
d_78 = Drift(l=0.4323, eid='D_78')
d_79 = Drift(l=0.18665, eid='D_79')
id_85978457_ = Drift(l=5.79965, eid='ID_85978457_')
id_85388498_ = Drift(l=0.35787, eid='ID_85388498_')
d_86 = Drift(l=0.07278, eid='D_86')
d_87 = Drift(l=0.0717, eid='D_87')
d_88 = Drift(l=0.2973, eid='D_88')
id_83891648_ = Drift(l=0.33287, eid='ID_83891648_')
id_54441912_ = Drift(l=0.025, eid='ID_54441912_')
id_17219612_ = Drift(l=0.57565, eid='ID_17219612_')
id_53856480_ = Drift(l=1.025, eid='ID_53856480_')
id_85095730_ = Drift(l=1.45, eid='ID_85095730_')
d_276 = Drift(l=0.524, eid='D_276')
id_42404360_ = Drift(l=0.04015, eid='ID_42404360_')

# Quadrupoles
qk_2027_tl = Quadrupole(l=1.0552, k1=-0.09035960007960576, eid='QK.2027.TL')
qf_2042_tl = Quadrupole(l=0.5321, k1=0.17919084760007517, eid='QF.2042.TL')
qk_2057_tl = Quadrupole(l=1.0552, k1=-0.09035960007960576, eid='QK.2057.TL')
qf_2072_t2 = Quadrupole(l=0.5321, k1=0.17919084760007517, eid='QF.2072.T2')
qf_2087_t2 = Quadrupole(l=0.5321, k1=-0.17919084760007517, eid='QF.2087.T2')
qf_2102_t2 = Quadrupole(l=0.5321, k1=0.17919084760007517, eid='QF.2102.T2')
qf_2117_t2 = Quadrupole(l=0.5321, k1=-0.17919084760007517, eid='QF.2117.T2')
qf_2132_t2 = Quadrupole(l=0.5321, k1=0.17919084760007517, eid='QF.2132.T2')
qh_2139_t2 = Quadrupole(l=1.0291, eid='QH.2139.T2')
qf_2145_t2 = Quadrupole(l=0.5321, k1=-0.17919084760007517, eid='QF.2145.T2')
qf_2162_t2 = Quadrupole(l=0.5321, k1=0.17919084760007517, eid='QF.2162.T2')
qf_2177_t2 = Quadrupole(l=0.5321, k1=-0.15508031819958654, eid='QF.2177.T2')
qf_2184_t2 = Quadrupole(l=0.5321, eid='QF.2184.T2')
qf_2192_t2 = Quadrupole(l=0.5321, k1=0.152454154899455, eid='QF.2192.T2')
qf_2207_t2 = Quadrupole(l=0.5321, k1=-0.12867168310092086, eid='QF.2207.T2')
qf_2218_t2 = Quadrupole(l=0.5321, k1=0.12061272159932343, eid='QF.2218.T2')
qa_2229_t2 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2229.T2')
qa_2235_t2 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2235.T2')
qa_2241_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2241.SA1')
qa_2247_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2247.SA1')
qa_2253_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2253.SA1')
qa_2259_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2259.SA1')
qa_2266_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2266.SA1')
qa_2272_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2272.SA1')
qa_2278_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2278.SA1')
qa_2284_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2284.SA1')
qa_2290_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2290.SA1')
qa_2296_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2296.SA1')
qa_2302_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2302.SA1')
qa_2308_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2308.SA1')
qa_2314_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2314.SA1')
qa_2320_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2320.SA1')
qa_2327_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2327.SA1')
qa_2333_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2333.SA1')
qa_2339_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2339.SA1')
qa_2345_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2345.SA1')
qa_2351_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2351.SA1')
qa_2357_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2357.SA1')
qa_2363_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2363.SA1')
qa_2369_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2369.SA1')
qa_2375_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2375.SA1')
qa_2381_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2381.SA1')
qa_2388_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2388.SA1')
qa_2394_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2394.SA1')
qa_2400_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2400.SA1')
qa_2406_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2406.SA1')
qa_2412_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2412.SA1')
qa_2418_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2418.SA1')
qa_2424_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2424.SA1')
qa_2430_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2430.SA1')
qa_2436_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2436.SA1')
qa_2442_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2442.SA1')
qa_2449_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2449.SA1')
qa_2455_sa1 = Quadrupole(l=0.1137, k1=0.5620550638962181, eid='QA.2455.SA1')
qa_2461_sa1 = Quadrupole(l=0.1137, k1=-0.5620550638962181, eid='QA.2461.SA1')

# SBends
bp_2197i_t2 = SBend(l=0.44, angle=-5e-07, e2=-5e-07, eid='BP.2197I.T2')
bp_2197ii_t2 = SBend(l=0.44, angle=1e-06, e1=5e-07, e2=-5e-07, eid='BP.2197II.T2')
bp_2198_t2 = SBend(l=0.44, angle=-1e-06, e1=-5e-07, e2=5e-07, eid='BP.2198.T2')
bp_2199i_t2 = SBend(l=0.44, angle=1e-06, e1=5e-07, e2=-5e-07, eid='BP.2199I.T2')
bp_2199ii_t2 = SBend(l=0.44, angle=-1e-06, e1=-5e-07, e2=5e-07, eid='BP.2199II.T2')
bp_2200_t2 = SBend(l=0.44, angle=1e-06, e1=5e-07, e2=-5e-07, eid='BP.2200.T2')
bp_2201_t2 = SBend(l=0.44, angle=-5e-07, e1=5e-07, eid='BP.2201.T2')
bs_2431_sa1 = SBend(l=0.3, eid='BS.2431.SA1')
bs_2432_sa1 = SBend(l=0.3, eid='BS.2432.SA1')
bs_2434_sa1 = SBend(l=0.3, eid='BS.2434.SA1')
bs_2435_sa1 = SBend(l=0.3, eid='BS.2435.SA1')

# RBends
bd_2079_t2 = RBend(l=1.0, eid='BD.2079.T2')
kl_2143_t2 = RBend(l=0.93, eid='KL.2143.T2')
kl_2193_t2 = RBend(l=0.93, eid='KL.2193.T2')
kl_2208_t2 = RBend(l=0.93, eid='KL.2208.T2')

# Hcors
cfx_2042_tl = Hcor(l=0.1, eid='CFX.2042.TL')
chx_2054_tl = Hcor(l=0.2, eid='CHX.2054.TL')
cfx_2072_t2 = Hcor(l=0.1, eid='CFX.2072.T2')
cfx_2102_t2 = Hcor(l=0.1, eid='CFX.2102.T2')
cfx_2133_t2 = Hcor(l=0.1, eid='CFX.2133.T2')
cmx_2140_t2 = Hcor(l=0.3, eid='CMX.2140.T2')
cfx_2162_t2 = Hcor(l=0.1, eid='CFX.2162.T2')
cmx_2178_t2 = Hcor(l=0.3, eid='CMX.2178.T2')
cfx_2192_t2 = Hcor(l=0.1, eid='CFX.2192.T2')
cbp_2197i_t2 = Hcor(eid='CBP.2197I.T2')
cbp_2201_t2 = Hcor(eid='CBP.2201.T2')
cfx_2219_t2 = Hcor(l=0.1, eid='CFX.2219.T2')
cex_2230_t2 = Hcor(l=0.1, eid='CEX.2230.T2')
cex_2234_t2 = Hcor(l=0.1, eid='CEX.2234.T2')
cax_2242_sa1 = Hcor(eid='CAX.2242.SA1')
cbx_2247_sa1 = Hcor(eid='CBX.2247.SA1')
cax_2248_sa1 = Hcor(eid='CAX.2248.SA1')
cbx_2253_sa1 = Hcor(eid='CBX.2253.SA1')
cax_2254_sa1 = Hcor(eid='CAX.2254.SA1')
cbx_2259_sa1 = Hcor(eid='CBX.2259.SA1')
cax_2260_sa1 = Hcor(eid='CAX.2260.SA1')
cbx_2265_sa1 = Hcor(eid='CBX.2265.SA1')
cax_2267_sa1 = Hcor(eid='CAX.2267.SA1')
cbx_2271_sa1 = Hcor(eid='CBX.2271.SA1')
cax_2273_sa1 = Hcor(eid='CAX.2273.SA1')
cbx_2277_sa1 = Hcor(eid='CBX.2277.SA1')
cax_2279_sa1 = Hcor(eid='CAX.2279.SA1')
cbx_2283_sa1 = Hcor(eid='CBX.2283.SA1')
cax_2285_sa1 = Hcor(eid='CAX.2285.SA1')
cbx_2289_sa1 = Hcor(eid='CBX.2289.SA1')
cax_2291_sa1 = Hcor(eid='CAX.2291.SA1')
cbx_2296_sa1 = Hcor(eid='CBX.2296.SA1')
cax_2297_sa1 = Hcor(eid='CAX.2297.SA1')
cbx_2302_sa1 = Hcor(eid='CBX.2302.SA1')
cax_2303_sa1 = Hcor(eid='CAX.2303.SA1')
cbx_2308_sa1 = Hcor(eid='CBX.2308.SA1')
cax_2309_sa1 = Hcor(eid='CAX.2309.SA1')
cbx_2314_sa1 = Hcor(eid='CBX.2314.SA1')
cax_2315_sa1 = Hcor(eid='CAX.2315.SA1')
cbx_2320_sa1 = Hcor(eid='CBX.2320.SA1')
cax_2321_sa1 = Hcor(eid='CAX.2321.SA1')
cbx_2326_sa1 = Hcor(eid='CBX.2326.SA1')
cax_2328_sa1 = Hcor(eid='CAX.2328.SA1')
cbx_2332_sa1 = Hcor(eid='CBX.2332.SA1')
cax_2334_sa1 = Hcor(eid='CAX.2334.SA1')
cbx_2338_sa1 = Hcor(eid='CBX.2338.SA1')
cax_2340_sa1 = Hcor(eid='CAX.2340.SA1')
cbx_2344_sa1 = Hcor(eid='CBX.2344.SA1')
cax_2346_sa1 = Hcor(eid='CAX.2346.SA1')
cbx_2350_sa1 = Hcor(eid='CBX.2350.SA1')
cax_2352_sa1 = Hcor(eid='CAX.2352.SA1')
cbx_2357_sa1 = Hcor(eid='CBX.2357.SA1')
cax_2358_sa1 = Hcor(eid='CAX.2358.SA1')
cbx_2363_sa1 = Hcor(eid='CBX.2363.SA1')
cax_2364_sa1 = Hcor(eid='CAX.2364.SA1')
cbx_2369_sa1 = Hcor(eid='CBX.2369.SA1')
cax_2370_sa1 = Hcor(eid='CAX.2370.SA1')
cbx_2375_sa1 = Hcor(eid='CBX.2375.SA1')
cax_2376_sa1 = Hcor(eid='CAX.2376.SA1')
cbx_2381_sa1 = Hcor(eid='CBX.2381.SA1')
cax_2382_sa1 = Hcor(eid='CAX.2382.SA1')
cbx_2387_sa1 = Hcor(eid='CBX.2387.SA1')
cax_2389_sa1 = Hcor(eid='CAX.2389.SA1')
cbx_2393_sa1 = Hcor(eid='CBX.2393.SA1')
cax_2395_sa1 = Hcor(eid='CAX.2395.SA1')
cbx_2399_sa1 = Hcor(eid='CBX.2399.SA1')
cax_2401_sa1 = Hcor(eid='CAX.2401.SA1')
cbx_2405_sa1 = Hcor(eid='CBX.2405.SA1')
cax_2407_sa1 = Hcor(eid='CAX.2407.SA1')
cbx_2411_sa1 = Hcor(eid='CBX.2411.SA1')
cax_2413_sa1 = Hcor(eid='CAX.2413.SA1')
cbx_2418_sa1 = Hcor(eid='CBX.2418.SA1')
cax_2419_sa1 = Hcor(eid='CAX.2419.SA1')
cbx_2424_sa1 = Hcor(eid='CBX.2424.SA1')
cax_2425_sa1 = Hcor(eid='CAX.2425.SA1')
cbx_2430_sa1 = Hcor(eid='CBX.2430.SA1')
cbs_2432_sa1 = Hcor(eid='CBS.2432.SA1')
cbs_2434_sa1 = Hcor(eid='CBS.2434.SA1')
cbs_2435_sa1 = Hcor(eid='CBS.2435.SA1')
cax_2437_sa1 = Hcor(eid='CAX.2437.SA1')
cbx_2442_sa1 = Hcor(eid='CBX.2442.SA1')
cax_2443_sa1 = Hcor(eid='CAX.2443.SA1')
cbx_2448_sa1 = Hcor(eid='CBX.2448.SA1')
cax_2450_sa1 = Hcor(eid='CAX.2450.SA1')
cbx_2454_sa1 = Hcor(eid='CBX.2454.SA1')
cax_2456_sa1 = Hcor(eid='CAX.2456.SA1')
cbx_2460_sa1 = Hcor(eid='CBX.2460.SA1')

# Vcors
chy_2054_tl = Vcor(l=0.2, eid='CHY.2054.TL')
cfy_2087_t2 = Vcor(l=0.1, eid='CFY.2087.T2')
cfy_2117_t2 = Vcor(l=0.1, eid='CFY.2117.T2')
cfy_2146_t2 = Vcor(l=0.1, eid='CFY.2146.T2')
cmy_2162_t2 = Vcor(l=0.3, eid='CMY.2162.T2')
cfy_2177_t2 = Vcor(l=0.1, eid='CFY.2177.T2')
cmy_2192_t2 = Vcor(l=0.3, eid='CMY.2192.T2')
cfy_2207_t2 = Vcor(l=0.1, eid='CFY.2207.T2')
cny_2229_t2 = Vcor(l=0.3, eid='CNY.2229.T2')
cny_2234_t2 = Vcor(l=0.3, eid='CNY.2234.T2')
cay_2242_sa1 = Vcor(eid='CAY.2242.SA1')
cby_2247_sa1 = Vcor(eid='CBY.2247.SA1')
cay_2248_sa1 = Vcor(eid='CAY.2248.SA1')
cby_2253_sa1 = Vcor(eid='CBY.2253.SA1')
cay_2254_sa1 = Vcor(eid='CAY.2254.SA1')
cby_2259_sa1 = Vcor(eid='CBY.2259.SA1')
cay_2260_sa1 = Vcor(eid='CAY.2260.SA1')
cby_2265_sa1 = Vcor(eid='CBY.2265.SA1')
cay_2267_sa1 = Vcor(eid='CAY.2267.SA1')
cby_2271_sa1 = Vcor(eid='CBY.2271.SA1')
cay_2273_sa1 = Vcor(eid='CAY.2273.SA1')
cby_2277_sa1 = Vcor(eid='CBY.2277.SA1')
cay_2279_sa1 = Vcor(eid='CAY.2279.SA1')
cby_2283_sa1 = Vcor(eid='CBY.2283.SA1')
cay_2285_sa1 = Vcor(eid='CAY.2285.SA1')
cby_2289_sa1 = Vcor(eid='CBY.2289.SA1')
cay_2291_sa1 = Vcor(eid='CAY.2291.SA1')
cby_2296_sa1 = Vcor(eid='CBY.2296.SA1')
cay_2297_sa1 = Vcor(eid='CAY.2297.SA1')
cby_2302_sa1 = Vcor(eid='CBY.2302.SA1')
cay_2303_sa1 = Vcor(eid='CAY.2303.SA1')
cby_2308_sa1 = Vcor(eid='CBY.2308.SA1')
cay_2309_sa1 = Vcor(eid='CAY.2309.SA1')
cby_2314_sa1 = Vcor(eid='CBY.2314.SA1')
cay_2315_sa1 = Vcor(eid='CAY.2315.SA1')
cby_2320_sa1 = Vcor(eid='CBY.2320.SA1')
cay_2321_sa1 = Vcor(eid='CAY.2321.SA1')
cby_2326_sa1 = Vcor(eid='CBY.2326.SA1')
cay_2328_sa1 = Vcor(eid='CAY.2328.SA1')
cby_2332_sa1 = Vcor(eid='CBY.2332.SA1')
cay_2334_sa1 = Vcor(eid='CAY.2334.SA1')
cby_2338_sa1 = Vcor(eid='CBY.2338.SA1')
cay_2340_sa1 = Vcor(eid='CAY.2340.SA1')
cby_2344_sa1 = Vcor(eid='CBY.2344.SA1')
cay_2346_sa1 = Vcor(eid='CAY.2346.SA1')
cby_2350_sa1 = Vcor(eid='CBY.2350.SA1')
cay_2352_sa1 = Vcor(eid='CAY.2352.SA1')
cby_2357_sa1 = Vcor(eid='CBY.2357.SA1')
cay_2358_sa1 = Vcor(eid='CAY.2358.SA1')
cby_2363_sa1 = Vcor(eid='CBY.2363.SA1')
cay_2364_sa1 = Vcor(eid='CAY.2364.SA1')
cby_2369_sa1 = Vcor(eid='CBY.2369.SA1')
cay_2370_sa1 = Vcor(eid='CAY.2370.SA1')
cby_2375_sa1 = Vcor(eid='CBY.2375.SA1')
cay_2376_sa1 = Vcor(eid='CAY.2376.SA1')
cby_2381_sa1 = Vcor(eid='CBY.2381.SA1')
cay_2382_sa1 = Vcor(eid='CAY.2382.SA1')
cby_2387_sa1 = Vcor(eid='CBY.2387.SA1')
cay_2389_sa1 = Vcor(eid='CAY.2389.SA1')
cby_2393_sa1 = Vcor(eid='CBY.2393.SA1')
cay_2395_sa1 = Vcor(eid='CAY.2395.SA1')
cby_2399_sa1 = Vcor(eid='CBY.2399.SA1')
cay_2401_sa1 = Vcor(eid='CAY.2401.SA1')
cby_2405_sa1 = Vcor(eid='CBY.2405.SA1')
cay_2407_sa1 = Vcor(eid='CAY.2407.SA1')
cby_2411_sa1 = Vcor(eid='CBY.2411.SA1')
cay_2413_sa1 = Vcor(eid='CAY.2413.SA1')
cby_2418_sa1 = Vcor(eid='CBY.2418.SA1')
cay_2419_sa1 = Vcor(eid='CAY.2419.SA1')
cby_2424_sa1 = Vcor(eid='CBY.2424.SA1')
cay_2425_sa1 = Vcor(eid='CAY.2425.SA1')
cby_2430_sa1 = Vcor(eid='CBY.2430.SA1')
cay_2437_sa1 = Vcor(eid='CAY.2437.SA1')
cby_2442_sa1 = Vcor(eid='CBY.2442.SA1')
cay_2443_sa1 = Vcor(eid='CAY.2443.SA1')
cby_2448_sa1 = Vcor(eid='CBY.2448.SA1')
cay_2450_sa1 = Vcor(eid='CAY.2450.SA1')
cby_2454_sa1 = Vcor(eid='CBY.2454.SA1')
cay_2456_sa1 = Vcor(eid='CAY.2456.SA1')
cby_2460_sa1 = Vcor(eid='CBY.2460.SA1')

# Undulators
u40s_2232_t2 = Undulator(lperiod=0.04, nperiods=3.0, eid='U40S.2232.T2')
u40_2244_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2244.SA1')
u40_2250_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2250.SA1')
u40_2256_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2256.SA1')
u40_2262_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2262.SA1')
u40_2269_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2269.SA1')
u40_2275_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2275.SA1')
u40_2281_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2281.SA1')
u40_2287_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2287.SA1')
u40_2293_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2293.SA1')
u40_2299_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2299.SA1')
u40_2305_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2305.SA1')
u40_2311_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2311.SA1')
u40_2317_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2317.SA1')
u40_2323_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2323.SA1')
u40_2330_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2330.SA1')
u40_2336_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2336.SA1')
u40_2342_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2342.SA1')
u40_2348_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2348.SA1')
u40_2354_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2354.SA1')
u40_2360_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2360.SA1')
u40_2366_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2366.SA1')
u40_2372_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2372.SA1')
u40_2378_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2378.SA1')
u40_2384_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2384.SA1')
u40_2391_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2391.SA1')
u40_2397_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2397.SA1')
u40_2403_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2403.SA1')
u40_2409_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2409.SA1')
u40_2415_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2415.SA1')
u40_2421_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2421.SA1')
u40_2427_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2427.SA1')
u40_2439_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2439.SA1')
u40_2445_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2445.SA1')
u40_2452_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2452.SA1')
u40_2458_sa1 = Undulator(lperiod=0.04, nperiods=125.0, eid='U40.2458.SA1')

# Monitors
bpma_2041_tl = Monitor(eid='BPMA.2041.TL')
bpma_2054_tl = Monitor(eid='BPMA.2054.TL')
bpma_2071_t2 = Monitor(eid='BPMA.2071.T2')
bpma_2086_t2 = Monitor(eid='BPMA.2086.T2')
bpma_2101_t2 = Monitor(eid='BPMA.2101.T2')
bpma_2116_t2 = Monitor(eid='BPMA.2116.T2')
bpma_2132_t2 = Monitor(eid='BPMA.2132.T2')
bpma_2145_t2 = Monitor(eid='BPMA.2145.T2')
bpma_2161_t2 = Monitor(eid='BPMA.2161.T2')
bpma_2176_t2 = Monitor(eid='BPMA.2176.T2')
bpma_2191_t2 = Monitor(eid='BPMA.2191.T2')
bpma_2206_t2 = Monitor(eid='BPMA.2206.T2')
bpma_2218_t2 = Monitor(eid='BPMA.2218.T2')
bpma_2223_t2 = Monitor(eid='BPMA.2223.T2')
bpme_2229_t2 = Monitor(eid='BPME.2229.T2')
bpme_2235_t2 = Monitor(eid='BPME.2235.T2')
bpme_2241_sa1 = Monitor(eid='BPME.2241.SA1')
bpme_2247_sa1 = Monitor(eid='BPME.2247.SA1')
bpme_2253_sa1 = Monitor(eid='BPME.2253.SA1')
bpme_2259_sa1 = Monitor(eid='BPME.2259.SA1')
bpme_2265_sa1 = Monitor(eid='BPME.2265.SA1')
bpme_2271_sa1 = Monitor(eid='BPME.2271.SA1')
bpme_2278_sa1 = Monitor(eid='BPME.2278.SA1')
bpme_2284_sa1 = Monitor(eid='BPME.2284.SA1')
bpme_2290_sa1 = Monitor(eid='BPME.2290.SA1')
bpme_2296_sa1 = Monitor(eid='BPME.2296.SA1')
bpme_2302_sa1 = Monitor(eid='BPME.2302.SA1')
bpme_2308_sa1 = Monitor(eid='BPME.2308.SA1')
bpme_2314_sa1 = Monitor(eid='BPME.2314.SA1')
bpme_2320_sa1 = Monitor(eid='BPME.2320.SA1')
bpme_2326_sa1 = Monitor(eid='BPME.2326.SA1')
bpme_2332_sa1 = Monitor(eid='BPME.2332.SA1')
bpme_2339_sa1 = Monitor(eid='BPME.2339.SA1')
bpme_2345_sa1 = Monitor(eid='BPME.2345.SA1')
bpme_2351_sa1 = Monitor(eid='BPME.2351.SA1')
bpme_2357_sa1 = Monitor(eid='BPME.2357.SA1')
bpme_2363_sa1 = Monitor(eid='BPME.2363.SA1')
bpme_2369_sa1 = Monitor(eid='BPME.2369.SA1')
bpme_2375_sa1 = Monitor(eid='BPME.2375.SA1')
bpme_2381_sa1 = Monitor(eid='BPME.2381.SA1')
bpme_2387_sa1 = Monitor(eid='BPME.2387.SA1')
bpme_2393_sa1 = Monitor(eid='BPME.2393.SA1')
bpme_2400_sa1 = Monitor(eid='BPME.2400.SA1')
bpme_2406_sa1 = Monitor(eid='BPME.2406.SA1')
bpme_2412_sa1 = Monitor(eid='BPME.2412.SA1')
bpme_2418_sa1 = Monitor(eid='BPME.2418.SA1')
bpme_2424_sa1 = Monitor(eid='BPME.2424.SA1')
bpme_2430_sa1 = Monitor(eid='BPME.2430.SA1')
bpme_2436_sa1 = Monitor(eid='BPME.2436.SA1')
bpme_2442_sa1 = Monitor(eid='BPME.2442.SA1')
bpme_2448_sa1 = Monitor(eid='BPME.2448.SA1')
bpme_2454_sa1 = Monitor(eid='BPME.2454.SA1')
bpme_2461_sa1 = Monitor(eid='BPME.2461.SA1')

# Markers
vcst40t10_2229_t2 = Marker(eid='VCST40T10.2229.T2')
match_2248_sa1 = Marker(eid='MATCH.2248.SA1')
ensec_2461_sa1 = Marker(eid='ENSEC.2461.SA1')

# Lattice
cell = (id_35371550_, qk_2027_tl, d_3, bpma_2041_tl, d_4, qf_2042_tl, d_5, cfx_2042_tl, d_6,
chx_2054_tl, d_7, chy_2054_tl, d_8, bpma_2054_tl, d_9, qk_2057_tl, id_23741906_, bpma_2071_t2, d_12,
qf_2072_t2, d_13, cfx_2072_t2, d_14, bd_2079_t2, d_15, bpma_2086_t2, d_12, qf_2087_t2, d_13,
cfy_2087_t2, d_18, bpma_2101_t2, d_12, qf_2102_t2, d_13, cfx_2102_t2, d_21, bpma_2116_t2, d_12,
qf_2117_t2, d_13, cfy_2117_t2, id_6624594_, bpma_2132_t2, d_12, qf_2132_t2, d_13, cfx_2133_t2, d_28,
qh_2139_t2, d_29, cmx_2140_t2, d_30, kl_2143_t2, d_31, bpma_2145_t2, d_12, qf_2145_t2, d_13,
cfy_2146_t2, d_34, bpma_2161_t2, d_12, qf_2162_t2, d_13, cfx_2162_t2, d_37, cmy_2162_t2, id_52103571_,
bpma_2176_t2, d_12, qf_2177_t2, d_13, cfy_2177_t2, d_42, cmx_2178_t2, d_43, qf_2184_t2, d_44,
bpma_2191_t2, d_12, qf_2192_t2, d_13, cfx_2192_t2, d_47, cmy_2192_t2, d_48, kl_2193_t2, id_38309706_,
bp_2197i_t2, cbp_2197i_t2, d_52, bp_2197ii_t2, d_52, bp_2198_t2, d_52, bp_2199i_t2, d_52, bp_2199ii_t2,
d_52, bp_2200_t2, d_52, bp_2201_t2, cbp_2201_t2, id_32688183_, bpma_2206_t2, d_12, qf_2207_t2, d_13,
cfy_2207_t2, d_63, kl_2208_t2, id_71170601_, bpma_2218_t2, d_12, qf_2218_t2, d_13, cfx_2219_t2, id_62721345_,
bpma_2223_t2, d_70, vcst40t10_2229_t2, d_71, bpme_2229_t2, d_72, qa_2229_t2, d_73, cny_2229_t2, d_74,
cex_2230_t2, d_75, u40s_2232_t2, d_76, cny_2234_t2, d_77, cex_2234_t2, d_78, bpme_2235_t2, d_79,
qa_2235_t2, id_85978457_, bpme_2241_sa1, d_79, qa_2241_sa1, id_85388498_, cax_2242_sa1, cay_2242_sa1, d_86, u40_2244_sa1,
d_87, cbx_2247_sa1, cby_2247_sa1, d_88, bpme_2247_sa1, d_79, qa_2247_sa1, id_83891648_, match_2248_sa1, id_54441912_,
cax_2248_sa1, cay_2248_sa1, d_86, u40_2250_sa1, d_87, cbx_2253_sa1, cby_2253_sa1, d_88, bpme_2253_sa1, d_79,
qa_2253_sa1, id_85388498_, cax_2254_sa1, cay_2254_sa1, d_86, u40_2256_sa1, d_87, cbx_2259_sa1, cby_2259_sa1, d_88,
bpme_2259_sa1, d_79, qa_2259_sa1, id_85388498_, cax_2260_sa1, cay_2260_sa1, d_86, u40_2262_sa1, d_87, cbx_2265_sa1,
cby_2265_sa1, d_88, bpme_2265_sa1, d_79, qa_2266_sa1, id_85388498_, cax_2267_sa1, cay_2267_sa1, d_86, u40_2269_sa1,
d_87, cbx_2271_sa1, cby_2271_sa1, d_88, bpme_2271_sa1, d_79, qa_2272_sa1, id_85388498_, cax_2273_sa1, cay_2273_sa1,
d_86, u40_2275_sa1, d_87, cbx_2277_sa1, cby_2277_sa1, d_88, bpme_2278_sa1, d_79, qa_2278_sa1, id_85388498_,
cax_2279_sa1, cay_2279_sa1, d_86, u40_2281_sa1, d_87, cbx_2283_sa1, cby_2283_sa1, d_88, bpme_2284_sa1, d_79,
qa_2284_sa1, id_85388498_, cax_2285_sa1, cay_2285_sa1, d_86, u40_2287_sa1, d_87, cbx_2289_sa1, cby_2289_sa1, d_88,
bpme_2290_sa1, d_79, qa_2290_sa1, id_85388498_, cax_2291_sa1, cay_2291_sa1, d_86, u40_2293_sa1, d_87, cbx_2296_sa1,
cby_2296_sa1, d_88, bpme_2296_sa1, d_79, qa_2296_sa1, id_85388498_, cax_2297_sa1, cay_2297_sa1, d_86, u40_2299_sa1,
d_87, cbx_2302_sa1, cby_2302_sa1, d_88, bpme_2302_sa1, d_79, qa_2302_sa1, id_85388498_, cax_2303_sa1, cay_2303_sa1,
d_86, u40_2305_sa1, d_87, cbx_2308_sa1, cby_2308_sa1, d_88, bpme_2308_sa1, d_79, qa_2308_sa1, id_85388498_,
cax_2309_sa1, cay_2309_sa1, d_86, u40_2311_sa1, d_87, cbx_2314_sa1, cby_2314_sa1, d_88, bpme_2314_sa1, d_79,
qa_2314_sa1, id_85388498_, cax_2315_sa1, cay_2315_sa1, d_86, u40_2317_sa1, d_87, cbx_2320_sa1, cby_2320_sa1, d_88,
bpme_2320_sa1, d_79, qa_2320_sa1, id_85388498_, cax_2321_sa1, cay_2321_sa1, d_86, u40_2323_sa1, d_87, cbx_2326_sa1,
cby_2326_sa1, d_88, bpme_2326_sa1, d_79, qa_2327_sa1, id_85388498_, cax_2328_sa1, cay_2328_sa1, d_86, u40_2330_sa1,
d_87, cbx_2332_sa1, cby_2332_sa1, d_88, bpme_2332_sa1, d_79, qa_2333_sa1, id_85388498_, cax_2334_sa1, cay_2334_sa1,
d_86, u40_2336_sa1, d_87, cbx_2338_sa1, cby_2338_sa1, d_88, bpme_2339_sa1, d_79, qa_2339_sa1, id_85388498_,
cax_2340_sa1, cay_2340_sa1, d_86, u40_2342_sa1, d_87, cbx_2344_sa1, cby_2344_sa1, d_88, bpme_2345_sa1, d_79,
qa_2345_sa1, id_85388498_, cax_2346_sa1, cay_2346_sa1, d_86, u40_2348_sa1, d_87, cbx_2350_sa1, cby_2350_sa1, d_88,
bpme_2351_sa1, d_79, qa_2351_sa1, id_85388498_, cax_2352_sa1, cay_2352_sa1, d_86, u40_2354_sa1, d_87, cbx_2357_sa1,
cby_2357_sa1, d_88, bpme_2357_sa1, d_79, qa_2357_sa1, id_85388498_, cax_2358_sa1, cay_2358_sa1, d_86, u40_2360_sa1,
d_87, cbx_2363_sa1, cby_2363_sa1, d_88, bpme_2363_sa1, d_79, qa_2363_sa1, id_85388498_, cax_2364_sa1, cay_2364_sa1,
d_86, u40_2366_sa1, d_87, cbx_2369_sa1, cby_2369_sa1, d_88, bpme_2369_sa1, d_79, qa_2369_sa1, id_85388498_,
cax_2370_sa1, cay_2370_sa1, d_86, u40_2372_sa1, d_87, cbx_2375_sa1, cby_2375_sa1, d_88, bpme_2375_sa1, d_79,
qa_2375_sa1, id_85388498_, cax_2376_sa1, cay_2376_sa1, d_86, u40_2378_sa1, d_87, cbx_2381_sa1, cby_2381_sa1, d_88,
bpme_2381_sa1, d_79, qa_2381_sa1, id_85388498_, cax_2382_sa1, cay_2382_sa1, d_86, u40_2384_sa1, d_87, cbx_2387_sa1,
cby_2387_sa1, d_88, bpme_2387_sa1, d_79, qa_2388_sa1, id_85388498_, cax_2389_sa1, cay_2389_sa1, d_86, u40_2391_sa1,
d_87, cbx_2393_sa1, cby_2393_sa1, d_88, bpme_2393_sa1, d_79, qa_2394_sa1, id_85388498_, cax_2395_sa1, cay_2395_sa1,
d_86, u40_2397_sa1, d_87, cbx_2399_sa1, cby_2399_sa1, d_88, bpme_2400_sa1, d_79, qa_2400_sa1, id_85388498_,
cax_2401_sa1, cay_2401_sa1, d_86, u40_2403_sa1, d_87, cbx_2405_sa1, cby_2405_sa1, d_88, bpme_2406_sa1, d_79,
qa_2406_sa1, id_85388498_, cax_2407_sa1, cay_2407_sa1, d_86, u40_2409_sa1, d_87, cbx_2411_sa1, cby_2411_sa1, d_88,
bpme_2412_sa1, d_79, qa_2412_sa1, id_85388498_, cax_2413_sa1, cay_2413_sa1, d_86, u40_2415_sa1, d_87, cbx_2418_sa1,
cby_2418_sa1, d_88, bpme_2418_sa1, d_79, qa_2418_sa1, id_85388498_, cax_2419_sa1, cay_2419_sa1, d_86, u40_2421_sa1,
d_87, cbx_2424_sa1, cby_2424_sa1, d_88, bpme_2424_sa1, d_79, qa_2424_sa1, id_85388498_, cax_2425_sa1, cay_2425_sa1,
d_86, u40_2427_sa1, d_87, cbx_2430_sa1, cby_2430_sa1, d_88, bpme_2430_sa1, d_79, qa_2430_sa1, id_17219612_,
bs_2431_sa1, id_53856480_, bs_2432_sa1, cbs_2432_sa1, id_85095730_, bs_2434_sa1, cbs_2434_sa1, id_53856480_, bs_2435_sa1, cbs_2435_sa1,
d_276, bpme_2436_sa1, d_79, qa_2436_sa1, id_85388498_, cax_2437_sa1, cay_2437_sa1, d_86, u40_2439_sa1, d_87,
cbx_2442_sa1, cby_2442_sa1, d_88, bpme_2442_sa1, d_79, qa_2442_sa1, id_85388498_, cax_2443_sa1, cay_2443_sa1, d_86,
u40_2445_sa1, d_87, cbx_2448_sa1, cby_2448_sa1, d_88, bpme_2448_sa1, d_79, qa_2449_sa1, id_85388498_, cax_2450_sa1,
cay_2450_sa1, d_86, u40_2452_sa1, d_87, cbx_2454_sa1, cby_2454_sa1, d_88, bpme_2454_sa1, d_79, qa_2455_sa1,
id_85388498_, cax_2456_sa1, cay_2456_sa1, d_86, u40_2458_sa1, d_87, cbx_2460_sa1, cby_2460_sa1, d_88, bpme_2461_sa1,
d_79, qa_2461_sa1, id_42404360_, ensec_2461_sa1)

# power supplies 

#  
qk_2027_tl.ps_id = 'QK.2.TL'
qf_2042_tl.ps_id = 'QF.1.TL'
qk_2057_tl.ps_id = 'QK.2.TL'
qf_2072_t2.ps_id = 'QF.1.T2'
qf_2087_t2.ps_id = 'QF.1.T2'
qf_2102_t2.ps_id = 'QF.1.T2'
qf_2117_t2.ps_id = 'QF.7.T2'
qf_2132_t2.ps_id = 'QF.8.T2'
qh_2139_t2.ps_id = 'QH.9.T2'
qf_2145_t2.ps_id = 'QF.10.T2'
qf_2162_t2.ps_id = 'QF.11.T2'
qf_2177_t2.ps_id = 'QF.3.T2'
qf_2184_t2.ps_id = 'QF.12.T2'
qf_2192_t2.ps_id = 'QF.4.T2'
qf_2207_t2.ps_id = 'QF.5.T2'
qf_2218_t2.ps_id = 'QF.6.T2'
qa_2229_t2.ps_id = 'QA.1.T2'
qa_2235_t2.ps_id = 'QA.2.T2'
qa_2241_sa1.ps_id = 'QA.1.SA1'
qa_2247_sa1.ps_id = 'QA.2.SA1'
qa_2253_sa1.ps_id = 'QA.1.SA1'
qa_2259_sa1.ps_id = 'QA.2.SA1'
qa_2266_sa1.ps_id = 'QA.1.SA1'
qa_2272_sa1.ps_id = 'QA.2.SA1'
qa_2278_sa1.ps_id = 'QA.1.SA1'
qa_2284_sa1.ps_id = 'QA.2.SA1'
qa_2290_sa1.ps_id = 'QA.1.SA1'
qa_2296_sa1.ps_id = 'QA.2.SA1'
qa_2302_sa1.ps_id = 'QA.1.SA1'
qa_2308_sa1.ps_id = 'QA.2.SA1'
qa_2314_sa1.ps_id = 'QA.1.SA1'
qa_2320_sa1.ps_id = 'QA.2.SA1'
qa_2327_sa1.ps_id = 'QA.1.SA1'
qa_2333_sa1.ps_id = 'QA.2.SA1'
qa_2339_sa1.ps_id = 'QA.1.SA1'
qa_2345_sa1.ps_id = 'QA.2.SA1'
qa_2351_sa1.ps_id = 'QA.1.SA1'
qa_2357_sa1.ps_id = 'QA.2.SA1'
qa_2363_sa1.ps_id = 'QA.1.SA1'
qa_2369_sa1.ps_id = 'QA.2.SA1'
qa_2375_sa1.ps_id = 'QA.1.SA1'
qa_2381_sa1.ps_id = 'QA.2.SA1'
qa_2388_sa1.ps_id = 'QA.1.SA1'
qa_2394_sa1.ps_id = 'QA.2.SA1'
qa_2400_sa1.ps_id = 'QA.1.SA1'
qa_2406_sa1.ps_id = 'QA.2.SA1'
qa_2412_sa1.ps_id = 'QA.1.SA1'
qa_2418_sa1.ps_id = 'QA.2.SA1'
qa_2424_sa1.ps_id = 'QA.1.SA1'
qa_2430_sa1.ps_id = 'QA.2.SA1'
qa_2436_sa1.ps_id = 'QA.1.SA1'
qa_2442_sa1.ps_id = 'QA.2.SA1'
qa_2449_sa1.ps_id = 'QA.1.SA1'
qa_2455_sa1.ps_id = 'QA.2.SA1'
qa_2461_sa1.ps_id = 'QA.1.SA1'

#  

#  

#  

#  
bd_2079_t2.ps_id = 'BD.10.T2'
kl_2143_t2.ps_id = 'KL.2143.T2'
kl_2193_t2.ps_id = 'KL.2193.T2'
bp_2197i_t2.ps_id = 'BP.1.T2'
bp_2197ii_t2.ps_id = 'BP.2.T2'
bp_2198_t2.ps_id = 'BP.2.T2'
bp_2199i_t2.ps_id = 'BP.2.T2'
bp_2199ii_t2.ps_id = 'BP.2.T2'
bp_2200_t2.ps_id = 'BP.2.T2'
bp_2201_t2.ps_id = 'BP.1.T2'
kl_2208_t2.ps_id = 'KL.2208.T2'
bs_2431_sa1.ps_id = 'BS.2.SA1'
bs_2432_sa1.ps_id = 'BS.2.SA1'
bs_2434_sa1.ps_id = 'BS.2.SA1'
bs_2435_sa1.ps_id = 'BS.2.SA1'
