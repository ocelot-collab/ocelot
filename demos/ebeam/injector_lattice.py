from ocelot import * 
tws = Twiss()
tws.E = 0.005
tws.beta_x  = 53.35971898
tws.beta_y  = 58.31307349
tws.alpha_x = 17.3149486
tws.alpha_y = 19.09296961
tws.s = 23.2
#tws = Twiss()
#tws.E = 0.005000000
#tws.beta_x  = 55.7887190242
#tws.beta_y  = 55.7887190242
#tws.alpha_x = 18.185436973
#tws.alpha_y = 18.185436973

# tws.beta_x =62.59953597
# tws.beta_y =99.94695119
# tws.alpha_x =  20.45567528
# tws.alpha_y =32.18973694


# Drifts
d_1 = Drift(l=0.276, eid='D_1')
d_2 = Drift(l=0.216, eid='D_2')
d_3 = Drift(l=0.311, eid='D_3')
d_4 = Drift(l=0.047, eid='D_4')
d_5 = Drift(l=0.788, eid='D_5')
d_6 = Drift(l=0.313, eid='D_6')
d_7 = Drift(l=0.421, eid='D_7')
d_8 = Drift(l=0.8992, eid='D_8')
d_9 = Drift(l=0.3459, eid='D_9')
d_16 = Drift(l=0.2043, eid='D_16')
d_17 = Drift(l=0.0432, eid='D_17')
d_19 = Drift(l=0.085, eid='D_19')
d_20 = Drift(l=0.679, eid='D_20')
d_21 = Drift(l=0.1282, eid='D_21')
d_23 = Drift(l=0.202, eid='D_23')
d_24 = Drift(l=0.262, eid='D_24')
d_31 = Drift(l=2.4914, eid='D_31')
d_32 = Drift(l=0.33305, eid='D_32')
d_33 = Drift(l=0.20765, eid='D_33')
d_34 = Drift(l=0.43315, eid='D_34')
d_35 = Drift(l=0.132315, eid='D_35')
d_36 = Drift(l=0.100827, eid='D_36')
d_37 = Drift(l=0.284415, eid='D_37')
d_38 = Drift(l=0.17575, eid='D_38')
d_39 = Drift(l=0.194, eid='D_39')
d_41 = Drift(l=0.120165, eid='D_41')
d_42 = Drift(l=0.100828, eid='D_42')
d_44 = Drift(l=0.38964, eid='D_44')
d_45 = Drift(l=0.303, eid='D_45')
d_46 = Drift(l=0.1015, eid='D_46')
d_47 = Drift(l=0.27965, eid='D_47')
d_48 = Drift(l=0.05115, eid='D_48')
d_50 = Drift(l=0.7143, eid='D_50')
d_51 = Drift(l=0.75615, eid='D_51')
d_52 = Drift(l=0.175, eid='D_52')
d_53 = Drift(l=0.15, eid='D_53')
d_54 = Drift(l=0.13115, eid='D_54')
d_60 = Drift(l=0.275, eid='D_60')
d_61 = Drift(l=0.18115, eid='D_61')
d_62 = Drift(l=0.20115, eid='D_62')
d_63 = Drift(l=0.38, eid='D_63')
d_64 = Drift(l=0.28115, eid='D_64')
d_65 = Drift(l=0.26285, eid='D_65')
d_66 = Drift(l=0.3643, eid='D_66')
d_68 = Drift(l=0.23515, eid='D_68')
d_69 = Drift(l=0.53115, eid='D_69')

# Quadrupoles
qln_23_i1 = Quadrupole(l=0.05, eid='QLN.23.I1')
qls_23_i1 = Quadrupole(l=0.05, tilt=0.785398163, eid='QLS.23.I1')
q_37_i1 = Quadrupole(l=0.2136, k1=-1.448361582, eid='Q.37.I1')
q_38_i1 = Quadrupole(l=0.2136, k1=1.4535982140, eid='Q.38.I1')
qi_46_i1 = Quadrupole(l=0.2377, k1=-0.2247807673, eid='QI.46.I1')
qi_47_i1 = Quadrupole(l=0.2377, k1=0.6927621982, eid='QI.47.I1')
qi_50_i1 = Quadrupole(l=0.2377, k1=-0.86466232944, eid='QI.50.I1')
qi_52_i1 = Quadrupole(l=0.2377, k1=-0.3522076137989062, eid='QI.52.I1')
qi_53_i1 = Quadrupole(l=0.2377, k1=2.1047941859486747, eid='QI.53.I1')
qi_54_i1 = Quadrupole(l=0.2377, k1=0.7943661064366849, eid='QI.54.I1')
qi_55_i1 = Quadrupole(l=0.2377, k1=-2.6383360778291967, eid='QI.55.I1')
qi_57_i1 = Quadrupole(l=0.2377, k1=3.278062023979807, eid='QI.57.I1')
qi_59_i1 = Quadrupole(l=0.2377, k1=-2.6383360778291967, eid='QI.59.I1')
qi_60_i1 = Quadrupole(l=0.2377, k1=1.9778194619267986, eid='QI.60.I1')
qi_61_i1 = Quadrupole(l=0.2377, k1=0.8708619869583508, eid='QI.61.I1')

# SBends
bl_48i_i1 = SBend(l=0.2, angle=-0.099484, e2=-0.099484, eid='BL.48I.I1')
bl_48ii_i1 = SBend(l=0.2, angle=0.099484, e1=0.099484, eid='BL.48II.I1')
bl_50i_i1 = SBend(l=0.2, angle=0.099484, e2=0.099484, eid='BL.50I.I1')
bl_50ii_i1 = SBend(l=0.2, angle=-0.099484, e1=-0.099484, eid='BL.50II.I1')

# Hcors
ckx_23_i1 = Hcor(l=0.025, eid='CKX.23.I1')
ckx_24_i1 = Hcor(l=0.025, eid='CKX.24.I1')
ckx_25_i1 = Hcor(l=0.025, eid='CKX.25.I1')
cx_37_i1 = Hcor(eid='CX.37.I1')
cx_39_i1 = Hcor(eid='CX.39.I1')
cix_51_i1 = Hcor(l=0.1, eid='CIX.51.I1')
cix_57_i1 = Hcor(l=0.1, eid='CIX.57.I1')

# Vcors
cky_23_i1 = Vcor(l=0.025, eid='CKY.23.I1')
cky_24_i1 = Vcor(l=0.025, eid='CKY.24.I1')
cky_25_i1 = Vcor(l=0.025, eid='CKY.25.I1')
cy_37_i1 = Vcor(eid='CY.37.I1')
cy_39_i1 = Vcor(eid='CY.39.I1')
ciy_51_i1 = Vcor(l=0.1, eid='CIY.51.I1')
ciy_55_i1 = Vcor(l=0.1, eid='CIY.55.I1')
ciy_58_i1 = Vcor(l=0.1, eid='CIY.58.I1')

# Undulators
Kx= 1.294
Kx=1.315
undu_49_i1 = Undulator(lperiod=0.074, nperiods=10.0, Kx=Kx, eid='UNDU.49.I1')

# Cavitys
dip_kick = 0
c_a1_1_1_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, vx_up=(-5.6813e-05+1.0751e-05j)*dip_kick, vy_up=(-4.1091e-05+5.739e-07j)*dip_kick, vxx_up=(0.00099943-0.00081401j), vxy_up=(0.0034065-0.0004146j), vx_down=(-2.4014e-05+1.2492e-05j)*dip_kick, vy_down=(3.6481e-05+7.9888e-06j)*dip_kick, vxx_down=(-0.004057-0.0001369j), vxy_down=(0.0029243-1.2891e-05j), eid='C.A1.1.1.I1')
c_a1_1_2_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, vx_up=(-5.6813e-05+1.0751e-05j)*dip_kick, vy_up=(-4.1091e-05+5.739e-07j)*dip_kick, vxx_up=(0.00099943-0.00081401j), vxy_up=(0.0034065-0.0004146j), vx_down=(-2.4014e-05+1.2492e-05j)*dip_kick, vy_down=(3.6481e-05+7.9888e-06j)*dip_kick, vxx_down=(-0.004057-0.0001369j), vxy_down=(0.0029243-1.2891e-05j), eid='C.A1.1.2.I1')
c_a1_1_3_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, vx_up=(-5.6813e-05+1.0751e-05j)*dip_kick, vy_up=(-4.1091e-05+5.739e-07j)*dip_kick, vxx_up=(0.00099943-0.00081401j), vxy_up=(0.0034065-0.0004146j), vx_down=(-2.4014e-05+1.2492e-05j)*dip_kick, vy_down=(3.6481e-05+7.9888e-06j)*dip_kick, vxx_down=(-0.004057-0.0001369j), vxy_down=(0.0029243-1.2891e-05j), eid='C.A1.1.3.I1')
c_a1_1_4_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, vx_up=(-5.6813e-05+1.0751e-05j)*dip_kick, vy_up=(-4.1091e-05+5.739e-07j)*dip_kick, vxx_up=(0.00099943-0.00081401j), vxy_up=(0.0034065-0.0004146j), vx_down=(-2.4014e-05+1.2492e-05j)*dip_kick, vy_down=(3.6481e-05+7.9888e-06j)*dip_kick, vxx_down=(-0.004057-0.0001369j), vxy_down=(0.0029243-1.2891e-05j), eid='C.A1.1.4.I1')
c_a1_1_5_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, vx_up=(-5.6813e-05+1.0751e-05j)*dip_kick, vy_up=(-4.1091e-05+5.739e-07j)*dip_kick, vxx_up=(0.00099943-0.00081401j), vxy_up=(0.0034065-0.0004146j), vx_down=(-2.4014e-05+1.2492e-05j)*dip_kick, vy_down=(3.6481e-05+7.9888e-06j)*dip_kick, vxx_down=(-0.004057-0.0001369j), vxy_down=(0.0029243-1.2891e-05j), eid='C.A1.1.5.I1')
c_a1_1_6_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, vx_up=(-5.6813e-05+1.0751e-05j)*dip_kick, vy_up=(-4.1091e-05+5.739e-07j)*dip_kick, vxx_up=(0.00099943-0.00081401j), vxy_up=(0.0034065-0.0004146j), vx_down=(-2.4014e-05+1.2492e-05j)*dip_kick, vy_down=(3.6481e-05+7.9888e-06j)*dip_kick, vxx_down=(-0.004057-0.0001369j), vxy_down=(0.0029243-1.2891e-05j), eid='C.A1.1.6.I1')
c_a1_1_7_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, vx_up=(-5.6813e-05+1.0751e-05j)*dip_kick, vy_up=(-4.1091e-05+5.739e-07j)*dip_kick, vxx_up=(0.00099943-0.00081401j), vxy_up=(0.0034065-0.0004146j), vx_down=(-2.4014e-05+1.2492e-05j)*dip_kick, vy_down=(3.6481e-05+7.9888e-06j)*dip_kick, vxx_down=(-0.004057-0.0001369j), vxy_down=(0.0029243-1.2891e-05j), eid='C.A1.1.7.I1')
c_a1_1_8_i1 = Cavity(l=1.0377, v=0.018125, freq=1300000000.0, vx_up=(-5.6813e-05+1.0751e-05j)*dip_kick, vy_up=(-4.1091e-05+5.739e-07j)*dip_kick, vxx_up=(0.00099943-0.00081401j), vxy_up=(0.0034065-0.0004146j), vx_down=(-2.4014e-05+1.2492e-05j)*dip_kick, vy_down=(3.6481e-05+7.9888e-06j)*dip_kick, vxx_down=(-0.004057-0.0001369j), vxy_down=(0.0029243-1.2891e-05j), eid='C.A1.1.8.I1')

c3_ah1_1_1_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, vx_up=dip_kick*(-0.00057076-1.3166e-05j), vy_up=dip_kick*(-3.5079e-05+0.00012636j), vxx_up=(-0.026045-0.042918j), vxy_up=(0.0055553-0.023455j), vx_down=dip_kick*(-8.8766e-05-0.00024852j), vy_down=dip_kick*(2.9889e-05+0.00014486j), vxx_down=(-0.0050593-0.013491j), vxy_down=(0.0051488+0.024771j), eid='C3.AH1.1.1.I1')
c3_ah1_1_2_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, vx_up=dip_kick*(0.00057076+1.3166e-05j), vy_up= dip_kick*(3.5079e-05-0.00012636j), vxx_up=(-0.026045-0.042918j), vxy_up=(0.0055553-0.023455j), vx_down= dip_kick*(8.8766e-05+0.00024852j), vy_down= dip_kick*(-2.9889e-05-0.00014486j), vxx_down=(-0.0050593-0.013491j), vxy_down=(0.0051488+0.024771j), eid='C3.AH1.1.2.I1')
c3_ah1_1_3_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, vx_up=dip_kick*(-0.00057076-1.3166e-05j), vy_up=dip_kick*(-3.5079e-05+0.00012636j), vxx_up=(-0.026045-0.042918j), vxy_up=(0.0055553-0.023455j), vx_down=dip_kick*(-8.8766e-05-0.00024852j), vy_down=dip_kick*(2.9889e-05+0.00014486j), vxx_down=(-0.0050593-0.013491j), vxy_down=(0.0051488+0.024771j), eid='C3.AH1.1.3.I1')
c3_ah1_1_4_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, vx_up=dip_kick*(0.00057076+1.3166e-05j), vy_up= dip_kick*(3.5079e-05-0.00012636j), vxx_up=(-0.026045-0.042918j), vxy_up=(0.0055553-0.023455j), vx_down= dip_kick*(8.8766e-05+0.00024852j), vy_down= dip_kick*(-2.9889e-05-0.00014486j), vxx_down=(-0.0050593-0.013491j), vxy_down=(0.0051488+0.024771j), eid='C3.AH1.1.4.I1')
c3_ah1_1_5_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, vx_up=dip_kick*(-0.00057076-1.3166e-05j), vy_up=dip_kick*(-3.5079e-05+0.00012636j), vxx_up=(-0.026045-0.042918j), vxy_up=(0.0055553-0.023455j), vx_down=dip_kick*(-8.8766e-05-0.00024852j), vy_down=dip_kick*(2.9889e-05+0.00014486j), vxx_down=(-0.0050593-0.013491j), vxy_down=(0.0051488+0.024771j), eid='C3.AH1.1.5.I1')
c3_ah1_1_6_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, vx_up=dip_kick*(0.00057076+1.3166e-05j), vy_up= dip_kick*(3.5079e-05-0.00012636j), vxx_up=(-0.026045-0.042918j), vxy_up=(0.0055553-0.023455j), vx_down= dip_kick*(8.8766e-05+0.00024852j), vy_down= dip_kick*(-2.9889e-05-0.00014486j), vxx_down=(-0.0050593-0.013491j), vxy_down=(0.0051488+0.024771j), eid='C3.AH1.1.6.I1')
c3_ah1_1_7_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, vx_up=dip_kick*(-0.00057076-1.3166e-05j), vy_up=dip_kick*(-3.5079e-05+0.00012636j), vxx_up=(-0.026045-0.042918j), vxy_up=(0.0055553-0.023455j), vx_down=dip_kick*(-8.8766e-05-0.00024852j), vy_down=dip_kick*(2.9889e-05+0.00014486j), vxx_down=(-0.0050593-0.013491j), vxy_down=(0.0051488+0.024771j), eid='C3.AH1.1.7.I1')
c3_ah1_1_8_i1 = Cavity(l=0.346, v=0.0025, phi=180.0, freq=3900000000.0, vx_up=dip_kick*(0.00057076+1.3166e-05j), vy_up= dip_kick*(3.5079e-05-0.00012636j), vxx_up=(-0.026045-0.042918j), vxy_up=(0.0055553-0.023455j), vx_down= dip_kick*(8.8766e-05+0.00024852j), vy_down= dip_kick*(-2.9889e-05-0.00014486j), vxx_down=(-0.0050593-0.013491j), vxy_down=(0.0051488+0.024771j), eid='C3.AH1.1.8.I1')

# TDCavitys
tdsa_52_i1 = TDCavity(l=0.7, freq=2800000000.0, eid='TDSA.52.I1')

# Solenoids
solb_23_i1 = Solenoid(eid='SOLB.23.I1')

# Monitors
bpmg_24_i1 = Monitor(eid='BPMG.24.I1')
bpmg_25i_i1 = Monitor(eid='BPMG.25I.I1')
bpmc_38i_i1 = Monitor(eid='BPMC.38I.I1')
bpmr_38ii_i1 = Monitor(eid='BPMR.38II.I1')
bpmf_47_i1 = Monitor(eid='BPMF.47.I1')
bpmf_48_i1 = Monitor(eid='BPMF.48.I1')
bpmf_52_i1 = Monitor(eid='BPMF.52.I1')
bpma_55_i1 = Monitor(eid='BPMA.55.I1')
bpma_57_i1 = Monitor(eid='BPMA.57.I1')
bpma_59_i1 = Monitor(eid='BPMA.59.I1')
bpma_60_i1 = Monitor(eid='BPMA.60.I1')
bpmatest_61_i1 = Drift(eid='BPMATEST.61.I1')

# Markers
stsec_23_i1 = Marker(eid='STSEC.23.I1')
tora_25_i1 = Marker(eid='TORA.25.I1')
match_37_i1 = Marker(eid='MATCH.37.I1')
tora_46_i1 = Marker(eid='TORA.46.I1')
otrl_48_i1 = Marker(eid='OTRL.48.I1')
otrl_50_i1 = Marker(eid='OTRL.50.I1')
match_55_i1 = Marker(eid='MATCH.55.I1')
otrc_55_i1 = Marker(eid='OTRC.55.I1')
otrc_56_i1 = Marker(eid='OTRC.56.I1')
otrc_58_i1 = Marker(eid='OTRC.58.I1')
otrc_59_i1 = Marker(eid='OTRC.59.I1')
tora_60_i1 = Marker(eid='TORA.60.I1')
stsub_62_i1 = Marker(eid='STSUB.62.I1')
lh_start = Marker()
lh_stop = Marker()

d_8_1 = Drift(l=0.5776, eid='D_12')
start_sim = Marker(eid="START_SIM")
d_8_2 = Drift(l=d_8.l - d_8_1.l)
d_8_n = (d_8_1, start_sim, d_8_2)

a1_sim_stop = Marker()
a1_1_stop = Marker()

d_35_1 = Drift(l=0.08115, eid='D_46')
stlat_47_i1 = Marker()
d_35_2 = Drift(l=d_35.l - d_35_1.l)

d_35_n = (d_35_1, stlat_47_i1, d_35_2)
# Lattice 
cell = (stsec_23_i1, d_1, solb_23_i1, qln_23_i1, qls_23_i1, d_2, ckx_23_i1, cky_23_i1, d_3, 
ckx_24_i1, cky_24_i1, d_4, bpmg_24_i1, d_5, tora_25_i1, d_6, bpmg_25i_i1, d_7, ckx_25_i1, 
cky_25_i1, d_8_n, c_a1_1_1_i1, a1_1_stop, d_9, c_a1_1_2_i1, d_9, c_a1_1_3_i1, d_9, c_a1_1_4_i1, d_9,
c_a1_1_5_i1, d_9, c_a1_1_6_i1, d_9, c_a1_1_7_i1, d_9, c_a1_1_8_i1, a1_sim_stop, d_16, match_37_i1, d_17,
q_37_i1, d_17, cx_37_i1, cy_37_i1, d_19, bpmc_38i_i1, d_20, bpmr_38ii_i1, d_21, q_38_i1, 
d_17, cx_39_i1, cy_39_i1, d_23, c3_ah1_1_1_i1, d_24, c3_ah1_1_2_i1, d_24, c3_ah1_1_3_i1, d_24, 
c3_ah1_1_4_i1, d_24, c3_ah1_1_5_i1, d_24, c3_ah1_1_6_i1, d_24, c3_ah1_1_7_i1, d_24, c3_ah1_1_8_i1, d_31, 
tora_46_i1, d_32, qi_46_i1, d_33, bpmf_47_i1, d_34, qi_47_i1, d_35_n, bl_48i_i1, d_36,
bl_48ii_i1, d_37, bpmf_48_i1, d_38, otrl_48_i1, d_39, lh_start, undu_49_i1,lh_stop, d_39, otrl_50_i1, d_41,
bl_50i_i1, d_42, bl_50ii_i1, d_35, qi_50_i1, d_44, ciy_51_i1, d_45, cix_51_i1, d_46, 
bpmf_52_i1, d_47, qi_52_i1, d_48, tdsa_52_i1, d_48, qi_53_i1, d_50, qi_54_i1, d_51, 
match_55_i1, otrc_55_i1, d_52, ciy_55_i1, d_53, bpma_55_i1, d_54, qi_55_i1, d_51, otrc_56_i1, 
d_52, cix_57_i1, d_53, bpma_57_i1, d_54, qi_57_i1, d_51, otrc_58_i1, d_60, ciy_58_i1, 
d_61, qi_59_i1, d_62, bpma_59_i1, d_63, otrc_59_i1, d_64, qi_60_i1, d_65, tora_60_i1, 
d_66, bpma_60_i1, d_53, bpmatest_61_i1, d_68, qi_61_i1, d_69, stsub_62_i1)

# power supplies 

#  
qln_23_i1.ps_id = 'QLN.23.I1'
qls_23_i1.ps_id = 'QLS.23.I1'
q_37_i1.ps_id = 'Q.A1.1.I1'
q_38_i1.ps_id = 'Q.AH1.1.I1'
qi_46_i1.ps_id = 'QI.1.I1'
qi_47_i1.ps_id = 'QI.2.I1'
qi_50_i1.ps_id = 'QI.3.I1'
qi_52_i1.ps_id = 'QI.4.I1'
qi_53_i1.ps_id = 'QI.5.I1'
qi_54_i1.ps_id = 'QI.6.I1'
qi_55_i1.ps_id = 'QI.7.I1'
qi_57_i1.ps_id = 'QI.8.I1'
qi_59_i1.ps_id = 'QI.9.I1'
qi_60_i1.ps_id = 'QI.11.I1'
qi_61_i1.ps_id = 'QI.12.I1'

#  

#  

#  
c_a1_1_1_i1.ps_id = 'C.A1.I1'
c_a1_1_2_i1.ps_id = 'C.A1.I1'
c_a1_1_3_i1.ps_id = 'C.A1.I1'
c_a1_1_4_i1.ps_id = 'C.A1.I1'
c_a1_1_5_i1.ps_id = 'C.A1.I1'
c_a1_1_6_i1.ps_id = 'C.A1.I1'
c_a1_1_7_i1.ps_id = 'C.A1.I1'
c_a1_1_8_i1.ps_id = 'C.A1.I1'
c3_ah1_1_1_i1.ps_id = 'C3.AH1.I1'
c3_ah1_1_2_i1.ps_id = 'C3.AH1.I1'
c3_ah1_1_3_i1.ps_id = 'C3.AH1.I1'
c3_ah1_1_4_i1.ps_id = 'C3.AH1.I1'
c3_ah1_1_5_i1.ps_id = 'C3.AH1.I1'
c3_ah1_1_6_i1.ps_id = 'C3.AH1.I1'
c3_ah1_1_7_i1.ps_id = 'C3.AH1.I1'
c3_ah1_1_8_i1.ps_id = 'C3.AH1.I1'

#  
bl_48i_i1.ps_id = 'BL.1.I1'
bl_48ii_i1.ps_id = 'BL.1.I1'
bl_50i_i1.ps_id = 'BL.3.I1'
bl_50ii_i1.ps_id = 'BL.4.I1'
