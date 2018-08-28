from ocelot import *

tws = Twiss()
tws.beta_x  = 10.506745988156398
tws.beta_y  = 42.02704133328497
tws.alpha_x = 0.676085898798978
tws.alpha_y = -2.1449229692237783
tws.E       = 14

# drifts 
d0 = Drift(l=1.472401, eid='D0')
d1 = Drift(l=14.206372, eid='D1')
d2 = Drift(l=14.206351999999999, eid='D2')
d3 = Drift(l=14.206391, eid='D3')
d4 = Drift(l=14.467899999999998, eid='D4')
d5 = Drift(l=14.4679, eid='D5')
d6 = Drift(l=14.468, eid='D6')
d7 = Drift(l=15.1679, eid='D7')
d8 = Drift(l=12.5679, eid='D8')
d9 = Drift(l=15.6679, eid='D9')
d10 = Drift(l=11.2679, eid='D10')
d11 = Drift(l=10.31115, eid='D11')
d12 = Drift(l=2.87065, eid='D12')
d13 = Drift(l=2.9956500000000004, eid='D13')
d14 = Drift(l=5.9863, eid='D14')
d15 = Drift(l=0.43065000000000003, eid='D15')
d16 = Drift(l=0.55565, eid='D16')
d17 = Drift(l=0.04015, eid='D17')

# quadrupoles 
qk_2027_tl = Quadrupole(l=1.0552, k1=-0.0903596001, tilt=0.0, eid='QK.2027.TL')
qf_2042_tl = Quadrupole(l=0.5321, k1=0.1791908476, tilt=0.0, eid='QF.2042.TL')
qf_2072_t2 = Quadrupole(l=0.5321, k1=0.1791908476, tilt=0.0, eid='QF.2072.T2')
qf_2087_t2 = Quadrupole(l=0.5321, k1=-0.1791908476, tilt=0.0, eid='QF.2087.T2')
qf_2102_t2 = Quadrupole(l=0.5321, k1=0.1791908476, tilt=0.0, eid='QF.2102.T2')
qf_2117_t2 = Quadrupole(l=0.5321, k1=-0.1791908476, tilt=0.0, eid='QF.2117.T2')
qf_2132_t2 = Quadrupole(l=0.5321, k1=0.1791908476, tilt=0.0, eid='QF.2132.T2')
qf_2145_t2 = Quadrupole(l=0.5321, k1=-0.1791908476, tilt=0.0, eid='QF.2145.T2')
qf_2162_t2 = Quadrupole(l=0.5321, k1=0.1791908476, tilt=0.0, eid='QF.2162.T2')
qf_2177_t2 = Quadrupole(l=0.5321, k1=-0.1525567278, tilt=0.0, eid='QF.2177.T2')
qf_2192_t2 = Quadrupole(l=0.5321, k1=0.1521721569, tilt=0.0, eid='QF.2192.T2')
qf_2207_t2 = Quadrupole(l=0.5321, k1=-0.1259982434, tilt=0.0, eid='QF.2207.T2')
qf_2218_t2 = Quadrupole(l=0.5321, k1=0.1191993195, tilt=0.0, eid='QF.2218.T2')
qa_2229_t2 = Quadrupole(l=0.1137, k1=-0.5445788788, tilt=0.0, eid='QA.2229.T2')
qa_2235_t2 = Quadrupole(l=0.1137, k1=0.5501622372, tilt=0.0, eid='QA.2235.T2')
qa_2241_sa1 = Quadrupole(l=0.1137, k1=-0.5445788788, tilt=0.0, eid='QA.2241.SA1')
qa_2247_sa1 = Quadrupole(l=0.1137, k1=0.5501622372, tilt=0.0, eid='QA.2247.SA1')

qa_2247_sa1_h = Quadrupole(l=0.1137/2, k1=0.5501622372, tilt=0.0, eid='QA.2247.SA1')
fodo_match = Marker()
# bending magnets 

# correctors 

# markers 

# monitor 

# sextupoles 

# octupole 

# undulator
K = 3
u40s_2232_t2 = Undulator(lperiod=0.04, nperiods=3, Kx=K, Ky=0.0, eid='U40S.2232.T2')
u40_2250_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2250.SA1')
u40_2256_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2256.SA1')
u40_2262_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2262.SA1')
u40_2269_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2269.SA1')
u40_2275_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2275.SA1')
u40_2281_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2281.SA1')
u40_2287_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2287.SA1')
u40_2293_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2293.SA1')
u40_2299_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2299.SA1')
u40_2305_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2305.SA1')
u40_2311_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2311.SA1')
u40_2317_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2317.SA1')
u40_2323_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2323.SA1')
u40_2330_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2330.SA1')
u40_2336_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2336.SA1')
u40_2342_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2342.SA1')
u40_2348_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2348.SA1')
u40_2354_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2354.SA1')
u40_2360_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2360.SA1')
u40_2366_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2366.SA1')
u40_2372_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2372.SA1')
u40_2378_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2378.SA1')
u40_2384_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2384.SA1')
u40_2391_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2391.SA1')
u40_2397_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2397.SA1')
u40_2403_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2403.SA1')
u40_2409_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2409.SA1')
u40_2415_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2415.SA1')
u40_2421_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2421.SA1')
u40_2427_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2427.SA1')
u40_2433_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2433.SA1')
u40_2439_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2439.SA1')
u40_2445_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2445.SA1')
u40_2452_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2452.SA1')
u40_2458_sa1 = Undulator(lperiod=0.04, nperiods=125, Kx=K, Ky=0.0, eid='U40.2458.SA1')

# cavity 

# tdcavity 

# UnknowElement 

# Matrices 

# Solenoids 
start = Marker()
# lattice 
cell = (start, d0, qk_2027_tl, d1, qf_2042_tl, d2, qk_2027_tl, d3,
qf_2072_t2, d4, qf_2087_t2, d5, qf_2102_t2, d6, qf_2117_t2, d7, 
qf_2132_t2, d8, qf_2145_t2, d9, qf_2162_t2, d5, qf_2177_t2, d5, 
qf_2192_t2, d5, qf_2207_t2, d10, qf_2218_t2, d11, qa_2229_t2, d12, 
u40s_2232_t2, d13, qa_2235_t2, d14, qa_2241_sa1, d14,qa_2247_sa1_h, fodo_match, qa_2247_sa1_h, #qa_2247_sa1,
        d15,
u40_2250_sa1, d16, qa_2241_sa1, d15, u40_2256_sa1, d16, qa_2247_sa1, d15, 
u40_2262_sa1, d16, qa_2241_sa1, d15, u40_2269_sa1, d16, qa_2247_sa1, d15, 
u40_2275_sa1, d16, qa_2241_sa1, d15, u40_2281_sa1, d16, qa_2247_sa1, d15, 
u40_2287_sa1, d16, qa_2241_sa1, d15, u40_2293_sa1, d16, qa_2247_sa1, d15, 
u40_2299_sa1, d16, qa_2241_sa1, d15, u40_2305_sa1, d16, qa_2247_sa1, d15, 
u40_2311_sa1, d16, qa_2241_sa1, d15, u40_2317_sa1, d16, qa_2247_sa1, d15, 
u40_2323_sa1, d16, qa_2241_sa1, d15, u40_2330_sa1, d16, qa_2247_sa1, d15, 
u40_2336_sa1, d16, qa_2241_sa1, d15, u40_2342_sa1, d16, qa_2247_sa1, d15, 
u40_2348_sa1, d16, qa_2241_sa1, d15, u40_2354_sa1, d16, qa_2247_sa1, d15, 
u40_2360_sa1, d16, qa_2241_sa1, d15, u40_2366_sa1, d16, qa_2247_sa1, d15, 
u40_2372_sa1, d16, qa_2241_sa1, d15, u40_2378_sa1, d16, qa_2247_sa1, d15, 
u40_2384_sa1, d16, qa_2241_sa1, d15, u40_2391_sa1, d16, qa_2247_sa1, d15, 
u40_2397_sa1, d16, qa_2241_sa1, d15, u40_2403_sa1, d16, qa_2247_sa1, d15, 
u40_2409_sa1, d16, qa_2241_sa1, d15, u40_2415_sa1, d16, qa_2247_sa1, d15, 
u40_2421_sa1, d16, qa_2241_sa1, d15, u40_2427_sa1, d16, qa_2247_sa1, d15, 
u40_2433_sa1, d16, qa_2241_sa1, d15, u40_2439_sa1, d16, qa_2247_sa1, d15, 
u40_2445_sa1, d16, qa_2241_sa1, d15, u40_2452_sa1, d16, qa_2247_sa1, d15, 
u40_2458_sa1, d16, qa_2241_sa1, d17)
# power supplies 

#  
qk_2027_tl.ps_id = 'QK.2.TL'
qf_2042_tl.ps_id = 'QF.1.TL'
qf_2072_t2.ps_id = 'QF.1.T2'
qf_2087_t2.ps_id = 'QF.1.T2'
qf_2102_t2.ps_id = 'QF.1.T2'
qf_2117_t2.ps_id = 'QF.1.T2'
qf_2132_t2.ps_id = 'QF.1.T2'
qf_2145_t2.ps_id = 'QF.2.T2'
qf_2162_t2.ps_id = 'QF.2.T2'
qf_2177_t2.ps_id = 'QF.3.T2'
qf_2192_t2.ps_id = 'QF.4.T2'
qf_2207_t2.ps_id = 'QF.5.T2'
qf_2218_t2.ps_id = 'QF.6.T2'
qa_2229_t2.ps_id = 'QA.1.T2'
qa_2235_t2.ps_id = 'QA.2.T2'
qa_2241_sa1.ps_id = 'QA.1.SA1'
qa_2247_sa1.ps_id = 'QA.2.SA1'

#  

#  

#  

#  
