#import sys
#sys.path.append("/Volumes/Promise RAID/UserFolders/zagor_xxl/ocelot")

import matplotlib
matplotlib.use('TkAgg')
from accelerator.s2e_sections.sections import *
from ocelot.utils.section_track import *
from ocelot.gui.accelerator import *
import time
from ext4s2e import *
from ocelot.adaptors.csrtrack2ocelot import *

#data_dir = "N:/4all/xxl/zagor/mpy_xxl"
data_dir = "/Users/zagor/dataxxl/ocelot"
#data_dir = "/Volumes/Promise RAID/UserFolders/zagor_xxl/dataxxl/ocelot"


FirstOrderOptics=False
all_sections = [A1, AH1, LH, DL,  BC0, L1, BC1, L2, BC2, L3, CL1, CL2, CL3, STN10]#, SASE1, T4, SASE3, T4D]

######################### design optics
tws0 = Twiss()
tws0.E = 0.005
tws0.beta_x  = 0.286527307369
tws0.beta_y  = 0.286527307369
tws0.alpha_x = -0.838833736086
tws0.alpha_y = -0.838833736086
section_lat = SectionLattice(sequence=all_sections, tws0=tws0, data_dir=data_dir)

sections = [DL]

if FirstOrderOptics:
    for s in sections:
        #section_lat.dict_sections[s].method.global_method = SecondTM
        section_lat.dict_sections[s].method.global_method = TransferMap
        section_lat.dict_sections[s].lattice.update_transfer_maps()

start = time.time()

config = {
        DL: {"match": False, "SC": False,"CSR": False, "wake": False},
    }

####################### TRACKING ###############################

#p_array = load_particle_array(data_dir + "/particles/section_LH.npz")
emit0=3.930758660269919e-09;
#p_array=generate_parray(sigma_x=np.sqrt(emit0), sigma_px=np.sqrt(emit0), sigma_y=None, sigma_py=None,
#                    sigma_tau=0, sigma_p=0, tau_p_cor=0, charge=5e-9, nparticles=6000, energy=0.13,
#                    tue_trunc=None)
#p_array.p()[2000:4000] = p_array.p()[2000:4000]-0.05
#p_array.p()[4000:6000] = p_array.p()[4000:6000]+0.05
#s00 = copy(p_array.s)
#p_array=csrtrackBeam2particleArray('/Users/zagor/dataxxl/works/Dogleg/N2c_dogleg/in/in.fmt1', orient="V")

p_array=generate_parray(sigma_x=np.sqrt(emit0), sigma_px=np.sqrt(emit0), sigma_y=None, sigma_py=None,
                    sigma_tau=0, sigma_p=0, tau_p_cor=0, charge=5e-9, nparticles=3, energy=0.13,
                    tue_trunc=None)
p_array.x()[:]=0.0
p_array.y()[:]=0.0
p_array.px()[:]=0.0
p_array.py()[:]=0.0
p_array.tau()[:]=0.0
p_array.p()[0:3]=-0.05


s00= 50.198195000000126
print('s00=',s00)

p_array = section_lat.track_sections(sections=sections, p_array=p_array, config=config, force_ext_p_array=True)

method=MethodTM()
method.global_method = TransferMap
lat = MagneticLattice(section_lat.elem_seq[3])
R = lattice_transfer_map(lat, energy = 130*1e+3)
r36_DL=R[2,5]
r16_DL=R[0,5]
print("r36_DL=",r36_DL,"   r16_DL=",r16_DL)

## collect tws for all sections
seq_global = []
tws_track_global = []
L = 0
for s in sections:
    sec = section_lat.dict_sections[s]
    seq_global.append(sec.lattice.sequence)
    for tws in sec.tws_track:
        tws.s += L
    tws_track_global = np.append(tws_track_global, sec.tws_track)
    L += sec.lattice.totalLen

## independent analysis
S = [tw.s+3.2 for tw in section_lat.tws]
BetaX = [tw.beta_x for tw in section_lat.tws]
BetaY = [tw.beta_y for tw in section_lat.tws]
AlphaX = [tw.alpha_x for tw in section_lat.tws]
AlphaY = [tw.alpha_y for tw in section_lat.tws]
GammaX = [tw.gamma_x for tw in section_lat.tws]
GammaY = [tw.gamma_y for tw in section_lat.tws]
E = [tw.E for tw in section_lat.tws]


S_tr = [tw.s+s00 for tw in tws_track_global]
BetaX_tr = [tw.beta_x for tw in tws_track_global]
BetaY_tr = [tw.beta_y for tw in tws_track_global]
AlphaX_tr = [tw.alpha_x for tw in tws_track_global]
AlphaY_tr = [tw.alpha_y for tw in tws_track_global]
GammaX_tr = [tw.gamma_x for tw in tws_track_global]
GammaY_tr = [tw.gamma_y for tw in tws_track_global]
X_tr = [tw.x for tw in tws_track_global]
Y_tr = [tw.y for tw in tws_track_global]
XX_tr = [tw.xx for tw in tws_track_global]
YY_tr = [tw.yy for tw in tws_track_global]
XY_tr = [tw.xy for tw in tws_track_global]



E_tr=[tw.E for tw in tws_track_global]


#Plots
s0 = 45; s1 = 75
plt.figure()
plt.plot(S_tr,X_tr,S_tr,Y_tr)
plt.axis([s0, s1,-0.0001,0.0001])
plt.figure()
plt.subplot(321)
plt.plot(S,BetaX,label='design');plt.plot(S_tr,BetaX_tr,label='beam')
plt.legend(); plt.ylabel(r"$\beta_x$[m]")
plt.axis([s0, s1,-1,40])
plt.subplot(322)
plt.plot(S,BetaY,S_tr,BetaY_tr); plt.axis([s0, s1,-1,40])
plt.subplot(323)
plt.plot(S,AlphaX,S_tr,AlphaX_tr)
plt.xlabel("s[m]");plt.ylabel(r"$\alpha_x$")
plt.axis([s0, s1,-10,10])
plt.subplot(324)
plt.plot(S,AlphaY,S_tr,AlphaY_tr)
plt.xlabel("s[m]")
plt.axis([s0, s1,-10,10])

plt.subplot(325)
plt.plot(S,GammaX,S_tr,GammaX_tr)
plt.xlabel("s[m]");plt.ylabel(r"$\gamma_x$")
plt.axis([s0, s1,0,7])
plt.subplot(326)
plt.plot(S,GammaY,S_tr,GammaY_tr)
plt.xlabel("s[m]")
plt.axis([s0, s1,0,7])

n_out = len(S_tr)
out = np.zeros((n_out, 13))
out[:, 0] = S_tr
out[:, 1] = AlphaX_tr
out[:, 2] = BetaX_tr
out[:, 3] = GammaX_tr
out[:, 4] = AlphaY_tr
out[:, 5] = BetaY_tr
out[:, 6] = GammaY_tr
out[:, 7] = E_tr
out[:, 8] = X_tr
out[:, 9] = Y_tr
out[:, 10] = XX_tr
out[:, 11] = YY_tr
out[:, 12] = XY_tr
np.savetxt("Optics_01.txt", out)

plt.show()


