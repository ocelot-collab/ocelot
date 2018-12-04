#import sys
#sys.path.append("/Volumes/Promise RAID/UserFolders/zagor_xxl/ocelot")

import matplotlib
matplotlib.use('TkAgg')
from accelerator.s2e_sections.sections import *
from ocelot.utils.section_track import *
from ocelot.gui.accelerator import *
import time
from ext4s2e import *

#data_dir = "N:/4all/xxl/zagor/mpy_xxl"
data_dir = "/Users/zagor/dataxxl/ocelot"
#data_dir = "/Volumes/Promise RAID/UserFolders/zagor_xxl/dataxxl/ocelot"

N_iter = 1
SaveRF = None #"RF.txt"
LoadRF = None #"RF_250_5.txt"
Restart=False
N_part_reduce = 10


# 500pC 5kA
E10=130; E20=700;  E30=2400; E40=14000
R3=700; R2=0; bounds=(-1, 1)
#r1=4.1218; r2= 8.3934; r3=14.7982;C10= 3; C20= 7; C30=192/(C10*C20);R3=1000;

# 500pC 10kA
#r1=4.1218; r2= 8.3934; r3=14.3639;C10= 3; C20= 7; C30=2*192/(C10*C20);R3=1000;

# 250pC 5kA
r1=4.1218; r2= 8.3934; r3=14.4111;C10= 3; C20= 7; C30=345/(C10*C20);R3=700;

CheckOptics=True
FirstOrderOptics=True

match_exec=False
smooth_exec=False
wake_exec=False
SC_exec=False
CSR_exec=False



c=299792458; grad=pi/180
f = 1.3e9; k=2*pi*f/c
all_sections = [A1, AH1, LH, DL,  BC0, L1, BC1, L2, BC2, L3, CL1, CL2, CL3, STN10]#, SASE1, T4, SASE3, T4D]

######################### design optics
tws0 = Twiss()
tws0.E = 0.005
tws0.beta_x  = 0.286527307369
tws0.beta_y  = 0.286527307369
tws0.alpha_x = -0.838833736086
tws0.alpha_y = -0.838833736086
section_lat = SectionLattice(sequence=all_sections, tws0=tws0, data_dir=data_dir)

#lat = MagneticLattice(section_lat.elem_seq)
#plot_opt_func(lat, section_lat.tws, top_plot=["Dx","Dy"],legend=False)
#plt.show()

if CheckOptics:
    r1=3.509851167058248;    r2=9.275657479513360;    r3=10.979999976659300
    V11 = 150 - 6.5;    V13 = 20;    V21 = 700 - 130;    V31 = 2400 - 700;
    fi11 = 0;    fi13 = 180 * grad;    fi21 = 0;    fi31 = 0;
########################### tracking

#sections = [A1, AH1, LH, DL,  BC0, L1, BC1, L2, BC2, L3, CL1, CL2, CL3, STN10]
#sections = [A1, AH1, LH, DL,  BC0, L1, BC1, L2, BC2,L3]
sections = [DL]

if FirstOrderOptics:
    for s in sections:
        #section_lat.dict_sections[s].method.global_method = SecondTM
        section_lat.dict_sections[s].method.global_method = TransferMap
        section_lat.dict_sections[s].lattice.update_transfer_maps()

start = time.time()

#################################### Start S2E ################################
p_array0 = load_particle_array(filename=data_dir+"/particles/N10_out_250pC.ast")
#p_array0 = load_particle_array(filename=data_dir+"/particles/N1_out_3p20_500pC.ast")
p_array = ParticleArray(n=round(200000/N_part_reduce))
p_array.x()[:]=p_array0.x()[::N_part_reduce]
p_array.y()[:]=p_array0.y()[::N_part_reduce]
p_array.px()[:]=p_array0.px()[::N_part_reduce]
p_array.py()[:]=p_array0.py()[::N_part_reduce]
p_array.tau()[:]=p_array0.tau()[::N_part_reduce]
p_array.p()[:]=p_array0.p()[::N_part_reduce]
p_array.q_array[:]=p_array0.q_array[::N_part_reduce]*N_part_reduce
p_array.s=3.2
p_array.E=p_array0.E
del p_array0

inds = np.argsort(p_array.tau()[:])
inds0 = np.argsort(inds)
ref_i0=inds0[0]
p_array.x()[:]=p_array.x()[inds]
p_array.y()[:]=p_array.y()[inds]
p_array.px()[:]=p_array.px()[inds]
p_array.py()[:]=p_array.py()[inds]
p_array.tau()[:]=p_array.tau()[inds]
p_array.p()[:]=p_array.p()[inds]
p_array.q_array[:]=p_array.q_array[inds]
print("ref_i0=",ref_i0, " ref_tau=",p_array.tau()[ref_i0], " ref_p=",p_array.p()[ref_i0])

p_array.y()[:] -= np.mean(p_array.y())
p_array.x()[:] -= np.mean(p_array.x())
p_array.px()[:] -= np.mean(p_array.px())
p_array.py()[:] -= np.mean(p_array.py())
Q=np.sum(p_array.q_array)

if CheckOptics:
    p_array.rparticles[5,:]=0
    p_array.rparticles[4,:]=0
else:
    show_e_beam(p_array, nparts_in_slice=10000, nfig=1, inverse_tau=False, title="", figsize=(10, 8))
save_particle_array(data_dir+"/particles/gun.npz", p_array)

gamref = p_array.E / m_e_GeV
betaref = np.sqrt(1 - gamref ** -2)
s0 = p_array.s
P = p_array.rparticles.view()
Np = int(P.size / 6)
xp = exact_xxstg_2_xp_mad(P, gamref)
E00=p_array.E*1e3
P0=np.zeros((Np,2))
P0[:, 0] = xp[:, 2]
P0[:, 1] = (1+P[5,:]*betaref)*E00


sig0=np.std(P0[:,0])
z0=np.mean(P0[:,0])
inds=np.argwhere(np.logical_and(P0[:,0]>z0+sig0*bounds[0], P0[:,0]<z0+sig0*bounds[1]));
inds=inds.reshape(inds.shape[0])

config = {
    BC0: {"rho": r1},
    BC1: {"rho": r2},
    BC2: {"rho": r3},
}

section_lat.update_sections(sections, config=config)

f0=[1/C10, 1/(C20*C10), 1/(C10*C20*C30), E10, E20, E30, R2, R3]
if LoadRF:
    RFpars=np.loadtxt(LoadRF)
    if Restart:
        fn=np.loadtxt('fn.txt')
    else:
        fn = f0.copy()
    V11=RFpars[0,0]; fi11=RFpars[0,1]*grad
    V13 = RFpars[1, 0]; fi13 = RFpars[1, 1] * grad
    V21 = RFpars[2, 0]; fi21 = RFpars[2, 1] * grad
    V31=RFpars[3,0]; fi31=RFpars[3,1]*grad
else:
    if not CheckOptics:
        fn = f0.copy()
        fi11, V11, fi13, V13, fi21, V21, fi31, V31 = FindParametersThreeBCs(P0[inds, :],E00,
            fn[3], fn[4], fn[5], r56_1c, t566_1c, u5666_1c, r56_2, t566_2, u5666_2,
            r56_3, t566_3, u5666_3, f, fn[0], fn[1], fn[2], fn[6], fn[7])


Res=np.zeros((N_iter,8))
for i in range(N_iter):
    config = {
        A1: {"phi": fi11/grad, "v": V11/8*1e-3,
            "SC": SC_exec, "smooth": True, "wake": wake_exec},
        AH1: {"phi": fi13/grad, "v": V13/8*1e-3,
            "match": True, "SC": SC_exec, "wake": wake_exec},
        LH: {"SC": SC_exec,"CSR": CSR_exec, "wake": wake_exec},
        DL: {"match": match_exec, "SC": SC_exec,"CSR": CSR_exec, "wake": wake_exec},
        BC0: {"rho": r1,
            "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        L1: {"phi": fi21/grad, "v": V21/32*1e-3,
            "SC": SC_exec, "wake": wake_exec, "smooth": smooth_exec},
        BC1: {"rho": r2,
            "match": match_exec, "SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        L2: {"phi": fi31/grad, "v": V31/96*1e-3,
             "SC": SC_exec, "wake": wake_exec,"smooth": smooth_exec},
        BC2: {"rho": r3,
            "match": match_exec,"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        L3: {"SC": SC_exec, "wake": wake_exec, "smooth": smooth_exec},
        CL1: {"match": match_exec,"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        CL2: {"match": match_exec,"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        CL3: {"match": match_exec,"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        STN10: {"match": match_exec,"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        SASE1: {"match": match_exec,"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec},
        T4: {"match": match_exec,"SC": SC_exec, "CSR": CSR_exec, "wake": wake_exec}
    }

####################### TRACKING ###############################

    #p_array = load_particle_array(data_dir+"/particles/gun.npz")
    p_array = load_particle_array(data_dir + "/particles/section_LH.npz")
    p_array.p()[:] = p_array.p()[:]+0.1
    s00 = copy(p_array.s)
    print('s00=',s00)


    p_array = section_lat.track_sections(sections=sections, p_array=p_array, config=config, force_ext_p_array=True)



    if not CheckOptics:
        Z1, Z2, Z3, E1, E2, E3, p2, p3, SigZ=FindGlobalCompressionParameters(bounds,ref_i0, data_dir)
        f1 = [Z1, Z2, Z3, E1, E2, E3, p2, p3]
        Res[i]=f1
        for j in range(8):
            fn[j] = fn[j] + 1*(f0[j] - f1[j])
        fi11, V11, fi13, V13, fi21, V21, fi31, V31 = FindParametersThreeBCs(P0[inds, :],E00,
            fn[3], fn[4], fn[5], r56_1c, t566_1c, u5666_1c, r56_2, t566_2, u5666_2,
            r56_3, t566_3, u5666_3, f, fn[0], fn[1], fn[2], fn[6], fn[7])
        print("\nC1=", 1/Z1, " C2=",1/Z2," C3=",1/Z3)
        print("E1=",E1," E2=", E2, " E3=",E3)
        print("p2=",p2," p3=", p3)

RFpars=[[V11, fi11/grad],[V13, fi13/grad], [V21, fi21/grad], [V31, fi31/grad]]
print("RFpars:\n", RFpars)
if not CheckOptics:
    print("Std_Z:\n", SigZ)
    if SaveRF:
        np.savetxt(SaveRF,RFpars)
        np.savetxt("fn.txt",fn)
    plt.figure()
    plt.plot(1/Res[:,2])


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
EmitX_tr = [tw.emit_x for tw in tws_track_global]
EmitY_tr = [tw.emit_y for tw in tws_track_global]
E_tr=[tw.E for tw in tws_track_global]
sig_tau = np.sqrt([tw.tautau for tw in tws_track_global])
if CheckOptics:
    sig_tau[:]=1
I = c*Q/np.sqrt(2*pi)/sig_tau
Sx=np.zeros(len(S_tr))
Sy=np.zeros(len(S_tr))
for i in range(len(S_tr)):
    gamma=E_tr[i-1]/m_e_GeV
    EmitX_tr[i-1]=EmitX_tr[i-1]*gamma*1e6
    EmitY_tr[i-1]=EmitY_tr[i-1]*gamma*1e6
    Sx[i-1] = I[i-1]  * BetaX_tr[i-1] / (1.7045e+04 * gamma * gamma * EmitX_tr[i-1]*1e-6)
    Sy[i - 1] = I[i - 1] * BetaY_tr[i - 1] / (1.7045e+04 * gamma * gamma * EmitY_tr[i - 1] * 1e-6)


#Plots
s0 = 45; s1 = 75
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

plt.figure()
plt.subplot(211)
plt.plot(S_tr,EmitX_tr,S_tr,EmitY_tr)
plt.grid(True); plt.title("Emittances"); plt.xlabel("s[m]"); plt.ylabel(r"$\epsilon_x,\epsilon_y [\mu m]$")
plt.axis([s0, s1,0,1.7])
plt.subplot(212)
plt.plot(S,E,S_tr,E_tr)
plt.grid(True); plt.title("Energy"); plt.xlabel("s[m]"); plt.ylabel(r"E [GeV]$")
plt.axis([s0, s1, 0, 3])

plt.figure()
plt.subplot(211)
plt.plot(S_tr,sig_tau*1e3)
plt.grid(False); #plt.title("Sigma_tau");
plt.ylabel("$\sigma_z$ [mm]")
plt.axis([s0, s1,0,2.1])
plt.subplot(212)
plt.plot(S_tr,Sx,S_tr,Sy)
plt.plot(S_tr,Sx, "b", label=r"$k_x^{sc}$")
plt.plot(S_tr,Sy, "r", label=r"$k_y^{sc}$")
plt.legend()
plt.grid(False);
plt.xlabel("z[m]")
plt.axis([s0, s1,0,1.25])

#show_e_beam(p_array, nparts_in_slice=10000, smooth_param=0.03, nfig=13, title="")
n_out = len(S_tr)
out = np.zeros((n_out, 8))
out[:, 0] = S_tr
out[:, 0] = out[:, 0]
out[:, 1] = AlphaX_tr
out[:, 2] = BetaX_tr
out[:, 3] = GammaX_tr
out[:, 4] = AlphaY_tr
out[:, 5] = BetaY_tr
out[:, 6] = GammaY_tr
out[:, 7] = E_tr
np.savetxt("Optics_00.txt", out)

plt.show()


