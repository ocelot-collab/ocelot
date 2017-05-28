
import numpy as np
from ocelot.common.globals import *
from ocelot.cpbd.beam import ParticleArray
from scipy import interpolate


def FindRefParticle(PD, col, weight, m=6):

    for i in range(0, m):
        x0 = np.mean(PD[:, i])
        PD[:, i] = PD[:, i] - x0
        sig0 = np.std(PD[:, i])
        PD[:, i] = np.abs(PD[:, i])
        if sig0>0:
            PD[:, i] = PD[:, i]/sig0

    PD[:, col] = PD[:, col]*weight
    [r, i0] = np.min(np.sum(np.transpose(PD)))
    return i0, r


def sortrows(x, col):
    return x[x[:, col].argsort()]


def SaveAstraParticles(outfile, PD, Q, findref=False, z0=None):
    #%SaveAstraParticles(outfile,PD,Q,findref,z0)
    n = len(PD[:,0])
    if z0 == None:
        findref=False
    if findref:
        i0, r = FindRefParticle(PD, 3, 10)
        P0 = PD[0,:]
        PD[0, :] = PD[i0, :]
        PD[i0, :] = P0
        PD[2:n, :] = sortrows(PD[1:n, :], 3)

    PD[1:n, 3-1] = PD[1:n, 3-1] - PD[0, 3-1]
    PD[1:n, 6-1] = PD[1:n, 6-1] - PD[0, 6-1] #substract reference particle
    PD1 = np.zeros((n, 10))
    PD1[:, 0:6] = PD[:, 0:6]
    PD1[:, 7] = -Q/n
    PD1[:, 8] = 1
    PD1[:, 9] = 5

    if z0 !=None:
        PD1[0,2] = z0
    np.save(outfile,PD1)


def load_Astra_particles(filename):
    PD = np.loadtxt(filename)
    n = len(PD[:, 0])
    Q = np.abs(np.sum(PD[:,7]))
    PD[1:n, 2] = PD[1:n, 2] + PD[0, 2]
    PD[1:n, 5] = PD[1:n, 5] + PD[0, 7] # add reference particle
    PD1 = PD[:, 0:6]
    return PD1, Q

def csrtrackBeam2particleArray(filename, orient="H"):
    #H z x y pz px py -> x y z px py pz
    #V z y x pz py px -> x y -z px py -pz
    PD = np.loadtxt(filename)
    #PD = load(infile)
    n = np.shape(PD)[0] - 1 #length(PD(:,1))- 1
    #print(n)
    PD1 = np.zeros((n, 6))
    Q = np.sum(PD[1:, 6])*1e9
    t0 = PD[0, 0]
    if orient=='H':
       PD1[:, 1-1] = PD[1:, 2-1]
       PD1[:, 2-1] = PD[1:, 3-1]
       PD1[:, 3-1] = PD[1:, 1-1]
       PD1[:, 4-1] = PD[1:, 5-1]
       PD1[:, 5-1] = PD[1:, 6-1]
       PD1[:, 6-1] = PD[1:, 4-1]
    else:
       PD1[:, 1-1] = -PD[1:, 3-1]
       PD1[:, 2-1] =  PD[1:, 2-1]
       PD1[:, 3-1] =  PD[1:, 1-1]
       PD1[:, 4-1] = -PD[1:, 6-1]
       PD1[:, 5-1] =  PD[1:, 5-1]
       PD1[:, 6-1] =  PD[1:, 4-1]
    #print("CSR", PD1[0, :])
    for i in range(6):
        PD1[1:n, i] = PD1[1:n, i] + PD1[0, i]

    p_ref = np.sqrt(PD1[0, 3]**2 + PD1[0, 4]**2 + PD1[0, 5]**2)
    px = PD1[:, 3]/ p_ref
    py = PD1[:, 4] / p_ref
    Eref = np.sqrt(m_e_eV ** 2 + p_ref ** 2)
    pe = (np.sqrt(m_e_eV**2 + (PD1[:, 3]**2 + PD1[:, 4]**2 + PD1[:, 5]**2)) - Eref) / p_ref

    p_array = ParticleArray(n)
    p_array.rparticles[0] = PD1[:, 0] -    0*PD1[0, 0]
    p_array.rparticles[2] = PD1[:, 1] -   0*PD1[0, 1]
    p_array.rparticles[4] = -(PD1[:, 2] - PD1[0, 2])
    p_array.rparticles[1] = px[:]
    p_array.rparticles[3] = py[:]
    p_array.rparticles[5] = pe[:]


    p_array.q_array[:] = PD[1:, 6]
    p_array.s = PD1[0, 2]
    p_array.E = Eref*1e-9
    return p_array

#def xyz2ParticleArray():









if __name__ == "__main__":
    import matplotlib.pyplot  as plt
    filename = "C:\\Users\\tomins\\Documents\\Dropbox\\DESY\\repository\\For_Sergey\\test_ocelot_bc2\\N5_BC3\\out\\out.fmt1"
    pd, Q, t = load_CSRtrack_particles(filename, orient="H")
    sig_2 = np.std(pd[:, 2])
    print(sig_2*0.003)
    pd[:, 2] = pd[:, 2] - np.mean(pd[:, 2])
    B = s_to_cur(pd[:, 2], 0.003*sig_2, Q*1e-9, v=299792458)

    print(B[:10, 0])
    print(B[:10, 1])
    plt.plot(B[:,0],B[:, 1],'b')
    plt.show()