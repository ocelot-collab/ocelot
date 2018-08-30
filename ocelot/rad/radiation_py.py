__author__ = 'Sergey Tomin'
"""
can read different types of files. By default, mag_file is in [mm] .
for python version, only vertical component of magnetic field (By) is taken into account.
In order to overcome this limitation, someone have to change function radiation_py.field_map2field_func(z, By).
Sergey Tomin 04.11.2016.
"""

from scipy import interpolate
from ocelot.cpbd.elements import *
from ocelot.cpbd.track import *
from ocelot.rad.spline_py import *
from ocelot.common.globals import *
import copy

class Motion:
    pass

def bspline(x, y, x_new):
    tck = interpolate.splrep(x, y, s=0)
    ynew = interpolate.splev(x_new, tck, der=0)
    return ynew

def integ_beta2(x, y):
    A,B,C,D, Z = cspline_coef(x, y)
    b2 = 0.
    beta2 = [0.]
    for i in range(len(x)-1):
        h = x[i+1] - x[i]
        a= A[i]
        b = B[i]
        c = C[i]
        d = D[i]
        #print h, a,b,c,d
        b2 += h*(d*d + h*(c*d + h*(1./3.* (c*c + 2*b*d) + h*(0.5*(b*c + a*d) + h*(0.2*(b*b + 2*a*c) + h*(1./3.*a*b + (a*a*h)/7.))))))
        beta2.append( b2)
        #print beta2
    return np.array(beta2)

def x2xgaus(X):
    """
    transform coordinates for gauss integration
    | | | | -> | x x x . x x x . x x x |
    | - coordinate
    x - new coordinate
    . - removed coordinate |

    :param X: array
    :return: new array
    """
    sqrt35 = 0.5*np.sqrt(3./5.)
    xnew = [X[0]]
    h = X[1] - X[0]
    xgaus = np.array([0.5-sqrt35, 0.5-sqrt35 + sqrt35, 0.5-sqrt35 + sqrt35 + sqrt35])*h
    for x in X[:-1]:
        xnew = np.append(xnew, x + xgaus)
    xnew = np.append(xnew, X[-1])
    return xnew

#from scipy.integrate import simps
def traj2motion(traj):
    motion = Motion()
    motion.x = traj[0::9]
    motion.y = traj[2::9]
    motion.z = traj[4::9]
    motion.bx = traj[1::9]
    motion.by = traj[3::9]
    motion.bz = traj[5::9]
    motion.Bx = traj[6::9]
    motion.By = traj[7::9]
    motion.Bz = traj[8::9]
    #new_motion = Motion()

    motion.z = motion.z.flatten()

    Z = x2xgaus(motion.z)

    motion.x = bspline(motion.z, motion.x, Z)*1000.
    motion.y = bspline(motion.z, motion.y, Z)*1000.
    #print "inegr = ", simps(motion.bx.flatten()**2, motion.z*1000)
    Ibx2 = integ_beta2(motion.z*1000., motion.bx)
    Iby2 = integ_beta2(motion.z*1000., motion.by)

    motion.XbetaI2 = bspline(motion.z, Ibx2, Z)
    motion.YbetaI2 = bspline(motion.z, Iby2, Z)

    motion.bx = bspline(motion.z, motion.bx, Z)
    motion.by = bspline(motion.z, motion.by, Z)

    motion.Bx = bspline(motion.z, motion.Bx, Z)
    motion.By = bspline(motion.z, motion.By, Z)
    motion.z = Z*1000.
    #plt.plot(motion.bx.flatten()**2)
    #plt.show()
    return motion


def und_field(x, y, z, lperiod, Kx, nperiods=None):
    if nperiods is None:
        ph_shift = 0
        z_coef = 1
    else:
        ph_shift = np.pi / 2
        z_coef = (0.25 * np.heaviside(z, 0) + 0.5 * np.heaviside(z - lperiod / 2., 0) + 0.25 * np.heaviside(z - lperiod,0)
                  - 0.25 * np.heaviside(z - (nperiods - 1) * lperiod, 0) - 0.5 * np.heaviside(
                    z - (nperiods - 0.5) * lperiod, 0)
                  - 0.25 * np.heaviside(z - nperiods * lperiod, 0))
    kx = 0.
    kz = 2*pi/lperiod
    ky = np.sqrt(kz*kz - kx*kx)
    c = speed_of_light
    m0 = m_e_eV
    B0 = Kx*m0*kz/c
    k1 = -B0*kx/ky
    k2 = -B0*kz/ky

    kx_x = kx*x
    ky_y = ky*y
    kz_z = kz*z
    cosx = np.cos(kx_x)
    sinhy = np.sinh(ky_y)
    cosz = np.cos(kz_z + ph_shift)*z_coef
    Bx = k1*np.sin(kx_x)*sinhy*cosz #// here kx is only real
    By = B0*cosx*np.cosh(ky_y)*cosz
    Bz = k2*cosx*sinhy*np.sin(kz_z)
    return (Bx, By, Bz)



def energy_loss_und(energy, Kx, lperiod, L, energy_loss = False):
    if energy_loss:
        k = 4.*pi*pi/3.*ro_e/m_e_GeV
        #print "k = ", k
        U = k*energy**2*Kx**2*L/lperiod**2
    else:
        U = 0.
    return U

def sigma_gamma_quat(energy, Kx, lperiod, L):
    """
    rate of energy diffusion

    :param energy: electron beam energy
    :param Kx: undulator parameter
    :param lperiod: undulator period
    :param L: length
    :return: sigma_gamma/gamma
    """
    gamma = energy/m_e_GeV
    lambda_compt = 2.4263102389e-12 #m
    lambda_compt_r = lambda_compt/2./pi
    f = lambda K: 1.2 + 1./(K + 1.33*K*K + 0.4*K**3)
    delta_Eq2 = 56.*pi**3/15.*lambda_compt_r*ro_e*gamma**4/lperiod**3*Kx**3*f(Kx)*L
    sigma_Eq = np.sqrt(delta_Eq2/(gamma*gamma))
    return sigma_Eq


def quantum_diffusion(energy, Kx, lperiod, L, quantum_diff = False):
    if quantum_diff:
        # gamma = energy/m_e_GeV
        # lambda_compt = 2.4263102389e-12 # h_eV_s/m_e_eV*speed_of_light
        # lambda_compt_r = lambda_compt/2./pi
        # f = lambda K: 1.2 + 1./(K + 1.33*K*K + 0.4*K**3)
        # delta_Eq2 = 56.*pi**3/15.*lambda_compt_r*ro_e*gamma**4/lperiod**3*Kx**3*f(Kx)*L
        sigma_Eq = sigma_gamma_quat(energy, Kx, lperiod, L) # sqrt(delta_Eq2/(gamma*gamma))
        # print "sigma_q = ", sigma_Eq, energy
        U = sigma_Eq*np.random.randn()*energy
    else:
        U = 0.
    return U


def field_map2field_func(z, By):
    tck = interpolate.splrep(z, By, k=3)
    func = lambda x, y, z: (0, interpolate.splev(z, tck, der=0), 0)
    return func


def track4rad(beam, lat, energy_loss=False, quantum_diff=False, accuracy=1):
    """
    Function calculates the electron trajectory

    :param beam: Beam class
    :param lat: MagneticLattice class
    :param energy_loss: False, flag to calculate energy loss
    :param quantum_diff: False, flag to calculate quantum diffusion
    :param accuracy: 1, accuracy
    :return: U, E; U - list of u, u is 9xN array (6 coordinates and 3 mag field), E - list of energies
    """
    energy = beam.E
    #Y0 = [beam.x, beam.xp, beam.y, beam.yp, 0, 0]
    p = Particle(x=beam.x, px=beam.xp, y=beam.yp, py=beam.yp, E=beam.E)
    L = 0.
    U = []
    E = []
    #n = 0
    non_u = []
    #K = 4
    for elem in lat.sequence:
        if elem.l == 0:
            continue
        if elem.__class__ != Undulator:
            non_u.append(elem)
            U0 = 0.
        else:
            if len(non_u) != 0:
                lat_el = MagneticLattice(non_u)
                if lat_el.totalLen != 0:
                    navi = Navigator(lat)
                    u = []
                    N = 500
                    for z in np.linspace(L, lat_el.totalLen + L, num=N):
                        h = lat_el.totalLen/(N)
                        tracking_step(lat_el, [p], h, navi)
                        #print p.s
                        ui = [p.x, p.px, p.y, p.py, z, np.sqrt(1. - p.px*p.px - p.py *p.py), 0., 0., 0.]
                        u.extend(ui)
                    U.append(np.array(u))
                    E.append(energy)
                L += lat_el.totalLen
            non_u = []



            U0 = energy_loss_und(energy, elem.Kx, elem.lperiod, elem.l, energy_loss)
            Uq = quantum_diffusion(energy, elem.Kx, elem.lperiod, elem.l, quantum_diff)
            #print U0, Uq
            U0 = U0 + Uq
            #U0 = 8889.68503061*1e-9
            #print "U0 = ", U0*1e9," eV"
            #if energy_loss  == False:
            #    U0 = 0.
            mag_length = elem.l
            try:
                mag_field = elem.mag_field
            except:
                if len(elem.field_map.z_arr) != 0:
                    #print("Field_map exist! Creating mag_field(x, y, z)")
                    unit_coef = 0.001 if elem.field_map.units == "mm" else 1
                    mag_length = elem.field_map.l*unit_coef
                    z_array = (elem.field_map.z_arr - elem.field_map.z_arr[0])*unit_coef
                    mag_field = field_map2field_func(z=z_array, By=elem.field_map.By_arr)
                else:
                    #print("Standard undulator field")
                    mag_field = lambda x, y, z: und_field(x, y, z, elem.lperiod, elem.Kx, nperiods=elem.nperiods)
            N = int((mag_length*1500 + 100)*accuracy)
            p_array = ParticleArray()
            p_array.list2array([p])
            u = rk_track_in_field(p_array.rparticles, mag_length, N, energy, mag_field, s_start=0)
            p.x = u[-9, 0]
            p.px =u[-8, 0]
            p.y = u[-7, 0]
            p.py =u[-6, 0]
            s = u[-5, 0]
            u[4::9] += L
            L += s
            U.append(u)

            E.append(energy)
        energy = energy - U0
        #print energy
    return U, E


def gintegrator(Xscr, Yscr, Erad, motion, screen, n, n_end, gamma, half_step):
    """

    :param Xscr:
    :param Yscr:
    :param Erad:
    :param motion:
    :param screen:
    :param n:
    :param n_end:
    :param gamma:
    :param half_step:
    :return:
    """
    Q = speed_of_light/m_e_eV/1000 #0.5866740802042227#; // e/mc = (mm*T)^-1
    hc = h_eV_s*speed_of_light*1000# 1.239841874330e-3 # // mm
    k2q3 = 1.1547005383792517#;//  = 2./sqrt(3)
    gamma2 = gamma*gamma
    w = [0.5555555555555556*half_step, 0.8888888888888889*half_step, 0.5555555555555556*half_step]
    LenPntrConst = screen.Distance - motion.z[0]#; // I have to pay attention to this
    phaseConst = np.pi*Erad/(gamma2*hc)
    for p in range(3):# // Gauss integration
        i = n*3 + p + 1
        #radConstAdd = w[p]*Q*k2q3*(screen.Distance - motion.z[0] - 0*screen.Zstart)
        XX = motion.x[i]
        YY = motion.y[i]
        ZZ = motion.z[i]
        BetX = motion.bx[i]
        BetY = motion.by[i]
        IbetX2 = motion.XbetaI2[i]
        IbetY2 = motion.YbetaI2[i]
        Bx = motion.Bx[i]
        By = motion.By[i]
        LenPntrZ = screen.Distance - ZZ

        prX = Xscr - XX #//for pointer nx(z)
        prY = Yscr - YY #//for pointer ny(z)
        nx = prX/LenPntrZ
        ny = prY/LenPntrZ
        tx = gamma*(nx - BetX)
        ty = gamma*(ny - BetY)
        tx2 = tx*tx
        ty2 = ty*ty
        tyx = 2.*tx*ty

        radConst = w[p]*Q*k2q3*(screen.Distance)/LenPntrZ/((1. + tx2 + ty2)*(1. + tx2 + ty2))
        
        radX = radConst*(By*(1. - tx2 + ty2) + Bx*tyx - 2.*tx/Q/LenPntrZ)#/*sigma*/
        radY = -radConst*(Bx*(1. + tx2 - ty2) + By*tyx + 2.*ty/Q/LenPntrZ)#;/*pi*/


        prXconst = Xscr - motion.x[0]
        prYconst = Yscr - motion.y[0]
        phaseConstIn = (prXconst*prXconst + prYconst*prYconst)/LenPntrConst
        phaseConstCur = (prX*prX + prY*prY)/LenPntrZ
        #// string below is for case direct accumulation
        #//double phase = screen->Phase[ypoint*xpoint*je + xpoint*jy + jx] + faseConst*(ZZ - motion->Z[0]  + gamma2*(IbetX2 + IbetY2 + phaseConstCur - phaseConstIn));

        phase = phaseConst*(ZZ - motion.z[0] + gamma2*(IbetX2 + IbetY2 + phaseConstCur - phaseConstIn)) + screen.arPhase

        cosf = np.cos(phase)
        sinf = np.sin(phase)
        EreX = radX*cosf #//(cosf *cos(fase0) - sinf*sin(fase0));
        EimX = radX*sinf #//(sinf *cos(fase0) + cosf*sin(fase0));
        EreY = radY*cosf
        EimY = radY*sinf

        screen.arReEx += EreX
        screen.arImEx += EimX
        screen.arReEy += EreY
        screen.arImEy += EimY
        if  i == n_end:# //(n == 5000 && p == 2)
            LenPntrZ = screen.Distance - motion.z[-1]
            prX = Xscr - motion.x[-1] # //for pointer nx(z)
            prY = Yscr - motion.y[-1] # //for pointer ny(z)
            IbetX2 = motion.XbetaI2[-1]
            IbetY2 = motion.YbetaI2[-1]
            phase = phaseConst*(motion.z[-1] - motion.z[0]  + gamma2*(IbetX2 + IbetY2 + prX*prX/LenPntrZ + prY*prY/LenPntrZ - phaseConstIn))
            screen.arPhase = screen.arPhase + phase
    return screen

#import matplotlib.pyplot as plt

def radiation_py(gamma, traj, screen):
    """
    screen format     screen->ReEx[ypoint*xpoint*je + xpoint*jy + jx] += EreX;
    """

    motion = traj2motion(traj)

    size = len(motion.z)
    Nmotion = int((size + 1)/3)
    half_step = (motion.z[-1]-motion.z[0])/2./(Nmotion-1)

    n_end = len(motion.z)-2
    Xscr = np.linspace(screen.x_start, screen.x_start+screen.x_step*(screen.nx-1), num = screen.nx)
    Yscr = np.linspace(screen.y_start, screen.y_start+screen.y_step*(screen.ny-1), num = screen.ny)
    Yscr = Yscr.reshape((screen.ny, 1))
    Erad = np.linspace(screen.e_start, screen.e_start+screen.e_step*(screen.ne-1), num = screen.ne)
    Erad = Erad.reshape((screen.ne, 1))

    shape_array = [screen.ne, screen.ny, screen.nx]
    if 1 in shape_array:
        if screen.ny >1 and screen.ne>1:
            Yscr = Yscr.reshape((1, screen.ny))
        shape_array.remove(1)
        screen.arReEx = screen.arReEx.reshape(shape_array)
        screen.arImEx = screen.arImEx.reshape(shape_array)
        screen.arReEy = screen.arReEy.reshape(shape_array)
        screen.arImEy = screen.arImEy.reshape(shape_array)
        screen.arPhase = screen.arPhase.reshape(shape_array)

        for n in range(Nmotion-1):
            screen = gintegrator(Xscr, Yscr, Erad, motion, screen, n, n_end, gamma, half_step)

        screen.arReEx = screen.arReEx.flatten()
        screen.arImEx = screen.arImEx.flatten()
        screen.arReEy = screen.arReEy.flatten()
        screen.arImEy = screen.arImEy.flatten()
        screen.arPhase = screen.arPhase.flatten()
    else:
        print( "SR 3D calculation")
        arReEx = np.array([])
        arImEx = np.array([])
        arReEy = np.array([])
        arImEy = np.array([])
        arPhase = np.array([])

        for erad in Erad:
            screen.arReEx = np.zeros((screen.ny, screen.nx))
            screen.arImEx = np.zeros((screen.ny, screen.nx))
            screen.arReEy = np.zeros((screen.ny, screen.nx))
            screen.arImEy = np.zeros((screen.ny, screen.nx))
            screen.arPhase =np.zeros((screen.ny, screen.nx))

            for n in range(Nmotion-1):
                screen = gintegrator(Xscr, Yscr, erad, motion, screen, n, n_end, gamma, half_step)

            arReEx = np.append(arReEx,screen.arReEx.flatten())
            arImEx = np.append(arImEx,screen.arImEx.flatten())
            arReEy = np.append(arReEy,screen.arReEy.flatten())
            arImEy = np.append(arImEy,screen.arImEy.flatten())
            arPhase= np.append(arPhase,screen.arPhase.flatten())
        screen.arReEx = arReEx
        screen.arImEx = arImEx
        screen.arReEy = arReEy
        screen.arImEy = arImEy
        screen.arPhase = arPhase
    return 1


#from ocelot.gui import *
#import matplotlib.pyplot as plt
def calculate_radiation(lat, screen, beam, energy_loss=False, quantum_diff=False, accuracy=1):
    screen.update()
    b_current = beam.I*1000. # b_current - beam current must be in [mA], but beam.I in [A]
    energy = beam.E
    gamma = energy/m_e_GeV

    screen.nullify()
    U, E = track4rad(beam, lat, energy_loss=energy_loss, quantum_diff=quantum_diff, accuracy=accuracy)
    for u, e in zip(U, E):
        # plt.plot(u[4::9], u[0::9], "r")
        # plt.plot(u[4::9], u[2::9], "b")
        radiation_py(e/m_e_GeV, u, screen)
    # plt.show()
    screen.distPhoton(gamma, current=b_current)
    screen.Ef_electron = E[-1]
    screen.motion = U
    return screen


def track4rad_beam(p_array, lat, energy_loss=False, quantum_diff=False, accuracy=1):
    """
    Function calculates the electron trajectory

    :param beam: Beam class
    :param lat: MagneticLattice class
    :param energy_loss: False, flag to calculate energy loss
    :param quantum_diff: False, flag to calculate quantum diffusion
    :param accuracy: 1, accuracy
    :return: U, E; U - list of u, u is 9xN array (6 coordinats and 3 mag field), E - list of energies
    """
    energy = p_array.E
    #Y0 = [beam.x, beam.xp, beam.y, beam.yp, 0, 0]
    #p = Particle(x=beam.x, px=beam.xp, y=beam.yp, py=beam.yp, E=beam.E)
    L = 0.
    U = []
    E = []
    non_u = []
    for elem in lat.sequence:
        if elem.l == 0:
            continue
        if elem.__class__ != Undulator:
            non_u.append(elem)
            U0 = 0.
        else:
            if len(non_u) != 0:
                lat_el = MagneticLattice(non_u)
                if lat_el.totalLen != 0:
                    navi = Navigator(lat)

                    N = 500
                    u = np.zeros((N * 9, np.shape(p_array.rparticles)[1]))
                    for i, z in enumerate(np.linspace(L, lat_el.totalLen + L, num=N)):
                        h = (lat_el.totalLen)/(N)
                        tracking_step(lat_el, p_array, h, navi)
                        u[i * 9 + 0, :] = p_array.rparticles[0]
                        u[i * 9 + 1, :] = p_array.rparticles[1]
                        u[i * 9 + 2, :] = p_array.rparticles[2]
                        u[i * 9 + 3, :] = p_array.rparticles[3]
                        u[i * 9 + 4, :] = p_array.rparticles[4] + z
                        u[i * 9 + 5, :] = p_array.rparticles[5]
                        #ui = [p.x, p.px, p.y, p.py, z, np.sqrt(1. - p.px*p.px - p.py *p.py), 0., 0., 0.]
                        #u.extend(ui)

                    U.append(u)
                    E.append(energy)
                L += lat_el.totalLen
            non_u = []

            U0 = energy_loss_und(energy, elem.Kx, elem.lperiod, elem.l, energy_loss)
            Uq = quantum_diffusion(energy, elem.Kx, elem.lperiod, elem.l, quantum_diff)
            U0 = U0 + Uq

            mag_length = elem.l
            try:
                mag_field = elem.mag_field
            except:
                if len(elem.field_map.z_arr) != 0:
                    #print("Field_map exist! Creating mag_field(x, y, z)")
                    unit_coef = 0.001 if elem.field_map.units == "mm" else 1
                    mag_length = elem.field_map.l*unit_coef
                    z_array = (elem.field_map.z_arr - elem.field_map.z_arr[0])*unit_coef
                    mag_field = field_map2field_func(z=z_array, By=elem.field_map.By_arr)
                else:
                    #print("Standard undulator field")
                    mag_field = lambda x, y, z: und_field(x, y, z, elem.lperiod, elem.Kx, nperiods=elem.nperiods)
            N = int((mag_length*1500 + 100)*accuracy)
            u = rk_track_in_field(p_array.rparticles, mag_length, N, energy, mag_field)

            p_array.x()[:] = u[-9, :]
            p_array.px()[:] =u[-8, :]
            p_array.y()[:] = u[-7, :]
            p_array.py()[:] =u[-6, :]
            s = u[-5, 0]
            u[4::9] += L
            L += s
            U.append(u)
            E.append(energy)
        energy = energy - U0
    #for u in U:
    #    print("here", len(u[4::9, 0]))
    #    plt.plot(u[4::9, :], u[7::9, :])
    #plt.show()
    return U, E


def calculate_beam_radiation(lat, screen, p_array, energy_loss=False, quantum_diff=False, accuracy=1, freq=1):
    screen.update()
    b_current = p_array.q_array[0]*1000.*freq # b_current - beam current must be in [mA], but beam.I in [A]
    energy = p_array.E
    gamma = energy/m_e_GeV
    tau0 = np.copy(p_array.tau())
    p_array.tau()[:] = 0

    U, E = track4rad_beam(p_array, lat, energy_loss=energy_loss, quantum_diff=quantum_diff, accuracy=accuracy)

    screen.nullify()

    screen2 = copy.deepcopy(screen)
    for i in range(len(p_array.x())):
        screen_copy = copy.deepcopy(screen2)
        screen_copy.nullify()

        wlengthes = h_eV_s*speed_of_light/screen_copy.Eph
        screen_copy.arPhase[:] = tau0[i]/wlengthes*2*np.pi

        for u, e in zip(U, E):
            #shift_x = (max(u[0::9, i]) + min(u[0::9, i]))/2
            #u[0::9, i] -= shift_x
            # plt.plot(u[4::9, i], u[0::9, i])
            # plt.show()
            u[4::9, i] -= U[0][4, i]
            gamma = (1 + p_array.p()[i]) * e / m_e_GeV
            # print(i, gamma, p_array.p()[i], screen_copy.Distance)

            radiation_py(gamma, u[:, i], screen_copy)
            #screen_copy.distPhoton(gamma, current=b_current)
            #show_flux(screen_copy, unit="mrad", title=str(i), nfig=i)
            #plt.plot(screen_copy.arReEx**2 + screen_copy.arImEx**2, label="Re " + str(i))
            #plt.plot(screen_copy.arImEx**2, label="Im " + str(i))
            # plt.title("90 periods, 20 km, Eph = 8.5044 meV")
            # plt.plot(screen_copy.Xph/1000, screen_copy.arPhase % 2*np.pi, label=r"$\phi \% 2 \pi$")
            # plt.plot(screen_copy.Xph/1000, np.arctan(screen_copy.arImEx/screen_copy.arReEx), label=r"$\arctan\left(\frac{Im(E_x)}{Re(E_x)}\right)$")
            #plt.plot(screen_copy.Yph/1000, screen_copy.arReEy, label=r"$Re(E_x)$")
            # plt.plot(screen_copy.Xph/1000, screen_copy.arImEx, label=r"$Im(E_x)$")
            # plt.plot(screen_copy.Xph/1000, np.sqrt(screen_copy.arImEx**2 + screen_copy.arReEx**2), label=r"$\sqrt{Re(E_x)^2 + Im(E_x)^2}$")
            # plt.grid(True)
            # plt.xlabel("m")
            # plt.legend()
            # plt.show()
            screen.arReEx += screen_copy.arReEx
            screen.arImEx += screen_copy.arImEx
            screen.arReEy += screen_copy.arReEy
            screen.arImEy += screen_copy.arImEy
            # screen.arPhase[0:10]/2./pi
    #plt.plot(screen.arReEx ** 2 + screen.arImEx ** 2, label="Re " + str(i))

    #plt.legend()
    #plt.show()
    screen.distPhoton(gamma, current = b_current)
    screen.Ef_electron = E[-1]
    screen.motion = U
    return screen


def calculate_beam_radiation_check(lat, screen, p_array, energy_loss=False, quantum_diff=False, accuracy=1, freq=1):


    screen.update()
    b_current = p_array.q_array[0]*1000.*freq # b_current - beam current must be in [mA], but beam.I in [A]
    energy = p_array.E
    gamma = energy/m_e_GeV

    screen.nullify()
    #U, E = track4rad_beam(p_array, lat, energy_loss=energy_loss, quantum_diff=quantum_diff, accuracy=accuracy)
    #print("DONE tracking")
    for i in range(len(p_array.x())):
        print(i)
        screen_copy = copy.deepcopy(screen)
        screen_copy.nullify()
        lat.sequence[0].l += p_array.tau()[i]
        lat.update_transfer_maps()
        p_array_i = ParticleArray(n=1)
        p_array_i.s = p_array.s
        p_array_i.E = p_array.E
        p_array_i.q_array = p_array.q_array[i]
        p_array_i.rparticles[:, 0] = p_array.rparticles[:, i]
        p_array_i.rparticles[4] = 0
        U, E = track4rad_beam(p_array_i, lat, energy_loss=energy_loss, quantum_diff=quantum_diff, accuracy=accuracy)
        #for u in U:
        #    plt.plot(u[4::9, :], u[0::9, :])
        for u, e in zip(U, E):
            #for i in range(np.shape(u)[1]):

            #print "Energy = ", e, e/m_e_GeV
            #plot(u[4::9], u[0::9],".-")
            gamma = (1 + p_array.p()[0])*e/m_e_GeV
            #print(i, gamma)
            radiation_py(gamma, u[:, 0], screen_copy)
            #print("phase = ", i, screen_copy.arPhase)
            screen.arReEx += screen_copy.arReEx
            screen.arImEx += screen_copy.arImEx
            screen.arReEy += screen_copy.arReEy
            screen.arImEy += screen_copy.arImEy
            # screen.arPhase[0:10]/2./pi
    #plt.show()
    screen.distPhoton(gamma, current = b_current)
    screen.Ef_electron = E[-1]
    screen.motion = U
    return screen


"""
void sum_screens(Screen *screen, Screen screen_up)
{
    //Screen screen;
    int size = screen->eNstep*screen->xNstep*screen->yNstep;

    for(int i = 0; i<size; i++)
    {
        double sinfa = sin(screen->Phase[i]);
        double cosfa = cos(screen->Phase[i]);
        screen->ReEx[i] += screen_up.ReEx[i]*cosfa - screen_up.ImEx[i]*sinfa;
        screen->ImEx[i] += screen_up.ImEx[i]*cosfa + screen_up.ReEx[i]*sinfa;
        screen->ReEy[i] += screen_up.ReEy[i]*cosfa - screen_up.ImEy[i]*sinfa;
        screen->ImEy[i] += screen_up.ImEy[i]*cosfa + screen_up.ReEy[i]*sinfa;
        screen->Phase[i] += screen_up.Phase[i];
    }
}
"""
if __name__ == "__main__":
    quantum_diffusion(17.5, 4., 0.04, 200. ,quantum_diff=True)
    x = np.linspace(0, 1, 4)
    xnew = x2xgaus(x)
    print( x)
    print( xnew)
