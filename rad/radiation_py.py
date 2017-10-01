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
    return array(beta2)

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


def und_field(x, y, z, lperiod, Kx):
    kx = 0.
    kz = 2*pi/lperiod
    ky = np.sqrt(kz*kz - kx*kx)
    c = speed_of_light
    m0 = m_e_eV
    B0 = Kx*m0*kz/c
    k1 =  -B0*kx/ky
    k2 = -B0*kz/ky

    kx_x = kx*x
    ky_y = ky*y
    kz_z = kz*z
    cosx = np.cos(kx_x)
    sinhy = np.sinh(ky_y)
    cosz = np.cos(kz_z)
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
    sigma_Eq = sqrt(delta_Eq2/(gamma*gamma))
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
                #print elem.type, elem.l, L
                lat_el = MagneticLattice(non_u)
                if lat_el.totalLen != 0:
                    navi = Navigator()
                    u = []
                    N = 500
                    for z in linspace(L, lat_el.totalLen + L, num=N):
                        h = lat_el.totalLen/(N)
                        tracking_step(lat_el, [p], h, navi)
                        #print p.s
                        ui = [p.x, p.px, p.y, p.py, z, np.sqrt(1. - p.px*p.px - p.py *p.py), 0., 0., 0.]
                        u.extend(ui)
                    U.append(array(u))
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
                    mag_field = lambda x, y, z: und_field(x, y, z, elem.lperiod, elem.Kx)
            N = int((mag_length*1500 + 100)*accuracy)

            u = rk_track_in_field(array([p.x, p.px, p.y, p.py, 0, 0]), mag_length, N, energy, mag_field)
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

#@jit
def gintegrator(Xscr, Yscr, Erad, motion, screen, n, n_end, gamma, half_step, tmp):
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
    :param tmp:
    :return:
    """
    Q = 0.5866740802042227#; // (mm*T)^-1
    hc = 1.239841874330e-3 # // mm
    k2q3 = 1.1547005383792517#;//  = 2./sqrt(3)
    gamma2 = gamma*gamma
    w = [0.5555555555555556*half_step, 0.8888888888888889*half_step, 0.5555555555555556*half_step]
    LenPntrConst = screen.Distance - motion.z[0]#; // I have to pay attention to this
    #print screen.Distance
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
        #print screen.arPhase
        #phase = (phase/2./pi - floor(phase/2./pi))*2.*pi
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
            #print "1part = ", (gamma2*(IbetX2 + IbetY2 + prX*prX/LenPntrZ + prY*prY/LenPntrZ - phaseConstIn))
            #print "phasw = ", phase[600]/2./pi
            #print "delta = ", motion.z[-1] - motion.z[0]
            #print "phaseConst = ", phaseConst[0]
            #print gamma2, EreX[0]
            screen.arPhase = screen.arPhase + phase*(1 + 0*2*tmp*5e-7)
    return screen



def radiation_py(gamma, traj, screen, tmp):
    """
    screen format     screen->ReEx[ypoint*xpoint*je + xpoint*jy + jx] += EreX;
    """
    #start = time()
    motion = traj2motion(traj)
    #Z = motion.z
    #print "motion = ", time() - start


    #gamma = gamma

    size = len(motion.z)
    Nmotion = int((size + 1)/3)
    half_step = (motion.z[-1]-motion.z[0])/2./(Nmotion-1)
    #plot(motion.z, motion.x, "r.-")
    #show()
    n_end = len(motion.z)-2
    start3 = time()
    Xscr = np.linspace(screen.x_start, screen.x_start+screen.x_step*(screen.nx-1), num = screen.nx)
    Yscr = np.linspace(screen.y_start, screen.y_start+screen.y_step*(screen.ny-1), num = screen.ny)
    Yscr = Yscr.reshape((screen.ny, 1))
    Erad = np.linspace(screen.e_start, screen.e_start+screen.e_step*(screen.ne-1), num = screen.ne)
    Erad = Erad.reshape((screen.ne, 1))
    #print Xscr
    #print Yscr
    shape_array = [screen.ne, screen.ny, screen.nx]
    #print shape(Xscr), shape(Yscr), shape(Erad)
    if 1 in shape_array:
        #ind = shape_array.index(1)
        if screen.ny >1 and screen.ne>1:
            Yscr = Yscr.reshape((1, screen.ny))
        shape_array.remove(1)
        #print shape_array
        screen.arReEx = screen.arReEx.reshape(shape_array)
        screen.arImEx = screen.arImEx.reshape(shape_array)
        screen.arReEy = screen.arReEy.reshape(shape_array)
        screen.arImEy = screen.arImEy.reshape(shape_array)
        screen.arPhase = screen.arPhase.reshape(shape_array)
        #print screen.arPhase

        for n in range(Nmotion-1):
            #print "n = ", n
            screen = gintegrator(Xscr, Yscr, Erad, motion, screen, n, n_end, gamma, half_step, tmp)
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
                #print arPhase
        screen.arReEx = arReEx
        screen.arImEx = arImEx
        screen.arReEy = arReEy
        screen.arImEy = arImEy
        screen.arPhase = arPhase
    #print "reshape = ", time() - start3

    return 1


def calculate_radiation(lat, screen, beam, energy_loss=False, quantum_diff=False, accuracy=1):
    screen.update()
    b_current = beam.I*1000. # b_current - beam current must be in [mA], but beam.I in [A]
    energy = beam.E
    gamma = energy/m_e_GeV

    screen.nullify()
    U, E = track4rad(beam, lat, energy_loss=energy_loss, quantum_diff=quantum_diff, accuracy=accuracy)
    tmp = 0
    for u, e in zip(U, E):
        #print "Energy = ", e, e/m_e_GeV
        #plot(u[4::9], u[0::9],".-")
        radiation_py(e/m_e_GeV, u, screen, tmp)
        # screen.arPhase[0:10]/2./pi
        tmp += 1
    screen.distPhoton( gamma, current = b_current)
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
