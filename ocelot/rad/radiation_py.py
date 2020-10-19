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
import time
from ocelot.common.ocelog import *
import copy
from scipy.integrate import cumtrapz
import numbers

_logger = logging.getLogger(__name__)

try:
    import numba as nb

    nb_flag = True
except:
    _logger.info("radiation_py.py: module NUMBA is not installed. Install it to speed up calculation")
    nb_flag = False


class Motion:
    def __init__(self):
        self.x = []
        self.y = []
        self.z = []
        self.bx = []
        self.by = []
        self.bz = []
        self.Bx = []
        self.By = []
        self.Bz = []
        self.XbetaI2 = []
        self.YbetaI2 = []


class BeamTraject:
    """
    A class for storing and retrieving the coordinates of the n-th particle from the table with all trajectories.
    method: x(n=0) - array, horizontal coordinates of the n-th particle
    method: y(n=0) - array, vertical coordinates of the n-th particle
    method: xp(n=0) - array, x' = dx/dz coordinates of the n-th particle
    method: yp(n=0) - array, y' = dy/dz coordinates of the n-th particle
    method: z(n=0) - array, longitudinal coordinates in Cartesian coordinate system of the n-th particle
    method: s(n=0) - array,  longitudinal coordinates in moving coordinate system of the n-th particle
    """

    def __init__(self, beam_trajectories):
        self.U = beam_trajectories

    def n(self):
        return np.shape(self.U[0])[1]

    def check(self, n):
        if n > self.n() - 1:
            raise Exception('n > number of particles')

    def x(self, n=0):
        self.check(n)
        x_array = np.array([])
        for u in self.U:
            x_array = np.append(x_array, u[0::9, n])
        return x_array

    def y(self, n=0):
        self.check(n)
        y_array = np.array([])
        for u in self.U:
            y_array = np.append(y_array, u[2::9, n])
        return y_array

    def z(self, n=0):
        self.check(n)
        z_array = np.array([])
        for u in self.U:
            z_array = np.append(z_array, u[4::9, n])
        return z_array

    def xp(self, n=0):
        self.check(n)
        xp_array = np.array([])
        for u in self.U:
            xp_array = np.append(xp_array, u[1::9, n])
        return xp_array

    def yp(self, n=0):
        self.check(n)
        yp_array = np.array([])
        for u in self.U:
            yp_array = np.append(yp_array, u[3::9, n])
        return yp_array

    def p(self, n=0):
        self.check(n)
        p_array = np.array([])
        for u in self.U:
            p_array = np.append(p_array, u[5::9, n])
        return p_array

    def s(self, n=0):
        self.check(n)

        xp2 = self.xp(n) ** 2
        yp2 = self.yp(n) ** 2
        # zp = np.sqrt(1. / (1. + xp2 + yp2))
        s = cumtrapz(np.sqrt(1. + xp2 + yp2), self.z(n), initial=0)
        return s

    def p_array_end(self, p_array):

        s_fin = p_array.tau()

        for u in self.U:
            x1 = u[1::9, :]
            y1 = u[3::9, :]
            # dz = u[4 + 9::9, :] - u[4:-9:9, :]
            z = u[4::9, :]
            s_fin += np.trapz(np.sqrt(1 + x1 * x1 + y1 * y1), z, axis=0)
            # s_fin += np.sum(dz * np.sqrt(1 + x1 * x1 + y1 * y1), axis=0)

        N = int(np.shape(self.U[-1])[0] / 9)

        # ref_path is reference path of the particle with zero initial conditions
        # in sake of speed and simplicity we assume that ref_path is equal to path of the beam in average

        ref_path = np.mean(s_fin)

        p_array.rparticles[0, :] = self.U[-1][(N - 1) * 9 + 0, :]
        p_array.rparticles[1, :] = self.U[-1][(N - 1) * 9 + 1, :]
        p_array.rparticles[2, :] = self.U[-1][(N - 1) * 9 + 2, :]
        p_array.rparticles[3, :] = self.U[-1][(N - 1) * 9 + 3, :]
        p_array.rparticles[4, :] = ref_path - s_fin
        p_array.rparticles[5, :] = self.U[-1][(N - 1) * 9 + 5, :]


def bspline(x, y, x_new):
    tck = interpolate.splrep(x, y, s=0)
    ynew = interpolate.splev(x_new, tck, der=0)
    return ynew


def integ_beta2(x, y):
    A, B, C, D, Z = cspline_coef(x, y)
    b2 = 0.
    beta2 = [0.]
    for i in range(len(x) - 1):
        h = x[i + 1] - x[i]
        a = A[i]
        b = B[i]
        c = C[i]
        d = D[i]
        # print h, a,b,c,d
        b2 += h * (d * d + h * (c * d + h * (1. / 3. * (c * c + 2 * b * d) + h * (0.5 * (b * c + a * d) + h * (
                    0.2 * (b * b + 2 * a * c) + h * (1. / 3. * a * b + (a * a * h) / 7.))))))
        beta2.append(b2)
        # print beta2
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
    sqrt35 = 0.5 * np.sqrt(3. / 5.)
    xnew = [X[0]]
    h = X[1] - X[0]
    xgaus = np.array([0.5 - sqrt35, 0.5 - sqrt35 + sqrt35, 0.5 - sqrt35 + sqrt35 + sqrt35]) * h
    for x in X[:-1]:
        xnew = np.append(xnew, x + xgaus)
    xnew = np.append(xnew, X[-1])
    return xnew


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
    # new_motion = Motion()

    motion.z = motion.z.flatten()

    Z = x2xgaus(motion.z)

    motion.x = bspline(motion.z, motion.x, Z) * 1000.
    motion.y = bspline(motion.z, motion.y, Z) * 1000.
    # print "inegr = ", simps(motion.bx.flatten()**2, motion.z*1000)
    Ibx2 = integ_beta2(motion.z * 1000., motion.bx)
    Iby2 = integ_beta2(motion.z * 1000., motion.by)

    motion.XbetaI2 = bspline(motion.z, Ibx2, Z)
    motion.YbetaI2 = bspline(motion.z, Iby2, Z)

    motion.bx = bspline(motion.z, motion.bx, Z)
    motion.by = bspline(motion.z, motion.by, Z)

    motion.Bx = bspline(motion.z, motion.Bx, Z)
    motion.By = bspline(motion.z, motion.By, Z)
    motion.z = Z * 1000.

    return motion


def und_field_py(x, y, z, lperiod, Kx, nperiods=None):
    kx = 0.
    kz = 2 * pi / lperiod
    ky = np.sqrt(kz * kz + kx * kx)
    c = speed_of_light
    m0 = m_e_eV
    B0 = Kx * m0 * kz / c
    k1 = -B0 * kx / ky
    k2 = -B0 * kz / ky

    kx_x = kx * x
    ky_y = ky * y
    kz_z = kz * z

    cosz = np.cos(kz_z)

    if nperiods is not None:
        ph_shift = np.pi / 2.
        heaviside = lambda x: 0.5 * (np.sign(x) + 1)
        z_coef = (0.25 * heaviside(z) + 0.5 * heaviside(z - lperiod / 2.) + 0.25 * heaviside(z - lperiod)
                  - 0.25 * heaviside(z - (nperiods - 1) * lperiod) - 0.5 * heaviside(
                    z - (nperiods - 0.5) * lperiod)
                  - 0.25 * heaviside(z - nperiods * lperiod))
        cosz = np.cos(kz_z + ph_shift) * z_coef

    cosx = np.cos(kx_x)
    sinhy = np.sinh(ky_y)
    # cosz = np.cos(kz_z + ph_shift)*z_coef
    Bx = k1 * np.sin(kx_x) * sinhy * cosz  # // here kx is only real
    By = B0 * cosx * np.cosh(ky_y) * cosz
    Bz = k2 * cosx * sinhy * np.sin(kz_z)
    return (Bx, By, Bz)


und_field = und_field_py if not nb_flag else nb.jit(forceobj=False)(und_field_py)


def energy_loss_und(energy, Kx, lperiod, L, energy_loss=False):
    if energy_loss:
        k = 4. * pi * pi / 3. * ro_e / m_e_GeV
        U = k * energy ** 2 * Kx ** 2 * L / lperiod ** 2
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
    gamma = energy / m_e_GeV
    lambda_compt = 2.4263102389e-12  # m
    lambda_compt_r = lambda_compt / 2. / pi
    f = lambda K: 1.2 + 1. / (K + 1.33 * K * K + 0.4 * K ** 3)
    delta_Eq2 = 56. * pi ** 3 / 15. * lambda_compt_r * ro_e * gamma ** 4 / lperiod ** 3 * Kx ** 3 * f(Kx) * L
    sigma_Eq = np.sqrt(delta_Eq2 / (gamma * gamma))
    return sigma_Eq


def quantum_diffusion(energy, Kx, lperiod, L, quantum_diff=False):
    if quantum_diff:
        # gamma = energy/m_e_GeV
        # lambda_compt = 2.4263102389e-12 # h_eV_s/m_e_eV*speed_of_light
        # lambda_compt_r = lambda_compt/2./pi
        # f = lambda K: 1.2 + 1./(K + 1.33*K*K + 0.4*K**3)
        # delta_Eq2 = 56.*pi**3/15.*lambda_compt_r*ro_e*gamma**4/lperiod**3*Kx**3*f(Kx)*L
        sigma_Eq = sigma_gamma_quat(energy, Kx, lperiod, L)  # sqrt(delta_Eq2/(gamma*gamma))
        U = sigma_Eq * np.random.randn() * energy
    else:
        U = 0.
    return U


def field_map2field_func(z, By):
    tck = interpolate.splrep(z, By, k=3)
    func = lambda x, y, z: (0, interpolate.splev(z, tck, der=0), 0)
    return func


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
    Q = 0.5866740802042227  # speed_of_light/m_e_eV/1000  // e/mc = (mm*T)^-1
    hc = 1.239841874330e-3  # h_eV_s*speed_of_light*1000  // mm
    k2q3 = 1.1547005383792517  # ;//  = 2./sqrt(3)
    gamma2 = gamma * gamma
    w = [0.5555555555555556 * half_step, 0.8888888888888889 * half_step, 0.5555555555555556 * half_step]
    LenPntrConst = screen.Distance - motion.z[0]  # ; // I have to pay attention to this
    phaseConst = np.pi * Erad / (gamma2 * hc)
    for p in range(3):  # // Gauss integration
        i = n * 3 + p + 1
        # radConstAdd = w[p]*Q*k2q3*(screen.Distance - motion.z[0] - 0*screen.Zstart)
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

        prX = Xscr - XX  # //for pointer nx(z)
        prY = Yscr - YY  # //for pointer ny(z)
        nx = prX / LenPntrZ
        ny = prY / LenPntrZ
        tx = gamma * (nx - BetX)
        ty = gamma * (ny - BetY)
        tx2 = tx * tx
        ty2 = ty * ty
        tyx = 2. * tx * ty

        radConst = w[p] * Q * k2q3 * (screen.Distance) / LenPntrZ / ((1. + tx2 + ty2) * (1. + tx2 + ty2))

        radX = radConst * (By * (1. - tx2 + ty2) + Bx * tyx - 2. * tx / Q / LenPntrZ)  # /*sigma*/
        radY = -radConst * (Bx * (1. + tx2 - ty2) + By * tyx + 2. * ty / Q / LenPntrZ)  # ;/*pi*/

        prXconst = Xscr - motion.x[0]
        prYconst = Yscr - motion.y[0]
        phaseConstIn = (prXconst * prXconst + prYconst * prYconst) / LenPntrConst
        phaseConstCur = (prX * prX + prY * prY) / LenPntrZ
        # // string below is for case direct accumulation
        # //double phase = screen->Phase[ypoint*xpoint*je + xpoint*jy + jx] + faseConst*(ZZ - motion->Z[0]  + gamma2*(IbetX2 + IbetY2 + phaseConstCur - phaseConstIn));

        phase = phaseConst * (
                    ZZ - motion.z[0] + gamma2 * (IbetX2 + IbetY2 + phaseConstCur - phaseConstIn)) + screen.arPhase

        cosf = np.cos(phase)
        sinf = np.sin(phase)
        EreX = radX * cosf  # //(cosf *cos(fase0) - sinf*sin(fase0));
        EimX = radX * sinf  # //(sinf *cos(fase0) + cosf*sin(fase0));
        EreY = radY * cosf
        EimY = radY * sinf

        screen.arReEx += EreX
        screen.arImEx += EimX
        screen.arReEy += EreY
        screen.arImEy += EimY
        if i == n_end:  # //(n == 5000 && p == 2)
            LenPntrZ = screen.Distance - motion.z[-1]
            prX = Xscr - motion.x[-1]  # //for pointer nx(z)
            prY = Yscr - motion.y[-1]  # //for pointer ny(z)
            IbetX2 = motion.XbetaI2[-1]
            IbetY2 = motion.YbetaI2[-1]
            phase = phaseConst * (motion.z[-1] - motion.z[0] + gamma2 * (
                        IbetX2 + IbetY2 + prX * prX / LenPntrZ + prY * prY / LenPntrZ - phaseConstIn))
            screen.arPhase = screen.arPhase + phase
    return screen


def gintegrator_over_traj_py(Nmotion, Xscr, Yscr, Erad, n_end, gamma, half_step, Distance, x, y, z, bx, by,
                             XbetaI2, YbetaI2, Bx, By, arReEx, arImEx, arReEy, arImEy, arPhase):

    q = 0.5866740802042227  # speed_of_light/m_e_eV/1000  // e/mc = (mm*T)^-1
    hc = 1.239841874330e-3  # h_eV_s*speed_of_light*1000  // mm
    k2q3 = 1.1547005383792517  # 2./sqrt(3)
    gamma2 = gamma * gamma
    w = np.array([0.5555555555555556 * half_step, 0.8888888888888889 * half_step, 0.5555555555555556 * half_step])
    LenPntrConst = Distance - z[0]  # I have to pay attention to this
    phaseConst = np.pi * Erad / (gamma2 * hc)

    for n in range(Nmotion - 1):
        for p in range(3):  # // Gauss integration
            i = n * 3 + p + 1
            # radConstAdd = w[p]*Q*k2q3*(screen.Distance - motion.z[0] - 0*screen.Zstart)

            LenPntrZ = Distance - z[i]

            prX = Xscr - x[i]  # //for pointer nx(z)
            prY = Yscr - y[i]  # //for pointer ny(z)
            nx = prX / LenPntrZ
            ny = prY / LenPntrZ
            tx = gamma * (nx - bx[i])
            ty = gamma * (ny - by[i])
            tx2 = tx * tx
            ty2 = ty * ty
            tyx = 2. * tx * ty

            radConst = w[p] * q * k2q3 * Distance / LenPntrZ / ((1. + tx2 + ty2) * (1. + tx2 + ty2))

            radX = radConst * (By[i] * (1. - tx2 + ty2) + Bx[i] * tyx - 2. * tx / q / LenPntrZ)  # /*sigma*/
            radY = -radConst * (Bx[i] * (1. + tx2 - ty2) + By[i] * tyx + 2. * ty / q / LenPntrZ)  # ;/*pi*/

            prXconst = Xscr - x[0]
            prYconst = Yscr - y[0]
            phaseConstIn = (prXconst * prXconst + prYconst * prYconst) / LenPntrConst
            phaseConstCur = (prX * prX + prY * prY) / LenPntrZ
            # // string below is for case direct accumulation
            # //double phase = screen->Phase[ypoint*xpoint*je + xpoint*jy + jx] + faseConst*(ZZ - motion->Z[0]  + gamma2*(IbetX2 + IbetY2 + phaseConstCur - phaseConstIn));

            # here the constant phase shift was subtracted
            phase = phaseConst * (
                    z[i] - z[0] + gamma2 * (XbetaI2[i] + YbetaI2[i] + phaseConstCur - phaseConstIn)) + arPhase

            # phase = phaseConst * (z[i] - z[0] + gamma2 * (XbetaI2[i] + YbetaI2[i] + phaseConstCur)) + arPhase # + (LenPntrConst *2*np.pi * Erad/hc)%(2*np.pi)
            cosf = np.cos(phase)
            sinf = np.sin(phase)
            EreX = radX * cosf  # (cosf *cos(fase0) - sinf*sin(fase0));
            EimX = radX * sinf  # (sinf *cos(fase0) + cosf*sin(fase0));
            EreY = radY * cosf
            EimY = radY * sinf

            arReEx += EreX
            arImEx += EimX
            arReEy += EreY
            arImEy += EimY
            if i == n_end:
                LenPntrZ = Distance - z[-1]
                prX = Xscr - x[-1]  # //for pointer nx(z)
                prY = Yscr - y[-1]  # //for pointer ny(z)
                phase = phaseConst * (z[-1] - z[0] + gamma2 * (
                        XbetaI2[-1] + YbetaI2[-1] + prX * prX / LenPntrZ + prY * prY / LenPntrZ - phaseConstIn))
                # phase = phaseConst * (z[-1] - z[0] + gamma2 * (
                #        XbetaI2[-1] + YbetaI2[-1] + prX * prX / LenPntrZ + prY * prY / LenPntrZ )) # + (LenPntrConst *2*np.pi * Erad/hc)%(2*np.pi)
                arPhase += phase


gintegrator_over_traj = gintegrator_over_traj_py if not nb_flag else nb.jit(gintegrator_over_traj_py)


def wrap_gintegrator(Nmotion, Xscr, Yscr, Erad, motion, screen, n_end, gamma, half_step):
    Distance = screen.Distance
    x, y, z = motion.x, motion.y, motion.z
    bx, by = motion.bx, motion.by
    XbetaI2, YbetaI2 = motion.XbetaI2, motion.YbetaI2
    Bx, By = motion.Bx, motion.By

    gintegrator_over_traj(Nmotion, Xscr, Yscr, Erad, n_end, gamma, half_step, Distance,
                          x, y, z, bx, by, XbetaI2,
                          YbetaI2, Bx, By,
                          screen.arReEx, screen.arImEx, screen.arReEy, screen.arImEy, screen.arPhase)

    return screen


def radiation_py(gamma, traj, screen):
    """
    screen format     screen->ReEx[ypoint*xpoint*je + xpoint*jy + jx] += EreX;
    """

    motion = traj2motion(traj)

    size = len(motion.z)
    Nmotion = int((size + 1) / 3)
    half_step = (motion.z[-1] - motion.z[0]) / 2. / (Nmotion - 1)

    n_end = len(motion.z) - 2
    Xscr = np.linspace(screen.x_start, screen.x_start + screen.x_step * (screen.nx - 1), num=screen.nx)
    Yscr = np.linspace(screen.y_start, screen.y_start + screen.y_step * (screen.ny - 1), num=screen.ny)
    Yscr = Yscr.reshape((screen.ny, 1))
    Erad = np.linspace(screen.e_start, screen.e_start + screen.e_step * (screen.ne - 1), num=screen.ne)
    Erad = Erad.reshape((screen.ne, 1))

    shape_array = [screen.ne, screen.ny, screen.nx]
    if 1 in shape_array:
        if screen.ny > 1 and screen.ne > 1:
            Yscr = Yscr.reshape((1, screen.ny))
        shape_array.remove(1)
        screen.arReEx = screen.arReEx.reshape(shape_array)
        screen.arImEx = screen.arImEx.reshape(shape_array)
        screen.arReEy = screen.arReEy.reshape(shape_array)
        screen.arImEy = screen.arImEy.reshape(shape_array)
        screen.arPhase = screen.arPhase.reshape(shape_array)

        wrap_gintegrator(Nmotion, Xscr, Yscr, Erad, motion, screen, n_end, gamma, half_step)
        # print("phase", screen.arPhase )
        # for n in range(Nmotion-1):
        #     screen = gintegrator(Xscr, Yscr, Erad, motion, screen, n, n_end, gamma, half_step)
        screen.arReEx = screen.arReEx.flatten()
        screen.arImEx = screen.arImEx.flatten()
        screen.arReEy = screen.arReEy.flatten()
        screen.arImEy = screen.arImEy.flatten()
        screen.arPhase = screen.arPhase.flatten()
    else:
        print("SR 3D calculation")
        arReEx = np.array([])
        arImEx = np.array([])
        arReEy = np.array([])
        arImEy = np.array([])
        screen_segment = copy.deepcopy(screen)

        n_pl = screen.ny * screen.nx
        for i, erad in enumerate(Erad):
            arPhase = screen.arPhase[i * n_pl:(i + 1) * n_pl]
            screen_segment.arReEx = np.zeros((screen.ny, screen.nx))
            screen_segment.arImEx = np.zeros((screen.ny, screen.nx))
            screen_segment.arReEy = np.zeros((screen.ny, screen.nx))
            screen_segment.arImEy = np.zeros((screen.ny, screen.nx))

            screen_segment.arPhase = arPhase.reshape((screen.ny, screen.nx))

            # for n in range(Nmotion-1):
            #    screen_segment = gintegrator(Xscr, Yscr, erad, motion, screen_segment, n, n_end, gamma, half_step)
            wrap_gintegrator(Nmotion, Xscr, Yscr, erad, motion, screen_segment, n_end, gamma, half_step)

            arReEx = np.append(arReEx, screen_segment.arReEx.flatten())
            arImEx = np.append(arImEx, screen_segment.arImEx.flatten())
            arReEy = np.append(arReEy, screen_segment.arReEy.flatten())
            arImEy = np.append(arImEy, screen_segment.arImEy.flatten())
            screen.arPhase[i * n_pl:(i + 1) * n_pl] = screen_segment.arPhase.flatten()[:]
        screen.arReEx[:] += arReEx[:]
        screen.arImEx[:] += arImEx[:]
        screen.arReEy[:] += arReEy[:]
        screen.arImEy[:] += arImEy[:]
    return 1


def calculate_radiation(lat, screen, beam, energy_loss=False, quantum_diff=False, accuracy=1, end_poles=False):
    """
    Function to calculate radation from the electron beam.

    :param lat: MagneticLattice should include element Undulator
    :param screen: Screen class
    :param beam: Beam class, the radiation is calculated from one electron
    :param energy_loss: False, if True includes energy loss after each period
    :param quantum_diff: False, if True introduces random energy kick
    :param accuracy: 1, scale for trajectory points number
    :param end_poles: False, if True includes end poles with 1/4, -3/4, 1, ...
    :return:
    """

    screen.update()

    if beam.__class__ is Beam:
        p = Particle(x=beam.x, y=beam.y, px=beam.xp, py=beam.yp, E=beam.E)
        p_array = ParticleArray()
        p_array.list2array([p])

    # elif beam.__class__ is ParticleArray:
    #    b_current = beam.q_array[0] * 1000.
    #    p_array = beam

    else:
        raise TypeError("'beam' object must be Beam or ParticleArray class")

    if beam.I == 0:
        print("Beam charge or beam current is 0. Default current I=0.1 A is used")
        beam.I = 0.1  # A

    tau0 = np.copy(p_array.tau())
    p_array.tau()[:] = 0

    screen.nullify()
    start = time.time()
    U, E = track4rad_beam(p_array, lat, energy_loss=energy_loss, quantum_diff=quantum_diff, accuracy=accuracy,
                          end_poles=end_poles)
    # print("traj time exec:", time.time() - start)
    # plt.plot(U[0][4::9, :], U[0][::9, :])
    # plt.show()
    for i in range(p_array.n):
        # print("%i/%i" % (i, p_array.n))
        screen_copy = copy.deepcopy(screen)
        screen_copy.nullify()

        # wlengthes = h_eV_s*speed_of_light/screen_copy.Eph
        # screen_copy.arPhase[:] = tau0[i]/wlengthes*2*np.pi
        for u, e in zip(U, E):
            gamma = (1 + p_array.p()[i]) * e / m_e_GeV
            radiation_py(gamma, u[:, i], screen_copy)
        screen.arReEx += screen_copy.arReEx
        screen.arImEx += screen_copy.arImEx
        screen.arReEy += screen_copy.arReEy
        screen.arImEy += screen_copy.arImEy
        screen.arPhase += screen_copy.arPhase
    gamma_mean = (1 + np.mean(p_array.p())) * p_array.E / m_e_GeV
    screen.distPhoton(gamma_mean, current=beam.I)
    screen.Ef_electron = E[-1]
    screen.motion = U
    beam_traj = BeamTraject(beam_trajectories=U)
    beam_traj.p_array_end(p_array)
    screen.beam_traj = beam_traj

    # adding fast oscillating term to the phase
    motion = traj2motion(U[0][:, 0])
    screen.rebuild_efields(x0=motion.x[0], y0=motion.y[0], z0=motion.z[0])

    return screen


def coherent_radiation(lat, screen, p_array, energy_loss=False, quantum_diff=False, accuracy=1, end_poles=False,
                       verbose=True):
    """
    Function to calculate radiation from the electron beam.

    :param lat: MagneticLattice should include element Undulator
    :param screen: Screen class
    :param p_array: ParticleArray - the radiation is calculated for the each particles in the ParticleArray
                    and field components is summing up afterwards.
    :param energy_loss: False, if True includes energy loss after each period
    :param quantum_diff: False, if True introduces random energy kick
    :param accuracy: 1, scale for trajectory points number
    :param verbose: True, print progress
    :return:
    """

    screen.update()

    if p_array.__class__ is not ParticleArray:
        raise TypeError("'beam' object must be Beam or ParticleArray class")

    tau0 = np.copy(p_array.tau())
    p_array.tau()[:] = 0

    screen.nullify()
    U, E = track4rad_beam(p_array, lat, energy_loss=energy_loss, quantum_diff=quantum_diff, accuracy=accuracy,
                          end_poles=end_poles)
    # plt.plot(U[0][4::9, :], U[0][::9, :])
    # plt.show()
    for i in range(p_array.n):
        # print("%i/%i" % (i, p_array.n))
        screen_copy = copy.deepcopy(screen)
        screen_copy.nullify()

        wlengthes = h_eV_s * speed_of_light / screen_copy.Eph
        screen_copy.arPhase[:] = tau0[i] / wlengthes * 2 * np.pi
        for u, e in zip(U, E):
            gamma = (1 + p_array.p()[i]) * e / m_e_GeV

            radiation_py(gamma, u[:, i], screen_copy)
            # number of electrons in one macro particle
            n_e = p_array.q_array[i] / q_e

            screen.arReEx += screen_copy.arReEx * n_e * gamma
            screen.arImEx += screen_copy.arImEx * n_e * gamma
            screen.arReEy += screen_copy.arReEy * n_e * gamma
            screen.arImEy += screen_copy.arImEy * n_e * gamma
        if verbose:
            sys.stdout.write("\r" + "n: " + str(i) + " / " + str(p_array.n - 1))
            sys.stdout.flush()
    screen.coherent_photon_dist()
    motion = traj2motion(U[0][:, 0])
    screen.rebuild_efields(x0=motion.x[0], y0=motion.y[0], z0=motion.z[0])

    screen.Ef_electron = E[-1]
    screen.motion = U
    beam_traj = BeamTraject(beam_trajectories=U)
    beam_traj.p_array_end(p_array)
    screen.beam_traj = beam_traj

    return screen


def track4rad_beam(p_array, lat, energy_loss=False, quantum_diff=False, accuracy=1, end_poles=False):
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
    # Y0 = [beam.x, beam.xp, beam.y, beam.yp, 0, 0]
    # p = Particle(x=beam.x, px=beam.xp, y=beam.yp, py=beam.yp, E=beam.E)
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

                    N = int((lat_el.totalLen * 2000 + 150) * accuracy)
                    u = np.zeros((N * 9, np.shape(p_array.rparticles)[1]))
                    for i, z in enumerate(np.linspace(L, lat_el.totalLen + L, num=N)):
                        h = (lat_el.totalLen) / (N)
                        tracking_step(lat_el, p_array, h, navi)
                        u[i * 9 + 0, :] = p_array.rparticles[0]
                        u[i * 9 + 1, :] = p_array.rparticles[1]
                        u[i * 9 + 2, :] = p_array.rparticles[2]
                        u[i * 9 + 3, :] = p_array.rparticles[3]
                        u[i * 9 + 4, :] = p_array.rparticles[4] + z
                        u[i * 9 + 5, :] = p_array.rparticles[5]
                        # ui = [p.x, p.px, p.y, p.py, z, np.sqrt(1. - p.px*p.px - p.py *p.py), 0., 0., 0.]
                        # u.extend(ui)

                    U.append(u)
                    E.append(energy)
                L += lat_el.totalLen
            non_u = []

            U0 = energy_loss_und(energy, elem.Kx, elem.lperiod, elem.l, energy_loss)
            Uq = quantum_diffusion(energy, elem.Kx, elem.lperiod, elem.l, quantum_diff)
            U0 = U0 + Uq

            mag_length = elem.l
            if elem.mag_field is not None:
                mag_field = elem.mag_field
            else:
                if len(elem.field_map.z_arr) != 0:
                    # print("Field_map exist! Creating mag_field(x, y, z)")
                    unit_coef = 0.001 if elem.field_map.units == "mm" else 1
                    mag_length = elem.field_map.l * unit_coef
                    z_array = (elem.field_map.z_arr - elem.field_map.z_arr[0]) * unit_coef
                    mag_field = field_map2field_func(z=z_array, By=elem.field_map.By_arr)
                else:
                    # print("Standard undulator field")
                    if end_poles:
                        nperiods = elem.nperiods
                    else:
                        nperiods = None
                    mag_field = lambda x, y, z: und_field(x, y, z, elem.lperiod, elem.Kx, nperiods=nperiods)
            N = int((mag_length * 1500 + 100) * accuracy)
            if hasattr(elem, "npoints") and isinstance(elem.npoints, numbers.Number):
                N = int((elem.npoints + 100) * accuracy)
            u = rk_track_in_field(p_array.rparticles, mag_length, N, energy, mag_field)

            p_array.x()[:] = u[-9, :]
            p_array.px()[:] = u[-8, :]
            p_array.y()[:] = u[-7, :]
            p_array.py()[:] = u[-6, :]
            s = u[-5, 0]
            u[4::9] += L
            L += s
            U.append(u)
            E.append(energy)
        energy = energy - U0
    # for u in U:
    #     print("here", len(u[4::9, 0]))
    #     plt.plot(u[4::9, :], u[0::9, :])
    # plt.show()
    return U, E


if __name__ == "__main__":
    quantum_diffusion(17.5, 4., 0.04, 200., quantum_diff=True)
    x = np.linspace(0, 1, 4)
    xnew = x2xgaus(x)
    print(x)
    print(xnew)
