__author__ = 'Sergey Tomin'

import numpy as np
from ocelot.lib.genera.src.python.trajectory.spline_py import *
from pylab import *
from time import time
from scipy import interpolate
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
    """
    sqrt35 = 0.5*np.sqrt(3./5.)
    xnew = [X[0]]
    h = X[1] - X[0]
    xgaus = np.array([0.5-sqrt35, 0.5-sqrt35 + sqrt35, 0.5-sqrt35 + sqrt35 + sqrt35])*h
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
    #new_motion = Motion()

    motion.z = motion.z.flatten()

    Z = x2xgaus(motion.z)

    motion.x = bspline(motion.z, motion.x, Z)*1000.
    motion.y = bspline(motion.z, motion.y, Z)*1000.

    Ibx2 = integ_beta2(motion.z*1000., motion.bx)
    Iby2 = integ_beta2(motion.z*1000., motion.by)

    motion.XbetaI2 = bspline(motion.z, Ibx2, Z)
    motion.YbetaI2 = bspline(motion.z, Iby2, Z)

    motion.bx = bspline(motion.z, motion.bx, Z)
    motion.by = bspline(motion.z, motion.by, Z)


    motion.Bx = bspline(motion.z, motion.Bx, Z)
    motion.By = bspline(motion.z, motion.By, Z)
    motion.z = Z*1000.


    return motion

def gintegrator(Xscr, Yscr, Erad, motion, screen, n, n_end, gamma):
    Q = 0.5866740802042227#; // (mm*T)^-1
    hc = 1.239841874330e-3 # // mm
    k2q3 = 1.1547005383792517#;//  = 2./sqrt(3)
    gamma2 = gamma*gamma
    size = len(motion.z)
    Nmotion = (size + 1)/3
    h2 = (motion.z[-1]-motion.z[0])/2./(Nmotion-1)

    w = [0.5555555555555556*h2, 0.8888888888888889*h2, 0.5555555555555556*h2]
    LenPntrConst = screen.Distance - motion.z[0]#; // I have to pay attention to this
    for p in range(3):# // Gauss integration
        i = n*3 + p + 1
        radConstAdd = w[p]*Q*k2q3*(screen.Distance - screen.Zstart)
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
        #denominator = (1. + tx2 + ty2)*(1. + tx2 + ty2)
        radConst = radConstAdd/LenPntrZ/((1. + tx2 + ty2)*(1. + tx2 + ty2))
        radX = radConst*(By*(1. - tx2 + ty2) + Bx*tyx - 2.*tx/Q/LenPntrZ)#/*sigma*/
        radY =-radConst*(Bx*(1. + tx2 - ty2) + By*tyx + 2.*ty/Q/LenPntrZ)#;/*pi*/
        phaseConst = np.pi*Erad/(gamma2*hc)
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
            #j = n*3 + p + 2
            prX = Xscr - motion.x[-1] # //for pointer nx(z)
            prY = Yscr - motion.y[-1] # //for pointer ny(z)
            IbetX2 = motion.XbetaI2[-1]
            IbetY2 = motion.YbetaI2[-1]
            phase = phaseConst*(motion.z[-1] - motion.z[0]  + gamma2*(IbetX2 + IbetY2 + prX*prX/LenPntrZ + prY*prY/LenPntrZ - phaseConstIn))
            screen.arPhase += phase
    return screen



def radiation_py(gamma, traj, screen):
    """
    screen format     screen->ReEx[ypoint*xpoint*je + xpoint*jy + jx] += EreX;
    """
    #start = time()
    motion = traj2motion(traj)
    #Z = motion.z
    #print "motion = ", time() - start


    gamma = gamma

    size = len(motion.z)
    Nmotion = (size + 1)/3

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
        for n in range(Nmotion-1):
            screen = gintegrator(Xscr, Yscr, Erad, motion, screen, n, n_end, gamma)
        screen.arReEx = screen.arReEx.flatten()
        screen.arImEx = screen.arImEx.flatten()
        screen.arReEy = screen.arReEy.flatten()
        screen.arImEy = screen.arImEy.flatten()
        screen.arPhase = screen.arPhase.flatten()
    else:
        print "SR 3D calculation"
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
                screen = gintegrator(Xscr, Yscr, erad, motion, screen, n, n_end, gamma)

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
    print "reshape = ", time() - start3

    return 1



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

    x = np.linspace(0, 1, 4)
    xnew = x2xgaus(x)
    print x
    print xnew