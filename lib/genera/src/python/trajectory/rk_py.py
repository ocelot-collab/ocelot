__author__ = 'Sergey Tomin'
import numpy as np
import numexpr as ne


def fields(x,y,z, kx, ky, kz, B0):
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
    #Bx = ne.evaluate("k1*sin(kx*x)*sinhy*cosz")
    #By = ne.evaluate("B0*cosx*cosh(ky*y)*cosz")
    #Bz = ne.evaluate("k2*cosx*sinhy*sin(kz*z)")
    return Bx, By, Bz

def moments(bx, by, Bx, By, Bz, dzk):
    bx2 = bx*bx
    by2 = by*by
    bxy = bx*by
    #sq = 1 + (bx2 + by2)/2.
    sq = np.sqrt(1. + bx2 + by2)
    k = sq*dzk
    #sq = ne.evaluate("sqrt(1 + bx2 + by2)")
    mx = k*(by*Bz - By*(1.+bx2) + bxy*Bx)
    my = -k*(bx*Bz - Bx*(1.+by2) + bxy*By)
    #mx = ne.evaluate("sq*(by*Bz - By*(1+bx2) + bxy*Bx)*dzk")
    #my = ne.evaluate("-sq*(bx*Bz - Bx*(1+by2) + bxy*By)*dzk")
    return mx, my

def track_und_RK(y0, h,N, kz, kx, Kx, energy):
    gamma = energy*1957.
    c = 299792458
    m0 = 0.510998928*1e+6
    B0 = Kx*m0*kz/c
    ky = np.sqrt(kz*kz - kx*kx)
    charge = 1
    mass = 1 #in electron mass
    cmm = 299792458
    massElectron = 0.510998910e+6 #// rest mass of electron

    u = np.zeros((N*9,len(y0)/6))
    px = y0[1::6]
    py = y0[3::6]
    dz = h
    dGamma2 = 1. - 0.5/(gamma*gamma)
    pz = dGamma2 - (px*px + py*py)/2.
    k = charge*cmm/(massElectron*mass*gamma)
    u[0,:] = y0[0::6]
    u[1,:] = y0[1::6]
    u[2,:] = y0[2::6]
    u[3,:] = y0[3::6]
    u[4,:] = y0[4::6]
    u[5,:] = pz
    dzk = dz*k
    for i in range(N-1):
        X = u[i*9 + 0]
        Y = u[i*9 + 2]
        Z = u[i*9 + 4]

        bxconst = u[i*9 + 1]
        byconst = u[i*9 + 3]
        #bz = u[i*6 + 5]
        bx = bxconst
        by = byconst
        kx1 = bx*dz
        ky1 = by*dz
        Bx, By, Bz = fields(X, Y, Z, kx, ky, kz, B0)
        mx1, my1 = moments(bx, by, Bx, By, Bz, dzk)
        u[i*9 + 6] = Bx
        u[i*9 + 7] = By
        u[i*9 + 8] = Bz
        #K2
        bx = bxconst + mx1/2.
        by = byconst + my1/2.
        kx2 = bx*dz
        ky2 = by*dz
        Bx, By, Bz = fields(X + kx1/2., Y + ky1/2., Z + dz/2., kx, ky, kz, B0)
        mx2, my2 = moments(bx, by, Bx, By, Bz, dzk)
        # K3
        bx = bxconst + mx2/2.
        by = byconst + my2/2.
        kx3 = bx*dz
        ky3 = by*dz
        Bx, By, Bz = fields(X + kx2/2., Y + ky2/2., Z + dz/2., kx, ky, kz, B0)
        mx3, my3 = moments(bx, by, Bx, By, Bz, dzk)
        #K4
        Z_n = Z + dz
        bx = bxconst + mx3
        by = byconst + my3
        kx4 = bx*dz
        ky4 = by*dz
        Bx, By, Bz = fields(X + kx3, Y + ky3, Z_n, kx, ky, kz, B0)
        mx4, my4 = moments(bx, by, Bx, By, Bz, dzk)

        u[(i+1)*9 + 0] = X + 1/6.*(kx1 + 2.*(kx2 + kx3) + kx4)
        u[(i+1)*9 + 1] = bxconst + 1/6.*(mx1 + 2.*(mx2 + mx3) + mx4) #// conversion in mrad
        u[(i+1)*9 + 2] = Y + 1/6.*(ky1 + 2.*(ky2 + ky3) + ky4)
        u[(i+1)*9 + 3] = byconst + 1/6.*(my1 + 2.*(my2 + my3) + my4)
        u[(i+1)*9 + 4] = Z_n
        u[(i+1)*9 + 5] = dGamma2 - (u[(i+1)*9 + 1]*u[(i+1)*9 + 1] + u[(i+1)*9 + 3]*u[(i+1)*9 + 3])/2.

        #u[(i+1)*6 + 1] = betax
        #u[(i+1)*6 + 3] = betay
    u[(N-1)*9 + 6], u[(N-1)*9 + 7], u[(N-1)*9 + 8] = fields(u[(N-1)*9 + 0], u[(N-1)*9 + 2], u[(N-1)*9 + 4], kx, ky, kz, B0)
    return u

def rk_track_in_field(y0, h, N, energy, mag_field):
    gamma = energy*1957.
    #c = 299792458
    #m0 = 0.510998928*1e+6
    #B0 = Kx*m0*kz/c
    #ky = np.sqrt(kz*kz - kx*kx)
    charge = 1
    mass = 1 #in electron mass
    cmm = 299792458
    massElectron = 0.510998910e+6 #// rest mass of electron

    u = np.zeros((N*9,len(y0)/6))
    px = y0[1::6]
    py = y0[3::6]
    dz = h
    dGamma2 = 1. - 0.5/(gamma*gamma)
    pz = dGamma2 - (px*px + py*py)/2.
    k = charge*cmm/(massElectron*mass*gamma)
    u[0,:] = y0[0::6]
    u[1,:] = y0[1::6]
    u[2,:] = y0[2::6]
    u[3,:] = y0[3::6]
    u[4,:] = y0[4::6]
    u[5,:] = pz
    dzk = dz*k
    for i in range(N-1):
        X = u[i*9 + 0]
        Y = u[i*9 + 2]
        Z = u[i*9 + 4]

        bxconst = u[i*9 + 1]
        byconst = u[i*9 + 3]
        #bz = u[i*6 + 5]
        bx = bxconst
        by = byconst
        kx1 = bx*dz
        ky1 = by*dz
        Bx, By, Bz = mag_field(X, Y, Z)
        mx1, my1 = moments(bx, by, Bx, By, Bz, dzk)
        u[i*9 + 6] = Bx
        u[i*9 + 7] = By
        u[i*9 + 8] = Bz
        #K2
        bx = bxconst + mx1/2.
        by = byconst + my1/2.
        kx2 = bx*dz
        ky2 = by*dz
        Bx, By, Bz = mag_field(X + kx1/2., Y + ky1/2., Z + dz/2.)
        mx2, my2 = moments(bx, by, Bx, By, Bz, dzk)
        # K3
        bx = bxconst + mx2/2.
        by = byconst + my2/2.
        kx3 = bx*dz
        ky3 = by*dz
        Bx, By, Bz = mag_field(X + kx2/2., Y + ky2/2., Z + dz/2.)
        mx3, my3 = moments(bx, by, Bx, By, Bz, dzk)
        #K4
        Z_n = Z + dz
        bx = bxconst + mx3
        by = byconst + my3
        kx4 = bx*dz
        ky4 = by*dz
        Bx, By, Bz = mag_field(X + kx3, Y + ky3, Z_n)
        mx4, my4 = moments(bx, by, Bx, By, Bz, dzk)

        u[(i+1)*9 + 0] = X + 1/6.*(kx1 + 2.*(kx2 + kx3) + kx4)
        u[(i+1)*9 + 1] = bxconst + 1/6.*(mx1 + 2.*(mx2 + mx3) + mx4) #// conversion in mrad
        u[(i+1)*9 + 2] = Y + 1/6.*(ky1 + 2.*(ky2 + ky3) + ky4)
        u[(i+1)*9 + 3] = byconst + 1/6.*(my1 + 2.*(my2 + my3) + my4)
        u[(i+1)*9 + 4] = Z_n
        u[(i+1)*9 + 5] = dGamma2 - (u[(i+1)*9 + 1]*u[(i+1)*9 + 1] + u[(i+1)*9 + 3]*u[(i+1)*9 + 3])/2.

        #u[(i+1)*6 + 1] = betax
        #u[(i+1)*6 + 3] = betay
    u[(N-1)*9 + 6], u[(N-1)*9 + 7], u[(N-1)*9 + 8] = mag_field(u[(N-1)*9 + 0], u[(N-1)*9 + 2], u[(N-1)*9 + 4])
    return u