__author__ = 'Sergey Tomin'

import numpy as np
from scipy.integrate import odeint
from scipy import optimize
from ocelot.common.globals import m_e_GeV, m_e_eV, speed_of_light
from numpy.linalg import norm
from ocelot.common.ocelog import *

_logger = logging.getLogger(__name__)

try:
    import numba as nb
    nb_flag = True
except:
    _logger.info("high_order.py: module NUMBA is not installed. Install it to speed up calculation")
    nb_flag = False
    
__MAD__ = True

"""
differential equation:

t_nnn'' + kx**2*t_nnn = f_nnn

here index nnn means (1,2,3,4,5,6) = (x, x', y, y', s, delta)

# h = 1/r0 - curvature
# k1 = 1/Bro*dBy/dx - gradient
# k2 = 1/Bro*d^2(By)/dx^2 - sextupole component
# h1 - derivative of h
# h11 - derivative of h1
# cx1, sx1, cy1 ... - derivative of cx, sx, cy, ... and so on
# cx = cos(kx*s)
# sx = sin(kx*s)/kx
# dx = h/kx2*(1. - cx)
# cy = cos(ky*s)
# sy = sin(ky*s)/ky

Defenition:
Brown -> OCELOT
n = ky**2/h**2
beta = K2/h**3

# Green's function for X is Gx(t, tau) = 1/kx*sin(kx*(t - tau))

f111 = -(h**3 + K2 - 2*h*ky**2)*cx**2 + 1./2.*h*kx**4*sx*2 - kx**3*cx*sx*h'
f112 = - 2 * kx * (h**3 + K2 - 2 * h * ky**2) * cx * sx - h * kx**2 * cx * sx + (cx**2 - kx**2 * sx**2)*h'
f116 = (2*h^2 -ky^2)* cx + 2 *(h^3+K2-2 *h *ky^2) *cx* dx-h^2*kx^2 *sx^2 + (h *kx *cx *sx-dx* kx^2*sx)* h'
f122 = 1/2* h* cx^2 + h^3* (-1-K2/h^3+(2 ky^2)/h^2) *sx^2 + sx *sx* h'
f126 = h^2 *(2-ky^2/h^2) sx+h^2 *cx* sx + 2 h^3 (-1-K2/h^3+(2 ky^2)/h^2) *dx* sx + (cx dx +h sx^2) h'
f166 = -h + dx*h^2*(2-ky^2/h^2) + dx^2*h^3*(-1-K2/h^3+(2 ky^2)/h^2) + 1/2*h*(dx')^2 + dx*dx'* h'
f133 = -(1/2)* h*ky4*sy^2 - ky2*cy* sy*h' + 1/2*cy^2*(2*K2 - h*ky^2+h'') = |h' = 0, h'' = 0 | = -(1/2)* h*ky2 + cy**2*K2
f134 = h*ky2*cy*sy + (cy2 - sy2*ky2)*h' + cy*sy*(2*k2 - h*ky2 + h'')
f144 = -(1/2)*h*cy2 + cy*sy*h' + (sy2*(2*K2 - h*ky2 + h''))/2 =|h' = 0, h'' = 0 | =  -(1/2) *h +sy^2*K2

# Green's function for Y is Gy(t, tau) = 1/ky*sin(ky*(t - tau))

f313 = 2*(K2 - ky2*h)*cx*cy + h*kx2*ky2*sx*sy + (kx^2 *cy *sx-ky^2 *cx* sy) h'
f314 = -h kx^2 cy sx+2 h^3 (K2/h^3-ky^2/h^2) cx sy+(cx cy+kx^2 sx sy) h'
f323 = 2 h^3 (K2/h^3-ky^2/h^2) cy sx-h ky^2 cx sy+(-cx cy-ky^2 sx sy) h'
f324 = h cx cy+2 h^3 (K2/h^3-ky^2/h^2) sx sy + (cy sx-cx sy) h'
f336 = ky^2 cy+2 h^3 (K2/h^3-ky^2/h^2) dx cy-h^2 ky^2 sx sy-(h cy sx+ky^2 dx sy) h'
f346 = h^2 cy sx+ky^2 sy+2 h^3 (K2/h^3-ky^2/h^2) dx sy-(-dx cy+h sx sy) h'

Integration:

I111 = Gx * cx**2 = ((2. + cx)*dx)/(3*h)
I122 = Gx * sx**2 = dx**2/3./h**2
I112 = Gx * cx*sx = sx*dx/(3.*h)
I11  = Gx * cx    = s*sx/2.
I116 = Gx * cx*dx = h/kx2*(cx - cx**2) = h/kx2*(I11 - I111)
I12  = Gx * sx    = 1./2./kx2*(sx - s*cx/2.)
I126 = Gx * sx*dx = h/kx2*(sx - sx*cx) = h/kx2*(I12 - I112)
I10  = Gx         = dx/h
I16  = Gx * dx    = h/kx2*(dx/h - s*sx/2)
I166 = Gx * dx**2 = h2/kx4*(1 - 2*cx + cx**2) = h2/kx4*(I10 - 2*I11 + I111)
I144 = Gx * sy**2 = (sy2 - 2*dx/h)/(kx2 - 4*ky2)
I133 = Gx * cy**2 = dx/h + ky2*(2*dx/h - sy2)/(kx2 - 4*ky2)
I134 = Gx * cy*sy = (sy*cy - sx)/(kx2 - 4.*ky2)
I313 = Gy * cx*cy = (kx2*cy*dx/h - 2*ky2*sx*sy)/(kx2 - 4 *ky2)
I324 = Gy * sx*sy = (kx2*cy*dx/h - 2*ky2*sx*sy)/(kx2 - 4*ky2)
I314 = Gy * cx*sy = (2*cy*sx - (1 + cx)*sy)/(kx2 - 4*ky2)
I323 = Gy * sx*cy = ((1 - 2*ky2*dx/h)*sy - cy*sx)/ (kx2 - 4*ky2)
I33  = Gy * cy    = s*sy/2.
I336 = Gy * dx*cy = h/kx2*(I33 - I313)
I34  = Gy * sy    = (sy - s*cy)/(2*ky2)
I346 = Gy * dx*sy = h/kx2*(I34 - I314)
"""

def t_nnn_transport(L, h, k1, k2, energy=0):
    """
    here is used the following set of variables:
    x, dx/ds, y, dy/ds, delta_l, dp/p0
    :param L:
    :param angle:
    :param k1:
    :param k2:
    :return:
    """

    igamma2 = 0.
    if energy != 0:
        gamma = energy/m_e_GeV
        gamma2 = gamma*gamma
        igamma2 = 1./gamma2

    beta = np.sqrt(1. - igamma2)
    h2 = h*h
    h3 = h2*h
    kx2 = (k1 + h*h)
    ky2 = -k1
    kx4 = kx2*kx2
    ky4 = ky2*ky2
    kx = np.sqrt(kx2 + 0.j)
    ky = np.sqrt(ky2 + 0.j)
    cx = np.cos(kx*L).real
    sx = (np.sin(kx*L)/kx).real if kx != 0 else L
    cy = np.cos(ky*L).real

    sy = (np.sin(ky*L)/ky).real if ky != 0 else L

    sx2 = sx*sx
    sy2 = sy*sy
    L2 = L*L
    L3 = L2*L
    L4 = L3*L
    L5 = L4*L
    dx = h/kx2*(1. - cx) if kx != 0. else L*L*h/2.
    dx_h = (1. - cx)/kx2 if kx != 0. else L*L/2.

    # Integrals
    denom = kx2 - 4.*ky2
    I111 = 1./3.*(sx2 + dx_h)                          #  I111 = Gx * cx**2
    I122 = dx_h*dx_h/3.                                #  I122 = Gx * sx**2
    I112 = sx*dx_h/3.                                  #  I112 = Gx * cx*sx
    I11  = L*sx/2.                                     #  I11  = Gx * cx
    I10  = dx_h                                        #  I10  = Gx
    I33  = L*sy/2.                                     #  I33  = Gy * cy
    I34  = (sy - L*cy)/(2.*ky2) if ky !=0. else L3/6.  #  I34  = Gy * sy
    I211 = sx/3.*(1. + 2.*cx)
    I222 = 2.*dx_h*sx/3.
    I212 = 1./3.*(2*sx2 - dx_h)
    I21  = 1./2.*(L*cx + sx)
    I22  = I11
    I20  = sx
    I43  = 0.5*(L*cy + sy)
    I44  = I33

    #I5xx = h*Integrate(I1xx, dx)
    I512 = h*dx_h*dx_h/6
    I51  = L*dx/2.

    if kx != 0:
        I116 = h/kx2*(I11 - I111)                     #  I116 = Gx * cx*dx
        I12  = 0.5/kx2*(sx - L*cx)                    #  I12  = Gx * sx
        I126 = h/kx2*(I12 - I112)                     #  I126 = Gx * sx*dx
        I16  = h/kx2*(dx_h - L*sx/2.)                 #  I16  = Gx * dx
        I166 = h2/kx4*(I10 - 2*I11 + I111)            #  I166 = Gx * dx**2
        I216 = h/kx2*(I21 - I211)
        I226 = h/kx2*(I22 - I212)
        I26  = h /(2.*kx2)*(sx - L*cx)
        I266 = h2/kx4*(I20 - 2.*I21 + I211)

        I511 = h*(3.*L - 2.*sx - sx*cx)/(6.*kx2)
        I522 = h*(3.*L - 4*sx + sx*cx)/(6.*kx4)
        I516 = h/kx2*(I51 - I511)
        #I52  =  (dx - 4.*L*sx)/(2.*kx2)
        I52  =  (2.*dx - h*L*sx)/(2.*kx2)
        I526 = h/kx2*(I52 - I512)
        I50  = h*(L - sx)/kx2
        I566 = h2/kx4*(I50 - 2*I51 + I511)
        I56  = (h2*(L*(1. + cx) - 2.*sx))/(2.*kx4)

    else:
        I116 = h*L4/24.                               #  I116 = Gx * cx*dx
        I12  = L3/6.                                  #  I12  = Gx * sx
        I126 = h*L5/40.                               #  I126 = Gx * sx*dx
        I16  = h*L4/24.                               #  I16  = Gx * dx
        I166 = h2*L5*L/120.                           #  I166 = Gx * dx**2
        I216 = h*L3/6.
        I226 = h*L4/8.
        I26  = h*L3/6.
        I266 = h2*L5/20.

        I511 = h*L3/6.
        I522 = h*L5/60.
        I516 = h2*L5/120.
        I52  = h*L4/24.
        I526 = h2*L5*L/240.
        I50  = h*L3/6.
        I566 = h2*h*L5*L2/840.
        I56  = h2*L5/120.

    if kx != 0 and ky != 0:
        I144 = (sy2 - 2.*dx_h)/denom                         #  I144 = Gx * sy**2
        I133 = dx_h - ky2*(sy2 - 2.*dx_h)/denom              #  I133 = Gx * cy**2
        I134 = (sy*cy - sx)/denom                            #  I134 = Gx * cy*sy
        I313 = (kx2*cy*dx_h - 2.*ky2*sx*sy)/denom            #  I313 = Gy * cx*cy
        I324 = (2.*cy*dx_h - sx*sy)/denom                    #  I324 = Gy * sx*sy
        I314 = (2.*cy*sx - (1. + cx)*sy)/denom               #  I314 = Gy * cx*sy
        I323 = (sy - cy*sx - 2.*ky2*sy*dx_h)/denom           #  I323 = Gy * sx*cy = (2*ky2/kx2*(1 + cx)*sy - cy*sx)/denom + sy/kx2
        #derivative of Integrals
        I244 = 2.*(cy*sy - sx)/denom
        I233 = sx - 2.*ky2*(cy*sy - sx)/denom
        I234 = (kx2*dx_h - 2.*ky2*sy2)/denom
        I413 = ((kx2 - 2.*ky2)*cy*sx - ky2*sy*(1. + cx))/denom
        I424 = (cy*sx - cx*sy - 2.*ky2*sy*dx_h)/denom
        I414 = ((kx2 - 2.*ky2)*sx*sy - (1. - cx)*cy)/denom
        I423 = (cy*dx_h*(kx2 - 2*ky2) - ky2*sx*sy)/denom      #  I423 = I323' = ((2.*ky2)/kx2*(1 + cx)*cy - cx*cy - ky2*sx*sy)/denom + cy/kx2

    elif kx != 0 and ky == 0:
        I323 = (L - sx)/kx2
        I324 = 2.*(1. - cx)/kx4 - L*sx/kx2
        I314 = (2.*sx - L*(1. + cx))/kx2
        I313 = (1. - cx)/kx2
        I144 = (-2. + kx2*L2 + 2.*cx)/kx4
        I133 = (1. - cx)/kx2
        I134 = (L - sx)/kx2
        # derivative of Integrals
        I423 = (1. - cx)/kx2
        I424 = (sx - L*cx)/kx2
        #I414 = 2*cx/kx2 + sx
        I414 = (cx - 1.)/kx2 + L*sx
        I413 = sx
        I244 = (2.*L - 2.*sx)/kx2
        I233 = sx
        I234 = (1. - cx)/kx2
    else:
        I144 = L4/12.                                          #  I144 = Gx * sy**2
        I133 = L2/2.                                           #  I133 = Gx * cy**2
        I134 = L3/6.                                           #  I134 = Gx * cy*sy
        I313 = L2/2.                                           #  I313 = Gy * cx*cy
        I324 = L4/12.                                          #  I324 = Gy * sx*sy
        I314 = L3/6.                                           #  I314 = Gy * cx*sy
        I323 = L3/6.                                           #  I323 = Gy * sx*cy
        I244 = L3/3.
        I233 = L
        I234 = L2/2.
        I413 = L
        I424 = L3/3.
        I414 = L2/2.
        I423 = L2/2.

    # print I323, L, kx, ky
    if kx == 0 and ky != 0:
        I336 = (h*L*(3.*L*cy + (2.*ky2*L2 - 3.)*sy))/(24.*ky2)
        I346 = (h*((3. - 2.*ky2*L2)*L*cy + 3.*(ky2*L2 - 1.)*sy))/(24.*ky4)
        I436 = I346
        I446 = (h*L*(-3.*L*cy + (3. + 2.*ky2*L2)*sy))/(24.*ky2)

        I533 = (h*(3.*L + 2.*ky2*L3 - 3.*sy*cy))/(24.*ky2)
        I534 = (h*(L2 - sy2))/(8.*ky2)
        I544 = (h*(-3.*L + 2.*ky2*L3 + 3.*sy*cy))/(24.*ky4)

    elif kx == 0 and ky == 0:
        I336 = (h*L4)/24.
        I346 = (h*L5)/40.
        I436 = (h*L3)/6.
        I446 = (h*L4)/8.

        I533 = h*L3/6.
        I534 = h*L4/24.
        I544 = h*L5/60.

    else:
        I336 = h/kx2*(I33 - I313)                                  #  I336 = Gy * dx*cy
        I346 = h/kx2*(I34 - I314)                                  #  I346 = Gy * dx*sy
        I436 = h/kx2*(I43 - I413)
        I446 = h/kx2*(I44 - I414)

        I533 = (h*(denom*L - 2.*(denom + 2.*ky2)*sx + kx2*cy*sy))/(2.*denom*kx2)
        I534 = (h*sy2 - 2*dx)/(2*denom)
        I544 = (sy2 - 2*dx_h)/denom

    K2 = k2/2.
    coef1 = 2.*ky2*h - h3 - K2
    coef3 = 2.*h2 - ky2

    t111 =    coef1*I111 + h*kx4*I122/2.
    t112 = 2.*coef1*I112 - h*kx2*I112
    t116 = 2.*coef1*I116 + coef3*I11 - h2*kx2*I122
    t122 =    coef1*I122 + 0.5*h*I111
    t126 = 2.*coef1*I126 + coef3*I12 + h2*I112
    t166 =    coef1*I166 + coef3*I16 + 0.5*h3*I122 - h*I10
    t133 =       K2*I133 - ky2*h*I10/2.
    t134 =    2.*K2*I134
    t144 =       K2*I144 - h*I10/2.

    t211 =    coef1*I211 + h*kx4*I222/2.
    t212 = 2.*coef1*I212 - h*kx2*I212
    t216 = 2.*coef1*I216 + coef3*I21 - h2*kx2*I222
    #print(I216, I21, I222, t216)
    t222 =    coef1*I222 + 0.5*h*I211
    t226 = 2.*coef1*I226 + coef3*I22 + h2*I212
    t266 =    coef1*I266 + coef3*I26 + 0.5*h3*I222 - h*I20
    t233 =       K2*I233 - ky2*h*I20/2.
    t234 =    2.*K2*I234
    t244 =       K2*I244 - h*I20/2.

    coef2 = 2*(K2 - ky2*h)

    t313 = coef2*I313 + h*kx2*ky2*I324
    t314 = coef2*I314 - h*kx2*I323
    t323 = coef2*I323 - h*ky2*I314
    t324 = coef2*I324 + h*I313
    t336 = coef2*I336 + ky2*I33 - h2*ky2*I324
    t346 = coef2*I346 + h2*I323 + ky2*I34
    #print "I323 = ", I323, "L, h, k1, k2 = ", L, h, k1, k2, "  t314 = ", t314
    t413 = coef2*I413 + h*kx2*ky2*I424
    t414 = coef2*I414 - h*kx2*I423
    t423 = coef2*I423 - h*ky2*I414
    t424 = coef2*I424 + h*I413
    t436 = coef2*I436 - h2*ky2*I424 + ky2*I43
    t446 = coef2*I446 + h2*I423 + ky2*I44
    #print "I414=", I414, "  I423 = ", I423,  "   t414 = ", t414, L, kx, ky
    # Coordinates transformation from Curvilinear to a Cartesian
    cx_1 = -kx2*sx
    sx_1 = cx
    cy_1 = -ky2*sy
    sy_1 = cy
    dx_1 = h*sx
    T = np.zeros((6, 6, 6))
    T[0, 0, 0] = t111
    T[0, 0, 1] = t112 + h*sx
    T[0, 0, 5] = t116
    T[0, 1, 1] = t122
    T[0, 1, 5] = t126/beta
    T[0, 5, 5] = t166
    T[0, 2, 2] = t133
    T[0, 2, 3] = t134
    T[0, 3, 3] = t144

    T[1, 0, 0] = t211 - h*cx*cx_1
    T[1, 0, 1] = t212 + h*sx_1 - h*(sx*cx_1 + cx*sx_1)
    T[1, 0, 5] = t216 - h*(dx*cx_1 + cx*dx_1)
    T[1, 1, 1] = t222 - h*sx*sx_1
    T[1, 1, 5] = t226 - h*(sx*dx_1 + dx*sx_1)
    T[1, 5, 5] = t266 - dx*h*dx_1

    T[1, 2, 2] = t233
    T[1, 2, 3] = t234
    T[1, 3, 3] = t244

    T[2, 0, 2] = t313
    T[2, 0, 3] = t314 + h*sy
    T[2, 1, 2] = t323
    T[2, 1, 3] = t324
    T[2, 2, 5] = t336
    T[2, 3, 5] = t346/beta

    T[3, 0, 2] = t413 - h*cx*cy_1
    T[3, 0, 3] = t414 + (1 - cx)*h*sy_1
    T[3, 1, 2] = t423 - h*sx*cy_1
    T[3, 1, 3] = t424 - h*sx*sy_1
    T[3, 2, 5] = t436 - h*dx*cy_1
    T[3, 3, 5] = t446 - h*dx*sy_1

    """
    Path length difference
    linear = cx*h*x0 + h*sx*x0' + dx*h*dp;
    nonlinear = (h*T111 + 1/2*(cx')^2)*x0^2 + (h*T112 + cx'*sx')*x0*x0' + (h*T116 + cx'*dx')*x0*dp
                +(h*T122 + 1/2*(sx')^2)*(x0')^2 + (h*T126 + dx'*sx')*dp*x0' + (h*T166 + 1/2*(dx')^2)*dp^2
                +(h*T133 + 1/2*(cy')^2)*y0^2 + (h*T134 + cy'*sy')*y0*y0' + (h*T144 + 1/2*(sy')^2)*(y0')^2
    dl = Integrate(linear*ds + nonlinear*ds)
    """

    t511 =    coef1*I511 + h*kx4*I522/2.
    t512 = 2.*coef1*I512 - h*kx2*I512
    t516 = 2.*coef1*I516 + coef3*I51 - h2*kx2*I522
    t522 =    coef1*I522 + 0.5*h*I511
    t526 = 2.*coef1*I526 + coef3*I52 + h2*I512
    t566 =    coef1*I566 + coef3*I56 + 0.5*h3*I522 - h*I50
    t533 =       K2*I533 - ky2*h*I50/2.
    t534 =    2.*K2*I534
    t544 =       K2*I544 - h*I50/2.
    i566 = h2*(L - sx*cx)/(4.*kx2) if kx != 0 else h2*L3/6.

    T511 = t511 + 1/4.*kx2*(L - cx*sx)

    T512 = t512 - (1/2.)*kx2*sx2 + h*dx
    T516 = t516 + h*(sx*cx - L)/2.
    T522 = t522 + (L + sx*cx)/4.
    T526 = t526 + h*sx2/2.
    T566 = t566 + i566
    T533 = t533 + 1/4.*ky2*(L - sy*cy )
    T534 = t534 - 1/2.*ky2*sy2
    T544 = t544 + (L + sy*cy)/4.

    #beta = 1.
    #if energy != 0:
    #    gamma = energy/m_e_GeV
    #    gamma2 = gamma*gamma
    #    beta = np.sqrt(1. - 1./gamma2)
        #T[4, 5, 5] = 1.5*L/(beta*beta*gamma2)

    T[4, 0, 0] = T511/beta
    T[4, 0, 1] = (T512 + h*dx)/beta
    T[4, 0, 5] = T516/beta
    T[4, 1, 1] = T522/beta
    T[4, 1, 5] = T526/beta
    T[4, 5, 5] = T566/beta + 1.5*L/(beta*beta*beta)*igamma2
    T[4, 2, 2] = T533/beta
    T[4, 2, 3] = T534/beta
    T[4, 3, 3] = T544/beta


    """
    #print("T511 = ", T511)
    #print("T512 = ", T512)
    #print("T516 = ", T516)
    #print("T522 = ", T522)
    #print("T526 = ", T526)
    #print("T566 = ", T566)
    #print("T533 = ", T533)
    #print("T534 = ", T534)
    #print("T544 = ", T544)



    print("T111 = ", T[1-1, 1-1, 1-1])
    print("T112 = ", T[1-1, 1-1, 2-1])
    print("T116 = ", T[1-1, 1-1, 6-1])
    print("T122 = ", T[1-1, 2-1, 2-1])
    print("T126 = ", T[1-1, 2-1, 6-1])
    print("T166 = ", T[1-1, 6-1, 6-1])
    print("T133 = ", T[1-1, 3-1, 3-1])
    print("T134 = ", T[1-1, 3-1, 4-1])
    print("T144 = ", T[1-1, 4-1, 4-1])
    print("T211 = ", T[2-1, 1-1, 1-1])
    print("T212 = ", T[2-1, 1-1, 2-1])
    print("T216 = ", T[2-1, 1-1, 6-1])
    print("T222 = ", T[2-1, 2-1, 2-1])
    print("T226 = ", T[2-1, 2-1, 6-1])
    print("T266 = ", T[2-1, 6-1, 6-1])
    print("T233 = ", T[2-1, 3-1, 3-1])
    print("T234 = ", T[2-1, 3-1, 4-1])
    print("T244 = ", T[2-1, 4-1, 4-1])
    print("T313 = ", T[3-1, 1-1, 3-1])
    print("T314 = ", T[3-1, 1-1, 4-1])
    print("T323 = ", T[3-1, 2-1, 3-1])
    print("T324 = ", T[3-1, 2-1, 4-1])
    print("T336 = ", T[3-1, 3-1, 6-1])
    print("T346 = ", T[3-1, 4-1, 6-1])
    print("T413 = ", T[4-1, 1-1, 3-1])
    print("T414 = ", T[4-1, 1-1, 4-1])
    print("T423 = ", T[4-1, 2-1, 3-1])
    print("T424 = ", T[4-1, 2-1, 4-1])
    print("T436 = ", T[4-1, 3-1, 6-1])
    print("T446 = ", T[4-1, 4-1, 6-1])
    """
    return T


def t_nnn_mad(L, h, k1, k2, energy=0):
    """
    Universal function to return second order matrix.
    Here are used the following set of variables: x, px/p0, y, py/p0, tau, dE/(p0c).
    See https://nbviewer.jupyter.org/github/ocelot-collab/ocelot/blob/master/demos/ipython_tutorials/2_tracking.ipynb

    :param L: length in [m]
    :param h: curvature (1/r) of the bend in [1/m]
    :param k1: strength of quadrupole component in [1/m^2],
    :param k2: strength of sextupole component in [1/m^3],
    :param energy: energy in [GeV]
    :return: second order matrix, numpy array with shape = (6, 6, 6)
    """

    igamma2 = 0.
    if energy != 0:
        gamma = energy/m_e_GeV
        gamma2 = gamma*gamma
        igamma2 = 1./gamma2

    beta = np.sqrt(1. - igamma2)
    beta2 = beta * beta
    beta3 = beta2*beta
    h2 = h*h
    h3 = h2*h
    kx2 = (k1 + h*h)
    ky2 = -k1
    kx4 = kx2*kx2
    ky4 = ky2*ky2
    kx = np.sqrt(kx2 + 0.j)
    ky = np.sqrt(-k1 + 0.j)
    cx = np.cos(kx*L).real
    sx = (np.sin(kx*L)/kx).real if kx != 0 else L
    cy = np.cos(ky*L).real

    sy = (np.sin(ky*L)/ky).real if ky != 0 else L

    sx2 = sx*sx
    sy2 = sy*sy
    L2 = L*L
    L3 = L2*L
    L4 = L3*L
    L5 = L4*L
    dx = (1. - cx)/kx2 if kx != 0. else L*L/2.
    dx2 = dx * dx
    K1 = k1
    K2 = k2

    d2y = 0.5 * sy * sy
    s2y = sy * cy
    c2x = np.cos(2 * kx * L).real
    c2y = np.cos(2 * ky * L).real
    fx = (L - sx) / kx2 if kx2 != 0 else L3 / 6.
    f2y = (L - s2y) / ky2 if ky2 != 0 else L3 / 6
    J1 = (L - sx) / kx2 if kx2 != 0 else L3 / 6.
    J2 = (3. * L - 4. * sx + sx * cx) / (2 * kx4) if kx != 0 else L ** 5 / 20.
    J3 = (15 * L - 22.5 * sx + 9 * sx * cx - 1.5 * sx * cx * cx + kx2 * sx2 * sx) / (
                6 * kx4 * kx2) if kx != 0 else L ** 7 / 56.
    #J3 = (15 * L - 22 * sx + 9 * sx * cx - 2 * sx * cx * cx ) / (
    #        6 * kx4 * kx2) if kx != 0 else L ** 7 / 56.
    J_denom = kx2 - 4 * ky2
    Jc = (c2y - cx) / J_denom if J_denom != 0 else 0.5 * L2
    Js = (cy * sy - sx) / J_denom if J_denom != 0 else L3 / 6
    Jd = (d2y - dx) / J_denom if J_denom != 0 else L4 / 24
    Jf = (f2y - fx) / J_denom if J_denom != 0 else L5 / 120

    khk = K2 + 2*h*K1
    T = np.zeros((6,6,6))

    T111 = -1/6*khk*(sx2 + dx) - 0.5*h*kx2*sx2
    T112 = -1/6*khk*sx*dx + 0.5*h*sx*cx
    T122 = -1 / 6 * khk * dx2 + 0.5 * h * dx * cx
    T116 = -h / 12/beta * khk * (3*sx*J1 - dx2) + 0.5 * h2/beta * sx2 + 0.25/beta*K1*L*sx
    T126 = -h / 12 / beta * khk * (
                sx*dx2 - 2*cx*J2) + 0.25 * h2 / beta * (sx*dx + cx*J1) - 0.25 / beta * (sx + L*cx)
    T166 = -h2 / 6 / beta2 * khk * (
                dx2*dx - 2*sx*J2) + 0.5 * h3 / beta2 * sx*J1 - 0.5*h / beta2 * L*sx  - 0.5*h/(beta2)*igamma2*dx
    T133 = K1*K2*Jd + 0.5*(K2 + h * K1)*dx
    T134 = 0.5*K2*Js
    T144 = K2*Jd - 0.5*h*dx

    T211 = -1/6*khk*sx*(1 + 2*cx)
    T212 = -1 / 6 * khk * dx * (1 + 2 * cx)
    T222 = -1 / 3 * khk * sx*dx  - 0.5*h*sx
    T216 = -h / 12/beta * khk * (3*cx * J1 + sx*dx) - 0.25/beta * K1 * (sx - L*cx)
    T226 = -h / 12/beta * khk * (3*sx * J1 + dx2) + 0.25/beta * K1 *L*sx
    T266 = -h2 / 6 / beta2 * khk * (sx * dx2 - 2*cx*J2) - 0.5 *h/ beta2 * K1 * (cx*J1 - sx*dx) - 0.5*h/beta2*igamma2*sx
    T233 = K1*K2*Js + 0.5*(K2 + h*K1)*sx
    T234 = 0.5*K2*Jc
    T244 = K2 * Js - 0.5*h*sx

    T313 = 0.5*K2*(cy*Jc - 2*K1*sy*Js) + 0.5*h*K1*sx*sy
    T314 = 0.5 * K2 * (sy * Jc - 2* cy * Js) + 0.5 * h * sx * cy
    T323 = 0.5 * K2 * (cy * Js - 2 *K1* sy * Jd) + 0.5 * h *K1 * dx * sy
    T324 = 0.5 * K2 * (sy * Js - 2 * cy * Jd) + 0.5 * h * dx * cy
    T336 = 0.5 * h/beta * K2 * (cy * Jd - 2 * K1 * sy * Jf) + 0.5 * h2/beta * K1 * J1 * sy - 0.25/beta*K1*L*sy
    T346 = 0.5 * h / beta * K2 * (sy * Jd - 2 * cy * Jf) + 0.5 * h2 / beta * J1 * cy - 0.25 / beta * (sy + L * cy)

    T413 = 0.5*K1*K2*(2*cy*Js - sy*Jc) + 0.5*(K2 + h*K1)*sx*cy
    T414 = 0.5 * K2 * (2 * K1*sy * Js - cy * Jc) + 0.5 * (K2 + h * K1) * sx * sy
    T423 = 0.5*K1*K2*(2*cy*Jd - sy*Js) + 0.5*(K2 + h*K1)*dx*cy
    T424 = 0.5 * K2 * (2 * K1 * sy * Jd - cy * Js) + 0.5 * (K2 + h * K1) * dx * sy
    T436 = 0.5 * h/beta * K1 * K2 * (2 * cy * Jf - sy * Jd) + 0.5 *h/beta* (K2 + h * K1) * J1 * cy + 0.25/beta*K1*(sy - L*cy)
    T446 = 0.5 *h/beta * K2 * (2 * K1 * sy * Jf - cy * Jd) + 0.5 *h/beta* (K2 + h * K1) * J1 * sy - 0.25/beta*K1*L*sy

    T511 = h/12/beta*khk*(sx*dx + 3*J1) - 0.25/beta*K1*(L - sx*cx)
    T512 = h / 12 / beta * khk * dx2 + 0.25 / beta * K1 * sx2
    T522 = h/6/beta*khk*J2 - 0.5/beta*sx - 0.25/beta*K1*(J1 - sx*dx)

    T516 = h2/12/beta2*khk*(3*dx*J1 - 4*J2) + 0.25*h/beta2*K1*J1*(1 + cx) + 0.5*h/beta2*igamma2*sx

    T526 = h2 / 12 / beta2 * khk * (dx * dx2 - 2 *sx * J2) + 0.25 * h / beta2 * K1*sx * J1 + 0.5 * h / beta2 * igamma2 * dx
    T566 = h3 / 6 / beta3 * khk * (3*J3 - 2 *dx * J2) + h2/6 / beta3 * K1*(sx*dx2 - J2*(1 + 2*cx)) + 1.5 / beta3 * igamma2 * (h2*J1 - L)
    T533 = - h/beta*K1*K2*Jf - 0.5*h/beta*(K2 + h*K1)*J1 + 0.25/beta*K1*(L - cy*sy)
    T534 = - 0.5*h/beta*K2*Jd - 0.25/beta*K1*sy2
    T544 = -h/beta * K2*Jf + 0.5*h2/beta*J1  - 0.25/beta*(L + cy*sy)

    T[0, 0, 0] = T111
    T[0, 0, 1] = T112*2
    T[0, 1, 1] = T122
    T[0, 0, 5] = T116*2
    T[0, 1, 5] = T126*2
    T[0, 5, 5] = T166
    T[0, 2, 2] = T133
    T[0, 2, 3] = T134*2
    T[0, 3, 3] = T144

    T[1, 0, 0] = T211
    T[1, 0, 1] = T212*2
    T[1, 1, 1] = T222
    T[1, 0, 5] = T216*2
    T[1, 1, 5] = T226*2
    T[1, 5, 5] = T266
    T[1, 2, 2] = T233
    T[1, 2, 3] = T234*2
    T[1, 3, 3] = T244

    T[2, 0, 2] = T313*2
    T[2, 0, 3] = T314*2
    T[2, 1, 2] = T323*2
    T[2, 1, 3] = T324*2
    T[2, 2, 5] = T336*2
    T[2, 3, 5] = T346*2


    T[3, 0, 2] = T413*2
    T[3, 0, 3] = T414*2
    T[3, 1, 2] = T423*2
    T[3, 1, 3] = T424*2
    T[3, 2, 5] = T436*2
    T[3, 3, 5] = T446*2

    T[4, 0, 0] = -T511
    T[4, 0, 1] = -T512*2
    T[4, 1, 1] = -T522
    T[4, 0, 5] = -T516*2
    T[4, 1, 5] = -T526*2
    T[4, 5, 5] = -T566
    T[4, 2, 2] = -T533
    T[4, 2, 3] = -T534*2
    T[4, 3, 3] = -T544
    return T


def t_nnn(L, h, k1, k2, energy=0):
    if __MAD__:
        T = t_nnn_mad(L, h, k1, k2, energy)
    else:
        T = t_nnn_transport(L, h, k1, k2, energy)
    return T


def fringe_ent(h, k1, e, h_pole=0., gap=0., fint=0.):
    """
    Calculate first and second order matrices for edge focusing element.
    Note: the first order (R) is shown here in the sake of completeness but it is duplicated
            in r_matrix.py module and used from there in optics.py

    :param h: curvature (1/r) of the bend in [1/m]
    :param k1: strength of quadrupole component of the dipole in [1/m^2],
    :param e: the angle of inclination of the entrance face [rad],
    :param h_pole: the curvature (1/r) of the entrance face
    :param gap: the magnet gap [m], NOTE in MAD and ELEGANT: HGAP = gap/2
    :param fint: fringe field integral
    :return: R, T. First and second order matrices (6, 6), (6, 6, 6)
    """
    sec_e = 1./np.cos(e)
    sec_e2 = sec_e*sec_e
    sec_e3 = sec_e2*sec_e
    tan_e = np.tan(e)
    tan_e2 = tan_e*tan_e
    phi = fint*h*gap*sec_e*(1. + np.sin(e)**2)
    R = np.eye(6)
    R[1, 0] = h*tan_e
    R[3, 2] = -h*np.tan(e - phi)

    T = np.zeros((6,6,6))
    T[0, 0, 0] = -h/2.*tan_e2
    T[0, 2, 2] = h/2.*sec_e2
    T[1, 0, 0] = h/2.*h_pole*sec_e3 + k1*tan_e
    T[1, 0, 1] = h*tan_e2
    T[1, 0, 5] = -h*tan_e
    T[1, 2, 2] = (-k1 + h * h / 2. + h * h * tan_e2) * tan_e - h / 2. * h_pole * sec_e3
    T[1, 2, 3] = -h*tan_e2
    T[2, 0, 2] = h*tan_e2
    T[3, 0, 2] = -h*h_pole*sec_e3 - 2*k1*tan_e
    T[3, 0, 3] = -h*tan_e2
    T[3, 1, 2] = -h*sec_e2
    T[3, 2, 5] = h*tan_e - h*phi/np.cos(e - phi)**2
    # MADx
    if __MAD__:
        T[1, 2, 2] = (-k1 + h * h / 2. * (1 + h * h * sec_e2)) * tan_e - h / 2. * h_pole * sec_e3
        T[1, 0, 5] = 0
        T[3, 2, 5] = 0
    return R, T


def fringe_ext(h, k1, e, h_pole=0., gap=0., fint=0.):
    """
    Calculate first and second order matrices for edge focusing element.
    Note: the first order (R) is shown here in the sake of completeness but it is duplicated
            in r_matrix.py module and used from there in optics.py

    :param h: curvature (1/r) of the bend in [1/m]
    :param k1: strength of quadrupole component of the dipole in [1/m^2],
    :param e: the angle of inclination of the exit face [rad],
    :param h_pole: the curvature (1/r) of the exit face
    :param gap: the magnet gap [m], NOTE in MAD and ELEGANT: HGAP = gap/2
    :param fint: fringe field integral
    :return: R, T. First and second order matrices (6, 6), (6, 6, 6)
    """
    sec_e = 1./np.cos(e)
    sec_e2 = sec_e*sec_e
    sec_e3 = sec_e2*sec_e
    tan_e = np.tan(e)
    tan_e2 = tan_e*tan_e
    tan_e3 = tan_e2 * tan_e
    phi = fint*h*gap*sec_e*(1. + np.sin(e)**2)
    h2 = h*h
    R = np.eye(6)

    R[1, 0] = h*tan_e
    R[3, 2] = -h*np.tan(e - phi)

    T = np.zeros((6,6,6))
    T[0, 0, 0] = h/2.*tan_e2
    T[0, 2, 2] = -h/2.*sec_e2
    T[1, 0, 0] = h/2.*h_pole*sec_e3 - (-k1 + h*h/2.*tan_e2)*tan_e
    T[1, 0, 1] = -h*tan_e2
    T[1, 0, 5] = -h*tan_e
    T[1, 2, 2] = (-k1 - h*h/2.*tan_e2)*tan_e - h/2.*h_pole*sec_e3
    T[1, 2, 3] = h*tan_e2
    T[2, 0, 2] = -h*tan_e2
    T[3, 0, 2] = -h * h_pole * sec_e3 + (-k1 + h * h * sec_e2) * tan_e
    T[3, 0, 3] = h*tan_e2
    T[3, 1, 2] = h*sec_e2
    T[3, 2, 5] = h*tan_e - h*phi/np.cos(e - phi)**2
    # MAD
    if __MAD__:
        T[3, 0, 2] = -h * h_pole * sec_e3 + (-2 * k1 + h * h * sec_e2) * tan_e
        T[1, 0, 5] = 0
        T[3, 2, 5] = 0
    return R, T


def H23(vec_x, h, k1, k2, beta=1., g_inv=0.):
    """
    {x, px}, {y, py}, {sigma, ps}
    sigma = s - ct
    ps = (E - E0)/(p0*c)
    H2 = (px**2 + py**2)/2 + (h**2 + k1)*x**2/2 - (k1*y**2)/2 - (h*pt*x)/beta
    H3 = (h*x - ps/beta)*(px**2 + py**2)/2 + (2*h*k1 + k2)*(x**3)/6 - (h*k1 + k2)*(x*y**2)/2 - ps/(beta*gamma**2*(1. + beta))
    H23 = H2 + H3
    :param vec_x: [x, px, y, py, sigma, psigma]
    :param h: curvature
    :param k1: quadrupole strength
    :param k2: sextupole strength
    :param beta: = 1, velocity
    :param g_inv: 1/gamma, by default 0.
    :return: [x', px', y', py', sigma', psigma']
    """
    #print "H23: ", vec_x
    x = vec_x[0]
    px = vec_x[1]
    y = vec_x[2]
    py = vec_x[3]
    ps = vec_x[5]
    px2 = px*px
    py2 = py*py
    x2 = x*x
    y2 = y*y
    ps_beta = ps/beta
    k = 1. + h*x - ps_beta
    x1 = px*k
    px1 = -(h*h + k1)*x + h*ps_beta + (-h*(px2 + py2) - (2.*h*k1 + k2)*x2 + (h*k1 + k2)*y2)/2.
    y1 = py*k
    py1 = k1*y + (h*k1 + k2)*x*y
    sigma1 = -(h*x)/beta - ((px2 + py2)/(2.*beta)) - g_inv*g_inv/(beta*(1. + beta))
    return [x1, px1, y1, py1, sigma1, 0.]


def verlet(vec_x, step, h, k1, k2, beta=1., g_inv=0.):
    """
    q_{n+1} = q_{n} + h * dH(p_{n}, q_{n+1})/dp
    p_{n+1} = p_{n} - h * dH(p_{n}, q_{n+1})/dq
    """
    #vec_x0 = copy(vec_x)
    x =     vec_x[0]
    px =    vec_x[1]
    y =     vec_x[2]
    py =    vec_x[3]
    sigma = vec_x[4]
    ps =    vec_x[5]

    px2_py2 = px*px + py*py
    ps_beta = ps/beta
    x1 = (x + step*px*(1. - ps_beta))/(1. - step*h*px)
    y1 = y + step*py*(1. + h*x1 - ps_beta)
    vec_x[4] = sigma + step*(-h*x1/beta - px2_py2/(2.*beta) - g_inv*g_inv/(beta*(1. + beta)))
    px_d = -(h*h + k1 + (h*k1 + k2/2.)*x1)*x1 + h*ps_beta + (-h*px2_py2 + (h*k1 + k2)*y1*y1)/2.
    py_d = (k1 + (h*k1 + k2)*x1)*y1
    vec_x[1] = px + step*px_d
    vec_x[3] = py + step*py_d
    vec_x[0] = x1
    vec_x[2] = y1
    return vec_x

def verlet1O(vec_x, step, h, k1, k2, beta=1., g_inv=0.):
    """
    q_{n+1} = q_{n} + h * dH(p_{n}, q_{n+1})/dp
    p_{n+1} = p_{n} - h * dH(p_{n}, q_{n+1})/dq
    """
    x =     vec_x[0]
    px =    vec_x[1]
    y =     vec_x[2]
    py =    vec_x[3]
    sigma = vec_x[4]
    ps =    vec_x[5]
    #px2_py2 = px*px + py*py
    ps_beta = ps/beta

    x1 = x + step*px
    y1 = y + step*py
    vec_x[4] = sigma - step*h*x1/beta

    # derivatives
    #px_d = -(h*h + k1)*x1 + h*ps_beta + (-h*px2_py2 - (2.*h*k1 + k2)*x1*x1 + (h*k1 + k2)*y1*y1)/2.
    #py_d = k1*y1 + (h*k1 + k2)*x1*y1

    vec_x[1] = px - step*((h*h + k1)*x1 - h*ps_beta)
    vec_x[3] = py + step*k1*y1
    vec_x[0] = x1
    vec_x[2] = y1
    return vec_x

#import
#import pyximport; pyximport.install()
#import high_order


#from time import time
def sym_map(z, X, h, k1, k2, energy=0.):

    if h != 0. or k1 != 0. or k2 != 0.:
        step = 0.005
    else:
        step = z
    if step > z:
        step = z
    if step == 0.:
        return X
    n = int(z/step) + 1
    gamma = energy/m_e_GeV
    g_inv = 0.
    beta = 1.
    if gamma != 0.:
        g_inv = 1./gamma
        beta = np.sqrt(1. - g_inv*g_inv)
    z_array = np.linspace(0., z, num=n)
    step = z_array[1] - z_array[0]

    x =     X[0::6]
    px =    X[1::6]
    y =     X[2::6]
    py =    X[3::6]
    sigma = X[4::6]
    ps =    X[5::6]

    ps_beta = ps/beta
    #vec = [x, px, y, py, sigma, ps]
    c1 = h*h + k1
    c2 = h*k1 + k2/2.
    c3 = h*k1 + k2
    c4 = g_inv*g_inv/(beta*(1. + beta))
    #Y = high_order.cython_test([x, px, y, py, sigma, ps], step, h,k1, len(z_array) , (c1,c2,c3,c4), ps_beta, beta)
    #x, px, y, py, sigma, ps = Y
    for z in range(len(z_array) - 1):
        #vec = verlet(vec, step, h, k1, k2, beta=beta, g_inv=g_inv)
        px2_py2 = px*px + py*py
        x = (x + step*px*(1. - ps_beta))/(1. - step*h*px)
        y = y + step*py*(1. + h*x - ps_beta)
        sigma = sigma + step*(-h*x/beta - px2_py2/(2.*beta) - c4)

        px = px + step*(h*ps_beta + (-h*px2_py2 + c3*y*y)/2. - (c1 + c2*x)*x)
        py = py + step*(k1 + c3*x)*y

    X[0::6] = x[:] #vec[0][:]
    X[1::6] = px[:] #vec[1][:]
    X[2::6] = y[:] #vec[2][:]
    X[3::6] = py[:] #vec[3][:]
    X[4::6] = sigma[:] #vec[4][:]
    return X



"""


        I111 = 1./3.*(sx2 + dx_h)                                                    #  I111 = Gx * cx**2
        I122 = dx_h*dx_h/3.                                                          #  I122 = Gx * sx**2
        I112 = sx*dx_h/3.                                                            #  I112 = Gx * cx*sx
        I11  = L*sx/2.                                                               #  I11  = Gx * cx
        I116 = h/kx2*(I11 - I111)                  if kx != 0 else h*L4/24.          #  I116 = Gx * cx*dx
        I12  = 0.5/kx2*(sx - L*cx)                 if kx != 0 else L3/6.             #  I12  = Gx * sx
        I126 = h/kx2*(I12 - I112)                  if kx != 0 else h*L5/40.          #  I126 = Gx * sx*dx
        I10  = dx_h                                                                  #  I10  = Gx
        I16  = h/kx2*(dx_h - L*sx/2.)              if kx != 0 else h*L4/24.          #  I16  = Gx * dx
        I166 = h2/kx4*(I10 - 2*I11 + I111)         if kx != 0 else h2*L5*L/120.      #  I166 = Gx * dx**2
        I144 = (sy2 - 2.*dx_h)/denom               if non_drift else L4/12.          #  I144 = Gx * sy**2
        I133 = dx_h - ky2*(sy2 - 2.*dx_h)/denom    if non_drift else L2/2.           #  I133 = Gx * cy**2
        I134 = (sy*cy - sx)/denom                  if non_drift else L3/6.           #  I134 = Gx * cy*sy
        I313 = (kx2*cy*dx_h - 2.*ky2*sx*sy)/denom  if non_drift else L2/2.           #  I313 = Gy * cx*cy
        I324 = (2.*cy*dx_h - sx*sy)/denom          if non_drift else L4/12.          #  I324 = Gy * sx*sy
        I314 = (2.*cy*sx - (1. + cx)*sy)/denom     if non_drift else L3/6.           #  I314 = Gy * cx*sy
        I323 = (sy - cy*sx - 2.*ky2*sy*dx_h)/denom if non_drift else L3/6.           #  I323 = Gy * sx*cy = (2*ky2/kx2*(1 + cx)*sy - cy*sx)/denom + sy/kx2
        I33  = L*sy/2.                                                               #  I33  = Gy * cy
        I34  = (sy - L*cy)/(2.*ky2)                if ky !=0. else L3/6.             #  I34  = Gy * sy
        I336 = h/kx2*(I33 - I313)                                                    #  I336 = Gy * dx*cy
        I346 = h/kx2*(I34 - I314)                                                    #  I346 = Gy * dx*sy

        #derivative of Integrals
        I211 = sx/3.*(1. + 2.*cx)
        I222 = 2.*dx_h*sx/3.
        I212 = 1./3.*(2*sx2 - dx_h)
        I21  = 1./2.*(L*cx + sx)
        I216 = h/kx2*(I21 - I211)
        I22  = I11
        I226 = h/kx2*(I22 - I212)
        I20  = sx
        I26  = h /(2.*kx2)*(sx - L*cx)
        I266 = h2/kx4*(I20 - 2.*I21 + I211)
        I244 = 2.*(cy*sy - sx)/denom                           if non_drift else L3/3.
        I233 = sx + 2.*ky2*(cy*sy - sx)/denom                  if non_drift else L
        I234 = (kx2*dx_h - 2.*ky2*sy2)/denom                   if non_drift else L2/2.
        I413 = ((kx2 - 2.*ky2)*cy*sx - ky2*sy*(1. + cx))/denom if non_drift else L
        I424 = (cy*sx - cx*sy - 2.*ky2*sy*dx_h)/denom          if non_drift else L3/3.
        I414 = ((kx2 - 2.*ky2)*sx*sy - (1. - cx)*cy)/denom     if non_drift else L2/2.
        I423 = (cy*dx_h*(kx2 - 2*ky2) - ky2*sx*sy)/denom       if non_drift else L2/2.   #  I423 = I323' = ((2.*ky2)/kx2*(1 + cx)*cy - cx*cy - ky2*sx*sy)/denom + cy/kx2
        I43  = 0.5*(L*cy + sy)
        I436 = h/kx2*(I43 - I413)
        I44  = I33
        I446 = h/kx2*(I44 - I414)
"""


def moments_py(bx, by, Bx, By, Bz, dzk):
    """
    mx = v/(dz/dt)*e/p*(y'*Bz - By*(1.+x'^2) + x'*y'*Bx)*dz
    my = -v/(dz/dt)*e/p*(x'*Bz - Bx*(1.+y'^2) + x'*y'*By)*dz
    :param bx: x' = dx/dz
    :param by: y' = dy/dz
    :param Bx: horizontal component of mag field
    :param By: vertical component of mag field
    :param Bz: longitudinal component of mag field
    :param dzk: dz*e/p; dz - step, e - charge, p - momentum
    :return: mx, my
    """
    bx2 = bx*bx
    by2 = by*by
    bxy = bx*by
    sq = np.sqrt(1. + bx2 + by2)
    k = sq*dzk
    mx = k*(by*Bz - By*(1.+bx2) + bxy*Bx)
    my = -k*(bx*Bz - Bx*(1.+by2) + bxy*By)
    return mx, my

moments = moments_py if not nb_flag else nb.jit(moments_py)


def rk_track_in_field(y0, s_stop, N, energy, mag_field, s_start=0.):
    """
    Runge-Kutta solver of the exact trajectory equations in the fixed coordinate system {X, Y, Z]}.
    X, Y - transverse coordinates, Z - longitudinal.

    Trajectory equations:
    x'' = v/(dz/dt)*e/p*(y'*Bz - By*(1.+x'^2) + x'*y'*Bx)
    y'' = -v/(dz/dt)*e/p*(x'*Bz - Bx*(1.+y'^2) + x'*y'*By)
    v/(dz/dt) = sqrt(1 + y'^2 + x'^2)
    z' is not needed in the initial conditions and z' is not used for SR calculation
        therefore z' can be arbitrary defined in the returning array. z' = dE/(ps*c)

    :param y0: array n*6; initial coordinates of n particles
                [x0, x0'=dx/dz, y0, y0=dy/dz, 0, dE0/pc,
                x1, x1', y1, y1', 0, dE1/pc, ...
                xn, xn', yn, yn', 0, dEn/pc]
    :param s_stop: longitudinal coordinate of stop
    :param N: number of steps
    :param energy: energy of particle in [GeV]
    :param mag_field: function. Bx, By, Bz = mag_field(X, Y, Z)
    :param s_start: starting longitudinal coordinate
    :return: array with length N*9- coordinates and magnetic fields on the trajectory
            [x, x', y, y', z, dE/pc, Bx, By, Bz, ...
            xn, xn', yn, yn', zn, zn', Bxn, Byn, Bzn]
    """

    # pc_ref = np.sqrt(energy**2 - m_e_GeV**2)
    z = np.linspace(s_start, s_stop, num=N)
    h = z[1] - z[0]
    N = len(z)
    gamma0 = energy/m_e_GeV
    charge = 1
    mass = 1 #in electron mass
    u = np.zeros((N*9, np.shape(y0)[1]))
    dz = h
    beta0 = np.sqrt(1. - 1./(gamma0*gamma0))
    gammai = gamma0*(1 + y0[5]*beta0)
    k = charge*speed_of_light/(m_e_eV*mass*gammai)
    u[0, :] = y0[0]
    u[1, :] = y0[1]
    u[2, :] = y0[2]
    u[3, :] = y0[3]
    u[4, :] = z[0] #y0[4]
    u[5, :] = y0[5]
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

        u[(i + 1) * 9 + 0] = X + 1/6.*(kx1 + 2.*(kx2 + kx3) + kx4)
        u[(i + 1) * 9 + 1] = bxconst + 1/6.*(mx1 + 2.*(mx2 + mx3) + mx4) #// conversion in mrad
        u[(i + 1) * 9 + 2] = Y + 1/6.*(ky1 + 2.*(ky2 + ky3) + ky4)
        u[(i + 1) * 9 + 3] = byconst + 1/6.*(my1 + 2.*(my2 + my3) + my4)
        u[(i + 1) * 9 + 4] = Z_n #u[i*9 + 4] + dz*np.sqrt(1 + u[(i+1)*9 + 1]**2 + u[(i+1)*9 + 3]**2)
        u[(i + 1) * 9 + 5] = y0[5]
        # beta_z as 6-th coordinate
        #u[(i+1)*9 + 5] = betai/np.sqrt(1 + u[(i+1)*9 + 1]*u[(i+1)*9 + 1] + u[(i+1)*9 + 3]*u[(i+1)*9 + 3])#dGamma2 - (u[(i+1)*9 + 1]*u[(i+1)*9 + 1] + u[(i+1)*9 + 3]*u[(i+1)*9 + 3])/2.

    u[(N-1)*9 + 6], u[(N-1)*9 + 7], u[(N-1)*9 + 8] = mag_field(u[(N-1)*9 + 0], u[(N-1)*9 + 2], u[(N-1)*9 + 4])
    return u


def rk_field(rparticles, s_start, s_stop, N, energy, mag_field, long_dynamics=True):
    """
    Method to track particles through mag field

    :param rparticles: initial coordinates of the particles shape (6xN), N number of particles (e.g. ParticleArray.rparticles)
    :param s_start:
    :param s_stop:
    :param N: number of points on the trajectory
    :param energy: in GeV
    :param mag_field: must be function e.g. lambda x, y, z: (Bx, By, Bz)
    :param long_dynamics: True, if True, includes longitudinal dynamics, otherwise only transverse
    :return:
    """
    if s_start > s_stop:
        print("rk_field: s_start > s_stop. Setup s_start = 0")
        s_start = 0.

    ref_path = 0

    if long_dynamics:
        traj_ref = rk_track_in_field(np.array([[0], [0], [0], [0], [0], [0]]), s_stop, N, energy, mag_field, s_start=s_start)
        x1 = traj_ref[1+9::9]
        y1 = traj_ref[3+9::9]
        dz = traj_ref[4+9::9] - traj_ref[4:-9:9]
        ref_path = np.sum(dz*np.sqrt(1 + x1*x1 + y1*y1))

    traj_data = rk_track_in_field(rparticles, s_stop, N, energy, mag_field, s_start=s_start)

    z_fin = rparticles[4, :]

    if long_dynamics:
        x1 = traj_data[1+9::9, :]
        y1 = traj_data[3+9::9, :]
        dz = traj_data[4+9::9, :] - traj_data[4:-9:9, :]
        z_fin += traj_data[4, :] + np.sum(dz*np.sqrt(1 + x1*x1 + y1*y1), axis=0) - ref_path

    rparticles[0, :] = traj_data[(N-1)*9 + 0, :]
    rparticles[1, :] = traj_data[(N-1)*9 + 1, :]
    rparticles[2, :] = traj_data[(N-1)*9 + 2, :]
    rparticles[3, :] = traj_data[(N-1)*9 + 3, :]
    rparticles[4, :] = z_fin
    rparticles[5, :] = traj_data[(N-1)*9 + 5, :]
    return rparticles


def scipy_track_in_field(y0, l, N, energy, mag_field):# y0, l, N, energy, mag_field
    z = np.linspace(0.,l, num=N)
    def func(y, z, fields):
        gamma = energy/m_e_GeV
        charge = 1
        mass = 1 #in electron mass
        cmm = speed_of_light
        massElectron = m_e_eV #0.510998910e+6 #// rest mass of electron
        k = charge*cmm/(massElectron*mass*gamma)
        n = len(y)/4
        X = y[0:n]
        bx = y[n:2*n]
        Y = y[2*n:3*n]
        by = y[3*n:4*n]
        bx2 = bx*bx
        by2 = by*by
        bxy = bx*by
        #sq = 1 + (bx2 + by2)/2.
        sq = np.sqrt(1. + bx2 + by2)
        k = sq*k
        Bx, By, Bz = fields(X, Y, z)
        mx = k*(by*Bz - By*(1.+bx2) + bxy*Bx)
        my = -k*(bx*Bz - Bx*(1.+by2) + bxy*By)
        y = np.append(bx, (mx, by, my))
        return np.array(y) #[y[1],mx, y[3], my]

    n = len(y0)/6
    x = y0[0::6]
    bx = y0[1::6]
    y = y0[2::6]
    by = y0[3::6]
    Y0 = np.append(x, (bx, y, by))

    u = odeint(func, Y0, z, args = (mag_field,), rtol=1e-10, atol=1e-10)

    w = np.zeros((6*len(z), n))
    w[0::6] = u[:, 0:n]
    w[1::6] = u[:, n:2*n]
    w[2::6] = u[:, 2*n:3*n]
    w[3::6] = u[:, 3*n:4*n]
    N = len(z)
    sub = y0[4::6]
    sub = np.tile(sub, N).reshape((N, n))
    w[4::6] = sub[:]
    sub = y0[5::6]
    sub = np.tile(sub, N).reshape((N, n))
    w[5::6] = sub[:]
    return w


def rk_field_scipy(y0, l, N, energy, mag_field):
    #z = linspace(0, l, num=N)
    #h = z[1]-z[0]
    #N = len(z)
    traj_data = scipy_track_in_field(y0, l, N, energy, mag_field)
    #print np.shape(traj_data), np.shape(y0)
    #print traj_data
    y0[0::6] = traj_data[(N-1)*6 + 0,:]
    y0[1::6] = traj_data[(N-1)*6 + 1,:]
    y0[2::6] = traj_data[(N-1)*6 + 2,:]
    y0[3::6] = traj_data[(N-1)*6 + 3,:]
    #y0[4::6] = traj_data[(N-1)*9 + 4,:]
    #y0[5::6] = traj_data[(N-1)*9 + 5,:]
    return y0


def track_und_weave(y0, z, kz, kx ,Kx, energy):
    from scipy import weave
    gamma = energy/m_e_GeV
    #rho = gamma/Kx/kz
    c = speed_of_light
    m0 = m_e_eV
    B0 = Kx*m0*kz/c
    #y = array(y0)
    ky = np.sqrt(kz*kz - kx*kx)
    #print "B0 = ", B0
    #print "rho = ", rho
    #print "kz = ", kz
    #print "kx = ", kx
    #print "gamma = ", gamma
    h = z[1]-z[0]
    N = len(z)
    #Ax =  1/kz * np.cos(kx*y0[0])*np.cosh(ky*y0[2])*np.sin(kz*z[0])
    #Ay =  kx/(ky*kz) * np.sin(kx*y0[0])*np.sinh(ky*y0[2])*np.sin(kz*z[0])
    q = np.array([y0[0],y0[1] ,y0[2],y0[3], y0[4], y0[5], kx, ky, kz, h, N, B0, gamma*(1+y0[5])])
    #print N

    u = np.zeros(N*6)
    support_code = """
    extern "C" {
    void fields(double x, double y, double z, double kx, double ky, double kz, double B0, double *Bx, double *By, double *Bz)
    {
        double k1 =  -B0*kx/ky;
        double k2 = -B0*kz/ky;
        double kx_x = kx*x;
        double ky_y = ky*y;
        double kz_z = kz*z;
        double cosx = cos(kx_x);
        double sinhy = sinh(ky_y);
        double cosz = cos(kz_z);
        *Bx = k1*sin(kx_x)*sinhy*cosz ;
        *By = B0*cosx*cosh(ky_y)*cosz;
        *Bz = k2*cosx*sinhy*sin(kz_z);
    };
    void moments(double bx, double by, double Bx, double By, double Bz, double dzk, double *mx, double *my)
    {
        double bx2 = bx*bx;
        double by2 = by*by;
        double bxy = bx*by;
        //#sq = 1 + (bx2 + by2)/2.
        double k = sqrt(1. + bx2 + by2)*dzk;
        *mx = k*(by*Bz - By*(1.+bx2) + bxy*Bx);
        *my = -k*(bx*Bz - Bx*(1.+by2) + bxy*By);
    }
    }
    """
    code = """

    double charge = 1;
    double mass = 1; //in electron mass

    double X, Y, Z;
    double bx, by;
    double kx1, ky1;
    double kx2, ky2;
    double kx3, ky3;
    double kx4, ky4;
    double mx1, my1;
    double mx2, my2;
    double mx3, my3;
    double mx4, my4;
    double bxconst;
    double byconst;

    double cmm = 299792458;
    //double k;
    double massElectron = 0.510998910e+6; // rest mass of electron
    double dzk ;
    double Bx;
    double By;
    double Bz;


    double px = Q1(1);
    double py = Q1(3);
    double kx = Q1(6);
    double ky = Q1(7);
    double kz = Q1(8);

    double dz = Q1(9);
    int N = Q1(10);
    double B0 = Q1(11);
    double gamma = Q1(12);
    double dGamma2 = 1. - 0.5/(gamma*gamma);
    double pz = dGamma2 - (px*px + py*py)/2.;

    int i;

    double k = charge*cmm/(massElectron*mass*gamma);

    U1(0) = Q1(0);
    U1(1) = Q1(1);
    U1(2) = Q1(2);
    U1(3) = Q1(3);
    U1(4) = 0.;
    U1(5) = pz;
    dzk = dz*k;
    for(i = 0; i < N-1; i++)
    {


        X = U1(i*6 + 0);
        Y = U1(i*6 + 2);
        Z = U1(i*6 + 4);

        bxconst = U1(i*6 + 1);
        byconst = U1(i*6 + 3);

        bx = bxconst;
        by = byconst;
        kx1 = bx*dz;
        ky1 = by*dz;

        fields( X,  Y,  Z,  kx,  ky,  kz,  B0, &Bx, &By, &Bz);
        moments( bx,  by,  Bx,  By,  Bz,  dzk,  &mx1,  &my1);
        //K2
        bx = bxconst + mx1/2.;
        by = byconst + my1/2.;
        kx2 = bx*dz;
        ky2 = by*dz;
        fields( X + kx1/2.,  Y + ky1/2.,  Z + dz/2.,  kx,  ky,  kz,  B0, &Bx, &By, &Bz);
        moments( bx,  by,  Bx,  By,  Bz,  dzk,  &mx2,  &my2);
        // K3
        bx = bxconst + mx2/2.;
        by = byconst + my2/2.;
        kx3 = bx*dz;
        ky3 = by*dz;
        fields( X + kx2/2.,  Y + ky2/2.,  Z + dz/2.,  kx,  ky,  kz,  B0, &Bx, &By, &Bz);
        moments( bx,  by,  Bx,  By,  Bz,  dzk,  &mx3,  &my3);
        //K4
        bx = bxconst + mx3;
        by = byconst + my3;
        kx4 = bx*dz;
        ky4 = by*dz;
        fields( X + kx3,  Y + ky3,  Z + dz,  kx,  ky,  kz,  B0, &Bx, &By, &Bz);
        moments( bx,  by,  Bx,  By,  Bz,  dzk,  &mx4,  &my4);

        U1((i+1)*6 + 0) = X + 1/6.*(kx1 + 2.*kx2 + 2.*kx3 + kx4);
        U1((i+1)*6 + 1) = bxconst + 1/6.*(mx1 + 2.*mx2 + 2.*mx3 + mx4); // conversion in mrad
        U1((i+1)*6 + 2) = Y + 1/6.*(ky1 + 2.*ky2 + 2.*ky3 + ky4);
        U1((i+1)*6 + 3) = byconst + 1/6.*(my1 + 2.*my2 + 2.*my3 + my4);

        U1((i+1)*6 + 4) = Z + dz;
        U1((i+1)*6 + 5) = dGamma2 - (U1((i+1)*6 + 1)*U1((i+1)*6 + 1) + U1((i+1)*6 + 3)*U1((i+1)*6 + 3))/2.; //bz;


    }
    //err = B3D(field, motion->X[i+1], motion->Y[i+1], motion->Z[i+1] + Zshift, pBx, pBy, pBz);
    //std::cout<<" B = "<<megaPack->motion.Z[i+1] <<std::endl;
    //motion->Bx[i+1] = Bx;
    //motion->By[i+1] = By;
    //motion->Bz[i+1] = Bz;


    """
    weave.inline(code, ['u',"q"], support_code= support_code)
    x = u[::6]
    y = u[2::6]
    #print x
    #print y
    px = u[1::6]
    py = u[3::6]
    z = u[4::6]
    pz = u[5::6]
    return x, px,y, py, z, pz


def track_und_weave_openmp(u, l, N, kz, kx ,Kx, energy):
    from scipy import weave
    import numpy
    gamma = energy/m_e_GeV
    #rho = gamma/Kx/kz
    c = speed_of_light
    m0 = m_e_eV
    B0 = Kx*m0*kz/c

    ky = np.sqrt(kz*kz - kx*kx)
    z = np.linspace(0, l, num=N)
    h = z[1]-z[0]
    N = len(z)

    nparticles = len(u)/6
    q = np.array([ kx, ky, kz, h, N, B0, gamma, nparticles])
    code = """
    double charge = 1;
    double mass = 1; //in electron mass

    double cmm = 299792458;

    double massElectron = 0.510998910e+6; // rest mass of electron

    double kx = Q1(0);
    double ky = Q1(1);
    double kz = Q1(2);

    double dz = Q1(3);
    int N = Q1(4);
    double B0 = Q1(5);

    int npart = Q1(7);

    double sq;
    #pragma omp parallel for
    for(int n = 0; n < npart; n++)
    {
        double kx1, ky1;
        double kx2, ky2;
        double kx3, ky3;
        double kx4, ky4;
        double mx1, my1;
        double mx2, my2;
        double mx3, my3;
        double mx4, my4;
        double bx2;
        double by2;
        double gamma = Q1(6);//*(1+U1(n*6 + 5));
        double dGamma2 = 1. - 0.5/(gamma*gamma);
        double k = charge*cmm/(massElectron*mass*gamma);

        double X = U1(n*6 + 0);
        double bxconst = U1(n*6 + 1);
        double Y = U1(n*6 + 2);
        double byconst = U1(n*6 + 3);
        double Bx;
        double By;
        double Bz;
        //double pz = dGamma2 - (bxconst*bxconst + byconst*byconst)/2.;
        //bz = pz;
        double dzk = dz*k;
        double Z = 0.;

        for(int i = 0; i < N-1; i++)
        {
            double x = X;
            double y = Y;
            double z = Z;
            double bx = bxconst;
            double by = byconst;

            Bx = -B0*kx/ky*sin(kx*x)*sinh(ky*y)*cos(kz*z); // here kx is only real
            By = B0*cos(kx*x)*cosh(ky*y)*cos(kz*z);
            Bz = -B0*kz/ky*cos(kx*x)*sinh(ky*y)*sin(kz*z);

            //motion->XbetaI2[i] = motion->By[i];
            kx1 = bx*dz;
            ky1 = by*dz;
            bx2 = bx*bx;
            by2 = by*by;
            //sq = 1 + (bx2 + by2)/2.;
            sq = sqrt(1 + bx2 + by2);

            mx1 = sq*(by*Bz - By*(1+bx2) + bx*by*Bx)*dzk;
            my1 = -sq*(bx*Bz - Bx*(1+by2) + bx*by*By)*dzk;
            //K2
            //err = B3D(field, X + kx1/2., Y + ky1/2., Z + dz/2. + Zshift, pBx, pBy, pBz);
            x = X + kx1/2.;
            y = Y + ky1/2.;
            z = Z + dz/2.;
            Bx = -B0*kx/ky*sin(kx*x)*sinh(ky*y)*cos(kz*z); // here kx is only real
            By = B0*cos(kx*x)*cosh(ky*y)*cos(kz*z);
            Bz = -B0*kz/ky*cos(kx*x)*sinh(ky*y)*sin(kz*z);

            bx = bxconst + mx1/2.;
            by = byconst + my1/2.;

            kx2 = bx*dz;
            ky2 = by*dz;
            bx2 = bx*bx;
            by2 = by*by;
            //sq = 1 + (bx2 + by2)/2.;
            sq = sqrt(1 + bx*bx + by*by);
            mx2 = sq*(by*Bz - By*(1+bx2) + bx*by*Bx)*dzk;
            my2 = -sq*(bx*Bz - Bx*(1+by2) + bx*by*By)*dzk;
            // K3
            //err = B3D(field, X + kx2/2., Y + ky2/2., Z + dz/2. + Zshift, pBx, pBy, pBz);

            x = X + kx2/2.;
            y = Y + ky2/2.;
            z = Z + dz/2.;
            Bx = -B0*kx/ky*sin(kx*x)*sinh(ky*y)*cos(kz*z); // here kx is only real
            By = B0*cos(kx*x)*cosh(ky*y)*cos(kz*z);
            Bz = -B0*kz/ky*cos(kx*x)*sinh(ky*y)*sin(kz*z);

            bx = bxconst + mx2/2.;
            by = byconst + my2/2.;

            kx3 = bx*dz;
            ky3 = by*dz;
            bx2 = bx*bx;
            by2 = by*by;
            //sq = 1 + (bx2 + by2)/2.;
            sq = sqrt(1 + bx*bx + by*by);
            mx3 = sq*(by*Bz - By*(1+bx2) + bx*by*Bx)*dzk;
            my3 = -sq*(bx*Bz - Bx*(1+by2) + bx*by*By)*dzk;

            //K4
            //err = B3D(field, X + kx3, Y + ky3, Z + dz + Zshift,pBx, pBy, pBz);

            x = X + kx3;
            y = Y + ky3;
            z = Z + dz;
            Bx = -B0*kx/ky*sin(kx*x)*sinh(ky*y)*cos(kz*z); // here kx is only real
            By = B0*cos(kx*x)*cosh(ky*y)*cos(kz*z);
            Bz = -B0*kz/ky*cos(kx*x)*sinh(ky*y)*sin(kz*z);

            bx = bxconst + mx3;
            by = byconst + my3;


            kx4 = bx*dz;
            ky4 = by*dz;
            bx2 = bx*bx;
            by2 = by*by;
            //sq = 1 + (bx2 + by2)/2.;
            sq = sqrt(1 + bx*bx + by*by);
            mx4 = sq*(by*Bz - By*(1+bx2) + bx*by*Bx)*dzk;
            my4 = -sq*(bx*Bz - Bx*(1+by2) + bx*by*By)*dzk;

            X = X + 1/6.*(kx1 + 2.*kx2 + 2.*kx3 + kx4);
            bxconst = bxconst + 1/6.*(mx1 + 2.*mx2 + 2.*mx3 + mx4); // conversion in mrad
            Y= Y + 1/6.*(ky1 + 2.*ky2 + 2.*ky3 + ky4);
            byconst = byconst + 1/6.*(my1 + 2.*my2 + 2.*my3 + my4);

            Z = Z + dz;
            //bz = dGamma2 - (bxconst*bxconst +byconst*byconst)/2.; //bz;
        }
        U1(n*6 + 0) = X;
        U1(n*6 + 1) = bxconst;
        U1(n*6 + 2) = Y;
        U1(n*6 + 3) = byconst;
    }
    """
    weave.inline(code, ['u',"q"],
    extra_compile_args =['-O3 -fopenmp'],
    compiler = 'gcc',
    libraries=['gomp'],
    headers=['<omp.h>']
    )
    return u

def arcline( SREin, Delta_S, dS, R_vect ):
    """
    Martin Dohlus DESY, 2015

    :param SREin: trajectory. traj[0,:] - longitudinal coordinate (s),
                             traj[1,:],traj[2,:],traj[3,:] - rectangular coordinates, (x, y, z)\
                             traj[4,:],traj[5,:],traj[6,:] - tangential unit vectors, (x', y', z')
    :param Delta_S:
    :param dS:
    :param R_vect: radius
    :return:
    """

    epsilon = 1e-8

    sre0 = SREin[:, -1]
    N = int(max(1, np.round(Delta_S/dS)))
    dS = float(Delta_S)/N
    SRE2 = np.zeros((7, N))
    SRE2[0,:] = sre0[0] + np.arange(1, N+1)*dS

    R_vect_valid = False

    if True:
        R_vect_valid = np.all(np.isfinite(R_vect)) and not np.all(np.equal(R_vect, np.zeros((3, 1))))
        if R_vect_valid:
            R = norm(R_vect)
            n_vect = R_vect/R
            e1 = sre0[4:7]
            if np.abs(np.dot(n_vect, e1)) > epsilon:
                R_vect_valid = False
                print('*** error in arcline: invalid R_vect --> using Drift. Consider to use Runge-Kutta integrator')

    if not R_vect_valid:
        SRE2[2-1, :] = sre0[2-1] + sre0[5-1]*np.arange(1, N+1)*dS
        SRE2[3-1, :] = sre0[3-1] + sre0[6-1]*np.arange(1, N+1)*dS
        SRE2[4-1, :] = sre0[4-1] + sre0[7-1]*np.arange(1, N+1)*dS
        SRE2[5-1, :] = sre0[5-1]
        SRE2[6-1, :] = sre0[6-1]
        SRE2[7-1, :] = sre0[7-1]
    else:
        R = norm(R_vect)
        n_vect = R_vect/R
        e2 = np.array([[n_vect[1]*e1[2] - n_vect[2]*e1[1]],
            [n_vect[2]*e1[0] - n_vect[0]*e1[2]],
            [n_vect[0]*e1[1] - n_vect[1]*e1[0]]])
        si = np.sin(np.arange(1, N+1)*dS/R)
        co = np.cos(np.arange(1, N+1)*dS/R)
        omco = 2*np.sin(np.arange(1, N+1)*dS/R/2)**2
        SRE2[2-1, :] = sre0[1] + R*(e1[0]*si + e2[0]*omco)
        SRE2[3-1, :] = sre0[2] + R*(e1[1]*si + e2[1]*omco)
        SRE2[4-1, :] = sre0[3] + R*(e1[2]*si + e2[2]*omco)
        SRE2[5-1, :] = e1[0]*co + e2[0]*si
        SRE2[6-1, :] = e1[1]*co + e2[1]*si
        SRE2[7-1, :] = e1[2]*co + e2[2]*si

    SRE = np.append(SREin, SRE2, axis=1)

    return SRE
