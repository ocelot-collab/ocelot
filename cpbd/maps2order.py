__author__ = 'Sergey Tomin'

from numpy import cos, sin, sqrt, zeros, eye, tan

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


def t_nnn(L, h, k1, k2):
    """
    :param L:
    :param angle:
    :param k1:
    :param k2:
    :return:

    here is used the following set of variables:
    x, dx/ds, y, dy/ds,
    """

    h2 = h*h
    h3 = h2*h
    kx2 = (k1 + h*h)
    ky2 = -k1
    kx4 = kx2*kx2
    ky4 = ky2*ky2
    kx = sqrt(kx2 + 0.j)
    ky = sqrt(ky2 + 0.j)
    cx = cos(kx*L).real
    sx = (sin(kx*L)/kx).real if kx != 0 else L
    cy = cos(ky*L).real

    sy = (sin(ky*L)/ky).real if ky != 0 else L

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
        I233 = sx + 2.*ky2*(cy*sy - sx)/denom
        I234 = (kx2*dx_h - 2.*ky2*sy2)/denom
        I413 = ((kx2 - 2.*ky2)*cy*sx - ky2*sy*(1. + cx))/denom
        I424 = (cy*sx - cx*sy - 2.*ky2*sy*dx_h)/denom
        I414 = ((kx2 - 2.*ky2)*sx*sy - (1. - cx)*cy)/denom
        I423 = (cy*dx_h*(kx2 - 2*ky2) - ky2*sx*sy)/denom      #  I423 = I323' = ((2.*ky2)/kx2*(1 + cx)*cy - cx*cy - ky2*sx*sy)/denom + cy/kx2
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

    if kx == 0 and ky != 0:
        I336 = (h*L*(3.*L*cy + (2.*ky2*L2 - 3.)*sy))/(24.*ky2)
        I346 = (h*((3. - 2.*ky2*L2)*L*cy + 3.*(ky2*L2 - 1.)*sy))/(24.*ky4)
        I436 = I346
        I446 = (h*L*(-3.*L*cy + (3. + 2.*ky2*L2)*sy))/(24.*ky2)
    elif kx == 0 and ky == 0:
        I336 = (h*L4)/24.
        I346 = (h*L5)/40.
        I436 = (h*L3)/6.
        I446 = (h*L4)/8.
    else:
        I336 = h/kx2*(I33 - I313)                                  #  I336 = Gy * dx*cy
        I346 = h/kx2*(I34 - I314)                                  #  I346 = Gy * dx*sy
        I436 = h/kx2*(I43 - I413)
        I446 = h/kx2*(I44 - I414)

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

    t413 = coef2*I413 + h*kx2*ky2*I424
    t414 = coef2*I414 - h*kx2*I423
    t423 = coef2*I423 - h*ky2*I414
    t424 = coef2*I424 + h*I413
    t436 = coef2*I436 - h2*ky2*I424 + ky2*I43
    t446 = coef2*I446 + h2*I423 + ky2*I44

    # Coordinates transformation from Curvilinear to a Restangular
    cx_1 = -kx2*sx
    sx_1 = cx
    cy_1 = -ky2*sy
    sy_1 = cy
    dx_1 = h*sx
    T = zeros((6,6,6))
    T[0, 0, 0] = t111
    T[0, 0, 1] = t112 + h*sx
    T[0, 0, 5] = t116
    T[0, 1, 1] = t122
    T[0, 1, 5] = t126
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
    T[2, 3, 5] = t346

    T[3, 0, 2] = t413 - h*cx*cy_1
    T[3, 0, 3] = t414 + (1 - cx)*h*sy_1
    T[3, 1, 2] = t423 - h*sx*cy_1
    T[3, 1, 3] = t424 - h*sx*sy_1
    T[3, 2, 5] = t436 - h*dx*cy_1
    T[3, 3, 5] = t446 - h*dx*sy_1
    """
    print "t111 = ", t111
    print "t112 = ", t112
    print "t116 = ", t116
    print "t122 = ", t122
    print "t126 = ", t126
    print "t166 = ", t166
    print "t133 = ", t133
    print "t134 = ", t134
    print "t144 = ", t144
    print "t211 = ", t211
    print "t212 = ", t212
    print "t216 = ", t216
    print "t222 = ", t222
    print "t226 = ", t226
    print "t266 = ", t266
    print "t233 = ", t233
    print "t234 = ", t234
    print "t244 = ", t244
    print "t313 = ", t313
    print "t314 = ", t314
    print "t323 = ", t323
    print "t324 = ", t324
    print "t336 = ", t336
    print "t346 = ", t346
    print "t413 = ", t413
    print "t414 = ", t414
    print "t423 = ", t423
    print "t424 = ", t424
    print "t436 = ", t436
    print "t446 = ", t446
    """
    return T

def fringe_ent(h, k1,  e, h_pole = 0., gap = 0., fint = 0.):

    sec_e = 1./cos(e)
    sec_e2 = sec_e*sec_e
    sec_e3 = sec_e2*sec_e
    tan_e = tan(e)
    tan_e2 = tan_e*tan_e
    phi = fint*h*gap*sec_e*(1. + sin(e)**2)
    R = eye(6)
    R[1,0] = h*tan_e
    R[3,2] = -h*tan(e - phi)
    #print R

    T = zeros((6,6,6))
    T[0,0,0] = -h/2.*tan_e2
    T[0,2,2] = h/2.*sec_e2
    T[1,0,0] = h/2.*h_pole*sec_e3 + k1*tan_e
    T[1,0,1] = h*tan_e2
    T[1,0,5] = -h*tan_e
    T[1,2,2] = (-k1 + h*h/2. + h*h*tan_e2)*tan_e - h/2.*h_pole*sec_e3
    T[1,2,3] = -h*tan_e2
    T[2,0,2] = h*tan_e2
    T[3,0,2] = -h*h_pole*sec_e3 - 2*k1*tan_e
    T[3,0,3] = -h*tan_e2
    T[3,1,2] = -h*sec_e2
    T[3,2,5] = h*tan_e - h*phi/cos(e - phi)**2
    return R, T

def fringe_ext(h, k1,  e, h_pole = 0., gap = 0., fint = 0.):

    sec_e = 1./cos(e)
    sec_e2 = sec_e*sec_e
    sec_e3 = sec_e2*sec_e
    tan_e = tan(e)
    tan_e2 = tan_e*tan_e
    phi = fint*h*gap*sec_e*(1. + sin(e)**2)
    R = eye(6)

    R[1,0] = h*tan_e
    R[3,2] = -h*tan(e - phi)
    #print R

    T = zeros((6,6,6))
    T[0,0,0] = h/2.*tan_e2
    T[0,2,2] = -h/2.*sec_e2
    T[1,0,0] = h/2.*h_pole*sec_e3 - (-k1 + h*h/2.*tan_e2)*tan_e
    T[1,0,1] = -h*tan_e2
    T[1,0,5] = -h*tan_e
    T[1,2,2] = (-k1 - h*h/2.*tan_e2)*tan_e - h/2.*h_pole*sec_e3
    T[1,2,3] = h*tan_e2
    T[2,0,2] = -h*tan_e2
    T[3,0,2] = -h*h_pole*sec_e3 +(-k1 + h*h*sec_e2)*tan_e
    T[3,0,3] = h*tan_e2
    T[3,1,2] = h*sec_e2
    T[3,2,5] = h*tan_e - h*phi/cos(e - phi)**2
    return R, T


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



