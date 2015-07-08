__author__ = 'Sergey Tomin'

from numpy import cos, sin, sqrt

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
beta = k2/h**3

# Green's function for X is Gx(t, tau) = 1/kx*sin(kx*(t - tau))

f111 = -(h**3 + K2 - 2*h*ky**2)*cx**2 + 1./2.*h*kx**4*sx*2 - kx**3*cx*sx*h'
f112 = -h * kx**2 * cx * sx - 2 * kx * (h**3 + K2 - 2 * h * ky**2) * cx * sx + (cx**2 - kx**2 * sx**2)*h'
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


def t_nnn(L, angle, k1, k2, k3):
    h = 0.
    if L >0:
        h = angle/L
    else:
        exit("error: l <= 0")

    h2 = h*h
    h3 = h2*h
    kx2 = (k1 + h*h)
    ky2 = -k1
    kx4 = kx2*kx2
    ky4 = ky2*ky2
    kx = sqrt(kx2 + 0.j)
    ky = sqrt(ky2 + 0.j)
    cx = cos(kx*L)
    sx = sin(kx*L)/kx
    cy = cos(ky*L)
    sy = sin(ky*L)/ky
    sx2 = sx*sx
    sy2 = sy*sy

    dx = h/kx2*(1. - cx)

    # Integrals
    denom = kx2 - 4.*ky2
    I111 = 1./3.*(dx/h + sx2)                               #  I111 = Gx * cx**2
    I122 = dx**2/3./h2                                      #  I122 = Gx * sx**2
    I112 = sx*dx/(3.*h)                                     #  I112 = Gx * cx*sx
    I11  = L*sx/2.                                          #  I11  = Gx * cx
    I116 = h/kx2*(I11 - I111)                               #  I116 = Gx * cx*dx
    I12  = 0.5/kx2*(sx - L*cx/2.)                           #  I12  = Gx * sx
    I126 = h/kx2*(I12 - I112)                               #  I126 = Gx * sx*dx
    I10  = dx/h                                             #  I10  = Gx
    I16  = h/kx2*(dx/h - L*sx/2.)                           #  I16  = Gx * dx
    I166 = h2/kx4*(I10 - 2*I11 + I111)                      #  I166 = Gx * dx**2
    I144 = (sy2 - 2.*dx/h)/denom                            #  I144 = Gx * sy**2
    I133 = dx/h + ky2*(2.*dx/h - sy2)/denom                 #  I133 = Gx * cy**2
    I134 = (sy*cy - sx)/denom                               #  I134 = Gx * cy*sy
    I313 = (kx2*cy*dx/h - 2*ky2*sx*sy)/denom                #  I313 = Gy * cx*cy
    I324 = (kx2*cy*dx/h - 2*ky2*sx*sy)/denom                #  I324 = Gy * sx*sy
    I314 = (2*cy*sx - (1 + cx)*sy)/denom                    #  I314 = Gy * cx*sy
    I323 = ((1 - 2*ky2*dx/h)*sy - cy*sx)/ denom             #  I323 = Gy * sx*cy
    I33  = L*sy/2.                                          #  I33  = Gy * cy
    I336 = h/kx2*(I33 - I313)                               #  I336 = Gy * dx*cy
    I34  = (sy - L*cy)/(2*ky2)                              #  I34  = Gy * sy
    I346 = h/kx2*(I34 - I314)                               #  I346 = Gy * dx*sy

    #derivative of Integrals
    I211 = sx/3.*(1. + 2.*cx)


    coef1 = 2*ky2*h - h3 - k2/2

    t111 =    coef1*I111 + 0.5*h*kx4*I122
    t112 = 2.*coef1*I112 - h*kx2*I112
    t116 = 2.*coef1*I116 + (2*h2 - ky2)*I11 - h2*kx2*I122
    t122 =    coef1*I122 + 0.5*h*I111
    t126 = 2.*coef1*I126 + (2*h2 - ky2)*I12 + h2*I112
    t166 =    coef1*I166 - dx + I16*(2.*h2 - ky2) + 0.5*h3*I122
    t133 = (k2*I133 - ky2*dx)*0.5
    t134 = I134*k2
    t144 = 0.5*(k2*I144 - dx)

    coef2 = 2*(K2 - ky2*h)

    t313 = coef2*I313 + h*kx2*ky2*I324
    t314 = coef2*I314 - h*kx2*I323
    t323 = coef2*I323 - h*ky2*I314
    t324 = coef2*I324 + h*I313
    t336 = coef2*I336 + ky2*I33 - h2*ky2*I324
    t346 = coef2*I346 + h2*I323 + ky2*I34

def map(l, angle, k1, k2, k3):
    pass












