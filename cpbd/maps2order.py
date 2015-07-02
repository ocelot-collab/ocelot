__author__ = 'Sergey Tomin'

from numpy import cos, sin, sqrt

"""
differential equation:

t_nnn'' + kx**2*t_nnn = f_nnn

here index nnn means (1,2,3,4,5,6) = (x, x', y, y', s, delta)

# h = 1/r0 - curvature
# k1 = 1/Bro*dBy/dx - gradient
# k2 = 1/Bro*d**2By/dx**2 - sextupole component
# h1 - derivative of h
# h11 - derivative of h1
# cx1, sx1, cy1 ... - derivative of cx, sx, cy, ... and so on
# cx = cos(kx*s)
# sx = sin(kx*s)/kx
# dx = h/kx2*(1. - cx)
# cy = cos(ky*s)
# sy = sin(ky*s)/ky

# Green's function for X is Gx(t, tau) = 1/kx*sin(kx*(t - tau))

f111 = -(h**3 + K2 - 2*h*ky**2)*cx**2 + 1./2.*h*kx**4*sx*2 - kx**3*cx*sx*h'
f112 = -h * kx**2 * cx * sx - 2 * kx * (h**3 + K2 - 2 * h * ky**2) * cx * sx + (cx**2 - kx**2 * sx**2)*h'
f116 = (2*h^2 -ky^2)* cx + 2 *(h^3+K2-2 *h *ky^2) *cx* dx-h^2*kx^2 *sx^2 + (h *kx *cx *sx-dx* kx^2*sx)* h'
f122 = 1/2* h* cx^2 + h^3* (-1-K2/h^3+(2 ky^2)/h^2) *sx^2 + sx *sx* h'
f126 = h^2 *(2-ky^2/h^2) sx+h^2 *cx* sx + 2 h^3 (-1-K2/h^3+(2 ky^2)/h^2) *dx* sx + (cx dx +h sx^2) h'
F166 = 2*dx*h**2 - dx**2*h**3 - dx*(dx*k2 + k1**2) + h*(-1 + 2*dx**2*k1**2) + 1./2.*h*(dx1)**2 + dx*dx1*h1
F133 = -(1./2.)*h*(cy1)**2 + cy*cy1*h1 + 1./2.*cy**2*(2*k2 - h*k1**2 + h11)
F134 = cy1*(sy*h1 - h*sy1) + cy*(h1*sy1 + sy*(2*k2 - h*k1**2 + h11))
f346 = sy*h1*sy1 + 1./2.*(-h*(sy1)**2 + sy**2*(2*k2 - h*k1**2 + h11))

# Green's function for Y is Gy(t, tau) = 1/ky*sin(ky*(t - tau))

f313 = cx1*(h*cy1 - cy*h1) + cx*(2*cy*(k2 - h*k1**2) + cy1*h1)
f314 = cx1*(-sy*h1 + h*sy1) + cx*(2*(k2 - h*k1**2)*sy + h1*sy1)
f323 = cy1*(sx*h1 + h*sx1) + cy*(2*(k2 - h*k1**2)*sx - h1*sx1)
f324 = 2*(k2 - h*k1**2)*sx*sy + h*sx1*sy1 + h1*(-sy*sx1 + sx*sy1)
f336 = cy1*(h*dx1 + dx*h1)+cy*(2*dx*k2+k1**2 - 2*dx*h*k1**2 - dx1*h1)
f346 = (2*dx*k2 + k1**2 - 2*dx*h*k1**2)*sy + dx*h1*sy1 + dx1*(-sy*h1+h*sy1)

Integration:

I111 = G * cx**2 = ((2. + cx)*dx)/(3*h)
I122 = G * sx**2 = dx**2/3./h**2
I112 = G * cx*sx = sx*dx/(3.*h)
I11  = G * cx    = s*sx/2.
I116 = G * cx*dx = h/kx2*(cx - cx**2) = h/kx2*(I11 - I111)
I12  = G * sx    = 1./2./kx2*(sx - s*cx/2.)
I126 = G * sx*dx = h/kx2*(sx - sx*cx) = h/kx2*(I12 - I112)
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
    kx = sqrt(kx2 + 0.j)
    ky = sqrt(ky2 + 0.j)
    cx = cos(kx*L)
    sx = sin(kx*L)/kx
    cy = cos(ky*L)
    sy = sin(ky*L)/ky
    dx = h/kx2*(1. - cx)
    # Integrals
    I111 = ((2. + cx)*dx)/(3*h)
    I122 = dx**2/3./h**2
    I112 = sx*dx/(3.*h)
    I11  = L*sx/2.
    I116 = h/kx2*(I11 - I111)
    I12  = 0.5/kx2*(sx - L*cx/2.)
    I126 = h/kx2*(I12 - I112)
    coef1 = h3 + k2 - 2*h*ky2

    t111 = -coef1*I111 + 0.5*h*kx2*kx2*I122
    t112 = -2.*kx*coef1*I112 - h*kx2*I112
    t116 = 2.*coef1*I116 + (2*h2 - ky2)*I11 - h2*kx2*I122
    t122 = 0.5*h*I111 + coef1*I122
    t126 = h2*(2 - ky2/h2)*I12 + h2*I112+ 2*coef1*I126



def map(l, angle, k1, k2, k3):
    pass












