__author__ = 'Sergey Tomin'

import numpy as np


def derivat(x, y):

    N=len(x)-1
    deriv1 = 0.
    deriv2 = 0.
    for i in range(2):
        n = i*(N-3)
        h0 = x[n+1] - x[n]
        h1 = x[n+2] - x[n]
        h2 = x[n+3] - x[n]
        df = (y[n + 1] - y[n])/h0
        F0 = (y[n + 2] - y[n])/h1 - df
        F1 = (y[n + 3] - y[n])/h2 - df
        a11 = h1*h1 - h0*h0
        a12 = h1 - h0
        a21 = h2*h2 - h0*h0
        a22 = h2 - h0
        det = a11*a22 - a12*a21
        a0 = (a22*F0 - a12*F1)/det
        b0 = (-a21*F0 + a11*F1)/det
        c0 = df - a0*h0*h0 - b0*h0

        derivative = (3.*a0*(x[N*i]-x[n*i])*(x[N*i]-x[n*i]) + 2.*b0*(x[N*i]-x[n*i]))*i + c0
        if(i==0):
            deriv1 = derivative
        else:
            deriv2 = derivative

    return deriv1 , deriv2


def moment(x, y):

    n = len(x)
    N =len(x)-1

    Alfa = np.zeros(n)
    Beta = np.zeros(n)
    F = np.zeros(n)
    h0 = x[1] - x[0]

    deriv1, deriv2 = derivat(x, y)
    #deriv1 = 100

    #deriv1 = 0
    Alfa[0] = 0.
    Beta[0] = 0.
    f = 6.*((y[1]-y[0])/h0-deriv1)
    F[0] = f
    Alfa[1] = -h0/(2.*h0)
    Beta[1] = f/(2.*h0)

    for i in range(1,N):
        h1 = x[i+1] - x[i]
        F[i] = 6.*((y[i+1] - y[i])/h1 - (y[i] - y[i-1])/h0)
        k = h0*Alfa[i] + 2.*(h1+h0)
        Alfa[i+1] = -h1/k
        Beta[i+1] = (F[i] - h0*Beta[i])/k
        h0 = h1

    F[N] = 6.*(deriv2 - (y[N]-y[N-1])/h0)
    M = np.zeros(len(x))
    M[N] = (F[N] - h0*Beta[N])/(h0*Alfa[N] + 2.*h0)
    #M[N-1:0:-1] = Alfa[N:1:-1]*M[N:1:-1] + Beta[N:1:-1]
    for j in range(N-1, -1, -1):
        M[j] = Alfa[j+1]*M[j+1] + Beta[j+1]
        #print j, M[j]
    #M[0] = -0.830334664514

    return M



def cspline_coef(x, y):
    n = len(x)
    #N = len(arZ) - 1
    a = np.zeros(n-1)
    b = np.zeros(n-1)
    c = np.zeros(n-1)
    d = np.zeros(n-1)
    z = np.zeros(n-1)
    if(n<4):
        print( "too short massive. Length must be > 3")
        return 0
    M = moment(x, y)
    #print M
    for i in range(n - 1):
        h = x[i+1]-x[i]
        M0 = M[i]
        a[i] = (M[i+1] - M0)/(6.*h)
        b[i] = M0/2.
        c[i] = (y[i+1] - y[i])/h - M[i+1]*h/6. -  M0*h/3.
        d[i] = y[i]
        z[i] = x[i]
    return a,b,c,d,z

def cinterp(x, y, x_new):
    Y_new = []
    A,B,C,D,Z = cspline_coef(x, y)
    for s in x_new:
        #if s<=Z[0]:
        #    x0 = s - Z[0]
        #    y_new = A[0]*x0**3 + B[0]*x0*x0 + C[0]*x0 +D[0]
        #elif s>Z[-1]:
        #    x0 = s - Z[-1]
        #    y_new = A[-1]*x0**3 + B[-1]*x0*x0 + C[-1]*x0 +D[-1]
        #else:
        ind = np.where(s>=Z)[0][-1]
        #print ind, s , Z[ind]
        x0 = s - Z[ind]
        y_new = ((A[ind]*x0 + B[ind])*x0 + C[ind])*x0 +D[ind]


        Y_new.append(y_new)
    return np.array(Y_new)



if __name__ == "__main__":
    from pylab import *
    from time import time
    f = lambda x: sin(x)
    x = np.linspace(0., 5, 100)
    y = np.array([f(xi) for xi  in x])
    start = time()
    A,B,C,D,Z = cspline_coef(x, y)
    print( time() - start)
    #print D
    x_new = np.linspace(-1.1, 1.1, 1000)
    y_new = cinterp(x, y, x_new)
    plot(x, y, "r.-", x_new, y_new,"b")
    show()
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
            b2 += h*(d*d + h*(c*d + h*(1./3.* (c*c + 2*b*d) + h*(0.5*(b*c + a*d) + h*(0.2*(b**2 + 2*a*c) + h*(1./3.*a*b + (a*a*h)/7.))))))
            beta2.append( b2)

        return array(beta2)
    y2 = integ_beta2(x, y)
    print( len(y2), len(x), len(y) )
    plot(x, x/2. - 1./4.*sin(2*x), x, y2)
    show()
