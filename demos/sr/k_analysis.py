__author__ = 'Sergey Tomin'

import numpy as np


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]


def data_analysis(e_ph, flux, method="least"):

    if method == "least":
        coeffs = np.polyfit(x=e_ph, y=flux, deg=11)
        polynom = np.poly1d(coeffs)


        x = np.linspace(e_ph[0], e_ph[-1], num=100)
        pd = np.polyder(polynom, m=1)
        indx = np.argmax(np.abs(pd(x)))
        eph_c = x[indx]

        pd2 = np.polyder(polynom, m=2)
        p2_roots = np.roots(pd2)
        p2_roots = p2_roots[p2_roots[:].imag == 0]
        p2_roots = np.real(p2_roots)
        Eph_fin = find_nearest(p2_roots,eph_c)
        return Eph_fin, polynom

    elif method == "new method":
        pass

        #plt.plot(Etotal, total, "ro")
        #plt.plot(x, polynom(x))
