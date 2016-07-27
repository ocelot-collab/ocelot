__author__ = 'Sergey Tomin'

"""
author Mathias Vogt, DESY.
converted to Python by Sergey Tomin, E.XFEL, 2015.

"USAGE: tpk2i <type-string> <momentum/GeV> <k_i /m^-(i+1) or angle/deg or kick/mrad> \n");
"""

from ocelot import *
import numpy as np
from scipy.optimize import newton
def p1(a, x):
    return a[0] + a[1] * x

def p1p(a, x):
    return a[1]


def p1inv(a, y):
    return ( y - a[0] ) / a[1]


def p3(a, x):
    return a[0] + x*( a[1] + x*( a[2] + x* a[3] ) )

def p3p(a, x):
    return a[1] + x*( 2.0*a[2] + x* 3.0*a[3] )

def p3inv(a, y):
    tol=1.e-8
    x1 = p1inv( a, y )
    if a[2]==0.0 and a[3]==0.0:
        return x1
    func = lambda x: p3(a, x) - y
    x = newton(func, x1, tol=tol)
    return x

def p3inv_old(a, y): # // uses undamped Newton-Raphson
    tol=1.e-8
    x1 = p1inv( a, y )
    if a[2]==0.0 and a[3]==0.0:
        return x1
    x  = x1
    y0 = p3(a, x)
    y1 = p3p(a, x)
    r  = y0 - y
    r  =  r if r >= 0.0 else -r
    x  = x - (y0-y)/y1
    #//  printf("p3inv: y=%f y0=%f r=%f y1=%f x=%f x1=%f\n",y,y0,r,y1,x,x1);

    while( r > tol ):
        #print r, tol
        y0 = p3(    a, x )
        y1 = p3p(   a, x )
        r  = y0 - y
        r  = r if r >= 0.0 else -r
        x  = x - (y0-y)/y1
        #//    printf("p3inv: y=%f y0=%f r=%f y1=%f x=%f x1=%f\n",y,y0,r,y1,x,x1);
    return x


table ={
    #EL                             KM                      A0              A1             A2             A3
    "QA":     {"type": "Q", "EL": 0.18500, "KM": 1.000, "A": [         0.0,         0.3866,            0.0,           0.0 ]}, # //00 QA
    "QAX":    {"type": "Q", "EL": 0.29750, "KM": 1.000, "A": [     -0.0192,        0.42411,         7.2e-6,       -6.7e-6 ]}, # //01 QAX
    "TQA":    {"type": "Q", "EL": 0.27680, "KM": 0.976, "A": [       0.058,          0.128,            0.0,           0.0 ]}, # //02 TQA
    "TQA_U":  {"type": "Q", "EL": 0.27680, "KM": 0.976, "A": [      -0.058,          0.128,            0.0,           0.0 ]}, # //03 TQA_U
    "TQB":    {"type": "Q", "EL": 0.32860, "KM": 0.976, "A": [         0.0,      0.0826469,    2.87120e-06,  -1.72190e-08 ]}, # //04 TQB
    "TQB_U":  {"type": "Q", "EL": 0.32860, "KM": 0.976, "A": [         0.0,      0.0826469,   -2.87120e-06,  -1.72190e-08 ]}, # //05 TQB_U
    "TQD":    {"type": "Q", "EL": 0.30460, "KM": 1.000, "A": [7.736526e-02,   9.421388e-02,   2.641595e-05, -4.441297e-07 ]}, # //06 TQD
    "TQF":    {"type": "Q", "EL": 0.37280, "KM": 1.000, "A": [3.579171e-02,   7.023619e-02,  -9.550849e-07, -1.649632e-08 ]}, # //07 TQF
    "TQG":    {"type": "Q", "EL": 0.12716, "KM": 0.981, "A": [         0.0,       0.668062,            0.0,           0.0 ]}, # //08 TQG
    "TQG70":  {"type": "Q", "EL": 0.07600, "KM": 1.000, "A": [  -0.0587244,       0.260282,    3.61552e-05,  -1.45329e-07 ]}, # //09 TQG70 : new data for Leff
    "TQM":    {"type": "Q", "EL": 0.32860, "KM": 0.976, "A": [  0.04906110,     0.05349805,    1.39626e-05,  -3.44832e-08 ]}, # //10 TQM
    "QMN":    {"type": "Q", "EL": 0.33500, "KM": 1.000, "A": [        0.02,       0.069437,            0.0,           0.0 ]}, # //11 QMN
    "QC":     {"type": "Q", "EL": 1.02000, "KM": 1.000, "A": [    0.039929,       0.042478,     1.8915e-05,   -4.9054e-08 ]}, # //12 QC
    "QL":     {"type": "Q", "EL": 0.45000, "KM": 1.000, "A": [        0.00,     1.1666e-02,            0.0,           0.0 ]}, # //13 QL    : provisonal 25.04.2017 vogtm
    "QTS_E":  {"type": "Q", "EL": 0.09300, "KM": 0.759, "A": [        0.12,         0.1356,            0.0,           0.0 ]}, # //14 QTS_EXT
    "QTS_I":  {"type": "Q", "EL": 0.13100, "KM": 0.914, "A": [        0.12,         0.1356,            0.0,           0.0 ]}, # //15 QTS_INT
    "TDA":    {"type": "B", "EL": 0.44600, "KM": 1.000, "A": [   3.3070e-3,      2.4110e-3,      2.3900e-6,    -6.2599e-9 ]}, #//16 TDA
    "TDB":    {"type": "B", "EL": 0.22000, "KM": 1.000, "A": [   2.1930e-3,      1.3792e-3,      1.1715e-6,    -3.1038e-9 ]}, #//17 TDB
    "TDC":    {"type": "B", "EL": 1.20000, "KM": 1.000, "A": [  5.7368e-03,     4.1344e-03,     2.1649e-06,   -5.4703e-09 ]}, #//18 TDC
    "TDD":    {"type": "B", "EL": 0.50000, "KM": 1.000, "A": [      0.0031,         0.0042,        1.8e-06,      -2.4e-08 ]}, #//19 TDD
    "TDV":    {"type": "B", "EL": 0.44600, "KM": 1.000, "A": [      0.0011,         0.0026,            0.0,           0.0 ]}, #//20 TDV
    "TBZ":    {"type": "B", "EL": 1.00000, "KM": 1.000, "A": [   -0.004823,       0.002740,            0.0,           0.0 ]}, #//21 TBZ (FL2SEP) : ***
    "TSB":    {"type": "S", "EL": 0.21500, "KM": 1.000, "A": [         0.0,           5.82,            0.0,           0.0 ]}, #//22 TSB
    "TDDt":   {"type": "T", "EL": 0.50000, "KM": 1.000, "A": [         0.0,         0.0021,            0.0,           0.0 ]}, #//23 TDDt
    "TCA20":  {"type": "C", "EL": 0.14880, "KM": 1.000, "A": [    3.107e-3,       0.007404,      -3.204e-4,     -1.971e-3 ]}, #//24 TCA20  : D1IDUMP
    "TCA40":  {"type": "C", "EL": 0.10000, "KM": 1.000, "A": [         0.0,        0.09578,            0.0,           0.0 ]}, #//25 TCA40
    "TCA50":  {"type": "C", "EL": 0.10000, "KM": 1.000, "A": [         0.0,        0.08070,            0.0,           0.0 ]}, #//26 TCA50
    "TCA70":  {"type": "C", "EL": 0.10000, "KM": 1.000, "A": [         0.0,        0.06250,            0.0,           0.0 ]}, #//27 TCA70
    "TCA40S": {"type": "C", "EL": 0.02000, "KM": 1.000, "A": [         0.0,    5.577750e-2,            0.0,  -3.875895e-4 ]}, #//28 TCA40S  : ###
    "TCA50S": {"type": "C", "EL": 0.02000, "KM": 1.000, "A": [  0.000410877,      0.0493091,  -6.57403e-05,   -0.000318351]}, #//29 TCA50S  : ### JZ tpk2i ??
    "HGUN":   {"type": "C", "EL": 0.04500, "KM": 1.000, "A": [          0.0,  2.02222222e-4,           0.0,           0.0 ]}, #//30 HGUN
    "VGUN":   {"type": "C", "EL": 0.04500, "KM": 1.000, "A": [          0.0,  1.88888889e-4,           0.0,           0.0 ]}, #//31 VGUN
    "CV":     {"type": "C", "EL": 0.20000, "KM": 1.000, "A": [2.59106514e-4,  0.12904764512,1.030998316e-2,  -2.817702e-3 ]}, #//32 CV
    "CA":     {"type": "C", "EL": 0.15000, "KM": 1.000, "A": [          0.0,       0.000950,           0.0,           0.0 ]}, #//33 CA
    "CAX":    {"type": "C", "EL": 0.29750, "KM": 1.000, "A": [          0.0,       0.000568,           0.0,           0.0 ]}  #//34 CAX
    #"TQG70": {"EL": 0.11000, "KM": 1.000, "A": [ 0.0543891,     0.254164,   3.13809e-05,  -1.34808e-07  ]}, //09 TQG70 : ***
    #{"TDA":{"EL": 0.44600, "KM": 1.000,   "A":[    0.0011,       0.0026,           0.0,           0.0 ]}, //16 TDA : Johann ? test ??
    #"TDB":{ "EL": 0.22000, "KM": 1.000,   "A":[   0.0000,   1.4013e-03,           0.0,           0.0 ]}, //17 TDB : provisonal 25.04.2017 vogtm
    #"TBZ":   {"EL": 1.00000, "KM": 1.000, "A":[  -0.0038568,   0.00262531,   2.87309e-07,  -5.07466e-09 ]}, //21 TBZ (FL2SEP) :
    #"TCA40S":{"EL": 0.02000, "KM":1.000, "A": [1.965654e-03, 3.136079e-01,  1.719241e-03,  -6.394868e-04]}, //28 TCA40S
    #"TCA50S":{"EL":0.02000, "KM": 1.000,  "A":[        0.0,       0.0981,           0.0,           0.0 ]}, //29 TCA50S
    #"TCA50S":{"EL":0.02000, "KM": 1.000,  "A":[        0.0,    0.0493091,           0.0,   -0.000318351]}, //29 TCA50S  : ### JZ tpi2k ??
}



# ***:  changed by jzemella's calib of MEA data (look into logbook doc magnets)
# ###:  (unipolar data  mirrored to get bipolar data)

def tpi2k(id, p, dI):
    p2brho = 10./2.99792458
    deg = pi/180.0
    if id not in table.keys():
        print( "type identifier is unknown")
        print( "known types:")
        print( list(table.keys())[:10])
        print( list(table.keys())[10:20])
        print( list(table.keys())[20:30])
        print( list(table.keys())[35:])
        return None
    magnet_type = table[id]["type"]
    a = table[id]["A"]
    el = table[id]["EL"]
    km = table[id]["KM"]
    dk0 = p3(a, dI ) * km/(p*p2brho)
    dk = 0.
    if magnet_type == "S" or magnet_type == "Q":
        #case kQ:
        dk = dk0
    elif magnet_type == "B":
        phi = 2.0*np.arcsin(0.5*el * dk0)
        dk = phi/deg
    elif magnet_type ==  "C":
        phi = 2.0*np.arcsin(0.5*el * dk0)
        dk = phi*1.e3
    elif magnet_type == "T":
        phi = 2.0*np.arcsin(0.5*el * dk0)
        dk = phi*1.e3
    return dk


def tpk2i(id, p, k):
    p2brho = 10./2.99792458
    deg = pi/180.0
    if id not in table.keys():
        print( "type identifier is unknown")
        print( "known types:")
        print( list(table.keys())[:10])
        print( list(table.keys())[10:20])
        print( list(table.keys())[20:30])
        print( list(table.keys())[35:])
        return None
    magnet_type = table[id]["type"]
    a = table[id]["A"]
    el = table[id]["EL"]
    km = table[id]["KM"]
    dk = 0.
    if magnet_type == "S" or magnet_type == "Q":
        #case kQ:
        dk = k
    elif magnet_type == "B":
        phi=k*deg
        dk=2.0*np.sin(0.5*phi)/el
    elif magnet_type ==  "C":
        phi=k*1.e-3
        dk=2.0*np.sin(0.5*phi)/el
    elif magnet_type == "T":
        phi=k*1.e-3
        dk=2.0*np.sin(0.5*phi)/el
    di = p3inv( a, p*p2brho*dk/km )
    return di

if __name__ == "__main__":
    I = 50.
    k2 = tpi2k(id = "TSB", p = 0.700, dI=I)
    I = tpk2i(id = "TSB", p = 0.700, k=k2)
    print ("p = 0.700, sextupole -> tpi2k: current = ", I, "A;", " strength = ", k2, "1/m**3")
    print ("p = 0.700, sextupole -> tpk2i: strength = ", k2, "1/m**3;", " current = ", I, "A")
    #print (tpk2i(id = "TCA40", p = 0.7, k=-12.3061297779))
