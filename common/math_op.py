'''
statistical analysis functions, fitting, optimization and the like
'''

import numpy as np
from numpy import cos, sin, tan, sqrt, log, exp, sum

def peaks(x, y, n=0):
    '''
    
    '''
    maxs = {}

    if len((np.where(y == y.max()))[0]) == 1:
        for i in np.arange(1, len(x)-1):
            if (y[i] - y[i-1]) > 0 and (y[i+1] - y[i]) < 0:
                maxs[y[i]] = x[i]
    else:
        for i in np.arange(2, len(x)-1):
            if bool((y[i-1] - y[i-2]) > 0 and (y[i] - y[i-1]) == 0 and (y[i+1] - y[i]) < 0):
                maxs[y[i]] = x[i]

    vals  = sorted(maxs.keys())
    f1 = []
    f2 = []
    for v in reversed(vals):
        f1.append(maxs[v])
        f2.append(v)
            
    if n>0:
        return np.array(f1[0:n]), np.array(f2[0:n])
    else:
        return np.array(f1), np.array(f2)

def gs_search(f, bracket, tol=1.e-6, nmax=50):
    '''
    golden section search
    '''
    n = 0
    a,b,c = bracket
    cur_tol = c-a
    
    w = 0.38197
    
    while True:
        if n > nmax: break
        if cur_tol < tol: break
        
        if c-b > b-a:
            x = b + (c-b)*w
            if f(b)<f(x):
                c = x
            else:
                a = b
                b = x
        else:
            x = b - (b-a)*w
            if f(b)<f(x):
                a = x
            else:
                c = b
                b = x
        
        n += 1
        cur_tol = c-a
    
    return np.array([a,b,c]), cur_tol

def fit_gauss_2d(x,y,F):
    
    n1 = len(x)
    n2 = len(y)
    
    mu1 = 0
    mu2 = 0
    
    for i in range( n1 ): mu1 += x[i]*sum(F[i,:])
    for i in range( n2 ): mu2 += y[i]*sum(F[:,i])
    mu2 /= sum(F[:,:])
    mu1 /= sum(F[:,:])
    
    sig1 = 0
    
    for i in range( n1 ): sig1 += (x[i]-mu1)**2*sum(F[i,:])
    sig1 /= sum(F[:,:])
    
    sig2 = 0
    
    for i in range( n2 ): sig2 += (y[i]-mu2)**2*sum(F[:,i])
    sig2 /= sum(F[:,:])

    sig12 = 0

    for i in range( n1 ):
        for j in range( n2 ):
            sig12 += (x[i]-mu1)*(y[j]-mu2)*F[i,j]
    
    sig12 /= sum(F[:,:])

    rho = sig12 / sqrt(sig1*sig2)
    
    return mu1, mu2, sqrt(sig1), sqrt(sig2), rho

def fit_gauss_1d(x,F):
    
    n1 = len(x)
    
    mu1 = 0
    
    for i in range( n1 ): mu1 += x[i]*F[i]
    mu1 /= sum(F[:])
    
    sig1 = 0
    
    for i in range( n1 ): sig1 += (x[i]-mu1)**2*F[i]
    sig1 /= sum(F[:])
    
    return mu1, sqrt(sig1)

def fwhm(x,F):
    m = np.max(F) / 2.0
    ups = []
    downs = []
    for i in range(len(x)-1):
        if F[i] <  m and F[i+1] > m:
            ups.append(i)
        if F[i] >=  m and F[i+1] < m:
            downs.append(i)
            
    #print ups, downs
    return x[downs[-1]] - x[ups[0]]

def stats(outputs):
    
    ''' return mean, std, median and extreme (farthest from mean) of a time series '''
    
    omean = np.zeros_like(outputs[0])

    for o in outputs:
        omean += np.array(o)
    omean = omean / len(outputs)
        
    difm = np.linalg.norm(omean, 2)
    difw = 0.0
    imed = 0
    iworst = 0
    std = np.zeros_like(outputs[0])
    
    for i in np.arange( len(outputs) ):
        
        difn = np.linalg.norm(omean - np.array(outputs[i]), 2)
        
        if difn < difm:
            imed = i
            difm = difn

        if difn > difw:
            iworst = i
            difw = difn
            
        std += (omean - np.array(outputs[i]))**2
                
    std = np.sqrt(std / len(outputs)) 
         
    #print 'median id=', imed
    return omean , std, outputs[imed], outputs[iworst], imed, iworst


def find_saturation(power, z, n_smooth=5):
    p = np.diff(np.log10(power))

    u = np.convolve(p, np.ones(n_smooth) / float(n_smooth), mode='same')
    um = np.max(u)
    
    ii = 0
    
    for i in range(len(u)):
        if u[i] < 0.0 * um and z[i] > 10: 
            ii = i
            break
    
    #plt.plot(g.z[1:], u, lw=3)
    #plt.plot(g.z[ii+1], p[ii], 'rd')

    #plt.plot(g.z, power, lw=3)
    #plt.plot(z[ii+1], np.log10(power[ii]), 'rd')

    return z[ii+1], ii+1

def find_nearest(array, value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
    
def n_moment(x, counts, c, n):
    if np.sum(counts)==0:
        return 0
    else:
        return (np.sum((x-c)**n*counts) / np.sum(counts))**(1./n)
        
def std_moment(x, counts):
    mean=n_moment(x, counts, 0, 1)
    return n_moment(x, counts, mean, 2)
    
def bin_array(array,bin_size):
    #bining the array by averaging values within bins with size bin_size pixels
    if bin_size > len(array):
        return np.mean(array)
    elif bin_size == 1:
        return array
    else:
        new_shape = (array.shape[0] // bin_size) * bin_size
        array_av = array[:new_shape]
        array_av = array_av.reshape(int(new_shape/bin_size), bin_size)
        array_av = np.mean(array_av, axis=1)
        return array_av

def bin_scale(scale,bin_size):
    if bin_size > len(scale):
        return np.array([0])
    elif bin_size == 1:
        return scale
    else:
        hbin = np.int(bin_size/2) #halfbin (to pick bin centers)
        new_shape = (scale.shape[0] // bin_size) * bin_size
        new_scale = scale[hbin : new_shape+hbin] [ :: bin_size]
        return new_scale

def index_of(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx