# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 21:04:59 2016

@author: Mitch
"""

import numpy as np
import pylab
from scipy.stats import norm
import matplotlib.pyplot as plt

def plotGPmodel(GP, xmin=-15, xmax=15, cov=False):
    mu = []
    xvals = []
    for i in range(1000):
        x = xmin + i*(xmax - xmin) / 1000.0
        xvals.append(x)
        mu.append(GP.predict(np.array(x,ndmin=2))[int(cov)])

    plt.scatter(xvals,mu)
    plt.show()

def vify(x, i):
    return np.array(x.iloc[i,:],ndmin=2)

def scat(y1,y2,fname=None):
    plt.scatter(range(len(y1)),y1,s=10,lw=0)
    plt.scatter(range(len(y1)),y2,s=10,lw=0,c=u'r')
    if(fname is not None):
        plt.savefig(fname,dpi=128)
    else:
        plt.show()

def BVplot(GP, f, xmin=-2.5, xmax=2.5, fname=None):
    xx = [xmin + (i/100.0)*(xmax-xmin) for i in range(100)]
    xxx = [xmin + (i/1000.0)*(xmax-xmin) for i in range(1000)]
    vec = np.reshape(np.array(xx),(100,1))
    vecx = np.reshape(np.array(xxx),(1000,1))
    pred,var = GP.predict(vec)
    var = np.reshape(np.diag(var),pred.shape)
    BVs = GP.BV
    fy = f(vecx)
    plt.errorbar(xx,pred,yerr=2*np.sqrt(var))#,c=u'b')
    plt.scatter(xxx,fy,s=10,lw=1,c=u'k')
    plt.scatter(BVs, [0 for x in BVs], s=30, c=u'r')
    plt.xlabel('x')
    plt.ylabel('y')
    if(fname is None):
        plt.show()
    else:
        plt.savefig(fname, dpi=128)

def errplot(y1,y2,fname=None,xlabel='Iteration number', ylabel='',mean=True):
    numiter = len(y1[0])
    if(not mean):
        plt.errorbar(range(numiter),np.median(y1,axis=0),yerr=np.std(y1,axis=0)/np.sqrt(len(y1)))
        plt.errorbar(range(numiter),np.median(y2,axis=0),yerr=np.std(y2,axis=0)/np.sqrt(len(y2)),c=u'r')
    else:
        plt.errorbar(range(numiter),np.mean(y1,axis=0),yerr=np.std(y1,axis=0)/np.sqrt(len(y1)))
        plt.errorbar(range(numiter),np.mean(y2,axis=0),yerr=np.std(y2,axis=0)/np.sqrt(len(y2)),c=u'r')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if(fname is not None):
        plt.savefig(fname,dpi=128)
    else:
        plt.show()

def GPheatmap(model, base, dim1, dim2, range1, range2, size=50, type='pred', ax=None, y_best=None, fname=None):
    # model: a GP that can predict
    # base: a vector to create a 2-D slice from
    # dim1: the index of the x-dimension to be sliced
    # dim2: the index of the y-dimension to be sliced
    # range1,range2: tuple (min,max) of plotting range for x,y
    # size: resolution of the heatmap
    # Var: untested, plots GP variance instead of predictions

    min1,max1 = range1
    min2,max2 = range2
    inc1 = (1.0 / size) * (max1 - min1)
    inc2 = (1.0 / size) * (max2 - min2)
    data = np.zeros(shape=(size*size,len(base)))
    for i in range(size):
        for j in range(size):
            base[dim1] = min1 + j * inc1
            base[dim2] = min2 + i * inc2
            data[i*size + j] = base

    pred,var = model.predict(data)
    if(type=='EI'):
        diff = pred - y_best
        var_vec = np.sqrt(np.reshape(np.diag(var),(size*size,1)))
        Z = diff / var_vec
        EI = diff * norm.cdf(Z) + var_vec * norm.pdf(Z)
        EI = np.nan_to_num(EI)
        mat = np.reshape(EI,(size,size))
    else:
        Mp = np.reshape(pred,(size,size))
        mat = Mp
        if(type=='var'):
            Mv = np.reshape(np.diag(var),(size,size))
            mat = Mv

    if(ax is None):
        plotter = pylab
        xt = pylab.xticks
        yt = pylab.yticks
    else:
        plotter = ax
        xt = ax.set_xticks
        yt = ax.set_yticks

    img = plotter.pcolor(mat)

    plt.colorbar(img,ax=ax)
    xlabels = [str(np.around(min1 + j * inc1 * int(size / 5.0),decimals=2)) for j in range(6)]
    ylabels = [str(np.around(min2 + i * inc2 * int(size / 5.0),decimals=2)) for i in range(6)]
    xt([i*int(size / 5.0) for i in range(6)], xlabels)
    yt([i*int(size / 5.0) for i in range(6)], ylabels)

    if(ax is not None):
        ax.set_xticklabels(xlabels)
        ax.set_yticklabels(ylabels)
        return

    if(fname is None):
        pylab.show()
    else:
        pylab.savefig(fname,dpi=128)

def regrets(y1,y2):
    y1 = np.array(y1)
    y2 = np.array(y2)
    maxs = np.reshape(np.max(np.concatenate((y1,y2),axis=1),axis=1),(y1.shape[0],1))
    return (maxs - y1, maxs - y2)

def rregrets(r1,r2,orig=True):
    if(orig):
        r1,r2 = regrets(r1,r2)
    sum1 = np.reshape(np.sum(r1,axis=1),(r1.shape[0],1))
    sum2 = np.reshape(np.sum(r2,axis=1),(r2.shape[0],1))
    rr1 = np.zeros(shape=r1.shape)
    rr2 = np.zeros(shape=r2.shape)
    rr1[:,[0]] = sum1
    rr2[:,[0]] = sum2
    for i in range(1,rr1.shape[1]):
        rr1[:,[i]] = rr1[:,[i-1]] - r1[:,[i-1]]
        rr2[:,[i]] = rr2[:,[i-1]] - r2[:,[i-1]]

    return rr1,rr2
