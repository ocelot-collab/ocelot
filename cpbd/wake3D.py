'''
Created on 17.05.2016
@author: Igor Zagorodnov
'''

from ocelot import *
from ocelot.adaptors import *
from ocelot.adaptors.astra2ocelot import *

def Filter(x,NF):
    Ns=x.shape[0]
    for i in range(NF):
        x[1:Ns]=(x[1:Ns]+x[0:Ns-1])*0.5
        x[0:Ns-1]=(x[1:Ns]+x[0:Ns-1])*0.5
    return x

def Der(x,y):
    #numerical derivative
    n=x.shape[0]
    dy = np.zeros(n)
    dy[1:n-1]=(y[2:n]-y[0:n-2])/(x[2:n]-x[0:n-2])
    dy[0]=(y[1]-y[0])/(x[1]-x[0])
    dy[n-1]=(y[n-1]-y[n-2])/(x[n-1]-x[n-2])
    return dy

def Int1(x,y):
    n=x.shape[0]
    Y=np.zeros(n)
    for i in range (1,n):
        Y[i]=Y[i-1]+0.5*(y(i)+y(i-1))*(x(i)-x(i-1))
    return Y

def Int1h(h,y):
    n=y.shape[0]	
    Y=np.zeros(n)
    # slow, switch to vector operations to be done
    for i in range(1,n):
        Y[i]=Y[i-1]+0.5*(y[i]+y[i-1])
    Y=Y*h
    return Y

def s2current(P0,q,Ns,NF,v):
    """ I=s2current(P0,q,Ns,NF)
        P0 - s-vector
        q - charge-vector
        Ns- number of sampling points
        NF- filter order 
        v - mean velocity
        """
    s0=np.min(P0) 
    s1=np.max(P0)
    NF2=int(np.floor(NF/2))
    Ns=Ns+2*NF2
    Np=P0.shape[0]
    ds=(s1-s0)/(Ns-2-2*NF2)
    s=s0+np.arange(-NF2,Ns-NF2)*ds

    # here we need a fast 1D linear interpolation of charges on the grid
    # in sc_old.py we use a fast 3D "near-point" interpolation
    # we need a stand-alone module with 1D,2D,3D parricles-to-grid functions
    Ip=(P0-s0)/ds
    I0=np.floor(Ip)
    dI0=Ip-I0
    I0=I0+NF2
    Ro=np.zeros(Ns)
    # slow, switch to vector operations to be done
    for i in range(Np):
        i0=int(I0[i])
        di0=dI0[i]
        Ro[i0]=Ro[i0]+(1-di0)*q[i]
        Ro[i0+1]=Ro[i0+1]+di0*q[i]
    if NF>0:
        Filter(Ro,NF)
    I=np.zeros([Ns,2])
    I[:,0]=s
    I[:,1]=Ro*v/ds
    return I

def convolution(xu,u,xw,w):
    #convolution of equally spaced functions
    hx=xu[1]-xu[0]
    wc=np.convolve(u,w)*hx
    nw=w.shape[0]
    nu=u.shape[0]
    x0=xu[0]+xw[0]
    xc=x0+np.arange(nw+nu)*hx
    return xc,wc

def wakeconvolution(xb,bunch,xw,wake):
    #convolution of unequally spaced functions
    #bunch defines the parameters
    nb=xb.shape[0]
    xwi=np.zeros(nb)
    xwi=xb-xb[0]
    wake1=np.interp(xwi,xw,wake,0,0)
    wake1[0]=wake1[0]*0.5
    xW,Wake=convolution(xb,bunch,xwi,wake1)
    return xW[0:nb],Wake[0:nb]

def AddWake (I,T):
    """[x, W] =AddWake (I,T)
        T - wake table in V/C, W in V
        (R,L,Cinv,nm,W0,N0,W1,N1)"""
    R,L,Cinv,nm,W0,N0,W1,N1=T
    c=299792458
    x=I[:,0]
    bunch=I[:,1]
    if L!=0 or N1>0: 
        d1_bunch=Der(x,bunch)
    nb=x.shape[0]
    W=np.zeros(nb) 
    if N0>0: 
        x,ww=wakeconvolution(x,bunch,W0[:,0],W0[:,1])  
        W=W-ww[0:nb]/c
    if N1>0: 
        x,ww=wakeconvolution(x,d1_bunch,W1[:,0],W1[:,1])
        W=W-ww[0:nb]
    if R!=0: 
        W=W-bunch*R
    if L!=0: 
        W=W-d1_bunch*L*c
    if Cinv!=0:
      int_bunch=Int1(x,bunch) 
      W=W-int_bunch*Cinv/c 
    return x,W

def LoadWakeTable (wakeFile):
    # return T- table of wakes coefs, H- matrix of the coefs place in T
    W=np.loadtxt(wakeFile)
    # head format %Nt 0 %N0 N1 %R L %C nm
    H=np.zeros([5,5])
    Nt=int(W[0,0])
    T=[]
    ind=0
    for i in range(Nt):
        ind=ind+1
        N0=int(W[ind,0])
        N1=int(W[ind,1])
        R=W[ind+1,0]
        L=W[ind+1,1]
        Cinv=W[ind+2,0]
        nm=int(W[ind+2,1])
        n=np.floor(nm/10) 
        m=nm-n*10
        H[n,m]=i
        ind=ind+2
        if N0>0:
            W0=np.zeros([N0,2]) 
            W0[0:N0,:]=W[ind+1:ind+N0+1,:] 
            ind=ind+N0 
        else: 
            W0=0
        if N1>0:
            W1=np.zeros([N1,2])
            W1[0:N1,:]=W[ind+1:ind+N1+1,:]
            ind=ind+N1 
        else: 
            W1=0
        T=T+[(R,L,Cinv,nm,W0,N0,W1,N1)]
    return (T,H)

def AddTotalWake (X,Y,Z,q,TH,Ns,NF):
    #function [Px Py Pz I00]=AddTotalWake (P,q,wakeFile,Ns,NF)
    T,H=TH
    c=299792458;
    #Z=-Z
    Np=X.shape[0]
    X2=X**2
    Y2=Y**2 
    XY=X*Y
    #generalized currents;
    I00=s2current(Z,q,Ns,NF,c)
    Nw=I00.shape[0]
    if (H[0,2]>0)or(H[2,3]>0)or(H[2,4]>0):
        qn=q*Y
        I01=s2current(Z,qn,Ns,NF,c)
    if (H[0,1]>0)or(H[1,3]>0)or(H[1,4]>0):
        qn=q*X
        I10=s2current(Z,qn,Ns,NF,c)
    if H[1,2]>0:
        qn=q*XY
        I11=s2current(Z,qn,Ns,NF,c)
    if H[1,1]>0:
        qn=q*(X2-Y2)
        I20_02=s2current(Z,qn,Ns,NF,c)
    #longitudinal wake
    #mn=0
    x, Wz =AddWake (I00,T[int(H[0,0])])
    if H[0,1]>0:
        x, w =AddWake(I10,T[int(H[0,1])]) 
        Wz=Wz+w;
    if H[0,2]>0:
        x, w =AddWake(I01,T[int(H[0,2])]) 
        Wz=Wz+w
    if H[1,1]>0:
        x, w =AddWake(I20_02,T[int(H[1,1])]) 
        Wz=Wz+w;
    if H[1,2]>0:
        x, w =AddWake(I11,T[int(H[1,2])]) 
        Wz=Wz+2*w
    Pz=np.interp(Z,x,Wz,0,0)
    Py=np.zeros(Np)
    Px=np.zeros(Np)
    #mn=01
    Wz[0:Nw]=0
    Wy=np.zeros(Nw)
    if H[0,4]>0:
        x, w =AddWake(I00,T[int(H[0,4])])
        Wz=Wz+w
        Wy=Wy+w
    if H[1,4]>0:
        x, w =AddWake(I10,T[int(H[1,4])]) 
        Wz=Wz+2*w
        Wy=Wy+2*w
    if H[2,4]>0:
        x, w =AddWake(I01,T[int(H[2,4])])
        Wz=Wz+2*w
        Wy=Wy+2*w
    Pz=Pz+np.interp(Z,x,Wz,0,0)*Y
    h=x[1]-x[0]
    Wy=-Int1h(h,Wy)
    Py=Py+np.interp(Z,x,Wy,0,0)
    #mn=10
    Wz[0:Nw]=0
    Wx=np.zeros(Nw)
    if H[0,3]>0:
        x, w =AddWake(I00,T[int(H[0,3])]) 
        Wz=Wz+w
        Wx=Wx+w
    if H[1,3]>0:
        x, w =AddWake(I10,T[int(H[1,3])]) 
        Wz=Wz+2*w
        Wx=Wx+2*w
    if H[2,3]>0:
        x, w =AddWake(I01,T[int(H[2,3])]) 
        Wz=Wz+2*w
        Wx=Wx+2*w
    Wx=-Int1h(h,Wx)
    Pz=Pz+np.interp(Z,x,Wz,0,0)*X
    Px=Px+np.interp(Z,x,Wx,0,0)
    #mn=11
    if H[3,4]>0:
        x, w =AddWake(I00,T[int(H[3,4])]); 
        Wx=-2*Int1h(h,w)
        p=np.interp(Z,x,Wx,0,0)
        Px=Px+p*Y
        Py=Py+p*X
        Pz=Pz+2*np.interp(Z,x,w,0,0)*XY
    #mn=02,20
    if H[3,3]>0:
        x, w =AddWake(I00,T[int(H[3,3])]) 
        Pz=Pz+np.interp(Z,x,w,0,0)*(X2-Y2)
        Wx=-2*Int1h(h,w)
        p=np.interp(Z,x,Wx,0,0)
        Px=Px+p*X
        Py=Py-p*Y
    I00[:,0]=-I00[:,0]
    #Z=-Z
    return Px,Py,Pz,I00
