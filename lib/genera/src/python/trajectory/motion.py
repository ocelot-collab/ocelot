__author__ = 'Sergey Tomin'
import numpy as np
from ctypes import *
#from xframework.cpbd.beam import Particle
def Py2C(array):
    arr_type =  c_double*len(array)
    c_array = arr_type(*array)
    return c_array


class Motion:
    def __init__(self, N = 0):
        self.N = N
        #if not bRough:
            #self.N = 3*N - 1

        #print "Motion N = ", N
        self.memory_motion = np.zeros(self.N*11)

        self.X = self.memory_motion[0:self.N]
        self.Y = self.memory_motion[self.N:2*self.N]
        self.Z = self.memory_motion[2*self.N:3*self.N]

        self.Xbeta = self.memory_motion[3*self.N:4*self.N]
        self.Ybeta = self.memory_motion[4*self.N:5*self.N]
        self.Zbeta = self.memory_motion[5*self.N:6*self.N]

        self.Bx = self.memory_motion[6*self.N:7*self.N]
        self.By = self.memory_motion[7*self.N:8*self.N]
        self.Bz = self.memory_motion[8*self.N:9*self.N]
        #self.Btest = np.zeros(self.N)

        self.XbetaI2 = self.memory_motion[9*self.N:10*self.N]
        self.YbetaI2 = self.memory_motion[10*self.N:11*self.N]
        self.ZI = np.zeros(self.N)

        #self.bRough = bRough
    
    def append(self, motion, Zshift = None):
        self.X = np.append(self.X, motion.X)
        self.Y = np.append(self.Y, motion.Y)
        #print "test ", self.Z, len(self.Z)
        if Zshift == None:
            Zshft = 0
        else:
            Zshft = Zshift
        self.Z = np.append(self.Z, motion.Z[:] + Zshft)
        self.Xbeta = np.append(self.Xbeta, motion.Xbeta)
        self.Ybeta = np.append(self.Ybeta, motion.Ybeta)
        self.Zbeta = np.append(self.Zbeta, motion.Zbeta)
        self.Bx = np.append(self.Bx, motion.Bx)
        self.By = np.append(self.By, motion.By)
        self.Bz = np.append(self.Bz, motion.Bz)

        self.XbetaI2 = np.append(self.XbetaI2, motion.XbetaI2)
        self.YbetaI2 = np.append(self.YbetaI2, motion.YbetaI2)
        
    def defineV(self):
        V = np.matrix([[self.X[-1]], [self.Xbeta[-1]], [1.], 
                        [self.Y[-1]], [self.Ybeta[-1]], [1.]])
        return V  
        
    def defineInitCond(self, gamma):
        IC = np.zeros(6)
        IC[0] = self.X[-1]
        IC[1] = self.Y[-1] #*1e-3
        IC[2] = self.Z[-1]
        IC[3] = self.Xbeta[-1]
        IC[4] = self.Ybeta[-1] #*1e-3
        IC[5] = gamma
        return IC      


    def Py2C(self):
        mtn = np.hstack((self.X, self.Y, self.Z, 
                        self.Xbeta, self.Ybeta, self.Zbeta, 
                        self.Bx, self.By, self.Bz, self.XbetaI2, self.YbetaI2))        
        return mtn

    def gaussIntegration(self):
        if(self.bRough>0):
            print (" gaussBetaSquare ERROR - bRough is not correct")
            return 0;

        size = (len(self.Z) + 1)/3;
        BX2 = BY2 = 0.;

        h2 = (self.Z[-1]-self.Z[0])/2./(size-1);
        w = [0.5555555555555556*h2, 0.8888888888888888*h2, 0.5555555555555556*h2];

        for i in range(size-1): #(i=0; i<size-1;i++)

            for j in range(3): #(j=0;j<3; j++)

                n = i*3 + j +1;
                fx = self.Bx[n];
                fy = self.By[n];
                BX2 += fx*fx*w[j];
                BY2 += fy*fy*w[j];

            #self.XbetaI2[i+1] =  integ1;
            #self.YbetaI2[i+1] =  integ2;
            #self.ZI[i+1] = self.ZI[i] + h2*2.;
        #self.citerpForZdef()
        return BX2, BY2;
        
    
