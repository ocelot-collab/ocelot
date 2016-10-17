    # -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 21:39:18 2015

@author: Mitch
"""

from ocelot.optimizer.GP.GP import SPGP_train, SPGP_predict

class SPGP(object):
    def __init__(self):
        self.X = []
        self.Y = []
        self.xb = []
        self.m = 0
        self.dim = 0
        self.hyps = ()
        
    def fit(self,X,Y,m):
        self.X = X
        self.Y = Y
        self.m = m
        self.dim = X.shape[1]
        (self.xb, self.hyps) = SPGP_train(X,Y,m)
        
    # too slow in practice
    def update(self, x_new, y_new):
        self.X.loc[len(self.X.index)] = x_new[0]
        self.Y[len(self.Y)] = y_new
        (self.xb, self.hyps) = SPGP_train(self.X,self.Y,min(self.m,self.X.shape[0]))
        
    
    def predict(self,X):
        return SPGP_predict(self.X, self.Y, self.xb, X, self.hyps)
        
    
    