# -*- coding: utf-8 -*-
"""
Contains simple interfaces for the Bayes optimization class.

Each interface must have the getState and setX methods as used below.
"""

import numpy as np

# Basically a dummy interface that just holds x-values
class TestInterface(object):
    def __init__(self, x_init, y_init=0):
        self.x = np.array(x_init,ndmin=2)
        self.y = y_init
        
    def getState(self):
        return self.x, self.y
            
    def setX(self, x_new):
        self.x = x_new
        
# an interface that evaluates a function to give y-values
class fint(object):
    def __init__(self, x_init, noise=0):
        self.x = np.array(x_init,ndmin=2)
        self.y = -1
        self.noise = noise
        
    def getState(self):
        return self.x, self.f(self.x)
        
    def f(self, x):
        res = 20 - np.reshape(np.sum((x - 1)**2,axis=1) * np.sum((x + 1)**2,axis=1),x.shape) + x
        res[np.abs(x) > 3.0] = 0
        return res
        
    def setX(self, x_new):
        self.x = x_new
        
# uses a GP model's predictions to give y-values
class GPint(object):
    def __init__(self, x_init, model):
        self.model = model
        self.x = x_init
        self.y = model.predict(np.array(x_init,ndmin=2))
        
    def getState(self):
        return self.x, self.model.predict(np.array(self.x,ndmin=2))[0]
        
    def setX(self, x_new):
        self.x = x_new
        
    
