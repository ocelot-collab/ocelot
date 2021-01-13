from copy import copy
from math import factorial

import numpy as np

from ocelot.cpbd.transformations.transfer_map import TransferMap


class MultipoleTM(TransferMap):
    def __init__(self, kn):
        TransferMap.__init__(self)
        self.kn = kn

    def kick(self, X, kn):
        p = -kn[0] * X[5] + 0j
        for n in range(1, len(kn)):
            p += kn[n] * (X[0] + 1j * X[2]) ** n / factorial(n)
        X[1] = X[1] - np.real(p)
        X[3] = X[3] + np.imag(p)
        X[4] = X[4] - kn[0] * X[0]
        return X

    def map_function(self, delta_length=None, length=None):
        return lambda X, energy: self.kick(X, self.kn)

    @classmethod
    def create_from_element(cls, element, params=None):
        return cls(kn=element.kn)
