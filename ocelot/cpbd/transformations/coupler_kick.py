from copy import copy
import logging

import numpy as np

from ocelot.cpbd.transformations.transfer_map import TransferMap

_logger = logging.getLogger(__name__)


class CouplerKickTM(TransferMap):
    def __init__(self, v=0, freq=0., phi=0., vx=0., vy=0.):
        TransferMap.__init__(self)
        self.v = v
        self.freq = freq
        self.phi = phi
        self.vx = vx
        self.vy = vy
        self.B_z = lambda z, energy: self.kick_b(self.v, self.phi, energy)

    def kick_b(self, v, phi, energy):
        phi = phi * np.pi / 180.
        dxp = (self.vx * v * np.exp(1j * phi)).real / energy
        dyp = (self.vy * v * np.exp(1j * phi)).real / energy
        b = np.array([[0.], [dxp], [0.], [dyp], [0.], [0.]])
        return b

    def kick(self, X, v, phi, energy):
        _logger.debug('invoking coupler kick zero order')
        b = self.kick_b(v, phi, energy)
        X1 = np.dot(self.R(energy), X)
        X1 = np.add(X1, b)
        X[:] = X1[:]
        return X

    def map_function(self, delta_length=None, length=None):
        return lambda X, energy: self.kick(X, self.v, self.phi, energy)

    @classmethod
    def create_from_element(cls, element, params=None):
        return cls(v=element.v, freq=element.freq, phi=element.phi, vx=element.vx, vy=element.vy)
