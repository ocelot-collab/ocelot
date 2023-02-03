import logging

import numpy as np

from ocelot.cpbd.transformations.transfer_map import TransferMap

_logger = logging.getLogger(__name__)


class PulseTM(TransferMap):
    def __init__(self, kn):
        TransferMap.__init__(self)

    def mul_parray(self, rparticles, energy=0.):
        n = len(rparticles)
        # if 'pulse' in self.__dict__:
        _logger.debug('TD transfer map')
        if n > 6: _logger.debug(
            'warning: time-dependent transfer maps not implemented for an array. Using 1st particle value')
        if n > 6: _logger.debug('warning: time-dependent transfer maps not implemented for steps inside element')
        tau = rparticles[4]
        dxp = self.pulse.kick_x(tau)
        dyp = self.pulse.kick_y(tau)
        _logger.debug('kick ' + str(dxp) + ' ' + str(dyp))
        b = np.array([0.0, dxp, 0.0, dyp, 0., 0.])
        # a = np.add(np.transpose(dot(self.R(energy), np.transpose(particles.reshape(int(n / 6), 6)))), b).reshape(n)
        a = np.add(np.dot(self.R(energy), rparticles), b)
        rparticles[:] = a[:]
        _logger.debug('return trajectory, array ' + str(len(rparticles)))
        return rparticles
