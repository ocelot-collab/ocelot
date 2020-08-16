__author__ = 'Svitozar Serkez'

from ocelot.rad.fel import FelParameters, calculateFelParameters, beam2fel
# from ocelot.adaptors.genesis import *
# from ocelot.rad.undulator_params import Ephoton2K
from ocelot.cpbd.beam import parray2beam
from ocelot.cpbd.track import update_effective_beta
from ocelot.cpbd.elements import *
import numpy as np
from ocelot.common.globals import *
from copy import deepcopy
import logging

_logger = logging.getLogger(__name__)

def beamlat2fel(beam, lat, smear_m=None, method='mxie', qf=0):
    
    _logger.info('estimating fel from beam and lat')
    
    beam_pk = beam[beam.idx_max()]
    
    E_beam = beam_pk.E
    
    indx_u = np.where([i.__class__ == Undulator for i in lat.sequence])[0]

    und = lat.sequence[indx_u[0]]
    l_period = und.lperiod
    # und.Kx = Ephoton2K(E_photon, und.lperiod, E_beam)
    K_peak = np.max([und.Kx, und.Ky])
    if und.Kx != und.Ky:
        iwityp = 0 # planar undulator
    elif und.Kx == und.Ky:
        iwityp = 1 # helical undulator
    else:
        raise ValueError('unknown undulator: neither planar nor helical, estimation not applicable')
    
    if smear_m is None:
        tcoh = beam2fel(beam_pk, l_period, K_peak, iwityp=iwityp).tcoh()
        smear_m = tcoh * speed_of_light * 4 #smear window

    beam_tmp = deepcopy(beam)
    beam_tmp.smear(smear_m)
    update_effective_beta(beam_tmp, lat)
    beam_tmp.I[beam_tmp.I==0] = np.nan
    fel = beam2fel(beam_tmp, l_period, K_peak, method='mxie', qf=qf)
    
    return fel


def parraylat2fel(parray, lat, smear=1e-7):
    
    _logger.info('estimating fel from parray and lat')
    
    beam = parray2beam(parray, step = 2 * smear)
    
    return beamlat2fel(beam, lat, smear_m = smear)
    
