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

def beamlat2fel(beam, lat, smear_m=1e-6):
    
    beam_pk = beam[beam.idx_max()]
    
    E_beam = beam_pk.E
    
    indx_u = np.where([i.__class__ == Undulator for i in lat.sequence])[0]
    
    und = lat.sequence[indx_u[0]]
    l_period = und.lperiod
    # und.Kx = Ephoton2K(E_photon, und.lperiod, E_beam)
    K_peak = und.Kx
    
    if smear_m is None:
        tcoh = beam2fel(beam_pk, l_period, K_peak).tcoh()
        smear_m = tcoh * speed_of_light * 2 #smear window

    beam_tmp = deepcopy(beam)
    beam_tmp.smear(smear_m)
    update_effective_beta(beam_tmp, lat)
    fel = beam2fel(beam_tmp, l_period, K_peak)
    
    return fel


def parraylat2fel(parray, lat, step=1e-7):

    beam = parray2beam(parray, step = 2 * step)
    
    return beamlat2fel(beam, lat, smear_m = step)
    
