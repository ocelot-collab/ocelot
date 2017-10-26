__author__ = 'Svitozar Serkez'

# from ocelot import *
from ocelot.rad.fel import FelParameters, calculateFelParameters, beam2fel
# from ocelot.adaptors.genesis import *
# from ocelot.rad.undulator_params import Ephoton2K
from ocelot.cpbd.beam import array2beam
from ocelot.cpbd.track import update_effective_beta
from ocelot.cpbd.elements import *
import numpy as np
from ocelot.common.globals import *


def beamlat2fel(beam, lat):
    
    beam_pk = beam[beam.idx_max()]
    
    E_beam = beam_pk.E
    
    indx_u = np.where([i.__class__ == Undulator for i in lat.sequence])[0]
    
    und = lat.sequence[indx_u[0]]
    l_period = und.lperiod
    # und.Kx = Ephoton2K(E_photon, und.lperiod, E_beam)
    K_peak = und.Kx
    
    tcoh = beam2fel(beam_pk, l_period, K_peak).tcoh()
    
    sw = tcoh * speed_of_light * 2 #smear window
    beam.smear(sw)
    
    update_effective_beta(beam, lat)
    fel = beam2fel(beam, l_period, K_peak)
    
    return fel


def parraylat2fel(parray, lat, step = 1e-7):

    beam = array2beam(parray, step = step)
    
    return beamlat2fel(beam, lat)
    
