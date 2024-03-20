#!/usr/bin/env python3

# Christoph Lechner, European XFEL, 2024-02-29
#
# This demonstrates the newly introduced functionality to
# manipulate "GENESIS 1.3" v4 lattices during the generation phase.
#
# In this demo, the same lattice is written twice:
# . one time with explicit 'aw' values (default behavior)
# . one time with references to a sequence (this a new feature added
#   by CL, available in stable "GENESIS 1.3", v4 releases in early 2024)

# import from Ocelot main modules and functions
import ocelot
from ocelot.adaptors.genesis4 import *
from ocelot.gui.beam_plot import *
from ocelot.utils.xfel_utils import create_fel_beamline

from copy import deepcopy



# Demo for G4 adaptor lattice callback:
# Generates undulator beamline lattice with aw values from sequence.
# This is based on modified code from OCELOT genesis4 adaptor
def my_cb_latline(lat, element, eleidx, element_prefix):
    # print(f'callback idx={eleidx}, data obj={str(element)}')

    ### code from OCELOT genesis4 adaptor ###
    if isinstance(element,Undulator):
        element_name = element_prefix + 'UND'
        if element.Kx == 0 or element.Ky == 0:
            is_helical = 'false'
        elif element.Kx == element.Ky:
            is_helical = 'true'
        else:
            _logger.warning(
                'undulator element with different non-zero Kx and Ky; not implemented; setting to helical')
            is_helical = 'true'
        aw = np.sqrt((element.Kx ** 2 + element.Ky ** 2) / 2)
        # TODO: implement ax,ay,kx,ky,gradx,grady
        if hasattr(element,'awref'):
            # aw value is reference to sequence defined in G4 main input file
            # -> aw parameter of instance of class Undulator has no effect!
            s = '{:}: UNDULATOR = {{lambdau = {:}, nwig = {:}, aw = @{:}, helical = {:}}};'.format(element_name,
                                                                                                  element.lperiod,
                                                                                                  element.nperiods,
                                                                                                  element.awref, is_helical)
        else:
            s = '{:}: UNDULATOR = {{lambdau = {:}, nwig = {:}, aw = {:.6f}, helical = {:}}};'.format(element_name,
                                                                                                     element.lperiod,
                                                                                                     element.nperiods,
                                                                                                     aw, is_helical)
        return (element_name,s)

    if isinstance(element,Marker) and hasattr(element,'verbatim'):
        # Another possible application might be a Marker that carries some additional information that is written verbatim to the lattice file.
        # Use case: lattice elements not supported yet; comments
        # Note that the callback function currently cannot not prevent the generation of the corresponding lattice elements (the assumption is that every element in the lattice sequence has its corresponding lattice element in the file), but as workaround one can generate one string containing multiple lines (with \n between)
        pass # no implemented yet

    # It is not one of the special cases -> signal that no special handling required
    return (None,None)
    pass

#########################
### MAIN part of demo ###
#########################
# matching procedure verified (same k1 values obtained) for:
# . standard beamline length (without override_und_N parameter)
# . override_und_N=2

# If you want the standard lattice
# sase_lat_pkg = create_fel_beamline('sase1')

'''
Shorter lattice with SASE1-type FODO (uses functionality added by CL to OCELOT in Feb 2024)
New feature added in March-2024: Request complete final half-cell
-> this allows to match two-undulator-cell lattice
'''
sase_lat_pkg = create_fel_beamline('sase1', override_und_N=2, final_halfcell_complete=True)

# just in case you want to see the lattice before matching
#sase_lat, _, _ = sase_lat_pkg
#print(str(sase_lat))

# set up electron optics, undulator aw value
betaavg = 20 # m
Ephot = 9000 # eV
beam = generate_beam(
    E=14, dE=0.003,
    I=5000, shape='flattop', l_beam=1e-6, l_window=1.5e-6,
    emit_n=0.4e-6, beta=betaavg, nslice=2)
prepare_el_optics(beam, sase_lat_pkg, E_photon=Ephot, beta_av=betaavg)

sase_lat, _, _ = sase_lat_pkg
sase_lat2 = deepcopy(sase_lat)
#print(str(sase_lat))

# Add custom property to the instances of class Undulator.
# In the callback function they trigger special functionality as the
# lattice is written out.
for ele in sase_lat.sequence:
    if isinstance(ele,Undulator):
        ele.awref='aw.seq' # add new property to class instance
for ele in sase_lat2.sequence:
    if isinstance(ele,Undulator):
        ele.awref='aw.seq2' # add new property to class instance

### Write lattices ###
# standard lattice
write_gen4_lat({'LX':sase_lat,'LY':sase_lat2}, filepath='./lat_noref.lat')
# lattice with callback function in effect, with aw values replaced by references to sequence
write_gen4_lat({'LX':sase_lat,'LY':sase_lat2}, filepath='./lat_ref.lat', cb_latline=my_cb_latline)


