# C. Lechner, European XFEL, 2024-03-24

import ocelot
from ocelot.adaptors.genesis4 import *
from ocelot.gui.beam_plot import *
from ocelot.utils.xfel_utils import create_fel_beamline

def my_cb_latline(lat, element, eleidx, element_prefix):
	# put marker for this test
    element.cb_invoked=True
    
    # standard behavior of G4 adaptor
    return (None,None)


# Remark: there are 30+ PendingDeprecationWarning messages, related to genesis.py and the use of np.matrix
def test_g4adaptor_writelat_cb(tmpdir):
    '''
    Code is from demo_lat_callback.py
    '''


    '''
    Shorter lattice with SASE1-type FODO (uses functionality added by CL to OCELOT in Feb 2024)
    New feature added in March-2024: Request complete final half-cell
    -> this allows to match two-undulator-cell lattice
    '''
    sase_lat_pkg = create_fel_beamline('sase1', override_und_N=2, final_halfcell_complete=True)

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

    # working directory for the test
    # https://docs.pytest.org/en/6.2.x/tmpdir.html
    wd = tmpdir
    import os.path

    ### Write lattices ###
    # standard lattice
    write_gen4_lat({'LX':sase_lat,'LY':sase_lat2}, filepath=os.path.abspath(wd.join('./lat_no_cb.lat')))
    # lattice with callback function in effect, with aw values replaced by references to sequence
    write_gen4_lat({'LX':sase_lat,'LY':sase_lat2}, filepath=os.path.abspath(wd.join('./lat_with_cb.lat')), cb_latline=my_cb_latline)
    
    # Arriving here, we know that function 'write_gen4_lat' accepted the 'cb_latline' argument

    # check for presence of marker put in callback function
    assert sase_lat.sequence[0].cb_invoked==True
