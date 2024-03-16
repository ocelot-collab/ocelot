# C. Lechner, European XFEL, 20240316

from ocelot import *
from ocelot.utils.xfel_utils import create_fel_beamline

def helper_count_undulators(pkg):
    test_lat, _, _ = pkg
    ucntr=0
    for ele in test_lat.sequence:
        if isinstance(ele,Undulator):
            ucntr+=1
    return(ucntr)

def test_create_fel_beamline_Nund_default():
    # remark CL 20240316: undulator counts also pass test against commit id 16425bd (date: 2023-11-03)
    test_lat_pkg = create_fel_beamline('sase1')
    ucntr = helper_count_undulators(test_lat_pkg)
    assert 37==ucntr
    #
    test_lat_pkg = create_fel_beamline('sase2')
    ucntr = helper_count_undulators(test_lat_pkg)
    assert 37==ucntr
    #
    test_lat_pkg = create_fel_beamline('sase3')
    ucntr = helper_count_undulators(test_lat_pkg)
    assert 21==ucntr

# Test that we always get the requested number of undulators #
def test_create_fel_beamline_Nund():
    for Nundreq in range(2,38):
        test_lat_pkg = create_fel_beamline('sase1', override_und_N=Nundreq)
        # check if number of undulators is as expected
        ucntr = helper_count_undulators(test_lat_pkg)
        assert Nundreq==ucntr
