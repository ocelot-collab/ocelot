# C. Lechner, European XFEL, 20240508
#
# Test case moved here from ocelot/rad/fel.py

from ocelot.rad.fel import FelParameters,calculateFelParameters
import numpy as np
from copy import deepcopy

'''
CHECK with Xie paper parameters

CL, Nov-2022: Probably
M. Xie, "Exact and variational solutions of 3D eigenmodes in high gain FELs", Nucl. Instruments Methods Phys. Res. Sect. A Accel. Spectrometers, Detect. Assoc. Equip., vol. 445, no. 1--3, pp. 59--66, 2000, DOI https://doi.org/10.1016/S0168-9002(00)00114-5
p. 63, left column
'''
def test_fel_MX():
    # inp = GenesisInput()
    inp = FelParameters()
    inp.curpeak = 3400
    inp.xlamd = 0.03
    inp.iwityp = 0 # 0=planar undulator
    inp.gamma0 = 28000
    inp.delgam = inp.gamma0 * 2e-4
    inp.betax = 18
    inp.betay = 18
    inp.emitx=1.5e-6
    inp.emity=1.5e-6
    inp.xlamd=0.03
    inp.aw0 = 3.7/np.sqrt(2)
    #
    p = calculateFelParameters(inp)
    #print(p.xie_lscale,'new')
    #p.lg1
    #p.rho1
    print('Comparison with Ming Xie paper')
    print('etad:     expected 0.0367, got {}'.format(p.xie_etad))
    print('etae:     expected 0.739,  got {}'.format(p.xie_etae))
    print('etagamma: expected 0.248,  got {}'.format(p.xie_etagamma))

    # Score the test result by comparing with values given in Xie paper
    # -> if needed, the evaluation of the test can be modified
    reldev = np.array([p.xie_etad/0.0367-1, p.xie_etae/0.739-1, p.xie_etagamma/0.248-1])
    relscore = np.sum(reldev*reldev)
    print('score={:e}'.format(relscore))
    maxrelscore = 3e-6 # I get relscore=1.74e-6
    ok = (relscore<maxrelscore)
    assert(ok)

    '''
    Check the values of some key parameters to protect against
    unwanted side-effects changing these results.
    Note that we compare with values computed with OCELOT source
    code version git commit id=eb3fbdf (date: 2024-05-08).
    '''
    def check_small_reldev(a,b, maxreldev=1e-6):
        reldev=a/b-1
        assert(abs(reldev)<maxreldev) # if condition=False -> test fails
    
    check_small_reldev(p.lg3,  5.73242744516371)
    check_small_reldev(p.rho3, 0.0002404430325092614)
