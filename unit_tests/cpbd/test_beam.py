import pytest
import numpy as np

from ocelot.cpbd.beam import (cov_matrix_from_twiss,
                              cov_matrix_to_parray,
                              optics_from_moments,
                              ParticleArray,
                              Twiss,
                              twiss_iterable_to_df
                              )

# ex=1.858178000891106e-11,
# ey=1.858178000891106e-11,
# sigma_tau=50e-15,
# sigma_p=0.0001,
# alpha_x=-0.5542165316,
# beta_x=9.701136465,
# alpha_y=2.310858304,
# beta_y=46.95602673,
# dx=0.01497010119,
# dpx=0.00697330389,
# dy=-0.02816641057,
# dpy=0.002391521083


def test_cov_matrix_to_parray():
    # Using entrance to TD20 as test case as it has dispersion in both planes.
    mean = [0, 0, 0, 0, 0, 0]
    cov = cov_matrix_from_twiss(ex=1.858178000891106e-11,
                                ey=1.858178000891106e-11,
                                sigma_tau=50e-15,
                                sigma_p=0.0001,
                                alpha_x=-0.5542165316,
                                beta_x=9.701136465,
                                alpha_y=2.310858304,
                                beta_y=46.95602673,
                                dx=0.01497010119,
                                dpx=0.00697330389,
                                dy=-0.02816641057,
                                dpy=0.002391521083
                                )
    energy = 16.5
    charge = 0.2e-9
    # Hand verified for big N that this results in beam array with correct
    # moments.  Here just generate one
    np.random.seed(0)
    parray = cov_matrix_to_parray(mean, cov, 16.5, 0.2e-9, 1)

    twiss = optics_from_moments(mean, cov)

    assert np.isclose(parray.rparticles, np.array([[ 1.06185830e-05],
                                                   [-3.56949798e-06],
                                                   [-6.30814792e-06],
                                                   [ 1.30761605e-06],
                                                   [ 4.88638940e-14],
                                                   [-1.76951338e-04]])).all()
    assert parray.E == energy
    assert sum(parray.q_array) == charge

def test_optics_from_moments():
    pass

@pytest.fixture
def a_twiss_dictionary():
    optics_dict = {'emit_x': 0.34763928534510036,
                   'emit_y': 0.6790349802898668,
                   'emit_xn': 0.10046549579167263,
                   'emit_yn': 0.22230695589625382,
                   'eigemit_1': 0.4896032635908605,
                   'eigemit_2': 0.4621735761528035,
                   'beta_x': 0.1288947746254585,
                   'beta_y': 0.09036354359497034,
                   'alpha_x': 0.9597776593035088,
                   'alpha_y': 0.8390564421996655,
                   'gamma_x': 0.14540121338731515,
                   'gamma_y': 0.7244156996636956,
                   'Dx': 0.010293384598673794,
                   'Dy': 0.8162183638949975,
                   'Dxp': 0.4672056535332062,
                   'Dyp': 0.2994689154547381,
                   'mux': 0.9993182745796386,
                   'muy': 0.6446622051396393,
                   'E': 0.5095571020715003,
                   's': 0.4799985227868193,
                   'q': 0.6522847419197526,
                   'x': 0.2852771630754366,
                   'y': 0.5087154088422419,
                   'p': 0.6810071075764375,
                   'tau': 0.067116684235901,
                   'xp': 0.140082529580558,
                   'yp': 0.2847354314054312,
                   'xx': 0.48184506797105153,
                   'xpx': 0.6798247854790131,
                   'pxpx': 0.26775563153047577,
                   'yy': 0.2375805309768111,
                   'ypy': 0.8878538454820061,
                   'pypy': 0.33197187208200896,
                   'tautau': 0.7690410387155361,
                   'xy': 0.5916561916755728,
                   'pxpy': 0.790476675553196,
                   'xpy': 0.2717175021299518,
                   'ypx': 0.6684978344712186,
                   "id": "Hello my friends!"
                   }
    return optics_dict

@pytest.fixture
def a_twiss_instance(a_twiss_dictionary):
    result = Twiss()
    for key, value in a_twiss_dictionary.items():
        if not hasattr(result, key):
            raise AttributeError(f"Twiss instance has no {key} attribute.")
        setattr(result, key, value)
    return result

def test_Twiss_to_series(a_twiss_instance, a_twiss_dictionary):
    # Random numbers between 0 and 1 just for fun.
    twiss = a_twiss_instance
    assert dict(twiss.to_series()) == a_twiss_dictionary

def test_twiss_iterable_to_df(a_twiss_instance):
    twisses = [a_twiss_instance, a_twiss_instance]

    twiss_df = twiss_iterable_to_df(twisses)
    assert dict(twiss_df.iloc[0]) == dict(twisses[0].to_series())
    assert dict(twiss_df.iloc[1]) == dict(twisses[1].to_series())


def test_particle_array_total_charge():
    parray = ParticleArray(2)
    charge0 = 1e-14
    charge1 = 5e-14
    parray.q_array[0] = charge0
    parray.q_array[1] = charge1

    assert parray.total_charge == charge0 + charge1
