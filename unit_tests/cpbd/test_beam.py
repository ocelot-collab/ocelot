import numpy as np

from ocelot.cpbd.beam import (cov_matrix_from_twiss,
                              cov_matrix_to_parray,
                              optics_from_moments
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
