__author__ = "Nikita Kuklev"

import time
from copy import copy

from unit_tests.params import *
from edrift_conf import *


def test_lattice_transfer_map(lattice):
    """R matrix with regular drift"""
    dr_matrix = lattice.sequence[1].transfer_map.R(0)
    d_expected = d_matrix
    result = check_matrix(dr_matrix, d_expected, TOL, assert_info=' r_matrix - ')
    assert check_result(result)

    r_matrix = lattice_transfer_map(lattice, 0.0)
    r_expected = m_matrix @ d_matrix @ d_matrix @ m_matrix
    result = check_matrix(r_matrix, r_expected, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_lattice_transfer_map_exact(lattice_exact):
    """R matrix with exact drift"""
    lattice = lattice_exact

    d_matrix = lattice.sequence[1].transfer_map.R(0)
    d_expected = d_matrix
    result = check_matrix(d_matrix, d_expected, TOL, assert_info=' r_matrix - ')
    assert check_result(result)

    r_matrix = lattice_transfer_map(lattice, 0.0)
    r_expected = m_matrix @ d_matrix @ d_matrix @ m_matrix
    result = check_matrix(r_matrix, r_expected, TOL, assert_info=' r_matrix - ')
    assert check_result(result)


def test_lattice_twiss_exact(lattice_exact):
    """Twiss with exact drift"""
    lattice = lattice_exact
    tws = twiss(lattice)

    result1 = check_value(tws[1].beta_x, beta_edge)
    result2 = check_value(tws[2].beta_x, beta_star)
    result3 = check_value(tws[1].alpha_x, -alpha_edge)
    result4 = check_value(tws[2].alpha_x, 0.0, tolerance_type='absolute')
    assert check_result([result1, result2, result3, result4])

    result1 = check_value(tws[1].beta_y, beta_edge)
    result2 = check_value(tws[2].beta_y, beta_star)
    result3 = check_value(tws[1].alpha_y, -alpha_edge)
    result4 = check_value(tws[2].alpha_y, 0.0, tolerance_type='absolute')
    assert check_result([result1, result2, result3, result4])


def test_lattice_setup():
    """ Check proper override of element default method """
    d = EDrift(l=l / 2, eid='D')  # defaults to method=1
    m = Matrix(eid='M')
    m.r = m_matrix.copy()
    mt = MethodTM(params={'edrift_method': 3})
    lat = MagneticLattice((m, d, d, m), method=mt)
    r1 = check_value(lat.sequence[1].method, 1) # Still 1
    r2 = check_value(lat.sequence[1].transfer_map.method, 3) # But map has 3
    assert check_result([r1, r2])


@pytest.mark.parametrize('config', [0, 1, 2, 3])
def test_lattice_track(drift_exact, config):
    """
    Its hard to test tracking directly, so compare to MAD-X, which uses exact drifts and similar coordinates
    Note that floating point epsilon is ~2.22e-16, need optimization to reduce error accumulation

    0 - on axis, all 0
    1 - on axis, energy offset
    2 - with transverse offsets
    3 - with all offsets

    For reference, MADX script is:
        beam,particle=electron,energy=0.1,radiate=false;
        d: DRIFT, L=2.0;
        iota: LINE=(d);
        use,sequence=iota;
        track,dump,damp=false,quantum=false,file="trackmadx";
            start, pt=-0.001;
            start, x=-0.3, y=0.7, px=0.1, py=-0.25;
            start, x=0.10, y=-0.2, px=-0.15, py=-0.07, t=-0.01, pt=0.03;
            run,turns=2,ffile=1;
        endtrack;
        stop;
    """
    lattice = drift_exact
    if config == 0:
        particle = Particle()
        ans0 = ans1 = ans2 = (0, 0, 0, 0, 0, 0, 0)
    elif config == 1:
        particle = Particle(p=-0.001, E=0.1)
        ans0 = (0, 0, 0, 0, 0, -0.001, 0.1)
        ans1 = (0, 0, 0, 0, -5.230379184e-08, -0.001, 0.1)
        ans2 = (0, 0, 0, 0, -1.046075837e-07, -0.001, 0.1)
    elif config == 2:
        particle = Particle(x=-0.3, px=0.1, y=0.7, py=-0.25, p=0.0, E=0.1)
        ans0 = (-0.3, 0.1, 0.7, -0.25, 0, -2.220446049e-16, 0.1)
        ans1 = (-0.09233034734, 0.1, 0.1808258683, -0.25, -0.07669752797, -2.220446049e-16, 0.1)
        ans2 = (0.1153393053, 0.1, -0.3383482633, -0.25, -0.1533950559, -2.220446049e-16, 0.1)
    elif config == 3:
        particle = Particle(x=0.1, px=-0.15, y=-0.2, py=-0.07, tau=0.01, p=0.03, E=0.1)
        ans0 = (0.1, -0.15, -0.2, -0.07, -0.01, 0.03, 0.1)
        ans1 = (-0.195097717, -0.15, -0.3377122679, -0.07, -0.0363372301, 0.03, 0.1)
        ans2 = (-0.490195434, -0.15, -0.4754245359, -0.07, -0.06267446021, 0.03, 0.1)
    else:
        raise Exception('?')

    p = particle
    navi = Navigator(lattice)
    dz = 0.01
    P = [copy(p)]
    for iii in range(int(lattice.totalLen / dz)):
        tracking_step(lattice, p, dz=dz, navi=navi)
        P.append(copy(p))
    navi2 = Navigator(lattice)
    tracking_step(lattice, p, dz=lattice.totalLen, navi=navi2)
    P.append(copy(p))

    p0, p1, p2 = P[0], P[-2], P[-1]
    p0.tau *= -1
    p1.tau *= -1
    p2.tau *= -1

    vars = ['x', 'px', 'y', 'py', 'tau', 'p', 'E']

    for p, ans in zip((p0, p1, p2), (ans0, ans1, ans2)):
        results = []
        for i, v in enumerate(vars):
            results.append(check_value(getattr(p, v), ans[i], tolerance=1e-10, tolerance_type='absolute'))
        assert check_result(results)


def setup_module(module):
    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### EDRIFT START ###\n\n')
    f.close()


def teardown_module(module):
    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write('### EDRIFT END ###\n\n\n')
    f.close()


def setup_function(function):
    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write(function.__name__)
    f.close()

    pytest.t_start = time.time()


def teardown_function(function):
    f = open(pytest.TEST_RESULTS_FILE, 'a')
    f.write(' execution time is ' + '{:.3f}'.format(time.time() - pytest.t_start) + ' sec\n\n')
    f.close()
