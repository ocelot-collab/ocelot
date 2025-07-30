import multiprocessing
from copy import deepcopy

import numpy as np
from scipy.optimize import *

from ocelot.cpbd.beam_params import radiation_integrals
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.beam import Particle
from ocelot.cpbd.elements import *
from ocelot.cpbd.beam import get_envelope
from ocelot.cpbd.track import track
from ocelot.cpbd.optics import lattice_transfer_map, twiss, Twiss
from ocelot.cpbd.elements.optic_element import OpticElement
from ocelot.cpbd.tm_utils import SecondOrderMult


def weights_default(val):
    if val == 'periodic':
        return 10000001.0
    if val == 'total_len':
        return 10000001.0
    if val in ['Dx', 'Dy']:
        return 10000002.0
    if val in ['Dxp', 'Dyp']:
        return 10000003.0
    if val == 'tau':
        return 10000004.0
    if val == 'i5':
        return 1.e14
    if val == 'negative_length':
        return 1.5e6
    if val in ['alpha_x', 'alpha_y']:
        return 100007.0
    if val in ['mux', 'muy']:
        return 10000006.0
    if val in ['beta_x', 'beta_y']:
        return 100007.0
    return 0.0001


def match(lat, constr, vars, tw, verbose=True, max_iter=1000, method='simplex', weights=weights_default,
          vary_bend_angle=False, min_i5=False, tol=1e-5):
    """
    Function to match twiss parameters. To find periodic solution for a lattice use MagneticLattice.periodic_twiss(tws)

    :param lat: MagneticLattice
    :param constr: dictionary, constrains. Example:
            'periodic':True - means the "match" function tries to find periodic solution at the ends of lattice:
                constr = {elem1:{'beta_x':15, 'beta_y':2}, 'periodic':True}

            "hard" constrains on the end of elements (elem1, elem2):
                constr = {elem1:{'alpha_x':5, 'beta_y':5}, elem2:{'Dx':0 'Dyp':0, 'alpha_x':5, 'beta_y':5}}

            or mixture of "soft" and "hard" constrains:
                constr = {elem1:{'alpha_x':[">", 5], 'beta_y':5}, elem2:{'Dx':0 'Dyp':0, 'alpha_x':5, 'beta_y':[">", 5]}}

                in case one wants to constrain absolute value of variable, the constrains can be:
                constr = {elem1:{'alpha_x':["a>", 5], "alpha_y": ["a<", 1]}}
                        - which means np.abs(alpha_x) > 5 and np.abs(alpha_y) < 1

            in case one needs global control on beta function, the constrains can be written following way.
                constr = {elem1:{'alpha_x':5, 'beta_y':5}, 'global': {'beta_x': ['>', 10]}}

            Experimental constrain (CAN BE DISABLED or CHANGED AT ANY MOMENT)
                constr = {"delta": {ELEM1: ["muy", 0],  ELEM2: ["muy", 0], "val":  3*np.pi/2, "weight": 100007}}
                        - try to satisfy: tws.muy at ELEM2 - tws.muy at ELEM1 == 'val'
                        - difference between ELEM1 and ELEM2 of twiss parameter "muy" (can be any) == "val"
                        - ELEM1: ["muy", 0] - 0 here means nothing it is just place to store value of muy
                        - pay attantion to order of elements. since val is sensitive to sign. SO element

    :param vars: list of elements e.g. vars = [QF, QD] or it can be initial twiss parameters:
                vars = [[tws0, 'beta_x'], [tws0, 'beta_y'], [tws0, 'alpha_x'], [tws0, 'alpha_y']].
                A tuple of quadrupoles can be passed as a variable to constrain their strengths
                to the same value, e.g., vars = [(QF, QD)].
                A dictionary with quadrupoles as keys and their relative strengths as values
                can be used for more flexible constraints, e.g., vars = [{QF: 1.0, QD: -1.0}],
                which constrains QD and QF to have equal strengths with opposite signs.
    :param tw: initial Twiss
    :param verbose: allow print output of minimization procedure
    :param max_iter:
    :param method: string, available 'simplex', 'cg', 'bfgs'
    :param weights: function returns weights, for example
                    def weights_default(val):
                        if val == 'periodic': return 10000001.0
                        if val == 'total_len': return 10000001.0
                        if val in ['Dx', 'Dy', 'Dxp', 'Dyp']: return 10000002.0
                        if val in ['alpha_x', 'alpha_y']: return 100007.0
                        if val in ['mux', 'muy']: return 10000006.0
                        if val in ['beta_x', 'beta_y']: return 100007.0
                        return 0.0001
    :param vary_bend_angle: False, allow to vary "angle" of the dipoles instead of the focusing strength "k1"
    :param min_i5: minimization of the radiation integral I5. Can be useful for storage rings.
    :param tol: tolerance default 1e-5
    :return: result
    """

    # tw = deepcopy(tw0)

    def errf(x):

        tw_loc = deepcopy(tw)
        tw0 = deepcopy(tw)


        # parameter to be varied is determined by variable class

        for i in range(len(vars)):
            if isinstance(vars[i], Drift):
                if x[i] < 0:
                    # print('negative length in match')
                    return weights('negative_length')

                vars[i].l = x[i]
            if isinstance(vars[i], Quadrupole):
                vars[i].k1 = x[i]
            if isinstance(vars[i], Solenoid):
                vars[i].k = x[i]
            if isinstance(vars[i], (RBend, SBend, Bend)):
                if vary_bend_angle:
                    vars[i].angle = x[i]
                else:
                    vars[i].k1 = x[i]
            if isinstance(vars[i], list):
                if isinstance(vars[i][0], Twiss) and isinstance(vars[i][1], str):
                    k = vars[i][1]
                    setattr(tw_loc, k, x[i])
            if isinstance(vars[i], tuple):
            # all quads strength in tuple varied simultaneously
                for v in vars[i]:
                    v.k1 = x[i]
            if isinstance(vars[i], dict):
            # all quads strength in dict keys varied simultaneously
            # with coupling parameters given as values.
                for q in vars[i].keys():
                    q.k1 = vars[i][q] * x[i]

        err = 0.0
        if "periodic" in constr.keys():
            if constr["periodic"]:
                tw_loc = lat.periodic_twiss(tw_loc)
                tw0 = deepcopy(tw_loc)
                if tw_loc is None:
                    print("########")
                    return weights('periodic')

        # save reference points where equality is asked

        ref_hsh = {}  # penalties on two-point inequalities
        for e in constr.keys():
            if e == 'periodic':
                continue
            if e == 'total_len':
                continue
            for k in constr[e].keys():
                if isinstance(constr[e][k], list):
                    if constr[e][k][0] == '->':
                        # print 'creating reference to', constr[e][k][1].id
                        ref_hsh[constr[e][k][1]] = {k: 0.0}
        # evaluating global and point penalties

        tw_loc.s = 0

        for e in lat.sequence:
            for tm in e.first_order_tms:
                tw_loc = tm * tw_loc  # apply transfer map

                # --- Global constraints ---
                if 'global' in constr:
                    for c, rule in constr['global'].items():
                        if isinstance(rule, list):
                            op, v1 = rule[0], rule[1]
                            val = getattr(tw_loc, c)

                            if op == '<' and val > v1:
                                err += weights(c) * (val - v1) ** 2
                            elif op == '>' and val < v1:
                                err += weights(c) * (val - v1) ** 2

                # --- Delta constraint update ---
                if 'delta' in constr and e in constr['delta']:
                    tw_k = constr['delta'][e][0]
                    constr['delta'][e][1] = getattr(tw_loc, tw_k)

                # --- Update reference hash if needed ---
                if e in ref_hsh:
                    ref_hsh[e] = deepcopy(tw_loc)

                # --- Local constraints ---
                if e in constr:
                    for k, rule in constr[e].items():
                        val = getattr(tw_loc, k)

                        if isinstance(rule, list):
                            op = rule[0]
                            v1 = rule[1]

                            if op == '<' and val > v1:
                                err += weights(k) * (val - v1) ** 2
                            elif op == '>' and val < v1:
                                err += weights(k) * (val - v1) ** 2
                            elif op == 'a<' and abs(val) > v1:
                                err += weights(k) * (abs(val) - v1) ** 2
                            elif op == 'a>' and abs(val) < v1:
                                err += weights(k) * (abs(val) - v1) ** 2
                            elif op == '->':
                                try:
                                    dv1 = float(rule[2]) if len(rule) > 2 else 0.0
                                    ref_val = getattr(ref_hsh[v1], k)
                                    err += (val - (ref_val + dv1)) ** 2
                                    if val < v1:
                                        err += (val - v1) ** 2
                                except Exception as ex:
                                    print(f'Constraint error: rval should precede lval in lattice ({ex})')

                            if val < 0:
                                err += (val - v1) ** 2

                        elif isinstance(rule, str):
                            # handle symbolic constraints if any
                            pass
                        else:
                            # direct comparison
                            ref_val = rule
                            err += weights(k) * (ref_val - val) ** 2
        if "total_len" in constr.keys():
            total_len = constr["total_len"]
            err = err + weights('total_len') * (tw_loc.s - total_len) ** 2

        if 'delta' in constr.keys():
            delta_dict = constr['delta']
            elems = []
            for e in delta_dict.keys():
                if isinstance(e, OpticElement):
                    elems.append(e)
            delta_err = delta_dict["weight"] * (delta_dict[elems[1]][1] - delta_dict[elems[0]][1] - delta_dict["val"])**2
            err = err + delta_err

        if min_i5:
            # evaluating integral parameters

            I1, I2, I3, I4, I5 = radiation_integrals(lat, tw0, nsuperperiod=1)
            err += I5 * weights('i5')

            Je = 2 + I4 / I2
            Jx = 1 - I4 / I2
            Jy = 1

            if Je < 0 or Jx < 0 or Jy < 0:
                err = 100000.0

        # c1, c2 = natural_chromaticity(lat, tw0)
        # err += ( c1**2 + c2**2) * 1.e-6

        if verbose:
            print('iteration error:', err)
        return err


    # list of arguments determined based on the variable class

    x = [0.0] * len(vars)
    for i in range(len(vars)):
        if vars[i].__class__ == list:
            if vars[i][0].__class__ == Twiss and vars[i][1].__class__ == str:
                k = vars[i][1]
                x[i] = getattr(tw, k)
                # if k in ['beta_x', 'beta_y']:
                #     x[i] = 10
                # else:
                #     x[i] = 0.0
        if vars[i].__class__ == tuple:
            x[i] = vars[i][0].k1
        if vars[i].__class__ == dict:
            q = list(vars[i].keys())[0]
            x[i] = q.k1 / vars[i][q]
        if vars[i].__class__ == Quadrupole:
            x[i] = vars[i].k1
        if vars[i].__class__ == Drift:
            x[i] = vars[i].l
        if vars[i].__class__ in [RBend, SBend, Bend]:
            # TODO: need a way to vary 2 aattributes of a class
            if vary_bend_angle:
                x[i] = vars[i].angle
            else:
                x[i] = vars[i].k1

    print("initial value: x = ", x)
    if method == 'simplex':
        res = fmin(errf, x, xtol=tol, maxiter=max_iter, maxfun=max_iter)
    if method == 'cg':
        res = fmin_cg(errf, x, gtol=tol, epsilon=1.e-5, maxiter=max_iter)
    if method == 'bfgs':
        res = fmin_bfgs(errf, x, gtol=tol, epsilon=1.e-5, maxiter=max_iter)


    # if initial twiss was varied set the twiss argument object to resulting value

    for i in range(len(vars)):
        if vars[i].__class__ == list:
            if vars[i][0].__class__ == Twiss and vars[i][1].__class__ == str:
                k = vars[i][1]
                setattr(tw, k, res[i])

    # update MagneticLattice total length in case a Drift length was in list of variables

    return res


def weights_default(val):
    if val == 'periodic':
        return 1
    if val == 'total_len':
        return 1
    if val == 'Dx':
        return 1
    if val == 'Dxp':
        return 1
    if val == 'tau':
        return 1
    if val == 'i5':
        return 1
    if val == 'negative_length':
        return 1
    if val in ['alpha_x', 'alpha_y']:
        return 1000
    if val in ['mux', 'muy']:
        return 1
    if val in ['beta_x', 'beta_y']:
        return 1000
    return 1


def match_beam(lat, constr, vars, p_array, navi, verbose=True, max_iter=1000, method='simplex', weights=weights_default,
               vary_bend_angle=False, min_i5=False, bounds=None, slice=None):
    """
    Function to match twiss parameters

    :param lat: MagneticLattice
    :param constr: dict in format {elem1:{'beta_x':15, 'beta_y':2}, 'periodic':True} try to find periodic solution or
                constr = {elem1:{'alpha_x':5, 'beta_y':5}, elem2:{'Dx':0 'Dyp':0, 'alpha_x':5, 'beta_y':5} and so on.
    :param vars: list of elements e.g. vars = [QF, QD]
    :param p_array: initial ParticleArray
    :param navi: Navigator with added PhysProcess if needed
    :param verbose: allow print output of minimization procedure
    :param max_iter:
    :param method: string, available 'simplex', 'cg', 'bfgs'
    :param weights: function returns weights, for example
                    def weights_default(val):
                        if val == 'periodic': return 10000001.0
                        if val == 'total_len': return 10000001.0
                        if val in ['Dx', 'Dy', 'Dxp', 'Dyp']: return 10000002.0
                        if val in ['alpha_x', 'alpha_y']: return 100007.0
                        if val in ['mux', 'muy']: return 10000006.0
                        if val in ['beta_x', 'beta_y']: return 100007.0
                        return 0.0001
    :param vary_bend_angle: False, allow to vary "angle" of the dipoles instead of the focusing strength "k1"
    :param min_i5: minimization of the radiation integral I5. Can be useful for storage rings.
    :return: result
    """

    # tw = deepcopy(tw0)

    run_number = 1
    def errf(x):
        p_array0 = deepcopy(p_array)
        tws = get_envelope(p_array0, bounds=bounds, slice=slice)
        tw_loc = deepcopy(tws)
        tw0 = deepcopy(tws)
        nonlocal run_number


        # parameter to be varied is determined by variable class

        for i in range(len(vars)):
            if vars[i].__class__ == Drift:
                if x[i] < 0:
                    # print('negative length in match')
                    return weights('negative_length')
                    pass
                vars[i].l = x[i]
            if vars[i].__class__ == Quadrupole:
                vars[i].k1 = x[i]
            if vars[i].__class__ == Solenoid:
                vars[i].k = x[i]
            if vars[i].__class__ in [RBend, SBend, Bend]:
                if vary_bend_angle:
                    vars[i].angle = x[i]
                else:
                    vars[i].k1 = x[i]
            if vars[i].__class__ == list:
                if vars[i][0].__class__ == Twiss and vars[i][1].__class__ == str:
                    k = vars[i][1]
                    setattr(tw_loc, k, x[i])
            if vars[i].__class__ == tuple:  # all quads strength in tuple varied simultaneously
                for v in vars[i]:
                    v.k1 = x[i]

        err = 0.0
        if "periodic" in constr.keys():
            if constr["periodic"] is True:
                tw_loc = lat.periodic_twiss(tw_loc)
                tw0 = deepcopy(tw_loc)
                if tw_loc is None:
                    print("########")
                    return weights('periodic')

        # save reference points where equality is asked

        ref_hsh = {}  # penalties on two-point inequalities

        for e in constr.keys():
            if e == 'periodic':
                continue
            if e == 'total_len':
                continue
            for k in constr[e].keys():
                if constr[e][k].__class__ == list:
                    if constr[e][k][0] == '->':
                        # print 'creating reference to', constr[e][k][1].id
                        ref_hsh[constr[e][k][1]] = {k: 0.0}

        # print 'references:', ref_hsh.keys()

        # evaluating global and point penalties

        # tw_loc.s = 0
        # print("start = ", get_envelope(p_array0))
        navi.go_to_start()
        tws_list, p_array0 = track(lat, p_array0, navi, print_progress=False, bounds=bounds, slice=slice)
        s = np.array([tw.s for tw in tws_list])
        # print("stop = ", tws_list[-1])
        L = 0.
        for e in lat.sequence:
            indx = (np.abs(s - L)).argmin()
            L += e.l
            tw_loc = tws_list[indx]

            # --- Global constraints ---
            if 'global' in constr:
                for c, rule in constr['global'].items():
                    if isinstance(rule, list):
                        op, v1 = rule[0], rule[1]
                        val = getattr(tw_loc, c)

                        if op == '<' and val > v1:
                            err += weights(c) * (val - v1) ** 2
                        elif op == '>' and val < v1:
                            err += weights(c) * (val - v1) ** 2

            # --- Store reference ---
            if e in ref_hsh:
                ref_hsh[e] = deepcopy(tw_loc)

            # --- Local constraints ---
            if e in constr:
                for k, rule in constr[e].items():
                    val = getattr(tw_loc, k)

                    if isinstance(rule, list):
                        op = rule[0]
                        v1 = rule[1]

                        if op == '<' and val > v1:
                            err += weights(k) * (val - v1) ** 2
                        elif op == '>' and val < v1:
                            err += weights(k) * (val - v1) ** 2
                        elif op == 'a<' and abs(val) > v1:
                            err += weights(k) * (abs(val) - v1) ** 2
                        elif op == 'a>' and abs(val) < v1:
                            err += weights(k) * (abs(val) - v1) ** 2
                        elif op == '->':
                            try:
                                dv1 = float(rule[2]) if len(rule) > 2 else 0.0
                                ref_val = getattr(ref_hsh[v1], k)
                                err += (val - (ref_val + dv1)) ** 2
                                if val < v1:
                                    err += (val - v1) ** 2
                            except Exception as ex:
                                print(f'Constraint error: rval should precede lval in lattice ({ex})')

                        if val < 0:
                            err += (val - v1) ** 2

                    else:
                        err += weights(k) * (rule - val) ** 2
        for v in vars:
            print(v.id, v.k1)
        if "total_len" in constr.keys():
            total_len = constr["total_len"]
            err = err + weights('total_len') * (tw_loc.s - total_len) ** 2

        if min_i5:
            # evaluating integral parameters

            I1, I2, I3, I4, I5 = radiation_integrals(lat, tw0, nsuperperiod=1)
            err += I5 * weights('i5')

            Je = 2 + I4 / I2
            Jx = 1 - I4 / I2
            Jy = 1

            if Je < 0 or Jx < 0 or Jy < 0:
                err = 100000.0

        # c1, c2 = natural_chromaticity(lat, tw0)
        # err += ( c1**2 + c2**2) * 1.e-6

        if verbose:
            print(f"Run Number: {run_number}/{max_iter}, iteration error: {err}" )
        run_number +=1

        return err


    # list of arguments determined based on the variable class

    x = [0.0] * len(vars)
    for i in range(len(vars)):
        if vars[i].__class__ == list:
            if vars[i][0].__class__ == Twiss and vars[i][1].__class__ == str:
                k = vars[i][1]
                if k in ['beta_x', 'beta_y']:
                    x[i] = 10.0
                else:
                    x[i] = 0.0
        if vars[i].__class__ == tuple:
            x[i] = vars[i][0].k1
        if vars[i].__class__ == Quadrupole:
            x[i] = vars[i].k1
        if vars[i].__class__ == Drift:
            x[i] = vars[i].l
        if vars[i].__class__ in [RBend, SBend, Bend]:
            # TODO: need a way to vary 2 aattributes of a class
            if vary_bend_angle:
                x[i] = vars[i].angle
            else:
                x[i] = vars[i].k1

    print("initial value: x = ", x)
    if method == 'simplex':
        res = fmin(errf, x, xtol=1e-3, maxiter=max_iter, maxfun=max_iter)
    if method == 'cg':
        res = fmin_cg(errf, x, gtol=1.e-5, epsilon=1.e-5, maxiter=max_iter)
    if method == 'bfgs':
        res = fmin_bfgs(errf, x, gtol=1.e-5, epsilon=1.e-5, maxiter=max_iter)
    if method == 'powell':
        res = minimize(errf, x, method='Powell', tol=1.e-5, options={"maxiter": max_iter})
    if method == "diff_evolution":
        workers = multiprocessing.cpu_count()
        bounds = []
        for xi in x:
            bounds.append((-5, 5))
        res = differential_evolution(errf, bounds, maxiter=max_iter, workers=1)

    # if initial twiss was varied set the twiss argument object to resulting value

    # for i in range(len(vars)):
    #    if vars[i].__class__ == list:
    #        if vars[i][0].__class__ == Twiss and vars[i][1].__class__ == str:
    #            k = vars[i][1]
    #            tw.__dict__[k] = res[i]


    # update MagneticLattice total length in case a Drift length was in list of variables


    return res


def match_matrix(lat, beam, varz, target_matrix):
    def error_func(x):

        for i in range(len(varz)):
            if varz[i].__class__ == Quadrupole:
                varz[i].k1 = x[i]
                varz[i].transfer_map = lat.method.create_tm(varz[i])  # create_transfer_map(varz[i])

        R = lattice_transfer_map(lat, beam.E)[0:2, 0:2]
        # print
        # R
        err = np.linalg.norm(np.abs(R - target_matrix) ** 2)

        # print
        # 'iteration error: ', err
        return err

    x = [0.0] * len(varz)

    for i in range(len(varz)):
        if varz[i].__class__ == Quadrupole:
            x[i] = varz[i].k1

    print("initial value: x = ", x)

    fmin(error_func, x, xtol=1e-8, maxiter=20000, maxfun=20000)


def match_tunes(lat, tw0, quads, nu_x, nu_y, ncells=1, max_iter=1000, tol=1e-5, print_proc=0):
    print("matching start .... ")
    end = Monitor(eid="end")
    lat = MagneticLattice(lat.sequence + [end])
    # tw0.E = lat.energy
    tws = twiss(lat, tw0, nPoints=None)

    nu_x_old = tws[-1].mux / 2 / np.pi * ncells
    nu_y_old = tws[-1].muy / 2 / np.pi * ncells
    # print nu_y, nu_y_old
    strengths1 = [p.k1 for p in quads]

    constr = {end: {'mux': 2 * np.pi * nu_x / ncells, 'muy': 2. * np.pi * nu_y / ncells}, 'periodic': True}
    # print constr
    vars = quads

    match(lat, constr, vars, tws[0], max_iter=max_iter, tol=tol)
    for i, q in enumerate(quads):
        print(q.id, ".k1: before: ", strengths1[i], "  after: ", q.k1)
    lat = MagneticLattice(lat.sequence[:-1])
    tws = twiss(lat, tw0, nPoints=None)
    print("nu_x: before: ", nu_x_old, "after: ", tws[-1].mux / 2 / np.pi * ncells)
    print("nu_y: before: ", nu_y_old, "after: ", tws[-1].muy / 2 / np.pi * ncells)
    print("matching end.")
    return lat


def closed_orbit(lattice, eps_xy=1.e-7, eps_angle=1.e-7, energy=0):
    __author__ = 'Sergey Tomin'

    """
    Searching of initial coordinates (p0) by iteration method.
    For initial conditions p uses exact solution of equation p = M*p + B
    :param lattice: class MagneticLattice
    :param eps_xy: tolerance on coordinates of beam in the start and end of lattice
    :param eps_angle: tolerance on the angles of beam in the start and end of lattice
    :return: class Particle
    """
    R = lattice_transfer_map(lattice, energy)
    smult = SecondOrderMult()

    ME = np.eye(4) - R[:4, :4]
    P = np.dot(np.linalg.inv(ME), lattice.B[:4])

    def errf(x):
        X = np.array([[x[0]], [x[1]], [x[2]], [x[3]], [0], [0]])
        smult.numpy_apply(X, R, lattice.T)
        X += lattice.B
        err = np.sum([1000 * (X[i, 0] - x[i]) ** 2 for i in range(4)])
        return err

    res = fmin(errf, P, xtol=1e-8, maxiter=2e3, maxfun=2.e3)

    return Particle(x=res[0], px=res[1], y=res[2], py=res[3])


def inverse_lattice(lat, phase=180):
    """
    Inverts a MagneticLattice for backtracking purposes.

    :param lat: MagneticLattice - the lattice to be inverted.
    :return: MagneticLattice - the inverted lattice.
    """
    if abs(phase) != 180:
        raise ValueError("phase must be +180 or -180")
    lat_inv = MagneticLattice(lat.sequence[::-1], method=lat.method)
    for elem in lat_inv.sequence:
        if isinstance(elem, Cavity):
            elem.phi += phase  # Adjust cavity phase for inversion
    return lat_inv


def inverse_twiss(tws):
    """
    Inverts the Twiss parameters for backtracking.

    :param tws: Twiss - Twiss object whose parameters need to be inverted.
    :return: Twiss - Twiss object with inverted alpha_x and alpha_y.
    """
    tws_inv = Twiss(tws)
    tws_inv.alpha_x *= -1
    tws_inv.alpha_y *= -1
    return tws_inv


def beam_matching(parray, navi, tws_end, vars, iter=3):
    """
    Matches the beam to desired Twiss parameters at the end of the lattice using
    an iterative approach based on linear tracking.

    :param parray: ParticleArray - particle array at the beginning of the lattice.
    :param navi: Navigator - navigator for the lattice.
    :param tws_end: Twiss - desired Twiss parameters at the end of the lattice.
    :param vars: list - list of adjustable variables (e.g., quadrupoles [q1, q2, q3, ...]).
    :param iter: int - number of iterations for the matching process (default: 3).
    :return: list - optimized variables after matching.
    """
    # Initial envelope based on the particle array
    tws_tmp = get_envelope(parray, bounds=[-2, 2])
    lat = navi.lat
    constr = {
        lat.sequence[-1]: {
            "beta_x": tws_end.beta_x,
            "beta_y": tws_end.beta_y,
            "alpha_x": tws_end.alpha_x,
            "alpha_y": tws_end.alpha_y,
        }
    }

    for i in range(iter):
        # Create a deep copy of the particle array for tracking
        parray_track = deepcopy(parray)

        # Perform linear matching to optimize variables
        res = match(lat, constr, vars, tws_tmp, verbose=False, max_iter=1000, tol=1e-5)

        # Track the beam through the lattice
        navi.reset_position()
        navi.lat = lat
        tws_track, parray_track = track(lat, p_array=parray_track, navi=navi)

        # Calculate Twiss parameters at the end of the lattice
        tws_end_tmp = get_envelope(parray_track, bounds=[-2, 2])

        # Backtrack Twiss parameters using the inverted lattice
        lat_inv = inverse_lattice(lat, phase=+180)
        tws_end_inv = inverse_twiss(tws_end_tmp)
        tws = twiss(lat_inv, tws_end_inv)
        tws_tmp = inverse_twiss(tws[-1])
        lat = inverse_lattice(lat_inv, phase=-180)

    return vars