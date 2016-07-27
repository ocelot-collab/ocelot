from scipy.optimize import *

from ocelot.cpbd.e_beam_params import radiation_integrals
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.cpbd.optics import *



def weights_default(val):
    if val == 'periodic': return 10000001.0
    if val == 'total_len': return 10000001.0
    if val == 'Dx': return 10000002.0
    if val == 'Dxp': return 10000003.0
    if val == 'tau': return 10000004.0
    if val == 'i5': return 1.e14
    if val == 'negative_length': return 1.5e6
    if val in ['alpha_x', 'alpha_y']: return 100005.0
    if val in ['mux', 'muy']: return 10000006.0
    if val in ['beta_x', 'beta_y']: return 100007.0
    return 0.0001


def match(lat, constr, vars, tw, verbose=True, max_iter=1000, method='simplex', weights=weights_default,
          vary_bend_angle=False, min_i5=False):
    '''
    matching stuff
    '''
    # tw = deepcopy(tw0)

    def errf(x):

        tw_loc = deepcopy(tw)
        tw0 = deepcopy(tw)

        '''
        parameter to be varied is determined by variable class
        '''
        for i in range(len(vars)):
            if vars[i].__class__ == Drift:
                if x[i] < 0:
                    # print('negative length in match')
                    return weights('negative_length')
                    pass
                vars[i].l = x[i]
                vars[i].transfer_map = lat.method.create_tm(vars[i])
            if vars[i].__class__ == Quadrupole:
                vars[i].k1 = x[i]
                vars[i].transfer_map = lat.method.create_tm(vars[i])
            if vars[i].__class__ in [RBend, SBend, Bend]:
                if vary_bend_angle:
                    vars[i].angle = x[i]
                else:
                    vars[i].k1 = x[i]
                vars[i].transfer_map = lat.method.create_tm(vars[i])
            if vars[i].__class__ == list:
                if vars[i][0].__class__ == Twiss and vars[i][1].__class__ == str:
                    k = vars[i][1]
                    tw_loc.__dict__[k] = x[i]
            if vars[i].__class__ == tuple: # all quads strength in tuple varied simultaneously 
                for v in vars[i]:
                    v.k1 = x[i]
                    v.transfer_map = lat.method.create_tm(v)
                     

        err = 0.0
        if "periodic" in constr.keys():
            if constr["periodic"] == True:
                tw_loc = periodic_twiss(tw_loc, lattice_transfer_map(lat, tw.E))
                tw0 = deepcopy(tw_loc)
                if tw_loc == None:
                    print("########")
                    return weights('periodic')


        # save reference points where equality is asked

        ref_hsh = {}  # penalties on two-point inequalities

        for e in constr.keys():
            if e == 'periodic': continue
            if e == 'total_len': continue
            for k in constr[e].keys():
                if constr[e][k].__class__ == list:
                    if constr[e][k][0] == '->':
                        # print 'creating reference to', constr[e][k][1].id
                        ref_hsh[constr[e][k][1]] = {k: 0.0}

        # print 'references:', ref_hsh.keys()

        ''' evaluating global and point penalties
        '''
                        
        tw_loc.s = 0
                        
        for e in lat.sequence:
            #print e.id,  e in constr.keys()
            # print tw_loc
            #print '--->'
            tw_loc = e.transfer_map * tw_loc
            # print tw_loc

            if 'global' in constr.keys():
                # print 'there is a global constraint', constr['global'].keys()
                for c in constr['global'].keys():
                    if constr['global'][c].__class__ == list:
                        # print 'list'
                        v1 = constr['global'][c][1]
                        if constr['global'][c][0] == '<':
                            if tw_loc.__dict__[c] > v1:
                                err = err + weights(k) * (tw_loc.__dict__[c] - v1) ** 2
                        if constr['global'][c][0] == '>':
                            # print '> constr'
                            if tw_loc.__dict__[c] < v1:
                                err = err + weights(k) * (tw_loc.__dict__[c] - v1) ** 2

            if e in ref_hsh.keys():
                # print 'saving twiss for', e.id
                ref_hsh[e] = deepcopy(tw_loc)

            if e in constr.keys():

                for k in constr[e].keys():
                    #print k
                    if constr[e][k].__class__ == list:                    
                        v1 = constr[e][k][1]

                        if constr[e][k][0] == '<':
                            if tw_loc.__dict__[k] > v1:
                                err = err + weights(k) * (tw_loc.__dict__[k] - v1) ** 2
                        if constr[e][k][0] == '>':
                            if tw_loc.__dict__[k] < v1:
                                err = err + weights(k) * (tw_loc.__dict__[k] - v1) ** 2
                        if constr[e][k][0] == 'a<':
                            if np.abs(tw_loc.__dict__[k]) > v1:
                                err = err + weights(k) * (tw_loc.__dict__[k] - v1) ** 2
                        if constr[e][k][0] == 'a>':
                            if np.abs(tw_loc.__dict__[k]) < v1:
                                err = err + weights(k) * (tw_loc.__dict__[k] - v1) ** 2


                        if constr[e][k][0] == '->':
                            try:
                                #print 'huh', k, e.id, float(constr[e][k][2])
                                
                                if len(constr[e][k]) > 2:
                                    dv1 = float(constr[e][k][2])
                                else:
                                    dv1 = 0.0
                                #print 'weiter'
                                err += (tw_loc.__dict__[k] - (ref_hsh[v1].__dict__[k]  + dv1 ) ) ** 2

                                if tw_loc.__dict__[k] < v1:
                                    err = err + (tw_loc.__dict__[k] - v1) ** 2
                            except:
                                print('constraint error: rval should precede lval in lattice')

                        if tw_loc.__dict__[k] < 0:
                            # print 'negative constr (???)'
                            err += (tw_loc.__dict__[k] - v1) ** 2

                    else:
                        # print "safaf", constr[e][k] , tw_loc.__dict__[k], k, e.id, x
                        err = err + weights(k) * (constr[e][k] - tw_loc.__dict__[k]) ** 2
                        # print err
        
        if "total_len" in constr.keys():
            total_len = constr["periodic"]
            err = err + weights('total_len')*(tw_loc.s - total_len)**2
        

        if min_i5:
            ''' evaluating integral parameters
            '''
            I1, I2, I3, I4, I5 = radiation_integrals(lat, tw0, nsuperperiod=1)
            err += I5 * weights('i5')

            Je = 2 + I4 / I2
            Jx = 1 - I4 / I2
            Jy = 1

            if Je < 0 or Jx < 0 or Jy < 0: err = 100000.0

        # c1, c2 = natural_chromaticity(lat, tw0)
        # err += ( c1**2 + c2**2) * 1.e-6

        if verbose:
            print('iteration error:', err)
        return err

    '''
    list of arguments determined based on the variable class
    '''
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
    if method == 'simplex': res = fmin(errf, x, xtol=1e-4, maxiter=max_iter, maxfun=max_iter)
    if method == 'cg': res = fmin_cg(errf, x, gtol=1.e-5, epsilon=1.e-5, maxiter=max_iter)
    if method == 'bfgs': res = fmin_bfgs(errf, x, gtol=1.e-5, epsilon=1.e-5, maxiter=max_iter)

    '''
    if initial twiss was varied set the twiss argument object to resulting value
    '''
    for i in range(len(vars)):
        if vars[i].__class__ == list:
            if vars[i][0].__class__ == Twiss and vars[i][1].__class__ == str:
                k = vars[i][1]
                tw.__dict__[k] = res[i]
    return res


def match_matrix(lat, beam, varz, target_matrix):
    def error_func(x):

        for i in range(len(varz)):
            if varz[i].__class__ == Quadrupole:
                varz[i].k1 = x[i]
                varz[i].transfer_map = lat.method.create_tm(varz[i]) # create_transfer_map(varz[i])

        R = lattice_transfer_map(lat, beam.E)[0:2, 0:2]
        #print
        #R
        err = np.linalg.norm(np.abs(R - target_matrix) ** 2)

        #print
        #'iteration error: ', err
        return err

    x = [0.0] * len(varz)

    for i in range(len(varz)):
        if varz[i].__class__ == Quadrupole:
            x[i] = varz[i].k1

    print("initial value: x = ", x)

    fmin(error_func, x, xtol=1e-8, maxiter=20000, maxfun=20000)


def match_tunes(lat, tw0, quads, nu_x, nu_y, ncells=1, print_proc=0):
    print("matching start .... ")
    end = Monitor(eid="end")
    lat = MagneticLattice(lat.sequence + [end])
    # tw0.E = lat.energy
    tws = twiss(lat, tw0, nPoints=None)

    nu_x_old = tws[-1].mux / 2 / pi * ncells
    nu_y_old = tws[-1].muy / 2 / pi * ncells
    # print nu_y, nu_y_old
    strengths1 = [p.k1 for p in quads]

    constr = {end: {'mux': 2 * pi * nu_x / ncells, 'muy': 2. * pi * nu_y / ncells}, 'periodic': True}
    # print constr
    vars = quads

    match(lat, constr, vars, tws[0])
    for i, q in enumerate(quads):
        print(q.id, ".k1: before: ", strengths1[i], "  after: ", q.k1)
    lat = MagneticLattice(lat.sequence[:-1])
    tws = twiss(lat, tw0, nPoints=None)
    print("nu_x: before: ", nu_x_old, "after: ", tws[-1].mux / 2 / pi * ncells)
    print("nu_y: before: ", nu_y_old, "after: ", tws[-1].muy / 2 / pi * ncells)
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
    #R = lattice_transfer_map(lattice, energy)
    navi = Navigator()
    t_maps = get_map(lattice, lattice.totalLen, navi)

    tm0 = TransferMap()
    for tm in t_maps:
        tm0 = tm * tm0
        #if tm.order != 2:
        #    tm0 = tm * tm0
        #else:
        #    sex = TransferMap()
        #    sex.R[0, 1] = tm.length
        #    sex.R[2, 3] = tm.length
        #    tm0 = sex * tm0

    R = tm0.R(energy)[:4, :4]

    ME = eye(4) - R
    P = dot(inv(ME), tm0.B(energy)[:4])

    def errf(x):

        p = Particle(x=x[0], px=x[1], y=x[2], py=x[3])
        for tm in t_maps:
            p = tm * p
        err = 1000. * (p.x - x[0]) ** 2 + 1000. * (p.px - x[1]) ** 2 + 1000. * (p.y - x[2]) ** 2 + 1000. * (p.py - x[
            3]) ** 2

        return err

    res = fmin(errf, P, xtol=1e-8, maxiter=2.e3, maxfun=2.e3)

    return Particle(x=res[0], px=res[1], y=res[2], py=res[3])