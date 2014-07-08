from ocelot.cpbd.optics import *
from ocelot.cpbd.elements import *

from copy import deepcopy
from scipy.optimize import *


def closed_orbit(lattice, eps_xy = 1.e-7, eps_angle = 1.e-7):
    __author__ = 'Sergey Tomin'
    
    """
    Searching of initial coordinates (p0) by iteration method.
    For initial conditions p uses exact solution of equation p = M*p + B
    :param lattice: class MagneticLattice
    :param eps_xy: tolerance on coordinates of beam in the start and end of lattice
    :param eps_angle: tolerance on the angles of beam in the start and end of lattice
    :return: class Particle
    """
    navi = Navigator(lattice)
    t_maps = get_map(lattice, lattice.totalLen, navi)

    tm0 = TransferMap()
    for tm in t_maps:
        if tm.order!=2:
            tm0 = tm*tm0
        else:
            sex = TransferMap()
            sex.R[0,1] = tm.length
            sex.R[2,3] = tm.length
            tm0 = sex*tm0

    R = tm0.R[:4,:4]

    ME = eye(4) - R
    P = dot(inv(ME), tm0.B[:4])

    p0 = Particle(x = P[0], px = P[1], y = P[2], py = P[3])

    for i in xrange(100):

        p = p0
        for tm in t_maps:
            p = tm*p
        d_x = p.x - p0.x
        d_y = p.y - p0.y

        d_angle_x = p.px - p0.px
        d_angle_y = p.py - p0.py
        if (abs(d_angle_x) <= eps_angle and
            abs(d_angle_y) <= eps_angle and
            abs(d_x) <= eps_xy and
            abs(d_y) <= eps_xy):

            print "closed orbit is found! Number iterations = ", i
            return p0
        n = 2.
        p0.px += d_angle_x/n
        p0.py += d_angle_y/n
        p0.x += d_x/n
        p0.y += d_y/n
    print "closed orbit is not found!"

    return p0


def weights(val):
    if val == 'Dx': return 10.0
    if val == 'Dxp': return 10.0
    
    return 0.0001

def match(lat, constr, vars, tw):
    
    #tw = deepcopy(tw0)
    
    def errf(x):
        
        #print 'iteration'
        tw_loc = deepcopy(tw)

        for i in xrange(len(vars)):
            if vars[i].__class__ == Quadrupole:
                vars[i].k1 = x[i]
                vars[i].transfer_map = create_transfer_map(vars[i], energy = tw.E)
            if vars[i].__class__ == list:
                if vars[i][0].__class__ == Twiss and  vars[i][1].__class__ == str:
                    k = vars[i][1]
                    tw_loc.__dict__[k] = x[i]
        
        err = 0.0

        for e in lat.sequence:
            #print e.id
            #print e.transfer_map
            #print tw
            #print 'k', e.transfer_map*tw
            tw_loc = e.transfer_map*tw_loc
    
            if e.id in constr.keys():
                #print e.id
                #print tw_loc.Dx
                #print constr[e.id]['Dx']
                for k in constr[e.id].keys():
                    if constr[e.id][k].__class__ == list:
                        print 'list'   
                        v1 = constr[e.id][k][1]
                        if constr[e.id][k][0] == '<':
                            print '< constr'                             
                            if tw_loc.__dict__[k] > v1:
                                err = err + (tw_loc.__dict__[k] - v1)**2
                        if constr[e.id][k] == '>':
                            if tw_loc.__dict__[k] < v1:
                                err = err + (tw_loc.__dict__[k] - v1)**2

                        if tw_loc.__dict__[k] < 0:
                            err += (tw_loc.__dict__[k] - v1)**2
                            
                    else:
                        err = err + weights(k) * (constr[e.id][k] - tw_loc.__dict__[k])**2
        print err
        return err

    
    x = [0.0]*len(vars)
    for i in xrange(len(vars)):
        if vars[i].__class__ == list:
            if vars[i][0].__class__ == Twiss and  vars[i][1].__class__ == str:
                k = vars[i][1]
                if k in ['beta_x', 'beta_y']:
                    x[i] = 10.0
                else:
                    x[i] = 0.0

    
    res = fmin(errf,x,xtol=1e-8, maxiter=2.e3, maxfun=2.e3)

    for i in xrange(len(vars)):
        if vars[i].__class__ == list:
            if vars[i][0].__class__ == Twiss and  vars[i][1].__class__ == str:
                k = vars[i][1]
                tw.__dict__[k] = res[i]
    
    #print res
