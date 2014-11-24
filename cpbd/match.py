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

    def errf(x):

        p = Particle(x = x[0], px = x[1], y = x[2], py = x[3])
        for tm in t_maps:
            p = tm*p
        err = 1000.*(p.x - x[0])**2 + 1000.*(p.px - x[1])**2 + 1000.*(p.y - x[2])**2 + 1000.*(p.py - x[3])**2

        return err

    res = fmin(errf,P,xtol=1e-8, maxiter=2.e3, maxfun=2.e3)

    return Particle(x = res[0], px = res[1], y = res[2], py = res[3])


def weights(val):
    if val == 'Dx': return 100.0
    if val == 'Dxp': return 100.0
    if val == 'tau': return 10000000.0
    
    return 0.0001

def match(lat, constr, vars, tw, print_proc = 1, max_iter=1000):
    
    #tw = deepcopy(tw0)
    
    def errf(x):
        
        #print 'iteration'
        tw_loc = deepcopy(tw)
        for i in xrange(len(vars)):
            if vars[i].__class__ == Quadrupole:
                vars[i].k1 = x[i]
                vars[i].transfer_map = create_transfer_map(vars[i], energy = tw.E)
            if vars[i].__class__ in [RBend, SBend, Bend]:
                vars[i].angle = x[i]
                vars[i].transfer_map = create_transfer_map(vars[i], energy = tw.E)
            if vars[i].__class__ == list:
                if vars[i][0].__class__ == Twiss and  vars[i][1].__class__ == str:
                    k = vars[i][1]
                    tw_loc.__dict__[k] = x[i]
        
        err = 0.0
        if "periodic" in constr.keys():
            if constr["periodic"] == True:
                tw_loc = periodic_solution(tw_loc, lattice_transfer_map(lat).R)
                if tw_loc == None:
                    return 1.0
        
        # save reference points where equality is asked   
        
        ref_hsh = {} # penalties on two-point inequalities
             
        for e in constr.keys():
            for k in constr[e].keys():
                if constr[e][k].__class__ == list:
                    if constr[e][k][0] == '->':
                        #print 'creating reference to', constr[e][k][1].id
                        ref_hsh[constr[e][k][1]] = {k:0.0}
                
        #print 'references:', ref_hsh.keys()
                
        #print "after: ", tw_loc.beta_x
        for e in lat.sequence:
            #print e.id, constr.keys()
            tw_loc = e.transfer_map*tw_loc
            
            if 'global' in constr.keys():
                #print 'there is a global constraint', constr['global'].keys()  
                for c in constr['global'].keys():
                    if constr['global'][c].__class__ == list:
                        #print 'list'   
                        v1 = constr['global'][c][1]
                        if constr['global'][c][0] == '<':
                            #print '< constr'                             
                            if tw_loc.__dict__[c] > v1:
                                err = err + (tw_loc.__dict__[c] - v1)**2
                        if constr['global'][c][0] == '>':
                            #print '> constr' 
                            if tw_loc.__dict__[c] < v1:
                                err = err + (tw_loc.__dict__[c] - v1)**2

            
            if e in ref_hsh.keys():
                #print 'saving twiss for', e.id 
                ref_hsh[e] = deepcopy(tw_loc)
            
            if e in constr.keys():
                #print 'match point', e.id
                #print tw_loc.Dx
                #print constr[e.id]['Dx']
                                
                for k in constr[e].keys():
                    #print 'evaluation debug: constr ', e.id, k
                                        
                    if constr[e][k].__class__ == list:
                        #print 'list'   
                        v1 = constr[e][k][1]
                        if constr[e][k][0] == '<':
                            #print '< constr'                             
                            if tw_loc.__dict__[k] > v1:
                                err = err + (tw_loc.__dict__[k] - v1)**2
                        if constr[e][k][0] == '>':
                            #print '> constr' 
                            if tw_loc.__dict__[k] < v1:
                                err = err + (tw_loc.__dict__[k] - v1)**2

                        if constr[e][k][0] == '->':
                            #print '-> constr'
                            try:
                                #print 'ref constraint:'
                                #print k
                                #print v1 
                                #print ref_hsh[v1]
                                #print ref_hsh[v1].__dict__.keys()
                                #print ref_hsh[v1].__dict__[k], tw_loc.__dict__[k]

                                err += (tw_loc.__dict__[k] - ref_hsh[v1].__dict__[k])**2

                                if tw_loc.__dict__[k] < v1:
                                    err = err + (tw_loc.__dict__[k] - v1)**2
                            except:
                                print 'constraint error: rval should precede lval in lattice'

                        if tw_loc.__dict__[k] < 0:
                            #print 'negative constr (???)'
                            err += (tw_loc.__dict__[k] - v1)**2
                            
                    else:
                        #print "safaf", constr[e.id][k] , tw_loc.__dict__[k], k
                        err = err + weights(k) * (constr[e][k] - tw_loc.__dict__[k])**2
                        #print err

        if print_proc == 1:
            print 'iteration error:', err
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
        if vars[i].__class__ == Quadrupole:
            x[i] = vars[i].k1
        if vars[i].__class__ in [RBend, SBend, Bend] :
            x[i] = vars[i].angle

    print "initial value: x = ", x
    res  = fmin(errf,x,xtol=1e-8, maxiter=max_iter, maxfun=max_iter)
    #print res
    for i in xrange(len(vars)):
        if vars[i].__class__ == list:
            if vars[i][0].__class__ == Twiss and  vars[i][1].__class__ == str:
                k = vars[i][1]
                tw.__dict__[k] = res[i]
    return res
    
    #print res

def math_tunes(lat, tw0, quads,  nu_x, nu_y, ncells= 1, print_proc = 0):
    print "matching start .... "
    end = Monitor(id = "end")
    lat = MagneticLattice(lat.sequence + [end], energy = lat.energy)
    tws=twiss(lat, tw0,nPoints=None)

    nu_x_old = tws[-1].mux/2/pi * ncells
    nu_y_old = tws[-1].muy/2/pi * ncells

    strengths1 = map(lambda p: p.k1, quads)

    constr = {'end':{'mux':2*pi*nu_x/ncells, 'muy':2.*pi*nu_y/ncells},'periodic':True}
    vars = quads

    match(lat, constr, vars, tws[0], print_proc = print_proc)
    for i, q in enumerate(quads):
        print q.id, ".k1: before: ",strengths1[i], "  after: ", q.k1
    lat = MagneticLattice(lat.sequence[:-1], energy = lat.energy)
    tws=twiss(lat, tw0,nPoints=None)
    print "nu_x: before: ", nu_x_old, "after: ", tws[-1].mux/2/pi * ncells
    print "nu_y: before: ", nu_y_old, "after: ", tws[-1].muy/2/pi * ncells
    print "matching end."
    return lat