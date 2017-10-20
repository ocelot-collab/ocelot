from ocelot import MagneticLattice

__author__ = 'Sergey Tomin'
from ocelot.lib.genera.src.python.trajectory.undulator import und_trace
from numpy import append
from ocelot.cpbd.optics import *
from ocelot.cpbd.elements import *
from ocelot.lib.genera.src.python.trajectory.motion import Motion
from scipy.interpolate import splrep, splev


def split_lat(lattice):
    cells = []
    cell = []
    for elem in lattice.sequence:
        if elem.__class__ == Undulator:
            if len(cell)>0:
                cells.append(MagneticLattice(cell))
            cells.append(MagneticLattice(elem))
            cell = []
        else:
            cell.append(elem)
            #if len(cell)>0:
                #cells.append(MagneticLattice(cell, energy = lattice.energy))
            #cells.append(MagneticLattice([elem], energy = lattice.energy))
            #cell = []
    if len(cell)>0:
        cells.append(MagneticLattice(cell))

    return cells

def particle_end(motion, p_start):
    m = motion
    p = Particle(x=m.X[-1]/1000., y=m.Y[-1]/1000., px=m.Xbeta[-1], py=m.Ybeta[-1], s=m.Z[-1]/1000.,
                 p=p_start.p,  tau=p_start.tau, E=p_start.E)
    return p


def particles2motion(part_list):
    #print map(lambda p:p.s, part_list)
    motion = Motion()
    motion.X = map(lambda p:p.x, part_list)
    motion.Y = map(lambda p:p.y, part_list)
    motion.Z = map(lambda p:p.s, part_list)
    motion.Xbeta= map(lambda p:p.px, part_list)
    motion.Ybeta = map(lambda p:p.py, part_list)
    return motion


def choice_dz(lat):
    return lat.totalLen/10000.


def integration_beta2(s, px, py):

    px2 = px*px
    py2 = py*py
    h = 2.*(s[1]-s[0])
    Ix = 0.
    Iy = 0.
    list_Ix = [Ix]
    list_Iy = [Iy]
    Is = s[0]
    for i in xrange(0,len(px2)-2,2):
        Ix += (px2[i] + 4.*px2[i+1] + px2[i+2])*h/6.
        Iy += (py2[i] + 4.*py2[i+1] + py2[i+2])*h/6.
        list_Ix = append(list_Ix, Ix)
        list_Iy = append(list_Iy, Iy)
        Is = append(Is, s[i+2])

    return Is, list_Ix, list_Iy

def coordinate_transform(x, y, xnew):

    tck = splrep(x,y)
    ynew = splev(xnew,tck,der=0)
    return ynew

def x2x(X):
    sqrt35 = 0.5*sqrt(3./5.)

    xnew = [X[0]]
    h = X[1] - X[0]
    xgaus = array([0.5-sqrt35, 0.5-sqrt35 + sqrt35, 0.5-sqrt35 + sqrt35 + sqrt35])*h
    #print h, xgaus
    for x in X:
        for xg in xgaus:
            xnew = append(xnew, x + xg)
    return xnew

def plist2arrays(part_list):
    x = array(map(lambda f: f.x, part_list))
    y = array(map(lambda f: f.y, part_list))
    s = array(map(lambda f: f.s, part_list))
    px = array(map(lambda f: f.px, part_list))
    py = array(map(lambda f: f.py, part_list))
    pz = speed_of_light - (px*px + py*py)*0.5
    return x,y,s,px,py,pz

def trace4radiation(lat,particle0, accuracy = 1):
    lattices = split_lat(lat)

    motions = []
    particle = particle0
    for lat in lattices:
        z = 0.
        #part_list = []
        if lat.sequence[0].__class__ == Undulator:
            undulator = lat.sequence[0]
            #print "energy = ", particle.E
            motion = und_trace(undulator, particle, energy = particle.E, bRough = 0, n_trajectory_points = None, accuracy = accuracy)
            particle = particle_end(motion, particle)
        else:
            dz = 0.05/accuracy #choise_dz(lat)
            part_list = trace_obj(lat, particle, nPoints = int(lat.totalLen/dz)*2+1)

            x,y,s,px,py,pz = plist2arrays(part_list)
            #print "sadfas", x,y,s,px,py,pz
            Is, Ipx2, Ipy2 = integration_beta2(s, px, py)

            s_new = x2x(Is)
            x = coordinate_transform(x = s, y = x, xnew = s_new)
            y = coordinate_transform(x = s, y = y, xnew = s_new)
            px = coordinate_transform(x = s, y = px, xnew = s_new)
            py = coordinate_transform(x = s, y = py, xnew = s_new)
            pz = coordinate_transform(x = s, y = pz, xnew = s_new)
            Ipx2 = coordinate_transform(x = Is, y = Ipx2, xnew = s_new)
            Ipy2 = coordinate_transform(x = Is, y = Ipy2, xnew = s_new)

            motion = Motion(N = len(s_new))

            motion.X[:] = x[:]*1000.
            motion.Y[:] = y[:]*1000.
            motion.Z[:] = s_new[:]*1000.
            motion.Xbeta[:] = px[:]
            motion.Ybeta[:] = py[:]
            motion.Zbeta[:] = pz[:]
            motion.XbetaI2[:] = Ipx2[:]
            motion.YbetaI2[:] = Ipy2[:]

            #plt.plot(motion.memory_motion[2*motion.N:3*motion.N], motion.memory_motion[9*motion.N:10*motion.N], "r.-", s_new, Ipy2,"b.-")
            #plt.show()
            #motion = prepare_gauss(motion)
            #print "z = ", z, dz
            particle = particle_end(motion, particle)
        motions.append(motion)
    return motions


