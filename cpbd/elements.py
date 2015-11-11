'''
definition of magnetic lattice
linear dimensions in [m]
'''

from ocelot.cpbd.field_map import FieldMap
from ocelot.cpbd.optics import create_transfer_map
from ocelot.common.globals import *
import numpy as np
from numpy import cos, sin


flatten = lambda *n: (e for a in n
                            for e in (flatten(*a) if isinstance(a, (tuple, list)) else (a,)))

'''
TODO: rename Element class into something else to avoid confusion 
'''
class Element:
    
    ''' Element is a beamline building element with an arbitrary set of parameters attached '''
    
    def __init__(self, id = None):
        self.id = id
        if id == None:
            #print 'setting automatic id'
            self.id = "ID_{0}_".format(np.random.randint(100000000))
        self.type = "none"
        self.tilt = 0.       # TILT: The roll angle about the longitudinal axis
                             # (default: 0 rad, i.e. a normal quadrupole).
                             # A positive angle represents a clockwise rotation.
                             # A TILT=pi/4 turns a positive normal quadrupole into a negative skew quadrupole.

        self.angle = 0.
        self.k1 = 0.
        self.k2 = 0.
        self.dx = 0.
        self.dy = 0.
        self.dtilt = 0.
        self.params = {}
    
    def __hash__(self):
        return hash( (self.id, self.type) )

    def __eq__(self, other):
        try:
            return (self.id, type) == (other.id, type)
        except:
            return False
    

# to mark locations of bpms and other diagnostics
class Monitor(Element):
        
    def __init__(self, l=0.0, id = None):
        Element.__init__(self, id)
        self.type = "monitor"
        self.l = l

class Marker(Element):
    def __init__(self, id = None):
        Element.__init__(self, id)
        self.type = "marker"
        self.l = 0.

class Quadrupole(Element):
    """
    k1 - strength of quadrupole lens in [1/m^2],
    l - length of lens in [m].
    """
    def __init__(self, l=0, k1=0, k2=0., tilt=0, id=None):
        Element.__init__(self, id)
        self.type = "quadrupole"
        self.l = l
        self.k1 = k1
        self.k2 = k2
        self.tilt = tilt

class Sextupole(Element):
    """
    m - strength of sextupole lens in [1/m^3],
    l - length of lens in [m].
    """
    def __init__(self, l=0, k2=0., ms=0., id=None, tilt=0):
        """
        k2 is sextupole strength
        ms = k2*l
        """
        Element.__init__(self, id)
        self.type = "sextupole"
        if l != 0 and ms != 0:
            if k2 == 0:
                k2 = ms/l
                ms = 0.
            else:
                ms = 0.

        self.l = l
        self.k2 = k2
        self.ms = ms
        self.tilt = tilt

class Octupole(Element):
    """
    m - strength of sextupole lens in [1/m^3],
    l - length of lens in [m].
    """
    def __init__(self, l=0, k3=None, id=None, tilt=0.):
        """
        k2 is sextupole strength
        moct = k3*l
        """
        Element.__init__(self, id)
        self.type = "octupole"
        self.l = l
        self.k3 = k3
        self.tilt = tilt
        self.moct = None

class Drift(Element):
    """
    l - length of lens in [m]
    """
    def __init__(self, l=0, id = None):
        Element.__init__(self, id)
        self.type = "drift"
        self.l = l
        self.tilt = 0.

class Bend(Element):
    def __init__(self, l, angle=0., k1 = 0., k2 = 0., tilt=0.0, e1 = 0., e2 = 0.,
                 gap = 0, h_pole1 = 0., h_pole2 = 0., fint = 0., fintx=0., id = None):
        Element.__init__(self, id)
        self.type = "bend"
        self.l = l
        self.angle = angle
        self.k1 = k1
        self.k2 = k2
        self.e1 = e1
        self.e2 = e2
        self.gap = gap
        self.h_pole1 = h_pole1
        self.h_pole2 = h_pole2
        self.fint1 = fint
        self.fint2 = fint
        if fintx>0:
            self.fint2 = fintx
        self.tilt = tilt

class Edge(Bend):
    def __init__(self, l=0, angle=0.0, k1 = 0, edge = 0.,
                 tilt=0.0, dtilt = 0.0, dx = 0.0, dy = 0.0,
                 h_pole = 0., gap = 0., fint = 0., pos = 1, id = None):
        Element.__init__(self, id)
        self.type = "edge"
        if l!=0.:
            self.h = angle/l
        else:
            self.h = 0
        self.l = 0.
        self.k1 = k1
        self.h_pole = h_pole
        self.gap = gap
        self.fint = fint
        self.edge = edge
        self.dx = dx
        self.dy = dy
        self.dtilt = dtilt
        self.tilt = tilt
        self.pos = pos
        """
        m = matrix(eye(6))
        y = tan(angle)*hx
        m1 = y*cos(tilt)*cos(tilt) - y*sin(tilt)*sin(tilt)
        m2 = y*sin(2*tilt)
        m[1,0] = m1; m[1,2] = m2
        m[3,1] = m2; m[3,2] = -m1
        return m
        """
class SBend(Bend):
    """
    sector bending magnet,
    l - length of magnet,
    angle - angle of bend in [rad],
    k - quadrupole strength in [1/m^2].
    """
    def __init__(self, l=0, angle=0.0,k1 = 0.0, k2 = 0., e1 = 0.0, e2 = 0.0, tilt=0.0,
                 gap = 0, h_pole1 = 0., h_pole2 = 0., fint = 0., fintx=0., id = None):
        Bend.__init__(self, l, angle=angle, k1=k1, k2=k2, e1=e1, e2=e2,
                      gap=gap, h_pole1=h_pole1, h_pole2=h_pole2, fint=fint, id=id)
        self.type = "sbend"
        self.l = l
        self.angle = angle
        self.k1 = k1
        self.k2 = k2
        self.tilt = tilt
        self.gap = gap
        self.h_pole1 = h_pole1
        self.h_pole2 = h_pole2
        self.fint1 = fint
        self.fint2 = fint
        if fintx>0:
            self.fint2 = fintx
        # for future

        #self.e1 = 0     # The rotation angle for the entrance pole face (default: 0 rad).
        #self.e2 = 0     # The rotation angle for the exit pole face (default: 0 rad).

class RBend(Bend):
    """
    rectangular bending magnet,
    l - length of magnet,
    angle - angle of bend in [rad],
    k - quadrupole strength in [1/m^2].
    """
    def __init__(self, l=0, angle=0,tilt=0, k1 = 0, k2 = 0.,  e1 = None, e2 = None,
                 gap=0, h_pole1=0., h_pole2=0., fint=0., fintx=0., id=None):
        if e1 == None:
            e1 = angle/2.
        else:
            e1 += angle/2.
        if e2 == None:
            e2 = angle/2.
        else:
            e1 += angle/2.
        Bend.__init__(self, l, angle=angle, e1=e1, e2=e2, k1=k1, k2=k2,
                      gap=gap, h_pole1=h_pole1, h_pole2=h_pole2, fint=fint, fintx=fintx, id=id)
        self.type = "rbend"
        self.l = l
        self.angle = angle
        self.k1 = k1
        self.k2 = k2
        self.tilt = tilt
        self.gap = gap
        self.h_pole1 = h_pole1
        self.h_pole2 = h_pole2
        self.fint1 = fint
        self.fint2 = fint
        if fintx > 0:
            self.fint2 = fintx

class Hcor(RBend):
    def __init__(self,l = 0, angle = 0, id = None):
        RBend.__init__(self, l=l, angle=angle, id = id)
        self.type = "hcor"
        self.l = l
        self.angle = angle
        self.tilt = 0

class Vcor(RBend):
    def __init__(self,l = 0, angle = 0, id = None):
        RBend.__init__(self, l=l, angle=angle, id = id)
        self.type = "vcor"
        self.l = l
        self.angle = angle
        self.tilt = pi/2.

class Undulator(Element):
    """
    lperiod - undulator period in [m];\n
    nperiod - number of periods;\n
    Kx - undulator paramenter for vertical field; \n
    Ky - undulator parameter for horizantal field;\n
    field_file_path - absolute path to magnetic field data;\n
    id - name of undulator. 
    """
    def __init__(self, lperiod, nperiods, Kx, Ky=0, field_file=None, id=None):
        Element.__init__(self, id)
        self.type = "undulator"
        self.lperiod = lperiod
        self.nperiods = nperiods
        self.l = lperiod * nperiods
        self.Kx = Kx
        self.Ky = Ky
        self.solver = "linear"  # can be "lin" is liear matrix,  "sym" - symplectic method and "rk" is Runge-Kutta
        self.phase = 0 # phase between Bx and By + pi/4 (spiral undulator)
        
        self.ax = -1              # width of undulator, when ax is negative undulator width is infinite
                                  # I need it for analytic description of undulator 
        
        self.field_file = field_file
        self.field_map = FieldMap(self.field_file)
        self.v_angle = 0.
        self.h_angle = 0.
        #self.processing()  # here we can check all data and here we can load magnetic map from file
                            # and if error will appear then it will be on stage of forming structure
                            # more over I suggest to store all data (including field map) in this class
                            
    def validate(self):
        pass
                            
        # maybe we will do two functions
        # 1. load data and check magnetic map
        # 2. check all input data (lperiod nperiod ...). domething like this we must do for all elements.

        # what do you think about ending poles? We can do several options
        # a) 1/2,-1,1,... -1,1/2
        # b) 1/2,-1,1,... -1,1,-1/2
        # c) 1/4,-3/4,1,-1... -1,3/4,-1/4   I need to check it.

class Cavity(Element):
    '''
    RF cavity
    v - voltage [V/m]
    f - frequency [GHz]
    '''
    def __init__(self, l, delta_e=0.0, freq=0.0, phi=0.0, id=None, volt=0., volterr=0.):
        Element.__init__(self, id)
        self.type = "cavity"
        self.l = l
        self.v = volt*1e-9   #in GV
        self.delta_e = delta_e
        self.f = freq
        self.phi = phi*np.pi/180.
        self.E = 0
        self.volterr = volterr


class Solenoid(Element):
    '''
    Solenoid
    '''
    def __init__(self, l, k = 0., id = None):
        Element.__init__(self, id)
        self.type = "solenoid"
        self.k = k # B0/(2B*rho)
        self.l = l

class Multipole(Element):
    """
    kn - list of strengths
    """
    def __init__(self, kn=0., id=None):

        Element.__init__(self, id)
        self.type = "multipole"
        kn = np.array([kn]).flatten()
        if len(kn) < 2:
            self.kn = np.append(kn, 0.)
        else:
            self.kn = kn
        self.n = len(self.kn)
        self.l = 0.


class Matrix(Element):
    def __init__(self, l = 0., rm11 = 0., rm12=0., rm13 = 0., rm21=0.0, rm22=0.0, rm33 = 0., rm34 = 0., rm43 = 0., rm44 = 0., id = None):
        Element.__init__(self, id)
        self.type = "matrix"
        self.l = l
        self.rm11 = rm11
        self.rm12 = rm12
        self.rm13 = rm13
        self.rm21 = rm21
        self.rm22 = rm22
        self.rm33 = rm33
        self.rm34 = rm34
        self.rm43 = rm43
        self.rm44 = rm44

class UnknownElement(Element):
    """
    l - length of lens in [m]
    """
    def __init__(self, l=0, kick = 0,xsize = 0, ysize = 0, volt = 0, lag = 0, harmon = 0, refer = 0,vkick = 0,hkick = 0, id = None):
        Element.__init__(self, id)
        self.type = "drift"
        self.l = l

class Sequence:
    def __init__(self, l = 0, refer = 0):
        self.l = l

class MagneticLattice:
    def __init__(self, sequence, start=None, stop=None):
        #self.energy = energy
        self.sequence = list(flatten(sequence))

        try:
            if start != None: id1 = self.sequence.index(start)
            else: id1 = 0
            if stop != None:
                id2 = self.sequence.index(stop) + 1
                self.sequence = self.sequence[id1:id2]
            else:
                self.sequence = self.sequence[id1:]
        except:
            print 'cannot construct sequence, element not found'
            raise

        #self.transferMaps = {}
        # create transfer map and calculate lattice length
        self.totalLen = 0
        if not self.check_edges():
            self.add_edges()
        self.update_transfer_maps()

        self.__hash__ = {}
        #print 'creating hash'
        for e in self.sequence:
            #print e
            self.__hash__[e] = e
    
    def __getitem__(self, el):
        try:
            return self.__hash__[el]
        except:
            return None

    def check_edges(self):
        """
        if there are edges on the ends of dipoles return True, else False
        """
        for i in range(len(self.sequence)-2):
            prob_edge1 = self.sequence[i]
            elem = self.sequence[i+1]
            prob_edge2 = self.sequence[i+2]
            if elem.type in ["bend", "sbend", "rbend"]: # , "hcor", "vcor"
                if prob_edge1.type != "edge" and prob_edge2 != "edge":
                    #print elem.type, prob_edge1.type, prob_edge2.type
                    return False
        return True

    def add_edges(self):
        n = 0
        for i in range(len(self.sequence)):
            elem = self.sequence[n]
            if elem.type in ["bend", "sbend", "rbend"] and elem.l != 0.: # , "hcor", "vcor"

                e_name = elem.id

                if elem.id == None:
                    e_name = "b_" + str(i)

                e1 = Edge(l=elem.l, angle=elem.angle, k1=elem.k1, edge=elem.e1, tilt=elem.tilt, dtilt=elem.dtilt,
                          dx=elem.dx, dy=elem.dy, h_pole=elem.h_pole1, gap=elem.gap, fint=elem.fint1, pos=1,
                          id=e_name + "_e1")

                self.sequence.insert(n, e1)

                e2 = Edge(l=elem.l, angle=elem.angle, k1=elem.k1, edge=elem.e2, tilt=elem.tilt, dtilt=elem.dtilt,
                          dx=elem.dx, dy=elem.dy, h_pole=elem.h_pole2, gap=elem.gap, fint=elem.fint2, pos=2,
                          id=e_name + "_e2")

                self.sequence.insert(n+2, e2)
                n += 2
            n +=1

    def update_transfer_maps(self):
        #E = self.energy
        self.totalLen = 0
        for element in self.sequence:
            if element.type == "undulator":
                if element.field_file != None:
                    element.l = element.field_map.l * element.field_map.field_file_rep
                    if element.field_map.units =="mm":
                        element.l = element.l*0.001
            self.totalLen += element.l
            """
            try:
                element.transfer_map = create_transfer_map(element, energy=element.E, track_acceleration=track_acceleration)
            except:
                element.transfer_map = create_transfer_map(element, energy=self.energy, track_acceleration=track_acceleration)
            """
            #if element.type == "cavity":
            #    E += element.delta_e
            #print "init = ", E

            element.transfer_map = create_transfer_map(element)
        return self

    def printElements(self):
        print( '\nLattice\n')
        for e in self.sequence:
            print('-->',  e.id, '[', e.l, ']')


def survey(lat, ang = 0.0, x0=0, z0=0):
    x = []
    z = []
    for e in lat.sequence:
        x.append(x0)
        z.append(z0)
        if e.__class__ in [Bend, SBend, RBend]:
            ang += e.angle
        x0 += e.l*cos(ang)  
        z0 += e.l*sin(ang)
    return x, z



if __name__ == "__main__":

    fm1 = FieldMap(field_file="center.dat", format = "tabular")
    fm2 = FieldMap(field_file="epu49cen.dat", format = "flat")
    print( fm1.z_arr)
    print(fm2.z_arr )
