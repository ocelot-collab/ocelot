'''
optics elements
'''

from numpy import sin, cos, pi, sqrt, log, array, random, sign
#from numpy.linalg import norm
import numpy as np


def debug(*args):
    print(' '.join( map(str, args) ))

def warn(*args):
    print(' '.join( map(str, args) ))

def info(*args):
    print(' '.join( map(str, args) ))


class OptElement:
    
    ''' OptElement is an optics building element with an arbitrary set of parameters attached '''
    
    def __init__(self, r, roll_ang, pitch_ang, yaw_ang, size, id):
        '''
        the roll pitch and yaw angles are not the standard aircraft notation
        yaw axis is always Y
        alternatively, one can specify the normal and z_ax explicitly after element is initialized
        see code for definition.
        '''
        self.id = id
        if id == None:
            #print 'setting automatic id'
            self.id = "ID_{0}_".format(np.random.randint(100000000))
        self.type = "none"
        self.r = np.array(r)
        self.size = size

        # outwards (from mirror surface) looking normal
        # default -- in y direction (elements in xz plane)
        no = np.matrix([[0., 1., 0.]]).T
        
        z_ax = np.matrix([[0., 0., 1.]]).T
        
        self.pitch_ang = pitch_ang
        M1 = np.matrix([[1., 0., 0.],
                        [0., cos(pitch_ang), sin(pitch_ang)],
                        [0., -sin(pitch_ang), cos(pitch_ang)]])
        #if pitch_ang 
        
        no = M1 * no
        z_ax = M1 * z_ax


        self.yaw_ang = yaw_ang
        M2 = np.matrix([[cos(yaw_ang), 0.,  sin(yaw_ang)],
                        [0.,           1.,  0.],
                        [-sin(yaw_ang),0.,  cos(yaw_ang)]])
        #if pitch_ang 
        
        no = M2 * no
        z_ax = M2 * z_ax


        
        self.roll_ang = roll_ang

        c = cos(roll_ang)
        s = sin(roll_ang)
        x, y, z = np.array(z_ax.T).reshape(-1) #no[0][0], no[0][1], no[0][2]
        #print x, y, z
        #print no
        
        M3 = np.matrix([[c + x**2*(1. - c), x*y*(1.-c) - z*s, x*z*(1.-c) + y*s],
                        [y*x*(1.-c) + z*s, c + y**2*(1.-c), y*z*(1.-c) - x*s],
                        [z*x*(1.-c) - y*s, z*y*(1.-c) + x*s , c + z**2*(1. - c)]])
        
        
        #print M3
        
        no = M3 * no
        z_ax = M3 * z_ax
        
        self.no = np.array(no.T).reshape(-1)
        self.z_ax = z_ax     


            
    def __hash__(self):
        return hash( (self.id, self.type) )

    def __eq__(self, other):
        try:
            return (self.id, type) == (other.id, type)
        except:
            return False


class Aperture:
    """
    all optical element sizes are for circular aperture
    add Aperture element in front to model more complex shapes/apertures
    """
    def __init__(self, r=[0,0,0], d=[0,0], no = [0,0,1], size=[1,1,0.1], id=None, type = "circular"): 
        self.r = r
        self.no = no
        self.d = d
        self.size = size
        self.type = type  # circular, rectangular, custom 
        self.id = id
    
    def slit(self,x,y):
        pass

class Lense:
    def __init__(self, r=[0,0,0], no=[0,0,1], s1=0.0, s2=0.0, D=1.0, n=0.99, id=""):
        self.n = n        # refractive index
        self.r = r        # position 
        self.s1 = s1      # sagitta (left/right)
        self.s2 = s2      # 
        self.D = D        # diameter
        self.size = [D/2.0,D/2.0,(s1+s2)/2.0]
        self.f = 1.0
        self.no = no
        self.id = id
        
    def c_fun(self,r):
        c = 1.0

        b1 = (r[2] >= self.r[2] - self.s1 + self.s1 / (self.D/2.0)**2 * r[1]**2 )
        b2 = (r[2] <= self.r[2] + self.s2 - self.s2 / (self.D/2.0)**2 * r[1]**2 )
        #print b1, b2
        if b1 and b2:
            #print 'LENSE'
            return 1.0 / self.n
        else: 
            return c


class Mirror(OptElement):
    """
    plane mirror
    """
    def __init__(self, r=[0,0,0], roll_ang=0.0, pitch_ang=0.0, yaw_ang=0.0, size=[1.0,1.0,0.1], id=None):
        OptElement.__init__(self, r, roll_ang, pitch_ang, yaw_ang, size, id)
        

class EllipticMirror(OptElement):
    def __init__(self, r=[0,0,0], roll_ang=0.0, pitch_ang=0.0, yaw_ang=0.0, size=[1.0,1.0,0.2], a=[1,1], id=""):
        OptElement.__init__(self, r, roll_ang, pitch_ang, yaw_ang, size, id)
        self.a = np.array(a) # major/minor axis of the ellipse


class ParabolicMirror:
    def __init__(self, r=[0,0,0], no=[0,0,1], a=[1,1], size=[1.0,1.0,0.1]):
        self.r = r
        self.no = no   # outwards (from mirror surface) looking normal
        self.a = a     # coefficient of the parabola z=a[0]*x**2 + a[1]*y**2    
        self.size = size

class Grating:
    def __init__(self, r=[0,0,0], no=[0,1,0], size=[1.0,1.0,0.1], d=0, id=""):
        self.r = r
        self.no = no
        self.size=size
        self.d = d
        self.id = id

class Crystal:
    def __init__(self, r=[0,0,0], no = [0,0,1], size=[1,1,0.1], id=""): 
        self.r = r
        self.no = no
        self.size = size
        self.id = id
    
    def slit(self,x,y):
        pass

class Filter_freq:
    def __init__(self, r=[0,0,0], no = [0,0,1], size=[1,1,1], id="", filename=""): 
        self.r = r
        self.no = no
        self.size = size
        self.id = id
        self.filename = filename

class OptDrift:
    '''
    Free space
    '''
    def __init__(self, r=[0,0,0], size=[0,0,0], id=""):
        self.r = r
        self.size = size
        self.id = id

class Detector:
    def __init__(self, r=[0,0,0], size=[1.,1.,0], no=[0,0,1], nx=101, ny=101, id=""):
        self.r = r
        self.size = size
        self.nx = nx
        self.ny = ny
        self.no = no
        self.matrix = np.zeros([self.nx,self.ny], dtype=np.double)
        self.id = id

    def hit(self, r):
        i1 = (self.size[0] + r[0])*self.nx / (2*self.size[0])
        i2 = (self.size[1] + r[1])*self.nx / (2*self.size[1])
        if (0 <= i1) and (i1 < self.nx) and (0 <= i2) and (i2 < self.ny):
            self.matrix[i1, i2] += 1
    
    def clear(self):
        for i1 in xrange(self.nx):
            for i2 in xrange(self.ny):
                self.matrix[i1, i2] = 0
                
                
class Geometry():
    def __init__(self, geo=[]):
        self.geo = []
        
        for o in geo:
            
            inserted = False
            for i in range(len(self.geo)):
                if o.r[2] < self.geo[i].r[2]:
                    self.geo.insert(i, o)
                    inserted = True
                    break
 
            if not inserted: self.geo.append(o)
                    
        self.create_lookup_tables()

    def __call__(self):
        for obj in self.geo: yield obj

    def __getitem__(self, idx):
        return self.geo[idx]
    
    def __setitem__(self, idx, val):
        self.geo[idx] = val
        
    def find(self, id):
        for obj in self.geo:
            if obj.id == id:
                return obj
        return None
        
    def create_lookup_tables(self):
        self.c_table = []
        
        c = 1.0   # speed of light
        
        for o in self():
            z_min = o.r[2] - o.size[2]
            z_max = o.r[2] + o.size[2]
            c_fun = lambda r: c
                        
            if o.__class__ == Lense:
                c_fun = o.c_fun 
                            
            # append free space between elements
            if len(self.c_table)>0 and z_min > self.c_table[-1][1]:
                self.c_table.append([self.c_table[-1][1], z_min, lambda r: c])
            
            self.c_table.append([z_min, z_max, c_fun])
            
        for o in self.c_table:
            print('lookup table', o)
            
    def get_c(self, r):
        for i in range( len(self.c_table) ):
            
            if r[2] >= self.c_table[i][0] and r[2] <= self.c_table[i][1]:
                #print r[2], self.c_table[i], self.c_table[i][2](r)
                return  self.c_table[i][2](r)
            
        print('warning: outside of geometry. Defaulting to c=1 ', r)
        return 1.0
    
    
# optical properties
    
class Material:
    def __init__(self, n=1.0):
        self.n = n
        
        
MaterialDb = {}

MaterialDb["valuum"] = Material(n=1.0)
MaterialDb["air"] = Material(n=1.000277)
MaterialDb["water"] = Material(n=1.333)
MaterialDb["diamond"] = Material(n=2.419)

# glass supposed to have 1.50 - 1.90
MaterialDb["glass"] = Material(n=1.50)

MaterialDb["silicon"] = Material(n=3.96)



