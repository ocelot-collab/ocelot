from xframework.cpbd.optics import *

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


class Constraints:
    def __init__(self):
        self.twiss_dict = {}
        self.map_dict = {}

    def __getattr__(self, name):
        if name not in self.__dict__.keys():
            return 0
        else:
            return self.__dict__[name]



def match(lat, constr, params):
    return lat

def test():
    constr = Constraints()
    
    # conatrsint type 1
    constr.twiss_dict['qf1'] = Twiss(beta_x = 15., beta_y = 12)
    constr.twiss_dict['qf1'] = Twiss(beta_x = 15., alpha_x = 0)
    
    # conatrsint type 2
    constr.map_dict['qf1','qf2'] = {'r45':0, 'r12':1}
    
    params = [lat['qf.k1'], lat['qd.k1']]
    
if __name__ == "__main__":
    test()
    

#lat2 = match(lat, constr, params)

#print lat['qf.k1'], lat['qd.k1']
    