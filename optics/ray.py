'''
ray optics
'''

from numpy import *
from ocelot.optics.elements import *

import numpy as np


intersection_tol = 1.e-6

class Ray(object):
    def __init__(self,r0=[0,0,0], k=[0,0,1], lamb = 2.0):
        self.r0 = [np.array(r0)]
        self.k = [np.array(k)]
        self.lamb = lamb
        self.s = [1]
        self.c = 3.e8
        self.obj = [OptDrift()]
        
    @property
    def w(self):
        """I'm the 'x' property."""
        print("getter of w called")
        return (2.*pi * self.c) / self.lamb

    @w.setter
    def w(self, value):
        print("setter of w called" )
        self.lamb = (2.*pi*self.c) / self.value

def find_intersections(ray, geo):
    """
    find the first intersection point of a ray with geometry
    """
    s = np.inf
    obj = None
    r_loc = np.array([0,0])
    no = None
    for o in geo():
        debug('checking intersection:', o.id, o.r, o.no)
        nk = np.dot(o.no, ray.k[-1])
        nr = np.dot(o.no, o.r - ray.r0[-1])
                
        #print nr, nk
        if nr*nk > 0:
            #TODO: check that intersection is on aperture
            s_int= nr/nk #nr/nk is path length to intersection along the ray
            if s_int < s and s_int > intersection_tol: 
                no = o.no
                
                r_int = ray.r0[-1] + s_int * ray.k[-1]
                debug('r_int=', r_int)                                
                # check intersection with elliptic 'aperture'                
                r_loc = r_int - o.r
                debug('r_loc unrotated=', r_loc)
                
                
                phi = np.arccos(o.no[2]/ np.linalg.norm(o.no))
                r_loc[1] = r_loc[1] * cos(phi) + r_loc[2] * sin(phi) 
                r_loc[2] = r_loc[2] * cos(phi) - r_loc[1] * sin(phi)
                                
                debug('r_loc=', r_loc, 'size=',o.size)
                
                # correct intersection for curved elements
                if o.__class__ == EllipticMirror:
                    # note that a[0] is the major axis
                    
                    #r_loc[0] = r_loc[0] * cos(o.roll) + r_loc[1] * sin(o.roll) 
                    #r_loc[1] = r_loc[1] * cos(o.roll) - r_loc[0] * sin(o.roll)

                    debug('r_loc=', r_loc, 'size=',o.size)

                    kz = ray.k[-1][2]                    
                    ky = ray.k[-1][1]

                    rz = r_int[2]
                    ry = r_int[1] - o.a[1]
                    
                    
                    az = o.a[0]
                    ay = o.a[1]
                    
                    
                    #debug('angle=', np.arctan(ry/rz) / pi)
                    
                    a_ = kz**2/ az**2 + ky**2/ ay**2
                    b_ = -2*(kz*rz / az**2 + ky*ry / ay**2)
                    c_ = rz **2/ az**2 + ry**2/ ay**2 - 1.
                    d_ = b_**2 - 4*a_*c_
               
                    s1 = (- b_ + np.sqrt(d_) ) / (2.*a_)
                    s2 = (- b_ - np.sqrt(d_) ) / (2.*a_)
                    
                    s_cor = np.min([s1,s2])
                    
                    #debug('D=', d_, 's12=',s1,s2, s_cor)
                    #debug( (rz - s_cor*kz)**2 / az**2 + (ry - s_cor*ky)**2 / ay**2 )
                    #debug( (rz )**2 / az**2 + (ry )**2 / ay**2 )
                    
                    debug('s_old=', s_int)
                    s_int = s_int - s_cor
                    debug('s_new=', s_int)
                    
                    r_int = r_int - s_cor * ray.k[-1]
                    
                    r_loc = r_int - o.r
                    #r_loc[1] = r_loc[1] * cos(phi)
                    #r_loc[2] = r_loc[2] * sin(phi)
                    debug('r_loc_new=', r_int, r_loc)
     
                    ang = arctan2(1./az*r_loc[2], 1./ay*(-r_loc[1] + ay))
                    
                    #debug(r_loc[2], r_loc[1] - ay)
                    debug('ellipse angle=', ang)                    
                    debug('local coord:', az*sin(ang), -ay*cos(ang) + ay)
                    
                    no = np.array([0, cos(ang),-ay/az*sin(ang)]) / np.sqrt(ay**2/az**2*sin(ang)**2 + cos(ang)**2 ) 
            
                    debug('no=',no)
                    debug(o.no)
     
                
                
                if (r_loc[0]/o.size[0])**2 + (r_loc[1]/o.size[1])**2 <= 1:
                    s = s_int
                    obj = o
                    debug('fits aperture')
                else:
                    debug('fits aperture not')
                
    
    return s, obj, r_loc, no


def refl_matrix(no):
    x, y, z = no
    M = np.matrix([[-1. + 2.*x**2, 2.*x*y, 2.*x*z],
                   [2.*y*x, -1. + y**2*(1.+1), y*z*(1.+1)],
                   [2.*z*x, 2.*z*y , -1. + 2.*z**2]])
    return M


def trace(ray, geo):
    """
    tracing the ray, starting from last segment
    """
    n_reflect = 0
    n_reflect_max = 4
    
    
    while n_reflect < n_reflect_max:
        
        debug('ray at: ', ray.r0[-1])
        
        # ray length to intersection
        s, obj, r_loc, no = find_intersections(ray, geo)
        
        if s == np.inf:
            info('ray leaves geometry, terminating')
            break

        
        debug('intersection: s=', s, 'obj:',  obj.id, 'normal',  no)
                
        #propagate to boundary
        ray.s[-1] = s 
        r0_new = ray.r0[-1] + ray.k[-1] * ray.s[-1]
        
        # reflect
        if obj.__class__ == Mirror:
            
            debug('reflecting off', obj.id)

            k_new = np.asarray(np.dot( refl_matrix(obj.no), -ray.k[-1]))[0]
            
            debug(ray.k[-1], '--->', k_new)
            s_new = 1
            ray.r0.append(r0_new)
            ray.k.append(k_new)
            ray.s.append(s_new)
            
            n_reflect += 1

        elif obj.__class__ == EllipticMirror:
            
            debug('reflecting off', obj.id)                        
            debug('no=',no,'k=',ray.k[-1])
            
            '''
            cs = np.dot(no, ray.k[-1]) / ( np.linalg.norm(no) * np.linalg.norm(ray.k[-1]) )
            debug('cos=',cs)
            
            if np.abs(cs) > 1:
                print 'warning, reflection angle adjustment by ', cs + 1.0
                if cs > 1: cs = 1.0
                else: cs = -1.0
            
            phi = np.arccos( cs )
            
            
            debug('ray/normal angle', phi / pi ,'pi')
            
            
            sgn = np.dot([1,0,0],np.cross(no, ray.k[-1]))
            if np.linalg.norm(sgn) < 1.e-9: 
                sgn = sgn / np.linalg.norm(sgn)
            else:
                sgn = 1.0
            
            debug('sgn=',sgn)
            
            phi = (2*phi - pi) * sgn
            
            
            debug('e:rotating by:', phi / pi, 'pi')
            
        
            M = np.matrix([[1, 0, 0],
                           [0, cos(phi), sin(phi)],
                           [0, -sin(phi), cos(phi)]])

            k_new = np.asarray(np.dot(M, ray.k[-1]))[0]
            '''
            
            k_new = np.asarray(np.dot( refl_matrix(no), -ray.k[-1]))[0]
                        
            debug(ray.k[-1], '--->', k_new)
            s_new = 1
            ray.r0.append(r0_new)
            ray.k.append(k_new)
            ray.s.append(s_new)
            
            n_reflect += 1
            
        elif obj.__class__ == Grating:
            
            debug('reflecting off', obj.id)
            debug(np.dot(obj.no, ray.k[-1]) / ( np.linalg.norm(obj.no) * np.linalg.norm(ray.k[-1]) ))
            phi = np.arccos( np.dot(obj.no, ray.k[-1]) / ( np.linalg.norm(obj.no) * np.linalg.norm(ray.k[-1]) ) )
            
            
            debug('ray/normal angle', phi / pi ,'pi')
            
            
            sgn = np.dot([1,0,0],np.cross(obj.no, ray.k[-1]))
            
            phi = (2*phi - pi) * sgn * (1+ 0.1 * ray.lamb)
            
            
            debug('rotating by:', phi / pi, 'pi')
            
        
            M = np.matrix([[1, 0, 0],
                           [0, cos(phi), sin(phi)],
                           [0, -sin(phi), cos(phi)]])

            k_new = np.asarray(np.dot(M, ray.k[-1]))[0]
            
            #print '###',ray.k[-1].shape, '###',obj.no
            #k_new = rotate_pi(ray.k[-1], obj.no )
            
            #k_new - 
            debug('k_new--->',k_new)
            s_new = 1
            ray.r0.append(r0_new)
            ray.k.append(k_new)
            ray.s.append(s_new)
            
            n_reflect += 1
            
        elif obj.__class__ == Aperture:
            if (r_loc[0] / obj.d[0])**2 + (r_loc[1] / obj.d[1])**2 > 1:
                debug('ray stopped at aperture')
                break
            else:
                r0_new = r0_new  +  ray.k[-1]*intersection_tol * 2
                k_new = ray.k[-1]
                s_new = ray.s[-1]
                n_reflect += 1
                ray.r0.append(r0_new)
                ray.k.append(k_new)
                ray.s.append(s_new)
        
        elif obj.__class__ == Crystal:
            r0_new = r0_new  +  ray.k[-1]*intersection_tol * 2
            k_new = ray.k[-1]
            s_new = ray.s[-1]
            n_reflect += 1
            ray.r0.append(r0_new)
            ray.k.append(k_new)
            ray.s.append(s_new)
        elif obj.__class__ == Lense:
            debug('tracing thru lense, f=',obj.f, ' [m]')
            r0_new = r0_new  +  ray.k[-1]*intersection_tol * 2
            k_new = np.array([ray.k[-1][0]- r_loc[0]*ray.k[-1][2] / obj.f, ray.k[-1][1] - r_loc[1]*ray.k[-1][2] / obj.f, ray.k[-1][2] ])
            s_new = ray.s[-1]
            n_reflect += 1
            ray.r0.append(r0_new)
            ray.k.append(k_new)
            ray.s.append(s_new)
        elif obj.__class__ == Detector:
            debug('detector hit')
            obj.hit(r_loc)
            k_new = ray.k[-1]
            s_new = ray.s[-1]
            n_reflect += 1
            ray.r0.append(r0_new)
            ray.k.append(k_new)
            ray.s.append(s_new)                
        else:
            warn('no propagator available, optics element:', obj)
            r0_new = r0_new  +  ray.k[-1]*intersection_tol * 2
            k_new = ray.k[-1]
            s_new = ray.s[-1]
            n_reflect += 1
            ray.r0.append(r0_new)
            ray.k.append(k_new)
            ray.s.append(s_new)

        ray.obj.append(obj)
        


