'''
ray optics
'''

from numpy import *
from elements import *

import numpy as np


intersection_tol = 1.e-6

class Ray(object):
    def __init__(self,r0=[0,0,0], k=[0,0,1], lamb = 2.0):
        self.r0 = [np.array(r0)]
        self.k = [np.array(k)]
        self.lamb = lamb
        self.s = [1]
        self.c = 3.e8
        self.obj = [Drift()]
        
    @property
    def w(self):
        """I'm the 'x' property."""
        print "getter of w called"
        return (2.*pi * self.c) / self.lamb

    @w.setter
    def w(self, value):
        print "setter of w called" 
        self.lamb = (2.*pi*self.c) / self.value

def find_intersections(ray, geo):
    """
    find the first intersection point of a ray with geometry
    """
    s = np.inf
    obj = None
    r_loc = np.array([0,0])
    for o in geo():
        print 'determining intersection with', o, o.r, o.no
        nk = np.dot(o.no, ray.k[-1])
        nr = np.dot(o.no, o.r - ray.r0[-1])
        #print nr, nk
        if nr*nk > 0:
            #TODO: check that intersection is on aperture
            s_int= nr/nk #nr/nk is path length to intersection along the ray
            if s_int < s and s_int > intersection_tol: 
                
                
                r_int = ray.r0[-1] + s_int * ray.k[-1]
                
                # check intersection with elliptic 'aperture'
                
                r_loc = r_int - o.r
                phi = np.arccos(o.no[2]/ np.linalg.norm(o.no))
                r_loc[1] = r_loc[1] * cos(phi)
                r_loc[2] = r_loc[2] * sin(phi)
                print 'r_loc=', r_loc, 'size=',o.size
                if (r_loc[0]/o.size[0])**2 + (r_loc[1]/o.size[1])**2 <= 1:
                    s = s_int
                    obj = o
                    print 'fits aperture'
                else:
                    print 'fits aperture not'
                
    
    return s, obj, r_loc

def trace(ray, geo):
    """
    tracing the ray, starting from last segment
    """
    n_reflect = 0
    n_reflect_max = 4
    
    
    while n_reflect < n_reflect_max:
        
        debug('current ray at: ', ray.r0[-1])
        
        # ray length to intersection
        s, obj, r_loc = find_intersections(ray, geo)
        debug('found intersection', s, obj)
        
        if s == np.inf:
            info('ray leaves geometry, terminating')
            break
        
        #propagate to boundary
        ray.s[-1] = s 
        r0_new = ray.r0[-1] + ray.k[-1] * ray.s[-1]
        
        # reflect
        if obj.__class__ == Mirror:
            
            debug('reflecting off', obj.id)
            debug(np.dot(obj.no, ray.k[-1]) / ( np.linalg.norm(obj.no) * np.linalg.norm(ray.k[-1]) ))
            phi = np.arccos( np.dot(obj.no, ray.k[-1]) / ( np.linalg.norm(obj.no) * np.linalg.norm(ray.k[-1]) ) )
            
            
            debug('ray/normal angle', phi / pi ,'pi')
            
            
            sgn = np.dot([1,0,0],np.cross(obj.no, ray.k[-1]))
            
            phi = (2*phi - pi) * sgn
            
            
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
        


