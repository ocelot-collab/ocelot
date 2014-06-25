'''
user interface for viewing/editing electron optics layouts
'''

import sys, os, csv


import matplotlib
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
import matplotlib.path as mpath

'''
try:
    from PyQt4.QtCore import *
    from PyQt4.QtGui import *
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
except:
    print 'WARNING: Qt not installed, some graphics may not work properly'
'''

from ocelot.cpbd.optics import *
import numpy as np


def plot_lattice(lat, axis, alpha=1.0, params={'kmax':2.0, 'ang_max':0.5e-2}):
    axis.grid(True)
    

    pos = 0.0
    offs = np.array([0,0.0])
     
    ang_max = params['ang_max']           # max dipole strength in lattice
    min_dipole_height = 0.1  # dipole of any strength will show at least this strength 

    kmax = params['kmax']
    min_quad_height = 0.1


    rendered_seq = []
    rendered_len = 0
    total_len = 0
    
    
    for e in lat.sequence:
        
        if e.type in ['bend','sbend', 'rbend', 'quadrupole', 'undulator', 'drift', 'monitor','hcor','vcor', 'cavity','edge']:
            e.s = total_len
            rendered_seq.append(e)
            rendered_len += e.l
        total_len += e.l
    
    for e in rendered_seq:
        dx = e.l 
        
        if e.type in ['bend', 'sbend', 'rbend','hcor','vcor']:
            axis.add_patch( mpatches.Rectangle(offs+np.array([pos, 0.0]), dx, 
                                               np.sign(-e.angle) * min_dipole_height - e.angle / ang_max * (1- min_dipole_height), 
                                               color='#0099FF', alpha = alpha))

        if e.type in ['quadrupole']:
                        
            if e.k1 >= 0:
                axis.add_patch( mpatches.Ellipse(offs+np.array([pos, 0.0]), dx, min_quad_height + abs(e.k1/kmax)*2, color='green', alpha=alpha) )
            else:
                Path = mpath.Path
                h = abs(e.k1/kmax) + min_quad_height / 2

                verts = np.array([
                                  (dx, h),
                                  (-dx, h),
                                  (-dx/4, 0),
                                  (-dx, -h),
                                  (dx, -h),
                                  (dx/4, 0),
                                  (dx, h)
                                  ])
                
                codes = [Path.MOVETO, Path.LINETO, Path.CURVE3, Path.LINETO, Path.LINETO, Path.CURVE3, Path.CURVE3]
                
                
                path = mpath.Path(verts+offs+np.array([pos, 0.0]), codes)
                patch = mpatches.PathPatch(path, color='green', alpha=alpha)
                axis.add_patch(patch)

            
        
        if e.type in ['undulator']:
            nper = 16
            dxs = dx / nper / 2.0
            
            height = 0.7
            gap = 0.1
            if e.Kx**2 + e.Ky**2 < 1.e-6:
                gap = 0.2
            
            for iseg in np.arange(0,nper,2):                                
                axis.add_patch( mpatches.Rectangle(offs+np.array([pos + iseg * dxs, gap]), dxs, height, color='red', alpha = alpha))
                axis.add_patch( mpatches.Rectangle(offs+np.array([pos + iseg * dxs, -height-gap]), dxs, height, color='blue', alpha = alpha) )
                axis.add_patch( mpatches.Rectangle(offs+np.array([pos + (iseg+1) * dxs, gap]), dxs, height, color='blue', alpha = alpha) )
                axis.add_patch( mpatches.Rectangle(offs+np.array([pos + (iseg+1) * dxs, -height-gap]), dxs, height, color='red', alpha = alpha  ) )

        if e.type in ['cavity']:
            nper = 16
            dxs = dx / nper / 2.0
            
            height = 0.7
            gap = 0.1
            axis.add_patch( mpatches.Ellipse(offs+np.array([pos + dx/2, 0.0]), 
                            dx/3, 0.5, 
                            color='#FF0033', alpha = alpha))

            axis.add_patch( mpatches.Ellipse(offs+np.array([pos + dx/6, 0.0]), 
                            dx/3, 0.5, 
                            color='#FF0033', alpha = alpha))

            axis.add_patch( mpatches.Ellipse(offs+np.array([pos + 5*dx/6, 0.0]), 
                            dx/3, 0.5, 
                            color='#FF0033', alpha = alpha))

     
        if e.type in ['monitor']:
            Path = mpath.Path
            h = 0.005
            offs_mon = np.array([0.0, 0.03])

            verts = np.array([
                              (0, h),
                              (0, -h),
                              (h, 0),
                              (-h, 0)
                              ])
            
            codes = [Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO]
            
            
            verts += offs_mon
            
            path = mpath.Path(verts+offs+np.array([pos, 0.0]), codes)
            patch = mpatches.PathPatch(path, color='black', lw=2, alpha = alpha*0.5)
            axis.add_patch(patch)
            axis.add_patch( mpatches.Circle(offs+offs_mon+np.array([pos, 0.0]), h/2, color='black', alpha = alpha*0.5))

        pos += dx
    
    
    lw = 0.005
    axis.add_patch( mpatches.Rectangle(offs - (0,lw/2), 1.0, lw, color='black', alpha = alpha*0.4) )

 
    axis.set_ylim(-1,1)
    axis.set_xticks([])
    axis.set_yticks([])
    