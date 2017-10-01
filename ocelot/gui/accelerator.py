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
import matplotlib.pyplot as plt
from ocelot.cpbd.optics import *
import numpy as np
from copy import deepcopy

import matplotlib.font_manager as font_manager


def plot_lattice(lat, axis, alpha=1.0, params={'kmax':2.0, 'ang_max':0.5e-2}, s_start = 0.0):
    axis.grid(True)
    

    pos = s_start
    offs = np.array([0,0.0])
     
    ang_max = params['ang_max']           # max dipole strength in lattice
    min_dipole_height = 0.1  # dipole of any strength will show at least this strength 

    min_solenoid_height = 0.1
    sol_max = 0.1

    kmax = params['kmax']
    min_quad_height = 0.1


    rendered_seq = []
    rendered_len = 0.0
    total_len = 0.0
    
    
    for e in lat.sequence:
        
        #if e.type in ['bend','sbend', 'rbend', 'quadrupole', 'undulator', 'drift', 'monitor','hcor','vcor', 'cavity','edge', 'solenoid']:
        if e.__class__ in [Bend, SBend, RBend, Quadrupole, Undulator, Drift, Monitor, Hcor, Vcor, Cavity, Edge, Solenoid]:
            e.s = total_len
            rendered_seq.append(e)
            rendered_len += e.l
        total_len += e.l
    
    for e in rendered_seq:
        dx = e.l 
        
        if e.__class__ in  [Bend, SBend, RBend,Hcor, Vcor]:
            axis.add_patch( mpatches.Rectangle(offs+np.array([pos, 0.0]), dx, 
                                               np.sign(-e.angle) * min_dipole_height - e.angle / ang_max * (1- min_dipole_height), 
                                               color='#0099FF', alpha = alpha))

        if e.__class__ in [Solenoid]:
            axis.add_patch( mpatches.Rectangle(offs+np.array([pos, 0.0]), dx, 
                                               np.sign(-e.k) * min_solenoid_height - e.k / sol_max * (1- min_solenoid_height), 
                                               color='#FF99FF', alpha = alpha))


        if e.__class__ in [Quadrupole]:
                        
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

            
        
        if e.__class__ in [Undulator]:
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

        if e.__class__ in [Cavity]:
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

     
        if e.__class__ in [Monitor]:
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

 
    #axis.set_ylim(-1,1)
    #axis.set_xticks([])
    #axis.set_yticks([])

# Tomin Sergey functions


def elem_cord(lat):
    quad = np.array([[0,0]])
    bend = np.array([[0,0]])
    sext = np.array([[0,0]])
    multi = np.array([[0,0]])
    c = []
    corr = np.array([[0,0]])
    mons = np.array([[0,0]])
    cav = np.array([[0,0]])
    mat = np.array([[0,0]])
    und = np.array([[0,0]])
    drft = np.array([[0,0]])
    L = 0
    #for elem in lat.sequence:
    #if elem.type == "drift" and elem.l == 0:
    #lat.sequence.remove(elem)
    
    for elem in lat.sequence:
        dL = 0.
        if elem.l == 0:
            dL = 0.03
        temp = np.array([[L - dL, 0]])
        temp = np.append(temp, [[L - dL, 1]], axis = 0)
        temp = np.append(temp, [[L+elem.l + dL, 1]], axis = 0)
        temp = np.append(temp, [[L+elem.l + dL, 0]], axis = 0)
        #temp = temp.reshape((4,2))
        
        if elem.__class__ == Quadrupole:
            k1 = elem.k1
            quad = np.append(quad, [[L, 0]], axis=0)
            quad = np.append(quad, [[L, k1]], axis=0)
            #quad = np.append(quad, [[L+elem.l/2, k1*(1 + 0.2 * np.sign(elem.k1)) ]], axis=0)
            quad = np.append(quad, [[L+elem.l, k1]],axis=0)
            quad = np.append(quad, [[L+elem.l, 0]],axis=0)

        elif elem.__class__ == Cavity:
            k1 = 1.
            cav = np.append(cav, [[L, 0]], axis=0)
            cav = np.append(cav, [[L, k1]], axis=0)
            cav = np.append(cav, [[L+elem.l, k1]],axis=0)
            cav = np.append(cav, [[L+elem.l, 0]],axis=0)

        elif elem.__class__  == Drift:
            k1 = 1.
            drft = np.append(drft, [[L, 0]], axis=0)
            drft = np.append(drft, [[L, 0]], axis=0)
            drft = np.append(drft, [[L+elem.l, 0]],axis=0)
            drft = np.append(drft, [[L+elem.l, 0]],axis=0)

        elif elem.__class__ == Matrix:
            #k1 = 1.
            mat = np.append(mat, [[L, 0]], axis=0)
            #mat = np.append(mat, [[L, k1]], axis=0)
            #mat = np.append(mat, [[L+elem.l, k1]],axis=0)
            mat = np.append(mat, [[L+elem.l, 0]],axis=0)

        elif elem.__class__ == Undulator:
            k1 = 1.
            und = np.append(und, [[L, 0]], axis=0)
            und = np.append(und, [[L, k1]], axis=0)
            und = np.append(und, [[L+elem.l, k1]],axis=0)
            und = np.append(und, [[L+elem.l, 0]],axis=0)

        elif elem.__class__ in [SBend, RBend, Bend]:
            if elem.l == 0:
                h = 0
            else:
                h = elem.angle/elem.l
            temp[:,1] = temp[:,1]*h
            bend = np.append(bend, temp, axis = 0)
        
        elif elem.__class__ == Sextupole:

            temp[:,1] = temp[:,1]*elem.k2
            sext = np.append(sext, temp, axis = 0)

        elif elem.__class__ == Multipole:
            if sum(abs(elem.kn)) != 0:
                temp[:,1] = temp[:,1]*sum(elem.kn)/sum(abs(elem.kn))
            else:
                temp[:,1] = temp[:,1]*0.
            multi = np.append(multi, temp, axis=0)
        
        elif elem.__class__ in [Hcor, Vcor]:
            temp[:,1] = temp[:,1]#*abs(elem.angle)
            corr = np.append(corr, temp, axis=0)
        
        elif elem.__class__ in [Monitor]:
            temp[:,1] = temp[:,1]#*abs(elem.angle)
            mons = np.append(mons, temp, axis=0)
        #c.append((L,  elem.l+0.03))
        L += elem.l
    if len(quad) != 1:
        quad[:,1] = quad[:,1]/max(quad[:,1])
    if len(bend) != 1:
        if max(bend[:,1]) == 0:
            bend[:,1] = 0
        else:
            bend[:,1] = bend[:,1]/max(bend[:,1])
    if len(sext) != 1 and max(sext[:,1] != 0):
        sext[:,1] = sext[:,1]/max(np.abs(sext[:,1]))
    if len(corr) != 1 and max(corr[:,1] != 0):
        corr[:,1] = corr[:,1]/max(corr[:,1])
    #if len(corr) != 1 and max(mons[:,1] != 0):
    #mons[:,1] = mons[:,1]/max(mons[:,1])
    return quad, bend, sext, corr, mons, cav, mat, und, multi, drft


dict_plot = {Quadrupole: {"scale": 0.7, "color": "r",            "edgecolor": "r",          "label": "quad"},
             Sextupole:  {"scale": 0.5, "color": "g",            "edgecolor": "g",          "label": "sext"},
             Octupole:   {"scale": 0.5, "color": "g",            "edgecolor": "g",          "label": "sext"},
             Cavity:     {"scale": 0.7, "color": "orange",       "edgecolor": "lightgreen", "label": "cav"},
             Bend:       {"scale": 0.7, "color": "lightskyblue", "edgecolor": "k",          "label": "bend"},
             RBend:      {"scale": 0.7, "color": "lightskyblue", "edgecolor": "k",          "label": "bend"},
             SBend:      {"scale": 0.7, "color": "lightskyblue", "edgecolor": "k",          "label": "bend"},
             Matrix:     {"scale": 0.7, "color": "pink",         "edgecolor": "k",          "label": "mat"},
             Multipole:  {"scale": 0.7, "color": "g",            "edgecolor": "k",          "label": "mult"},
             Undulator:  {"scale": 0.7, "color": "pink",         "edgecolor": "k",          "label": "und"},
             Monitor:    {"scale": 0.5, "color": "orange",       "edgecolor": "orange",     "label": "mon"},
             Hcor:       {"scale": 0.7, "color": "c",            "edgecolor": "c",          "label": "cor"},
             Vcor:       {"scale": 0.7, "color": "c",            "edgecolor": "c",          "label": "cor"},
             Drift:      {"scale": 0.,  "color": "k",            "edgecolor": "k",          "label": ""},
             Marker:     {"scale": 0.,  "color": "k",            "edgecolor": "k",          "label": "mark"},
             Edge:       {"scale": 0.,  "color": "k",            "edgecolor": "k",          "label": ""},
             Solenoid:   {"scale": 0.7, "color": "g",            "edgecolor": "g",          "label": "sol"},
             UnknownElement:{"scale": 0.7, "color": "g",            "edgecolor": "g",          "label": "unk"},
             }


def new_plot_elems(fig, ax, lat, s_point = 0, nturns = 1, y_lim = None,y_scale = 1, legend = True):
    dict_copy=deepcopy(dict_plot)
    alpha = 1
    ax.set_ylim((-1,1.5))
    if y_lim != None:
        ax.set_ylim(y_lim)
    points_with_annotation = []
    L = 0.
    q = []
    b = []
    c = []
    s = []
    u = []
    rf = []
    m = []
    for elem in lat.sequence:
        if elem.__class__ == Quadrupole:
            q.append(elem.k1)
        elif elem.__class__ in [Bend, RBend, SBend]:
            b.append(elem.angle)
        elif elem.__class__ in [Hcor, Vcor]:
            c.append(elem.angle)
        elif elem.__class__ == Sextupole:
            s.append(elem.k2)
        elif elem.__class__ == Undulator:
            u.append(elem.Kx + elem.Ky)
        elif elem.__class__ == Cavity:
            rf.append(elem.v )
        elif elem.__class__ == Multipole:
            m.append(sum(np.abs(elem.kn)))
    q_max = np.max(np.abs(q))if len(q) !=0 else 0
    b_max = np.max(np.abs(b))if len(b) !=0 else 0
    s_max = np.max(np.abs(s))if len(s) !=0 else 0
    c_max = np.max(np.abs(c))if len(c) !=0 else 0
    u_max = np.max(np.abs(u))if len(u) !=0 else 0
    rf_max = np.max(np.abs(rf))if len(rf) !=0 else 0
    m_max = np.max(m) if len(m) !=0 else 0
    ncols = np.sign(len(q)) + np.sign(len(b)) + np.sign(len(s)) + np.sign(len(c)) + np.sign(len(u)) + np.sign(len(rf))+ np.sign(len(m))

    labels_dict = {}
    for elem in dict_copy.keys():
        labels_dict[elem] = dict_copy[elem]["label"]
    for elem in lat.sequence:
        if elem.__class__ in [Marker, Edge]:
            L +=elem.l
            continue
        l = elem.l
        if l == 0:
            l = 0.03
        #type = elem.type
        scale = dict_copy[elem.__class__]["scale"]
        color = dict_copy[elem.__class__]["color"]
        label = dict_copy[elem.__class__]["label"]
        ecolor = dict_copy[elem.__class__]["edgecolor"]
        ampl = 1
        s_coord = np.array([L + elem.l/2. - l/2., L + elem.l/2. - l/2., L + elem.l/2. +l/2., L + elem.l/2. +l/2., L + elem.l/2.- l/2.]) + s_point
        if elem.__class__ == Quadrupole:
            ampl = elem.k1/q_max if q_max != 0 else 1
            point, = ax.fill(s_coord,  (np.array([-1, 1, 1, -1, -1])+1)*ampl*scale*y_scale, color, edgecolor=ecolor,
                             alpha = alpha, label=dict_copy[elem.__class__]["label"])
            dict_copy[elem.__class__]["label"] = ""

        elif elem.__class__ in [Bend, RBend, SBend]:
            ampl = elem.angle/b_max if b_max != 0 else 1
            point, = ax.fill(s_coord, (np.array([-1, 1, 1, -1, -1])+1)*ampl*scale*y_scale, color,
                             alpha = alpha, label=dict_copy[elem.__class__]["label"])
            dict_copy[elem.__class__]["label"] = ""

        elif elem.__class__ in [Hcor, Vcor]:

            ampl = elem.angle/c_max if c_max != 0 else 0.5
            #print c_max, elem.angle, ampl
            if elem.angle == 0:
                ampl=0.5
                point, = ax.fill(s_coord, (np.array([-1, 1, 1, -1, -1]))*ampl*scale*y_scale, "lightcyan",  edgecolor="k",
                             alpha = 0.5, label=dict_copy[elem.__class__]["label"])
            else:
                point, = ax.fill(s_coord, (np.array([-1, 1, 1, -1, -1])+1)*ampl*scale*y_scale, color,  edgecolor=ecolor,
                             alpha = alpha, label=dict_copy[elem.__class__]["label"])
            dict_copy[Hcor]["label"] = ""
            dict_copy[Vcor]["label"] = ""

        elif elem.__class__ == Sextupole:
            ampl = (elem.k2)/s_max if s_max != 0 else 1
            point, = ax.fill(s_coord, (np.array([-1, 1, 1, -1, -1])+1)*ampl*scale*y_scale, color,
                             alpha = alpha, label=dict_copy[elem.__class__]["label"])
            dict_copy[elem.__class__]["label"] = ""

        elif elem.__class__ == Cavity:
            ampl = 1 # elem.v/rf_max if rf_max != 0 else 0.5
            point, = ax.fill(s_coord, np.array([-1, 1, 1, -1, -1])*ampl*scale*y_scale, color,
                             alpha = alpha, edgecolor = "lightgreen", label=dict_copy[elem.__class__]["label"])
            dict_copy[elem.__class__]["label"] = ""

        elif elem.__class__ == Undulator:
            ampl = elem.Kx/u_max if u_max != 0 else 0.5
            point, = ax.fill(s_coord, np.array([-1, 1, 1, -1, -1])*ampl*scale*y_scale, color,
                             alpha = alpha, label=dict_copy[elem.__class__]["label"])
            dict_copy[elem.__class__]["label"] = ""

        elif elem.__class__ == Multipole:
            ampl = sum(elem.kn)/m_max if u_max != 0 else 0.5
            point, = ax.fill(s_coord, np.array([-1, 1, 1, -1, -1])*ampl*scale*y_scale, color,
                             alpha = alpha, label=dict_copy[elem.__class__]["label"])
            dict_copy[elem.__class__]["label"] = ""

        else:
            point, = ax.fill(s_coord, np.array([-1, 1, 1, -1, -1])*ampl*scale*y_scale, color, edgecolor=ecolor,
                             alpha = alpha)
        annotation = ax.annotate(elem.__class__.__name__+": " + elem.id,
            xy=(L+l/2., 0), #xycoords='data',
            #xytext=(i + 1, i), textcoords='data',
            horizontalalignment="left",
            arrowprops=dict(arrowstyle="simple", connectionstyle="arc3,rad=+0.2"),
            bbox=dict(boxstyle="round", facecolor="w", edgecolor="0.5", alpha=0.9),
                                 fontsize=16
            )
        # by default, disable the annotation visibility
        annotation.set_visible(False)
        L +=elem.l
        points_with_annotation.append([point, annotation])

    def on_move(event):

        visibility_changed = False
        for point, annotation in points_with_annotation:
            should_be_visible = (point.contains(event)[0] == True)

            if should_be_visible != annotation.get_visible():
                visibility_changed = True
                annotation.set_visible(should_be_visible)

        if visibility_changed:
            plt.draw()

    on_move_id = fig.canvas.mpl_connect('motion_notify_event', on_move)
    if legend:
        ax.legend(loc='upper center', ncol=ncols, shadow=False, prop=font_manager.FontProperties(size=15))


def plot_elems(ax, lat, s_point = 0, nturns = 1, y_lim = None,y_scale = 1, legend = True):
    quad, bend, sext, corr, mons, cav, mat, und, multi, drft = elem_cord(lat)
    #print len(quad), len(bend), len(sext), len(corr ),len( mons), len( cav)
    #print cav
    alpha = 1
    ax.set_ylim((-1,1.5))
    if y_lim != None:
        ax.set_ylim(y_lim)
    n = 0
    for i in range(nturns):
        n = 0
        if len(quad)>1:
            ax.fill(quad[:,0]+i*lat.totalLen + s_point, quad[:,1]*y_scale*0.8, "r", alpha = alpha, label = "quad")
            n += 1
        if len(bend)>1:
            ax.fill(bend[:,0]+i*lat.totalLen + s_point, bend[:,1]*y_scale*0.7, "lightskyblue", alpha = alpha, label = "bend")
            n += 1
        if len(sext)>1:
            ax.fill(sext[:,0]+i*lat.totalLen + s_point, sext[:,1]*y_scale*0.8, "green",  edgecolor = "green", alpha = alpha, label = "sext")
            n += 1
        if len(multi)>1:
            ax.fill(multi[:,0]+i*lat.totalLen + s_point, multi[:,1]*y_scale*0.8, "green",  edgecolor = "green", alpha = alpha, label = "multi")
            n += 1
        if len(corr)>1:
            ax.fill(corr[:,0]+i*lat.totalLen + s_point, corr[:,1]*y_scale*0.7, "b", edgecolor = "b", alpha = alpha, label = "corr")
            n += 1
        if len(mons)>1:
            ax.fill(mons[:,0]+i*lat.totalLen + s_point, mons[:,1]*y_scale*0.7, "orange", edgecolor = "orange", alpha = alpha, label = "bpm")
            n += 1
        if len(cav)>1:
            ax.fill(cav[:,0]+i*lat.totalLen + s_point, cav[:,1]*y_scale*0.7, "orange", edgecolor = "lightgreen", alpha = alpha, label = "cav")
            #print cav[:,0]+i*lat.totalLen, cav[:,1]*y_scale*0.7, i*lat.totalLen
            n += 1
        if len(mat)>1:
            ax.fill(mat[:,0]+i*lat.totalLen + s_point, mat[:,1]*y_scale*0.7, "pink", edgecolor = "lightgreen", alpha = alpha, label = "mat")
            n += 1
        if len(und)>1:
            ax.fill(und[:,0]+i*lat.totalLen + s_point, und[:,1]*y_scale*0.7, "pink", edgecolor = "lightgreen", alpha = alpha, label = "und")
            n += 1
        if len(drft)>1:
            ax.fill(drft[:,0]+i*lat.totalLen + s_point, drft[:,1]*y_scale*0.7, "k")
            n += 1
    #ax.broken_barh(s , (y0, yw*1.3), facecolors='green', edgecolor = "green", alpha = alpha, label = "Sext")
    #ax.broken_barh(c , (y0, yw), facecolors='b',edgecolor = "b", alpha = alpha)
    if legend:
        ax.legend(loc='upper center', ncol=n, shadow=False, prop=font_manager.FontProperties(size=15))

def plot_disp(ax,tws, top_plot, font_size):
    S = [p.s for p in tws]#map(lambda p:p.s, tws)
    d_Ftop = []
    Fmin = []
    Fmax = []
    for elem in top_plot:
        #print(elem, tws.__dict__[elem] )
        Ftop = [p.__dict__[elem] for p in tws]
        #for f in Ftop:
        #    print(f)
        #print (max(Ftop))
        Fmin.append(min(Ftop))
        Fmax.append(max(Ftop))
        top_label = r"$"+elem+"$"
        ax.plot(S, Ftop, lw = 2, label=top_label)
        d_Ftop.append( max(Ftop) - min(Ftop))
    d_F = max(d_Ftop)
    if d_F == 0:
        d_Dx = 1
        ax.set_ylim(( min(Fmin)-d_Dx*0.1, max(Fmax)+d_Dx*0.1))
    if top_plot[0] == "E":
        top_ylabel = r"$"+"/".join(top_plot) +"$"+ ", GeV"
    else:
        top_ylabel = r"$"+"/".join(top_plot) +"$"+ ", m"

    yticks = ax.get_yticks()
    yticks = yticks[2::2]
    ax.set_yticks(yticks)
    #for i, label in enumerate(ax.get_yticklabels()):
    #    if i == 0 or i == 1:
    #        label.set_visible(False)
    ax.set_ylabel(top_ylabel)
    
    #ax.plot(S, Dx,'black', lw = 2, label=lable)
    leg2 = ax.legend(loc='upper right', shadow=False, fancybox=True, prop=font_manager.FontProperties(size=font_size))
    leg2.get_frame().set_alpha(0.2)




def plot_betas(ax, S, beta_x, beta_y, font_size):
    ax.set_ylabel(r"$\beta_{x,y}$, m")
    ax.plot(S, beta_x,'b', lw = 2, label=r"$\beta_{x}$")
    ax.plot(S, beta_y,'r', lw = 2, label=r"$\beta_{y}$")
    leg = ax.legend(loc='upper right', shadow=False, fancybox=True, prop=font_manager.FontProperties(size=font_size))
    leg.get_frame().set_alpha(0.2)


def plot_opt_func(lat, tws, top_plot=["Dx"], legend=True, fig_name=None, grid=True, font_size=18):
    """
    function for plotting: lattice (bottom section), vertical and horizontal beta-functions (middle section),
    other parameters (top section) such as "Dx", "Dy", "E", "mux", "muy", "alpha_x", "alpha_y", "gamma_x", "gamma_y"
    lat - MagneticLattice,
    tws - list if Twiss objects,
    top_plot=["Dx"] - parameters which displayed in top section. Example top_plot=["Dx", "Dy", "alpha_x"]
    legend=True - displaying legend of element types in bottom section,
    fig_name=None - name of figure,
    grid=True - grid
    font_size=18 - font size.
    """
    if fig_name == None:
        fig = plt.figure()
    else:
        fig = plt.figure(fig_name)

    plt.rc('axes', grid=grid)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    left, width = 0.1, 0.85
    
    rect1 = [left, 0.63, width, 0.3]
    rect2 = [left, 0.18, width, 0.45]
    rect3 = [left, 0.05, width, 0.13]

    ax_top = fig.add_axes(rect1)
    ax_b = fig.add_axes(rect2, sharex=ax_top)  #left, bottom, width, height
    ax_el = fig.add_axes(rect3, sharex=ax_top)
    for ax in ax_b, ax_el, ax_top:
        if ax!=ax_el:
            for label in ax.get_xticklabels():
                label.set_visible(False)
    
    ax_b.grid(grid)
    ax_top.grid(grid)
    ax_el.set_yticks([])
    ax_el.grid(grid)

    fig.subplots_adjust(hspace=0)
    beta_x = [p.beta_x for p in tws] # list(map(lambda p:p.beta_x, tws))
    beta_y = [p.beta_y for p in tws] #list(map(lambda p:p.beta_y, tws))
    S = [p.s for p in tws] #list(map(lambda p:p.s, tws))
    #plt.plot(S, beta_x)

    plt.xlim(S[0], S[-1])

    plot_disp(ax_top,tws, top_plot, font_size)

    plot_betas(ax_b, S, beta_x, beta_y, font_size)
    #plot_elems(ax_el, lat, s_point = S[0], legend = legend, y_scale=0.8) # plot elements
    new_plot_elems(fig, ax_el, lat, s_point = S[0], legend = legend, y_scale=0.8)



def plot_xy(ax, S, X, Y, font_size):
    ax.set_ylabel(r"$X, Y$, m")
    ax.plot(S, X,'r', lw = 2, label=r"$X$")
    ax.plot(S, Y,'b', lw = 2, label=r"$Y$")
    leg = ax.legend(loc='upper right', shadow=True, fancybox=True, prop=font_manager.FontProperties(size=font_size))
    leg.get_frame().set_alpha(0.5)


def body_trajectory(fig, ax_xy, ax_el, lat, plist):
    X = [p.x for p in plist]
    Y = [p.y for p in plist]
    S = [p.s for p in plist]
    
    font_size = 16
    
    for ax in ax_xy, ax_el:
        if ax!=ax_el:
            for label in ax.get_xticklabels():
                label.set_visible(False)
    
    ax_xy.grid(True)
    ax_el.set_yticks([])
    ax_el.grid(True)
    #plt.xlim(S[0], S[-1])
    
    fig.subplots_adjust(hspace=0)
    
    plot_xy(ax_xy, S, X, Y, font_size)

    plot_elems(ax_el, lat, nturns = 1, legend = False) # plot elements

"""
def plot_current(p_array, charge, num_bins = 200):
    z = p_array.particles[4::6]
    hist, bin_edges = np.histogram(z, bins=num_bins)
    delta_Z = max(z) - min(z)
    delta_z = delta_Z/num_bins
    t_bins = delta_z/speed_of_light
    print "Imax = ", max(hist)*charge/t_bins
    hist = np.append(hist, hist[-1])
    plt.plot(bin_edges, hist*charge/t_bins)
    plt.grid(True)
    plt.title("current")
    plt.xlabel("s, m")
    plt.ylabel("I, A")
    plt.show()
"""

def plot_trajectory(lat, list_particles):
    fig = plt.figure()
    plt.rc('axes', grid=True)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    left, width = 0.1, 0.85
    rect2 = [left, 0.2, width, 0.7]
    rect3 = [left, 0.05, width, 0.15]
    
    ax_xy = fig.add_axes(rect2)  #left, bottom, width, height
    ax_el = fig.add_axes(rect3, sharex=ax_xy)
    
    body_trajectory(fig, ax_xy, ax_el, lat, list_particles)
    plt.show()


def plot_API(lat, legend=True, fig_name=1):
    """
    Function creates a picture with lattice on the bottom part of the picture and top part of the picture can be
    plot arbitrary lines.

    :param lat: MagneticLattice
    :param legend: True, description of the elements, if False it is switched legend off
    :return: fig, ax
    """
    fig = plt.figure(fig_name)
    plt.rc('axes', grid=True)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    left, width = 0.1, 0.85
    rect2 = [left, 0.2, width, 0.7]
    rect3 = [left, 0.05, width, 0.15]

    ax_xy = fig.add_axes(rect2)  #left, bottom, width, height
    ax_el = fig.add_axes(rect3, sharex=ax_xy)

    font_size = 16

    for ax in ax_xy, ax_el:
        if ax!=ax_el:
            for label in ax.get_xticklabels():
                label.set_visible(False)

    ax_xy.grid(True)
    ax_el.set_yticks([])
    ax_el.grid(True)
    #plt.xlim(S[0], S[-1])

    fig.subplots_adjust(hspace=0)

    #plot_xy(ax_xy, S, X, Y, font_size)

    #plot_elems(ax_el, lat, nturns = 1, legend = True) # plot elements
    new_plot_elems(fig, ax_el, lat, nturns = 1, legend = legend)
    return fig, ax_xy


def plot_traj_pulse(lat, list_particles, list_particles2, U1, U2):
    fig = plt.figure()
    plt.rc('axes', grid=True)
    plt.rc('grid', color='0.75', linestyle='-', linewidth=0.5)
    left, width = 0.1, 0.85
    rect2 = [left, 0.2, width, 0.7]
    rect3 = [left, 0.05, width, 0.15]
    
    ax_xy = fig.add_axes(rect2)  #left, bottom, width, height
    ax_el = fig.add_axes(rect3, sharex=ax_xy)
    S = map(lambda p:p.s, list_particles)
    ax_xy.plot(S, U1,'g', lw = 1, label=r"$U1$")
    ax_xy.plot(S, U2,'b', lw = 1, label=r"$U2$")
    #print map(lambda p:p.s, list_particles2)
    ax_xy.plot(map(lambda p:p.s, list_particles2), map(lambda p:p.x, list_particles2),'b', lw = 2, label=r"$p2$")
    body_trajectory(fig, ax_xy, ax_el, lat, list_particles)
    #body_trajectory(fig, ax_xy, ax_el, lat, list_particles2)
    plt.show()



def plot_elem_disp(lat, tws):
    Dx = map(lambda p:p.Dx, tws)
    S = map(lambda p:p.s, tws)
    fig, ax = plt.subplots(1,sharex=True)
    ax.grid(False)
    ax.set_yticks([])
    plt.xlim(tws[0].s, tws[-1].s)
    lim = max(Dx)*1.1
    plot_disp(ax, S, Dx, font_size = 16)
    plot_elems(ax, lat, y_lim = None, y_scale = lim*0.6, legend = False) # plot elements
    ax.set_ylim((0,lim))
    plt.show()


def resonance(Qx, Qy, order = 5):
    ORD = order
    qx1, qy1 = 0,0
    qx2, qy2 = 2,2
    X = []
    Y = []
    Order = []
    params = []
    #Qx = 0.22534
    #Qy = 0.301
    for order in range(1, ORD + 1):
        n = np.arange(-order, order+1)
        m = np.array([order - abs(n), - order + abs(n)]).flatten()
        n = np.tile(n, 2)
        ziped = []
        for n,m in zip(n,m):
            if (n,m) in ziped:
                continue
            ziped.append((n,m))
        #ziped =  zip(n,m)
        #print ziped
        for p in range(-50, 50):
            for n, m in ziped:
                #print p, n,m
                if m != 0:
                    x = [qx1, qx2]
                    y = [(p - n*qx1)/float(m), (p - n*qx2)/float(m)]
                else:
                    x = [p/float(n), p/float(n)]
                    y = [qy1, qy2]
                params.append([n,m,p])
                X.append(x)
                Y.append(y)
                Order.append(order)
    return X,Y,Order,params

def plot_resonance_diag(ax, Qx, Qy, order):
    X,Y,Order,params = resonance(Qx, Qy, order)
    indsort = np.argsort(Order)
    #print Order
    #print len(indsort), len(X)
    X = np.array(X)
    #print X[indsort]
    #fig = plt.figure()
    #ax1 = fig.add_subplot(111)
    """
        def onpick3(event):
        ind = event.ind
        print 'onpick3 scatter:', ind #, np.take(x, ind), np.take(y, ind)
        
        ax1.plot(X[:], Y[:], picker=True)
        plt.xlim(0,1)
        plt.ylim(0,1)
        """
    for i, order in enumerate(Order):
        if order == 1:
            color = "k"
            lw = 3
        elif order == 2:
            color ="r"
            lw = 2
        elif order == 3:
            color = "b"
            lw = 2
        elif order == 4:
            color = "g"
            lw = 0.6
        else:# order == 5:
            color = "c"
            lw = 0.3
        #print array(X[i])+Qx
        #print array([i])+Qy
        plt.plot(array(X[i])+Qx, array(Y[i])+Qy, color, lw = lw, picker=True)
        plt.xlim(Qx,Qx+1)
        plt.ylim(Qy,Qy+1)
        #plt.xticks(x, labels, rotation='vertical')


def show_da(out_da, x_array, y_array, title=""):
    from matplotlib import pyplot as plt
    from numpy import linspace, max, min
    #print "time execution = ", time() - start , " s"
    nx = len(x_array)
    ny = len(y_array)
    #print(nx, ny, len(out_da))
    out_da = out_da.reshape(ny,nx)
    xmin, xmax, ymin, ymax = min(x_array), max(x_array), min(y_array), max(y_array)
    #plt.subplot(111, axisbg='darkslategray')
    extent = xmin, xmax, ymin, ymax
    #print extent
    #plt.savetxt("da.txt", da)
    plt.figure(figsize=(10, 7))
    fig1 = plt.contour(out_da, linewidths=2,extent = extent)#, colors = 'r')
    #fig1 = plt.contourf(out_da, 20,cmap=plt.cm.rainbow,extent = extent)#, colors = 'r')
    #plt.axis_bgcolor("#bdb76b")
    plt.grid(True)
    plt.title(title)
    plt.xlabel("X, m")
    plt.ylabel("Y, m")
    cb = plt.colorbar()
    cb.set_label('Nturns')
    #cb.ax.set_yticklabels(map(str, linspace(min(out_da), max(out_da), 5) ))

    #plt.savefig('da_error_'+str(int(np.random.rand()*100))+'.png')
    plt.show()

def show_mu(contour_da, mux, muy, x_array, y_array, zones = None ):
    from matplotlib import pyplot as plt

    nx = len(x_array)
    ny = len(y_array)
    t= linspace(0,3.14, num = 100)
    contour_da = contour_da.reshape(ny,nx)
    mux = mux.reshape(ny,nx)
    muy = muy.reshape(ny,nx)
    xmin, xmax, ymin, ymax = min(x_array), max(x_array), min(y_array), max(y_array)
    plt.figure(1,figsize=(10, 7)) #axisbg='darkslategray'
    extent = xmin, xmax, ymin, ymax

    my_cmap = plt.cm.Paired
    #my_cmap.set_under('w')
    #norm = mlb.colors.Normalize(vmin=-0.005, vmax=max(mux))
    fig1 = plt.contour(contour_da, 1,extent = extent, linewidths=2,colors='k')#, colors = 'r')
    fig1 = plt.contourf(mux,40, cmap=my_cmap, extent = extent)#, colors = 'r')
    cb = plt.colorbar(cmap=my_cmap)
    fig1 = plt.contourf(mux,10, levels=[-1,-.0001], colors='w',extent = extent)
    if zones != None:
        x_zone = zones[0]
        y_zone = zones[1]
        plt.plot(x_zone*cos(t), y_zone*sin(t), "g", lw = 2)
        plt.plot(2*x_zone*cos(t), 2*y_zone*sin(t), "b", lw = 2)
        plt.plot(3*x_zone*cos(t), 3*y_zone*sin(t), "r", lw = 2)
        plt.plot(4*x_zone*cos(t), 4*y_zone*sin(t), "y", lw = 2)
    plt.grid(True)
    #plt.figure(figsize=(10, 7))
    plt.xlabel("X, m")
    plt.ylabel("Y, m")
    cb.set_label('Qx')
    plt.figure(2,figsize=(10, 7))

    fig1 = plt.contour(contour_da, 1,extent = extent, linewidths=2,colors='k')#, colors = 'r')
    fig1 = plt.contourf(muy,40, cmap=my_cmap, extent = extent)#, colors = 'r')
    if zones != None:
        x_zone = zones[0]
        y_zone = zones[1]
        plt.plot(x_zone*cos(t), y_zone*sin(t), "g", lw = 2)
        plt.plot(2*x_zone*cos(t), 2*y_zone*sin(t), "b", lw = 2)
        plt.plot(3*x_zone*cos(t), 3*y_zone*sin(t), "r", lw = 2)
        plt.plot(4*x_zone*cos(t), 4*y_zone*sin(t), "y", lw = 2)
    #x = np.linspace(-, 0.01, 0.0001)
    #plt.plot()
    cb = plt.colorbar(cmap=my_cmap)
    fig1 = plt.contourf(muy,10, levels=[-1,-.0001], colors='w',extent = extent)
    plt.xlabel("X, m")
    plt.ylabel("Y, m")
    plt.grid(True)
    cb.set_label('Qy')
    plt.show()
