from PyQt5 import QtGui
import pyqtgraph as pg

from ocelot.cpbd.elements import *
from app.parser import *
from app.ui_forms.widgets.tune_panel import *


class LatticeElementsPics():

    def __init__(self, lattice):
        
        self.scale = 1.0
        self.h_drift = 0.01

        self.lattice = lattice
        self.elements_pics = []

    
    def __del__(self):
        pass


    def calc_lattice_elements_pics(self):
        '''prepare lattice elements to draw them on PG plot'''

        elements = {'bends': [], 'cors': [], 'quads': [], 'sexts': [], 'octs': [], 'mults': [], 'cavs': [], 'unduls': [], 'monits': [], 'drifts': []}
        e_max =    {'bends':0.0, 'cors':1.0, 'quads':0.0, 'sexts':0.0, 'octs':0.0, 'mults':0.0, 'cavs':1.0, 'unduls':0.0, 'monits':1.0, 'drifts':1.0}
        colors = {'bends':(135,206,250), 'cors':(0,255,55), 'quads':(255,0,0), 'sexts':(0,128,0), 'octs':(0,128,0), 'mults':(0,128,0), 'cavs':(255,165,0), 'unduls':(255,192,203), 'monits':(255,165,0), 'drifts':(255,255,255)}
        color_white = QtGui.QColor(255, 255, 255)

        pos = 0.0
        for elem in self.lattice.lattice.sequence:

            if elem.l != 0.0:
                l = elem.l
                d_pos = 0.0
            else:
                l = 0.001
                d_pos = l / 2.0
            
            if elem.__class__ == Edge:
                pass
            
            elif elem.__class__ in [Bend, RBend, SBend]:
                field = elem.angle / elem.l if elem.l != 0.0 else 1.0
                h = np.abs(field)
                y = (field - h) * 0.5
                elements['bends'].append([h, pos-d_pos, y, l, elem])
                if e_max['bends'] < np.abs(field):
                    e_max['bends'] = np.abs(field)
            
            elif elem.__class__ in [Hcor, Vcor]:
                h = np.abs(elem.angle)
                y = (elem.angle - h) * 0.5
                elements['cors'].append([h, pos-d_pos, y, l, elem])
                if e_max['cors'] < np.abs(elem.angle):
                    e_max['cors'] = np.abs(elem.angle)
            
            elif elem.__class__ == Quadrupole:
                h = np.abs(elem.k1)
                y = (elem.k1 - h) * 0.5
                elements['quads'].append([h, pos-d_pos, y, l, elem])
                if e_max['quads'] < np.abs(elem.k1):
                    e_max['quads'] = np.abs(elem.k1)
            
            elif elem.__class__ == Sextupole:
                h = np.abs(elem.k2)
                y = (elem.k2 - h) * 0.5
                elements['sexts'].append([h, pos-d_pos, y, l, elem])
                if e_max['sexts'] < np.abs(elem.k2):
                    e_max['sexts'] = np.abs(elem.k2)
            
            elif elem.__class__ == Octupole:
                h = np.abs(elem.k3)
                y = (elem.k3 - h) * 0.5
                elements['octs'].append([h, pos-d_pos, y, l, elem])
                if e_max['octs'] < np.abs(elem.k3):
                    e_max['octs'] = np.abs(elem.k3)
            
            elif elem.__class__ == Multipole:
                h = np.abs(elem.kn)
                y = (elem.kn - h) * 0.5
                elements['mults'].append([h, pos-d_pos, y, l, elem])
                if e_max['mults'] < np.abs(elem.kn):
                    e_max['mults'] = np.abs(elem.kn)
            
            elif elem.__class__ == Cavity:
                elements['cavs'].append([2.0, pos-d_pos, -1.0-self.h_drift, l, elem])
            
            elif elem.__class__ == Undulator:
                h = np.abs(elem.Kx)
                y = (elem.Kx - h) * 0.5
                elements['unduls'].append([h, pos-d_pos, y, l, elem])
                if e_max['unduls'] < np.abs(elem.Kx + elem.Ky):
                    e_max['unduls'] = np.abs(elem.Kx + elem.Ky)
            
            elif elem.__class__ == Monitor:
                elements['monits'].append([0.1, pos-d_pos, -0.05, l, elem])
            
            else:
                elements['drifts'].append([2.0*self.h_drift, pos-d_pos, -self.h_drift, l, elem])

            pos += elem.l

        self.elements_pics = []
        zero_h = []

        for elem_class in elements:
            if elem_class in ['drifts', 'monits']:
                continue
            for elem in elements[elem_class]:
                dy = self.h_drift if elem[2] == 0.0 else -self.h_drift
                y = elem[2]/e_max[elem_class] - dy
                h = elem[0]/e_max[elem_class]
                if h == 0.0:
                    zero_h.append(elem)
                    continue
                item = callbackElement(elem[1], y*self.scale, elem[3], h*self.scale)
                item.element = elem[4]
                item.setAcceptHoverEvents(True)
                item.color = colors[elem_class]
                item.setPen(pg.mkPen(colors[elem_class], width=0.0))
                item.setBrush(pg.mkBrush(colors[elem_class]))
                self.elements_pics.append(item)
        
        # add white line
        item = callbackElement(0.0, -self.h_drift*self.scale, self.lattice.lattice.totalLen, 2.0*self.h_drift*self.scale)
        item.color = color_white
        item.setPen(pg.mkPen(color_white, width=0.0))
        item.setBrush(pg.mkBrush(color_white))
        self.elements_pics.append(item)

        # add drifts
        for elem in elements['drifts']:
            item = callbackElement(elem[1], elem[2]*self.scale, elem[3], elem[0]*self.scale)
            item.element = elem[4]
            item.setAcceptHoverEvents(True)
            item.color = colors['drifts']
            item.setPen(pg.mkPen(colors['drifts'], width=0.0))
            item.setBrush(pg.mkBrush(colors['drifts']))
            self.elements_pics.append(item)

        # add monitors
        for elem in elements['monits']:
            item = callbackElement(elem[1], elem[2]*self.scale, elem[3], elem[0]*self.scale)
            item.element = elem[4]
            item.setAcceptHoverEvents(True)
            item.color = colors['monits']
            item.setPen(pg.mkPen(colors['monits'], width=0.0))
            item.setBrush(pg.mkBrush(colors['monits']))
            self.elements_pics.append(item)

        # add zero strength elements
        for elem in zero_h:
            item = callbackElement(elem[1], -self.h_drift*self.scale, elem[3], 2.0*self.h_drift*self.scale)
            item.element = elem[4]
            item.setAcceptHoverEvents(True)
            item.color = colors['drifts']
            item.setPen(pg.mkPen(colors['drifts'], width=0.0))
            item.setBrush(pg.mkBrush(colors['drifts']))
            self.elements_pics.append(item)
        

class callbackElement(QtGui.QGraphicsRectItem):
    '''Element call-back'''

    def add_tune_block(self):
        pass


    def mouseDoubleClickEvent(self, event):
        
        if not self.element.is_tuneable:
            return

        brush = pg.mkBrush(self.color)
        pg.QtGui.QGraphicsRectItem.setBrush(self, brush)
        
        self.add_tune_block(self.element)

    
    def hoverMoveEvent(self, event):
        pass


    def hoverEnterEvent(self, event):
        
        message = 'Element: ' + self.element.id
        
        if self.element.is_tuneable:
            message += '  (double click to select element)'

            color = QtGui.QColor(255, 250, 250)
            pg.QtGui.QGraphicsRectItem.setBrush(self, color)
        
        self.statusBar.showMessage(message)


    def hoverLeaveEvent(self, event):
        
        self.statusBar.showMessage('')

        brush = pg.mkBrush(self.color)
        pg.QtGui.QGraphicsRectItem.setBrush(self, brush)
