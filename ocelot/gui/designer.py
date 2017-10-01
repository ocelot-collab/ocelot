from __future__ import print_function
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.console
import numpy as np
from pyqtgraph.dockarea import *
import pyqtgraph.parametertree.parameterTypes as pTypes
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType


from ocelot import *


import sys

sys.path.append('/home/iagapov/workspace/petra/p4/rb/ocelot/')

lattice_file = '/home/iagapov/workspace/petra/p4/rb/ocelot/undcell_v2.py'


'''
def load_file(lattice_file):

    #global lat, tws, lattice_names, old_quad_values, __knobs__, __lattice__, __lattice_list__, __tws__
    
    exec( open(lattice_file))

    print 'local:', locals()[tws]
    print 'global:',globals()[tws]
    
    print tws
    print __lattice__
    

    
    
    lat = MagneticLattice(__lattice__)
    tws = twiss(lat, __tws__, nPoints = 100)
    
    
    lattice_names = sorted(__lattice_list__.keys())
    
    old_quad_values = {}
    if __knobs__ != None:
        for v in sorted( globals().keys() ) : 
            try: 
                if globals()[v].__class__ == Quadrupole: 
                    old_quad_values[globals()[v]] = globals()[v].k1
                    print 'stored value for', globals()[v].id
            except:
                pass

__periodic__ = False
n_twiss_points = 100
    



load_file(lattice_file)
'''
exec( open(lattice_file))


lat = MagneticLattice(__lattice__)
tws = twiss(lat, __tws__, nPoints = 100)


lattice_names = sorted(__lattice_list__.keys())

old_quad_values = {}
if __knobs__ != None:
    for v in sorted( globals().keys() ) : 
        try: 
            if globals()[v].__class__ == Quadrupole: 
                old_quad_values[globals()[v]] = globals()[v].k1
                print ('stored value for', globals()[v].id)
        except:
            pass

__periodic__ = False
n_twiss_points = 100




app = QtGui.QApplication([])
win = QtGui.QMainWindow()
area = DockArea()
win.setCentralWidget(area)
win.resize(1000,500)
win.setWindowTitle('Lattice designer')

## Create docks, place them into the window one at a time.
## Note that size arguments are only a suggestion; docks will still have to
## fill the entire dock area and obey the limits of their internal widgets.
dock1 = Dock("Dock1", size=(1, 1))     ## give this dock the minimum possible size
dock2 = Dock("Dock2 - Console", size=(500,300), closable=True)
dock4 = Dock("Twiss -- dispersion", size=(500,200))
dock5 = Dock("Dock5 - Image", size=(500,200))
d6 = Dock("Twiss - beta", size=(500,200))
area.addDock(dock1, 'left')      ## place dock1 at left edge of dock area (it will fill the whole space since there are no other docks yet)
area.addDock(dock2, 'right')     ## place dock2 at right edge of dock area
area.addDock(dock4, 'right')     ## place dock4 at right edge of dock area
area.addDock(dock5, 'left', dock1)  ## place dock5 at left edge of dock1
area.addDock(d6, 'top', dock4)   ## place dock5 at top edge of dock4

## Test ability to move docks programatically after they have been placed
area.moveDock(dock4, 'top', dock2)     ## move dock4 to top edge of dock2
area.moveDock(d6, 'above', dock4)   ## move d6 to stack on top of dock4
area.moveDock(dock5, 'top', dock2)     ## move dock5 to top edge of dock2


## Add widgets into each dock

## first dock gets save/restore buttons
w1 = pg.LayoutWidget()
label = QtGui.QLabel(""" -- DockArea Example -- 
Place tree widget here
""")
saveBtn = QtGui.QPushButton('Save')
loadBtn = QtGui.QPushButton('Load')
loadBtn.setEnabled(True)
w1.addWidget(label, row=0, col=0)
w1.addWidget(saveBtn, row=1, col=0)
w1.addWidget(loadBtn, row=2, col=0)
dock1.addWidget(w1)
state = None

def save():
    global state
    state = area.saveState()
    loadBtn.setEnabled(True)
    
    print ('saving something...')
    
def load():
    #global state
    #area.restoreState(state)
    filename = QtGui.QFileDialog.getOpenFileName(win, 'Open File', '/')
    print (filename)
    
    print('loading {}'.format(filename))
    
    
saveBtn.clicked.connect(save)
loadBtn.clicked.connect(load)


w2 = pg.console.ConsoleWidget()
dock2.addWidget(w2)


w4 = pg.PlotWidget(title="Dock 4 plot")
w4.plot([t.s for t in tws],[t.Dx for t in tws], pen=(0, 0, 255))
w4.plot([t.s for t in tws],[t.Dxp for t in tws], pen=(255, 0, 0))

dock4.addWidget(w4)



params = [
    {'name': 'Quadrupoles', 'type': 'group', 'children': []},
    {'name': 'Drifts', 'type': 'group', 'children': []},
    {'name': 'Bends', 'type': 'group', 'children': []},
    {'name': 'BendsFocusing', 'type': 'group', 'children': []},
    {'name': 'Twiss', 'type': 'group', 'children': [
        {'name': 'Lattices', 'type': 'list', 'values': lattice_names, 'value': lattice_names[0]},
        {'name': 'Periodic', 'type': 'bool', 'value': False},
        {'name': 'npoints', 'type': 'int', 'values': n_twiss_points, 'value': 100},]},
    {'name': 'Knobs', 'type': 'group', 'children': []},
    {'name': 'ShowParameters', 'type': 'group', 'children': [
        {'name': 'Save State', 'type': 'action'},
        {'name': 'Restore State', 'type': 'action', 'children': [
            {'name': 'Add missing items', 'type': 'bool', 'value': True},
            {'name': 'Remove extra items', 'type': 'bool', 'value': True},
        ]},
    ]},
]

for v in sorted( globals().keys() ) : 
    try: 
        if globals()[v].__class__ == Quadrupole: params[0]['children'].append({'name': v, 'type': 'float', 'value': globals()[v].k1, 'step': 0.1})
        if globals()[v].__class__ == Drift: params[1]['children'].append({'name': v, 'type': 'float', 'value': globals()[v].l, 'step': 0.1})
        if globals()[v].__class__ in (SBend, RBend, Bend): params[2]['children'].append({'name': v, 'type': 'float', 'value': globals()[v].angle, 'step': 0.1})
    except:
        pass
    
if __tws__ != None:
    params[4]['children'].append({'name': 'beta_x', 'type': 'float', 'value': __tws__.beta_x, 'step': 0.1})
    params[4]['children'].append({'name': 'alpha_x', 'type': 'float', 'value': __tws__.alpha_x, 'step': 0.1})
    params[4]['children'].append({'name': 'beta_y', 'type': 'float', 'value': __tws__.beta_y, 'step': 0.1})
    params[4]['children'].append({'name': 'alpha_y', 'type': 'float', 'value': __tws__.alpha_y, 'step': 0.1})

if __knobs__ != None:
    for k in __knobs__.keys():
        params[5]['children'].append({'name': k, 'type': 'float', 'value': 0.0, 'step': 0.1})


p = Parameter.create(name='params', type='group', children=params)


w5 = ParameterTree()
w5.setParameters(p, showTop=False)
w5.setWindowTitle('pyqtgraph example: Parameter Tree')


## If anything changes in the tree, print a message
def change(param, changes):
    global __periodic__, __tws__, lat, n_twiss_points
    global __knobs__
    print("tree changes:")
    for param, change, data in changes:
        path = p.childPath(param)
        if path is not None:
            print (path)
            childName = '.'.join(path)
            parName = path[1]
            groupName = path[0]
        else:
            childName = param.name()
        print('  parameter (full): %s'% childName)
        print('  parameter: %s'% parName)
        print('  change:    %s'% change)
        print('  data:      %s'% str(data))
        print('  ----------')
        
        #print globals()
        
        if groupName == 'Quadrupoles':
            globals()[parName].k1 = float(data)
        if groupName == 'Drifts':
            globals()[parName].l = float(data)
        if groupName == 'Bends':
            globals()[parName].angle = float(data)
        if groupName == 'Knobs':
            knob = __knobs__[parName]
            for e in knob.keys():
                e.k1 = old_quad_values[e] + float(data) * knob[e]
        if groupName == 'Twiss':
            if parName == 'Periodic':
                __periodic__ = data
            if parName == 'npoints':
                n_twiss_points = int(data)
            if parName in ('beta_x','alpha_x','beta_y','alpha_y','Dx','Dxp','Dy','Dyp'):
                __tws__.__dict__[parName] = float(data)
            if parName == 'Lattices':
                __lattice__, __tws__ =  __lattice_list__[str(data)]
                lat = MagneticLattice(__lattice__)



        lat.update_transfer_maps()
        
        if __periodic__:
            tws = twiss(lat, None, nPoints = n_twiss_points)
        else:
            tws = twiss(lat, __tws__, nPoints = n_twiss_points)
            
        w6.clear()
        w4.clear()
        
        w6.plot([t.s for t in tws],[t.beta_x for t in tws], pen=(255, 0, 0))
        w6.plot([t.s for t in tws],[t.beta_y for t in tws], pen=(0, 255, 0))
        w4.plot([t.s for t in tws],[t.Dx for t in tws], pen=(0, 0, 255))
        w4.plot([t.s for t in tws],[t.Dxp for t in tws], pen=(255, 0, 0))

        label.setText( 'Twiss Start:\n{}\nTwiss End:\n{}'.format(tws[0], tws[-1]) )

        
    
p.sigTreeStateChanged.connect(change)


def valueChanging(param, value):
    print("Value changing (not finalized):", param, value)
    
# Too lazy for recursion:
for child in p.children():
    child.sigValueChanging.connect(valueChanging)
    for ch2 in child.children():
        ch2.sigValueChanging.connect(valueChanging)



#w5 = pg.ImageView()
#w5.setImage(np.random.normal(size=(100,100)))
dock5.addWidget(w5)



w6 = pg.PlotWidget(title=u"\u03b2 x,y") # beta
w6.plot([t.s for t in tws],[t.beta_x for t in tws], pen=(255, 0, 0))
w6.plot([t.s for t in tws],[t.beta_y for t in tws], pen=(0, 255, 0))

d6.addWidget(w6)



win.show()



## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
