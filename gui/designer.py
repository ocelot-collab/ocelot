import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.console
import numpy as np
from pyqtgraph.dockarea import *
import pyqtgraph.parametertree.parameterTypes as pTypes
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType


from ocelot import *


import sys
lattice_file = '/home/iagapov/workspace/petra/p4/rb/cell2rf.py'

sys.path.append('/home/iagapov/workspace/petra/p4/rb/')

exec( open(lattice_file))
lat = MagneticLattice(__lattice__)
tws = twiss(lat, __tws__, nPoints = 100)


app = QtGui.QApplication([])
win = QtGui.QMainWindow()
area = DockArea()
win.setCentralWidget(area)
win.resize(1000,500)
win.setWindowTitle('Lattice designer')

## Create docks, place them into the window one at a time.
## Note that size arguments are only a suggestion; docks will still have to
## fill the entire dock area and obey the limits of their internal widgets.
d1 = Dock("Dock1", size=(1, 1))     ## give this dock the minimum possible size
d2 = Dock("Dock2 - Console", size=(500,300), closable=True)
d4 = Dock("Twiss -- dispersion", size=(500,200))
d5 = Dock("Dock5 - Image", size=(500,200))
d6 = Dock("Twiss - beta", size=(500,200))
area.addDock(d1, 'left')      ## place d1 at left edge of dock area (it will fill the whole space since there are no other docks yet)
area.addDock(d2, 'right')     ## place d2 at right edge of dock area
area.addDock(d4, 'right')     ## place d4 at right edge of dock area
area.addDock(d5, 'left', d1)  ## place d5 at left edge of d1
area.addDock(d6, 'top', d4)   ## place d5 at top edge of d4

## Test ability to move docks programatically after they have been placed
area.moveDock(d4, 'top', d2)     ## move d4 to top edge of d2
area.moveDock(d6, 'above', d4)   ## move d6 to stack on top of d4
area.moveDock(d5, 'top', d2)     ## move d5 to top edge of d2


## Add widgets into each dock

## first dock gets save/restore buttons
w1 = pg.LayoutWidget()
label = QtGui.QLabel(""" -- DockArea Example -- 
Place tree widget here
""")
saveBtn = QtGui.QPushButton('Save dock state')
restoreBtn = QtGui.QPushButton('Restore dock state')
restoreBtn.setEnabled(False)
w1.addWidget(label, row=0, col=0)
w1.addWidget(saveBtn, row=1, col=0)
w1.addWidget(restoreBtn, row=2, col=0)
d1.addWidget(w1)
state = None
def save():
    global state
    state = area.saveState()
    restoreBtn.setEnabled(True)
    
    print 'saving something...'
    
def load():
    global state
    area.restoreState(state)
    
    print 'loading something...'
    
    
saveBtn.clicked.connect(save)
restoreBtn.clicked.connect(load)


w2 = pg.console.ConsoleWidget()
d2.addWidget(w2)


w4 = pg.PlotWidget(title="Dock 4 plot")
w4.plot([t.s for t in tws],[t.Dx for t in tws], pen=(0, 0, 255))

d4.addWidget(w4)



params = [
    {'name': 'Quadrupoles', 'type': 'group', 'children': []},
    {'name': 'Drifts', 'type': 'group', 'children': []},
    {'name': 'Bends', 'type': 'group', 'children': []},
    {'name': 'BendsFocusing', 'type': 'group', 'children': []},
    {'name': 'Twiss', 'type': 'group', 'children': [
        {'name': 'Save State', 'type': 'action'},
        {'name': 'Restore State', 'type': 'action', 'children': [
            {'name': 'Add missing items', 'type': 'bool', 'value': True},
            {'name': 'Remove extra items', 'type': 'bool', 'value': True},
        ]},
    ]},
]

for v in sorted( globals().keys() ) : 
    try: 
        if v == 'd1': print 'd1'
        if globals()[v].__class__ == Quadrupole: params[0]['children'].append({'name': v, 'type': 'float', 'value': globals()[v].k1, 'step': 0.1})
        if globals()[v].__class__ == Drift: params[1]['children'].append({'name': v, 'type': 'float', 'value': globals()[v].l, 'step': 0.1})
        if globals()[v].__class__ in (SBend, RBend, Bend): params[2]['children'].append({'name': v, 'type': 'float', 'value': globals()[v].angle, 'step': 0.1})
    except:
        pass

p = Parameter.create(name='params', type='group', children=params)


w5 = ParameterTree()
w5.setParameters(p, showTop=False)
w5.setWindowTitle('pyqtgraph example: Parameter Tree')


## If anything changes in the tree, print a message
def change(param, changes):
    print("tree changes:")
    for param, change, data in changes:
        path = p.childPath(param)
        if path is not None:
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


        lat.update_transfer_maps()
        tws = twiss(lat, __tws__, nPoints = 100)
        
        w6.clear()
        w4.clear()
        
        w6.plot([t.s for t in tws],[t.beta_x for t in tws], pen=(255, 0, 0))
        w6.plot([t.s for t in tws],[t.beta_y for t in tws], pen=(0, 255, 0))
        w4.plot([t.s for t in tws],[t.Dx for t in tws], pen=(0, 0, 255))

        
    
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
d5.addWidget(w5)



w6 = pg.PlotWidget(title="Dock 6 plot")
w6.plot([t.s for t in tws],[t.beta_x for t in tws], pen=(255, 0, 0))
w6.plot([t.s for t in tws],[t.beta_y for t in tws], pen=(0, 255, 0))

d6.addWidget(w6)



win.show()



## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
