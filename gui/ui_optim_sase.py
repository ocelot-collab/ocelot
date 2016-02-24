__author__ = 'sergey'
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
from pyqtgraph.dockarea import *

import pyqtgraph.parametertree.parameterTypes as pTypes
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
import numpy as np



class Ui_MainWindow(object):
    def setupUi(self, MainWindow, params):
        #app = QtGui.QApplication([])
        #self.win = QtGui.QMainWindow()
        self.area = DockArea()
        MainWindow.setCentralWidget(self.area)
        MainWindow.resize(1200,600)
        MainWindow.setWindowTitle('pyqtgraph example: dockarea')

        ## Create docks, place them into the window one at a time.
        ## Note that size arguments are only a suggestion; docks will still have to
        ## fill the entire dock area and obey the limits of their internal widgets.
        self.orb_fig = Dock("Orbit", size=(400, 300))     ## give this dock the minimum possible size
        self.sase_fig = Dock("SASE", size=(400,300), closable=True)
        self.blm_fig = Dock("BLM", size=(400,200))
        self.seq_cntr = Dock("Sequence", size=(150,200))
        self.sase_cntr = Dock("Controls", size=(50,200))
        self.orb_cntr = Dock("orb contr.", size=(400,100))


        self.area.addDock(self.orb_fig, 'left')      ## place d1 at left edge of dock area (it will fill the whole space since there are no other docks yet)
        self.area.addDock(self.sase_fig, 'right')     ## place d2 at right edge of dock area
        self.area.addDock(self.blm_fig, 'bottom', self.sase_fig)## place d3 at bottom edge of d1
        self.area.addDock(self.sase_cntr, 'right')  ## place d5 at left edge of d1
        self.area.addDock(self.seq_cntr, 'left', self.sase_cntr)     ## place d4 at right edge of dock area
        self.area.addDock(self.orb_cntr, 'bottom', self.orb_fig)

        ## Test ability to move docks programatically after they have been placed
        self.area.moveDock(self.sase_fig, 'bottom', self.orb_cntr)     ## move d4 to top edge of d2
        self.area.moveDock(self.blm_fig, 'bottom', self.sase_fig)   ## move d6 to stack on top of d4



        ## Add widgets into each dock

        ## first dock gets save/restore buttons
        self.t = ParameterTree()
        if params != None:
            self.p = Parameter.create(name='params', type='group', children=params)
            self.t.setParameters(self.p, showTop=False)

        self.t.setWindowTitle('pyqtgraph example: Parameter Tree')
        self.seq_cntr.addWidget(self.t)


        self.seq = pg.LayoutWidget()
        self.label = QtGui.QLabel("""sequence control""")
        self.add_seq_btn = QtGui.QPushButton('Add seq')
        self.save_seq_btn = QtGui.QPushButton('save seqs')
        self.load_seq_btn = QtGui.QPushButton('load seqs')
        #self.restoreBtn.setEnabled(False)
        self.seq.addWidget(self.label, row=0, col=0)
        self.seq.addWidget(self.add_seq_btn, row=1, col=0)
        self.seq.addWidget(self.save_seq_btn, row=2, col=0)
        self.seq.addWidget(self.load_seq_btn, row=3, col=0)
        self.seq_cntr.addWidget(self.seq)





        # Orbit graphics
        #self.orbit = pg.PlotWidget(title="Orbit")
        self.orbit = pg.GraphicsWindow(title="Basic plotting examples")
        #self.w1.plot(np.random.normal(size=100))
        self.orb_fig.addWidget(self.orbit)



        #BLM graphics
        ## Hide title bar on dock 3
        #d3.hideTitleBar()
        self.blm = pg.PlotWidget(title="BLM")
        #self.blm.plot(np.random.normal(size=100))
        self.blm_fig.addWidget(self.blm)


        #SASE graphics
        self.sase = pg.PlotWidget(title="SASE")
        #self.sase.plot(np.random.normal(size=100))
        self.sase_fig.addWidget(self.sase)

        #controls
        self.w5 = pg.LayoutWidget()
        self.start_btm = QtGui.QPushButton('start')
        self.restore_btn = QtGui.QPushButton('restore')
        self.stop_btn = QtGui.QPushButton('stop')
        #self.stop_btn = QtGui.QPushButton('stop')
        self.w5.addWidget(self.start_btm, row=0, col=0)
        self.w5.addWidget(self.restore_btn, row=1, col=0)
        self.w5.addWidget(self.stop_btn, row=2, col=0)
        self.sase_cntr.addWidget(self.w5)


        # orbit controls
        self.w6 = pg.LayoutWidget()
        self.ref_btm = QtGui.QPushButton('ref orbit')
        self.save_btn = QtGui.QPushButton('save')
        self.load_btn = QtGui.QPushButton('load')
        self.w6.addWidget(self.ref_btm, row=0, col=0)
        self.w6.addWidget(self.save_btn, row=0, col=1)
        self.w6.addWidget(self.load_btn, row=0, col=2)
        self.orb_cntr.addWidget(self.w6)
        #return win



class Ui_ChildWindow(object):
    def setupUi(self, MainWindow, params):
        #app = QtGui.QApplication([])
        #self.win = QtGui.QMainWindow()
        self.area = DockArea()
        MainWindow.setCentralWidget(self.area)
        MainWindow.resize(500, 700)
        MainWindow.setWindowTitle('pyqtgraph example: dockarea')


        self.seq_cntr = Dock("Sequence", size=(150,200))
        self.area.addDock(self.seq_cntr, 'left')

        ## first dock gets save/restore buttons
        self.t = ParameterTree()


        if params != None:
            self.p_child = Parameter.create(name='params', type='group', children=params)
            self.t.setParameters(self.p_child, showTop=False)

        self.t.setWindowTitle('pyqtgraph example: Parameter Tree')
        self.seq_cntr.addWidget(self.t)


        self.seq = pg.LayoutWidget()
        self.label = QtGui.QLabel(""" -- DockArea Example --""")
        self.saveBtn = QtGui.QPushButton('Save dock state')
        self.restoreBtn = QtGui.QPushButton('Restore dock state')
        self.restoreBtn.setEnabled(False)
        self.seq.addWidget(self.label, row=0, col=0)
        self.seq.addWidget(self.saveBtn, row=1, col=0)
        self.seq.addWidget(self.restoreBtn, row=2, col=0)
        self.seq_cntr.addWidget(self.seq)
