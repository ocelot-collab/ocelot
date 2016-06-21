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
        MainWindow.setWindowTitle('SASE optimization')

        ## Create docks, place them into the window one at a time.
        ## Note that size arguments are only a suggestion; docks will still have to
        ## fill the entire dock area and obey the limits of their internal widgets.
        #self.orb_fig = Dock("Orbit", size=(400, 300))     ## give this dock the minimum possible size
        #self.orb_fig.float()
        #self.sase_fig = Dock("SASE", size=(400,300), closable=True)
        self.sase_fig = Dock("SASE", size=(400,300))
        self.blm_fig = Dock("BLM", size=(400,200))
        self.seq_cntr = Dock("Sequence", size=(150,200))
        self.sase_cntr = Dock("Controls", size=(150,200))
        #self.orb_cntr = Dock("orb contr.", size=(400,100))
        self.cur_fig = Dock("Settings", size=(400,300))
        self.logger = Dock("Logger", size=(100,300))

        self.area.addDock(self.cur_fig, 'left')

        #self.area.addDock(self.orb_fig, 'above', self.cur_fig)      ## place d1 at left edge of dock area (it will fill the whole space since there are no other docks yet)
        self.area.addDock(self.sase_fig, 'bottom', self.cur_fig)     ## place d2 at right edge of dock area

        self.area.addDock(self.blm_fig, 'top', self.sase_fig)## place d3 at bottom edge of d1
        self.area.addDock(self.sase_cntr, 'right')  ## place d5 at left edge of d1
        self.area.addDock(self.seq_cntr, 'left', self.sase_cntr)     ## place d4 at right edge of dock area
        #self.area.addDock(self.orb_cntr, 'bottom', self.orb_fig)

        ## Test ability to move docks programatically after they have been placed
        #self.area.moveDock(self.sase_fig, 'bottom', self.orb_fig)     ## move d4 to top edge of d2
        #self.area.moveDock(self.blm_fig, 'bottom', self.sase_fig)   ## move d6 to stack on top of d4
        self.area.addDock(self.logger, 'bottom', self.sase_fig)
        self.area.moveDock(self.blm_fig, 'above', self.logger)

        ## Add widgets into each dock

        #add Logger
        self.log_lab = QtGui.QTextBrowser()
        #self.log_lab.verticalScrollBar().setValue(self.log_lab.verticalScrollBar().maximum())
        self.logger.addWidget(self.log_lab)

        ## first dock gets save/restore buttons
        self.t = ParameterTree()
        if params != None:
            self.p = Parameter.create(name='params', type='group', children=params)
            self.t.setParameters(self.p, showTop=False)

        self.t.setWindowTitle('SASE optimization')
        self.seq_cntr.addWidget(self.t)


        self.seq = pg.LayoutWidget()
        self.label = QtGui.QLabel("""sequence control""")
        self.add_seq_btn = QtGui.QPushButton('Add Action')
        self.save_seq_btn = QtGui.QPushButton('Save seqs')
        self.load_seq_btn = QtGui.QPushButton('Load seqs')
        #self.restoreBtn.setEnabled(False)
        self.seq.addWidget(self.label, row=0, col=0)
        self.seq.addWidget(self.add_seq_btn, row=1, col=0)
        self.seq.addWidget(self.save_seq_btn, row=2, col=0)
        self.seq.addWidget(self.load_seq_btn, row=3, col=0)
        self.seq_cntr.addWidget(self.seq)


        #Currents graphics
        self.t_cur_cntr = ParameterTree()
        #param = [{'name': 'Devices', 'type': 'list', 'values': {}, 'value': 0}]
        param = []
        self.p_cur_cntr = Parameter.create(name='control', type='group', children=param)
        self.t_cur_cntr.setParameters(self.p_cur_cntr, showTop=False)

        self.current = pg.PlotWidget(title="Settings")

        self.cur_fig.addWidget(self.current, row=0, col=0)
        self.cur_fig.addWidget(self.t_cur_cntr, row=0, col=1)



        #BLM graphics
        ## Hide title bar on dock 3
        #d3.hideTitleBar()
        self.blm = pg.PlotWidget(title="BLM")
        self.blm_fig.addWidget(self.blm)


        #SASE graphics

        self.sase = pg.PlotWidget(title="SASE")
        self.sase_fig.addWidget(self.sase)

        #controls
        self.w5 = pg.LayoutWidget()
        self.start_opt_btm = QtGui.QPushButton('start')

        params = [
            {'name': 'Basic opt. parameters', 'type': 'group', 'children': [
        {'name': 'debug', 'type': 'bool', 'value': False},
        {'name': 'logging', 'type': 'bool', 'value': True},
        {'name': 'log file', 'type': 'str', 'value': 'test.log'},
        {'name': 'timeout', 'type': 'float', 'value': 0.5, 'step': 0.1, 'limits': (0, 10)},
        #{'name': 'SASE det.', 'type': 'list', 'values': {'mcp': 'mcp', 'gmd slow':'gmd_fl1_slow', 'bkr':'bkr', 'default':'gmd_default'}, 'value': "default"},
        {'name': 'detector', 'type': 'list', 'values': ['gmd_default', 'mcp', 'gmd_fl1_slow', 'bkr'], 'value': "gmd_default"} ]}
        ]

        self.t_cntr = ParameterTree()
        self.p_cntr = Parameter.create(name='control', type='group', children=params)
        self.t_cntr.setParameters(self.p_cntr, showTop=False)


        self.restore_cur_btn = QtGui.QPushButton('Restore')
        #self.restore_cur_btn.setEnabled(False)
        self.setmax_opt_btn = QtGui.QPushButton('Set currents for max SASE')
        #self.setmax_opt_btn.setEnabled(False)
        self.stop_opt_btn = QtGui.QPushButton('Stop')
        self.clear_disp_btn = QtGui.QPushButton('Clear display')

        self.save_machine_btn = QtGui.QPushButton('Save new tuning')
        #self.save_machine_btn.setEnabled(False)
        #self.stop_btn = QtGui.QPushButton('stop')
        #self.w5.addWidget(self.start_opt_btm, row=0, col=0)
        self.w5.addWidget(self.stop_opt_btn, row=0, col=0)
        self.w5.addWidget(self.restore_cur_btn, row=1, col=0)
        #self.w5.addWidget(self.setmax_opt_btn, row=2, col=0)

        self.w5.addWidget(self.clear_disp_btn, row=4, col=0)

        #self.w5.addWidget(self.debug_opt_chk, row=3, col=0)
        #self.w5.addWidget(self.log_opt_chk, row=4, col=0)
        self.w5.addWidget(self.t_cntr, row=5, col=0)

        #self.w5.addWidget(QtGui.QLabel("""machine settings"""), row=6, col=0)
        #self.w5.addWidget(self.save_machine_btn, row=7, col=0)
        self.sase_cntr.addWidget(self.w5)

        # Orbit graphics
        #self.orbit = pg.PlotWidget(title="Orbit")
        #self.orbit = pg.GraphicsWindow(title="Orbit")
        #self.orb_fig.addWidget(self.orbit)

        # orbit controls
        #self.w6 = pg.LayoutWidget()
        #self.ref_btm = QtGui.QPushButton('Set ref. orbit')
        #self.save_btn = QtGui.QPushButton('Save')
        #self.save_btn.setEnabled(False)
        #self.load_btn = QtGui.QPushButton('Load')
        #self.load_btn.setEnabled(False)
        #self.w6.addWidget(self.ref_btm, row=0, col=0)
        #self.w6.addWidget(self.save_btn, row=0, col=1)
        #self.w6.addWidget(self.load_btn, row=0, col=2)
        #self.orb_cntr.addWidget(self.w6)
        #self.orb_fig.addWidget(self.w6)
        #return win



class Ui_ChildWindow(object):
    def setupUi(self, MainWindow, params):
        #app = QtGui.QApplication([])
        #self.win = QtGui.QMainWindow()
        self.area = DockArea()
        MainWindow.setCentralWidget(self.area)
        MainWindow.resize(500, 700)
        MainWindow.setWindowTitle('Action construct')


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
        self.label = QtGui.QLabel("""Controls""")
        self.saveBtn = QtGui.QPushButton('Add Action')
        self.restoreBtn = QtGui.QPushButton('Modify table')
        self.restoreBtn.setEnabled(False)
        self.seq.addWidget(self.label, row=0, col=0)
        self.seq.addWidget(self.saveBtn, row=1, col=0)
        self.seq.addWidget(self.restoreBtn, row=2, col=0)
        self.seq_cntr.addWidget(self.seq)
