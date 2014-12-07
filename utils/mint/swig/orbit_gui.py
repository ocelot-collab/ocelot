#from __future__ import unicode_literals
import random
import pickle

from PyQt4 import QtGui, QtCore

from numpy import arange, sin, pi, abs, max, sum, mean, array
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import ocelot.utils.mint.mint as mint
import ocelot.utils.mint.swig.dcs as dcs
import os, sys

progname = os.path.basename(sys.argv[0])

class MyData:
    def __init__(self):
        self.i = 0


bpm_names = ['2UND1', '2UND2', '2UND4', '2UND5', '2UND6']
blm_names = ['1L.UND1','1R.UND1', '1L.UND2', '1R.UND2', '1R.UND3', '1R.UND4', '1R.UND5']

gmd_channel = 'TTF2.DAQ/PHFLUX/OUT33/VAL' # ??


bpms = []

for bpm_name in bpm_names:
	bpm = dcs.BPM("TTF2.DIAG/BPM/" + bpm_name)
	bpms.append(bpm)


blms = []

		
for blm_name in blm_names:
	blm = dcs.Device("TTF2.DIAG/BLM/" + blm_name)
	dcs.get_device_info(blm)
	blm_channel = blm.id + '/CH00.TD'
	print 'blm info:', blm.id, bpm.z_pos
	blms.append(blm)



class MyMplCanvas(FigureCanvas):
    """Ultimately, this is a QWidget (as well as a FigureCanvasAgg, etc.)."""
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        # We want the axes cleared every time plot() is called
        self.axes.hold(False)

        self.compute_initial_figure()

        #
        FigureCanvas.__init__(self, fig)
        self.setParent(parent)

        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def compute_initial_figure(self):
        pass


class MyStaticMplCanvas(MyMplCanvas):
    """Simple canvas with a sine plot."""
    def compute_initial_figure(self):
        t = arange(0.0, 3.0, 0.01)
        s = sin(2*pi*t)
        self.axes.plot(t, s)


class XYCanvas(MyMplCanvas):
    """A canvas that updates itself every second with a new plot."""
    def __init__(self, *args, **kwargs):
        MyMplCanvas.__init__(self, *args, **kwargs)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(1000)

    def compute_initial_figure(self):
        self.axes.bar(range(14), range(14), width=0.2, color='b')
	self.axes.set_ylim((-4,4))

    def update_figure(self):
        
	for bpm in bpms:
		dcs.get_bpm_val(bpm)
		print 'bpm read:', bpm.id, bpm.x, bpm.y
		
        x = [bpm.x for bpm in bpms]
	y = [bpm.y for bpm in bpms]
        z = [bpm.z_pos for bpm in bpms]
	
	self.axes.bar(z, x, width=0.2, color='b', alpha=0.5)
	self.axes.hold(True)
	self.axes.bar(z, y, width=0.2, color='g', alpha=0.5)
	self.axes.hold(False)
	self.axes.set_ylim((-4,4))
            
        self.axes.set_title('Orbit')
        self.draw()

class BLMCanvas(MyMplCanvas):
    """A canvas that updates itself every second with a new plot."""
    def __init__(self, *args, **kwargs):
        MyMplCanvas.__init__(self, *args, **kwargs)
	self.plane = 'x'
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(1000)

    def compute_initial_figure(self):
        self.axes.bar(range(14), range(14), width=0.2, color='b')
	#self.axes.set_ylim((-4,4))

    def update_figure(self):

        z = [blm.z_pos for blm in blms]   
	loss = {}
	
	for blm in blms:
		blm_channel = blm.id + '/CH00.TD'
		print 'blm info:', blm.id, blm.z_pos
		h = array(dcs.get_device_td(blm_channel))
		loss[blm] = abs( sum(h) )
		print 'blm read:', blm.id, loss[blm]
		
        y = [loss[blm] for blm in blms]

	
	self.axes.bar(z, y, width=0.2, color='b')
            
        self.axes.set_title('Beam loss')
        self.draw()


class SASECanvas(MyMplCanvas):
    """A canvas that updates itself every second with a new plot."""
    def __init__(self, *args, **kwargs):
        MyMplCanvas.__init__(self, *args, **kwargs)
	self.data = []
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(1000)

    def compute_initial_figure(self):
        self.axes.plot(range(14), range(14), color='b')
	#self.axes.set_ylim((-4,4))

    def update_figure(self):
	val = dcs.get_device_val(gmd_channel)
	self.data.append(val)
	self.axes.plot(xrange(len(self.data)), self.data, color='b')
            
        self.axes.set_title('SASE')
        self.draw()

	
class ApplicationWindow(QtGui.QMainWindow):
    def __init__(self):
                
        QtGui.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("application main window")

        self.file_menu = QtGui.QMenu('&File', self)
        self.file_menu.addAction('&Quit', self.fileQuit,
                                 QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.help_menu = QtGui.QMenu('&Help', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)

        self.main_widget = QtGui.QWidget(self)

        l = QtGui.QVBoxLayout(self.main_widget)
        dc1 = XYCanvas(self.main_widget, width=5, height=4, dpi=100)
        l.addWidget(dc1)

        dc2 = BLMCanvas(self.main_widget, width=5, height=4, dpi=100)
        l.addWidget(dc2)

        dc3 = SASECanvas(self.main_widget, width=5, height=4, dpi=100)
        l.addWidget(dc3)


        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        
    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()
    

qApp = QtGui.QApplication(sys.argv)

aw = ApplicationWindow()
aw.setWindowTitle("%s" % progname)
aw.show()
sys.exit(qApp.exec_())
#qApp.exec_()
