#from __future__ import unicode_literals
import sys
import os
import random
import pickle
import numpy as np

from sim_info import SimInfo, RunInfo

    
from PyQt4 import QtGui, QtCore

from numpy import arange, sin, pi
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from ocelot.common.math import fwhm

progname = os.path.basename(sys.argv[0])
data_file = sys.argv[1]


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


class RunStatCanvas(MyMplCanvas):
    """A canvas that updates itself every second with a new plot."""
    def __init__(self, *args, **kwargs):
        MyMplCanvas.__init__(self, *args, **kwargs)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(100)

    def compute_initial_figure(self):
        self.axes.bar([0,], [0,], width=0.2, color='b')

    def update_figure(self):
        print ('drawing', self.data.run_ids)

        stages = []
        for r_id in self.data.run_ids:
            stages.append(self.data.runs[r_id].stage)

        rects = self.axes.bar(self.data.run_ids, self.data.max_powers, width=0.2, color='b', alpha=0.5)

        for i in range(len(rects)):
            rect = rects[i]
            height = rect.get_height()
            self.axes.text(rect.get_x()+rect.get_width()/2., 0.50*height, 'st. %d'%int(stages[i]), ha='center', va='bottom')


        self.axes.set_title('run statistics')
        self.draw()

class RunPulseCanvas(MyMplCanvas):
    def __init__(self, *args, **kwargs):
        MyMplCanvas.__init__(self, *args, **kwargs)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(100)

    def compute_initial_figure(self):
        self.axes.plot([0,], [0,], color='r', alpha=0.6)

    def update_figure(self):

        print ('drawing', self.data.show_run_id)
        try:
            print ('have runs', self.data.runs.keys())
            r = self.data.runs[self.data.show_run_id]

            #self.axes.set_yscale('log')
            #self.axes.hold(True)
            p1, = self.axes.plot(r.t, r.power, color='r', alpha=0.6)

            try:
                self.axes.hold(True)
                p2, = self.axes.plot(r.t, r.power_ref, color='black')
                self.axes.hold(False)
            except:
                print ('run has no power_ref')
                self.axes.hold(False)

            self.axes.legend([p1], ['run ' + str(self.data.show_run_id)])            
            #self.axes.set_yscale('log')
            
            self.draw()
        except Exception as e:
            print ('run not ready...')
            raise(e)

class RunSpecCanvas(MyMplCanvas):
    def __init__(self, *args, **kwargs):
        MyMplCanvas.__init__(self, *args, **kwargs)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(100)

    def compute_initial_figure(self):
        self.axes.plot([0,], [0,], color='r')

    def update_figure(self):

        try:
            print ('drawing', self.data.show_run_id)
            print ('have runs', self.data.runs.keys())
            r = self.data.runs[self.data.show_run_id]
            
            
            si = np.zeros_like( np.abs(r.spec))
            start_id = 0
            for i in range(start_id,self.data.show_run_id +1):
                print ('averaging run', i)
                si += np.roll( np.abs(self.data.runs[i].spec)**2, len(r.freq_ev)/2) / int(self.data.show_run_id +1 - start_id)
            

            x = np.roll(r.freq_ev, len(r.freq_ev)/2)#[ 4*len(r.spec)/10:6*len(r.spec)/10]
            y = np.roll( np.abs(r.spec)**2, len(r.freq_ev)/2)#[ 4*len(r.spec)/10:6*len(r.spec)/10]



            p1, = self.axes.plot(x, y, color='r')
            self.axes.hold(True)
            p2, = self.axes.plot(x, si, color='b')
            self.axes.hold(False)

            spec_width = fwhm(x, y)
            print ('bw = ', spec_width)
            print (np.max(np.abs(x)))

            self.axes.legend([p1], ['run ' + str(self.data.show_run_id) + ' ' + str(spec_width) + ' ev'])            
            
            self.draw()
        except Exception as e:
            print ('run not ready...')
            raise(e)


class RunPowerCanvas(MyMplCanvas):
    def __init__(self, *args, **kwargs):
        MyMplCanvas.__init__(self, *args, **kwargs)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        timer.start(100)

    def compute_initial_figure(self):
        self.axes.plot([0,], [0,], color='r')

    def update_figure(self):
        #run_id = np.max(self.data.runs.keys())
        #run_id = 9

        run_id = self.data.show_run_id

        print ('power_z: drawing', run_id)
        try:
            print ('have runs', self.data.runs.keys())
            r = self.data.runs[run_id]
            self.axes.plot(r.z, r.power_z, color='b')
            self.axes.set_title('run ' + str(run_id))
            self.draw()
        except:
            print ('run not ready...')

class ApplicationWindow(QtGui.QMainWindow):
    def __init__(self):
        
        self.data = SimInfo()
        self.data.run_ids = []
        self.data.max_powers = []
        
        QtGui.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("application main window")

        self.file_menu = QtGui.QMenu('&File', self)
        self.file_menu.addAction('&Quit', self.fileQuit, QtCore.Qt.CTRL + QtCore.Qt.Key_Q)
        self.menuBar().addMenu(self.file_menu)

        self.help_menu = QtGui.QMenu('&Help', self)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)

        self.main_widget = QtGui.QWidget(self)

        l = QtGui.QVBoxLayout(self.main_widget)
        dc = RunStatCanvas(self.main_widget, width=5, height=4, dpi=100)
        dc2 = RunPulseCanvas(self.main_widget, width=5, height=4, dpi=100)
        dc3 = RunPowerCanvas(self.main_widget, width=5, height=4, dpi=100)
        dc4 = RunSpecCanvas(self.main_widget, width=5, height=4, dpi=100)
        dc.data = self.data
        dc2.data = self.data
        dc3.data = self.data
        dc4.data = self.data
        l.addWidget(dc)
        l.addWidget(dc2)
        l.addWidget(dc3)
        l.addWidget(dc4)

        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        #self.statusBar().showMessage("All hail matplotlib!", 2000)
                
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_data)
        timer.start(1000)

        

    def fileQuit(self):
        self.close()

    def closeEvent(self, ce):
        self.fileQuit()


    def update_data(self):
        print ('updating...')
        data = pickle.load(open(data_file, 'rb'))
        self.data.run_ids = data.runs.keys()
        self.data.max_powers = map(lambda i: data.runs[i].max_power, self.data.run_ids)
        self.data.runs = data.runs

        self.data.stages = []
        for r_id in self.data.run_ids:
            self.data.stages.append(self.data.runs[r_id].stage)
        
        self.data.stage = np.max(self.data.stages)
        self.data.show_run_id = np.min(self.data.runs.keys())


        for r_id in range(len(self.data.run_ids)):
            run_id = self.data.run_ids[r_id]
            if run_id > self.data.show_run_id and self.data.stages[r_id] == self.data.stage: self.data.show_run_id = run_id

        #self.data.show_run_id = 21
        print (self.data.run_ids)
    

qApp = QtGui.QApplication(sys.argv)

aw = ApplicationWindow()
aw.setWindowTitle("%s" % progname)
aw.show()
sys.exit(qApp.exec_())
#qApp.exec_()
