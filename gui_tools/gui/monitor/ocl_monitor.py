"""
S.Tomin PyQtGraph monitor
"""

import sys
from PyQt5.QtWidgets import QCheckBox,QPushButton, QHBoxLayout, QHeaderView, QApplication,QMenu, QWidget, QAction, QTableWidget, QTableWidgetItem, QDoubleSpinBox

from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot
from gui.monitor.ui_ocl_monitor import *
import numpy as np
import pyqtgraph as pg
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

class OclMonitor(QWidget):
    def __init__(self, parent=None):
        super().__init__()
        self.title = 'Twiss'
        self.ui = Ui_Widget()
        self.ui.setupUi(self)
        #make the custom table for middle click
        #self.ui.tableWidget.setParent(None) #remove old table
        #self.ui.tableWidget = CustomTableWidget(self) # make new widget
        #self.ui.gridLayout_2.addWidget(self.ui.tableWidget, 0, 0)
        #self.createTable()
        self.add_plot()
        self.loadStyleSheet()

    def loadStyleSheet(self):
        """ Load in the dark theme style sheet. """
        try:
            self.cssfile = "style.css"
            with open(self.cssfile, "r") as f:
                self.setStyleSheet(f.read())
        except IOError:
            print('No style sheet found!')



    def zoom_signal(self):
        pass
        # s = self.plot1.viewRange()[0][0]
        # s_pos = np.array([q.s_pos for q in self.quads])
        #s_pos = np.array([q.s_pos for q in self.quads]) + self.lat_zi
        #s_up = self.plot1.viewRange()[0][0]
        #s_down = self.plot1.viewRange()[0][1]
        #s_up = s_up if s_up <= s_pos[-1] else s_pos[-1]
        #s_down = s_down if s_down >= s_pos[0] else s_pos[0]
        #indexes = np.arange(np.argwhere(s_pos >= s_up)[0][0], np.argwhere(s_pos <= s_down)[-1][0] + 1)
        #mask = np.ones(len(self.quads), np.bool)
        #mask[indexes] = 0
        #self.quads = np.array(self.quads)
        #[q.ui.set_hide(hide=False) for q in self.quads[indexes]]
        #[q.ui.set_hide(hide=True) for q in self.quads[mask]]

    def update_plot(self, s, bx, by, dx, dy):
        # Line
        s = np.array(s) #+ self.lat_zi
        self.beta_x.setData(x=s, y=bx)
        self.beta_y.setData(x=s, y=by)
        self.plot1.update()
        #self.plot1.setYRange(-5, 200)
        self.plot2.update()
        self.Dx.setData(x=s, y=dx)
        self.Dy.setData(x=s, y=dy)
        self.plot3.update()

    def update_plot_track(self, s, bx, by, E):
        # Line
        s = np.array(s) #+ self.lat_zi
        self.beta_x.setData(x=s, y=bx)
        self.beta_y.setData(x=s, y=by)
        self.plot1.update()
        #self.plot1.setYRange(-5, 200)
        #self.plot2.update()
        self.plot3.removeItem(self.Dx)
        self.plot3.removeItem(self.Dy)

        color = QtGui.QColor(0, 255, 255)
        pen = pg.mkPen(color, width=3)
        self.E = pg.PlotCurveItem(x=[], y=[], pen=pen, name='E', antialias=True)
        self.plot3.addItem(self.E)
        self.E.setData(x=s, y=E)
        self.plot3.update()

    def add_plot(self):
        win = pg.GraphicsLayoutWidget()
        self.plot3 = win.addPlot(row=0, col=0)
        win.ci.layout.setRowMaximumHeight(0, 200)
        self.plot3.showGrid(1, 1, 1)
        self.plot1 = win.addPlot(row=1, col=0)
        self.plot3.setXLink(self.plot1)
        self.plot1.showGrid(1, 1, 1)
        self.plot1.getAxis('left').enableAutoSIPrefix(enable=False)  # stop the auto unit scaling on y axes
        layout = QtGui.QGridLayout()
        layout.setContentsMargins(0,0,0,0)
        self.ui.w_monitor.setLayout(layout)
        layout.addWidget(win, 0, 0)
        self.plot1.setAutoVisible(y=True)
        self.plot1.addLegend()
        color = QtGui.QColor(0, 255, 255)
        pen = pg.mkPen(color, width=3)
        self.beta_x = pg.PlotCurveItem(x=[], y=[], pen=pen, name='beta_x', antialias=True)
        self.plot1.addItem(self.beta_x)
        pen = pg.mkPen(color, width=1)
        self.beta_x_des = pg.PlotCurveItem(x=[], y=[], pen=pen, name='beta_x', antialias=True)
        self.plot1.addItem(self.beta_x_des)
        color = QtGui.QColor(255, 0, 0)
        pen = pg.mkPen(color, width=3)
        self.beta_y = pg.PlotCurveItem(x=[], y=[], pen=pen, name='beta_y', antialias=True)
        self.plot1.addItem(self.beta_y)
        color = QtGui.QColor(255, 0, 0)
        pen = pg.mkPen(color, width=1)
        self.beta_y_des = pg.PlotCurveItem(x=[], y=[], pen=pen, name='beta_y', antialias=True)
        self.plot1.addItem(self.beta_y_des)
        self.plot2 = win.addPlot(row=2, col=0)
        win.ci.layout.setRowMaximumHeight(2, 150)
        self.plot2.setXLink(self.plot1)
        self.plot2.showGrid(1, 1, 1)
        self.plot3.addLegend()
        color = QtGui.QColor(0, 255, 255)
        pen = pg.mkPen(color, width=3)
        self.Dx = pg.PlotCurveItem(x=[], y=[], pen=pen, name='Dx', antialias=True)
        self.plot3.addItem(self.Dx)
        color = QtGui.QColor(255, 0, 0)
        pen = pg.mkPen(color, width=3)
        self.Dy = pg.PlotCurveItem(x=[], y=[], pen=pen, name='Dy', antialias=True)
        self.plot3.addItem(self.Dy)
        self.plot2.sigRangeChanged.connect(self.zoom_signal)

    #@pyqtSlot()
    #def on_click(self):
    #    print("\n")
    #    for currentQTableWidgetItem in self.ui.tableWidget.selectedItems():
    #        print(currentQTableWidgetItem.row(), currentQTableWidgetItem.column(), currentQTableWidgetItem.text())


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = OclMonitor()
    window.show()
    sys.exit(app.exec_())