
from PyQt5 import QtGui
import pyqtgraph as pg

from ocelot.cpbd.optics import *
from app.forms.widgets.draw_lattice_elements import *


class LatticePlot():

    def __init__(self, MainWindow, ElementPanel):
        self.mw = MainWindow
        self.ep = ElementPanel
        self.scale = 0.2
        self.init_plot()


    def init_plot(self):
        
        pg.setConfigOptions(antialias=True)
        self.lattice_plot = pg.GraphicsLayoutWidget()

        self.plot_lattice = self.lattice_plot.addPlot(row=0, col=0)
        self.plot_lattice.showGrid(x=False, y=False)
        self.plot_lattice.setYRange(-self.scale*1.5, self.scale*1.5)
        self.plot_lattice.hideAxis('left')
        self.plot_lattice.setMenuEnabled(enableMenu=False)

        de = DrawElements(self.mw, self.ep)
        data4 = de.draw_lattice_elements(self.mw.lattice)
        for r in data4:
            self.plot_lattice.addItem(r)
