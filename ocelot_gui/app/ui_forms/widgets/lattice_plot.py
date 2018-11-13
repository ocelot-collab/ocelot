from PyQt5 import QtGui
import pyqtgraph as pg


class LatticePlot():

    def __init__(self):
        
        self.lattice_curves = []
        self.scale = 0.2
        
        pg.setConfigOptions(antialias=True)
        self.lattice_plot = pg.GraphicsLayoutWidget()

        self.plot_lattice = self.lattice_plot.addPlot(row=3, col=0)
        self.plot_lattice.showGrid(x=False, y=False)
        self.plot_lattice.setYRange(-self.scale*1.5, self.scale*1.5)
        self.plot_lattice.hideAxis('left')
        self.plot_lattice.setMenuEnabled(enableMenu=False)
    