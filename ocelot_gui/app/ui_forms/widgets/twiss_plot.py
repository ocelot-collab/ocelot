from PyQt5 import QtGui
import pyqtgraph as pg


class TwissPlot():

    def __init__(self):
        
        pg.setConfigOptions(antialias=True)
        self.twiss_plot = pg.GraphicsLayoutWidget()

        # Switch to using white background and black foreground
        #pg.setConfigOption('background', 'w')
        #pg.setConfigOption('foreground', 'k')

        self.plot_disp_x = self.twiss_plot.addPlot(row=0, col=0)
        self.plot_disp_x.showGrid(x=True, y=True)

        self.plot_beta = self.twiss_plot.addPlot(row=1, col=0)
        self.plot_beta.showGrid(x=True, y=True)

        self.plot_lattice = self.twiss_plot.addPlot(row=3, col=0)
        self.plot_lattice.showGrid(x=False, y=False)
        #self.plot_lattice.hideAxis('left')
        self.plot_lattice.setMenuEnabled(enableMenu=False)

        self.plot_disp_x.setXLink(self.plot_lattice)
        self.plot_disp_x.addLegend()
        self.plot_beta.setXLink(self.plot_lattice)
        self.plot_beta.addLegend()

        color_blue = QtGui.QColor(0, 0, 255)
        color_red = QtGui.QColor(255, 0, 0)
        color_aqua = QtGui.QColor(0, 255, 255)
        
        pen_blue = pg.mkPen(color_blue, width=2)
        pen_red = pg.mkPen(color_red, width=2)
        pen_aqua = pg.mkPen(color_aqua, width=2)

        self.curv1 = self.plot_disp_x.plot(pen=pen_aqua, name='Dx')
        self.curv2 = self.plot_beta.plot(pen=pen_aqua, name='betaX')
        self.curv3 = self.plot_beta.plot(pen=pen_red, name='betaY')
        