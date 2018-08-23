
from PyQt5 import QtGui
import pyqtgraph as pg

from ocelot.cpbd.optics import *
from app.forms.widgets.draw_lattice_elements import *

class TwissPlot():

    def __init__(self, MainWindow, ElementPanel):
        self.mw = MainWindow
        self.ep = ElementPanel
        self.lattice_curves = []
        self.scale = 0.2
        self.init_plot()


    def init_plot(self):
        
        pg.setConfigOptions(antialias=True)
        self.tws_plot = pg.GraphicsLayoutWidget()

        self.plot_disp_x = self.tws_plot.addPlot(row=0, col=0)
        self.plot_disp_x.showGrid(x=True, y=True)

        self.plot_beta = self.tws_plot.addPlot(row=1, col=0)
        self.plot_beta.showGrid(x=True, y=True)

        self.plot_lattice = self.tws_plot.addPlot(row=3, col=0)
        self.plot_lattice.showGrid(x=False, y=False)
        self.plot_lattice.setYRange(-self.scale*1.5, self.scale*1.5)
        #self.plot_lattice.hideAxis('left')
        #self.plot_lattice.setMenuEnabled(enableMenu=False)

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
        self.curv2 = self.plot_beta.plot(pen=pen_blue, name='betaX')
        self.curv3 = self.plot_beta.plot(pen=pen_red, name='betaY')


    def update_plot(self, tws):

        if tws is None:
            return

        for cur in self.lattice_curves:
            self.plot_lattice.removeItem(cur)

        s = [p.s for p in tws]
        beta_x = [p.beta_x for p in tws]
        beta_y = [p.beta_y for p in tws]
        disp_x = [p.Dx for p in tws]
        
        self.curv1.setData(s, disp_x)
        self.curv2.setData(s, beta_x)
        self.curv3.setData(s, beta_y)

        de = DrawElements(self.mw, self.ep)
        data4 = de.draw_lattice_elements(self.mw.work_lattice)
        for r in data4:
            self.lattice_curves.append(r)
            self.plot_lattice.addItem(r)   
