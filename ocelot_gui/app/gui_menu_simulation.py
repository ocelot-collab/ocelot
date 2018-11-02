"""Menu Simulation description"""

from app.forms.twiss_monitor import *
from app.forms.matching import *


class GUIMenuSimulation():

    def __init__(self, MainWindow):
        self.mw = MainWindow

    
    def calc_parameters(self):

        self.mw.init_central_widget()

        twiss_monitor = TwissMonitor(self.mw)
        self.mw.layout.addWidget(twiss_monitor, 0, 0)


    def matching(self):
        
        self.mw.init_central_widget()

        match = Matching(self.mw)
        self.mw.layout.addWidget(match, 0, 0)
