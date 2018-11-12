"""Menu Simulation description"""

from app.action_simulation_twiss import *
from app.action_simulation_parameters import *
from app.action_simulation_matching import *


class GUIMenuSimulation():

    def __init__(self, MainWindow):
        self.mw = MainWindow

    
    def calc_twiss(self):

        self.mw.init_central_widget()

        twiss_monitor = TwissMonitor(self.mw)
        self.mw.layout.addWidget(twiss_monitor, 0, 0)


    def calc_params(self):

        #self.mw.init_central_widget()
        
        self.params = CalcParameters(self.mw)
        self.params.show()
        
        
    def matching(self):
        
        self.mw.init_central_widget()

        matching = Matching(self.mw)
        self.mw.layout.addWidget(matching, 0, 0)
