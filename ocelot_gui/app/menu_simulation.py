"""Menu Simulation description"""

from app.form_simulation_twiss import *


class GUIMenuSimulation():

    def __init__(self, MainWindow):
        self.mw = MainWindow

    
    def __del__(self):
        pass
        
        
    def calc_twiss(self):

        self.mw.clean_central_widget()

        twiss_monitor = TwissMonitor(self.mw)
        self.mw.layout.addWidget(twiss_monitor, 0, 0)
        