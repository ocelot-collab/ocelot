'''
definition of a sampling screen
'''

class Screen:
    def __init__(self):
        # position of screen center
        self.x = 0.0
        self.y = 0.0
        self.z = 0.0
        
        # Screen sizes and resolution
        self.size_x = 1.0 # half-size horizontal
        self.size_y = 1.0 # half-size vertical
        self.nx = 1 # points per horizontal direction
        self.ny = 1 # points per vertical direction
        
        # parameters of energy options
        self.start_energy = 100.0
        self.end_energy = 10000.0
        self.num_energy = 1000

        # half of angle aperture.  Angle relative to undulator axis for emittance influence on spectrum
        self.theta_x = 0 # rad
        self.theta_y = 0 # rad