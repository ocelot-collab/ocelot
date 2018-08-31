"""Menu Edit description"""

from app.forms.edit_field import *


class GUIMenuEdit():

    def __init__(self, MainWindow):
        self.mw = MainWindow


    def edit_lattice(self):
        
        self.mw.init_central_widget()
        
        elements_editor = EditField(self.mw)
        self.mw.layout.addWidget(elements_editor, 0, 0)
