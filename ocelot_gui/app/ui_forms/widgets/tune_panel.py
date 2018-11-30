
from app.ui_forms.widgets.tune_panel_ui import *


class TunePanel(QtWidgets.QWidget):

    def __init__(self):
        super().__init__()

        # init user interface
        self.ui = Ui_Form_TunePanel()
        self.ui.setupUi(self)

        self.ui.factor.editingFinished.connect(self.change_factor)


    def change_factor(self):
        step = self.get_factor()
        self.ui.curval.setSingleStep(step)


    def set_name(self, name):
        self.ui.check.setText(name)

    
    def set_oldval(self, val):
        self.ui.oldval.setText(str(round(val, 6)))

    
    def set_curval(self, val, block_signal=False):
        
        self.ui.curval.blockSignals(block_signal)
        self.ui.curval.setValue(val)
        self.ui.curval.blockSignals(False)

        if round(val, 5) != round(self.get_oldval(), 5):
            self.ui.curval.setStyleSheet("background-color:LightCoral;")
        else:
            self.ui.curval.setStyleSheet("background-color:")


    def get_curval(self):
        return self.ui.curval.value()

    
    def get_oldval(self):
        try:
            val = float(self.ui.oldval.text())
        except:
            val = 0.0
        return val


    def get_factor(self):
        return self.ui.factor.value()


    def is_selected(self):
        return self.ui.check.isChecked()
 