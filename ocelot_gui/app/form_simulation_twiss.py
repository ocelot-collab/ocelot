"""Menu Twiss calculation description"""

from copy import deepcopy

from app.parser import *

from app.ui_forms.sim_twiss import *

from app.functions.calc_tws import *
from app.functions.calc_lattice_element_pics import *

from app.ui_forms.widgets.twiss_plot import *
from app.ui_forms.widgets.tune_panel import *
from app.ui_forms.widgets.ebeam_table import *
from app.ui_forms.widgets.matching import *


class TwissMonitor(QtWidgets.QWidget):

    def __init__(self, MainWindow):
        super().__init__()
        
        # init user interface
        self.ui = Ui_Form_Twiss()
        self.ui.setupUi(self)

        self.mw = MainWindow

        # lattice and other parameters containers
        self.lattice = None
        self.tws = None
        self.elements_pics = None

        # elements are chosen to use for tuning and matching
        self.tune_elements = {}

        # list of opend chiled windows
        self.child_windows = {}
        self.child_windows['params'] = None
        self.child_windows['matching'] = None
               
        # init twiss plot
        self.twiss_plot = TwissPlot()
        self.ui.tws_plot_widget.setLayout(QtWidgets.QVBoxLayout())
        self.ui.tws_plot_widget.layout().addWidget(self.twiss_plot.twiss_plot)
        
        # init tune panel widgets list
        self.ui.widgetArea.setAlignment(QtCore.Qt.AlignTop)
        self.ui.scrollAreaWidgetContents.setLayout(self.ui.widgetArea)
        
        # init button actions
        self.ui.btn1.clicked.connect(self.action_update_lattice)
        self.ui.btn2.clicked.connect(self.action_reset_lattice)
        
        self.ui.plus_button.step = 1.0
        self.ui.minus_button.step = -1.0
        self.ui.plus_button.clicked.connect(self.action_value_up_down)
        self.ui.minus_button.clicked.connect(self.action_value_up_down)

        self.ui.calc_params.clicked.connect(self.action_open_params_window)
        self.ui.calc_matching.clicked.connect(self.action_open_matching_window)
        self.ui.edit_tws_step.editingFinished.connect(self.change_tws_step)
        
        # init lattice and calc twiss functions
        self.init_twiss()


    def deleteLater(self):

        # this is just for fix then calling TwissMonitor.__del__
        # without deleting elem_pic.add_tune_block function TwissMonitor.__del__ don't call
        for elem_pic in self.elements_pics.elements_pics:
            elem_pic.statusBar = None
            elem_pic.add_tune_block = None

        # close all opened child windows
        for type, window in self.child_windows.items():
            if window:
                window.close()
        
        # delete winget TwissMonitor
        super().deleteLater()


    def __del__(self):
        pass


    def init_twiss(self):
        
        self.lattice = deepcopy(self.mw.lattice)
        self.lattice.init_lattice()
        
        self.twiss_plot.plot_lattice.setYRange(-1.5, 1.5, padding=0)
        self.twiss_plot.plot_lattice.setXRange(0.0, self.lattice.lattice.totalLen, padding=0)

        self.tws = TwissParameters(self.lattice)
        self.elements_pics = LatticeElementsPics(self.lattice)

        self.update_tws_plot()


    def action_open_params_window(self):

        self.child_windows['params'] = EbeamWindow(self.lattice, self.tws)
        self.child_windows['params'].closeEvent = self.action_close_params_window
        self.child_windows['params'].show()


    def action_close_params_window(self, event=None):

        self.child_windows['params'].close()
        self.child_windows['params'] = None


    def action_open_matching_window(self):

        self.child_windows['matching'] = MatchingWindow(self.lattice, self.tune_elements)
        self.child_windows['matching'].closeEvent = self.action_close_matching_window
        self.child_windows['matching'].update_tws_plot = self.update_tws_plot
        self.child_windows['matching'].show()

    
    def action_close_matching_window(self, event=None):

        self.child_windows['matching'].close()
        self.child_windows['matching'] = None

    
    def update_tws_plot(self):

        # calc elements pics
        self.elements_pics.calc_lattice_elements_pics()

        # rescale lattice elements plot
        axX = self.twiss_plot.plot_lattice.getAxis('bottom')
        axY = self.twiss_plot.plot_lattice.getAxis('left')

        self.twiss_plot.plot_lattice.setYRange(axY.range[0], axY.range[1], padding=0)
        self.twiss_plot.plot_lattice.setXRange(axX.range[0], axX.range[1], padding=0)

        # redraw elements pics
        self.twiss_plot.plot_lattice.clear()

        for elem_pic in self.elements_pics.elements_pics:
            
            # add hoverMouseEvent functionallity
            elem_pic.statusBar = self.mw.statusBar
            elem_pic.add_tune_block = self.add_tune_block

            self.twiss_plot.plot_lattice.addItem(elem_pic)

        # calc twiss functions
        self.tws.calc_twiss()

        if self.tws.tws is None:
            self.twiss_plot.curv1.setData([], [])
            self.twiss_plot.curv2.setData([], [])
            self.twiss_plot.curv3.setData([], [])
            return

        # redraw twiss functions
        s = [p.s for p in self.tws.tws]
        beta_x = [p.beta_x for p in self.tws.tws]
        beta_y = [p.beta_y for p in self.tws.tws]
        disp_x = [p.Dx for p in self.tws.tws]
        
        self.twiss_plot.curv1.setData(s, disp_x)
        self.twiss_plot.curv2.setData(s, beta_x)
        self.twiss_plot.curv3.setData(s, beta_y)
    
    
    def change_tws_step(self):

        new_step = self.ui.edit_tws_step.value()
            
        if new_step != self.tws.tws_step:
            self.tws.tws_step = new_step
            self.update_tws_plot()


    def action_reset_lattice(self):

        # update tune blocks
        for elem_id in self.tune_elements:
            for type, panel in self.tune_elements[elem_id].items():

                value = self.mw.lattice.elements[elem_id].__dict__[type]
                panel.set_oldval(value)
                panel.set_curval(value, True)

        # reset lattice
        for elem_id, element in self.lattice.elements.items():
            if element.is_tuneable:
                
                for type in element.tune_params:
                    element.__dict__[type] = self.mw.lattice.elements[elem_id].__dict__[type]

        # reset Twiss parameters
        self.lattice.tws0 = deepcopy(self.mw.lattice.tws0)
        self.lattice.lattice.update_transfer_maps()

        self.update_tws_plot()
        

    def action_update_lattice(self):

        # update tune blocks
        for elem_id in self.tune_elements:
            for type, panel in self.tune_elements[elem_id].items():

                value = self.lattice.elements[elem_id].__dict__[type]
                panel.set_oldval(value)
                panel.set_curval(value, True)
                
        # update lattice
        for elem_id, element in self.lattice.elements.items():
            if element.is_tuneable:
                
                for type in element.tune_params:
                    self.mw.lattice.elements[elem_id].__dict__[type] = element.__dict__[type]

        self.mw.lattice.lattice.update_transfer_maps()

    
    def add_tune_block(self, element):

        p = Parser()
        elem_id = p.convert_elem_id(element.id)
        
        if not element.is_tuneable:
            return
        
        if elem_id not in self.tune_elements:
            self.tune_elements[elem_id] = {}
        
        for i, type in enumerate(element.tune_params):

            if type in self.tune_elements[elem_id]:
                continue

            tune_panel = TunePanel()
            tune_panel.set_name(element.id + '.' + type)
            tune_panel.set_oldval(self.mw.lattice.elements[elem_id].__dict__[type])
            tune_panel.set_curval(element.__dict__[type], block_signal=True)

            tune_panel.ui.curval.params = [elem_id, type]
            tune_panel.ui.curval.valueChanged.connect(self.change_value)

            tune_panel.ui.xbutton.params = [elem_id, type]
            tune_panel.ui.xbutton.clicked.connect(self.del_tune_block)

            self.tune_elements[elem_id][type] = tune_panel
            self.ui.widgetArea.addWidget(tune_panel)

    
    def del_tune_block(self, elem_params=None):

        if not elem_params:
            target = self.sender()
            elem_id = target.params[0]
            type = target.params[1]
        else:
            elem_id = elem_params[0]
            type = elem_params[1]

        # delet widget from the main form
        tune_panel = self.tune_elements[elem_id][type]
        self.ui.widgetArea.removeWidget(tune_panel)
        tune_panel.setParent(None)

        # delete element from array
        del self.tune_elements[elem_id][type]
        if self.tune_elements[elem_id] == {}:
            del self.tune_elements[elem_id]
    
    
    def change_value(self, d, elem_params=None):

        if not elem_params:
            target = self.sender()
            elem_id = target.params[0]
            type = target.params[1]
        else:
            elem_id = elem_params[0]
            type = elem_params[1]

        self.lattice.elements[elem_id].__dict__[type] = d
        self.lattice.lattice.update_transfer_maps()

        # just only for update background color while mouse sliding
        self.tune_elements[elem_id][type].set_curval(d, block_signal=True)
        
        self.update_tws_plot()
       
        
    def action_value_up_down(self, step=None):
    
        if not step:
            target = self.sender()
            step = target.step

        d_val = self.ui.edit_delta.value()

        for elem_id in self.tune_elements:
            for type, panel in self.tune_elements[elem_id].items():
                
                if panel.is_selected():

                    factor = panel.get_factor()
                    value = panel.get_curval()
                    new_value = value + step * factor * d_val

                    panel.set_curval(new_value, block_signal=True)
                    self.lattice.elements[elem_id].__dict__[type] = new_value
        
        self.lattice.lattice.update_transfer_maps()
        
        self.update_tws_plot()
