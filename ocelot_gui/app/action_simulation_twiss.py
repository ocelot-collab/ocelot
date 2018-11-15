"""Menu Twiss calculation description"""

from copy import deepcopy
from functools import partial

from ocelot.cpbd.optics import *
from app.parser import *

from app.ui_forms.sim_twiss import *

from app.ui_forms.widgets.twiss_plot import *
from app.ui_forms.widgets.tune_panel_top import *
from app.ui_forms.widgets.tune_panel_bottom import *
from app.ui_forms.widgets.draw_lattice_elements import *


class TwissMonitor(QtWidgets.QWidget):

    def __init__(self, MainWindow):
        super().__init__()
        
        self.tune_elements = {}
        self.mw = MainWindow
        
        # init user interface
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        
        # init twiss plot
        self.twiss_plot = TwissPlot()
        self.ui.tws_plot_widget.setLayout(QtWidgets.QVBoxLayout())
        self.ui.tws_plot_widget.layout().addWidget(self.twiss_plot.twiss_plot)
        
        # init tune panel
        self.ui.widgetArea.setAlignment(QtCore.Qt.AlignTop)
        self.ui.scrollAreaWidgetContents.setLayout(self.ui.widgetArea)
        
        # init button actions
        self.ui.btn1.clicked.connect(self.action_update_lattice)
        self.ui.btn2.clicked.connect(self.action_reset_lattice)
        self.ui.edit_tws_step.editingFinished.connect(self.change_tws_step)
        self.ui.plus_button.clicked.connect(partial(self.action_value_up_down, 1.0))
        self.ui.minus_button.clicked.connect(partial(self.action_value_up_down, -1.0))
        
        # init lattice and calc twiss functions
        self.tws_split = 0.0
        self.tws_split_old = 0.0
        self.init_twiss()


    def init_twiss(self):
        
        self.mw.work_lattice = deepcopy(self.mw.lattice)
        self.mw.work_lattice.init_lattice()
        
        self.calc_twiss()
        self.update_plot()
        

    def calc_twiss(self):

        if self.mw.work_lattice.periodic_solution:
            tws0 = periodic_twiss(Twiss(), lattice_transfer_map(self.mw.work_lattice.lattice, self.mw.work_lattice.tws0.E))
        else:
            tws0 = self.mw.work_lattice.tws0

        if tws0 is not None:
            n = None
            if self.tws_split != 0.0:
                n = int(self.mw.work_lattice.lattice.totalLen / self.tws_split) + 1
            self.tws = twiss(self.mw.work_lattice.lattice, tws0=tws0, nPoints=n)
        else:
            self.tws = None

    
    def update_plot(self):
        
        if self.tws is None:
            return

        for cur in self.twiss_plot.lattice_curves:
            self.twiss_plot.plot_lattice.removeItem(cur)

        s = [p.s for p in self.tws]
        beta_x = [p.beta_x for p in self.tws]
        beta_y = [p.beta_y for p in self.tws]
        disp_x = [p.Dx for p in self.tws]
        
        self.twiss_plot.curv1.setData(s, disp_x)
        self.twiss_plot.curv2.setData(s, beta_x)
        self.twiss_plot.curv3.setData(s, beta_y)

        de = DrawElements(self.mw, self.add_tune_block, self.mw.tunable_elements)
        data4 = de.draw_lattice_elements(self.mw.work_lattice)
        for r in data4:
            self.twiss_plot.lattice_curves.append(r)
            self.twiss_plot.plot_lattice.addItem(r)     
    
    
    def change_tws_step(self):

        try:
            self.tws_split = float(self.ui.edit_tws_step.text())
        except:
            self.ui.edit_tws_step.setText('0.0')
            self.tws_split = 0.0
            
        if self.tws_split_old != self.tws_split:
            self.tws_split_old = self.tws_split
            
            self.calc_twiss()
            self.update_plot()


    def action_reset_lattice(self):
    
        # delete tune blocks
        ids = self.get_all_tune_blocks()
        for id in ids:
            self.del_tune_block(id)

        # reset lattice
        self.mw.work_lattice = deepcopy(self.mw.lattice)
        self.mw.work_lattice.init_lattice()

        self.calc_twiss()
        self.update_plot()


    def action_update_lattice(self):
        
        # delete tune blocks
        ids = self.get_all_tune_blocks()
        for id in ids:
            self.del_tune_block(id)

        # update lattice
        self.mw.lattice = deepcopy(self.mw.work_lattice)

    
    def add_tune_block(self, element_data):
        
        p = Parser()
        elem_id = p.convert_elem_id(element_data.id)
        eclass = element_data.__class__.__name__
        
        if elem_id in self.tune_elements or eclass not in self.mw.tunable_elements:
            return

        self.tune_elements[elem_id] = {}
        self.tune_elements[elem_id]['element'] = element_data
        self.tune_elements[elem_id]['frame'] = []
        
        for i, type in enumerate(self.mw.tunable_elements[eclass]):
 
            # init panel
            form_tune_panel = QtWidgets.QWidget()
            
            if i == 0:
                tp_ui = Ui_Form_TunePanel_Top()
            else:
                tp_ui = Ui_Form_TunePanel_Bottom()
                
            tp_ui.setupUi(form_tune_panel)
                
            # init actions
            tp_ui.curval.valueChanged.connect(partial(self.change_value, [elem_id, type]))
            tp_ui.factor.editingFinished.connect(partial(self.change_factor, [elem_id, type]))
            if i == 0:
                tp_ui.xbutton.clicked.connect(partial(self.del_tune_block, elem_id))
            
            # fields filling out
            tp_ui.check.setText(element_data.id + '.' + type)
            
            tp_ui.curval.setValue(element_data.__dict__[type])
            tp_ui.curval.setSingleStep(float(tp_ui.factor.text()))
            
            tp_ui.oldval.setText(str(round(self.mw.lattice.elements[elem_id].__dict__[type], 6)))

            # adding panel
            self.ui.widgetArea.addWidget(form_tune_panel)
            
            self.tune_elements[elem_id][type] = {}
            self.tune_elements[elem_id][type]['cb'] = tp_ui.check
            self.tune_elements[elem_id][type]['factor'] = tp_ui.factor
            self.tune_elements[elem_id][type]['val'] = tp_ui.curval
            self.tune_elements[elem_id]['frame'].append(form_tune_panel)
    
    
    def del_tune_block(self, id):

        for frame in self.tune_elements[id]['frame']:
            self.ui.widgetArea.removeWidget(frame)
            frame.setParent(None)
        
        del self.tune_elements[id]
        

    def get_all_tune_blocks(self):

        ids = []
        for id in self.tune_elements:
            ids.append(id)
        
        return ids

    
    def change_factor(self, id):
        
        step = float(self.tune_elements[id[0]][id[1]]['factor'].text())
        self.tune_elements[id[0]][id[1]]['val'].setSingleStep(step)

    
    def change_value(self, id, d):

        self.tune_elements[id[0]]['element'].__dict__[id[1]] = float(d)
        self.mw.work_lattice.lattice.update_transfer_maps()
        
        self.calc_twiss()
        self.update_plot()
       
        
    def action_value_up_down(self, step):
    
        d_val = float(self.ui.edit_delta.text())

        values = []
        ids = self.get_all_tune_blocks()
        for id in ids:
            for key in self.tune_elements[id]:
                if isinstance(self.tune_elements[id][key], dict) and self.tune_elements[id][key]['cb'].checkState():
                    factor = float(self.tune_elements[id][key]['factor'].text())
                    value = float(self.tune_elements[id][key]['val'].value())
                    new_value = value + step * factor * d_val
                    values.append([id, key, new_value])

                    self.tune_elements[id][key]['val'].blockSignals(True)
                    self.tune_elements[id][key]['val'].setValue(new_value)
                    self.tune_elements[id][key]['val'].blockSignals(False)

        for val in values:
            self.tune_elements[val[0]]['element'].__dict__[val[1]] = val[2]

        self.mw.work_lattice.lattice.update_transfer_maps()
        
        self.calc_twiss()
        self.update_plot()
    