"""Menu Matching description"""

from copy import deepcopy
from functools import partial
from time import sleep

from ocelot.cpbd.optics import *
from app.parser import *

from app.ui_forms.sim_matching import *

from app.ui_forms.widgets.lattice_plot import *
from app.ui_forms.widgets.matching_panel import *
from app.ui_forms.widgets.draw_lattice_elements import *


class Matching(QtWidgets.QWidget):

    def __init__(self, MainWindow):
        super().__init__()
        
        self.vars_elements = {}
        self.mw = MainWindow
        self.tws_labels = ['beta_x', 'beta_y', 'alpha_x', 'alpha_y', 'Dx', 'Dxp']
        self.matching_dane = False
        
        # init user interface
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        
        # init lattice plot
        self.lattice_plot = LatticePlot()
        self.ui.lattice_plot_widget.setLayout(QtWidgets.QVBoxLayout())
        self.ui.lattice_plot_widget.layout().addWidget(self.lattice_plot.lattice_plot)
        self.update_plot()
        
        # init element matching panel
        self.ui.widgetArea.setAlignment(QtCore.Qt.AlignTop)
        self.ui.scrollAreaWidgetContents.setLayout(self.ui.widgetArea)
        
        # init button actions
        self.ui.btn1.clicked.connect(self.action_apply_results)
        self.ui.btn2.clicked.connect(self.action_matching)
        
        # filling out init twiss parameters
        self.prepare_columbs()
        
    
    def update_plot(self):
        
        self.mw.lattice.init_lattice()
        
        de = DrawElements(self.mw, self.add_tune_block, self.mw.matchable_elements)
        data4 = de.draw_lattice_elements(self.mw.lattice)
        for r in data4:
            self.lattice_plot.plot_lattice.addItem(r)

    
    def prepare_columbs(self):
        
        for label in self.tws_labels:
            val = self.mw.lattice.tws0.__dict__[label]
            if val != 0.0:
                self.ui.__dict__['comboBox_'+label+'_in'].setCurrentIndex(1)
                self.ui.__dict__['lineEdit_'+label+'_in'].setText(str(round(val, 6)))
    
    
    def add_tune_block(self, element_data):

        elem_id = convert_elem_id(element_data.id)
        eclass = element_data.__class__.__name__
        
        if elem_id in self.vars_elements or eclass not in self.mw.matchable_elements:
            return

        self.vars_elements[elem_id] = {}
        self.vars_elements[elem_id]['element'] = element_data
        
        type = self.mw.matchable_elements[eclass]
 
        # init panel
        form_vars_panel = QtWidgets.QWidget()
        tp_ui = Ui_Form_MatchingPanel()   
        tp_ui.setupUi(form_vars_panel)
            
        # init actions
        tp_ui.xbutton.clicked.connect(partial(self.del_tune_block, elem_id))
        
        # fields filling out
        tp_ui.check.setText(element_data.id + '.' + type)
        tp_ui.oldval.setText(str(round(self.mw.lattice.elements[elem_id].__dict__[type], 6)))

        # adding panel
        self.ui.widgetArea.addWidget(form_vars_panel)
        
        self.vars_elements[elem_id]['cb'] = tp_ui.check
        self.vars_elements[elem_id]['frame'] = form_vars_panel
        
    
    def del_tune_block(self, id):

        self.ui.widgetArea.removeWidget(self.vars_elements[id]['frame'])
        self.vars_elements[id]['frame'].setParent(None)
        
        del self.vars_elements[id]

        
    def get_all_tune_blocks(self):

        ids = []
        for id in self.vars_elements:
            ids.append(id)
        
        return ids
        
        
    def action_matching(self):
    
        self.ui.results.setText('')
        QtWidgets.QApplication.processEvents()
        
        self.mw.work_lattice = deepcopy(self.mw.lattice)
        
        m_start = Monitor(eid="matching_start")
        m_end = Monitor(eid="matching_end")
        
        # collect variables
        vars = []
        for elem_id in self.vars_elements:
            if self.vars_elements[elem_id]['cb'].checkState():
                vars.append(self.mw.work_lattice.elements[elem_id])

        if vars == []:
            self.mw.error_window('Attention', 'Select matching variables', 'Warning')
            return
            
        # collect constraints
        constr = {}

        # check periodic constraint
        if self.ui.periodic_solution.checkState():
            constr['periodic'] = True
        
        # check global constraints and constraints at lattice start and end
        g, g1, g2 = {}, {}, {}
        for label in self.tws_labels:

            val = self.ui.__dict__['lineEdit_'+label+'_all'].text()
            val1 = self.ui.__dict__['lineEdit_'+label+'_in'].text()
            val2 = self.ui.__dict__['lineEdit_'+label+'_out'].text()

            cond = str(self.ui.__dict__['comboBox_'+label+'_all'].currentText())
            cond1 = str(self.ui.__dict__['comboBox_'+label+'_in'].currentText())
            cond2 = str(self.ui.__dict__['comboBox_'+label+'_out'].currentText())
            
            if val != '' and cond != '':
                try:
                    val = float(val)
                    g[label] = [cond, val]
                except:
                    pass
                
            if val1 != '' and cond1 != '':
                try:
                    val1 = float(val1)
                    if self.ui.periodic_solution.checkState():
                        g1[label] = val1 if cond1 == '=' else [cond1, val1]
                    elif cond1 == '=':
                        self.mw.work_lattice.tws0.__dict__[label] = val1
                        g1[label] = val1
                except:
                    pass
            
            if val2 != '' and cond2 != '':
                try:
                    val2 = float(val2)
                    g2[label] = val2 if cond2 == '=' else [cond2, val2]
                except:
                    pass

        if g != {}:
            constr['global'] = g
        if g1 != {}:
            constr[m_start] = g1
        if g2 != {}:
            constr[m_end] = g2

        if constr == {}:
            self.mw.error_window('Attention', 'Set matching constraints', 'Warning')
            return
        
        # prepare lattice
        self.mw.work_lattice.elements[m_start.id] = m_start
        self.mw.work_lattice.elements[m_end.id] = m_end

        self.mw.work_lattice.cell = (m_start,) + self.mw.work_lattice.cell + (m_end,)
        self.mw.work_lattice.init_lattice()

        # matching
        self.ui.results.setText('Matching ...')
        QtWidgets.QApplication.processEvents()
        result = match(self.mw.work_lattice.lattice, constr, vars, self.mw.work_lattice.tws0)
            
        # delete matching monitors
        self.mw.work_lattice.cell = self.mw.work_lattice.cell[1:-1]
        del self.mw.work_lattice.elements[m_start.id]
        del self.mw.work_lattice.elements[m_end.id]
        
        lines = 'New variable values:\n'
        for i, v in enumerate(vars):

            eclass = v.__class__.__name__
            type = self.mw.matchable_elements[eclass]
            lines += str(v.id) + '.' + str(type) + ' = ' + str(result[i]) + '\n'

        self.ui.results.setText(lines)
        self.matching_dane = True
        
        
    def action_apply_results(self):

        if not self.matching_dane:
            return
            
        self.ui.results.setText('')
        QtWidgets.QApplication.processEvents()
            
        self.matching_dane = False
        
        # delete vars blocks
        ids = self.get_all_tune_blocks()
        for id in ids:
            self.del_tune_block(id)
        
        # update lattice
        self.mw.lattice = deepcopy(self.mw.work_lattice)
        self.ui.results.setText('Lattice updated')
        