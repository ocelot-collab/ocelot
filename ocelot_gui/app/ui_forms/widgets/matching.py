"""Matching description"""

from ocelot.cpbd.elements import *
from ocelot.cpbd.match import *

from app.ui_forms.widgets.matching_ui import *


class MatchingWindow(QtWidgets.QWidget):

    def __init__(self, lattice, tune_elements):
        super().__init__()
        
        # init user interface
        self.ui = Ui_Form_Matching()
        self.ui.setupUi(self)
        
        self.ui.tws_labels = ['beta_x', 'beta_y', 'alpha_x', 'alpha_y', 'Dx', 'Dxp']

        self.lattice = lattice
        self.tune_elements = tune_elements

        #self.mw = None
        self.matching_done = False
        self.matching_vars = []
        
        # init button actions
        self.ui.btn1.clicked.connect(self.action_reset_results)
        self.ui.btn2.clicked.connect(self.action_matching)
        
        # filling out init twiss parameters
        self.prepare_columbs()

    
    def __del__(self):
        pass


    def update_tws_plot(self):
        pass


    def prepare_columbs(self):

        for label in self.ui.tws_labels:
            val = self.lattice.tws0.__dict__[label]
            if val != 0.0:
                self.ui.__dict__['comboBox_'+label+'_in'].setCurrentIndex(1)
                self.ui.__dict__['lineEdit_'+label+'_in'].setText(str(round(val, 6)))

        if self.lattice.periodic_solution:
            self.ui.periodic_solution.toggle()


    def action_matching(self):
    
        self.ui.old_results.setText('')
        self.ui.new_results.setText('')
        QtWidgets.QApplication.processEvents()
        
        m_start = Monitor(eid='matching_start')
        m_end = Monitor(eid='matching_end')
        
        # collect variables
        string = self.prepare_matching_vars()

        self.ui.old_results.setText(string)
        QtWidgets.QApplication.processEvents()

        if self.matching_vars == []:
            self.ui.new_results.setText('Attention\nSelect matchable variables to do matching')
            QtWidgets.QApplication.processEvents()
            return

        # collect constraints
        constr = self.collect_constraints(monitors=[m_start, m_end])
        
        # prepare lattice
        self.lattice.elements[m_start.id] = m_start
        self.lattice.elements[m_end.id] = m_end

        self.lattice.cell = (m_start,) + self.lattice.cell + (m_end,)
        self.lattice.init_lattice()

        # matching
        string += 'Matching ...'
        self.ui.old_results.setText(string)
        QtWidgets.QApplication.processEvents()
        
        matching_vars = [line['element'] for line in self.matching_vars]
        result = match(self.lattice.lattice, constr, matching_vars, self.lattice.tws0)
            
        # delete matching monitors
        self.lattice.cell = self.lattice.cell[1:-1]
        self.lattice.init_lattice()
        del self.lattice.elements[m_start.id]
        del self.lattice.elements[m_end.id]
        
        # prerare output results
        old_lines = 'Old variable values:\n'
        new_lines = 'New variable values:\n'

        for line in self.matching_vars:
            old_lines += str(line['element'].id) + '.' + str(line['type']) + ' = ' + str(line['init_val']) + '\n'
            new_lines += str(line['element'].id) + '.' + str(line['type']) + ' = ' + str(line['element'].__dict__[line['type']]) + '\n'
        
        self.ui.old_results.setText(old_lines)
        self.ui.new_results.setText(new_lines)
        
        # update tune blocks
        for line in self.matching_vars:
            self.tune_elements[line['elem_id']][line['type']].set_curval(line['element'].__dict__[line['type']], block_signal=True)

        # replot twiss functions
        self.matching_done = True
        self.update_tws_plot()


    def collect_constraints(self, monitors):

        constr = {}

        # check periodic constraint
        if self.ui.periodic_solution.checkState():
            constr['periodic'] = True
        
        # check global constraints and constraints at lattice start and end
        g, g1, g2 = {}, {}, {}
        for label in self.ui.tws_labels:

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
                        self.lattice.tws0.__dict__[label] = val1
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
            constr[monitors[0]] = g1
        if g2 != {}:
            constr[monitors[1]] = g2

        return constr


    def prepare_matching_vars(self):

        string = ''
        self.matching_vars = []

        for elem_id in self.tune_elements:

            element = self.lattice.elements[elem_id]
            elem_class = element.__class__.__name__

            if not element.is_matchable:
                continue

            for type, panel in self.tune_elements[elem_id].items():
                
                if panel.is_selected():

                    if type == self.lattice.matchable_elements[elem_class]:
                        self.matching_vars.append({'element':element, 'elem_id':elem_id, 'type':type, 'init_val':element.__dict__[type]})
                    else:
                        string += 'Variable \'' + type + '\' of element \'' + element.id + '\' can\'t be used for matching\n'

        return string


    def action_reset_results(self):

        if not self.matching_done:
            return
        
        self.matching_done = False

        self.ui.old_results.setText('')
        self.ui.new_results.setText('')
        QtWidgets.QApplication.processEvents()
        
        # reset elements and update tune blocks
        for line in self.matching_vars:
            line['element'].__dict__[line['type']] = line['init_val']
            self.tune_elements[line['elem_id']][line['type']].set_curval(line['init_val'], block_signal=True)
        
        self.lattice.lattice.update_transfer_maps()

        # replot twiss functions
        self.update_tws_plot()
