from PyQt5 import QtGui
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget
from functools import partial


class TunePanel():

    def __init__(self, MainWindow, tune_elements, funcs):
        self.mw = MainWindow
        self.tune_elements = tune_elements
        self.change_factor = funcs[0]
        self.change_value = funcs[1]
        self.reset_lattice = funcs[2]
        self.update_lattice = funcs[3]
        self.change_value_up_down = funcs[4]
        self.init_panel()


    def init_panel(self):
        
        self.tune_panel = QtGui.QVBoxLayout()
        self.tune_panel.setAlignment(Qt.AlignTop)

        # title
        label_title = QtGui.QLabel('Tuning elements')
        self.tune_panel.addWidget(label_title)

        # group edit
        box_edit = QtGui.QHBoxLayout()
        box_edit.setAlignment(Qt.AlignCenter)

        plus_button = QtGui.QPushButton("+")
        plus_button.clicked.connect(partial(self.action_value_up_down, 1.0))
        plus_button.setFixedWidth(30)

        self.edit_delta = QtGui.QLineEdit()
        self.edit_delta.setFixedWidth(80)
        self.edit_delta.setAlignment(Qt.AlignCenter)
        self.edit_delta.setText('0.0')

        minus_button = QtGui.QPushButton("-")
        minus_button.clicked.connect(partial(self.action_value_up_down, -1.0))
        minus_button.setFixedWidth(30)

        box_edit.addWidget(plus_button)
        box_edit.addWidget(self.edit_delta)
        box_edit.addWidget(minus_button)
                
        self.tune_panel.addLayout(box_edit)

        widget = QWidget()
        widget.setStyleSheet(""".QWidget {background-color: rgb(255, 255, 255);}""")

        self.box_tune = QtGui.QVBoxLayout()
        self.box_tune.setContentsMargins(0, 0, 0, 0)
        self.box_tune.setSizeConstraint(QtGui.QLayout.SetMinAndMaxSize)
        widget.setLayout(self.box_tune)

        scroll = QtGui.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setMinimumWidth(260)
        scroll.setMinimumHeight(330)
        scroll.setAlignment(Qt.AlignTop)
        scroll.setContentsMargins(0, 0, 0, 0)
        scroll.setWidget(widget)

        self.tune_panel.addWidget(scroll)

        # update / reset buttons
        self.tune_panel.addStretch(1)
        box_buttons = QtGui.QHBoxLayout()

        update_button = QtGui.QPushButton("Update")
        update_button.clicked.connect(self.action_update_lattice)
        reset_button = QtGui.QPushButton("Reset")
        reset_button.clicked.connect(self.action_reset_lattice)

        box_buttons.addWidget(update_button)
        box_buttons.addWidget(reset_button)
        self.tune_panel.addLayout(box_buttons)


    def add_tune_block(self, element_data):

        element = element_data[0]
        var_name = element_data[1]

        eclass = element.__class__.__name__

        if element.id in self.tune_elements or eclass not in self.mw.tunable_elements:
            return

        self.tune_elements[element.id] = {}
        self.tune_elements[element.id]['element'] = element

        # element tuning group
        vbox_tune = QtGui.QVBoxLayout()

        frame_height = 0
        x_button = True
        for type in self.mw.tunable_elements[eclass]:

            frame_height += 70

            self.tune_elements[element.id][type] = {}

            # element tuning group (line 1)
            hbox1 = QtGui.QHBoxLayout()
            hbox1.setAlignment(Qt.AlignLeft)

            hbox_check = QtGui.QCheckBox()
            self.tune_elements[element.id][type]['cb'] = hbox_check

            label_f = QtGui.QLabel('factor:')
            label_f.setAlignment(Qt.AlignRight | Qt.AlignVCenter)

            hbox_factor = QtGui.QLineEdit()
            hbox_factor.setFixedWidth(50)
            hbox_factor.setText('0.01')
            hbox_factor.textChanged.connect(partial(self.change_factor, [element.id, type]))
            self.tune_elements[element.id][type]['factor'] = hbox_factor

            elem_name = QtGui.QLabel(element.id + '.' + type)

            hbox_xbutton = QtGui.QPushButton("x")
            hbox_xbutton.setFixedWidth(20)
            hbox_xbutton.clicked.connect(partial(self.del_tune_block, element.id))

            hbox1.addWidget(hbox_check)
            hbox1.addWidget(elem_name)
            hbox1.addStretch(1)
            hbox1.addWidget(label_f)
            hbox1.addWidget(hbox_factor)
            hbox1.addStretch(1)

            if x_button:
                hbox1.addWidget(hbox_xbutton)
                x_button = False

            vbox_tune.addLayout(hbox1)

            # element tuning group (line 2)
            hbox2 = QtGui.QHBoxLayout()
            hbox2.setAlignment(Qt.AlignLeft)

            label_c = QtGui.QLabel('value:')
            label_c.setAlignment(Qt.AlignRight)

            hbox_curval = QtGui.QDoubleSpinBox()
            hbox_curval.setFixedWidth(90)
            hbox_curval.setRange(-1.0e16, 1.0e16)
            hbox_curval.setSingleStep(float(hbox_factor.text()))
            hbox_curval.setDecimals(6)
            hbox_curval.setAlignment(Qt.AlignRight)
            hbox_curval.setValue(element.__dict__[type])
            hbox_curval.valueChanged.connect(partial(self.change_value, [element.id, type]))
            self.tune_elements[element.id][type]['val'] = hbox_curval

            label_o = QtGui.QLabel('old value:')
            label_o.setAlignment(Qt.AlignRight)
            
            hbox_oldval = QtGui.QLineEdit()
            hbox_oldval.setFixedWidth(75)
            hbox_oldval.setReadOnly(True)
            hbox_oldval.setAlignment(Qt.AlignLeft)
            hbox_oldval.setText(str(round(self.mw.lattice.elements[var_name].__dict__[type], 6)))
            
            hbox2.addWidget(label_c)
            hbox2.addWidget(hbox_curval)
            hbox2.addStretch(1)
            hbox2.addWidget(label_o)
            hbox2.addWidget(hbox_oldval)

            vbox_tune.addLayout(hbox2)

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)
        frame.setFrameShadow(QtGui.QFrame.Plain)
        frame.setLayout(vbox_tune)
        frame.setMaximumHeight(frame_height)

        self.box_tune.addWidget(frame)
        self.tune_elements[element.id]['frame'] = frame


    def del_tune_block(self, id):

        self.box_tune.removeWidget(self.tune_elements[id]['frame'])
        self.tune_elements[id]['frame'].setParent(None)
        del self.tune_elements[id]


    def _get_all_tune_blocks(self):

        ids = []
        for id in self.tune_elements:
            ids.append(id)
        
        return ids


    def action_update_lattice(self):
        
        ids = self._get_all_tune_blocks()
        for id in ids:
            self.del_tune_block(id)
        
        self.update_lattice()

    
    def action_reset_lattice(self):

        ids = self._get_all_tune_blocks()
        for id in ids:
            self.del_tune_block(id)

        self.reset_lattice()


    def action_value_up_down(self, step):
        
        d_val = float(self.edit_delta.text())

        values = []
        ids = self._get_all_tune_blocks()
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
        
        self.change_value_up_down(values)
