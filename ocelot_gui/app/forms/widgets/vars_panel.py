from PyQt5 import QtGui
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget
from functools import partial


class VarsPanel():

    def __init__(self, MainWindow, vars_elements):
        self.mw = MainWindow
        self.vars_elements = vars_elements
        self.init_panel()


    def init_panel(self):
        
        self.vars_panel = QtGui.QVBoxLayout()
        self.vars_panel.setAlignment(Qt.AlignTop)

        # title
        label_title = QtGui.QLabel('Matching variables')
        self.vars_panel.addWidget(label_title)

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
        scroll.setMinimumHeight(450)
        scroll.setAlignment(Qt.AlignTop)
        scroll.setContentsMargins(0, 0, 0, 0)
        scroll.setWidget(widget)

        self.vars_panel.addWidget(scroll)


    def add_tune_block(self, element_data):

        element = element_data[0]
        var_name = element_data[1]

        eclass = element.__class__.__name__

        if element.id in self.vars_elements or eclass not in self.mw.matchable_elements:
            return

        self.vars_elements[element.id] = {}
        self.vars_elements[element.id]['element'] = var_name

        # element tuning group
        vbox_tune = QtGui.QVBoxLayout()
        frame_height = 40

        type = self.mw.matchable_elements[eclass][0]

        # element tuning group (line 1)
        hbox1 = QtGui.QHBoxLayout()
        hbox1.setAlignment(Qt.AlignLeft)

        label_v = QtGui.QLabel('value:')
        label_v.setAlignment(Qt.AlignRight)
        
        hbox_val = QtGui.QLineEdit()
        hbox_val.setFixedWidth(75)
        hbox_val.setReadOnly(True)
        hbox_val.setAlignment(Qt.AlignLeft)
        hbox_val.setText(str(round(self.mw.lattice.elements[var_name].__dict__[type], 6)))

        elem_name = QtGui.QLabel(element.id + '.' + type)

        hbox_xbutton = QtGui.QPushButton("x")
        hbox_xbutton.setFixedWidth(20)
        hbox_xbutton.clicked.connect(partial(self.del_tune_block, element.id))

        hbox1.addWidget(elem_name)
        hbox1.addStretch(1)
        hbox1.addWidget(label_v)
        hbox1.addWidget(hbox_val)
        hbox1.addStretch(1)

        hbox1.addWidget(hbox_xbutton)

        vbox_tune.addLayout(hbox1)

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)
        frame.setFrameShadow(QtGui.QFrame.Plain)
        frame.setLayout(vbox_tune)
        frame.setMaximumHeight(frame_height)

        self.box_tune.addWidget(frame)
        self.vars_elements[element.id]['frame'] = frame


    def del_tune_block(self, id):

        self.box_tune.removeWidget(self.vars_elements[id]['frame'])
        self.vars_elements[id]['frame'].setParent(None)
        del self.vars_elements[id]
