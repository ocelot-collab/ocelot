from ocelot.gui import ui_optim_sase
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import sys
import numpy as np
from ocelot.utils.mint.flash1_interface_pydoocs3 import *
from ocelot.utils.mint.machine_setup import *
from desy.flash.lattices.lattice_rf_red import *
from ocelot.utils.mint.mint import TestInterface
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
import pickle


def generate_tree_params(lat):
    devices = []
    cors = {'name': 'correctors', 'type': 'group', 'children':[]}
    #cors['children'] = cor_places
    quads = {'name': 'quadrupoles', 'type': 'group', 'children':[]}
    #quad_chld = []
    cav = {}
    section = {'STARTUBC2':"BC2", 'STARTACC2':"ACC23", 'STARTUBC3':"BC3", 'STARTACC4':"ACC4-7", 'ENDACC7':"DOGLEG",
               'STARTSMATCH1':"MATCH", 'STARTUND':"UND"}

    name_seq = "GUN-ACC39"
    sub_cor_seq = {'name': name_seq, 'type': 'group','expanded': False, 'children':[]}
    sub_cor_chld = []

    sub_quad_seq = {'name': name_seq, 'type': 'group','expanded': False, 'children':[]}
    sub_quad_chld = []

    print("test = ", sub_cor_seq['children'])
    for elem in lat.sequence:
        if elem.__class__ == Marker and elem.id in section.keys():

            name_seq = section[elem.id]
            #correctors
            sub_cor_seq['children']= sub_cor_chld
            cors['children'].append(sub_cor_seq)
            sub_cor_seq = {'name': name_seq, 'type': 'group','expanded': False,'children':[]}
            sub_cor_chld = []

            #quadrupoles
            sub_quad_seq['children'] = sub_quad_chld
            quads['children'].append(sub_quad_seq)
            sub_quad_seq = {'name': name_seq, 'type': 'group','expanded': False,'children':[]}
            sub_quad_chld = []

        #if elem.id == 'STARTUBC2':
        #    break

        if elem.__class__ in [Hcor, Vcor]:
            tmp = {}
            tmp["name"] = elem.id
            tmp["type"] = "bool"
            tmp["value"] = False
            sub_cor_chld.append(tmp)
        elif elem.__class__ in [Quadrupole]:
            if elem.id in [p['name'] for p in sub_quad_chld]:
                continue
            tmp = {}
            tmp["name"] = elem.id
            tmp["type"] = "bool"
            tmp["value"] = False
            sub_quad_chld.append(tmp)

    sub_cor_seq['children'] = sub_cor_chld
    cors['children'].append(sub_cor_seq)

    sub_quad_seq['children'] = sub_quad_chld
    quads['children'].append(sub_quad_seq)
    devices.append(cors)
    devices.append(quads)
    #print("cors = ", cors)
    #print("deevice = ", devices)
    return devices


class ExampleApp(QtGui.QMainWindow, ui_optim_sase.Ui_MainWindow):
    def __init__(self, high_level_mi, parent=None, params=None, devices=None):
        super(ExampleApp, self).__init__(parent)
        self.setupUi(self, params=params)

        with open("test.seq", 'rb') as f:
            state = pickle.load(f)
            self.p.restoreState(state)
        #print(state)
        self.n_seqs = len(self.p.children())

        #print(self.sequence)
        self.new_seq = []
        self.hlmi = high_level_mi

        self.save_seq_btn.clicked.connect(self.save_sequences)
        self.load_seq_btn.clicked.connect(self.load_sequences)


        self.save_btn.clicked.connect(self.orbit2file)
        self.load_btn.clicked.connect(self.file2orb)
        self.ref_btm.clicked.connect(self.ref_orbit)

        self.add_seq_btn.clicked.connect(self.create_child)
        #self.restore_btn.clicked.connect(self.create_tree)
        #self.start_btm.clicked.connect(self.update_sase)
        self.window2 = None
        self.devices = devices

        self.sase.setDownsampling(mode='peak')
        self.sase.setClipToView(True)
        self.curve_sase_fast = self.sase.plot()
        self.curve_sase_slow = self.sase.plot()

        self.curve_blm = self.blm.plot()
        self.blm.setYRange(0, 1)
        #self.orbit.setDownsampling(mode='peak')
        #self.orbit.setClipToView(True)
        self.p_x = self.orbit.addPlot(title="X", row=0, col=0)
        self.p_y = self.orbit.addPlot(title="Y", row=1, col=0)
        self.curve_orb_x = self.p_x.plot()
        self.curve_orb_x_ref = self.p_x.plot()
        self.curve_orb_y = self.p_y.plot()
        self.curve_orb_y_ref = self.p_y.plot()
        self.data_fast = np.empty(100)
        self.ptr1 = 0
        self.data_slow = np.empty(100)
        self.ptr2 = 0
        self.x, self.y = self.hlmi.read_bpms()
        self.x_ref = np.zeros(len(self.x))
        self.y_ref = np.zeros(len(self.y))
        self.p.sigTreeStateChanged.connect(self.change)
        #self.tmp = None


    def save_sequences(self):
        fileName = QtGui.QFileDialog.getSaveFileName(self, 'Save sequences')
        state = self.p.saveState()
        if len(self.p.children()) == 0:
            state = {'type': 'group', 'expanded': True, 'readonly': False, 'value': None, 'renamable': False, 'removable': False, 'enabled': True, 'visible': True, 'children': []}
        if fileName:
            print( fileName)
            with open(fileName, 'wb') as f:
                print(self.p.saveState())
                pickle.dump(state, f)


    def load_sequences(self):
        fileName = QtGui.QFileDialog.getOpenFileName(self, 'Load file')
        if fileName:
            print(fileName)
            with open(fileName, 'rb') as f:
                state = pickle.load(f)
                print(state)
                self.p.restoreState(state)# = Parameter.create(name='params', type='group', children=params)


    def change(self, param, changes):
        self.sequence = []
        for p in self.p:
            new_seq = {}
            print(p.name())
            new_seq['name'] = p.name()
            for sub in p:
                if sub.opts["name"] == "order":
                    new_seq['order'] = sub.opts["value"]
                new_seq['devices'] = []
                new_seq['tol'] = []
                for ch in sub:
                    #print(ch.opts["name"])
                    new_seq['devices'].append(ch.opts["name"])
                    new_seq['tol'].append(ch.opts["value"])

            self.sequence.append(new_seq)
        print(self.sequence)


    def create_tree(self):
        self.n_seqs += 1
        print(" len p.children", len(self.p.children()))
        new_chld= []


        for name in self.new_seq:
            tmp = {}
            tmp['name'] = name
            tmp['type'] = 'float'
            tmp['value'] = 10
            new_chld.append(tmp)
        #new_seq['children'] = new_chld

        new_seq = {'name': 'seq'+str(self.n_seqs), 'type': 'group', 'removable': True, 'children': [
                    {'name': 'order', 'type': 'int', 'value': 0},
                        {'name': 'devices', 'type': 'group', 'children': new_chld},
        ]}

        #self.sequence.append(new_seq)
        #print (self.sequence)
        #print(params)
        self.p.addChild(new_seq)# = Parameter.create(name='params', type='group', children=self.sequence)
        self.t.setParameters(self.p, showTop=False)
        #print(self.tmp)



    def orbit2file(self):
        fileName = QtGui.QFileDialog.getSaveFileName(self, 'Dialog Title')
        if fileName:
            print( fileName)

    def file2orb(self):
        fileName = QtGui.QFileDialog.getOpenFileName(self, 'Load file')
        if fileName:
            print(fileName)

    def ref_orbit(self):
        self.x_ref = self.x
        self.y_ref = self.y

    def update_sase(self):
        #self.hlmi.read_sase()
        sase_fast, sase_slow = self.hlmi.read_sase()

        self.data_fast[self.ptr1] = sase_fast #np.random.normal()
        self.ptr1 += 1

        self.data_slow[self.ptr2] = sase_slow #np.random.normal()
        self.ptr2 += 1
        if self.ptr1 >= self.data_fast.shape[0]:
            tmp1 = self.data_fast
            tmp2 = self.data_slow
            self.data_fast = np.empty(self.data_fast.shape[0] * 2)
            self.data_fast[:tmp1.shape[0]] = tmp1

            self.data_slow = np.empty(self.data_slow.shape[0] * 2)
            self.data_slow[:tmp2.shape[0]] = tmp2
        self.curve_sase_fast.setData(self.data_fast[:self.ptr1])

        self.curve_sase_slow.setData(self.data_slow[:self.ptr2], pen='r')
        #self.sase.plot(np.random.normal(size=100))

    def update_orbit(self):
        self.x, self.y = self.hlmi.read_bpms()
        x = self.x - self.x_ref # np.random.normal(size=100)
        #x_ref = np.random.normal(size=100)*0.0001
        y = self.y - self.y_ref # np.random.normal(size=100)
        #y_ref = np.random.normal(size=100)*0.0001
        self.curve_orb_x.setData(x*1000., pen='r', symbol='o', symbolPen='r', symbolBrush=0.5, name='new')
        #self.curve_orb_x_ref.setData(self.x_ref*1000., pen='r', symbol='o', symbolPen='g', symbolBrush=0.5, name='ref')

        self.curve_orb_y.setData(y*1000., pen='r', symbol='o', symbolPen='r', symbolBrush=0.5, name='new')
        #self.curve_orb_y_ref.setData(self.y_ref*1000., pen='r', symbol='o', symbolPen='g', symbolBrush=0.5, name='ref')

    def update_blm(self):
        alarm_vals = self.hlmi.mi.get_alarms()
        print(alarm_vals)
        self.curve_blm.setData(range(len(alarm_vals)+1), alarm_vals, stepMode=True, fillLevel=0, brush=(0,0,255,150))

    def create_child(self):
        # here put the code that creates the new window and shows it.
        if self.window2 is None:
            self.window2 = Form2(parent=self, params=self.devices)
        self.window2.show()

class Form2(QtGui.QMainWindow, ui_optim_sase.Ui_ChildWindow):
    def __init__(self, parent=None, params=None, sequence=[]):
        super(Form2, self).__init__(parent)
        self.setupUi(self, params=params)
        self.state = self.p_child.saveState()
        self.saveBtn.clicked.connect(self.create_seq)
        #print("self.x =", self.x)
        self.sequence = sequence
        self.p_child.sigTreeStateChanged.connect(self.change)
        #self.p_child.param('Save/Restore functionality', 'Save State').sigActivated.connect(self.save)
        #self.p_child.param('Save/Restore functionality', 'Restore State').sigActivated.connect(self.restore)

    def change(self, param, changes):
        """
        print("tree changes:")
        for param, change, data in changes:
            path = self.p_child.childPath(param)
            if path is not None:
                childName = '.'.join(path)
            else:
                childName = param.name()
            print('  parameter: %s'% childName)
            print('  change:    %s'% change)
            print('  data:      %s'% str(data))
            print('  ----------')
            #if data == True:
        """
        self.sequence = []
        for p in self.p_child:
            print(p)
            for sub in p:
                for ch in sub:
                    #print(ch.opts)
                    if ch.opts["value"] == True:
                        self.sequence.append(ch.opts["name"])
        #print(self.sequence)

    def create_seq(self):
        self.parent().new_seq = self.sequence
        self.parent().create_tree()
        self.restore()


    def save(self):
        #global state
        self.state = self.p_child.saveState()
        #return self
        #print(state)

    def restore(self):
        #global state
        #add = p['Save/Restore functionality', 'Restore State', 'Add missing items']
        #rem = p['Save/Restore functionality', 'Restore State', 'Remove extra items']
        self.p_child.restoreState(self.state)
        #self.p_child.restoreState(self.state, addChildren=add, removeChildren=rem)


def main():
    mi = FLASH1MachineInterface()
    mi = TestInterface()
    dp = FLASH1DeviceProperties()
    lat = MagneticLattice(lattice)
    hlmi = HighLevelInterface(lat, mi, dp)
    app = QtGui.QApplication(sys.argv)
    params = [

    ]



    devices = generate_tree_params(lat)

    form = ExampleApp(hlmi, params=params, devices=devices)
    timer = pg.QtCore.QTimer()
    timer.timeout.connect(form.update_sase)
    timer.start(100)

    timer1 = pg.QtCore.QTimer()
    timer1.timeout.connect(form.update_orbit)
    timer1.start(1000)

    timer2 = pg.QtCore.QTimer()
    timer2.timeout.connect(form.update_blm)
    timer2.start(100)

    form.show()
    app.exec_()


if __name__ == '__main__':
    main()