from ocelot.gui import ui_optim_sase
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import sys
import numpy as np
from ocelot.utils.mint.flash1_interface_pydoocs3 import *
from ocelot.utils.mint.machine_setup import *
from desy.flash.lattices.lattice_rf_red import *
#from ocelot.utils.mint.mint import TestInterface
from pyqtgraph.parametertree import Parameter, ParameterTree, ParameterItem, registerParameterType
import pickle
from ocelot.utils.mint.mint import Optimizer, Action
#from ocelot.utils.mint.flash1_interface_pydoocs3 import FLASH1MachineInterface, FLASH1DeviceProperties
#from ocelot.utils.mint.flash_tune_gui import *
from time import sleep


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

    #print("test = ", sub_cor_seq['children'])
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
    return devices


def tree2seq(tree):
    sequence = []
    for p in tree:

        new_seq = {}
        new_seq['name'] = p.name()
        new_seq['method'] = "simplex"
        new_seq['func'] = "max_sase"
        new_seq['maxiter'] = None
        new_seq['order'] = p.opts["value"]
        new_seq['devices'] = []
        new_seq['tol'] = []
        new_seq['values'] = []
        new_seq["type_devs"] = []

        for sub in p:
            if sub.opts["type"] == "action":
                continue
            type_devs, dev = sub.opts["name"].split("/")
            new_seq['devices'].append(dev)
            new_seq['type_devs'].append(type_devs)
            new_seq['values'].append(sub.opts["value"])
            for ch in sub:
                limits = [0, 0]
                I = float(sub.opts["value"])
                tol = ch.opts["value"]

                #print(I, tol)
                limits[0] = I - tol
                limits[1] = I + tol
                """
                lim = [float(s) for s in tol.split(',')]
                if len(lim) == 1:
                    limits[0] = I*(1. - np.sign(I)*lim[0]/100.)
                    limits[1] = I*(1. + np.sign(I)*lim[0]/100.)
                else:
                    limits[0] = lim[0]
                    limits[1] = lim[1]
                """
                new_seq['tol'].append(limits)
        sequence.append(new_seq)
    return sequence


def seq2tree(tree, seq):

    for act in seq:
        new_chld= []
        for i, dev in enumerate(act['devices']):
            tol = act['tol'][i]
            tol = (tol[1] - tol[0])/2.
            #print(tol)
            #tol_str = [str(np.around(x, 3)) for x in tol]
            name = act["type_devs"][i]+'/'+dev
            tmp = {'name': name, 'type': "str", "value": "0",'readonly': True, 'expanded': False, 'children': []}
            tmp['children'].append({'name': 'tol., A', 'type': "float", "value": tol, 'step': 0.01})
            new_chld.append(tmp)

        new_chld.append({'name': 'start Action', 'type': 'action'})

        new_seq = {'name': act['name'], 'type': 'int', 'limits': (0, 20),  'value': act['order'], 'removable': True, 'children': new_chld}
        tree.addChild(new_seq)  # = Parameter.create(name='params', type='group', children=self.sequence)
    return tree


def devices2def_seq(devs, type_devs):

    def_dict = {}
    def_dict['name'] = "action"
    def_dict['order'] = 0
    def_dict['devices'] = devs
    def_dict["type_devs"] = type_devs

    def_dict['tol'] = []
    for dev in devs:
        def_dict['tol'].append([0.1, 0.3])
    def_dict['method'] = "simplex"
    def_dict['func'] = "max_sase"
    def_dict['maxiter'] = None
    def_seq = [def_dict]
    return def_seq


def optim_params(tree):
    opt_params = {}
    for p in tree:
        for sub in p:
            if sub.opts["name"] == "debug":
                opt_params["debug"] = sub.opts["value"]
            elif sub.opts["name"] == "logging":
                opt_params["logging"] = sub.opts["value"]
            elif sub.opts["name"] == 'log file':
                opt_params['log file'] = sub.opts["value"]
            elif sub.opts["name"] == 'timeout':
                opt_params['timeout'] = sub.opts["value"]
    return opt_params


class AThread(QtCore.QThread):
    def __init__(self, optimizer):
        QtCore.QThread.__init__(self)
        self.opt = optimizer

    def __del__(self):
        self.wait()

    def stop(self):
        print("qthread, stop")

    def run(self):
        self.opt.run(self.seq_dict, self.opt_params)
        """
        count = 0
        while count < 50:
            time.sleep(1)
            print ("Increasing")
            count += 1
        """




class OptimApp(QtGui.QMainWindow, ui_optim_sase.Ui_MainWindow):
    def __init__(self, parent=None, params=None, devices=None, optimizer=None):
        super(OptimApp, self).__init__(parent)
        self.setupUi(self, params=params)
        self.opt_thread = optimizer

        try:
            with open("default.seq", 'rb') as f:
                seq = pickle.load(f)

        except:
            seq = []
        #print(seq)
        seq2tree(tree=self.p, seq=seq)

        #self.p.restoreState(state)
        self.debug = False
        self.n_seqs = len(self.p.children())

        #print(self.sequence)
        self.new_seq = []
        #self.hlmi = high_level_mi

        self.start_opt_btm.clicked.connect(self.start_opt)
        self.restore_cur_btn.clicked.connect(self.restore)
        self.stop_opt_btn.clicked.connect(self.force_stop)
        self.setmax_opt_btn.clicked.connect(self.set_max_sase)
        self.save_machine_btn.clicked.connect(self.write_machine)

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
        #orbit= self.hlmi.read_bpms()
        orbit = self.opt_thread.opt.sop.read_bpms()
        self.x = np.array([z[2] for z in orbit])
        self.y = np.array([z[3] for z in orbit])
        self.s = np.array([z[1] for z in orbit])
        self.x_ref = np.zeros(len(self.x))
        self.y_ref = np.zeros(len(self.y))
        self.p.sigTreeStateChanged.connect(self.change)
        #self.tmp = None

        self.opt_thread.finished.connect(self.stop)
        self.work_seq = []
        self.curves_cur = [self.current.plot()]
        self.ndevs = 1
        self.pntr_cur = 0
        self.data = np.zeros((self.ndevs, 100))


    def force_stop(self):
        print("force stop")

        self.opt_thread.opt.isRunning=False
        #print("wasSaved:", self.opt_thread.opt.wasSaved)
        for i in range(300):
            sleep(0.01)
            if self.opt_thread.opt.wasSaved == True:
                break
            #print("inside", i)
        print("*** Params were saved!: ", self.opt_thread.opt.wasSaved)
        self.opt_thread.terminate()
        self.opt_thread.isRunning = False
        self.start_opt_btm.setEnabled(True)
        self.restore_cur_btn.setEnabled(True)

    def stop(self):
        print("stop")
        self.start_opt_btm.setEnabled(True)
        self.restore_cur_btn.setEnabled(True)

    def error_box(self):
        QtGui.QMessageBox.about(self, "Error box", "all Actions have zero order" )

    def start_opt(self):
        self.sequence = np.array(tree2seq(tree=self.p))
        order = []
        for act in self.sequence:
            order.append(act['order'])
        order = np.array(order)
        indx_sort = np.argsort(order)
        indx = np.nonzero(order[indx_sort])[0]
        indx = indx_sort[indx]
        if len(indx) == 0:
            self.error_box()
        else:
            self.work_seq = self.sequence[indx]
            self.optimization()

    def restore(self):
        if not self.opt_thread.opt.isRunning:
            self.opt_thread.opt.restore_currents()


    def optimization(self):

        if not self.opt_thread.opt.isRunning:
            self.create_tree_cur_contr()
            self.start_opt_btm.setEnabled(False)
            self.restore_cur_btn.setEnabled(False)
            opt_params = optim_params(self.p_cntr)

            self.opt_thread.seq_dict=self.work_seq
            self.opt_thread.opt_params=opt_params
            self.opt_thread.start()

            # currents drawing
            self.ndevs = 0
            for act in self.work_seq:
                self.ndevs += len(act["devices"])

            self.current.clear()
            self.data = np.zeros((self.ndevs*2+1, 100))
            self.current.addLegend()
            self.curves_cur = [self.current.plot(pen='r', name='set'), self.current.plot(pen='g', name='RBV')]

            self.pntr_cur = 0


    def save_sequences(self):
        fileName = QtGui.QFileDialog.getSaveFileName(self, 'Save sequences', filter ="txt (*.seq *.)")

        seq = tree2seq(tree=self.p)
        if fileName:
            print( fileName)
            with open(fileName, 'wb') as f:
                #print(self.p.saveState())
                pickle.dump(seq, f)


    def load_sequences(self):
        fileName = QtGui.QFileDialog.getOpenFileName(self, 'Load file', filter ="txt (*.seq *.)")
        if fileName:
            #print(fileName)
            with open(fileName, 'rb') as f:
                seq = pickle.load(f)
                self.p.clearChildren()
                seq2tree(tree=self.p, seq=seq)


    def change(self, param, changes):
        self.sequence = tree2seq(tree=self.p)

        if self.debug: print("tree changes:")
        for param, change, data in changes:
            #print(param,change)
            path = self.p.childPath(param)
            if self.debug:
                if path is not None:
                    childName = '.'.join(path)
                else:
                    childName = param.name()
                print('  parameter: %s'% childName)
                print('  change:    %s'% change)
                print('  data:      %s'% str(data))
                print('  ----------')

            if change == "activated" and not self.opt_thread.opt.isRunning:
                self.work_seq = []
                act_name = path[0]
                for act in self.sequence:
                    if act['name'] == act_name:
                        act['order'] = 1
                        self.work_seq.append(act)
                    else:
                        act['order'] = 0
                self.optimization()

    def create_tree(self):
        self.n_seqs += 1
        self.new_seq[0]["name"] = self.new_seq[0]["name"]+str(self.n_seqs )
        #print("new name = ", self.new_seq)
        seq2tree(tree=self.p, seq=self.new_seq)
        self.t.setParameters(self.p, showTop=False)


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

    def update_tree_currents(self):
        self.sequence = tree2seq(tree=self.p)
        dict_cur = {}
        for p, act in zip(self.p.children(), self.sequence):
            for child, dev in zip(p, act["devices"]):
                if child.opts["type"] == "action":
                    continue
                if dev in dict_cur.keys():
                    current = dict_cur[dev]
                else:
                    current = self.opt_thread.opt.mi.get_value(dev)
                    dict_cur[dev] = current
                child.setValue(str(current))

        self.update_current()

    def create_tree_cur_contr(self):
        devices = []
        for act in self.work_seq:
            for i, devname in enumerate(act["devices"]):
                devices.append(devname)

        self.p_cur_cntr.clearChildren()
        self.p_cur_cntr.addChild({'name': 'Devices', 'type': 'list', 'values': devices, 'value': 0})
        #self.p_cur_cntr.opts["values"] = devices

        self.t_cur_cntr.setParameters(self.p_cur_cntr, showTop=False)
        return devices

    def update_current(self):
        devices = []
        devname_sel = ''
        for p in self.p_cur_cntr:
            devname_sel = p.opts["value"]
            devices = p.opts["values"]
        #print(devices, devname_sel)
        n = 0
        dict_cur = {}
        for act in self.work_seq:
            for i, devname in enumerate(act["devices"]):
                if devname in dict_cur.keys():
                    current_RBS = dict_cur[devname]
                    surrent_set = current_RBS
                else:
                    current_RBS = self.opt_thread.opt.mi.get_value(devname)
                    dict_cur[devname] = current_RBS
                    surrent_set = current_RBS #self.opt_thread.opt.mi.get_value_ps(devname)
                self.data[n, self.pntr_cur] = surrent_set
                self.data[n+1, self.pntr_cur] = current_RBS
                self.data[-1, self.pntr_cur] = self.opt_thread.opt.mi.get_sase()
                n += 2

        self.pntr_cur += 1
        if self.pntr_cur >= self.data.shape[1]:
            tmp = self.data
            self.data = np.empty((self.ndevs*2+1, self.data.shape[1] * 2))
            self.data[:, :tmp.shape[1]] = tmp[:, :]


        for i, name in enumerate(devices):
            #print(name == devname_sel, name, devname_sel)
            if name == devname_sel:
                #print(i)
                self.curves_cur[0].setData(self.data[2*i, :self.pntr_cur])
                self.curves_cur[1].setData(self.data[2*i + 1, :self.pntr_cur])
        #for i, x in enumerate(self.data):
        #    self.curves_cur[i].setData(x[:self.pntr_cur])

    def set_max_sase(self):
        if not self.opt_thread.opt.isRunning:
            sase = self.data[-1,:]
            indx = np.argmax(sase)
            devices = []
            for p in self.p_cur_cntr:
                devices = p.opts["values"]

            for i, name in enumerate(devices):
                I = self.data[2*i+1, indx]
                print(name, "<-- ", I)
                self.opt_thread.opt.mi.set_value(name, I)


    def update_sase(self):

        sase_fast = self.opt_thread.opt.mi.get_sase()
        sase_slow = self.opt_thread.opt.mi.get_sase(detector="gmd_fl1_slow")
        if np.isnan(sase_slow):
            sase_slow = 0.

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

    def update_orbit(self):
        orbit = self.opt_thread.opt.sop.read_bpms()
        self.x = np.array([z[2] for z in orbit])
        self.y = np.array([z[3] for z in orbit])
        #print(self.x_ref)
        x = self.x - self.x_ref # np.random.normal(size=100)
        y = self.y - self.y_ref # np.random.normal(size=100)

        self.curve_orb_x.setData(self.s, x*1000., pen='r', symbol='o', symbolPen='r', symbolBrush=0.5, name='new')

        self.curve_orb_y.setData(self.s, y*1000., pen='r', symbol='o', symbolPen='r', symbolBrush=0.5, name='new')

    def update_blm(self):
        alarm_vals = self.opt_thread.opt.sop.mi.get_alarms()
        self.curve_blm.setData(range(len(alarm_vals)+1), alarm_vals, stepMode=True, fillLevel=0, brush=(0,0,255,150))

    def write_machine(self):
        self.opt_thread.opt.sop.new_tuning()


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
        self.sequence = sequence
        self.p_child.sigTreeStateChanged.connect(self.change)

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
        devices = []
        type_devices = []
        for p in self.p_child:
            for sub in p:
                for ch in sub:
                    if ch.opts["value"] == True:
                        devices.append(ch.opts["name"])
                        type_devices.append(p.opts["name"][:4])

        self.sequence = devices2def_seq(devices, type_devices)

    def create_seq(self):
        #print( "new action = ", self.sequence)
        self.parent().new_seq = self.sequence
        self.parent().create_tree()
        self.restore()

    def save(self):
        self.state = self.p_child.saveState()

    def restore(self):
        self.p_child.restoreState(self.state)


def main():
    mi = FLASH1MachineInterface()
    mi = TestInterface()
    dp = FLASH1DeviceProperties()

    lat = MagneticLattice(lattice)
    sop = SaveOptParams(mi, dp, lat)
    #hlmi = HighLevelInterface(lat, mi, dp)
    #opt = Optimizer(mi, dp)
    opt = Optimizer(mi, dp, sop)
    opt = AThread(optimizer=opt)
    #opt.finished.connect(app.exit)

    app = QtGui.QApplication(sys.argv)


    devices = generate_tree_params(lat)

    form = OptimApp(params=[], devices=devices, optimizer=opt)

    timer = pg.QtCore.QTimer()
    timer.timeout.connect(form.update_sase)
    timer.start(300)

    timer1 = pg.QtCore.QTimer()
    timer1.timeout.connect(form.update_orbit)
    timer1.start(1000)

    timer2 = pg.QtCore.QTimer()
    timer2.timeout.connect(form.update_blm)
    timer2.start(300)

    timer3 = pg.QtCore.QTimer()
    timer3.timeout.connect(form.update_tree_currents)
    timer3.start(2000)

    form.show()
    app.exec_()


if __name__ == '__main__':
    main()