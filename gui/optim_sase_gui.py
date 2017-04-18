import sys
sys.path.append("/home/ttflinac/ocelot/dev/ocelot/gui")
sys.path.append("/home/ttflinac/ocelot/releases/pyqtgraph-0.9.10/build/lib")
sys.path.append("/home/ttflinac/ocelot/dev")

from ocelot.gui import ui_optim_sase
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
from ocelot.mint.flash1_interface import *
#from ocelot.utils.mint.machine_setup import *
from ocelot.gui.flash_tree import *
import pickle
from ocelot.mint.mint import Optimizer, Action
from time import sleep

#pg.setConfigOption('background', 'w')
#pg.setConfigOption('foreground', 'k')

class Logger(object):
    def __init__(self, form, filename="default.log"):
        self.terminal = sys.stderr
        self.log = form.log_lab
        self.form = form
        #self.log.setText("test")
        self.log2file = open(filename, "a")

    def write(self, message):
        self.form.area.moveDock(self.form.logger, 'above', self.form.blm_fig)
        self.terminal.write(message)
        #self.log.append(message)
        self.log.moveCursor(QtGui.QTextCursor.End)
        self.log.insertPlainText( message )
        self.log2file.write(message)
    def flush(self):
        pass



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
        new_seq["pBPM"] = {}
        new_seq["pBPM"]['centers'] = [0, 0, 0, 0]  # centers of the photon beam [tun.x, tun.y, bda.x, bda.y]
        new_seq["pBPM"]['delta'] = [0.5, 0.5, 0.5, 0.5]  # zeros penalty for range = centers +- delta; delta=[tun.x, tun.y, bda.x, bda.y]
        new_seq["pBPM"]['Tunnel'] = False
        new_seq["pBPM"]['BDA'] = False
        for sub in p:
            if sub.opts["type"] == "action":
                continue
            if sub.opts["name"] == "max iter":
                new_seq['maxiter'] = sub.opts["value"]
                continue
            if sub.opts["name"] == "photon BPMs":
                for ch in sub:
                    new_seq["pBPM"][ch.opts["name"]] = ch.opts["value"]
                continue

            indx = sub.opts["name"].find("/")
            type_devs = sub.opts["name"][:indx]
            dev = sub.opts["name"][indx+1:]

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
                new_seq['tol'].append(limits)
        sequence.append(new_seq)
        #print(sequence)
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
            tmp['children'].append({'name': 'lim', 'type': "float", "value": tol, 'step': 0.01, 'tip': 'range: X +/- lim'})
            new_chld.append(tmp)

        if act['maxiter'] == None:
            maxiter = 50
        else:
            maxiter = act['maxiter']
        new_chld.append({'name': 'max iter', 'type': 'int', "value": maxiter})
        try:
            bBDA = act['pBPM']["BDA"]
            bTun = act['pBPM']["Tunnel"]
        except:
            bBDA = False
            bTun = False

        #new_chld.append({'name': 'photon BPMs', 'type': "group", 'expanded': False, 'children':
        #    [{'name': 'BDA', 'type': "bool", "value": bBDA},
        #     {'name': 'Tunnel', 'type': "bool", "value": bTun}]})

        new_chld.append({'name': 'start Action', 'type': 'action'})

        new_seq = {'name': act['name'], 'type': 'int', 'limits': (0, 20),  'value': act['order'], 'renamable': True,
                   'removable': True, 'children': new_chld}
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
    def_dict['pBPM'] = {'BDA': False, "Tunnel":False}
    def_seq = [def_dict]
    return def_seq


def optim_params(tree):
    opt_params = {}
    #print([p for p in tree], tree.opts["children"])
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
            elif sub.opts["name"] == 'detector':
                opt_params['detector'] = sub.opts["value"]
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
        self.opt.run(self.sequence)
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
        self.init_currents = []
        self.sequence = []
        opt_params = optim_params(self.p_cntr)
        self.detector = opt_params["detector"]
        #self.hlmi = high_level_mi

        #self.start_opt_btm.clicked.connect(self.start_opt)
        self.restore_cur_btn.clicked.connect(self.reverting_currents)
        self.stop_opt_btn.clicked.connect(self.force_stop)
        #self.setmax_opt_btn.clicked.connect(self.set_max_sase)
        self.save_machine_btn.clicked.connect(self.write_machine)
        self.clear_disp_btn.clicked.connect(self.clear_sase)
        self.save_seq_btn.clicked.connect(self.save_sequences)
        self.load_seq_btn.clicked.connect(self.load_sequences)


        self.add_seq_btn.clicked.connect(self.create_child)
        #self.restore_btn.clicked.connect(self.create_tree)
        #self.start_btm.clicked.connect(self.update_sase)
        self.window2 = None
        self.devices = devices

        self.sase.setDownsampling(mode='peak')
        self.sase.setClipToView(True)
        self.sase.addLegend()
        self.curve_sase_fast = self.sase.plot(name = "det")
        self.curve_sase_slow = self.sase.plot(name='slow')
        self.curve_penalty = self.sase.plot(name='pen')

        self.curve_blm = self.blm.plot()
        self.blm.setYRange(0, 1)


        self.data_sase = np.empty((3, 100))
        self.ptr_sase = 0

        self.p.sigTreeStateChanged.connect(self.change)
        self.p_cntr.sigTreeStateChanged.connect(self.change_detector)

        self.opt_thread.finished.connect(self.stop)
        self.work_seq = []
        self.curves_cur = [self.current.plot()]
        self.ndevs = 1
        self.pntr_cur = 0
        self.data = np.zeros((self.ndevs, 100))


    def force_stop(self):

        if self.opt_thread.opt.isRunning == False:
            #QtGui.QMessageBox.about(self, "Error box", "Action is nor running" )
            return 0
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

    #def start_opt(self):
    #    self.sequence = np.array(tree2seq(tree=self.p))
    #    order = []
    #    for act in self.sequence:
    #        order.append(act['order'])
    #    order = np.array(order)
    #    indx_sort = np.argsort(order)
    #    indx = np.nonzero(order[indx_sort])[0]
    #    indx = indx_sort[indx]
    #    if len(indx) == 0:
    #        self.error_box()
    #    else:
    #        self.work_seq = self.sequence[indx]
    #        self.optimization()

    def reverting_currents(self):
        if not self.opt_thread.opt.isRunning:

            for act_init_cur in self.init_currents:
                for x in act_init_cur:
                    if len(x) == 2:
                        devname = x[0]
                        current = x[1]
                        print('reverting', devname, '->', current)
                        self.opt_thread.opt.mi.set_value(devname, current)

    def save_init_currents(self):
        init_currents = []
        for act in self.work_seq:
            act_init_cur = []
            for i, devname in enumerate(act['devices']):
                act_init_cur.append([devname, self.opt_thread.opt.mi.get_value(devname)])
                init_currents.append(act_init_cur)
        return init_currents

    def set_limits(self):
        for act in self.work_seq:
            for i, devname in enumerate(act["devices"]):
                limits = act["tol"][i]
                self.opt_thread.opt.dp.set_limits(dev_name=devname, limits=limits)

    def create_sequence(self):

        sequence = []
        for act in self.work_seq:
            func = self.opt_thread.opt.max_sase
            if act["func"] == "min_orbit":
                func = self.opt_thread.opt.min_orbit

            args = [act["devices"], act["method"], {'maxiter': act["maxiter"], 'pBPM': act['pBPM']}]
            #print(args)
            action = Action(func=func, args=args)
            sequence.append(action)
        return sequence

    def set_optimizer_params(self, opt_params):
        self.opt_thread.opt.debug = opt_params["debug"]             #True
        self.opt_thread.opt.logging = opt_params["logging"]         #True
        self.opt_thread.opt.log_file = opt_params['log file']       #'test.log'
        self.opt_thread.opt.timeout = opt_params['timeout']         #1.2
        self.opt_thread.opt.detector = opt_params['detector']       #'gmd_default'
        #print("SET PARAMS", opt_params)

    def optimization(self):

        if not self.opt_thread.opt.isRunning:
            self.create_tree_cur_contr()
            self.start_opt_btm.setEnabled(False)
            self.restore_cur_btn.setEnabled(False)

            print(self.work_seq)
            opt_params = optim_params(self.p_cntr)
            self.set_optimizer_params(opt_params)
            self.set_limits()
            self.init_currents = self.save_init_currents()
            self.opt_thread.sequence = self.create_sequence()

            self.opt_thread.start()

            # currents drawing
            self.ndevs = 0
            for act in self.work_seq:
                self.ndevs += len(act["devices"])

            self.current.clear()
            self.data = np.zeros((self.ndevs*2+1, 100))
            self.current.addLegend()

            self.curves_cur = []
            c = ["r", "g", "b", "y", "w", 'm']
            for i, name in enumerate(self.work_seq[0]["devices"]):
                self.curves_cur.append(self.current.plot(pen=c[i%6],  name=name))
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

    def change_detector(self, param, changes):
        opt_params = optim_params(self.p_cntr)
        self.detector = opt_params["detector"]
        #for param, change, data in changes:

    def create_tree(self):
        self.n_seqs += 1
        self.new_seq[0]["name"] = self.new_seq[0]["name"]+str(self.n_seqs )
        #print("new name = ", self.new_seq)
        seq2tree(tree=self.p, seq=self.new_seq)
        self.t.setParameters(self.p, showTop=False)


    #def orbit2file(self):
    #    fileName = QtGui.QFileDialog.getSaveFileName(self, 'Dialog Title')
    #    if fileName:
    #        print( fileName)

    #def file2orb(self):
    #    fileName = QtGui.QFileDialog.getOpenFileName(self, 'Load file')
    #    if fileName:
    #        print(fileName)

    #def ref_orbit(self):
    #    self.x_ref = self.x
    #    self.y_ref = self.y

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
        #if self.opt_thread.opt.isRunning == True:
        self.update_current()

    def create_tree_cur_contr(self):
        #devices = []
        self.p_cur_cntr.clearChildren()
        for act in self.work_seq:
            for i, devname in enumerate(act["devices"]):
                self.p_cur_cntr.addChild({'name': devname, 'type': 'bool','value': False})
        self.t_cur_cntr.setParameters(self.p_cur_cntr, showTop=False)
        return #devices

    def update_current(self):
        n = 0
        dict_cur = {}

        devices = [p.opts["name"] for p in self.p_cur_cntr]
        devbools= [p.opts["value"] for p in self.p_cur_cntr]
        for i, name in enumerate(devices):

            if name in dict_cur.keys():
                current_RBS = dict_cur[name]
                surrent_set = current_RBS
            else:
                current_RBS = self.opt_thread.opt.mi.get_value(name)
                dict_cur[name] = current_RBS
                surrent_set = current_RBS #self.opt_thread.opt.mi.get_value_ps(devname)

            self.data[n, self.pntr_cur] = surrent_set
            self.data[n+1, self.pntr_cur] = current_RBS
            self.data[-1, self.pntr_cur] = self.opt_thread.opt.mi.get_sase(detector=self.detector)
            n += 2

            if devbools[i]:
                self.curves_cur[i].setData(self.data[2*i, :self.pntr_cur])
            else:
                self.curves_cur[i].clear()
        self.pntr_cur += 1
        if self.pntr_cur >= self.data.shape[1]:
            tmp = self.data
            self.data = np.empty((self.ndevs * 2 + 1, self.data.shape[1] * 2))
            self.data[:, :tmp.shape[1]] = tmp[:, :]

    #def set_max_sase(self):
    #    if not self.opt_thread.opt.isRunning:
    #        sase = self.data[-1,:]
    #        indx = np.argmax(sase)
    #        devices = []
    #        for p in self.p_cur_cntr:
    #            devices = p.opts["values"]
    #
    #        for i, name in enumerate(devices):
    #            I = self.data[2*i+1, indx]
    #            print(name, "<-- ", I)
    #            self.opt_thread.opt.mi.set_value(name, I)


    def update_sase(self):
        sase_fast = self.opt_thread.opt.mi.get_sase(detector=self.detector)
        sase_slow = self.opt_thread.opt.mi.get_sase(detector="gmd_fl1_slow")
        penalty = -self.opt_thread.opt.penalty
        if np.isnan(sase_slow):
            sase_slow = 0.

        self.data_sase[0, self.ptr_sase] = sase_fast #np.random.normal()
        self.data_sase[1, self.ptr_sase] = sase_slow
        self.data_sase[2, self.ptr_sase] = penalty
        self.ptr_sase += 1
        if self.ptr_sase >= self.data_sase.shape[1]:
            tmp1 = self.data_sase
            self.data_sase = np.empty((3, self.data_sase.shape[1] * 2))
            self.data_sase[:, :tmp1.shape[1]] = tmp1[:, :]

        self.curve_sase_fast.setData(self.data_sase[0, :self.ptr_sase])
        self.curve_sase_slow.setData(self.data_sase[1, :self.ptr_sase], width=3,  pen='r', name="sdf")
        self.curve_penalty.setData(self.data_sase[2, :self.ptr_sase], width=3,  pen='g')


    def clear_sase(self):
        self.data_sase = np.empty((3, 100))
        self.ptr_sase = 0


    #def update_orbit(self):
    #    orbit = self.opt_thread.opt.sop.read_bpms()
    #    self.x = np.array([z[2] for z in orbit])
    #    self.y = np.array([z[3] for z in orbit])
    #    #print(self.x_ref)
    #    x = self.x - self.x_ref # np.random.normal(size=100)
    #    y = self.y - self.y_ref # np.random.normal(size=100)
    #
    #    self.curve_orb_x.setData(self.s, x*1000., pen='r', symbol='o', symbolPen='r', symbolBrush=0.5, name='new')
    #
    #    self.curve_orb_y.setData(self.s, y*1000., pen='r', symbol='o', symbolPen='r', symbolBrush=0.5, name='new')

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
                if sub.opts["value"] == True:
                    devices.append(sub.opts["name"])
                    type_devices.append(p.opts["name"][:4])
                for ch in sub:
                    if ch.opts["value"] == True:
                        devices.append(ch.opts["name"])
                        type_devices.append(p.opts["name"][:4])

        self.sequence = devices2def_seq(devices, type_devices)

    def create_seq(self):
        print( "new action = ", self.sequence)
        if len(self.sequence) == 0:
            QtGui.QMessageBox.about(self, "Error box", "Nothing selected")
            return 0
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
    sop = SaveOptParams(mi, dp, lat, dbname='../../data/flash.db')

    #opt = Optimizer(mi, dp)
    opt = Optimizer(mi, dp, sop)
    opt = AThread(optimizer=opt)
    #opt.finished.connect(app.exit)

    app = QtGui.QApplication(sys.argv)


    devices = generate_tree_params(lat)

    form = OptimApp(params=[], devices=devices, optimizer=opt)
    sys.stderr = Logger(form)

    timer = pg.QtCore.QTimer()
    timer.timeout.connect(form.update_sase)
    timer.start(300)

    timer2 = pg.QtCore.QTimer()
    timer2.timeout.connect(form.update_blm)
    timer2.start(300)

    timer3 = pg.QtCore.QTimer()
    timer3.timeout.connect(form.update_tree_currents)
    timer3.start(1000)

    form.show()
    app.exec_()


if __name__ == '__main__':
    main()
