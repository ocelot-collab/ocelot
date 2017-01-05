
1. opt_objects.py has Device and TargetFunction objects. These object can use MachineInterface
(e.g. lcls_interface.py or xfel_interface.py) to get information from the control system.

simple example (it is not real example):

class Device(object):
    def __init__(self, eid=None):
        self.eid = eid
        self.mi = LCLSMachineInterface()

    def set_value(self, x):
        mi.set_value(self.eid, value)

    def get_value(self):
        return mi.get_value(self.eid)

d1 = Device("QUAD:LTU1:680:BCTRL")
d2 = Device("QUAD:LTU1:660:BCTRL")
d1.get_value()
d2.set_value(12)

Thus, we can avoid the constructions like mi.get_value(QUAD:LTU1:680:BCTRL") inside code. We have to initialize devices somewhere.
And it give us compatibility with different facility.

Other reason, we can use this philosophy to construct more abstract devices or "multiknob" like (this is also not real example :) ):

class Device_BC_R56(object):
    def __init__(self, eid=None):
        self.eid = eid
        self.mi = LCLSMachineInterface()

    def set_value(self, R56):
        ...
        # some calculation
        ...
        mi.set_value("Dipol_1", I1)
        mi.set_value("Dipol_2", I2)
        mi.set_value("Dipol_3", I3)
        mi.set_value("Dipol_4", I4)

    def get_value(self):
        I1 = mi.set_value("Dipol_1")
        I2 = mi.set_value("Dipol_2")
        I3 = mi.set_value("Dipol_3")
        I4 = mi.set_value("Dipol_4")
        ...
        # some calculation
        ...
        return R56



Structure of the Optimizer:

ocelot
    /optimizer

        /GP                    (if I remember correctly I changed GP to search minimum and not maximum)
            bayes_optimization.py
            OnlineGP.py

       /docs
           ….

       /mint
          mint.py              (here is different minimizers (Simplex/GaussProcess/.. ), Optimizer, OptimizerControl and MachineStatus)
          lcls_interface.py
          xfel_interface.py
          opt_objects.py       (here is Device and TargetFunction (or ObjectiveFunction) )
          obj_function.py      (experiment with reloadable from GUI Objective Function)

       /parameters
          default.json        (here is saved all GUI settings (update during closing of the GUI))
          hyperparameters.npy

       /resetpanel            (this is all about table where devices are placed)
          restpanel.py
          restpanelbox.py
          UIrestpanel.py

    generic_optim.py           (Main file)
    gui_main.py                (this is UIOcelotInterface_gen.py "child". I decided to keep in that file all things connected to GUI (almost all things))
    UIOcelotInterface_gen.py   (GUI is generated from UIOcelotInterface_gen.ui using pyui4.bat (for Windows))

Here you can see how to use optimizer with different Minimizers. Almost the same you can find in generic_optim.py self.start_scan().


def test():
    """
    test simplex method
    :return:
    """
    # init Devices
    d1 = TestDevice(eid="d1")
    d2 = TestDevice(eid="d2")
    d3 = TestDevice(eid="d3")

    def get_limits():
        return [-100, 100]
    # reload method get_limits

    d1.get_limits = get_limits
    d2.get_limits = get_limits
    d3.get_limits = get_limits

    devices = [d1, d2, d3]

    # init objective Function
    target = TestTarget()

    # init Optimizer
    opt = Optimizer()
    opt.timeout = 0

    # init GP
    minimizer1 = GaussProcess()
    minimizer1.seed_iter = 3
    minimizer1.max_iter = 300

    # init Simplex
    minimizer2 = Simplex()
    minimizer2.max_iter = 300

    # choose optimization method
    opt.minimizer = minimizer1
    # opt.minimizer = minimizer2

    # create action
    seq = [Action(func=opt.max_target_func, args=[ target, devices])]
    opt.seq = seq

    # start
    opt.eval(seq)

    # or
    # opt.start()


