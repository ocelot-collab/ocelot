__author__ = 'Sergey Tomin'

import cPickle
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

from ocelot.cpbd.orbit_correction import *


class SuperOrbit(QObject):
    def __init__(self):
        QObject.__init__(self)

    def open_file(self, path):
        exec( open(path ))
        self.errors = errors
        self.lattice = MagneticLattice(ring)
        self.orbit = Orbit()
        self.lat_errors = create_copy(self.lattice, nsuperperiods=1)
        self.orbit.lattice_analysis(self.lat_errors)


    def seed_errors(self):

        self.lat_errors,a = errors_seed(self.lat_errors, self.errors)
        #self.orbit = Orbit()
        self.orbit.lattice_analysis(self.lat_errors)
        print "len = ", len(self.orbit.bpms)
        self.emit(SIGNAL("resp_off()"))

    def read_orbit(self):
        self.orbit.read_virtual_orbit(self.lat_errors)
        self.send_data()

    def seed_and_read(self):
        self.seed_errors()
        self.read_orbit()

    def send_data(self):
        self.emit(SIGNAL("plot()"))
        #print self.z

    def find_resp(self):
        self.orbit.ring_response_matrix(self.lat_errors)
        #real_resp = measure_response_matrix(orbit, lat_errors)
        #orbit.resp = real_resp
        self.emit(SIGNAL("resp_ok()"))

    def correction(self):
        self.lat_errors = self.orbit.correction( self.lat_errors)
        max_coe = []
        for elem in self.lat_errors.sequence:
            if elem.type == "hcor" or elem.type == "vcor":
                #print elem.angle
                max_coe.append(elem.angle)
        print "max corrector angle = ", max(max_coe), min(max_coe)
        self.read_orbit()


class AppForm(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        self.create_menu()
        self.suOrbit = SuperOrbit()
        self.create_main_frame()
        self.lastFile = ''
        self.track = 0
        if os.path.exists("file.txt"):
            f = open("file.txt", 'r')
            self.lastFile = cPickle.load(f)
            self.suOrbit.open_file(self.lastFile)
            f.close()
            self.create_points()
        else:
            self.errors_button.setEnabled(False)
        #self.suOrbit = SuperOrbit()


    def create_main_frame(self):
        self.main_frame = QWidget()

        self.create_graphics()
        self.errors_button = QPushButton("seed errors")
        self.cor_button = QPushButton("correction")
        self.cor_button.setEnabled(False)
        self.resp_button = QPushButton("response matrix calculation")
        self.resp_button_meas = QPushButton("response matrix measurement")
        self.resp_button.setEnabled(False)
        self.resp_button_meas.setEnabled(False)
        self.cb_track = QCheckBox('Show trajectory')
        #cb.toggle()
        self.cb_track.stateChanged.connect(self.show_track)

        QObject.connect(self.errors_button, SIGNAL("clicked()"), self.suOrbit.seed_and_read )
        QObject.connect(self.suOrbit, SIGNAL("plot()"), self.draw_plot )
        QObject.connect(self.resp_button, SIGNAL("clicked()"), self.suOrbit.find_resp )
        QObject.connect(self.suOrbit, SIGNAL("resp_ok()"), self.switch_correction )
        QObject.connect(self.suOrbit, SIGNAL("resp_off()"), self.switch_off )
        QObject.connect(self.cor_button, SIGNAL("clicked()"), self.suOrbit.correction )
        # Layout with box sizers
        #
        grid = QGridLayout()
        grid.setSpacing(5)
        grid.addWidget(self.canvas, 0,0, 3,3)
        grid.addWidget(self.mpl_toolbar, 4,0, 1,3)
        grid.addWidget(self.resp_button, 5,0, 1,1)
        grid.addWidget(self.resp_button_meas, 6, 0, 1,1)
        grid.addWidget(self.cor_button, 5,1,1,1)
        grid.addWidget(self.errors_button, 5,2,1,1)
        grid.addWidget(self.cb_track, 6,1)

        self.main_frame.setLayout(grid)
        self.setCentralWidget(self.main_frame)


    def create_menu(self):
        self.file_menu = self.menuBar().addMenu("&File")
        load_file_action = self.create_action("&Open file", shortcut="Ctrl+S", slot=self.open_file, tip="Open input file")
        quit_action = self.create_action("&Quit", slot=self.close, shortcut="Ctrl+Q", tip="Close the application")
        self.add_actions(self.file_menu,(load_file_action, None, quit_action))
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About",shortcut='F1', slot=self.on_about, tip='About the demo')
        self.add_actions(self.help_menu, (about_action,))

    def create_action(  self, text, slot=None, shortcut=None, icon=None, tip=None, checkable=False, signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def open_file(self):
        file_choices = "*.py;*.inp"
        path = QFileDialog.getOpenFileNames(self, 'Open file', self.lastFile,file_choices)

        if len(path) != 0:
            self.path = unicode(path[0])
            f = open("file.txt", 'w')
            cPickle.dump(self.path, f)
            f.close()
            self.analysis()

    def on_about(self):
        msg = """ Orbit correction
        """
        QMessageBox.about(self, "About the demo", msg.strip())

    def analysis(self):
        self.suOrbit.open_file(self.path)
        self.create_points()

    def create_points(self):
        self.points_x = []
        self.points_y = []
        for i, bpm in enumerate(self.suOrbit.orbit.bpms):
            point_x, = self.axes.plot([], 'ro')
            point_y, = self.axes2.plot([], 'ro')
            self.points_x.append(point_x)
            self.points_y.append(point_y)

        self.errors_button.setEnabled(True)


    def create_graphics(self):
        self.setWindowTitle('Orbit correction')
        self.resize(800, 600)
        self.points_with_annotation = []
        self.fig = Figure()
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)
        self.fig.clear()
        self.fig.canvas.set_window_title('Horizontal orbit [m]')
        self.axes = self.fig.add_subplot(211)
        self.axes.set_title("Horizontal orbit [m]")
        self.axes.set_ylabel("X, [m]")
        self.axes.grid(True)

        self.fig.canvas.set_window_title('Vertical orbit [m]')
        self.axes2 = self.fig.add_subplot(212)
        self.axes2.set_title("Vertical orbit [m]")
        self.axes2.set_xlabel("S, [m]")
        self.axes2.set_ylabel("Y, [m]")
        self.axes2.grid(True)

        self.p1, = self.axes.plot([], 'ro-',lw=2.0)
        self.p2, = self.axes.plot([], 'bo-', lw=2.0)
        self.show_x, = self.axes.plot([], 'black', lw=1.0)
        self.axes.legend( [r'$old$', r'$new$'])
        self.p3, = self.axes2.plot([],'ro-', lw=2.0)
        self.p4, = self.axes2.plot([], 'bo-',lw=2.0)
        self.show_y, = self.axes2.plot([], 'black', lw=1.0)
        self.axes2.legend( [r'$old$', r'$new$'])
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)

        #self.canvas.mpl_connect('pick_event', self.line_picker)
        self.canvas.mpl_connect('motion_notify_event', self.on_move)

    def show_track(self):
        #print self.cb_track.checkState(), self.track
        if self.track and self.cb_track.checkState():
            self.suOrbit.orbit.calc_track(self.suOrbit.lat_errors)
            self.show_x.set_data(self.suOrbit.orbit.s_track, self.suOrbit.orbit.x_track)
            self.show_y.set_data(self.suOrbit.orbit.s_track, self.suOrbit.orbit.y_track)
        else:
            self.show_x.set_data([], [])
            self.show_y.set_data([], [])
        self.canvas.draw()

    def draw_plot(self):
        self.track = 1
        s = map(lambda b: b.s, self.suOrbit.orbit.bpms)
        x = map(lambda b: b.x, self.suOrbit.orbit.bpms)
        y = map(lambda b: b.y, self.suOrbit.orbit.bpms)


        if len(self.s) == 0:
            self.s = s
            self.x = x
            self.y = y


        self.p1.set_data(self.s,self.x)
        self.p3.set_data(self.s,self.y)
        self.p2.set_data(s,x)
        self.p4.set_data(s,y)
        self.axes.grid(True)
        self.axes2.grid(True)

        self.points_with_annotation = []
        max_x = max(self.x)

        for i, bpm in enumerate(self.suOrbit.orbit.bpms):
            point = self.points_x[i]
            point.set_data(self.s[i],self.x[i])

            annotation = self.axes.annotate(bpm.id, xy=(self.s[i], self.x[i]), xycoords='data',
            xytext=(self.s[i], self.x[i]+max_x/5), textcoords='data', horizontalalignment="upper",
            arrowprops=dict(arrowstyle="simple", connectionstyle="arc3,rad=-0.2"),
            bbox=dict(boxstyle="round", facecolor="w", edgecolor="0.5", alpha=0.9))

            # by default, disable the annotation visibility
            annotation.set_visible(False)
            self.points_with_annotation.append([point, annotation])
        max_y = max(self.y)

        for i, bpm in enumerate(self.suOrbit.orbit.bpms):
            point = self.points_y[i]
            point.set_data(self.s[i],self.y[i])

            annotation = self.axes2.annotate(bpm.id, xy=(self.s[i], self.y[i]), xycoords='data',
            xytext=(self.s[i], self.y[i]+max_y/5), textcoords='data', horizontalalignment="upper",
            arrowprops=dict(arrowstyle="simple", connectionstyle="arc3,rad=-0.2"),
            bbox=dict(boxstyle="round", facecolor="w", edgecolor="0.5", alpha=0.9))

            # by default, disable the annotation visibility
            annotation.set_visible(False)
            self.points_with_annotation.append([point, annotation])
        self.show_track()
        self.s = s
        self.x = x
        self.y = y

        self.axes.relim()
        self.axes.autoscale_view(True, True, True)
        self.axes2.relim()
        self.axes2.autoscale_view(True, True, True)
        self.canvas.draw()
        #self.s = copy.deepcopy(s)
        #self.x = copy.deepcopy(x)
        #self.y = copy.deepcopy(y)
        self.resp_button.setEnabled(True)

    def switch_correction(self):
        self.cor_button.setEnabled(True)


    def switch_off(self):
        self.cor_button.setEnabled(False)
        self.track = 0
        self.s = []
        self.x = []
        self.y = []

    def on_move(self, event):
        visibility_changed = False
        for point, annotation in self.points_with_annotation:
            should_be_visible = (point.contains(event)[0] == True)

            if should_be_visible != annotation.get_visible():
                visibility_changed = True
                annotation.set_visible(should_be_visible)

        if visibility_changed:
            self.canvas.draw()

    """
    def line_picker(self, event):

        if isinstance(event.artist, Line2D):
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind

            pickx = np.take(xdata, ind)
            picky = np.take(ydata, ind)
            #props = dict(ind=ind, pickx=pickx, picky=picky)
            print('onpick1 line:', zip(np.take(xdata, ind), np.take(ydata, ind)))
            msg = "You've clicked on a bar with coords:\n %s" % [pickx, picky]

            QMessageBox.information(self, "Click!", msg)
            #return props

    def on_pick(self, event):
        # The event received here is of the type
        # matplotlib.backend_bases.PickEvent
        #
        # It carries lots of information, of which we're using
        # only a small amount here.
        #
        #box_points = event.artist.get_bbox().get_points()
        print event.pickx, event.picky
        msg = "You've clicked on a bar with coords:\n %s" % event.pickx, event.picky

        QMessageBox.information(self, "Click!", msg)
    """
def main():
    app = QApplication(sys.argv)
    form = AppForm()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()