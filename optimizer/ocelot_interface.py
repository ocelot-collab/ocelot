#!/usr/local/lcls/package/python/current/bin/python
"""
Ocelot GUI, interface for running and testing accelerator optimization methods

This file primarily contains the code for the UI and GUI
The scanner classes are contained in an external file, scannerThreads.py
The resetpanel widget is also contained in a separate module, resetpanel

Tyler Cope, 2016
"""
#QT imports
from PyQt4.QtGui import QApplication, QFrame, QPixmap, QMessageBox
from PyQt4 import QtGui, QtCore

#normal imports
import numpy as np
#import epics
import sys
import os
import time
import pyqtgraph as pg

##logbook imports
#from re import sub
#from xml.etree import ElementTree
from shutil import copy
#from datetime import datetime
#import Image

#GUI layout file
from UIOcelotInterface import Ui_Form

#slac python toolbox imports
#from ocelot.optimizer.mint.lcls_interface import MatLog
from ocelot.optimizer.mint.lcls_interface import TestMatLog as MatLog
#local imports

import scanner_threads
from ocelot.optimizer.mint.opt_objects import *

#import taperThread

matlog = MatLog()
class OcelotInterfaceWindow(QFrame):
    """ Main class for the GUI application """
    def __init__(self):
        """
        Initialize the GUI and QT UI aspects of the application.

        Create the epicsGet class, to try and get around network errors.
        Initialize the scan parameters.
        Connect start and logbook buttons on the scan panel.
        Initialize the plotting.
        Make the timer object that updates GUI on clock cycle durring a scan.
        """
        # initialize

        QFrame.__init__(self)
        self.ui = Ui_Form()
        self.ui.setupUi(self)
        #init getter class to combat pyepics errors
        #self.epicsGet = mi.getter# ocelot.optimizer.mint.lcls_interface.epicsGet()

        #method to get defaults for all the scan parameters
        self.setScanParameters()

        #scan panel button connections
        self.ui.pushButton.clicked.connect(self.runScan)
        #self.ui.pushButton.clicked.connect(self.runScript)

        #logbooking
        #self.ui.pushButton_2.clicked.connect(lambda:self.logbook())
        self.ui.pushButton_2.clicked.connect(lambda:self.logTextVerbose())

        #launch heatmap button
        self.ui.pushButton_3.clicked.connect(self.launchHeatMap)

        #help and documentation launch button
        self.ui.pushButton_5.clicked.connect(lambda: os.system("firefox file:///usr/local/lcls/tools/python/toolbox/OcelotInterface/docs/build/html/index.html"))

        #ocelot edm panel for development
        self.ui.pushButton_4.clicked.connect(lambda: os.system("edm -x /home/physics/tcope/edm/ocelot_dev.edl &"))

        #Save path for data, default will put the data in the current matlab data directory
        #See data logging module 'matlog'
        self.save_path = 'default'

        #init plots
        self.addPlots()

        #object funciton selectinator (gdet)
        self.setObFunc()

        #load in the dark theme style sheet
        self.loadStyleSheet()

        #timer for plots, starts when scan starts
        self.multiPvTimer = QtCore.QTimer()
        self.multiPvTimer.timeout.connect(self.getPlotData)


        #taperThread.Taper(initPVs = True, mi=self.mi) #update taper PVs with initial taper parameters


    def loadStyleSheet(self):
        """ Sets the dark GUI theme from a css file."""
        try:
            self.cssfile = "style.css"
            with open(self.cssfile, "r") as f:
                self.setStyleSheet(f.read())
        except IOError:
            print ('No style sheet found!')

    def setListener(self,state):
        """
        Method to set epics flag inducating that this GUI is running.

        Args:
                state (bool): Bools to set the PV stats flag true or false
        """
        #watcher cud flag
        try:
            self.mi.caput("PHYS:ACR0:OCLT:OPTISRUNNING",state)
        except:
            print ("No watcher cud PV found!")
        #listener application flag
        try:
            self.mi.caput("SIOC:SYS0:ML03:AO702" ,state)
        except:
            print ("No listener PV found!")

        #sets the hostname env to another watcher cud PV
        try:
            opi = os.environ['HOSTNAME']
            self.mi.caput("SIOC:SYS0:ML00:CA999",opi)
        except:
            print ("No OPI enviroment variable found")

    def devmode(self):
        """
        Used to setup a development mode for quick testing.

        This method contains settings for a dev mode on GUI startup.

        Uses the following PVs as dev devices:
                SIOC:SYS0:ML00:CALCOUT997
                SIOC:SYS0:ML00:CALCOUT998
                SIOC:SYS0:ML00:CALCOUT999
                SIOC:SYS0:ML00:CALCOUT000

        Uses the following PV as an objective function:
                SIOC:SYS0:ML00:CALCOUT993

        Best used with the epics dev control panel fromt he GUIs options panel.
        """
        #select GP alg for testing
        self.ui.comboBox.setCurrentIndex(0)

        #set dev objective function
        self.objective_func_pv = "SIOC:SYS0:ML00:CALCOUT993"

        #faster timing
        self.trim_delay = 0.2 #fast trim time
        self.data_delay = 0.2 #fast delay time

        #select normalization for simplex
        self.norm_params_bool = True

        #GP settings
        self.GP_hyp_file = "parameters/hype3.npy"
        self.GP_seed_iters   = 5
        self.GP_simpelx_seed = True

        #set the save path to tmp instead of the lcls matlab data directory
        self.save_path = '/tmp/'

    def disableStuff(self):
        """
        Method to disable development functionality durring production use.
        """
        #optmizer select
        self.ui.label_5.setEnabled( True)
        self.ui.comboBox.setEnabled(True)#Select Optimizer Algorithm

        #Scipy Normalization
        #self.ui.checkBox_2.setEnabled(False)#Use normalization File
        #self.ui.lineEdit_7.setEnabled(False)#Text file name

        #GP Scanner
        self.ui.groupBox_2.setEnabled(True)#GP Scanner Setup group

        #buttons
        self.ui.pushButton_3.setEnabled(True)#GP 2D Heatmap
        self.ui.pushButton_4.setEnabled(True)#Dev Ocelot Panel

        #simplex amp coeff
        #self.ui.lineEdit_8.setEnabled(False)#Normalization Scaling Coefficient
        #self.ui.label_11.setEnabled(False)#Normilization Text



#==============================================================#
# -------------- Start code for scan options UI -------------- #
#==============================================================#



    def setScanParameters(self):
        """
        Initialize default parameters for a scan when the GUI starts up.

        Creates connection for parameter changes on options panel.
        """
        #set number of scan iterations for any scan
        #could make this an argument from the UI as well
        self.iters = 45

        #normalization amp coeff for scipy scanner
        #multiplicative factor
        self.norm_amp_coeff = 1.0
        self.ui.lineEdit_8.setText(str(self.norm_amp_coeff))
        self.ui.lineEdit_8.returnPressed.connect(self.setNormAmpCoeff)

        #set objection method (gdet or some other pv to optimize)
        self.objective_func_pv = "GDET:FEE1:241:ENRCHSTBR"
        self.ui.lineEdit.setText(str(self.objective_func_pv))
        self.ui.lineEdit.returnPressed.connect(self.setObFunc)

        #set bool optinon to use normalization from file
        self.norm_params_bool = False
        self.ui.checkBox_2.stateChanged.connect(self.setNormParamsBool)
        self.ui.checkBox_2.setCheckState(0)

        #set the normalizaiton file
        self.norm_params_file = "parameters/hype3.npy"
        self.ui.lineEdit_7.setText(self.norm_params_file)
        self.ui.lineEdit_7.returnPressed.connect(self.setNormParamsFile)

        #set trim delay
        self.trim_delay = 0.5
        self.ui.lineEdit_2.setText(str(self.trim_delay))
        self.ui.lineEdit_2.returnPressed.connect(self.setTrimDelay)

        #set data delay
        self.data_delay = 1.5
        self.ui.lineEdit_3.setText(str(self.data_delay))
        self.ui.lineEdit_3.returnPressed.connect(self.setDataDelay)

        #set GP Seed data file
        self.GP_seed_file = "parameters/simSeed.mat"
        self.ui.lineEdit_4.setText(str(self.GP_seed_file))
        self.ui.lineEdit_4.returnPressed.connect(self.setGpSeed)

        #set GP Hyperparameters from a file
        self.GP_hyp_file = "parameters/hype3.npy"
        self.ui.lineEdit_5.setText(str(self.GP_hyp_file))
        self.ui.lineEdit_5.returnPressed.connect(self.setGpHyps)

        #set number GP simplex seed Iterations
        #used in adition to the standard number of iterations defined above
        self.GP_seed_iters = 5
        self.ui.lineEdit_6.setText(str(self.GP_seed_iters))
        self.ui.lineEdit_6.returnPressed.connect(self.setGpSeedIters)

        #set the "use GP Simplex Seed" bool for the GP optimizer class
        self.GP_simpelx_seed = True
        self.ui.checkBox.stateChanged.connect(self.setGpSimplexSeed)
        self.ui.checkBox.setCheckState(2)

        #set cycle period from sum of trim and data delay
        self.ui.label_7.setText("Cycle Period = "+str(round(self.trim_delay+self.data_delay,4)))

        #initialize algorithm names for UI, add items to combobox
        self.name1 = "Nelder-Mead Simplex"
        self.name2 = "Gaussian Process"
        self.name3 = "Conjugate Gradient"
        self.name4 = "Powell's Method"
        self.ui.comboBox.addItem(self.name1)
        self.ui.comboBox.addItem(self.name2)
        #self.ui.comboBox.addItem(self.name3)
        #self.ui.comboBox.addItem(self.name4)

        #initialize GUI with simplex method
        self.name_current = "Nelder-Mead Simplex"

        #set up taper coef default
        self.taper_coeff = 1.0
        #self.ui.lineEdit_9.setText(str(1.0))
        #self.ui.lineEdit_9.returnPressed.connect(self.setTaperCoeff)
        #
        ##initialize taper PV names
        #self.taperPVs = taperThread.Taper(mi=self.mi).getPVs()

    def scanMethodSelect(self):
        """
        Sets scanner method from options panel combo box selection.

        This method executes from the runScan() method, when the UI "Start Scan" button is pressed.

        Returns:
                 Selected scanner object
                 These objects are contrained in the scannerThreads.py file
        """
        index = self.ui.comboBox.currentIndex()
        #simplex Method
        if index == 0:
            self.name_current = self.name1
            scanner = scanner_threads.OcelotScanner(
                            parent=self,
                            norm_params_bool = self.norm_params_bool, #normalize toggle
                            norm_params_file = self.norm_params_file, #param filename
                            norm_amp_coeff   = self.norm_amp_coeff,   #strength scale
                            taper_coeff = self.taper_coeff)         #taper strength scale

        #GP Method
        if index == 1:
            self.name_current = self.name2
            #print("##############")
            scanner = scanner_threads.GpScanner(
                            parent           = self,
                            seed_file        = self.GP_seed_file,     #saved .mat file to make GP from file
                            hyp_file         = self.GP_hyp_file,      #hyperparameter data filename
                            seedScanBool     = self.GP_simpelx_seed,  #simplex seed or .mat seed toggle
                            seedIters        = self.GP_seed_iters,    #iterations of simplex to seed GP
                            norm_params_bool = self.norm_params_bool, #normalize toggle for simplex
                            norm_params_file = self.norm_params_file, #filename to use for simplex normalization
                            norm_amp_coeff   = self.norm_amp_coeff,   #strength scale for simplex normalizaiton
             )
        #print("##############")
        # Conjugate Gradient
        if index == 2:
            self.name_current = self.name3
            scanner = scanner_threads.OcelotScanner(parent=self, method='cg')

        # Powells Method
        if index == 3:
            self.name_current = self.name4
            scanner = scanner_threads.OcelotScanner(parent=self,method='powell')

        print()
        print ("Selected Scanner Object", scanner)
        print()
        return scanner

    def setObFunc(self):
        """
        Method to select new objective function PV (GDET).

        Typically the gas detector, but it could be some other calc PV.
        """
        text = str(self.ui.lineEdit.text())
        #check for blank string that will break it
        if text == '':
            self.ui.lineEdit.setStyleSheet("color: red")
            return #exit

        # TODO: add method for checking state of the detector
        state = True # self.mi.state(text)#epics.PV(str(text),connection_timeout=0.1).get()
        print ("state")


        if state != None:
            self.objective_func_pv = text
            self.ui.lineEdit.setStyleSheet("color: rgb(85, 255, 0);")
            self.plot1.setLabel('left',text=text)
        else:
            self.ui.lineEdit.setStyleSheet("color: red")

    def setNormAmpCoeff(self):
        """Changes the scaling parameter for the simplex/scipy normalization."""
        try:
            self.norm_amp_coeff = float(self.ui.lineEdit_8.text())
            print ("Norm scaling coeff = ", self.norm_amp_coeff)
        except:
            self.ui.lineEdit_8.setText(str(self.norm_amp_coeff))
            print ("Bad float for norm amp coeff")

    def setNormParamsFile(self):
        """ Set the normalization parameter directory."""
        self.norm_params_file = str(self.ui.lineEdit_7.text())
        print (self.norm_params_file)

    def setNormParamsBool(self):
        """
        Set boolian on whether the scan should normalize input/output for devices.

        This only applies to the simplex/scipy minimize methods, not the GP.
        """
        if self.ui.checkBox_2.isChecked():
            self.norm_params_bool = True
        else:
            self.norm_params_bool = False
        print ("Norm Params == ",self.norm_params_bool)

    def setTrimDelay(self):
        """
        Select a new trim time for a device from GUI line edit.

        Scanner will wait this long before starting data acquisition.
        """
        try:
            self.trim_delay = float(self.ui.lineEdit_2.text())
            self.ui.label_7.setText("Cycle Period = "+str(round(self.trim_delay+self.data_delay,4)))
            print ("Trim delay =", self.trim_delay)
        except:
            self.ui.lineEdit_2.setText(str(self.trim_delay))
            print ("bad float for trim delay")

    def setDataDelay(self):
        """
        Select time for objective method data collection time.

        Scanner will wait this long to collect new data.
        """
        try:
            self.data_delay = float(self.ui.lineEdit_3.text())
            self.ui.label_7.setText("Cycle Period = "+str(round(self.trim_delay+self.data_delay,4)))
            print ("Data delay =",self.data_delay)
        except:
            self.ui.lineEdit_3.setText(str(self.data_delay))
            print ("bad float for data delay")

    def setGpSeed(self):
        """
        Set directory string to use for the GP scanner seed file.
        """
        self.GP_seed_file = str(self.ui.lineEdit_4.text())

    def setGpHyps(self):
        """
        Set directory string to use for the GP hyper parameters file.
        """
        self.GP_hyp_file = str(self.ui.lineEdit_5.text())
        print (self.GP_hyp_file)

    def setGpSeedIters(self):
        """
        Set number of iterations to run the simplex scan before GP runs.
        """
        try:
            self.GP_seed_iters = int(self.ui.lineEdit_6.text())
            print ("GP seed iterations = ",self.GP_seed_iters)
        except:
            self.ui.lineEdit_6.setText(str(self.GP_seed_iters))
            print ("bad float for data delay")

    def setGpSimplexSeed(self):
        """
        Sets the bool to run GP in a simplex seed mode.
        """
        if self.ui.checkBox.isChecked():
            self.GP_simpelx_seed = True
        else:
            self.GP_simpelx_seed = False
        print ("GP seed bool == ",self.GP_simpelx_seed)


    def setTaperCoeff(self):
        """Changes the scaling parameter for taper optimization."""
        try:
            self.taper_coeff = float(self.ui.lineEdit_9.text())
            print ("Taper scaling coeff = ", self.taper_coeff)
        except:
            self.ui.lineEdit_9.setText(str(self.taper_coeff))
            print ("Bad float for taper amp coeff")

#======================================================================#
# -------------- Start code for setting/updateing plots -------------- #
#======================================================================#



    def getPlotData(self):
        """
        Collects data and updates plot on every GUI clock cycle.
        """
        #get x,y obj func data from the machine interface
        try:
            #print("get plot data = ", self.thread.mi.data[self.thread.mi.detector])
            #print("get plot data = ", self.mi.data[self.mi.detector])
            #y = self.thread.mi.data[self.thread.mi.detector]
            y = self.objective_func.y
        except:
            self.scanFinished()

        #x = np.array(self.thread.mi.data['timestamps'])-self.scanStartTime
        x = np.array(self.objective_func.x) #- self.scanStartTime
        #print("X ==== ", np.array(self.thread.mi.data['timestamps']), x, self.scanStartTime)
        #set data to liek pg line object
        self.obj_func_line.setData(x=x, y=y)

        #plot data for all devices being scanned
        for dev in self.devices:
            y = np.array(dev.data)-self.multiPlotStarts[dev.eid]
            #print(self.scanStartTime, dev.time)
            x = np.array(dev.time) - self.scanStartTime
            line = self.multilines[dev.eid]
            line.setData(x=x, y=y)

    def addPlots(self):
        """
        Initializes the GUIs plots and labels on startup.
        """

        #setup plot 1 for obj func monitor
        self.plot1 = pg.PlotWidget(title = "Objective Function Monitor",labels={'left':str(self.objective_func_pv),'bottom':"Time (seconds)"})
        self.plot1.showGrid(1,1,1)
        self.plot1.getAxis('left').enableAutoSIPrefix(enable=False) # stop the auto unit scaling on y axes
        layout = QtGui.QGridLayout()
        self.ui.widget_2.setLayout(layout)
        layout.addWidget(self.plot1,0,0)

        #setup plot 2 for device monitor
        self.plot2 = pg.PlotWidget(title = "Device Monitor",labels={'left':"Device (Current - Start)",'bottom':"Time (seconds)"})
        self.plot2.showGrid(1,1,1)
        self.plot2.getAxis('left').enableAutoSIPrefix(enable=False) # stop the auto unit scaling on y axes
        layout = QtGui.QGridLayout()
        self.ui.widget_3.setLayout(layout)
        layout.addWidget(self.plot2,0,0)

        #legend for plot 2
        self.leg2 = customLegend(offset=(75,20))
        self.leg2.setParentItem(self.plot2.graphicsItem())

        #create the obj func line object
        color = QtGui.QColor(0,255,255)
        pen=pg.mkPen(color,width=3)
        self.obj_func_line = pg.PlotCurveItem(x=[],y=[],pen=pen,antialias=True)
        self.plot1.addItem(self.obj_func_line)

    def randColor(self):
        """
        Generate random line color for each device plotted.

        Returns:
                QColor object of a random color
        """
        hi = 255
        lo = 128
        c1 = np.random.randint(lo,hi)
        c2 = np.random.randint(lo,hi)
        c3 = np.random.randint(lo,hi)
        return QtGui.QColor(c1,c2,c3)

    def setUpMultiPlot(self, devices):
        """
        Reset plots when a new scan is started.
        """
        self.plot2.clear()
        self.multilines      = {}
        self.multiPvData     = {}
        self.multiPlotStarts = {}
        x = []
        y = []
        self.leg2.scene().removeItem(self.leg2)
        self.leg2 = customLegend(offset=(50,10))
        self.leg2.setParentItem(self.plot2.graphicsItem())

        default_colors = [QtGui.QColor(255,51,51),QtGui.QColor(51,255,51),QtGui.QColor(255,255,51),QtGui.QColor(178,102,255)]
        for i, dev in enumerate(devices):

            #set the first 4 devices to have the same default colors
            if i < 4:
                color = default_colors[i]
            else:
                color = self.randColor()

            pen=pg.mkPen(color,width=2)
            self.multilines[dev.eid]  = pg.PlotCurveItem(x,y,pen=pen,antialias=True,name=str(dev.eid))
            self.multiPvData[dev.eid] = []
            self.multiPlotStarts[dev.eid] = dev.get_value()
            self.plot2.addItem(self.multilines[dev.eid])
            self.leg2.addItem(self.multilines[dev.eid], dev.eid,color=str(color.name()))


    def launchHeatMap(self):
        """
        Launches script to display a GP heatmap of two PVs selected from table.

        Can only show data from the GUIs last scan.
        """
        pvnames = self.ui.widget.getPvsFromCbState()
        if len(pvnames) != 2:
            print ("Pick only 2 PVs for a slice!")
            return
        com = "python ../GP/analyze_script.py "+str(self.last_filename)+" "+pvnames[0]+" "+pvnames[1]+" &"
        print ('Heatmap command:',com)
        os.system(com)


#========================================================================#
# -------------- Start code for running optimization scan -------------- #
#========================================================================#

    def get_devices(self):
        # TODO: add new method for creation of devices
        self.devices = []
        for pv in self.pvs:
            self.devices.append(SLACDevice(eid=pv))
        return self.devices

    def get_object_func(self, objective_func_pv):
        # TODO: add new method for creation of objective function
        target = TestTarget(eid=objective_func_pv)
        return target

    def runScan(self):
        """
        Method to start the scanner thread object.
        """
        #pulls selected pvs from resetpanelbox widget
        self.pvs = self.ui.widget.getPvsFromCbState()
        self.devices = self.get_devices()

        self.scanStartTime = time.time()
        #self.taperPVSelected = not (set(self.pvs).intersection(self.taperPVs)) == set([]) #convoluted set method of checking if any taper PVs are being optimized
        #kill scan if user hits this button while the scan is already running
        if self.ui.pushButton.text() == "Stop scan":
            try:
                #kill taperThread if it exists
                self.taperThread.kill = True
            except:
                pass
            try:
                #kill a seeded scan
                self.thread.seedThread.opt.kill = True
                self.ui.pushButton.setEnabled(False)
            except:
                pass
            try:
                #kill a normal scan
                self.thread.opt.kill = True
            except:
                pass
            #update the GUI after killing scan
            self.scanFinished()
            return

        #set gui button and color
        self.ui.pushButton.setText("Stop scan")
        #self.ui.pushButton.setStyleSheet("color: rgb(85, 255, 127);")
        self.ui.pushButton.setStyleSheet("color: red")

        #Make new plots to show the scan data
        self.setUpMultiPlot(self.devices)
        self.multiPvTimer.start(100)

        #set listener flag to True
        index = self.ui.comboBox.currentIndex()
        if index == 0:
            self.setListener(1)
        elif index ==1:
            self.setListener(2)


        #start thread scan
        try:
            self.thread = self.scanMethodSelect()
            #print(self.thread.__class__)
            self.objective_func = self.get_object_func(self.objective_func_pv)
            self.thread.setup(self.devices, self.objective_func, iters=self.iters)
            self.thread.daemon = True
            self.thread.start()
            #Display error if the scanner object fails to inititialize and setup
        except Exception:
            print()
            print("Error starting up scan thread!")
            print()
            #print(str(e))
            self.scanFinished()
            return


    def scanFinished(self):
        """
        Reset the GUI after a scan is complete.
        """
        print ("Scan Finished")
        #set flag PV to zero
        self.setListener(0)
        #reset UI controls
        self.multiPvTimer.stop()
        try:
            tThread = self.taperThread
            msg = QMessageBox()
            msg.setText("Reset Taper?")
            msg.setWindowTitle("Reset Taper")
            msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
            retval = msg.exec_()
            if(retval == QMessageBox.Yes):
                tThread.taper.resetTaper()
            del(self.taperThread)
        except:
            pass

        self.ui.pushButton.setStyleSheet("color: rgb(85, 255, 127);")
        self.ui.pushButton.setText("Start scan")


#======================================================================#
# -------------- Start code for saveing/logbooking data -------------- #
#======================================================================#



    def saveData(self, data_new):
        """
        Save scan data to the physics matlab data directory.

        Uses module matlog to save data dict in machine interface file.
        """
        #get the first and last points for GDET gain
        self.detValStart = data_new[self.objective_func_pv][0]
        self.detValStop  = data_new[self.objective_func_pv][-1]

        #replace with matlab friendly strings
        for key in data_new:
            key2 = key.replace(":","_")
            data_new[key2] = data_new.pop(key)

        #extra into to add into the save file
        data_new["BEND_DMP1_400_BDES"]   = self.epicsGet.caget("BEND:DMP1:400:BDES")
        data_new["ScanAlgorithm"]        = str(self.name_current)      #string of the algorithm name
        data_new["ObjFuncPv"]            = str(self.objective_func_pv) #string identifing obj func pv
        data_new["NormAmpCoeff"]         = self.norm_amp_coeff

        #save data
        #import matlog
        self.last_filename = matlog.save("OcelotScan",data_new,path=self.save_path)


    def logbook(self,extra_log_text='default'):
        """
        Send a screenshot to the physics logbook.

        Args:
                extra_log_text (str): string to set if verbose text should be printed to logbook. 'default' prints only gain and algorithm
        """
        #Put an extra string into the logbook function
        log_text = "Gain ("+str(self.objective_func_pv)+"): "+str(round(self.detValStart,4))+" > "+str(round(self.detValStop,4))+"\nScan Method: "+self.name_current
        fileName, path = matlog.logbook(log_text, extra_log_text=extra_log_text)
        self.screenShot(fileName, 'png')
        copy(fileName+'.ps', path)
        copy(fileName+'.png', path)
        copy(fileName+'.xml', path)
        #if extra_log_text != 'default':
        #    log_text = log_text+'\n'+str(extra_log_text)
        #
        #curr_time = datetime.now()
        #timeString = curr_time.strftime("%Y-%m-%dT%H:%M:%S")
        #log_entry = ElementTree.Element(None)
        #severity  = ElementTree.SubElement(log_entry, 'severity')
        #location  = ElementTree.SubElement(log_entry, 'location')
        #keywords  = ElementTree.SubElement(log_entry, 'keywords')
        #time      = ElementTree.SubElement(log_entry, 'time')
        #isodate   = ElementTree.SubElement(log_entry, 'isodate')
        #log_user  = ElementTree.SubElement(log_entry, 'author')
        #category  = ElementTree.SubElement(log_entry, 'category')
        #title     = ElementTree.SubElement(log_entry, 'title')
        #metainfo  = ElementTree.SubElement(log_entry, 'metainfo')
        #imageFile = ElementTree.SubElement(log_entry, 'link')
        #imageFile.text = timeString + '-00.ps'
        #thumbnail = ElementTree.SubElement(log_entry, 'file')
        #thumbnail.text = timeString + "-00.png"
        #text      = ElementTree.SubElement(log_entry, 'text')
        #log_entry.attrib['type'] = "LOGENTRY"
        #category.text = "USERLOG"
        #location.text = "not set"
        #severity.text = "NONE"
        #keywords.text = "none"
        #time.text = curr_time.strftime("%H:%M:%S")
        #isodate.text =  curr_time.strftime("%Y-%m-%d")
        #metainfo.text = timeString + "-00.xml"
        #fileName = "/tmp/" + metainfo.text
        #fileName=fileName.rstrip(".xml")
        #log_user.text = " "
        #title.text = unicode("Ocelot Interface")
        #text.text = log_text
        #if text.text == "": text.text = " " # If field is truly empty, ElementTree leaves off tag entirely which causes logbook parser to fail
        #xmlFile = open(fileName+'.xml',"w")
        #rawString = ElementTree.tostring(log_entry, 'utf-8')
        #parsedString = sub(r'(?=<[^/].*>)','\n',rawString)
        #xmlString=parsedString[1:]
        #xmlFile.write(xmlString)
        #xmlFile.write("\n")  # Close with newline so cron job parses correctly
        #xmlFile.close()
        #self.screenShot(fileName,'png')
        #path = "/u1/lcls/physics/logbook/data/"
        #copy(fileName+'.ps', path)
        #copy(fileName+'.png', path)
        #copy(fileName+'.xml', path)

    def screenShot(self,filename,filetype):
        """
        Takes a screenshot of the whole gui window, saves png and ps images to file

        Args:
                fileName (str): Directory string of where to save the file
                filetype (str): String of the filetype to save
        """
        s = str(filename)+"."+str(filetype)
        p = QPixmap.grabWindow(self.winId())
        p.save(s, 'png')
        im = Image.open(s)
        im.save(s[:-4]+".ps")
        p = p.scaled(465,400)
        #save again a small image to use for the logbook thumbnail
        p.save(str(s), 'png')


#============================================================================#
# -------------- Code for script mode running and documenation --------------#
#============================================================================#


    def runScript(self):
        """
        Method to run a set of scripted commands throught the GUI.

        This is used in development/MD for quickly scanning many sets in sequence.
        """
        #kill scan if user hits this button while the scan is already running
        if self.ui.pushButton.text() == "Stop scan":
            try:
                #kill a seeded scan
                self.thread.seedThread.opt.kill = True
            except:
                pass
            try:
                self.thread.opt.kill = True
            except:
                pass

            #update the GUI after killing scan
            self.scanFinished()
            return

        #obj func PV
        #self.objective_func_pv = "SIOC:SYS0:ML00:CALCOUT993"

        #Timesetup
        self.trim_delay = 0.5 #fast trim time
        self.data_delay = 1.5 #fast delay time
        #self.trim_delay = 0.1 #fast trim time
        #self.data_delay = 0.1 #fast delay time

        #Iters
        self.iters = 45

        #GP MODE
        self.GP_seed_iters   = 5
        self.GP_simpelx_seed = True

        #select ocelot normalization for OcelotScanner
        self.norm_params_bool = True

        # ------- START NORM SCRIPT THE SCRIPT ------- #

        print ("STARTING SCRIPT RUN")

        self.norm_amp_coeff = 0.5
        self.scriptScan(index = 0)

        self.norm_amp_coeff = 1.0
        self.scriptScan(index = 0)

        self.norm_amp_coeff = 2.0
        self.scriptScan(index = 0)

        self.norm_amp_coeff = 3.0
        self.scriptScan(index = 0)


    def scriptScan(self,index):
        """
        Single scan to run in the script.

        The sequency of events:
        * Select optimizer type based on index arg
        * Execute runScan() method
        * Start the wait() method until thread finishes
        * Send data to logbook (data file saved separatly automaticaly)
        * Reset all devices to start

        Args:
                index (int): Integer to choose algorithm method from combobox
        """
        #select algorithm
        self.ui.comboBox.setCurrentIndex(index)
        #start the scan
        self.runScan()
        #wait for scan to finish
        self.wait()
        #extra log material
        self.logTextVerbose()
        #reset devices
        self.ui.widget.resetAll()
        time.sleep(2)

    def logTextVerbose(self):
        """
        Logbook method with extra info in text string>
        """
        e1 = "Iterations: "+str(self.iters)+"\n"
        e2 = "Trim delay: "+str(self.trim_delay)+"\n"
        e3 = "Data delay: "+str(self.data_delay)+"\n"
        e4 = "Using Scipy Normalization: "+str(self.norm_params_bool)+"\n"
        e5 = "Normalization Amp Coeff: "+str(self.norm_amp_coeff)+"\n"
        e6 = "Using Live Simplex Seed: "+str(self.GP_simpelx_seed)+"\n"
        e7 = "Iters of simplex Seed: "+str(self.GP_seed_iters)+"\n"

        extra_log_text = e1+e2+e3+e4+e5+e6+e7
        self.logbook(extra_log_text)

    def wait(self):
        """
        Method to wait for thread completion.

        When running in script mode for development, this method waits until the scan thread has finished to launch the next test scan. Without this, the GUI try to start may threads at the same time instead of ordered,
        """
        while 1:
            time.sleep(0.1)
            QApplication.processEvents()
            try:
                if not self.thread.is_alive():
                    print ('Breaking out')
                    break
            except:
                #thread object is removed when it finishes sometimes
                break
        print ("WAIT FINISHED")



#==========================================================================#
# -------------- Start code for reformating the plot legend -------------- #
#==========================================================================#



# Ignore most of thus stuff, only cosmetic for device plot

class customLegend(pg.LegendItem):
    """
    STUFF FOR PG CUSTOM LEGEND (subclassed from pyqtgraph).
    Class responsible for drawing a single item in a LegendItem (sans label).
    This may be subclassed to draw custom graphics in a Legend.
    """
    def __init__(self,size=None,offset=None):
        pg.LegendItem.__init__(self,size,offset)

    def addItem(self, item, name, color="CCFF00"):

        label = pg.LabelItem(name,color=color,size="6pt",bold=True)
        sample = None
        row = self.layout.rowCount()
        self.items.append((sample, label))
        self.layout.addItem(sample, row, 0)
        self.layout.addItem(label, row, 1)
        self.layout.setSpacing(0)

class ItemSample(pg.GraphicsWidget):
    """ MORE STUFF FOR CUSTOM LEGEND """

    ## Todo: make this more generic; let each item decide how it should be represented.
    def __init__(self, item):
        pg.GraphicsWidget.__init__(self)
        self.item = item

    def boundingRect(self):
        return QtCore.QRectF(0, 0, 20, 20)

    def paint(self, p, *args):
        #p.setRenderHint(p.Antialiasing)  # only if the data is antialiased.
        opts = self.item.opts

        if opts.get('fillLevel',None) is not None and opts.get('fillBrush',None) is not None:
            p.setBrush(fn.mkBrush(opts['fillBrush']))
            p.setPen(fn.mkPen(None))
            p.drawPolygon(QtGui.QPolygonF([QtCore.QPointF(2,18), QtCore.QPointF(18,2), QtCore.QPointF(18,18)]))

        if not isinstance(self.item, ScatterPlotItem):
            p.setPen(fn.mkPen(opts['pen']))
            p.drawLine(2, 18, 18, 2)

        symbol = opts.get('symbol', None)
        if symbol is not None:
            if isinstance(self.item, PlotDataItem):
                opts = self.item.scatter.opts

            pen = fn.mkPen(opts['pen'])
            brush = fn.mkBrush(opts['brush'])
            size = opts['size']

            p.translate(10,10)
            path = drawSymbol(p, symbol, size, pen, brush)



#==============================================================#
# --------------- main method for starting GUI --------------- #
#==============================================================#



def main():

    """
    Funciton to start up the main program.

    Slecting a PV parameter set:
    If launched from the command line will take an argument with the filename of a parameter file.
    If no argv[1] is provided, the default list in ./parameters/lclsparams is used.

    Development mode:
    If devmode == False - GUI defaults to normal parameter list, defaults to nelder mead simplex
    if devmode == True  - GUI uses 4 development matlab PVs and loaded settings in the method "devmode()"
    """

    #try to get a pv list file name from commandline arg
    #this goes into initializing the reset panel PVs that show up in the GUI
    try:
        pvs = sys.argv[1]   # arg filename of params
    except:
        pvs = 'parameters/lclsparams.txt'#default filename

    #make pyqt threadsafe
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)
    #create the application
    app    = QApplication(sys.argv)


    #dp = TestLCLSDeviceProperties()
    #mi = TestLCLSMachineInterface()


    window = OcelotInterfaceWindow()

    #setup development mode if devmode==True
    devmode = False
    #devmode = True
    if devmode:
        pvs = ["SIOC:SYS0:ML00:CALCOUT000",
               "SIOC:SYS0:ML00:CALCOUT999",
               "SIOC:SYS0:ML00:CALCOUT998",
               "SIOC:SYS0:ML00:CALCOUT997"]
        window.devmode()
    else:
        pass
        window.disableStuff()

    #Build the PV list from dev PVs or selected source
    window.ui.widget.getPvList(pvs)

    #set checkbot status
    if not devmode:
        window.ui.widget.uncheckBoxes()

    #show app
    window.setWindowIcon(QtGui.QIcon('ocelot.png'))
    window.show()

    #Build documentaiton if source files have changed
    # TODO: make more universal
    #os.system("cd ./docs && xterm -T 'Ocelot Doc Builder' -e 'bash checkDocBuild.sh' &")
    #exit script
    sys.exit(app.exec_())

if __name__ == "__main__":

    main()
