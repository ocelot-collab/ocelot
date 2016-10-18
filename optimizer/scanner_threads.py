"""
File to contain scanner classes to be used in the ocelot interface GUI

Tyler Cope, 2016
modified S. Tomin 2016
"""
import threading
import pandas as pd
import numpy as np
import time
import scipy.io
import csv
import copy
#import epics

#OCELOT files
from mint.mint import Optimizer, Action, Simplex
#from mint.lcls_interface import LCLSMachineInterface, LCLSDeviceProperties
from itertools import chain

#load in GP modules
from GP.GPtools import *
from GP.OnlineGP import OGP
from GP.SPGPmodel import SPGP
from GP.BasicInterfaces import fint
import GP.bayes_optimization as BOpt

#from taperThread import Taper

class OcelotScanner(threading.Thread):
    """
    This is a threaded class that runs the simplex and other scipy methods from ocelot.

    Uses the lcls machine interface and device properties class to interface with epics.
    """
    def __init__(

            self,
            parent=None,
            method='simplex',
            norm_params_file='parameters/normParams',
            norm_params_bool=False,
            norm_amp_coeff = 1.0,
            taper_coeff = 1.0,
            ):
        """
        Init method for the ocelot scanner.

        The init method takes args specific to the scanner type.
        For instance this scanner takes in information for normalization, while the GP takes for a seed.

        Args:
                parent: 'self' from the parent GUI that creates this object
                method (str): String for the optmization method (see docs for scipy.optimize.minimize package)
                norm_params_file (str): String for directory of normalization csv file
                norm_params_bool (bool): Bool to tell scanner if it should normalize using norm_params_file
                norm_amp_coeff (float): Multiplication scaling factor to the normalization coeffs
        """
        #threading setup
        super(OcelotScanner, self).__init__()
        self.stoprequest = threading.Event()
        self.parent      = parent

        #input parameters
        self.method           = method
        self.norm_params_file = norm_params_file
        self.norm_params_bool = norm_params_bool
        self.norm_amp_coeff   = norm_amp_coeff

    def setup(self, devices, objective_func, iters=45):
        """
        Basic setup of the ocelot scanner.

        The setup method is called from the GUI when a scan is started.
        This method should have the same arguments for each scanner.

        Args:
                devices (List[str]): List of the PVs to be used in scan
                objective_func (str): String for the PV for the value to maximize
                iters (int): Number of iterations for the scanner to run
                mi: MachineInterface object from the mint file
                dp: DeviceProperties oject fromt the mint file
        """
        #define pvlist
        self.devices = devices

        self.opt = Optimizer()
        self.opt.minimizer = Simplex()
        self.opt.minimizer.max_iter = iters

        #set up the timing delays
        self.opt.timeout = self.parent.trim_delay+self.parent.data_delay

        #sequence, we only use one sequence here
        #Ocelot was designed to launch multiple sequences of scan
        self.seq = [Action(func=self.opt.max_target_func, args=[objective_func, self.devices])]

    def run(self, guimode=True):
        """
        Method to start up the thread.

        In this class, the squence length is only one, using all the PVs from the setup funciton.
        If guimode is selected, when the scan is finished the GUI update clock stops and data is saved to file.

        Args:
                guimode (bool): Bool that tells parent GUI to reset UI and save data.
                If running a GP seeded scan, dont reset GUI after simpelx is done

        Returns:
                The kill status (bool) of optimizer object if not in GUI mode
                This is used in the case that this object is used in seeded mode for the GP optimizer
        """
        #Start command for the ocelot scanner
        self.opt.eval(self.seq)
        if guimode:
            pass
            #self.parent.saveData(self.mi.data)
            #self.parent.scanFinished()
        else:
            #Flips a flag in the optimizer to exit
            #There is probably a better way to do this, but it works for now
            return self.opt.kill


class GpScanner(threading.Thread):

    """
    Class to setup gaussian process bayes optimization scanner.

    Similar style to the scipy/simplex ocelot scanner class.
    Uses the GP bayes optimization package writen by Mitch McIntire
    """
    def __init__(
            self,
            parent           = None,
            seed_file        = 'parameters/simSeed.mat',
            hyp_file         = 'parameters/hyperparameters.npy',
            seedScanBool     = False,
            seedIters        = 5,
            norm_params_bool = False,
            norm_params_file = 'parameters/normParams',
            norm_amp_coeff   = 1.0,
    ):
        """
        Initilizaiton function for the GP scanner, contains specific parameters for the GP.

        This is the threaded implementation of the GP code, setup to interface with the LCLS control system.
        It can be used in several modes listed here:

                Normal mode:
                        Builds the GP model from seed_file.mat matrix.
                        Can use older Ocelot save files, provided matlab file contains data for all pvs.
                Seeded mode:
                        Builds the GP model from a live optimization seed scan.
                        In this mode a simplex or other scan run to gather fresh data.
                        After 'seedIters' of simplex, the data collected is used to build a GP.
                        Then GP scanner then runs for the rest of the scan.

        Both modes need information about the hyperparameters, found in the 'hyperparameter.npy' file.

        Args:
                parent: 'self' from the parent GUI that creates this object.
                seed_file (str): Directory string for a matlab file to use as seed data to build model.
                hyp_file (str): Directory string for a hyperparameter csv file.
                seedScanBool (bool): Bool to determine if the scan will use a simplex as seed data.
                seedIters (int): Number of iteration to the the simple seed scan run before swithing to GP.
                seedScaleFactor (float): Mutiplicitive factor passed to the simplex scanner to scale strength of initial simpelx seed scan.

        """
        #threading setup
        super(GpScanner, self).__init__()
        self.stoprequest = threading.Event()
        self.parent      = parent

        #filenames to load in
        self.seed_file = seed_file
        self.hyp_file  = hyp_file

        #stuff for simplex seeded scan
        self.seedScanBool = seedScanBool
        self.seedIters    = seedIters
        #load in parameters for simplex
        self.norm_params_bool = norm_params_bool
        self.norm_params_file = norm_params_file
        self.norm_amp_coeff = norm_amp_coeff

        #error checker
        # TODO: add object for ErrorCheck.
        #self.initErrorCheck()

    def setup(self, devices, objective_func, iters = 45):
        """
        Basic setup procedure for the GP scan objects, input args that are common to optimizer classes.

        Similar setup funtion with the same args as the ocelot scanner.
        Does not contain option to load in MachineInterface and DeviceProperties

        Could probably just write all this into the init function if you have motivation.

        Args:
                devices (List[str]): List of the PVs to be used in scan
                objective_func (str): String for the PV for the value to maximize
                iters (int): Number of iterations for the scanner to run
        """

        #load in the pvs to scan
        self.devices = devices
        self.objective_func = objective_func
        # for testing
        self.objective_func.devices = self.devices

        #number of iterations to scan (hardcode for now)
        self.iters  = iters

        #set new timing variables
        #self.mi.secs_to_ave = self.parent.data_delay
        self.total_delay = self.parent.trim_delay+self.parent.data_delay

        # -------------- GP config setup -------------- #

        #GP parameters
        self.numBV    = 30
        self.xi       = 0.01

        #no input bounds on GP selection for now
        bnds = None

        pvs = [dev.eid for dev in devices]

        hyp_params = BOpt.HyperParams(pvs=pvs, filename=self.hyp_file)
        hyps = hyp_params.loadHyperParams(objective_func.get_energy(), detector_stat_params=objective_func.get_stat_params())


        #init model
        dim = len(devices)

        self.model = OGP(dim, hyps, maxBV=self.numBV, weighted=False)

        #load model stuff
        filename = '/u1/lcls/matlab/data/2016/2016-04/2016-04-25/OcelotScan-2016-04-25-181811.mat'

        #if you would lake to load a model from file, use this function
        #self.model = self.loadModelParams(self.model,filename)


        #if this is not a simplex seeded scan, setup the seed data from mat file and build optimizer object\
        if not self.seedScanBool:
            #load seed data file
            s_data = self.loadSeedData(self.seed_file)
            #optimzer object
            self.opt = Optimizer()
            self.opt.timeout = self.total_delay
            minimizer = BOpt.BayesOpt(self.model, self.objective_func, acq_func='EI', xi = self.xi, bounds = bnds, prior_data=pd.DataFrame(s_data))
            #self.opt = BOpt.BayesOpt(self.model, self.interface, xi=self.xi, acq_func='EI', bounds=bnds, prior_data=pd.DataFrame(s_data))
            #bool to kill scanner thread from GUI
            minimizer.devices = devices
            minimizer.max_iter = iter
            self.opt.minimizer = minimizer
            self.seq = [Action(func=self.opt.max_target_func, args=[objective_func, self.devices])]
            self.opt.kill = False


    def run(self, guimode=True):
        """
        Method to start up the thread.

        In this class, the squence length is only one, using all the PVs from the setup funciton.
        If guimode is selected, when the scan is finished the GUI update clock stops and data is saved to file.

        Args:
                guimode (bool): Bool that tells parent GUI to reset UI and save data.
                If running a GP seeded scan, dont reset GUI after simpelx is done

        Returns:
                The kill status (bool) of optimizer object if not in GUI mode
                This is used in the case that this object is used in seeded mode for the GP optimizer
        """
        #Start command for the ocelot scanner
        if self.seedScanBool:
           if not self.runSeed():
               #If the ocelot scanner return False, user aborted in the simplex scan
               self.parent.ui.pushButton.setEnabled(True)
               return
        self.opt.eval(self.seq)
        if guimode:
            pass
            #self.parent.saveData(self.mi.data)
            #self.parent.scanFinished()
        else:
            #Flips a flag in the optimizer to exit
            #There is probably a better way to do this, but it works for now
            return self.opt.kill


        #self.opt.minimizer.terminate(self.opt.error_func)
        #print "Getting last data point"
        #time.sleep(3)
        #self.objective_func.get_penalty()
        #self.parent.getPlotData()
        #print ("Done getting data")

        #save data and update GUI on finish if not run as script
        #if (self.parent != None):

            #self.saveModel()
            #self.parent.saveData(self.mi.data)
            #self.parent.scanFinished()

        #Re-enable to let user start another scan
        try:
            self.parent.ui.pushButton.setEnabled(True)
        except:
            pass


    def runSeed(self):

        """
        Method to run if user has chose a to the live simplex seed scan mode.

        1). Creates a new simplex/scipy OcelotScanner object.
        2). Takes in the current mi and dp objects so that the saved data comes out in the same file
        3). Runs the simplex/scipy for 'self.seedIters' number of iterations
        4). Takes the data fromthe simplex/scipy scan and uses that to start the GP

        There is some wierd logic here to make sure that the both scanners exit if user hits abort
        Not ideal and could get fixed later

        Returns:
                False: If the user choose to abort durring the seed optimize
                True:  If the user did not abort
        """

        #set up ocelot scan, simplex default
        self.seedThread = OcelotScanner(
                        parent=self.parent,
                        method='simplex',
                        norm_params_bool=self.norm_params_bool,
                        norm_params_file = self.norm_params_file,
                        norm_amp_coeff=self.norm_amp_coeff)
        self.seedThread.setup(
                        self.devices,
                        self.objective_func,
                        self.seedIters, #add in a different value for number of iters here\
        )
        #use the seedScan interfaces for GP

        #
        #guimode=false -> dont try to save data or stop gui on scan completion
        status = self.seedThread.run(guimode=False)
        if status:
            #exit if seed thread was killed by user
            return False

        #use the simplex data saved in the mi object as a seed for the GP

        s_data = np.append(np.vstack(self.seedThread.opt.x_data), np.transpose(-np.array([self.seedThread.opt.y_data])), axis=1)

        #Make the GP optimizer object after acquiring seed data

        #self.opt = BOpt.BayesOpt(self.model, self.interface, xi=self.xi, acq_func='EI', bounds=None, prior_data=pd.DataFrame(s_data))
        #self.opt.kill = False
        #tell the mi that it should no longer normalize output to control system after simplex seed finishes
        #otherwise the mi class will try and caput crazy values
        #self.mi.norm_params_bool = False

        self.opt = Optimizer()
        self.opt.timeout = self.total_delay
        minimizer = BOpt.BayesOpt(self.model, self.objective_func, acq_func='EI', xi=self.xi,
                                  prior_data=pd.DataFrame(s_data))
        # bool to kill scanner thread from GUI
        minimizer.devices = self.devices
        minimizer.max_iter = self.iters
        self.opt.minimizer = minimizer
        self.seq = [Action(func=self.opt.max_target_func, args=[self.objective_func, self.devices])]
        self.opt.kill = False
        self.opt.minimizer.kill = False
        return True


    def loadSeedData(self, filename):
        """ Load in the seed data from a matlab ocelot scan file.

        Input file should formated like OcelotInterface file format.
        ie. the datasets that are saved into the matlab data folder.

        Pulls out the vectors of data from the save file.
        Sorts them into the same order as this scanner objects pv list.
        The GP wont work if the data is in the wrong order and loaded data is not ordered.

        Args:
                filename (str): String for the .mat file directory.

        Returns:
                Matrix of ordered data for GP. [ len(num_iterations) x len(num_devices) ]
        """
        print()
        dout = []
        if type(filename) == type(''):
            print "Loaded seed data from file:", filename
            # stupid messy formating to unest matlab format
            din = scipy.io.loadmat(str(filename))['data']
            names = np.array(din.dtype.names)
            for dev in self.devices:
                pv = dev.eid
                pv = pv.replace(":", "_")
                if pv in names:
                    x = din[pv].flatten()[0]
                    x = list(chain.from_iterable(x))
                    dout.append(x)

            # check if the right number of PV were pulled from the file
            if len(self.devices) != len(dout):
                print ("The seed data file device length unmatched with scan requested PVs!")
                print ('PV len         = ', len(self.devices))
                print ('Seed data len = ', len(dout))
                self.parent.scanFinished()

            # add in the y values
            ydata = din[self.objective_func.eid.replace(':', '_')].flatten()[0]
            dout.append(list(chain.from_iterable(ydata)))

        # If passing seed data from a seed scan
        else:
            print ("Loaded Seed Data from Seed Scan:")
            din = filename
            for dev in self.devices:
                pv = dev.eid
                if pv in din.keys():
                    dout.append(din[pv])
            dout.append(din[self.objective_func.eid])

        # transpose to format for the GP
        dout = np.array(dout).T

        # dout = dout.loc[~np.isnan(dout).any(axis=1),:]
        dout = dout[~np.isnan(dout).any(axis=1)]

        # prints for debug
        pvs = [dev.eid for dev in self.devices]
        print("[device_1, ..., device_N] detector")
        print (pvs, self.objective_func.eid)
        print (dout)
        print()

        return dout
