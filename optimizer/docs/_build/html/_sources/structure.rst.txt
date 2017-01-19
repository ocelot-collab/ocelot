Project Structure
=================

This section details information about the most important files and classes used int the project

OcelotInterface Main File
-------------------------
| This is the GUI main file. It contrains primarility all of the code for controling the UI.
| The GUI design follows a similar styles to other PyQT GUIs at SLAC. 
| For the most part the optimization algorthims are kept separate in a file mint/mint.py in the projjust 

:class:`optimizer.generic_optim.OcelotInterfaceWindow`

.. literalinclude:: ../../generic_optim.py
        :language: python
        :pyobject: OcelotInterfaceWindow.scan_method_select

The scan_method_select() method returns a threaded scanner object from mint/mint.py

:class:`optimizer.generic_optim.OcelotInterfaceWindow.scan_method_select`

.. literalinclude:: ../../generic_optim.py
        :language: python
        :pyobject: OcelotInterfaceWindow.start_scan
        :emphasize-lines: 6,37



| The start_scan() method starts when the 'Start scan' button is hit in the UI.
| A list of selected PV for the scan is called from the reset panel widget.
| The thread scanner object is returned from scan_method_select() and started. 
| :class:`optimizer.generic_optim.OcelotInterfaceWindow.start_scan`


XFEL Interface File
-------------------

This file is used a wrapper to translate between requests by the optimizer and the control system.
.. It also serves to hold the scan data, as well as perform normalization and unformalization for the optimzer class.

| There are two classes, the MachineInterface and DeviceProperties class.
| For the most part only the MachineInterface is used. The DeviceProperties is used now for getting limits from the GUI.

:class:`OcelotInterface.mint.xfel_interface`

**Getters and Setters**

Basic getter and setter methods using epics:
.. The SASE or objective function measurement is done in a separate method for averageing over the BSA waveform.

Getter function 

.. literalinclude:: ../../mint/xfel_interface.py
        :language: python
        :pyobject: XFELMachineInterface.get_value

Setter function

.. literalinclude:: ../../mint/xfel_interface.py
        :language: python
        :pyobject: XFELMachineInterface.set_value

.. *Saving Data*

.. Durring a scan the XFELMachineInterface is used to save data for evey step in the scan. The function *get_sase()* is used to trigger a save event. Everytime an optimizer object calls this funciton, data is saved for the setpoint of every devices and the objective funciton. :class:`OcelotInterface.mint.lcls_interface.get_sase`

.. When a scan is finished, the data is written to a file in the matlab data directory using a module "matlog.py" imported from the python toolbox.

.. **Normalization**

.. | If desired, a normalized value is passed to the optimizer based on an input mean and standrard deviation.
.. | If not normalizing the simplex, will use a default value of %5 difference from the devices current setpoint.

.. .. literalinclude:: ../../mint/lcls_interface.py
..        :language: python
..        :pyobject: LCLSMachineInterface.normalize

.. In order to caput optimizer output back to the control system, a similar unnormalize funcion is used.

.. Scanner Threads File
   --------------------

.. | This file contains scanner classes that are subclassed from a python thread object.
.. | When a new scan is started from the GUI, the scanner runs as a separate thread to avoid tying up the GUI
.. :class:`OcelotInterface.scannerThreads`

**Ocelot Optimizer**

| The Ocelot Optimizer uses the scipy optimizer fmin function to run an optimization.
| It has an argument to pass in different methods of optimizations that the scipy.optimize.fmin function can run

.. literalinclude:: ../../mint/mint.py
        :language: python
        :pyobject: Optimizer.__init__

**Simplex**

| The Ocelot scanner uses the scipy optimizer fmin function to run an optimization.
| It has an argument to pass in different methods of optimizations that the scipy.optimize.fmin function can run

.. literalinclude:: ../../mint/mint.py
        :language: python
        :pyobject: Simplex.__init__


**GP Scanner**

| This is the threaded object that runs the GP scanner.
| Below shows the initialization for the object, and required arguments.

.. literalinclude:: ../../mint/mint.py
        :language: python
        :pyobject: GaussProcess.__init__

Reset Panel Module
------------------

| The resetpanel is stand alone widget used to control the devices.
| More details in the Usage section of this document. 
| This is the origonal version of reset panel with checkboxes
:class:`optimizer.resetpanel.resetpanel`
| This version subclasses resetpanel and adds in the active checkbox column
:class:`optimizer.resetpanel.resetpanelbox`

.. epicsGet Class
.. --------------

.. | The epicsGet class if a very simple wrapper of the epics caget funciton.
.. | It's only purpose is to deal with some issues of network connectivity, and the pyepics.caget function returning None or NaN

.. .. literalinclude:: ../../epicsGet.py
        :language: python
        :pyobject: epicsGet


.. Parameters and File IO
.. ----------------------
..
.. GP File IO
.. __________
..
.. The program needs to be able to pull information from external files in order to work correctly, based on the input for parameter files and settings in the options panel.
..
.. Ocelot Normalization Parameters
.. _______________________________
..
.. The Ocelot scipy based scanners use the normalization parameters located in *./parameters/normParams*
..
.. These parameters are auto generated from the historical ranges of the devices.
.. Format: *PVNAME,MEAN,STD*
..
.. Load funciton :class:`OcelotInterface.scannerThreads.OcelotScanner.loadNormFile`
..
.. **Hyperparameters**
..
.. | The GP will always need hyperparameter information for all the devices you wish to scan.
.. | Hyperparameters are calculated from data in a hyperparameter file *./parameters/hyperparameters.npy*
.. | The file is a numpy biniary with a statistic information for each device, binned by machine L3 end energy in GeV.
.. | These are calculated outside the ocelot GUI project, in a python module named archiveRetrieve in the python toolbox.
.. | The format of this file is a nested dictionary object, where keys of the outer dict are beam energy in GeV, and inner are keys are the PVs.
.. | Load funciton :class:`OcelotInterface.scannerThreads.GpScanner.loadHyperParams`
..
..         | data[
..         |         {"3": {"PV1": [AVE,STD], "PV2": [AVE,STD], "PVN": [AVE,STD] }
..         |         {"4": {"PV1": [AVE,STD], "PV2": [AVE,STD], "PVN": [AVE,STD] }
..         |         ...
..         |         ...
..         |         {"16": {"PV1": [AVE,STD], "PV2": [AVE,STD], "PVN": [AVE,STD] }
..         | ]
..
.. | For example when the GUI starts up, it will chooses a dict of parameters based on the current L3 beam energy from *BEND:DMP1:400:BDES*
.. | Once this data in pulled from the file it is then used to calculate hyperparameters that are input into the GP
.. |
.. | Load funciton :class:`OcelotInterface.scannerThreads.GpScanner.calcLengthScaleHP`
.. | Load funciton :class:`OcelotInterface.scannerThreads.GpScanner.calcAmpCoeffHP`
.. | Load funciton :class:`OcelotInterface.scannerThreads.GpScanner.calcNoiseHP`
..
.. Understanding how these hyperparameters affct the optimizer is a bit tricky. What we have found is that typically the length scales determine the step size that GP will take when choosing a new point. The amplitude coefficent seems to determine how much the scanner will explore versus stay on a local maximum. The noise parameter determines how much of the obj func measured is accounted for by random fluctuation and noise from the beam. We have usually letf this constant, using the log of the GDET standrard deviation.
..
.. .. literalinclude:: ../../parameters/simHyps
..         :lines: 30-35
..
.. **Matlab Seed File**
..
.. | The GP scanner needs some information to build an initial model.
.. | One method of initializing this model is loading a matlab file with data from a previous scan.
.. | The function *loadSeedData* will load a previous save file from ocelot to build the GP model.
.. | The most cases this function will not be called to build the GP, but rather we will use the simpelx seeded method.
.. | This was mostly used in development when we wanted to load pre saved data sets intead of using live data.
..
.. Load funciton :class:`OcelotInterface.scannerThreads.GpScanner.loadSeedData`
..
.. .. image:: images/ocelot_savefile.png
..         :width: 600
..
.. *Example of data from a Ocelot Scan. This data is formated into a matrix and feed into the GP*
..
..
.. Calculating Parameters with the archiveRetrieve Module
.. ------------------------------------------------------
..
.. The mean and standard deviation used in the IO files above are generated in a script named energySeparation.py. More detailed informatino on this is located in the readme file within the module. The module is located in the python toolbox */usr/local/lcls/tools/python/toolbox*
..
.. | Here are the steps to get usefull information out of this module:
..
.. * 1). Check out the files from CVS into a directory (or use the existing one in /home/physics/tcope/cvswork/tools/python/toolbox/archiveRetrieve)
.. * 2). Edit the "parameter.py" file to add in the PVs and time ranges you want to get history data for
.. * 3). Execute the "generateDataSet.py" file. This generates raw data in the "./data" folder for the PVs specified.
.. * 4). Edit the "energySeparation.py" file to add in PVs you want to generate mean and std data for. Kind of redundant with step 2, but necessary.
.. * 5). Execute the energySeparation.py file to generate the hyperparameter.npy file. Then copy the file into the OcelotInterface/parameter folder.
..
.. Ocelot Developmeng and Test Mode
.. --------------------------------
..
.. This describes the process of running ocelot in development mode, in order to test new algorithms or debugging. To do this 4 dev matlab PVs are used as dummy devices, and te GUI reads in a linear objective function is place of the GDET.
..
.. .. literalinclude:: ../../OcelotInterface.py
..         :language: python
..         :pyobject: main
..
.. The main function contains a boolean 'devmode' to turn the development mode on and off. If devmode is true, the GUI loads up the dev device PVs and reads in the and objective function from the PV *"SIOC:SYS0:ML00:CALCOUT993"*. A script to start this output can be started from the blue button on the edm dev panel.
..
.. .. literalinclude:: ../../OcelotInterface.py
..         :language: python
..         :pyobject: OcelotInterfaceWindow.devmode
..
.. This is the function that is executed when you choose dev mode. You can change options here depending on what you are testing.
..
.. .. image:: images/ocelot_dev_panel.jpg
..         :width: 400
..
.. This is the edm development panel that controls the y funciton read in durring dev mode.
..
..