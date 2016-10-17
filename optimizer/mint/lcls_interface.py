"""
Machine interface file for the LCLS to ocelot optimizer

Tyler Cope, 2016
"""

import time
import sys

import numpy as np
try:
    import epics
except:
    pass
#logbook imports
from re import sub
from xml.etree import ElementTree

from datetime import datetime

sys.path.append("..")
#from ocelot.optimizer.taperThread import Taper
#import matlog


class MatLog:
    def __init__(self):
        pass

    def save(self, name, data_new, path):
        #name = "OcelotScan"
        last_filename = matlog.save(name, data_new, path=path)
        #last_filename = "test.txt"
        return last_filename

    def logbook(self, log_text, extra_log_text='default'):
        """
        Send a screenshot to the physics logbook.

        Args:
                extra_log_text (str): string to set if verbose text should be printed to logbook. 'default' prints only gain and algorithm
        """
        #Put an extra string into the logbook function
        #log_text = "Gain ("+str(self.objective_func_pv)+"): "+str(round(self.detValStart,4))+" > "+str(round(self.detValStop,4))+"\nScan Method: "+self.name_current
        if extra_log_text != 'default':
            log_text = log_text+'\n'+str(extra_log_text)

        curr_time = datetime.now()
        timeString = curr_time.strftime("%Y-%m-%dT%H:%M:%S")
        log_entry = ElementTree.Element(None)
        severity  = ElementTree.SubElement(log_entry, 'severity')
        location  = ElementTree.SubElement(log_entry, 'location')
        keywords  = ElementTree.SubElement(log_entry, 'keywords')
        time      = ElementTree.SubElement(log_entry, 'time')
        isodate   = ElementTree.SubElement(log_entry, 'isodate')
        log_user  = ElementTree.SubElement(log_entry, 'author')
        category  = ElementTree.SubElement(log_entry, 'category')
        title     = ElementTree.SubElement(log_entry, 'title')
        metainfo  = ElementTree.SubElement(log_entry, 'metainfo')
        imageFile = ElementTree.SubElement(log_entry, 'link')
        imageFile.text = timeString + '-00.ps'
        thumbnail = ElementTree.SubElement(log_entry, 'file')
        thumbnail.text = timeString + "-00.png"
        text      = ElementTree.SubElement(log_entry, 'text')
        log_entry.attrib['type'] = "LOGENTRY"
        category.text = "USERLOG"
        location.text = "not set"
        severity.text = "NONE"
        keywords.text = "none"
        time.text = curr_time.strftime("%H:%M:%S")
        isodate.text =  curr_time.strftime("%Y-%m-%d")
        metainfo.text = timeString + "-00.xml"
        fileName = "/tmp/" + metainfo.text
        fileName=fileName.rstrip(".xml")
        log_user.text = " "
        title.text = unicode("Ocelot Interface")
        text.text = log_text
        if text.text == "": text.text = " " # If field is truly empty, ElementTree leaves off tag entirely which causes logbook parser to fail
        xmlFile = open(fileName+'.xml',"w")
        rawString = ElementTree.tostring(log_entry, 'utf-8')
        parsedString = sub(r'(?=<[^/].*>)','\n',rawString)
        xmlString=parsedString[1:]
        xmlFile.write(xmlString)
        xmlFile.write("\n")  # Close with newline so cron job parses correctly
        xmlFile.close()
        #self.screenShot(fileName, 'png')
        path = "/u1/lcls/physics/logbook/data/"
        #copy(fileName+'.ps', path)
        #copy(fileName+'.png', path)
        #copy(fileName+'.xml', path)
        return fileName, path

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

class EpicsGet:

    def __init__(self):
        pass
        ''' Separate getter class to add logic for dealing with channel access return errors '''

    def caget(self,device_name):
        #need to do this while/try loop stuff because of CA errors
        #when channel acces trys to connect for the first time in a separate thread
        #seems to be some problem with pyepics
        data = None
        ct = 0
        #print("caget")
        while 1:
            #print("caget")
            try:
                 data = 0.6# epics.caget(device_name)
                 if data == None:
                      continue
                 return data
            except:
                print ("Error retriving ca data! Tyring to caget data again")
                time.sleep(.1)
            ct+=1
            if ct > 5:
                raise Exception("Too many caget trys ,exiting")
                return None


class LCLSMachineInterface:
    """ Start machine interface class """

    def __init__(self):
        """ Initialize parameters for the scanner class. """
        self.secs_to_ave = 2         #time to integrate gas detector
        self.getter = EpicsGet()     #getter class for channel access
        self.caput = epics.caput
        self.state = lambda text: epics.PV(str(text), connection_timeout=0.1).get()
        self.inputNormParams = None  #normalization parameters
        self.norm_params_bool= False #normalization parameters
        self.taperParams = None

    #=================================================================#
    # -------------- Original interface file functions -------------- #
    #=================================================================#

    def get_alarms(self):
        """ Does not need alarms for now, proabaly dont need this with LCLS MPS. """
        return [0.0]

    def get_sase(self, seconds=None):
        """
        Returns data for the ojective function (sase) from the selected detector PV.

        At lcls the repetition is  120Hz and the readout buf size is 2800.
        The last 120 entries correspond to pulse energies over past 1 second.

        Args:
                seconds (float): Variable input on how many seconds to average data

        Returns:
                Float of SASE or other detecor measurment
        """
        datain = self.getter.caget(self.detector)
        try: #try to average over and array input
            if seconds == None: #if a resquested seconds is passed
                dataout = np.mean(datain[-(self.secs_to_ave*120):])
                sigma   = np.std( datain[-(self.secs_to_ave*120):])
            else:
                dataout = np.mean(datain[-(seconds*120):])
                sigma   = np.std( datain[-(seconds*120):])
        except: #if average fails use the scaler input
            print ("Detector is not a waveform PV, using scalar value")
            dataout = datain
            sigma   = -1

        self.record_data(dataout, sigma)
        return dataout

    def get_value(self, device_name):
        """
        Getter function for lcls.

        Args:
                device_name (str): String of the pv name used in caget

        Returns:
                Data from caget, variable data type depending on PV
        """

        return self.getter.caget(str(device_name))

    def get_energy(self):
        return self.getter.caget("BEND:DMP1:400:BDES")

    def set_value(self, device_name, val):
        """
        Setter function for lcls.

        Args:
                device_name (str): String of the pv name used in caput
                val (variable): Value to caput to device, variable data type depending on PV
        """
        unnormed = self.unnormalize(device_name, val)
        return epics.caput(device_name, unnormed)


    def init_corrector_vals(self, correctors):
        """
        Gathers starting values for the correcters/devices.

        Args:
                correctors ([str]): list of pv names

        Returns:
                Float list of corrector start values
        """
        vals = [0.0]*len(correctors) #np.zeros(len(correctors))
        for i in range(len(correctors)):
            mag_channel = correctors[i]
            val = self.getter.caget(str(mag_channel))
            #normalize val
            val = self.normalize(correctors[i],val)

            vals[i] = val
        return vals

    #========================================================================#
    # --------------- New setup functions and data recording --------------- #
    #========================================================================#

    def setUpDetector(self,pvs,detector="GDET:FEE1:241:ENRCHSTBR"):
        """
        Initializes detector parameter to optimize.
        Usefull for switching desired parameter from GUI.

        Default PV is the gas detector: GDET:FEE1:241:ENRCHSTBR

        Args:
                pvs ([str]): List of PV names
                detector (str): String of the detector PV name, usually gas detector but can change
        """
        self.detector = detector
        self.setup_data_record(pvs) #reinit the data recording

    def setup_data_record(self,pvs):
        """
        Initializing blank arrays for data storage.

        Args:
                pvs ([str]): List of pv names
        """
        self.pvs = pvs
        self.data = {} #dict of all devices deing scanned
        self.data[self.detector] = [] #detector data array
        self.data['DetectorStd'] = [] #detector std array
        self.data['timestamps']  = [] #timestamp array
        for pv in pvs:
            self.data[pv] = []

    def record_data(self, gdet, sigma):
        """
        Get data for devices everytime the SASE is measured to save syncronous data.

        Args:
                gdet (str): String of the detector PV, usually gas detector
                simga (float): Float of the measurement standard deviation

        """
        self.data[self.detector].append(gdet)
        self.data['DetectorStd'].append(sigma)
        self.data['timestamps'].append(time.time())
        for pv in self.pvs:
            self.data[pv].append(self.getter.caget(pv))


    #=======================================================#
    # -------------- Normalization functions -------------- #
    #=======================================================#


    def normalize(self, corrector, x):
        """
        Transform to normalized data for optimizer input.

        Args:
                correcter: pv name of the devices
                x: the input x val to be normalized

        Returns:
                Float normalized value of x
        """
        if(Taper().isTaperPV(corrector)):
            mu = self.taperParams[corrector][0]
            sig = self.taperParams[corrector][1]
            y = (float(x) -mu)/(sig)
            print ("SCALED TAPER",y)
            return y

        if not self.norm_params_bool:
            return x
        mu  = self.inputNormParams[corrector][0]
        sig = self.inputNormParams[corrector][1]
        y = (float(x)-mu)/(sig)
        print ("NORMALIZED",y)
        return y

    def unnormalize(self,corrector,y):

        """
        Transform back to to machine units for optimizer output

        Args:
                correcter: pv name of the devices
                x: the input x val to be normalized

        Returns:
                Float un-normalized value of x
        """
        if(Taper(mi=self).isTaperPV(corrector)):
            mu = self.taperParams[corrector][0]
            sig = self.taperParams[corrector][1]
            x = float(y)*(sig)+mu
            print ("UNSCALED TAPER",x)
            return x

        if not self.norm_params_bool:
            return y
        mu  = self.inputNormParams[corrector][0]
        sig = self.inputNormParams[corrector][1]
        x = float(y)*(sig)+mu
        print ("UN-NORMALIZED",x)
        return x

    def undsOut(self, numK, undsNotUsedIndex):
        """
        check if undulators that should be in (i.e. not in undsNotUsedIndex) are out, returns those segments
        :param numK:
        :param undsNotUsedIndex:
        :return: outPVs
        """
        outPVs = []
        for i in range(1, numK+1):
            if(i-1 not in undsNotUsedIndex):
                statPV = 'USEG:UND1:' + str(i) + '50:LOCATIONSTAT'
                if(not(self.getter.caget(statPV) == 1)):
                    outPVs.append(i-1)
        return outPVs


    def stillMoving(self, movePVs):
        """
        check if undulator segments are moving
        :param movePVs:
        :return:
        """
        time.sleep(.1)
        for pv in movePVs:
            if(not(self.getter.caget(pv) == 1)):
                print(pv, ' moving')
                return True
        return False

    def withinTol(self, newTaper, despv):
        """
        check if k values are within tolerance of undulator
        :param newTaper:
        :param despv:
        :return:
        """
        for pv, newK in zip(despv, newTaper):
            hilim = self.getter.caget(pv + '.HOPR')
            lolim = self.getter.caget(pv + '.LOPR')
            if(newK > hilim or newK < lolim):
                return False
        return True


    def init_error_check(self):
        """
        Initialize PVs and setting used in the errorCheck method.
        """
        #setup pvs to check
        self.error_bcs      = "BCS:MCC0:1:BEAMPMSV"
        self.error_mps      = "SIOC:SYS0:ML00:CALCOUT989"
        self.error_gaurdian = "SIOC:SYS0:ML00:AO466"
        self.error_und_tmit = "BPMS:UND1:3290:TMITTH"

        #pv to bypass the error pause
        self.error_bypass  = "SIOC:SYS0:ML00:CALCOUT990"
        self.error_tripped = "SIOC:SYS0:ML00:CALCOUT991"

        #set the unlatch pv to zero
        #epics.caput(self.error_bypass, 0)
        #epics.caput(self.error_tripped, 0)

        self.set_value(self.error_bypass, 0)
        self.set_value(self.error_tripped, 0)

    def error_check(self):
        while 1:
            #check for bad state
            if self.get_value(self.error_bypass)     == 1:
                out_msg="Bypass flag is TRUE"
            elif self.get_value(self.error_bcs)      != 1:
                out_msg="BCS tripped"
            elif self.get_value(self.error_mps)      != 0:
                out_msg="MPS tripped"
            elif self.get_value(self.error_gaurdian) != 0:
                out_msg="Gaurdian tripped"
            elif self.get_value(self.error_und_tmit) < 5.0e7:
                out_msg="UND Tmit Low"
            else:
                out_msg='Everything Okay'

            #exit if the stop button is set
            if not self.get_value("SIOC:SYS0:ML03:AO702"):
                break

            #set the error check message
            #epics.caput ("SIOC:SYS0:ML00:CA000",out_msg)
            self.set_value("SIOC:SYS0:ML00:CA000",out_msg)
            print (out_msg)

            #break out if error check is bypassed
            if (out_msg=="Bypass flag is TRUE"):
                break

            #break out if everything is okay
            if (out_msg=="Everything Okay"):
                #epics.caput(self.error_tripped, 0)
                self.set_value(self.error_tripped, 0)
                break
            else:
                #epics.caput(self.error_tripped, 1)
                self.set_value(self.error_tripped, 1)
            time.sleep(0.1)

class LCLSDeviceProperties:

    """ Start the device properties class """

    def __init__(self):
        self.getter = EpicsGet()

    def get_limits(self, device,percent=0.25):
        """
        Function to get device limits.
        Executes on every iteration of the optimizer function evaluation.
        Currently does not work with the normalization scheme.
        Defaults to + 25% of the devices current values.

        Args:
                device (str): PV name of the device to get a limit for
                percent (float): Generates a limit based on the percent away from the devices current value
        """
        val = self.start_values[device]
        tol = abs(val*percent)
        lim_lo = val-tol
        lim_hi = val+tol
        limits = [lim_lo,lim_hi]
        #print device, 'LIMITS ->',limits
        #return limits
        #Dosnt work with normalizaiton, big limits
        return [-10000,10000]


    def get_start_values(self, devices,percent=0.25):
        """
        Function to initialize the starting values for get_limits methomethodd.

        Called from tuning file or GUI

        Args:
                devices ([str]): PV list of devices
                percent (float): Percent around the mean to generate limits

        """
        self.start_values={}
        self.norm_minmax={}
        for d in devices:
            val = self.getter.caget(str(d))
            self.start_values[str(d)] = val
            tol = abs(val*percent)
            lim_lo = val-tol
            lim_hi = val+tol
            limits = [lim_lo,lim_hi]
            self.norm_minmax[str(d)] = [lim_lo,lim_hi]


class TestLCLSMachineInterface:
    """ Start machine interface class """

    def __init__(self):
        """ Initialize parameters for the scanner class. """
        self.secs_to_ave = 2         #time to integrate gas detector
        self.getter = EpicsGet()     #getter class for channel access
        self.caput = lambda x1, x2: 0
        self.state = lambda text: epics.PV(str(text), connection_timeout=0.1).get()
        self.inputNormParams = None  #normalization parameters
        self.norm_params_bool= False #normalization parameters
        self.taperParams = None
        self.currents = {}
    #=================================================================#
    # -------------- Original interface file functions -------------- #
    #=================================================================#

    def get_alarms(self):
        """ Does not need alarms for now, proabaly dont need this with LCLS MPS. """
        return [0.0]

    def get_sase(self, seconds=None):
        #mus = np.random.rand(len(self.currents))
        if len(self.currents) == 0:
            return 0.
        values = np.array(self.currents.values()) - np.ones(len(self.currents))
        dataout = 1./(0.1 + np.sum((np.power(values, 2))))
        #dataout = np.random.rand()
        print("SASE = ", self.currents.values(), dataout, values,np.sum((np.power(values, 2))), self.currents)
        self.record_data(dataout,0)

        return dataout

    def get_value(self, pv):
        #print(self.currents, pv)
        try:
            return self.currents[pv]
        except:
            dataout = np.random.rand()
            self.currents[pv] = dataout
            return dataout

        #dataout = np.random.rand()
        ##dataout = self.getter.caget(pv)
        #return dataout

    def get_energy(self):
        return 10# self.getter.caget("BEND:DMP1:400:BDES")

    def set_value(self, device_name, val):
        #unnormed = self.unnormalize(device_name, val)
        self.currents[device_name] = val
        return 1

    def init_corrector_vals(self, correctors):
        vals = [0.0]*len(correctors) #np.zeros(len(correctors))
        for i, cor in enumerate(correctors):
            vals[i] = self.get_value(cor)
        return vals

    #========================================================================#
    # --------------- New setup functions and data recording --------------- #
    #========================================================================#

    def setUpDetector(self,pvs,detector="GDET:FEE1:241:ENRCHSTBR"):
        self.detector = detector
        self.setup_data_record(pvs) #reinit the data recording

    def setup_data_record(self,pvs):
        """
        Initializing blank arrays for data storage.

        Args:
                pvs ([str]): List of pv names
        """
        self.pvs = pvs
        self.data = {} #dict of all devices deing scanned
        self.data[self.detector] = [] #detector data array
        self.data['DetectorStd'] = [] #detector std array
        self.data['timestamps']  = [] #timestamp array
        for pv in pvs:
            self.data[pv] = []

    def record_data(self,gdet,sigma):
        """
        Get data for devices everytime the SASE is measured to save syncronous data.

        Args:
                gdet (str): String of the detector PV, usually gas detector
                simga (float): Float of the measurement standard deviation

        """
        self.data[self.detector].append(gdet)
        self.data['DetectorStd'].append(sigma)
        self.data['timestamps'].append(time.time())
        for pv in self.pvs:
            self.data[pv].append(self.get_value(pv))
        print("record data ", self.data['timestamps'])


    #=======================================================#
    # -------------- Normalization functions -------------- #
    #=======================================================#
    def init_error_check(self):
        """
        Initialize PVs and setting used in the errorCheck method.
        """
        #setup pvs to check
        self.error_bcs      = "BCS:MCC0:1:BEAMPMSV"
        self.error_mps      = "SIOC:SYS0:ML00:CALCOUT989"
        self.error_gaurdian = "SIOC:SYS0:ML00:AO466"
        self.error_und_tmit = "BPMS:UND1:3290:TMITTH"

        #pv to bypass the error pause
        self.error_bypass  = "SIOC:SYS0:ML00:CALCOUT990"
        self.error_tripped = "SIOC:SYS0:ML00:CALCOUT991"

        #set the unlatch pv to zero
        #epics.caput(self.error_bypass, 0)
        #epics.caput(self.error_tripped, 0)

        self.set_value(self.error_bypass, 0)
        self.set_value(self.error_tripped, 0)

    def normalize(self, corrector, x):
        #if(Taper().isTaperPV(corrector)):
        #    mu = self.taperParams[corrector][0]
        #    sig = self.taperParams[corrector][1]
        #    y = (float(x) -mu)/(sig)
        #    print "SCALED TAPER",y
        #    return y

        if not self.norm_params_bool:
            return x
        mu  = self.inputNormParams[corrector][0]
        sig = self.inputNormParams[corrector][1]
        y = (float(x)-mu)/(sig)
        print ("NORMALIZED",y)
        return y
    def unnormalize(self,corrector,y):
        #if(Taper().isTaperPV(corrector)):
        #    mu = self.taperParams[corrector][0]
        #    sig = self.taperParams[corrector][1]
        #    x = float(y)*(sig)+mu
        #    print "UNSCALED TAPER",x
        #    return x

        if not self.norm_params_bool:
            return y
        mu  = self.inputNormParams[corrector][0]
        sig = self.inputNormParams[corrector][1]
        x = float(y)*(sig)+mu
        print ("UN-NORMALIZED",x)
        return x

    def undsOut(self, numK, undsNotUsedIndex):
        """
        check if undulators that should be in (i.e. not in undsNotUsedIndex) are out, returns those segments
        :param numK:
        :param undsNotUsedIndex:
        :return: outPVs
        """
        outPVs = []
        for i in range(1, numK+1):
            if(i-1 not in undsNotUsedIndex):
                statPV = 'USEG:UND1:' + str(i) + '50:LOCATIONSTAT'
                getter_caget_statPV = 1
                if(not(getter_caget_statPV == 1)):
                    outPVs.append(i-1)
        return outPVs


    def stillMoving(self, movePVs):
        """
        check if undulator segments are moving
        :param movePVs:
        :return:
        """
        return False

    def withinTol(self, newTaper, despv):
        """
        check if k values are within tolerance of undulator
        :param newTaper:
        :param despv:
        :return:
        """
        return True

    def error_check(self):
        pass


class TestLCLSDeviceProperties:

    """ Start the device properties class """

    def __init__(self):
        self.getter = EpicsGet()

    def get_limits(self, device,percent=0.25):
        return [-10000,10000]


    def get_start_values(self, devices,percent=0.25):
        """
        Function to initialize the starting values for get_limits methomethodd.

        Called from tuning file or GUI

        Args:
                devices ([str]): PV list of devices
                percent (float): Percent around the mean to generate limits

        """
        self.start_values={}
        self.norm_minmax={}
        for d in devices:
            val = self.getter.caget(str(d))
            self.start_values[str(d)] = val
            tol = abs(val*percent)
            lim_lo = val-tol
            lim_hi = val+tol
            limits = [lim_lo,lim_hi]
            self.norm_minmax[str(d)] = [lim_lo,lim_hi]


class TestMatLog:
    def __init__(self):
        pass

    def save(self, name, data_new, path):
        #name = "OcelotScan"
        #last_filename = matlog.save(name, data_new, path=path)
        last_filename = "test.txt"
        return last_filename

    def logbook(self, log_text, extra_log_text='default'):
        """
        Send a screenshot to the physics logbook.

        Args:
                extra_log_text (str): string to set if verbose text should be printed to logbook. 'default' prints only gain and algorithm
        """
        #Put an extra string into the logbook function
        #log_text = "Gain ("+str(self.objective_func_pv)+"): "+str(round(self.detValStart,4))+" > "+str(round(self.detValStop,4))+"\nScan Method: "+self.name_current
        if extra_log_text != 'default':
            log_text = log_text+'\n'+str(extra_log_text)

        curr_time = datetime.now()
        timeString = curr_time.strftime("%Y-%m-%dT%H:%M:%S")
        log_entry = ElementTree.Element(None)
        severity  = ElementTree.SubElement(log_entry, 'severity')
        location  = ElementTree.SubElement(log_entry, 'location')
        keywords  = ElementTree.SubElement(log_entry, 'keywords')
        time      = ElementTree.SubElement(log_entry, 'time')
        isodate   = ElementTree.SubElement(log_entry, 'isodate')
        log_user  = ElementTree.SubElement(log_entry, 'author')
        category  = ElementTree.SubElement(log_entry, 'category')
        title     = ElementTree.SubElement(log_entry, 'title')
        metainfo  = ElementTree.SubElement(log_entry, 'metainfo')
        imageFile = ElementTree.SubElement(log_entry, 'link')
        imageFile.text = timeString + '-00.ps'
        thumbnail = ElementTree.SubElement(log_entry, 'file')
        thumbnail.text = timeString + "-00.png"
        text      = ElementTree.SubElement(log_entry, 'text')
        log_entry.attrib['type'] = "LOGENTRY"
        category.text = "USERLOG"
        location.text = "not set"
        severity.text = "NONE"
        keywords.text = "none"
        time.text = curr_time.strftime("%H:%M:%S")
        isodate.text =  curr_time.strftime("%Y-%m-%d")
        metainfo.text = timeString + "-00.xml"
        fileName = "/tmp/" + metainfo.text
        fileName=fileName.rstrip(".xml")
        log_user.text = " "
        title.text = unicode("Ocelot Interface")
        text.text = log_text
        if text.text == "": text.text = " " # If field is truly empty, ElementTree leaves off tag entirely which causes logbook parser to fail
        xmlFile = open(fileName+'.xml',"w")
        rawString = ElementTree.tostring(log_entry, 'utf-8')
        parsedString = sub(r'(?=<[^/].*>)','\n',rawString)
        xmlString=parsedString[1:]
        xmlFile.write(xmlString)
        xmlFile.write("\n")  # Close with newline so cron job parses correctly
        xmlFile.close()
        #self.screenShot(fileName, 'png')
        path = "/u1/lcls/physics/logbook/data/"
        #copy(fileName+'.ps', path)
        #copy(fileName+'.png', path)
        #copy(fileName+'.xml', path)
        return fileName, path


