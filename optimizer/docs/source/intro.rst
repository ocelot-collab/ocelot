Introduction
============

This is deep modification of SLAC version of the Ocelot GUI (Tyler Cope, SLAC, 2016) for the European XFEL facility.

Sergey Tomin, 2017.


What is OceotInterface?
-----------------------

OcelotInterface is a a python and PyQt based GUI for running and testing accelerator optimization methods.

This project stems from work by Ilya Agapov (DESY) on the OCELOT python package, an accelerator simulaiton and interface framework.
The OCELOT package is pirmarily used for x-ray simulation but it contains files for beam optimization, which we have used and expanded upon for accelerator tuning.
The main optimization file Ocelot.util.mint was used to test the scanner on LCLS by writing a machine interface wrapper file. Soon after a GUI interface was created to facilitate easy testing.

The goal is two fold:

* To provide a stable user interface to run optimization software in standard tuning procedures. 
* To provide a testing framework for new methods of optimization

The software is now used in standard tuning and saves data to the matlab physics data directory for analysis.


What can it do?
---------------

Currently the production GUI is used to run optimization scans using the Nelder-Mead simplex algorithm, to optimize the machine FEL or some other parameter. 
The interface provides a method to select multiple tuning devices for a scan, and quickly reset devices in even of a problem. 
The options panel tab allows a user to change settings for the scan.

The software is also in development to use a new Bayesian optimization method called the Gaussian Process.
The GP scanner uses the algorithm detailed in the papers below in the reference secton.
So far initial tests have been conducted using the GP scanner, but it is not yet ready for production use. 


Resources
---------

**OCELOT Info**

* `OCELOT GitHub                 <https://github.com/iagapov/ocelot>`_
* `Simplex algorithm wiki        <https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method>`_
* `IPAC16 Automated Tuning:      <http://accelconf.web.cern.ch/AccelConf/ipac2016/papers/wepoy036.pdf>`_

**Bayesian Optimization**

* `Gaussian Processes textbook:  <http://www.gaussianprocess.org/gpml/chapters/RW.pdf>`_
* `Bayesian Optimization Text:   <http://arxiv.org/pdf/1012.2599v1.pdf>`_
* `IPAC16 Bayesian Optimization: <http://accelconf.web.cern.ch/AccelConf/ipac2016/papers/wepow055.pdf>`_
