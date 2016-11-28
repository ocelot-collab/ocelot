.. ocelot documentation master file, created by
   sphinx-quickstart on Thu Jul  7 14:19:32 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

OCELOT Documentation
====================

Overview
--------
OCELOT is a framework for synchrotron light source and FEL design and operation.


Ocelot is a multiphysics simulation toolkit designed for studying FEL and storage ring based light sources.
Ocelot is written in Python. Its central concept is the writing of python's scripts for simulations with the usage of
Ocelot's modules and functions and the standard Python libraries.

Ocelot includes following main modules:

* **Charged particle beam dynamics module (CPBD)**
    - optics
    - tracking
    - matching
    - collective effects
        - Space Charge (true 3D Laplace solver)
        - CSR (Coherent Synchrotron Radiation) (1D model with arbitrary number of dipoles) (under development).
        - Wakefields (Taylor expansion up to second order for arbitrary geometry).
    - MOGA (Multi Objective Genetics Algorithm). (under development but we have already applied it for a storage ring `application <http://accelconf.web.cern.ch/AccelConf/ipac2016/papers/thpmb034.pdf>`_
* **Native module for spontaneous radiation calculation**
* **FEL calculations: interface to GENESIS and pre/post-processing**
* **Modules for online beam control and online optimization of accelerator performances.** `Work1 <http://accelconf.web.cern.ch/accelconf/IPAC2014/papers/mopro086.pdf>`_ , `work2 <https://jacowfs.jlab.org/conf/y15/ipac15/prepress/TUPWA037.PDF>`_ , `work3 <http://accelconf.web.cern.ch/AccelConf/ipac2016/papers/wepoy036.pdf>`_ .

Ocelot extensively  uses Python's `NumPy (Numerical Python <http://numpy.org>`_ and `SciPy (Scientific Python) <http://scipy.org>`_ libraries,
which enable efficient in-core numerical and scientific computation within Python and give you access to various mathematical and optimization techniques and algorithms.
To produce high quality figures Python's `matplotlib <http://matplotlib.org/index.html>`_ library is used.

It is an open source project and it is being developed by physicists from  `The European XFEL <http://www.xfel.eu/>`_ , `DESY <http://www.desy.de/>`_ (Germany), `NRC Kurchatov Institute <http://www.nrcki.ru/>`_ (Russia).


* **Ocelot user profile**


Ocelot is designed for researchers who want to have the flexibility that is given by high-level languages such as Matlab, Python (with Numpy and SciPy) or Mathematica.
However if someone needs a GUI  it can be developed using Python's libraries like a `PyQtGraph <http://www.pyqtgraph.org/>`_ or `PyQt <http://pyqt.sourceforge.net/Docs/PyQt4/>`_ .

Downloads
---------

Installation:
-------------
You have to download from GitHub zip file.

Unzip ocelot-master.zip to your working folder ../your_working_dir/.

Rename folder ../your_working_dir/ocelot-master to ../your_working_dir/ocelot.
Add ../your_working_dir/ to PYTHONPATH
Windows 7: go to Control Panel -> System and Security -> System -> Advance System Settings -> Environment Variables. and in User variables add ../your_working_dir/ to PYTHONPATH. If variable PYTHONPATH does not exist, create it
Variable name: PYTHONPATH
Variable value: ../your_working_dir/
Linux:
$ export PYTHONPATH=**../your_working_dir/**:$PYTHONPATH


Contents:
---------
.. toctree::
   :maxdepth: 2

   cpbd
   radiation
   optics
   mint
   adaptors




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

