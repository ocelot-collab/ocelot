Accelerator, radiation and x-ray optics simulation framework
# An Introduction to Ocelot

Ocelot is a multiphysics simulation toolkit designed for studying FEL and storage ring based light sources. Ocelot is written in Python. Its central concept is the writing of python's scripts for simulations with the usage of Ocelot's modules and functions and the standard Python libraries. 

Ocelot includes following main modules:
* **Charged particle beam dynamics module (CPBD)**
    - optics
    - tracking
    - matching
    - collective effects (description can be found [here](http://vrws.de/ipac2017/papers/wepab031.pdf) )
        - Space Charge (true 3D Laplace solver) 
        - CSR (Coherent Synchrotron Radiation) (1D model with arbitrary number of dipoles) (under development).
        - Wakefields (Taylor expansion up to second order for arbitrary geometry).
    - MOGA (Multi Objective Genetics Algorithm). (under development but we have already applied it for a storage ring [application](http://accelconf.web.cern.ch/AccelConf/ipac2016/papers/thpmb034.pdf))
* **Native module for spontaneous radiation calculation**
* **FEL calculations: interface to GENESIS and pre/post-processing**
* **Modules for online beam control and online optimization of accelerator performances.** [Work1](http://accelconf.web.cern.ch/accelconf/IPAC2014/papers/mopro086.pdf), [work2](https://jacowfs.jlab.org/conf/y15/ipac15/prepress/TUPWA037.PDF), [work3](http://accelconf.web.cern.ch/AccelConf/ipac2016/papers/wepoy036.pdf), [work4](https://arxiv.org/pdf/1704.02335.pdf).

Ocelot extensively  uses Python's [NumPy (Numerical Python)](http://numpy.org) and [SciPy (Scientific Python)](http://scipy.org) libraries, which enable efficient in-core numerical and scientific computation within Python and give you access to various mathematical and optimization techniques and algorithms. To produce high quality figures Python's [matplotlib](http://matplotlib.org/index.html) library is used.

It is an open source project and it is being developed by physicists from  [The European XFEL](http://www.xfel.eu/), [DESY](http://www.desy.de/) (Germany), [NRC Kurchatov Institute](http://www.nrcki.ru/) (Russia).

We still have no documentation but you can find a lot of examples in ocelot/demos/ 


## Ocelot user profile

Ocelot is designed for researchers who want to have the flexibility that is given by high-level languages such as Matlab, Python (with Numpy and SciPy) or Mathematica.
However if someone needs a GUI  it can be developed using Python's libraries like a [PyQtGraph](http://www.pyqtgraph.org/) or [PyQt](http://pyqt.sourceforge.net/Docs/PyQt4/). 

For example, you can see GUI for SASE optimization (uncomment and run next block)
