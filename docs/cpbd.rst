Charged Particle Beam Dynamics (CPBD) module
============================================

Overview
--------

Charged Particle Beam Dynamics module provides features for charged particle (electron) beam optics, including  
calculating and matching Twiss parameters, single-particle tracking as well as tracking with collective effects (CSR, space charge and wakefields)

Getting started
---------------

Import OCELOT

.. code-block:: python
   :emphasize-lines: 3,5

   from ocelot import *


Define a magnetic lattice

.. code-block:: python
   :emphasize-lines: 3,5

    q1 = Quadrupole(l = 0.3, k1 = 5)
    q2 = Quadrupole(l = 0.3, k1 = -5)
    d = Drift(l = 0.5)
    lat = MagneticLattice( (d, q1, d, q2, d, q1, d,q2,d) )


Use :py:func:`twiss` to find linear optics (Twiss) functions for given initial values

.. code-block:: python
   :emphasize-lines: 3,5

    tw0 = Twiss()
    tw0.beta_x = 5.
    tw0.alpha_x = -0.87
    tw0.beta_y = 2.1
    tw0.alpha_y = 0.96
    tws = twiss(lat, tw0)


Find periodic Twiss solution

.. code-block:: python
   :emphasize-lines: 3,5

    tws = twiss(lat)


Find periodic Twiss solution with given longitudinal resolution (500 points)

.. code-block:: python
   :emphasize-lines: 3,5

    tws = twiss(lat, nPoints=500)


Plot Twiss parameters

.. code-block:: python
   :emphasize-lines: 3,5

    from pylab import *
    plot([t.s for t in tws], [t.beta_x for t in tws])
    plot([t.s for t in tws], [t.beta_y for t in tws])

Plot Twiss parameters in the lattice display

.. code-block:: python
   :emphasize-lines: 3,5


    from ocelot.gui.accelerator import *
    plot_opt_func(lat, tws)
    show()
	
	
Linear optics functions
-----------------------
.. function:: twiss(lat [, nPoints=None])	
	
Matching
--------

.. function:: match(lattice, constarints, variables[, start=0])

   lattice a :py:class:`MagneticLattice` object

Tracking
--------



Elements
--------

.. class:: MagneticLattice

.. class:: Drift

.. class:: Quadrupole

.. class:: Bend
   same as SBend

.. class:: SBend

.. class:: RBend


Transfer maps
-------------------------------

Transfer maps define how the element map acts in tracking.  
The default transfer map attachment scheme is as follows:

* Drifts, Quadrupoles, and bends have first order transfer maps
* Sextupoles have a drift-kick-drift map 


API documentation
-----------------


.. autofunction:: ocelot.cpbd.optics.twiss

.. autofunction:: ocelot.cpbd.match.match
   

