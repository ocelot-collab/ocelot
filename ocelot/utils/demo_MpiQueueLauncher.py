#!/usr/bin/env python3.6

import sys, os
import matplotlib.pyplot as plt
from ocelot.utils.xfel_utils import *
from ocelot.cpbd.magnetic_lattice import MagneticLattice
from ocelot.gui.genesis_plot import *
import time
from ocelot.common.globals import *  #import of constants like "h_eV_s" and "speed_of_light"
from ocelot.common.py_func import *
# import logging
from ocelot import ocelog

ocelog.setLevel(logging.DEBUG)


# generate launcher template
launcher0 = MpiQueueLauncher()
launcher0.launch_wait=True	# when launch() is called from run_genesis, just submit job, don't wait for results
launcher0.jobpartition='maxwell'
# launcher0.jobpartition='maxwell,exfel-wp72' # <== how to specify multiple queues
launcher0.jobnodes=2
launcher0.jobntaskspernode=40
launcher0.jobscript_template='T_run_demo.sh'
launcher0.dir='.'  # when used with GENESIS v2 adaptor provided by ocelot, this is not needed
launcher0.jobscript_setup=False # not setting working directory for demo (currently setup procedure still has many GENESIS2 specifics)

# tidy up file that is used by GENESIS v2 job script to signal simulation complete
os.system('rm -f flag.finish')
##############
### DEMO 1 ###
##############
# First of all, let's wait for the job to finish
# This was designed as drop-in replacement of MpiLauncher(), with genesis v2 adapter in mind
my_launcher = deepcopy(launcher0)
print('Submitting first job. Only returns once job is complete...')
my_launcher.launch()
print('Job completed')


# tidy up file that is used by GENESIS v2 job script to signal simulation complete
os.system('rm -f flag.finish')
##############
### DEMO 2 ###
##############
another_launcher = deepcopy(launcher0)
# This time, we won't wait for job completion (using "submit" instead of "launch").
# Ideal for launching and managing multiple jobs.
another_launcher.launch_wait=False
another_launcher.submit()
print('Job was submitted as jobid={}'.format(another_launcher.jobid))
print('Waiting for job to finish...')
another_launcher.wait_for_completion(verbose=True) # wait for completion (call returns once job complete)
print('Job is complete')


# tidy up file that is used by GENESIS v2 job script to signal simulation complete
os.system('rm -f flag.finish')
##############
### DEMO 2 ###
##############
yet_another_launcher = deepcopy(launcher0)
# This time, we won't wait for job completion (using "submit" instead of "launch").
# Ideal for launching and managing multiple jobs.
yet_another_launcher.launch_wait=False
yet_another_launcher.submit()
print('Job was submitted as jobid={}'.format(yet_another_launcher.jobid))
print('Waiting for job to finish...')
while True:
    # is job finished?
    if yet_another_launcher.is_complete()==True:
        break
    # no => do some other work
    print('doing important stuff in a loop in the demo script, then asking again if job is complete')
    time.sleep(12) # <= "important stuff"

print('Job is complete')
