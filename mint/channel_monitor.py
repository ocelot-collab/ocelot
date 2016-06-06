'''
tuning using 4 corrector direct sase optimization
'''

from time import time

from ocelot.mint.flash1_interface import FLASH1MachineInterface

mi = FLASH1MachineInterface()

while True:
    print ("{0}\t{1}".format(time(), mi.get_sase(detector='gmd_fl1_slow')))