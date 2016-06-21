'''
utilities to add externaly calculated wakes to beam files
'''
from xframework.adaptors.genesis import *
try:
    import matplotlib.animation as anim
except:
    print 'animation not installed'
import numpy as np
import matplotlib.pyplot as plt


if len(sys.argv)>3:
    command = sys.argv[1]
    beamf = sys.argv[2]    
    outf = sys.argv[3]
else:
    beamf = '/home/iagapov/tmp/run_2/tmp.beam'
    command = 'current'
    outf = 'tmp.beam'


if command == 'current':
        beam = read_beam_file(beamf)
        beam.columns = ['ZPOS','CURPEAK']
        beam.I = beam.I[::-1]
        f=open(outf,'w')
        f.write(beam_file_str(beam))
        f.close()
        #print beam_file_str(beam) 

if command == 'add':

    if len(sys.argv)>2:
        wakef = sys.argv[3]
        outf = sys.argv[4]    
    else:
        beamf = '/home/iagapov/tmp/run_2/tmp.beam'
        command = 'current'

    buf = open(wakef).read().split('\n')
    wake = []
    for l in buf:
        d = l.split()
        if len(d)>1:
            wake.append(float(d[1]))

    print 'read ', len(wake), ' slice values'

    beam = read_beam_file(beamf)
    
    beam.eloss = wake.reverse()
    
    f=open(outf,'w')
    f.write(beam_file_str(beam))
    f.close()

    beam = read_beam_file(beamf)
