import sqlite3
import os, sys
import datetime, time

from ocelot.utils.db import PerfDB

db = PerfDB(dbname = "flash.db")
#db.new_tuning({'wl':13.6, 'charge':0.1,'comment':'test tuning'})
tunes = db.get_tunings()

print ('*** all tunings ***')

for t in tunes:
    print (t.time, t.id, t.charge, t.wl)
    
print (' *** selection ***')
for t in tunes:
        #if abs(t.wl - 13.6) < 0.01 and t.charge  < 0.4:
        print ('\ntime={} id={} c={} wl={}'.format(t.time, t.id, t.charge, t.wl))
        print ('machine parameters', [p[1:] for p in db.get_machine_parameters(t.id)])
        print ('   actions:')
        for a in db.get_actions(t.id):
            print ('   id = {}, sase {}-->{}'.format(a.id, a.sase_start, a.sase_end))
            print ('     ' + str([(p.par_name, p.start_value, p.end_value) for p in db.get_action_parameters(a.tuning_id, a.id)]))