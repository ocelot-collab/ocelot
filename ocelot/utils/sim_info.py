import numpy as np
import pickle
import time


class RunInfo:
    def __init__(self, id):
        self.id = id
        self.stage = 1
        self.power = 0.

class SimInfo:
    def __init__(self, n=10):
        self.x = np.linspace(0,1,n)
        self.y = np.sin(self.x)

        self.runs = {}
        
        
if __name__ == "__main__":
    o = SimInfo()
    id = 0
    while True:
        r = RunInfo(id)
        r.power = np.random.randint(100)
        o.runs[r.id] = r
        print ('saving run ', id)
        f_obj = open('dump.dat', 'wb')
        pickle.dump(o, f_obj)
        f_obj.close()
        id += 1
        time.sleep(2)

