from pylab import *
from time import sleep

class ActionGraphNode:
    def __init__(self, id):
        self.id = id
        self.next_options = []
    def next(self):
        if len(self.next_options) > 0 :
            opt_id  =  np.random.randint(0, len(self.next_options)) 
            return self.next_options [opt_id]
        else:
            return None 
    def apply(self):
        print 'applying ', self.id
        sleep(1)
        

if __name__ == "__main__":

    a = ActionGraphNode(id="a")
    b = ActionGraphNode(id="b")
    c = ActionGraphNode(id="c")        

    a.next_options = [b,c]
    b.next_options = [b,c]
    c.next_options = [b,c]

    print a.next().id

    
    state = a

    while True:
        state.apply()
        state = state.next() 

