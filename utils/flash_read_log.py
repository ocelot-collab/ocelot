import numpy as np
from datetime import datetime
#import dateutil
import matplotlib.pyplot as plt
import matplotlib.dates as md
import matplotlib
import sys
font = {'size'   : 20}
matplotlib.rc('font', **font)


if len(sys.argv)>1:
    filename = sys.argv[1]
else:
    filename = "test.txt"

f = open(filename, "r")
line = f.readline()
f.close()

names = line.split("\t")

def read_wrong(filename):
    f = open(filename, "r")
    lines = f.readlines()
    f.close() 
    
    #print lines
    L = []
    for line in lines:
        #print line.find("Read error 123"), line
        #print "end"
        
        if line.find("Read error 123")>=0:
            #line = line.replace("Read error 123", "")
            #L.append(line)
            continue
        elif line.find("Read error 100")>=0:
            #line = line.replace("Read error 100", "")
            #L.append(line)
            continue
        elif line.find("\n")==0:
            continue
        else:
            L.append(line)
        #print line
    new_name = filename.split(".")[0]+"_mod.txt"
    print new_name
    f2 = open(new_name, "w")
    #for l in L:
    #    print l
    f2.writelines(L)
    f2.close()  
    return new_name


#filename = read_wrong(filename)

data = np.loadtxt(filename, comments = "#")
inrv = 1

ncols = np.shape(data)[1]
timestamp = data[:, 0]

#data = data[:, 1]


dates = [datetime.fromtimestamp(t) for t in timestamp]

fig, ax = plt.subplots(figsize=(20,10))
#plt.xticks( rotation=25 )

xfmt = md.DateFormatter('%H:%M:%S')
ax.xaxis.set_major_formatter(xfmt)
ax.set_ylabel(r"$W, \mu J$")
pax2 = ax.plot( dates[::inrv], data[::inrv,-1],"r-", lw=2, label = "sase")
ax2 = ax.twinx()
pict = []
for i in range(1, ncols-1):
    x = data[:,i]
    shift = np.around(x[0], decimals=2)
    
    ax2.plot( dates[::inrv], x[::inrv] - shift, label = names[i]+ ": " + str(shift) )
    #pict.append(p)
fig.autofmt_xdate()


pict.append(pax2)
#ax2.legend(loc='upper center')
ax2.legend(loc=1)
ax2.grid(True)
ax2.set_ylabel(r"$I, A$")
#ax.legend(pict, [l.get_label() for l in pict])
ax.grid()
#plt.plot(data[::inrv],'.-')
print dates[-1]
ax.set_xlim([dates[0], datetime.fromtimestamp(timestamp[-1]+120)])
plt.savefig(filename.split(".")[0]+".png")

plt.show()


