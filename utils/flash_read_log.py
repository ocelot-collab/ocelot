import numpy as np
from datetime import datetime
#import dateutil
import matplotlib.pyplot as plt
import matplotlib.dates as md
import matplotlib
import sys
font = {'size'   : 20}
matplotlib.rc('font', **font)



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
    print (new_name)
    f2 = open(new_name, "w")
    #for l in L:
    #    print l
    f2.writelines(L)
    f2.close()  
    return new_name


def read_log(filename):
    f = open(filename, "r")
    line = f.readline()
    line = line.replace("#", "")
    line = line.replace(" ", "")
    line = line.replace("\n","")
    f.close()
    names = line.split("\t")
    data = np.loadtxt(filename, comments="#")

    dict_data = {}
    for name, col in zip(names, np.transpose(data)):
        dict_data[name] = col

    return dict_data


def plot_dict(dict_data, filename=None, interval=1, mode="%"):
    inrv = interval
    times = [datetime.fromtimestamp(t) for t in dict_data["time"]]
    devices = list(dict_data.keys())
    devices.remove("time")
    devices.remove("sase")
    #print devices
    fig, ax = plt.subplots(figsize=(20, 10))
    #fig, ax = plt.subplots()

    xfmt = md.DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_formatter(xfmt)
    ax.set_ylabel(r"$U,\, \mu J$")
    pax2 = ax.plot(times[::inrv],  dict_data["sase"][::inrv],"r-", lw=2, label = "sase")
    ax.grid()

    ax.set_xlim([times[0], datetime.fromtimestamp(dict_data["time"][-1])])
    #ax.legend(loc=1, framealpha=0.7)
    ax.legend(loc=1)
    ax2 = ax.twinx()
    pict = []
    for device in devices:
        x = dict_data[device]
        shift = np.around(x[0], decimals=2)
        if mode == "%":
            ax2.plot( times[::inrv], (x[::inrv] - x[0])/x[0], label = device)
        else:
            ax2.plot( times[::inrv], x[::inrv] - shift, label = device + ": " + str(shift) )
    fig.autofmt_xdate()
    pict.append(pax2)
    #ax2.legend(loc=4, framealpha=0.7)
    ax2.legend(loc=4)
    ax2.grid(True)
    if mode == "%":
        ax2.set_ylabel(r"$\Delta I/I$")
    else:
        ax2.set_ylabel(r"$I, A$")
    if filename != None:
        plt.savefig(filename.split(".")[0]+".png")
    #fig.set_size_inches(20, 10)
    #fig.savefig("optim_vect_pict.svg", format="svg")
    plt.show()

def plot_log(filename):
    dict_data = read_log(filename)
    plot_dict(dict_data,filename=filename, interval=1)


def rm_nonwork_devices(dict_data, threshold=0.01, debug=False, rm_devices=["",]):

    new_dict = {}
    for name in dict_data.keys():
        if name in rm_devices:
            continue
        x = dict_data[name][:int(len(dict_data["time"])*1)]
        delta = abs((max(x) - min(x))/max(x))
        if delta > threshold:
            new_dict[name] = x
        if name == "time":
            new_dict[name] = x

    if debug == True:
        #print( new_dict.keys() )
        plot_dict(new_dict)

        devices = list(new_dict.keys())
        devices.remove("sase")
        devices.remove("time")
        print( len(devices))
        for name in devices:
            n = len(dict_data[name])
            x = new_dict[name]
            delta = abs((max(x) - min(x))/max(x))
            print(name, "delta = ", delta, "  A0 = ", x[0])
            plt.plot((x-x[0])/x[0], label=name)
        plt.legend()
        plt.show()
    return new_dict


def save_new_dict(new_dict, filename):
    names = list(new_dict.keys())
    ncols = len(names)
    nrows = len(new_dict[names[0]])
    data = np.zeros((nrows, ncols))
    data[:,0] = new_dict["time"]
    data[:,-1] = new_dict["sase"]
    names.remove("sase")
    names.remove("time")
    header = "time"
    for i, name in enumerate(names):
        data[:, i+1] = new_dict[name]
        header = header+ "\t" + name
    header = header + "\t" + "sase"
    np.savetxt(filename, data, header=header)


if __name__ == "__main__":
    
    if len(sys.argv)>1:
        filename = sys.argv[1]
    else:
        filename = "opt_3.txt"
    #filename = read_wrong(filename)
    plot_log(filename)