import numpy as np
import datetime
import time
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
        if name == "time":
            col = col.astype(int)
        dict_data[name] = col

    return dict_data


def read_orbit_log(filename):
    f = open(filename, "r")
    line = f.readline()
    line = line.replace("#", "")
    line = line.replace(" ", "")
    line = line.replace("\n","")

    names = line.split("\t")
    names.remove("plane")
    n_n = len(names)
    #content = f.readlines()
    lines = [line.rstrip('\n') for line in f]
    f.close()

    x_plane = []
    y_plane = []
    for line in lines:
        data_line = line.split("\t")
        if data_line[-1] == "X":
            x = np.array(data_line[:-1]).astype(float)
            x_plane.append(x)
        else:
            y = np.array(data_line[:-1]).astype(float)
            y_plane.append(y)

    #x_plane = np.array(x_plane)
    x_plane = np.array([item for sublist in x_plane for item in sublist])
    nx = len(x_plane)
    x_plane.reshape((nx/n_n, n_n))
    dict_data_x = {}
    for name, col in zip(names, np.transpose(x_plane)):
        if name == "time":
            col = col.astype(int)
        dict_data_x[name] = col

    y_plane = np.array([item for sublist in y_plane for item in sublist])
    ny = len(y_plane)
    y_plane.reshape((ny/n_n, n_n))
    dict_data_y = {}
    for name, col in zip(names, np.transpose(y_plane)):
        if name == "time":
            col = col.astype(int)
        dict_data_y[name] = col

    return dict_data_x, dict_data_y


def plot_dict(dict_data, filename=None, interval=1, mode="%", grid=True):
    inrv = interval
    times = [datetime.datetime.fromtimestamp(t) for t in dict_data["time"]]
    #print (times)
    devices = list(dict_data.keys())
    devices.remove("time")
    devices.remove("sase")
    #print devices
    fig, ax = plt.subplots(figsize=(20, 10))
    #fig, ax = plt.subplots()

    xfmt = md.DateFormatter('%H:%M:%S')
    ax.xaxis.set_major_formatter(xfmt)
    ax.set_ylabel(r"U [$\mu$ J]")
    pax2 = ax.plot(times[::inrv],  dict_data["sase"][::inrv],"r-", lw=4, label = "sase")
    ax.grid(grid)

    ax.set_xlim([times[0], datetime.datetime.fromtimestamp(dict_data["time"][-1])])
    ax.legend(loc=1, framealpha=0.7)
    #ax.legend(loc=1)
    ax2 = ax.twinx()
    pict = []
    #devices=['V3DBC3', 'H3DBC3', 'H10SMATCH', 'V14SMATCH', 'H12SMATCH', 'V7SMATCH']
    lines = []
    for device in devices:
        x = dict_data[device]
        shift = np.around(x[0], decimals=2)
        if mode == "%":
            line, = ax2.plot( times[::inrv], (x[::inrv] - x[0])/x[0], lw = 2, label = device)
        else:
            line, = ax2.plot( times[::inrv], x[::inrv] - shift,lw=2,  label = device + ": " + str(shift) )
        lines.append(line)
    fig.autofmt_xdate()
    pict.append(pax2)
    ax2.legend(loc=4, framealpha=0.5)
    # Create a legend for the first line.
    # first_legend = plt.legend(handles=[lines[0], lines[1]], loc=3, framealpha=0.5)
    # Add the legend manually to the current Axes.
    # ax2.add_artist(first_legend)
    # ax2.legend(handles=[lines[2], lines[3], lines[4], lines[5]], loc=4,framealpha=0.5)
    #ax2.legend(loc=4)
    ax2.grid(grid)
    if mode == "%":
        ax2.set_ylabel(r"$\Delta I/I$")
    else:
        #ax2.set_ylabel(r"$I, A$")
        ax2.set_ylabel(mode)
    if filename != None:
        plt.savefig(filename.split(".")[0]+".png")
    #fig.set_size_inches(20, 10)
    fig.savefig("figure_2.eps", format="eps")
    plt.show()

def timslot_extract(dict_data, timeslot):
    #d0 = datetime.datetime.fromtimestamp(dict_data["time"][0])
    times = [datetime.datetime.fromtimestamp(t) for t in dict_data["time"]]
    start = timeslot[0].split(":")
    end = timeslot[1].split(":")
    print(start, end)
    ind_0 = 0
    ind_1 = -1
    for i, t in enumerate(times):
        if t.hour == int(start[0]) and t.minute == int(start[1]):
            ind_0 = i
        if t.hour == int(end[0]) and t.minute == int(end[1]):
            ind_1 = i
            break
    for name in dict_data.keys():
        dict_data[name] = dict_data[name][ind_0:ind_1]
    return dict_data

def plot_log(filename, devices=None, timeslot=None, mode="%", grid=True):

    dict_data = read_log(filename)
    if timeslot is not None:
        dict_data = timslot_extract(dict_data, timeslot)

    if devices is not None:
        dict_data_new = {}

        dict_data_new["sase"] = dict_data["sase"]
        dict_data_new["time"] = dict_data["time"]
        for name in devices:
            dict_data_new[name] = dict_data[name]
        dict_data = dict_data_new
    plot_dict(dict_data, filename=filename, interval=1, mode=mode, grid=grid)


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
        filename = "../../desy/flash/exp_files/new_opt_4.txt"
        #filename = "../../desy/flash/exp_files/optim.txt"
        #filename = "../../desy/flash/exp_files_1_02/cor_sase_shif.txt"
        #filename = "../../desy/flash/exp_files_1_02/cor_sase_shif.txt"
    #filename = read_wrong(filename)
    #print filename
    #, "bda_x", "bda_y"
    #"tun_x", "tun_y"
    #plot_log(filename, devices=["bda_x", "bda_y"], timeslot=["01:52", "02:02"], mode="units")
    plot_log(filename, devices=['V3DBC3', 'H3DBC3', 'H10SMATCH', 'H12SMATCH', 'V14SMATCH', 'V7SMATCH'], timeslot=None, grid=False)