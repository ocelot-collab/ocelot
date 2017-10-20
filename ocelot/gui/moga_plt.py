import matplotlib.pyplot as plt
import pickle
import datetime, time
import threading


class ParetoPlot():

    def __init__(self):
       
       self.work = False

       self.limx = None
       self.limy = None
       self.file = 'moga.dat'

       self.x = []
       self.y = []
       self.xnd = []
       self.ynd = []

       self.title = ''
       self.xlabel = 'x variable'
       self.ylabel = 'y variable'

       self.fig, self.ax = plt.subplots()
       self.lineF, = self.ax.plot(self.x, self.y, 'ro', ms=3)
       self.lineND, = self.ax.plot(self.xnd, self.ynd, 'bo', ms=4, alpha=0.75)
       self.ax.set_xlabel("f1")
       self.ax.set_xlabel("f2")
       self.ax.grid(True)

       self.update_time = 1000

    def get_data_once(self):
        t_start = datetime.datetime.now()

        try:
            with open(self.file, 'rb') as f:
                data_file = pickle.load(f)
        except Exception:
            time.sleep(0.1)
            return

        self.title = 'Iteration ' + str(data_file[0][0]) + ' from ' + str(data_file[0][1])

        self.x = []
        self.y = []
        for ind in data_file[1]:
            self.x.append(ind[0])
            self.y.append(ind[1])

        if self.limx == None:
            self.limx = [0, max(self.x)]
        if self.limy == None:
            self.limy = [0, max(self.y)]

        self.xnd = []
        self.ynd = []
        for ind in data_file[2]:
            self.xnd.append(ind[0])
            self.ynd.append(ind[1])

        sleep = (datetime.datetime.now().microsecond - t_start.microsecond) / 1000.0
        if sleep < 0.0: sleep += 1000.0
        sleep = (self.update_time - sleep) / 1000.0 - 0.001  # 0.001 - correction factor
        if sleep > 0.0: time.sleep(sleep)

    def get_data(self):
        
        while(self.work):
            self.get_data_once()


    def plot(self):

        plt.title(self.title)

        self.ax.set_xlim(self.limx)
        self.ax.set_ylim(self.limy)

        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)

        self.lineF.set_data(self.x, self.y)
        self.lineND.set_data(self.xnd, self.ynd)
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def show_once(self):

        #self.work = True

        self.get_data_once()

        self.plot()
        plt.show()
        #self.work = False

        #thread.join()


    def run(self):

        self.work = True

        thread = threading.Thread(target=self.get_data)
        thread.start()


        timer = self.fig.canvas.new_timer(interval=self.update_time)
        timer.add_callback(self.plot)
        timer.start()
        plt.show()
        self.work = False

        thread.join()


if __name__ == "__main__":
    pic = ParetoPlot()
    pic.limx = [0.0, 20.0]
    pic.limy = [0.0, 10.0e-3]
    pic.file ="../test/moga.dat"
    pic.run()
