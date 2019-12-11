"""
Launches simulations on various remote configurations
"""

import os
import time
import subprocess
from ocelot.common.ocelog import *

_logger = logging.getLogger(__name__)


def createId(prefix):
    t = time.localtime()
    return prefix + '_' + str(t.tm_year) + '_' + str(t.tm_mon) + '_' + str(t.tm_mday) + '_' + str(t.tm_hour) + '_' + str(t.tm_min) + '_' + str(t.tm_sec)


class Launcher:
    def __init__(self):
        pass
    
    def backup(self, dir='./'):
        newDir = createId(dir)
        os.system('cp -r ' + dir + ' ' + newDir)

    def cleanAllBut(self, dir='./',  keepFiles=[]):
        files = os.listdir(dir)
        
        for f in files:
            if f not in keepFiles:
                os.remove(dir + os.path.sep + f)
    
class LocalLauncher(Launcher):
    def __init__(self):
        self.command = ""
        self.dir = '.'
        self.id = ''

    def prepare(self):
        pass

    def outputfile(self,filename):
        self.filename = filename

    def launch(self):
        t1 = time.time()
        self.command = self.program
        print ('launching in ', self.dir)

        if hasattr(self, "filename"):
            self.output = open(self.filename,"w")
            self.output.write("running: "+self.command+"\n")
            self.output.write("working dir: "+self.dir+"\n")

        try:
            p = subprocess.Popen(self.command,
                                cwd=self.dir,
                                shell=True,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)

            #Real-time output:
            #not beautifull

            while True:
                p.poll()
                line = p.stdout.readline()
                if line:
                    print (line)
                    if hasattr(self, "output"):
                        self.output.write(line)
                    p.stdout.flush()

                if line == "":
                    if hasattr(self, "output"):
                        self.output.write("process finished\n")
                    break

            #Error handling
            stats = p.communicate()

            if stats[1] != "":
                print ("subprocess popen error:\n", stats[1])
                if hasattr(self, "filename"):
                    self.output.write("subprocess popen error:\n"+stats[1])
                exit()

            t2 = time.time()
            print ('execution time ', t2 - t1, ' sec')

        except (OSError) as err:
            print ("os error:",err.errno,err.strerror, "initial cmd is:\n", self.program)
            if hasattr(self, "output"):
                self.output.write("os error: "+err.errno+" "+err.strerror)
            exit()


    
    def collect(self, inputDir, outputDir, pattern = "*"):
        #os.system('cp ' + inputDir + pattern + ' ' + outputDir)
        pass
    

class MpiLauncher(Launcher):
    def __init__(self):
        self.id = ''
        self.program = ""
        self.argument = ""
        self.nproc = 0
        self.mpiParameters = ""
        
    def prepare(self):
        return False

    def launch(self):
        t1 = time.time()
        _logger.info('launching mpi job')
        command = 'mkdir -p '+ self.dir + '; ' + 'cd '+ self.dir + '; '+ "`which mpirun` " + str(self.mpiParameters) + " " + self.program + self.argument
        _logger.debug(ind_str + 'launcher command "{}"'.format(command))
        # print (command)
        os.system(command)
        t2 = time.time()
        _logger.debug(ind_str + 'execution time {:} sec'.format(t2 - t1))
            
    # collect data to directory
    def collect(self, inputDir, outputDir, pattern = "*"):
        #print 'collecting results'
        #os.system('cp ' + inputDir + pattern + ' ' + outputDir)
        pass

class NewLauncher(Launcher):
    '''
    for MOGA, tmp?
    '''
    def __init__(self):
        self.id = ''
        self.host = 'localhost'
        self.program = ""
        self.nproc = 0
        self.mpiParameters = ""
        
    def prepare(self):
        return False

    def launch(self):
        t1 = time.time()
        print ('launching job')
        command = 'mkdir -p '+ self.dir + '; ' + 'cd '+ self.dir + '; ' + self.program
        print (command)
        os.system(command)
        t2 = time.time()
        print ('execution time ', t2 - t1, ' sec')
            
    # collect data to directory
    def collect(self, inputDir, outputDir, pattern = "*"):
        #print 'collecting results'
        #os.system('cp ' + inputDir + pattern + ' ' + outputDir)
        pass
        
class PollableMpiLauncher(Launcher):
    def __init__(self):
        self.id = ''
        self.host = 'localhost'
        self.program = ""
        self.nproc = 0
        self.mpiParameters = ""

    def prepare(self):
        return False

    def launch(self):
        t1 = time.time()
        print ('launching mpi job')

        ##THIS LINE AND THE NEXT CHANGED BY GG
        # ## command = 'mkdir -p '+ self.dir + '; ' + 'cd '+ self.dir + '; '+ "mpirun " + str(self.mpiParameters) + " -n " + str(self.nproc) + " " + self.program
        command = 'mkdir -p '+ self.dir + '; ' + 'cd '+ self.dir + '; '+ "`which mpirun` " + str(self.mpiParameters) + " " + self.program
        print (command)
        os.system(command)
        t2 = time.time()
        print ('execution time ', t2 - t1, ' sec')
            
    # collect data to directory
    def collect(self, inputDir, outputDir, pattern = "*"):
        #print ('collecting results')
        #os.system('cp ' + inputDir + pattern + ' ' + outputDir)
        pass

class SshMpiLauncher(Launcher):
    def __init__(self):
        self.id = ''
        self.host = 'localhost'
        self.program = ""
        self.nproc = 0
        self.input = ["tmp.gen","tmp.cmd"]
        
    def prepare(self):
        
        os.system('mkdir -p job_sandbox')
        for f in self.input:
            os.system('cp ' + f + ' job_sandbox')
        
        command1 = "ssh " + self.port + " " + self.host + " 'mkdir -p " + self.dir + "'"
        command2 = "ssh " + self.port + " " + self.host + " 'rm -f " + self.dir + "/*'"
        command3 = "scp -r " + self.port.replace("-p","-P") + " job_sandbox/* " + self.host + ":" + self.dir

        print (command1)
        print (command2)
        print (command3)
        os.system(command1)
        os.system(command2)
        os.system(command3)

    
    def launch(self):
        t1 = time.time()
        print ('launching remote mpi job')
        command = "ssh " + self.port + " " + self.host + " 'mkdir -p " + self.dir + "; cd " + self.dir + "; mpirun " + str(self.mpiParameters) + " -n " + str(self.nproc) + " " + self.program + "'"
        print (command)
        os.system(command)
        t2 = time.time()
        print ('execution time ', t2 - t1, ' sec')

    def collect(self, inputDir, outputDir, pattern = "*"):
        print ('collecting results from remote host')
        os.system("mkdir -p " + outputDir)
        
        command1 = "ssh " + self.port + " " + self.host + " 'cd " + self.dir + "; python merge.py simulation.gout.slice >> simulation.gout'"
        command2 = "scp -r " + self.port.replace("-p","-P") + " " + self.host + ":" + inputDir + "/" + pattern + " " + outputDir 

        os.system(command1 + " ; " + command2)
        
            
