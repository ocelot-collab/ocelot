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
        if hasattr(self, 'command'):
            command = self.command #command override
        else:
            command = 'mkdir -p '+ self.dir + '; ' + 'cd '+ self.dir + '; '+ "`which mpirun` " + str(self.mpiParameters) + " " + self.program + " " + self.argument
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

######################################################################
#
# Current issues with MpiQueueLauncher:
# * Need to check that self.dir is valid directory with permissions
#   enabling write (or that directory can at least created under this name).
#   If needed, create directory not using 'mkdir' shell command but
#   using Python-embedded mkdir (should result in better error handling).
# * When update_job_state is called immediately after submitting the job,
#   it does not work: It prints the line
#   "unable to extract job status from 'sacct' output, exit code was: 0".

class MpiQueueLauncher(Launcher):
    def __init__(self):
        self.id = ''
        self.program = "" # old
        self.argument = "" # old
        self.dir = ""
        self.nproc = 0
        self.mpiParameters = "" # old

        self.jobpartition='maxwell'
        self.jobnodes=2
        self.jobntaskspernode=32
        self.jobname=None
        self.launch_wait = True
        ### params of running job
        self.jobid=-1        # assigned on job submission
        self.jobid_valid=False
        self.jobstate=''     # updated in 'update_job_state'
        self.jobtime=''      # updated in 'update_job_state'


    def prepare(self):
        return False


    def submit(self):
        _logger.info('preparing job')

        # deploy SLURM job script
        # generate output dir for stdout and stderr
        # delete any flag files that might still be there for an earlier run (for example, the previous stage in a series of stages)
        jobscript=os.path.dirname(__file__)+'/T_run_g2.sh'
        setupcmd='cp "'+jobscript+'" "'+self.dir +'/run_g2.sh"; mkdir -p '+self.dir+'/jobout; rm -f '+self.dir+'/flag.start '+self.dir+'/flag.finish'
        _logger.debug(ind_str + 'setup command "{}"'.format(setupcmd))
        os.system(setupcmd)
        #
        # assemble cmd line
        # !no quotes around strings, the submission of the args as list elements ensures that for instance strings with spaces in them are treated as single arguments!
        _logger.info('submitting job to MAXWELL cluster queue 1/2 (assembling cmd line)')
        args=['sbatch', '--parsable']
        if self.jobname is not None:
            args.append('--job-name={0:s}'.format(str(self.jobname)))
        args.append('--partition={0:s}'.format(self.jobpartition))
        args.append('--nodes={0:d}'.format(self.jobnodes))
        args.append('--ntasks-per-node={0:d}'.format(self.jobntaskspernode))
        args.append('./run_g2.sh')
        _logger.info('submitting job to MAXWELL cluster queue 2/2 (submitting)')
        _logger.debug('running command: {0:s}'.format(" ".join(args)))
        p=subprocess.Popen(args, cwd=self.dir, stdout=subprocess.PIPE)
        (out,err)=p.communicate()
        rc=p.returncode # get exit code of 'sbatch' call (but not used as of now)
        _logger.debug('exit code of command: {0}'.format(str(rc)))
        # print("output:", out)
        # As of 2021-Jan, success job submission on maxwell will return jobid in format: b'6759516\n'
        try:
            tstr=out.rstrip()
            self.jobid=int(tstr)  # retain job id for future reference
            self.jobid_valid=True
            _logger.info('job submitted: jobid={0}, sbatch returned status {1}'.format(str(self.jobid),str(rc)))
        except ValueError:
            # unexpected output (format) from sbatch (or on error: no output at all)
            self.jobid_valid=False
            _logger.info('unclear if job submission successful, sbatch returned status {0}'.format(str(rc)))
            return False # job submission not successful
        #
        return True # true => job submission successful

####
    def update_job_state(self):
        dbg=False
        if(self.jobid_valid==False):
            print('error: update_job_state: jobid is not valid. Was job submission successful?')
            return False

        p = subprocess.Popen(['sacct', '--parsable2', '--noheader', '--format', 'JobID,State,Partition,AllocCPUs,Elapsed', '-j', str(self.jobid)], stdout=subprocess.PIPE)

        # process output of 'sacct' command (typically multiple lines)
        gotit=False
        while (gotit==False):
            # example for typical format: l='6761134|COMPLETED|maxwell|80|00:02:48'
            l=p.stdout.readline()
            if not l:
                break

            l=l.decode()
            l=l.rstrip()
            print('update_job_state dbg: considering acct output line "{0}"'.format(l))
            ls=l.split('|')

            # expecting 5 fields
            nele=len(ls)
            if nele!=5:
                print('error: was expecting 5 fields from sacct output')
                break # return

            my_jobstate=ls[1]
            my_jobtime=ls[4]
            part=ls[2]
            # only one of the lines from 'sacct' cmd has partition info => this one we want, skip all others
            if not part:
                if dbg:
                    print('no partition id found, skipping to next line')
                continue
            gotit=True

        p.wait(timeout=10) # exit code is only available after termination of sub-process
        rc=p.returncode # get exit code of 'sacct' call (not used as of now)
        if gotit:
            if dbg:
                print('jobtime:'+my_jobtime+', jobstate:'+my_jobstate)
            self.jobtime=my_jobtime
            self.jobstate=my_jobstate
        else:
            print("unable to extract job status from 'sacct' output, exit code was: {0}".format(rc))


    def is_complete(self):
        # 1) query scheduler and update variables in class (status codes aren't evaluated yet)
        self.update_job_state()

        # 2) check the status file generated in the job script
        flagfile=self.dir+'/flag.finish'
        if (os.path.isfile(flagfile)==True):
            return True
        return False


    def wait_for_completion(self,verbose=False):
        flagfile=self.dir+'/flag.finish'
        _logger.info('waiting for job completion, polling file "{}"'.format(flagfile))
        while True:
            is_complete=self.is_complete()
            if verbose:
                print('job {0}: state={1}, runtime={2}, is_complete={3}'.format(self.jobid, self.jobstate, self.jobtime, is_complete))
            if(is_complete):
                break
            time.sleep(10) # <== won't wait if the job is already complete
####
    # Simple interface if you don't plan to run multiple job in parallel
    # Also to support code built around MpiLauncher()
    def launch(self):
        t1 = time.time()
        self.submit()
        if self.launch_wait==False:
            return
        self.wait_for_completion()
        t2 = time.time()
        _logger.debug(ind_str + 'execution time {:} sec'.format(t2 - t1))
            
    # collect data to directory
    def collect(self, inputDir, outputDir, pattern = "*"):
        #print 'collecting results'
        #os.system('cp ' + inputDir + pattern + ' ' + outputDir)
        pass

######################################################################


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
        
            
