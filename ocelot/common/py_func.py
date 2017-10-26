import inspect
import os

class bcolors:
    HEADER = '\033[95m'
    BOLD = "\033[1m"
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def list_attr(obj):
    '''
    lists all attributes of an object
    '''
    for a in dir(obj):
        print(bcolors.BOLD, a, bcolors.ENDC)
        print(getattr(obj,a))

def loc_class(obj):
    '''
    locates the file where class of the object is defined
    '''
    print(inspect.getfile(obj.__class__))

def background(command):
    '''
    start command as background process
    the argument shohuld preferably be a string in triple quotes '
    '''
    import subprocess
    imports = 'from ocelot.adaptors.genesis import *; from ocelot.gui.genesis_plot import *; from ocelot.utils.xfel_utils import *; '
    subprocess.Popen(["python3", "-c", imports + command])

def copy_this_script(scriptName, scriptPath, folderPath):
    cmd = 'cp ' + scriptPath + ' ' + folderPath + 'exec_' + scriptName + ' '
    os.system(cmd)

def filename_from_path(path_string):
    # return path_string[-path_string[::-1].find(os.path.sep)::]
    return path_string.split(os.path.sep)[-1]

def deep_sim_dir(path, **kwargs):
    print(kwargs)
    from collections import OrderedDict
    import os
    if 'struct' in kwargs:
        struct = OrderedDict(kwargs['struct'])
    else:
        struct = OrderedDict([
                       ('beam_filename','%s'),
                       ('E_beam','%.1fGeV'),
                       ('beamline_no','sase%i'),
                       ('E_photon','%.0feV')
                       ])
    
    for parameter, style in struct.items():
        for key, value in kwargs.items():
            
            if key == parameter:
                path += os.path.sep + (style %(value))
    
    if 'ending' in kwargs:
        path += kwargs['ending']
        
    path += os.path.sep
    return path