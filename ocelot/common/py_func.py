import inspect
import os

python_v = "python3"

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
    subprocess.Popen([python_v, "-c", imports + command])
    
def copy_this_script(file, folderPath):
    scriptName = os.path.basename(file)
    scriptPath = os.path.realpath(file)
    cmd = 'cp ' + scriptPath + ' ' + folderPath + 'exec_' + scriptName + ' '
    os.system(cmd)

def filename_from_path(path_string):
    tree = path_string.split(os.path.sep)
    return tree[-1]
    
def directory_from_path(path_string):
    tree = path_string.split(os.path.sep)
    return os.path.sep.join(tree[:-1]) + os.path.sep

def deep_sim_dir(path, **kwargs):
    from collections import OrderedDict
    import os
    pathsep = os.path.sep
    path_g = ''  #generated path
    
    if not path.endswith(pathsep):
        path += pathsep
    
    if 'struct' in kwargs:
        struct = OrderedDict(kwargs['struct'])
    else:
        struct = OrderedDict([
                       ('beam_filename','%s'),
                       ('E_beam','%.1fGeV'),
                       ('beamline','%s'),
                       ('E_photon','%.0feV')
                       ])

    separator = '__'
    if 'nested' in kwargs:
        if kwargs['nested'] == True:
            separator = pathsep
            
    if 'separator' in kwargs:
            separator = kwargs['separator']
    
    for parameter, style in struct.items():
        for key, value in kwargs.items():
            if key == parameter:
                path_g += separator + (style %(value))
    
    if 'ending' in kwargs:
        path_end = kwargs['ending']
    else:
        path_end = ''
    
    path = (path + pathsep + path_g + separator + path_end + pathsep).replace(pathsep + separator, pathsep)
    
    return path