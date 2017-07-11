import inspect

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