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

class Logger:
    '''
    use logger instead of print statements throughout the code
    for debug and other purposes
    a different implementation of same interface can be assigned to ocelot.logger if other features are needed
    '''
        
    def __init__(self):

        self.show_info = True
        self.show_warning = False
        self.show_debug = False
        self.file = None

    def info(self, txt):
        if self.show_info: print(txt)

    def warn(self, txt):
        if self.show_warning: print(bcolors.WARNING + txt + bcolors.ENDC)

    def debug(self, txt):
        if self.show_debug: print(txt)