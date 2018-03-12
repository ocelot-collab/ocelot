#!/usr/bin/env python

#  Copyright (c) 2017 Kurt Jacobson

#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:

#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.

#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.
import sys
from copy import copy
from logging import Formatter

MAPPING = {
    'DEBUG'   : '0;37', # white
    'INFO'    : '0;36', # cyan
    'WARNING' : '1;33', # yellow
    'ERROR'   : '0;31', # red
    'CRITICAL': '0;41', # white on red bg
}
#yellow: print('\033[33myellow\033[0m')
PREFIX = '\033['
SUFFIX = '\033[0m'

class ColoredFormatter(Formatter):

    def __init__(self, patern):
        Formatter.__init__(self, patern)

    def format(self, record):
        colored_record = copy(record)
        levelname = colored_record.levelname
        seq = MAPPING.get(levelname, 37) # default white
        colored_levelname = ('{0}{1}m{2}{3}') \
            .format(PREFIX, seq, levelname, SUFFIX)
        colored_record.levelname = colored_levelname
        return Formatter.format(self, colored_record)

import logging

# logging.basicConfig(stream=sys.stdout) #test


# _console_format = "[%(name)s][%(levelname)s]  %(message)s (%(filename)s:%(lineno)d)"
# _file_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

_console_format = "[%(levelname)s]  %(message)s [%(name)s] \033[37m(%(filename)s:%(lineno)d)\033[0m"
_file_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s - %(filename)s:%(lineno)d'
_logtofile = False

ocelog = logging.getLogger('ocelot')
ocelog.handlers=[]
# Add console handler

if True:
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    cf = ColoredFormatter(_console_format)
    ch.setFormatter(cf)
    ocelog.addHandler(ch)


if _logtofile:
    # Add file handler
    fh = logging.FileHandler('ocelot.log')
    fh.setLevel(logging.DEBUG)
    ff = logging.Formatter(_file_format)
    fh.setFormatter(ff)
    ocelog.addHandler(fh)
    
# Set log level
ocelog.setLevel(logging.INFO)







# # to be deprecated soon?
# from ocelot.common.py_func import bcolors

# class Logger:
    # '''
    # use logger instead of print statements throughout the code
    # for debug and other purposes
    # a different implementation of same interface can be assigned to ocelot.logger if other features are needed
    # '''
        
    # def __init__(self):

        # self.show_info = True
        # self.show_warning = False
        # self.show_debug = False
        # self.file = None

    # def info(self, txt):
        # if self.show_info: print(txt)

    # def warn(self, txt):
        # if self.show_warning: print(bcolors.WARNING + txt + bcolors.ENDC)

    # def debug(self, txt):
        # if self.show_debug: print(txt)



