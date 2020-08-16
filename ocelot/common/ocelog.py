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
from copy import copy
import logging
from logging import Formatter
import inspect

# import traceback
# len(traceback.extract_stack())

ind_str = ': '



_log_colored = True
_log_indented = True
_log_debugging = True


_MAPPING = {
'INFO'   : '0', # default
# 'INFO'   : '0;37', # white
'DEBUG'    : '0;36', # cyan
'WARNING' : '1;33', # yellow
'ERROR'   : '0;31', # red
'CRITICAL': '0;41', # white on red bg
}
#yellow: print('\033[33myellow\033[0m')
_PREFIX = '\033['
_SUFFIX = '\033[0m'

# Create basic OCELOT logger. All loggers with name ocelot.something.something. ... inherit/affected from basic logger
ocelog = logging.getLogger('ocelot')
#ocelog.propagate = 1    # do not affect root logger
ocelog.indent0 = len(inspect.stack())

class OcelogFormatter(Formatter):
    
    
    def __init__(self, patern):
        Formatter.__init__(self, patern)
        
    def format(self, record):
        fmt_orig = self._style._fmt
        ocelog_record = copy(record)
        
        if _log_colored:
            seq = _MAPPING.get(ocelog_record.levelname, '0') # default
            ocelog_record.msg = ('{0}{1}m{2}{3}').format(_PREFIX, seq, ocelog_record.msg, _SUFFIX)
            # ocelog_record.levelname = ('{0}{1}m{2}{3}').format(_PREFIX, seq, ocelog_record.levelname, _SUFFIX)
            
        if _log_indented:
            # print('stack ',len(inspect.stack()))
            # print('_indent0 ', _indent0)
            # if hasattr(ocelog, 'indent0'):
                # print('ocelog.indent0', ocelog.indent0)
            # print('logging.indent0_before = ' + str(logging.indent0))
            indent = len(inspect.stack())
            if indent < ocelog.indent0:
                ocelog.indent0 = indent
            # print('indent = ' + str(indent - logging.indent0))
            # print('logging.indent0_after = ' + str(logging.indent0))
            ind_space = ind_str * (indent - ocelog.indent0)
            ocelog_record.msg = ind_space + ocelog_record.msg
            
        if _log_debugging:
            if ocelog_record.levelname != 'INFO':
                self._style._fmt += ' \033[37m(%(filename)s:%(lineno)d)\033[0m'
        
        result = Formatter.format(self, ocelog_record)
        self._style._fmt = fmt_orig
        return result

def ocelog_indentate():
    ocelog.indent0 = len(inspect.stack())
    print('ocelog.indent0.init', ocelog.indent0)

# logging.basicConfig(stream=sys.stdout) #test


# _console_format = "[%(name)s][%(levelname)s]  %(message)s (%(filename)s:%(lineno)d)"
# _file_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

# _console_format = "[%(levelname)s]  %(message)s [%(name)s] \033[37m(%(filename)s:%(lineno)d)\033[0m"

# ocelog.console_format = "[%(levelname)-8s]  %(message)s (%(filename)s:%(lineno)d)  [%(name)s]" # full
# ocelog.console_format = "[%(levelname)-8s] %(message)s [%(name)s]" # with name
ocelog.console_format = "[%(levelname)-8s] %(message)s" # minimum
ocelog.file_format = '%(asctime)s - [%(levelname)-8s] - %(message)s - %(name)s - %(filename)s:%(lineno)d'

ocelog.handlers=[]
# Add console handler
if True:
    ch = logging.StreamHandler()
    # ch.setLevel(logging.DEBUG)
    cf = OcelogFormatter(ocelog.console_format)
    ch.setFormatter(cf)
    ocelog.addHandler(ch)


if False:
    # Add file handler
    fh = logging.FileHandler('ocelot.log')
    # fh.setLevel(logging.DEBUG)
    ff = logging.Formatter(ocelog.file_format)
    fh.setFormatter(ff)
    ocelog.addHandler(fh)
    
# Set log level
ocelog.setLevel(logging.INFO)