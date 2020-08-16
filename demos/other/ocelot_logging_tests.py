# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 12:22:32 2019

@author: sserkez

getLogger() returns a reference to a logger instance with the specified name if it is provided, or root if not. 
The names are period-separated hierarchical structures. Multiple calls to getLogger() with the same name will return a reference 
to the same logger object. Loggers that are further down in the hierarchical list are children of loggers higher up in the list. 
For example, given a logger with a name of foo, loggers with names of foo.bar, foo.bar.baz, and foo.bam are all descendants of foo.

Loggers have a concept of effective level. If a level is not explicitly set on a logger, the level of its parent is used instead 
as its effective level. If the parent has no explicit level set, its parent is examined, and so on - all ancestors are searched 
until an explicitly set level is found. The root logger always has an explicit level set (WARNING by default). 
When deciding whether to process an event, the effective level of the logger is used to determine whether the event is passed 
to the logger’s handlers.

Child loggers propagate messages up to the handlers associated with their ancestor loggers. 
Because of this, it is unnecessary to define and configure handlers for all the loggers an application uses. 
It is sufficient to configure handlers for a top-level logger and create child loggers as needed. 
(You can, however, turn off propagation by setting the propagate attribute of a logger to False.)

"""

import logging

import ocelot
from ocelot.common.ocelog import *

print(__name__)

#----------------------
# SETTING UP LOGGER1

_logger1 = logging.getLogger('logger0.FirstLoggerName')
#_logger1.propagate=0
#print('Logger1 handlers before cleaning: ',_logger1.handlers)
_logger1.handlers = []
#print('Logger1 handlers after cleaning: ',_logger1.handlers)
 

   
ch1 = logging.StreamHandler()
ch1.setLevel(logging.DEBUG)
#log_console_format = logging.Formatter("[%(levelname)-8s] %(message)s")
#cf1 = log_console_format
#ch1.setFormatter(cf1)
_logger1.addHandler(ch1) #added Console handler


fh1 = logging.FileHandler('testlogger1.log')
fh1.setLevel(logging.DEBUG)
log_file_format = '%(asctime)s - [%(levelname)-8s] - %(message)s - %(name)s - %(filename)s:%(lineno)d'
ff1 = logging.Formatter(log_file_format)
fh1.setFormatter(ff1)
_logger1.addHandler(fh1) #added File handler
print('Logger1 file handler will write to {}'.format(fh1.baseFilename))
print('Logger1 handlers after adding custom handlers: ',_logger1.handlers)

#----------------------
# SETTING UP LOGGER2

_logger2 = logging.getLogger('logger0.SecondLoggerName')
#_logger2.propagate=1
_logger2.handlers = []


ch2 = logging.StreamHandler()
ch2.setLevel(logging.DEBUG)
cf2 = logging.Formatter("[%(levelname)-8s] %(message)s")
ch2.setFormatter(cf2)
_logger2.addHandler(ch2)

fh2 = logging.FileHandler('testlogger2.log')
fh2.setLevel(logging.DEBUG)
log_file_format = '%(asctime)s - [%(levelname)-8s] - %(message)s - %(name)s - %(filename)s:%(lineno)d'
ff1 = logging.Formatter(log_file_format)
fh2.setFormatter(ff1)
_logger2.addHandler(fh2)
print('Logger2 file handler will write to {}'.format(fh2.baseFilename))

#----------------------
# SETTING UP LOGGER LEVELS


# logging levels: [DEBUG, INFO, WARNING, ERROR, CRITICAL, or a number, e.g. debug is 10, info is 12 etc..]
_logger1.setLevel(logging.INFO) 
_logger2.setLevel(logging.DEBUG)
ocelog.setLevel(logging.DEBUG) #ocelog is a global logger


#----------------------
# SAMPLE OUTPUT


_logger1.critical('logger_1_CRITICAL')
_logger1.error('logger_1_ERROR')
_logger1.warning('logger_1_WARNING')
_logger1.info('logger_1_INFO')
_logger1.debug('logger_1_DEBUG')


_logger2.critical('logger_2_CRITICAL')
_logger2.error('logger_2_ERROR')
_logger2.warning('logger_2_WARNING')
_logger2.info('logger_2_INFO')
_logger2.debug('logger_2_DEBUG')


ocelog.critical('ocelog_CRITICAL')
ocelog.error('ocelog_ERROR')
ocelog.warning('ocelog_WARNING')
ocelog.info('ocelog_INFO')
ocelog.debug('ocelog_DEBUG')
ocelog.log(5, 'log level 5')
ocelog.log(15, 'log level 15')
ocelog.log(25, 'log level 25')