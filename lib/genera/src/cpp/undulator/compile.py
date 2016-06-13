# for windows - important to add this file to PATH ..\Pre-built.2\lib\pthreadGC2.dll 
# otherwise xcode will not work. 
import subprocess
import shlex
import os
from sys import platform 


print "platform = ", platform, "  os.name = ", os.name
if os.name == "nt":
	cmd = 'c:/MinGW/bin/mingw32-make -f Makefile'
else:
	cmd = "make -f Makefile.ubu"
	if platform == "darwin":
		cmd = "make -f Makefile"
args = shlex.split(cmd)

p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
result = p.communicate()[0]
print result
cmd = "make clean"
if os.name == "nt":
	cmd = 'c:/MinGW/bin/mingw32-make -f Makefile'
args = shlex.split(cmd)