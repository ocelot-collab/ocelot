# for windows - important to add this file to PATH ..\Pre-built.2\lib\pthreadGC2.dll 
# otherwise xcode will not work. 
import subprocess
import shlex
import os
from sys import platform, path
import shutil


def remove_file(path_to_folder, filename):
    files = os.listdir(path_to_folder)
    if filename in files:
        os.remove(path_to_folder + "/" + filename)

def move_file(source_folder, filename, destin_folder):
    files = os.listdir(source_folder)
    if filename in files:
        shutil.move(source_folder+"/" + filename,destin_folder)

def compile(dir_path):
    os.chdir(dir_path)
    print "platform = ", platform, "  os.name = ", os.name
    if os.name == "nt":
        cmd = 'c:/MinGW/bin/mingw32-make -f Makefile'
    else:
        cmd = "make -f Makefile"

        if platform == "darwin":
            # just in case if mac OS
            cmd = "make -f Makefile"
    args = shlex.split(cmd)
    #print args
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    result = p.communicate()[0]
    print result

# clean
def clean_folder(dir_path):
    """
    it does not need
    """
    print dir_path
    os.chdir(dir_path)
    cmd = "make clean"
    if os.name == "nt":
        cmd = 'c:/MinGW/bin/mingw32-make -f Makefile'
    args = shlex.split(cmd)
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    result = p.communicate()[0]
    print result

# move lib to genera_libs
home_path = path[0]
#print "compile ", home_path
import ocelot
#print ocelot.__file__
import os
path_to_ocelot = os.path.dirname(ocelot.__file__)


#indx = path[0].find("siberia2")

#xcode_path = home_path[:indx]
gen_path = path_to_ocelot + "/lib/genera/"

libs_path = gen_path + "build/genera_libs"

cpp_path = gen_path + "src/cpp"

for dirname in ["undulator", "radiation", "convolution"]:
    dir_path = cpp_path + "/" + dirname

    compile(dir_path)

    filename = dirname +".so"

    remove_file(libs_path, filename)

    move_file(dir_path, filename, libs_path)

#shutil.move(dir_path+"/undulator.so",libs_path)

