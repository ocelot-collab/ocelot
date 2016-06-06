download swig
http://www.swig.org/download.html

dcs.i contsins interface definition
dcs.hh and dcs.cc are the library

generate python wrappers (alice_wrap.cxx, alice.py) 
> swig -c++ -python dcs.i

compile python wrappers (linux)
> g++ -I /usr/include/python2.7/ -shared -fPIC -o _dcs.so dcs_wrap.cxx dcs.cc

some usage examples in example.py

compile main version if needed
> g++ -o dcs dcs_main.cc dcs.cc
