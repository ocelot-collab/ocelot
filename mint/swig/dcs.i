/* File: dcs.i */
%module dcs

%{
#define SWIG_FILE_WITH_INIT
#include "dcs.hh"
%}

%include "dcs.hh"



%extend Func_1d {
    double __getitem__(int i) {
         return $self->f[i];
   }


    double __setitem__(int i, double v) {
         return $self->f[i] = v;
   }
 
    double __len__() {
         return $self->n;
   }

    
    char *__str__() {
       static char tmp[1024];
       sprintf(tmp,"[%g,%g,%g ...", $self->f[0],$self->f[1],$self->f[2]);
       return tmp;
   }
   
   Func_1d __add__(Func_1d *other) {
         Func_1d f2($self->n);
         for (int i=0; i < $self->n; i++) f2.f[i] = $self->f[i] + other->f[i];
         return f2;
    }
   

};

 %extend Orbit {
	 BPM __getitem__(int i) {
        return $self->bpms[i];
    }

 };

