#
# File:         Makefile 
#
#		source/server/ttf/daq_server
#

DOOCSROOT = ./doocs/
PROJECT_NAME = dcs


# to define DOOCSROOT as an absolute path

include $(DOOCSROOT)/$(DOOCSARCH)/DEFINEDOOCSROOT
include $(DOOCSROOT)/$(DOOCSARCH)/CONFIG


#__THIS IS CORRECT WAY
#OBJDIR = $(DOOCSROOT)/$(DOOCSARCH)/lib/$(PKGDIR)
OBJDIR = $(DOOCSROOT)/$(DOOCSARCH)

#__TEMPORARY HACK by Ast
#OBJDIR = $(DOOCSROOT)/bin/

#SRCDIR = $(DOOCSROOT)/source/$(PKGDIR)
SRCDIR = .

SOURCEFILES = \
        $(SRCDIR)/$(PROJECT_NAME).cc 

SOURCEOBJ = \
        $(OBJDIR)/$(PROJECT_NAME).o

ALLPROGS = \
	$(OBJDIR)/$(PROJECT_NAME)

CPPFLAGS += $(ROOTCFLAGS) -I$(DOOCSINC)



all: $(ALLPROGS)

$(OBJDIR)/.depend depend:
		@if [ ! -f $(OBJDIR) ] ; then \
		  echo ---------- create dir $(OBJDIR) --------------; \
		  mkdir -p $(OBJDIR) ; \
		fi
#		$(CCDEP) -r -I. -I$(ROOTSYS)/include > $(OBJDIR)/.depend_temp
		for i in $(SRCDIR)/*.cc ;do $(CCDEP) $$i ;done >> $(OBJDIR)/.depend_temp
		cat $(OBJDIR)/.depend_temp | sed -e "/:/s/^/\$$\(OBJDIR\)\//g" > $(OBJDIR)/.depend
		chmod g+w $(OBJDIR)/.depend*

include $(OBJDIR)/.depend


#$(OBJDIR)/$(PROJECT_NAME):    $(SOURCEOBJ)	
#		$(LINK.cc)  -I.  $(LDFLAGS) -o $(OBJDIR)/$(PROJECT_NAME) $(SOURCEOBJ) \
#		        $(LDLIBS)   -lDOOCSapi -lldap
#		@chmod g+w $(OBJDIR)/$(PROJECT_NAME)
#		@echo "----------------$(PROJECT_NAME) done---------------"


swig:
	swig -c++ -python dcs.i

lib:
	g++ $(CCFLAGS) $(CPPFLAGS) -L/usr/X11/lib -I /usr/include/python2.7 -I. -I./doocs/include -I/local/lib/include $(LDFLAGS) -shared -fPIC -o _dcs.so dcs_wrap.cxx dcs.cc $(LDLIBS)   -lDOOCSapi -lldap -lrt


clean:
	rm -f $(SOURCEOBJ)

#  -lnsl -ldl -lpthread -lm -lrt
