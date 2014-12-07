#
# File:         Makefile 
#
#		source/server/ttf/daq_server
#

DOOCSROOT = ../../
PROJECT_NAME = names_read


# to define DOOCSROOT as an absolute path

include $(DOOCSROOT)/$(DOOCSARCH)/DEFINEDOOCSROOT
include $(DOOCSROOT)/$(DOOCSARCH)/CONFIG


#__THIS IS CORRECT WAY
OBJDIR = $(DOOCSROOT)/$(DOOCSARCH)/obj/$(PKGDIR)

#__TEMPORARY HACK by Ast
#OBJDIR = $(DOOCSROOT)/bin/

SRCDIR = $(DOOCSROOT)/source/$(PKGDIR)

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


$(OBJDIR)/$(PROJECT_NAME):    $(SOURCEOBJ)	
		$(LINK.cc)  -I.  $(LDFLAGS) -o $(OBJDIR)/$(PROJECT_NAME) $(SOURCEOBJ) \
		        $(LDLIBS)   -lDOOCSapi -lldap
		@chmod g+w $(OBJDIR)/$(PROJECT_NAME)
		@echo "----------------$(PROJECT_NAME) done---------------"


clean:
	rm -f $(SOURCEOBJ)

#  -lnsl -ldl -lpthread -lm -lrt
