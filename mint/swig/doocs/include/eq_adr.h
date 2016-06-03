/* file eq_adr.h
//
// EqAdr class
//
// Kay Rehlich, -PVAK-
//
// last update:
// 13. Apr. 1993
*/

/*
 *
 * $Date: 2012-02-09 13:13:11 $
 * $Source: /doocs/doocssvr1/cvsroot/source/clientlib/eq_adr.h,v $
 * $Revision: 1.27 $
 * $State: Exp $
 *
 * $Log: eq_adr.h,v $
 * Revision 1.27  2012-02-09 13:13:11  arthura
 * *** empty log message ***
 *
 * Revision 1.26  2011/12/06 08:31:52  arthura
 * *** empty log message ***
 *
 * Revision 1.25  2010-01-06 08:15:26  arthura
 * new clientlib 5.8.1
 *
 * Revision 1.24  2009/08/10 06:55:41  arthura
 * *** empty log message ***
 *
 * Revision 1.23  2009/06/22 07:52:45  arthura
 * *** empty log message ***
 *
 * Revision 1.22  2008/05/07 10:36:14  arthura
 * extended name fields for RPC protocol data structures
 *
 * Revision 1.21  2008/04/17 07:23:20  arthura
 * *** empty log message ***
 *
 * Revision 1.20  2006/09/14 12:10:09  arthura
 * introduces a "cycling".
 *
 * Revision 1.19  2005/03/04 14:52:38  arthura
 * *** empty log message ***
 *
 * Revision 1.17  2004/02/03 10:19:33  arthura
 * fixed the bug in EqAdr class
 *
 * Revision 1.16  2003/11/07 16:19:48  arthura
 * modified EqAdr (eq_adr.cc, eq_adr.h) and eq_res (eq_res.h, eq_res.cc) classes
 *
 * Revision 1.15  2003/06/27 12:48:44  arthura
 * fixed the crash in rpc call
 *
 * Revision 1.13  2003/06/10 08:29:58  arthura
 * first MT clientlib
 *
 * Revision 1.12  2003/03/11 16:52:29  rehlich
 * new in eq_adr merge method: $f $d $l $p variables to insert the facility .. property in the result addr
 *
 * Revision 1.11  2000/03/07 11:19:25  grygiel
 * change toupper to To_Upper for GNU compatibility
 *
 * Revision 1.10  1999/09/21 16:26:56  rehlich
 * 1. TINE implementation (monitors emulated)
 *
 * Revision 1.9  1999/08/02 18:27:54  grygiel
 * changes for libc6 - GNU C Library
 *
 * Revision 1.8  1998/02/12 09:44:41  rehlich
 * longer EPICS timeout, longer delay if first connect failed
 *
 * Revision 1.7  1997/05/07 14:49:54  rehlich
 * destructor and set_name_string bug fixed
 *
 * Revision 1.6  1997/05/02 07:00:21  rehlich
 * environment vari NOEPICS, update of monitor in get call, bug fix in get_locations(eq_services)
 *
 * Revision 1.5  1996/01/11 15:00:24  rehlich
 * EqAdr bug in merge fixed
 *
 * Revision 1.4  1996/01/04 17:33:42  rehlich
 * *** empty log message ***
 *
 * Revision 1.3  1995/11/13 17:00:41  grygiel
 * Test of cvs
 *
 * Revision 1.2  1995/11/12 16:55:25  cobraadm
 * modificationas for new gnu compilers
 *
 * Revision 1.1.1.1  1995/11/03 12:23:39  grygiel
 * DOOCS sources
 *
 * Revision 1.1  1995/10/24 13:04:18  xgrygiel
 * Initial revision
 *
 *
 */

#ifndef eq_adr_h
#define eq_adr_h


#include "eq_types.h"
#include <strings.h>


#define PROTOCOL_LENGTH     8
#define HOSTNAME_LENGTH     32

#define ADDR_STRING_LENGTH  (FACILITY_STRING_LENGTH + \
                             DEVICE_STRING_LENGTH   + \
                             LOCATION_STRING_LENGTH + \
                             PROPERTY_STRING_LENGTH)
enum {

        ADDR,
        ADDR_FAC,
        ADDR_DEV,
        ADDR_LOC,
        ADDR_PROP,
        ADDR_HOST,
        ADDR_PROT,
        ADDR_DOM
};

class   eq_res;



class	EqAdr    {

private:
	    // holds address parts:

        char            prot  [PROTOCOL_LENGTH        + 1];
        char            host  [HOSTNAME_LENGTH        + 1];
        char	        fac   [FACILITY_STRING_LENGTH + 1];
        char            dev   [DEVICE_STRING_LENGTH   + 1];
        char            loc   [LOCATION_STRING_LENGTH + 1];
        char            prop  [PROPERTY_STRING_LENGTH + 1];

        char            alloc_nsp;     // object has allocated the nsp space
        void            *handle;       // communication object handle
        struct timeval  list_id;       // when handle was created
        EqNameString    *nsp;	       // points to 4 name parts
        int             cycle;         // cycle mask

        void            clean           (int = 0);
        void            parse           (void);
        int             stringlen       (char *);
        void            to_upper        (char *);
        void            set_handle      (void *, struct timeval * = 0);
        void            *get_handle     (struct timeval * = 0);
        char            *get_domain     (char *);

        friend class eq_res;

public:
                       EqAdr           (EqNameString *);
                       EqAdr           (void);
                       ~EqAdr          (void);

        EqAdr          &operator =     (EqAdr *);
        EqAdr          &operator =     (EqAdr &);

        EqNameString   *name_string    (void);
        char           *protocol       (void);
        char           *hostname       (void);
        char           *facility       (void);
        char           *device         (void);
        char           *location       (void);
        char           *property       (void);

        void           set_name_string (EqNameString *);
        void           set_protocol    (char *);
        void           set_host        (char *);
        void           set_facility    (char *);
        void           set_device      (char *);
        void           set_location    (char *);
        void           set_property    (char *);
        void           set             (char *, char *, char *, char *);

        void           show_adr        (char *, size_t);

        void           adr             (char *);
        void           merge           (char *);

        void           check_upper     (char); // convert to upper case
                                                // 0 - facility
                                                // 1 - device
                                                // 2 - location
                                                // 3 - property
        int            get_cycle       (void);
        int            get_limits      (int);
};



// -----------------------------------------------------------------------
//:DESCRIPTION
//
//.Class library to the address part of the RPC based communication.
//.An address is composed of 4 parts:
//.	1) facility name e.g. "HERA.VAC"
//.	2) device name e.g. "ION_PUMP"
//.	3) location name e.g. "WR123"
//.	4) property name e.g. "PRESSURE"
//.Normaly the facitlity and the device part
//.of the name specifies a server process. A server process is
//.a program on a CPU which implements the device specific code.
//.The translation from a facility name and device name to a CPU
//.name and server program number is implemented by a
//.call to a equipment name server.
//.The client side of the communication libraries are implemented by
//.one EqServerEntry class per server program which needs to be accessed.
//.All EqAdr's with the same server target are pointing to one EqServerEntry
//.( obj_p pointer). When the address is modified this pointer is
//.automatically recalculated if necessary.
//.
//.
//.EqAdr ( )
//.	Constructor.
//.~EqAdr();
//.	Destructor.
//.
//.char*	facility()
//.char*	device()  
//.char*	location()
//.char*	property()
//.	Return the 4 address parts as a character string.
//.

//.void adr ( char* name);
//.	To set the address. The name argument is an ASCII string of the form
//.	"facility/device/location/property".
//.
//.void merge ( char* name);
//.	Merges the address from the current address with the string "name"
//.	of the "merge" call. The existing field of the current address is
//.	used if the corresponding field in the name argument is empty.
//.	A filled field in name e.g. "//VERT/" overwrites the existing
//.	one and a string starting with a "." character is appended to
//.	the existing string.
//.	An adr("a/b//")  followed by a merge("//c/d") composes a resulting
//.	name of "a/b/c/d" followed by a merge("//.f/") gives "a/b/c.f/d".
//.	Wild cards like merge("f/d/*/p") are allowed and sets the address to
//.	"f/d/*/p".
//.     In addtion a variable substitution with a $ in name is implemented.
//.     $f is replaced with the facility name of this class, $d is the device name, 
//.     $l the location and $p for the property.
//.     Example: ("a/b/lll/ppp")  followed by a merge("x$fy//c$lxx/d$p$f") composes a 
//.     result of "xay/b/clllxx/dddda".
//!void	show_adr(char* dest, int length)
//.	Helper function to extract the address from this EqAdr instance.
//.	The result is written as an ASCII string into dest with a
//.	maximum length of length character to be filled into dest.
//.	The string Written into dest has the form "fac/dev/loc/prop".
//.	dest has to be provided by the caller.
//!EqAdr& operator = (EqAdr *a)
//.	The copy operator sets the left argument equal to the
//.	facility, device, location and property name of a.



//!Old calls (for compatibility):
//.void set ( char* facility, char* location, char* device, char* property);
//.void set_facility (char* facility);
//.void set_device   (char* device);
//.void set_location (char* location);
//.void set_property (char* property);
//.	To set all or the 4 parts individually of the address.
//.
//:EXAMPLE
//.EqAdr	ea, ea2;
//.
//.ea.adr ("TTF.VAC/MASS_SPECTR/VERT/RANGE");
//.
//.result = eq->get (&ea, &ed);	// reads "range"
//.
//.ea.merge ("///P");
//.
//.result = eq->get (&ea, &ed);	// reads "P" on the same device
//.
//.ae2 = ea;			// copy ea
//.
//.ea.merge ("///.HIST");
//.
//.result = eq->get (&ea, &ed);	// reads "P.HIST" on the same device
//.
//.ea2.merge ("///.EGU");	
//.result = eq->get (&ea2, &ed); // reads "P.EGU" on the same device
//.

//:BUGS
//.Limited wildcards. Only one "*" character is allowed in a device name.
//.And, the "*" must be the first or last character in a location or property string.


//:AUTHOR
//.KayRehlich   DESY  -PVAK-

#endif
