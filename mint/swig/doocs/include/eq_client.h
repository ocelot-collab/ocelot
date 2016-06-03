//! \file eq_client.h
//
//! \mainpage
//! 
//! \section intro Introduction
//!
//! Eq client classes
//!
//! The client classes are the interface to the RPC based communication.
//! The data is passed by the EqData class and the address is specified
//! by the EqAdr class.
//!
//! \author	K.Rehlich, O.Hensler, G.Grygiel, A.Agababyan
//

/*
 *
 * $Date: 2012-02-09 13:13:13 $
 * $Source: /doocs/doocssvr1/cvsroot/source/clientlib/eq_client.h,v $
 * $Revision: 1.40 $
 * $State: Exp $
 *
 * $Log: eq_client.h,v $
 * Revision 1.40  2012-02-09 13:13:13  arthura
 * *** empty log message ***
 *
 * Revision 1.39  2011/12/06 08:31:53  arthura
 * *** empty log message ***
 *
 * Revision 1.38  2011-01-31 15:21:54  arthura
 * *** empty log message ***
 *
 * Revision 1.37  2010/08/02 07:19:38  arthura
 * *** empty log message ***
 *
 * Revision 1.36  2010/05/10 09:01:56  arthura
 * *** empty log message ***
 *
 * Revision 1.35  2010/03/11 09:01:09  arthura
 * *** empty log message ***
 *
 * Revision 1.34  2009/11/17 09:20:18  arthura
 * new clientlib 5.8.0
 *
 * Revision 1.33  2009/08/10 06:55:43  arthura
 * *** empty log message ***
 *
 * Revision 1.32  2009/07/06 09:14:58  arthura
 * *** empty log message ***
 *
 * Revision 1.31  2008/08/19 09:12:50  arthura
 * *** empty log message ***
 *
 * Revision 1.30  2008/07/21 08:53:04  arthura
 * *** empty log message ***
 *
 * Revision 1.29  2008/05/16 10:41:35  arthura
 * *** empty log message ***
 *
 * Revision 1.28  2008/02/08 14:22:11  arthura
 * *** empty log message ***
 *
 * Revision 1.27  2005/04/29 07:18:53  arthura
 * *** empty log message ***
 *
 * Revision 1.26  2005/04/05 09:57:33  arthura
 * new clientlib 5.1.38
 *
 * Revision 1.25  2005/04/04 13:58:06  arthura
 * *** empty log message ***
 *
 * Revision 1.24  2005/03/18 16:16:18  arthura
 * *** empty log message ***
 *
 * Revision 1.23  2005/03/17 16:50:20  arthura
 * *** empty log message ***
 *
 * Revision 1.22  2005/01/19 15:27:11  arthura
 * *** empty log message ***
 *
 * Revision 1.21  2004/11/26 14:21:07  arthura
 * new clientlib 5.1.31
 *
 * Revision 1.20  2004/03/19 12:56:58  arthura
 * major change in TINE subsystem for monitoring.
 *
 * Revision 1.19  2003/07/31 09:38:53  arthura
 * avoid stop at doxygen
 *
 * Revision 1.18  2003/07/29 12:42:01  arthura
 * added 3 update rates for asynchronous mode
 *
 * Revision 1.17  2003/07/29 11:27:19  grygiel
 * Adaptations for the GNU C++ compiler
 *
 * Revision 1.16  2003/07/22 08:22:38  arthura
 * added disconnect method
 *
 * Revision 1.15  2003/07/21 07:05:50  ohensler
 * new clientlib 5.1.0
 *
 * Revision 1.14  2003/06/27 12:48:45  arthura
 * fixed the crash in rpc call
 *
 * Revision 1.13  2003/06/10 08:29:59  arthura
 * first MT clientlib
 *
 * Revision 1.12  2000/08/03 16:10:30  rehlich
 * bug fix:zero length hist SEGV; workaround: do not delete TINE monitors
 *
 * Revision 1.11  2000/06/28 07:46:13  rehlich
 * TINE bug fixes
 *
 * Revision 1.10  2000/06/19 15:10:12  ohensler
 * Added TINE async calls
 *
 * Revision 1.9  1999/09/21 16:26:59  rehlich
 * 1. TINE implementation (monitors emulated)
 *
 * Revision 1.8  1999/04/30 16:59:35  rehlich
 * includes shared memory
 *
 * Revision 1.6  1997/06/02 11:51:40  rehlich
 * next_server fixed, disconnect ENS after 120sec, disconnect Monitors with disconnect_all
 *
 * Revision 1.5  1997/02/10 14:50:00  rehlich
 * bug in names (MASK_LOCATION on)
 *
 * Revision 1.4  1997/01/23 13:31:14  rehlich
 * get/set_option new, EPICS monitor from ENS, new: * in location in get calls
 *
 * Revision 1.3  1996/11/07 10:16:28  rehlich
 * byte_array bug fix, static constructor for SunOS fixed
 *
 * Revision 1.2  1996/11/04 10:23:07  rehlich
 * pointer to data block restore, ENS disconnect, pointer to ServerEntry --> fixed
 *
 * Revision 1.1.1.1  1995/11/03 12:23:39  grygiel
 * DOOCS sources
 *
 * Revision 1.1  1995/10/24 13:04:18  xgrygiel
 * Initial revision
 *
 *
 */

#ifndef eq_client_h
#define eq_client_h

#include "eq_stl.h"

#include "mtp.h"
#include "eq_rpc.h"
#include "eq_adr.h"
#include "eq_data.h"
#include "eq_ens.h"
#include "eq_svr.h"
#include "eq_res.h"
#include "eq_shm.h"
#include "eq_sapi.h"

/*
         mask bits for the get_monitor () call:

         UPDATE_NORMAL       1
         UPDATE_FAST         2
         UPDATE_SLOW         4
         UPDATE_ONCE         8
         UPDATE_BY_SYSMASK   16
         
         UPDATE_ON_CHANGE
         UPDATE_CONNECTED

         scalar EqOption values for set_option () / get_option () calls:

         EQ_HOSTNAME
         EQ_PROTOCOL
         EQ_PROTOCOLID
         EQ_RETRY
         EQ_ONLINE
         EQ_TCP
         EQ_UDP
         EQ_AUTHTABLE
         EQ_AUTHMASK
         EQ_SVRNAME
         EQ_SVRMASK
         EQ_SVRSTATUS
         EQ_SVROPTION
         EQ_FILE_CHANNEL
         EQ_TIMEOUT
         EQ_UPDATEMASK
         EQ_SHM
         EQ_SHM_CREATE
         EQ_SHM_DELETE
         EQ_UPDATEMODE
         EQ_LIFETIME
         EQ_CALLTIME
         EQ_YIELD
         EQ_SBUFSIZE
         EQ_RBUFSIZE
         EQ_LIBNO
         EQ_DEBUG
         EQ_TINEVERS
};
*/



class	EqCall {

private:
        static void    *pdap;

        eq_ens         *ens_call;
        eq_shm         *shm_call;
        eq_res         *res_call;

        SEM_T          sem;
        EqData         data;

public:
                   EqCall         (void);
                   ~EqCall        (void);

        // new MT-safe API

        int        names          (EqAdr *, EqData *);
        int        get            (EqAdr *, EqData *, EqData *);
        int        set            (EqAdr *, EqData *, EqData *);
        int        get_option     (EqAdr *, EqData *, EqData *, EqOption);
        int        set_option     (EqAdr *, EqData *, EqData *, EqOption);
        int        get_monitor    (EqAdr *, EqData *, EqData *, int, float = 0.0/*, EQ_CB = EQ_IGNORE, void * = (void *) 0*/);
        int        clear_monitor  (EqAdr *);
        int        update         (int);
        int        servers        (void);
        int        monitors       (void);
        int        disconnect     (void);

        // old MT-unsafe API

        EqData	   *names         (EqAdr *);
        EqData	   *get           (EqAdr *, EqData *);
        EqData	   *set           (EqAdr *, EqData *);
        EqData	   *get_option    (EqAdr *, EqOption);
        EqData	   *set_option    (EqAdr *, EqOption, EqData * = 0);
        EqData	   *get_monitor   (EqAdr *, EqData *, int, float = 0.0);
};




/*
 Class library to do the actual transfer of the data from/to the devices.
 A transfer always needs an address defined by an object of the EqAdr 
 class and returns a data object of EqData class. Some calls also need
 a data object to be send to the device. The current implementation has
 3 communication methods - a get, a set and a names calls. The get call
 is used read device data, the set call sets device data and the names
 call returns a list with the available names in the communication network.
 All calls return a pointer to an EqData structure. This structure is
 freed on the next call of get/set/names. To keep the data it is necessary
 to copy the result or to create a further EqCall instance. "disconnect" 
 disconnects all servers from this application. Three additional calls are
 used to operate with monitors. A monitor buffers the data from one property
 (channel) in a local EqData structure. The data is updated by an update call.
 A first call of get_monitor reads the data from the device and creates 
 a monitor. Further calls get the data from the buffer without doing 
 a network call. This interface is prepared to use EPICS monitors also.
 The user program needs not to be changed when this feature is added later.

   		EqCall ()

	Constructor

   int         get   (EqAdr *, EqData *, EqData *)
   EqData      *get  (EqAdr *, EqData *)

	The get call sends a data and an address object to the device 
	and reads back a data object. The data send to the device is
        only usefull for some properties to specify a time range or
        to select a stored spectrum for instance.


   int         set   (EqAdr *, EqData *, EqData *)
   EqData      *set  (EqAdr *, EqData *)

	The set call sends a data and an address object to the device
        and gets a result data object returned. The data returned data
        contains the error messages of the call.

   int         names  (EqAdr *, EqData *)
   EqData      *names (EqAdr *)

	The names call is used to query names from the network. A wild card
        character "*" in one of the 4 components of the address in the names
        call selects the list which will be returned. The data type of the
        returned data is always DATA_A_USTR. Currently the "*" is the only
        supported wildcard character. A "*" in the facility, device, location
        or property field gives a list of the according facilities, devices,
        locations or properties. The list of names are in the ASCII part of
        the DATA_A_USTR structure. The list of properties contains the data
        type in the integer part of the DATA_A_USTR structure also. A description
        of the property is appended to the property name string and separated
        by a blank character.

   int         get_monitor  (EqAdr *, EqData *, EqData *, int, float = 0)
   EqData      *get_monitor (EqAdr *, EqData *, int, float = 0)

	Reads back the data from a buffer (async call). It creates a 
        monitor instance. Further calls to the same address read from the
        monitor instance only. All monitors are updated by the by the update
        call. The mask arguments specifies an update mask which must match
        one or more mask bits in the update call. The ref_rate arg is used
        to specify the max refresh reate of the device server in seconds
        (for TINE calls only, float parameter).

   int         update (int)

	Updates all monitors with the corresponding mask bit set.
	Convention for masks are: 1 = normal, 2 = fast, 4 = slow updates.

   void        clear_monitor (EqAdr *)

	Clears the monitor from the list of updated entries.

   int         get_option  (EqAdr *, EqData *, EqData *, enum EqOption)
   EqData      *get_option (EqAdr *, enum EqOption)

   int         set_option  (EqAdr *, EqData *, EqData *, enum EqOption)
   EqData      *set_option (EqAdr *, enum EqOption, EqData * = 0)

	Read and set optional parameters of the communication. It also
	shows most of the internal settings known to the API library.

	EqOption   read/write	Function

	EQ_HOSTNAME	    R	hostname of the server
	EQ_PROTOCOL	    R	protocol name of communication
	EQ_PROTOCOLID	R	protocol as integer id
	EQ_RETRY	    W	resets the retry counter and reconnects
	EQ_ONLINE 	    R	device is online flag from ENS
	EQ_TCP		    RW	use TCP protocol (DOOCS only) *)
	EQ_UDP		    RW	use UDP protocol (DOOCS only) *)
	EQ_AUTHTABLE	R	shows who has access permissions
	EQ_AUTHMASK	    R	shows authentication for this client
	EQ_SVRNAME	    R	name of the server
	EQ_SVRMASK	    R	server mask: part of name which points to the server *)
	EQ_SVRSTATUS	R	connection status (enum ConState)
	EQ_SVROPTION	R	available server options
	EQ_FILE_CHANNEL	R	file or channel name argument from ENS
	EQ_TIMEOUT	    RW	timeout parameter in calls *)
	EQ_UPDATEMASK	RW	priority of updates 
                        (UPDATE_NORMAL, UPDATE_FAST, UPDATE_SLOW)
	EQ_SHM		    RW	shows that the shm creator-process exists (R).
				        On (W) deletes its name when creator does not exist
	EQ_SHM_CREATE	W	creates the shared memory of the type specified in
				        EqData and the name specified in EqAdr
	EQ_SHM_DELETE	W	deletes the shared memory specified in EqAdr
	EQ_UPDATEMODE	RW	synchronous/asynchronous mode for get_monitor ()
                        (IIII type in EqData):
                         i1_data: update interval timer (sec)  (for FAST);
                         i2_data: update interval timer (usec) (for FAST);
                         i3_data: decrement count (for NORMAL);
                         i4_data: decrement count (for SLOW);
	EQ_LIFETIME	    RW	defines the duration of the existense of the
                        monitor created by get_monitor ();
	EQ_YIELD	    RW	enables/disables 100 ms pause between monitor updates


	*) not yet implemented


	EXAMPLE 1 reads a single item:

	#include	<eq_client.h>

	EqCall      eq;
	EqAdr		ea;
	EqData		src;
	EqData		dst;

    int         err;
    char        buf [256];

	ea.adr ("TTF.VAC/MASS_SPECTR/VERT/RANGE");

	err = eq.get (&ea, &src, &dst);

	The RANGE information from a mass spectrometer at position "VERT" is read and
	the resulting data is in "result". The error code of result is returned either 
    in  err  or by:

	int err_code = dst.error ();

	The data is:

	float f = dst.get_float ();


	EXAMPLE 2 writes a single item:

	ea.adr ("TTF.VAC/MASS_SPECTR/VERT/RANGE");

    src.set (3);
	err = eq.set (&ea, &src, &dst);

	int err_code = dst.error ();

*/

#endif
