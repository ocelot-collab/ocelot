/* file : eq_errors.h
//
//	Client ERROR codes
//
// text strings in : eq_errors.cc
//
// Kay Rehlich
// last update:
//  7. Oct. 1993
*/


/*
 *
 * $Date: 2012-02-10 12:03:52 $
 * $Source: /doocs/doocssvr1/cvsroot/source/clientlib/eq_errors.h,v $
 * $Revision: 1.16 $
 * $State: Exp $
 *
 * $Log: eq_errors.h,v $
 * Revision 1.16  2012-02-10 12:03:52  arthura
 * *** empty log message ***
 *
 * Revision 1.15  2012/02/09 13:13:15  arthura
 * *** empty log message ***
 *
 * Revision 1.14  2011/12/06 08:31:54  arthura
 * *** empty log message ***
 *
 * Revision 1.13  2009-08-10 06:55:45  arthura
 * *** empty log message ***
 *
 * Revision 1.12  2009/05/04 13:41:46  arthura
 * *** empty log message ***
 *
 * Revision 1.11  2008/02/08 14:22:12  arthura
 * *** empty log message ***
 *
 * Revision 1.10  2005/07/08 14:33:43  arthura
 * *** empty log message ***
 *
 * Revision 1.9  2005/04/11 14:49:28  arthura
 * *** empty log message ***
 *
 * Revision 1.8  2004/03/03 12:42:57  arthura
 * fixed the problem of timeouts and erroneous behaviour of names () call
 * when reading ENS database.
 *
 * Revision 1.7  2003/06/27 12:48:46  arthura
 * fixed the crash in rpc call
 *
 * Revision 1.5  2003/06/10 08:30:00  arthura
 * first MT clientlib
 *
 * Revision 1.4  2000/07/19 11:21:41  rehlich
 * bug fixes for TINE
 *
 * Revision 1.3  1999/09/21 16:27:06  rehlich
 * 1. TINE implementation (monitors emulated)
 *
 * Revision 1.2  1997/01/23 13:31:27  rehlich
 * get/set_option new, EPICS monitor from ENS, new: * in location in get calls
 *
 * Revision 1.1.1.1  1995/11/03 12:23:39  grygiel
 * DOOCS sources
 *
 * Revision 1.1  1995/10/24 13:04:18  xgrygiel
 * Initial revision
 *
 *
 */


#ifndef eq_errors_h
#define eq_errors_h

#define ERR_NONE             0
#define ERR_RPC	             100 
#define ERR_ILL_FACILITY     101 /* old name */
#define ERR_ILL_SERVICE      101 /* new name */
#define ERR_RPC_TO           102
#define ERR_ENS              103
#define ERR_ILL_MON          104
#define ERR_ILL_PROTOCOL     105
#define ERR_EPICS_CALL       106
#define ERR_EPICS_UNSUP      107
#define ERR_NO_EPICS         108
#define ERR_EPICS_UNSUP_DAT  109
#define ERR_UNSUP_DAT        110
#define ERR_DIFF_DAT         111
#define ERR_OPT	             112
#define ERR_RO_OPT           113
#define ERR_ILL_TYPE         114
#define ERR_STALE_DATA       115
#define ERR_OFFLINE          116
#define ERR_TMO	             117
#define ERR_NO_DATA          118
#define ERR_ENS_NO_DATA      119
#define ERR_FAULTY_CHANS     120
#define ERR_SHMEM            121
#define ERR_UNSUP_SERV 	     122

extern char     *str_no_error;
extern char     *str_unav_service;
extern char     *str_ill_service;
extern char     *str_timeout;
extern char     *str_enserror;
extern char     *str_monerror;
extern char     *str_ill_protocol;
extern char     *str_epics_call;
extern char     *str_epics_unsuprec;
extern char     *str_no_epics;
extern char     *str_epics_unsupdat;
extern char     *str_unsupdat;
extern char     *str_diffdat;
extern char     *str_ill_option;
extern char     *str_ro_option;
extern char     *str_ill_type;
extern char     *str_stale_data;
extern char     *str_offline;
extern char     *str_tmo;
extern char     *str_no_data;
extern char     *str_ens_no_data;
extern char     *str_faulty_chans;
extern char     *str_shmem;
extern char     *str_unsupserv;

#endif
