/* file eq_data.h
//
// EqData class
//
// Kay Rehlich, -PVAK-
//
// last update:
// 19. Mar. 1993
//  8. Aug. 1995	new data types
*/


/*
 *
 * $Date: 2012-10-15 11:59:43 $
 * $Source: /doocs/doocssvr1/cvsroot/source/clientlib/eq_data.h,v $
 * $Revision: 1.59 $
 * $State: Exp $
 *
 * $Log: eq_data.h,v $
 * Revision 1.59  2012-10-15 11:59:43  arthura
 * *** empty log message ***
 *
 * Revision 1.58  2012/06/11 07:55:38  arthura
 * *** empty log message ***
 *
 * Revision 1.57  2012/05/23 13:43:21  arthura
 * *** empty log message ***
 *
 * Revision 1.56  2012/02/09 13:13:14  arthura
 * *** empty log message ***
 *
 * Revision 1.55  2011/12/06 08:31:53  arthura
 * *** empty log message ***
 *
 * Revision 1.54  2011-08-15 11:31:13  arthura
 * introduced a mapping DOOCS-TINE for KEYVAL, removed some other types
 *
 * Revision 1.53  2011/07/28 08:05:45  arthura
 * added new data types: DATA_A_BOOL, DATA_KEYVAL
 *
 * Revision 1.52  2010/06/14 13:33:03  arthura
 * *** empty log message ***
 *
 * Revision 1.51  2010/03/11 09:01:11  arthura
 * *** empty log message ***
 *
 * Revision 1.50  2010/02/24 15:39:24  arthura
 * *** empty log message ***
 *
 * Revision 1.49  2010/01/06 08:15:27  arthura
 * new clientlib 5.8.1
 *
 * Revision 1.48  2009/11/18 14:17:23  arthura
 * new clientlib 5.8.0
 *
 * Revision 1.47  2009/11/10 12:06:52  arthura
 * *** empty log message ***
 *
 * Revision 1.46  2009/08/10 06:55:44  arthura
 * *** empty log message ***
 *
 * Revision 1.45  2009/06/22 13:55:30  arthura
 * *** empty log message ***
 *
 * Revision 1.44  2009/06/22 07:52:46  arthura
 * *** empty log message ***
 *
 * Revision 1.43  2008/06/30 09:58:32  arthura
 * *** empty log message ***
 *
 * Revision 1.42  2008/04/17 07:23:20  arthura
 * *** empty log message ***
 *
 * Revision 1.41  2007/11/14 11:55:13  arthura
 * *** empty log message ***
 *
 * Revision 1.40  2007/10/05 08:21:51  arthura
 * *** empty log message ***
 *
 * Revision 1.39  2007/08/24 09:58:43  arthura
 * added image data type
 *
 * Revision 1.38  2007/06/06 10:26:50  arthura
 * *** empty log message ***
 *
 * Revision 1.37  2007/04/24 10:32:38  arthura
 * new clientlib 5.5.0
 *
 * Revision 1.35  2006/09/14 12:10:10  arthura
 * introduces a "cycling".
 *
 * Revision 1.34  2006/03/30 13:06:51  arthura
 * *** empty log message ***
 *
 * Revision 1.33  2006/01/30 08:47:03  arthura
 * added a new data type DATA_MDA_FLOAT - multidimentional float array
 *
 * Revision 1.32  2005/10/20 15:40:53  arthura
 * *** empty log message ***
 *
 * Revision 1.31  2005/10/20 15:39:24  arthura
 * *** empty log message ***
 *
 * Revision 1.30  2005/10/20 15:36:19  arthura
 * *** empty log message ***
 *
 * Revision 1.29  2005/10/20 15:30:14  arthura
 * *** empty log message ***
 *
 * Revision 1.28  2005/10/20 15:11:03  arthura
 * *** empty log message ***
 *
 * Revision 1.27  2005/10/20 14:16:00  arthura
 * linux bug fix for EPICS monitors, complete description of EqData class
 *
 * Revision 1.26  2005/08/18 08:15:41  arthura
 * removed void get_xml (char *, int) member function
 *
 * Revision 1.25  2005/04/21 07:16:39  arthura
 * increased the size of the DATA_A_BYTE object to 8MB
 * added two methods to eq_data class to process DATA_A_BYTE
 * type: set_byte (len, ptr) and get_byte (*len, **ptr)
 *
 * Revision 1.24  2005/03/14 10:49:30  arthura
 * *** empty log message ***
 *
 * Revision 1.23  2005/03/01 13:46:50  arthura
 * Some changes in EqData class to fix bugs and improve the code
 * major changes in the internal structure of eq_res and eq_svr
 * classes to allow multisource addresses for DOOCS channels,
 * complete redesign of EPICS interface, redesign of TINE interface
 * to deal in event handler instead if the address with the ID.
 * Still there is a bug in the TINE subsystem releated to the
 * memory leak in TINE stack after calling destroy () method.
 * This bug fix comes later.
 *
 * Revision 1.22  2004/11/26 14:21:07  arthura
 * new clientlib 5.1.31
 *
 * Revision 1.21  2004/11/16 09:10:37  arthura
 * new clientlib 5.1.30
 *
 * Revision 1.20  2004/10/12 11:14:28  arthura
 * removed MT unsafe functions:
 * time_string (), get_time_string (), get_string (), get_string_arg (), get_string_extra ()
 * CVe: ----------------------------------------------------------------------
 *
 * Revision 1.19  2004/10/06 14:28:32  arthura
 * added the following methods to replace the old ones which are not MT safe:
 * char *time_string (void)
 * char *get_time_string (void)
 * char *get_time_string (int)
 * char *get_string (void)
 * char *get_string (int)
 * char *get_string_arg (void)
 * char *get_string_arg (int)
 * char *get_string_extra (int)
 *
 * modified the following files where these functions are used:
 * DOOCSapi.cc, eq_shm.cc, eq_res.cc, eq_svr.cc, lv_client.cc
 *
 * Revision 1.18  2004/09/29 12:42:49  arthura
 * *** empty log message ***
 *
 * Revision 1.17  2003/11/14 10:18:28  arthura
 * new clientlib 5.1.10
 *
 * Revision 1.16  2003/10/16 07:27:59  arthura
 * added a new data type: DATA_XYZS and DATA_A_XYZS
 * changed TIMEOUT value for RPC protocol from 25 sec to 8 sec
 * changed retry values of rpc threads from 10 to 5
 * added file pvak_types_xdr.cc
 * moved config class (config.h, config.cc) from clientlib to serverlib
 *
 * Revision 1.15  2003/06/27 12:48:45  arthura
 * fixed the crash in rpc call
 *
 * Revision 1.14  2003/06/10 09:19:56  ohensler
 * first clientlib MT without AList and MkvString
 *
 * Revision 1.13  2003/06/10 08:29:59  arthura
 * first MT clientlib
 *
 * Revision 1.12  2003/01/27 08:36:17  arthura
 * fixed the bug in EqData class
 * modified get_tds2 (int) and get_tds (int) functions in eqData class
 *
 * Revision 1.11  2001/10/23 12:02:36  ohensler
 * new clientlib 4.3.3
 *
 * Revision 1.10  2000/03/02 17:53:46  rehlich
 * bug fix: SEGV in free -> string termination of XDR string removed
 *
 * Revision 1.9  1999/07/14 12:48:07  rehlich
 * move math functions out of eq_data, casts for signals modified
 *
 * Revision 1.8  1999/04/30 16:59:40  rehlich
 * includes shared memory
 *
 * Revision 1.7  1998/09/11 15:18:15  rehlich
 * bug fix in spectrum comment with zero length
 *
 * Revision 1.6  1998/02/12 09:44:53  rehlich
 * longer EPICS timeout, longer delay if first connect failed
 *
 * Revision 1.5  1997/05/02 07:00:27  rehlich
 * environment vari NOEPICS, update of monitor in get call, bug fix in get_locations(eq_services)
 *
 * Revision 1.4  1997/01/23 13:31:22  rehlich
 * get/set_option new, EPICS monitor from ENS, new: * in location in get calls
 *
 * Revision 1.3  1996/11/04 10:23:12  rehlich
 * pointer to data block restore, ENS disconnect, pointer to ServerEntry --> fixed
 *
 * Revision 1.2  1996/05/15 08:06:57  rehlich
 * DATA_A_BYTE fromat changed, minor bugs corrected
 *
 * Revision 1.1.1.1  1995/11/03 12:23:39  grygiel
 * DOOCS sources
 *
 * Revision 1.1  1995/10/24 13:04:18  xgrygiel
 * Initial revision
 *
 *
 */

#ifndef eq_data_h
#define eq_data_h

#include  "eq_types.h"
#include  "math.h"

#define	AUTH_MARKER    0x80000000
#define	CYCLE_MARKER   0x40000000


/**
 * EqData is a container class which simplifies an access to specific
 * fields of data being sent to or received from servers througth its
 * member functions. At present this class supports the following
 * data types:
 * 
 *       \tScalar types:\t        Array types:\n
 * 
 *       \tDATA_NULL\t            DATA_A_INT\n
 *       \tDATA_BOOL\t            DATA_A_FLOAT\n
 *       \tDATA_INT\t             DATA_A_XY\n
 *       \tDATA_FLOAT\t           DATA_A_TDS\n
 *       \tDATA_STRING\t          DATA_A_XYZS\n
 *       \tDATA_STRING16\t        DATA_A_BYTE\n
 *       \tDATA_XY\t              DATA_A_USTR\n
 *       \tDATA_TDS\t
 *       \tDATA_IIII\t
 *       \tDATA_IFFF\n
 *       \tDATA_TTII\n
 *       \tDATA_SPECTRUM\n
 *       \tDATA_XML\n
 *       \tDATA_XYZS\n
 * 
 *  The exact structure of each mentioned data type can be found in
 *  the file pvak_types.h
 *
 *  All data provided to application programs by the DOOCS API are incapsulated
 *  in EqData class and the other way round, data provided by applications
 *  to the API to be trasferred to servers should be saved in EqData objects.
 *  EqData class is an interface between application and the API which hides
 *  all specifics of real data being transferred over the netwoks.
 *
 *  All member funtions if they have as a return value of type int
 *  then they return 1 in case of success or 0 in case of fail.
 *
 **/

class EqData {

private:
        size_t	      space;      // size of allocated data
        size_t	      preserve;   // preserve allocated data
        char         *buffer;    // pointer to allocated data for dbp
        EqDataBlock  *dbpo;      // pointer to allocated data block
        EqDataBlock  *dbp;       // pointer to actual data block

        void         assign_type       (int type);
        void         alloc_space       (size_t len);
        void         extend_space      (size_t len);
        void         alloc_array       (void);
        void         free_array        (void);
        int          check_error       (char *bufp, size_t len);
        char         *str_copy         (char *bp, size_t blen, const char *sp, size_t slen);
        float        *get_pos          (int *dp, int dims, MDA_FLOAT *mp);
        int          find_key          (char *keyp);

public:
	/**
	 *  Constructor. Internaly allocates EqDataBlock
	 *
	**/
                     EqData            (void);

	/**
	 *  Destructor. Frees all internally allocated memory
	 *
	**/
                     ~EqData           (void);

	/**
	 *  Constructor. EqDataBlock supplied externaly
	 *
	 * @param edbp  The external EqDataBlock 
	 * @param flag  1, prohibits class destructor to free EqDataBlock,\n
         *              0 allows to free all resources allocated
	**/
                     EqData            (EqDataBlock *edbp, int flag = 0);

	/**
	 *  Internal call to supply an external EqDataBlock
	 *
	 * @param edbp  The external EqDataBlock 
	**/
	void	        data_block     (EqDataBlock *edbp);

	/**
	 *  Internal call to return an EqDataBlock
	 *
	 * @returns   The pointer to internal EqDataBlock 
	**/
	EqDataBlock	*data_block       (void);


	/**
	 *  Initializes the object and resets data type to DATA_NULL
	 *
	**/
	void	        init              (void);

	/**
	 *  Copies content of EqData object, acts only for scalar types
	 *
	 * @param ed  Source EqData object
	**/
	void	        load_simple       (EqData *ed);

	/**
	 *  Copies content of EqData object including arrays
	 *
	 * @param ed  Source EqData object
	**/
	void		copy_from         (EqData *ed);

	/**
	 *  Appends content of EqData object including arrays
	 *
	 * @param ed  Source EqData object
	**/
	void		add_from          (EqData *ed);


	/**
	 *  Returns the time stamp of the data.
	 *
	 * @return       The time stamp 
	**/
	time_t		time              (void);

	/**
	 *  Returns the error code of the data.
	 *
	 * @return       The error 
	**/

	int		error             (void);

	/**
	 *  Returns the authentication mask of the data.
	 *
	 * @return       The authentication mask 
	**/
	int		auth_mask         (void);

	/**
	 *  Returns the cycle mask of the data.
	 *
	 * @return       The cycle mask 
	**/
	int		cycle_mask        (void);

	/**
	 *  Returns the data type.
	 *
	 * @return       The data type 
	**/

	int		type              (void);

	/**
	 *  Returns the data type as an ascii string
	 *
	 * @return       The data type string
	**/

	char		*type_string      (void);

	/**
	 *  Returns an ascii string for provided data type
	 *
	 * @param dtype  The provided data type
	 * @return       The data type string
	**/
	char		*type_string      (int dtype);

	/**
	 *  Returns the limit of the specified type.
	 *
	 * @param type   Any data type supported by DOOCS
	 * @return       The data type size 
	**/

	int		get_limits        (int type);

	/**
	 *  Returns the data length
         *
         *  Applicable to:\n
         *  DATA_IIII           - returns 4\n
         *  DATA_IFFF           - returns 4\n
         *  DATA_STRING         - returns length of string\n
         *  DATA_STRING16       - returns length of short string\n
         *  DATA_XML            - returns length of XML array\n
         *  DATA_SPECTRUM       - returns number of elements in spectrum array\n
         *  DATA_GSPECTRUM      - returns number of elements in grouped spectrum array\n
         *  DATA_IMAGE          - returns number of elements in image array\n
         *  DATA_A_SHORT        - returns number of elements in short array\n
         *  DATA_A_LONG         - returns number of elements in long long array\n
         *  DATA_A_BOOL         - returns number of elements in bool array\n
         *  DATA_A_INT          - returns number of elements in int array\n
         *  DATA_A_FLOAT        - returns number of elements in float array\n
         *  DATA_A_TDS          - returns number of elements in TDS array\n
         *  DATA_A_XY           - returns number of elements in XY array\n
         *  DATA_A_USTR         - returns number of elements in USTR array\n
         *  DATA_A_XYZS         - returns number of elements in XYZS array\n
         *  DATA_A_THUMBNAIL    - returns number of elements in THUMBNAIL array\n
         *  DATA_A_BYTE         - returns number of elements in BYTE array\n
         *  DATA_A_TS_BOOL      - returns number of elements in TS_BOOL array\n
         *  DATA_A_TS_INT       - returns number of elements in TS_INT array\n
         *  DATA_A_TS_FLOAT     - returns number of elements in TS_FLOAT array\n
         *  DATA_A_TS_DOUBLE    - returns number of elements in TS_DOUBLE array\n
         *  DATA_A_TS_STRING16  - returns number of elements in TS_STRING16 array\n
         *  DATA_A_TS_STRING    - returns number of elements in TS_STRING array\n
         *  DATA_A_TS_USTR      - returns number of elements in TS_USTR array\n
         *  DATA_A_TS_XML       - returns number of elements in TS_XML array\n
         *  DATA_A_TS_XY        - returns number of elements in TS_XY array\n
         *  DATA_A_TS_IIII      - returns number of elements in TS_IIII array\n
         *  DATA_A_TS_IFFF      - returns number of elements in TS_IFFF array\n
         *  DATA_A_TS_TTII      - returns number of elements in TS_TTII array\n
         *  DATA_A_TS_XYZS      - returns number of elements in TS_XYZS array\n
         *  DATA_A_TS_SPECTRUM  - returns number of elements in TS_SPECTRUM array\n
         *  DATA_A_TS_GSPECTRUM - returns number of elements in TS_GSPECTRUM array\n
         *  DATA_KEYVAL         - returns number of elements in KEYVAL table\n
         *  others              - returns 1\n
	 *
	 * @return       The data length
	**/
	int			length            (void);

	/**
	 *  Returns the length of array type data
         *
         *  Applicable to:\n
         *  DATA_XML            - returns length of XML array\n
         *  DATA_SPECTRUM       - returns number of elements in spectrum array\n
         *  DATA_GSPECTRUM      - returns number of elements in grouped spectrum array\n
         *  DATA_IMAGE          - returns number of elements in image array\n
         *  DATA_A_SHORT        - returns number of elements in short array\n
         *  DATA_A_LONG         - returns number of elements in long long array\n
         *  DATA_A_BOOL         - returns number of elements in bool array\n
         *  DATA_A_INT          - returns number of elements in int array\n
         *  DATA_A_FLOAT        - returns number of elements in float array\n
         *  DATA_A_TDS          - returns number of elements in TDS array\n
         *  DATA_A_XY           - returns number of elements in XY array\n
         *  DATA_A_USTR         - returns number of elements in USTR array\n
         *  DATA_A_XYZS         - returns number of elements in XYZS array\n
         *  DATA_A_THUMBNAIL    - returns number of elements in THUMBNAIL array\n
         *  DATA_A_BYTE         - returns number of elements in BYTE array\n
         *  DATA_A_TS_BOOL      - returns number of elements in TS_BOOL array\n
         *  DATA_A_TS_INT       - returns number of elements in TS_INT array\n
         *  DATA_A_TS_FLOAT     - returns number of elements in TS_FLOAT array\n
         *  DATA_A_TS_DOUBLE    - returns number of elements in TS_DOUBLE array\n
         *  DATA_A_TS_STRING16  - returns number of elements in TS_STRING16 array\n
         *  DATA_A_TS_STRING    - returns number of elements in TS_STRING array\n
         *  DATA_A_TS_USTR      - returns number of elements in TS_USTR array\n
         *  DATA_A_TS_XML       - returns number of elements in TS_XML array\n
         *  DATA_A_TS_XY        - returns number of elements in TS_XY array\n
         *  DATA_A_TS_IIII      - returns number of elements in TS_IIII array\n
         *  DATA_A_TS_IFFF      - returns number of elements in TS_IFFF array\n
         *  DATA_A_TS_TTII      - returns number of elements in TS_TTII array\n
         *  DATA_A_TS_XYZS      - returns number of elements in TS_XYZS array\n
         *  DATA_A_TS_SPECTRUM  - returns number of elements in TS_SPECTRUM array\n
         *  DATA_A_TS_GSPECTRUM - returns number of elements in TS_GSPECTRUM array\n
         *  DATA_KEYVAL         - returns number of elements in KEYVAL table\n
         *  others        - returns 0\n
	 *
	 * @return       The length of array type data
	**/
	int			array_length      (void);


	/**
	 *  Returns the length of string type data with index\n
         *
         *  Applicable to:\n
         *  DATA_STRING    - returns length of string\n
         *  DATA_STRING16  - returns length of short string\n
         *  DATA_SPECTRUM  - returns length of comment string of spectrum\n
         *  DATA_GSPECTRUM - returns length of comment string of grouped spectrum\n
         *  DATA_IMAGE     - returns length of comment string of image\n
         *  DATA_A_USTR    - returns length of string in USTR element\n
         *  DATA_A_XYZS    - returns length of string in XYZS element\n
         *  DATA_XML       - returns number of elements in XML array\n
	 *
	 * @param index  The index of data array element
	 * @return       The length of string type data
	**/
	int			string_length     (int index = 0);

	/**
	 *  Returns int array
         *
         *  Applicable to:\n
         *  DATA_STRING    - returns pointer to string\n
         *  DATA_TEXT      - returns pointer to text\n
         *  DATA_XML       - returns pointer to XML string\n
         *  DATA_A_BYTE    - returns pointer to byte array\n
         *  DATA_IMAGE     - returns pointer to image array\n
         *  others         - returns 0\n
	 *
	 * @return       The array pointer
	**/
	char			*get_char_array  (void);

	/**
	 *  Returns long long with index
         *
         *  Applicable to:\n
         *  DATA_A_LONG    - returns element of long long array\n
         *
	 * @param index  The index of data array element
	 * @return       The integer
	**/
	long long              get_long         (int index);

	/**
	 *  Returns long long array
         *
         *  Applicable to:\n
         *  DATA_A_LONG    - returns pointer to long long array\n
         *  others         - returns 0\n
	 *
	 * @return       The array pointer
	**/
	long long	       *get_long_array  (void);

	/**
	 *  Returns integer
         *
         *  Applicable to:\n
         *  DATA_BOOL   - returns d_bool\n
         *  DATA_INT    - returns d_int\n
         *  DATA_FLOAT  - returns d_float\n
         *  DATA_DOUBLE - returns d_double\n
	 *
	 * @return       The integer
	**/
	int			get_int           (void);

	/**
	 *  Returns integer with index
         *
         *  Applicable to:\n
         *  DATA_BOOL      - returns d_bool\n
         *  DATA_INT       - returns d_float\n
         *  DATA_FLOAT     - returns d_float\n
         *  DATA_DOUBLE    - returns double\n
         *  DATA_A_INT     - returns element of int array\n
         *  DATA_A_BOOL    - returns element of bool array\n
         *  DATA_A_SHORT   - returns element of short array\n
         *  DATA_A_FLOAT   - returns element of float array\n
         *  DATA_A_DOUBLE  - returns element of double array\n
         *  DATA_A_TDS     - returns tm element of TDS array\n
         *  DATA_A_BYTE    - returns element of BYTE array\n
         *  DATA_A_USTR    - returns the following fields of USTR [0] element\n
         *                   \ti1_data if index == 0\n
         *                   \tf1_data if index == 1\n
         *                   \tf2_data if index == 2\n
         *                   \ttm      if index == 3\n
         *  DATA_A_XYZS    - returns the following fields of XYZS [0] element\n
         *                   \tstatus if index == 0\n
         *                   \tx      if index == 1\n
         *                   \ty      if index == 2\n
         *                   \tz      if index == 3\n
         *  DATA_IIII      - returns\n
         *                   \ti1_data if index == 0\n
         *                   \ti2_data if index == 1\n
         *                   \ti3_data if index == 2\n
         *                   \ti4_data if index == 3\n
         *  DATA_IFFF      - returns\n
         *                   \ti1_data if index == 0\n
         *                   \tf1_data if index == 1\n
         *                   \tf2_data if index == 2\n
         *                   \tf3_data if index == 3\n
         *  DATA_MDA_FLOAT - returns d_float\n
         *  DATA_IMAGE     - returns d_image\n
	 *
	 * @param index  The index of data array element
	 * @return       The integer
	**/
	int			get_int           (int index);

	/**
	 *  Returns integer array
         *
         *  Applicable to:\n
         *  DATA_A_INT     - returns pointer to int array\n
         *  DATA_A_BOOL    - returns pointer to bool array\n
         *  DATA_A_SHORT   - returns pointer to short array\n
         *  DATA_A_LONG    - returns pointer to long long array\n
         *  others         - returns 0\n
	 *
	 * @return       The array pointer
	**/
	int			*get_int_array  (void);

	/**
	 *  Returns float
         *
         *  Applicable to:\n
         *  DATA_INT    - returns d_float\n
         *  DATA_FLOAT  - returns d_float\n
         *  DATA_DOUBLE - returns d_double\n
	 *
	 * @return       The float
	**/
	float			get_float         (void);

	/**
	 *  Returns float with index
         *
         *  Applicable to:\n
         *  DATA_SPECTRUM  - returns element of SPECTRUM float array\n
         *  DATA_GSPECTRUM - returns element of GSPECTRUM float array\n
         *  DATA_A_INT     - returns element of int array\n
         *  DATA_A_FLOAT   - returns element of float array\n
         *  DATA_A_DOUBLE  - returns element of double array\n
         *  DATA_A_TDS     - returns data element of TDS array\n
         *  DATA_A_BYTE    - returns element of BYTE array\n
         *  DATA_A_USTR    - returns the following fields of USTR [0] element\n
         *                   \ti1_data if index == 0\n
         *                   \tf1_data if index == 1\n
         *                   \tf2_data if index == 2\n
         *                   \ttm      if index == 3\n
         *  DATA_A_XYZS    - returns the following fields of XYZS [0] element\n
         *                   \tstatus if index == 0\n
         *                   \tx      if index == 1\n
         *                   \ty      if index == 2\n
         *                   \tz      if index == 3\n
         *  DATA_IIII      - returns\n
         *                   \ti1_data if index == 0\n
         *                   \ti2_data if index == 1\n
         *                   \ti3_data if index == 2\n
         *                   \ti4_data if index == 3\n
         *  DATA_IFFF      - returns\n
         *                   \ti1_data if index == 0\n
         *                   \tf1_data if index == 1\n
         *                   \tf2_data if index == 2\n
         *                   \tf3_data if index == 3\n
         *  DATA_MDA_FLOAT - returns d_float\n
	 *
	 * @param index  The index of data array element
	 * @return       The float
	**/
	float			get_float         (int index);

	/**
	 *  Returns float array
         *
         *  Applicable to:\n
         *  DATA_SPECTRUM  - returns pointer to SPECTRUM float array\n
         *  DATA_GSPECTRUM - returns pointer to GSPECTRUM float array\n
         *  DATA_A_FLOAT   - returns pointer to float array\n
         *  DATA_MDA_FLOAT - returns pointer to multidimentional float array\n
         *  others         - returns 0\n
	 *
	 * @return       The array pointer
	**/
	float			*get_float_array  (void);

	/**
	 *  Returns double
         *
         *  Applicable to:\n
         *  DATA_INT    - returns d_float\n
         *  DATA_FLOAT  - returns d_float\n
         *  DATA_DOUBLE - returns d_double\n
	 *
	 * @return       The double
	**/
	double			get_double         (void);

	/**
	 *  Returns double with index
         *
         *  Applicable to:\n
         *  DATA_SPECTRUM  - returns element of SPECTRUM float array\n
         *  DATA_GSPECTRUM - returns element of GSPECTRUM float array\n
         *  DATA_A_INT     - returns element of int array\n
         *  DATA_A_FLOAT   - returns element of float array\n
         *  DATA_A_DOUBLE  - returns element of double array\n
         *  DATA_A_TDS     - returns data element of TDS array\n
         *  DATA_A_BYTE    - returns element of BYTE array\n
         *  DATA_A_USTR    - returns the following fields of USTR [0] element\n
         *                   \ti1_data if index == 0\n
         *                   \tf1_data if index == 1\n
         *                   \tf2_data if index == 2\n
         *                   \ttm      if index == 3\n
         *  DATA_A_XYZS    - returns the following fields of XYZS [0] element\n
         *                   \tstatus if index == 0\n
         *                   \tx      if index == 1\n
         *                   \ty      if index == 2\n
         *                   \tz      if index == 3\n
         *  DATA_IIII      - returns\n
         *                   \ti1_data if index == 0\n
         *                   \ti2_data if index == 1\n
         *                   \ti3_data if index == 2\n
         *                   \ti4_data if index == 3\n
         *  DATA_IFFF      - returns\n
         *                   \ti1_data if index == 0\n
         *                   \tf1_data if index == 1\n
         *                   \tf2_data if index == 2\n
         *                   \tf3_data if index == 3\n
         *  DATA_MDA_FLOAT - returns d_float\n
	 *
	 * @param index  The index of data array element
	 * @return       The double
	**/
	double			get_double         (int index);

	/**
	 *  Returns double array
         *
         *  Applicable to:\n
         *  DATA_A_DOUBLE  - returns pointer to double array\n
         *  others         - returns 0\n
	 *
	 * @return       The array pointer
	**/
	double			*get_double_array  (void);

	/**
	 *  Returns time stamp as an ascii string
	 *
	 * @param bufp   The pointer to buffer where an ascii string will be stored
	 * @param len    The length of buffer
	 * @return       The pointer to string
	**/
	char			*time_string      (char *bufp, size_t len);

	/**
	 *  Returns time stamp field of data element as an ascii string
         *
         *  Applicable to:\n
         *  DATA_TDS       - returns tm element of TDS\n
         *  DATA_TTII      - returns tm1 and tm2 elements of TTII as a string\n
         *  DATA_SPECTRUM  - returns tm element of SPECTRUM array\n
         *  DATA_GSPECTRUM - returns tm element of GSPECTRUM array\n
         *  DATA_A_USTR    - returns tm element of USTR [0] array as a string\n
         *  others         - returns time stamp as an ascii string\n
	 *
	 * @param bufp   The pointer to buffer where an ascii string will be stored
	 * @param len    The length of buffer
	 * @return       The pointer to string
	**/
	char			*get_time_string  (char *bufp, size_t len);

	/**
	 *  Returns data as an ascii string
         *
         *  Applicable to:\n
         *  all data types\n
	 *
	 * @param bufp   The pointer to buffer where an ascii string will be stored
	 * @param len    The length of buffer
	 * @return       The pointer to string
	**/
	char			*get_string       (char *bufp, size_t len);

	/**
	 *  Returns string component of data
         *
         *  Applicable to:\n
         *  DATA_STRING    - returns string\n
         *  DATA_STRING16  - returns short string\n
	        *  DATA_SPECTRUM  - returns comment string of SPECTRUM\n
	        *  DATA_GSPECTRUM - returns comment string of GSPECTRUM\n
         *  DATA_IMAGE     - returns comment string of IMAGE\n
         *  DATA_A_USTR    - returns string component of USTR [0] element\n
         *  DATA_A_XYZS    - returns string component of XYZS [0] element\n
         *  others         - returns null string\n
	 *
	 * @param bufp   The pointer to buffer where an ascii string will be stored
	 * @param len    The length of buffer
	 * @return       The pointer to string
	**/
	char			*get_string_arg   (char *bufp, size_t len);

	/**
	 *  Returns data as an ascii string from specified position
         *  see EqData::get_string ()
	 *
	 * @param index  The position of data
	 * @param bufp   The pointer to buffer where an ascii string will be stored
	 * @param len    The length of buffer
	 * @return       The pointer to string
	**/
	char			*get_string       (int index, char *bufp, size_t len);

	/**
	 *  Returns time stamp field of data element as an ascii string from specified position
         *  see EqData::get_time_string ()
	 *
	 * @param index  The position of data
	 * @param bufp   The pointer to buffer where an ascii string will be stored
	 * @param len    The length of buffer
	 * @return       The pointer to string
	**/
	char			*get_time_string  (int index, char *bufp, size_t len);

	/**
	 *  Returns string part of data from specified position
         *  see EqData::get_string_arg ()
	 *
	 * @param index  The position of data
	 * @param bufp   The pointer to buffer where an ascii string will be stored
	 * @param len    The length of buffer
	 * @return       The pointer to string
	**/
	char			*get_string_arg   (int index, char *bufp, size_t len);

	/**
	 *  Returns data as a string from specified position
         *
         *  Applicable to:\n
         *  DATA_A_TDS\n
	 *
	 * @param index  The position of data
	 * @param bufp   The pointer to buffer where an ascii string will be stored
	 * @param len    The length of buffer
	 * @return       The pointer to string
	**/
	char			*get_string_extra (int index, char *bufp, size_t len);

	/**
	 *  Returns TDS data
         *
         *  Applicable to:\n
         *  DATA_TDS - returns pointer to TDS\n
         *  others   - returns 0\n
	 *
	 * @return       The pointer to data
	**/
	TDS			*get_tds          (void);

	/**
	 *  Returns TDS data from specified position
         *
         *  Applicable to:\n
         *  DATA_TDS - returns pointer to element in TDS array\n
         *  others   - returns 0\n
	 *
	 * @param index  The position of data
	 * @return       The pointer to data
	**/
	TDS			*get_tds          (int index);

	/**
	 *  Returns XY data
         *
         *  Applicable to:\n
         *  DATA_XY  - returns pointer to XY\n
         *  others   - returns 0\n
	 *
	 * @return       The pointer to data
	**/
	XY			*get_xy           (void);

	/**
	 *  Returns XY data from specified position
         *
         *  Applicable to:\n
         *  DATA_XY  - returns pointer to element in XY array\n
         *  others   - returns 0\n
	 *
	 * @param index  The position of data
	 * @return       The pointer to data
	**/
	XY			*get_xy           (int index);

	/**
	 *  Returns IIII data
         *
         *  Applicable to:\n
         *  DATA_IIII - returns pointer to IIII\n
         *  others    - returns 0\n
	 *
	 * @return       The pointer to data
	**/
	IIII			*get_iiii         (void);

	/**
	 *  Returns IFFF data
         *
         *  Applicable to:\n
         *  DATA_IFFF - returns pointer to IFFF\n
         *  others    - returns 0\n
	 *
	 * @return       The pointer to data
	**/
	IFFF			*get_ifff         (void);

	/**
	 *  Returns TTII data
         *
         *  Applicable to:\n
         *  DATA_TTII - returns pointer to TTII\n
         *  others    - returns 0\n
	 *
	 * @return       The pointer to data
	**/
	TTII			*get_ttii         (void);

	/**
	 *  Returns USTR data from specified position
         *
         *  Applicable to:\n
         *  DATA_USTR - returns pointer to element in USTR array\n
         *  others    - returns 0\n
	 *
	 * @param index  The position of data
	 * @return       The pointer to data
	**/
	USTR			*get_ustr         (int index);

	/**
	 *  Returns XYZS data from specified position
         *
         *  Applicable to:\n
         *  DATA_XYZS - returns pointer to element in XYZS array\n
         *  others    - returns 0\n
	 *
	 * @param index  The position of data
	 * @return       The pointer to data
	**/
	XYZS			*get_xyzs         (int index);

	/**
	 *  Returns SPECTRUM data
         *
         *  Applicable to:\n
         *  DATA_SPECTRUM - returns pointer to SPECTRUM\n
         *  others        - returns 0\n
	 *
	 * @return       The pointer to SPECTRUM data
	**/
	SPECTRUM		*get_spectrum     (void);

	/**
	 *  Returns GSPECTRUM data
	 *
	 *  Applicable to:\n
	 *  DATA_GSPECTRUM - returns pointer to GSPECTRUM\n
	 *  others         - returns 0\n
	 *
	 * @return       The pointer to GSPECTRUM data
	 **/
	GSPECTRUM		*get_gspectrum     (void);
	
	/**
	 *  Returns the comment string of SPECTRUM/GSPECTRUM data
	 *
	 *  Applicable to:\n
	 *  DATA_SPECTRUM  - returns comment string to SPECTRUM\n
	 *  DATA_GSPECTRUM - returns comment string to GSPECTRUM\n
     *  others         - returns 0\n
     *
	 * @param bufp   The pointer to buffer where a comment string will be stored
	 * @param len    The length of buffer
	 * @return       The pointer to string
	 **/
	char			*get_spec_comm    (char *bufp, size_t len);

	/**
	 *  Returns the identification of spectrum data
	 *
	 *  Applicable to:\n
	 *  DATA_SPECTRUM  - returns the value of status field of SPECTRUM\n
	 *  DATA_GSPECTRUM - returns the value of id field of GSPECTRUM\n
	 *  others         - returns 0\n
	 *
	 * @param idp    The pointer to buffer where an id will be stored
	 * @return       integer
	 **/
	int			get_spec_id     (u_int *idp);

	/**
	 *  Returns the statistics of grouped spectrum data
	 *
	 *  Applicable to:\n
	 *  DATA_GSPECTRUM - returns the statistics mask and block of GSPECTRUM\n
	 *  others         - returns 0\n
	 *
	 * @param stmskp The pointer to buffer where a statistics mask will be stored
	 * @param stp    The pointer to buffer where a statistics block will be stored
	 * @return       integer
	**/
	int			get_spec_statistics    (int *stmskp, float *stp);

	/**
	 *  Returns the number of groups and group increment in grouped spectrum data
	 *
	 *  Applicable to:\n
	 *  DATA_GSPECTRUM - returns the number of groupes and group increment of GSPECTRUM\n
	 *  others         - returns 0\n
	 *
	 * @param groupsp The pointer to buffer where the number of blocks in spectrum will be stored
	 * @param grpincp The pointer to buffer where the group increment parameter will be stored
	 * @return        The pointer to string
	 **/
	int			get_spec_groups    (int *groupsp, int *grpsize, float *grpincp);

	/**
	 *  Returns the comment field of XML data
	 *
	 *  Applicable to:\n
	 *  DATA_XML - returns pointer to XML array\n
	 *  others   - returns 0\n
	 *
	 * @return       The pointer to data
	**/
	char			*get_xml	       (void);

	/**
	 *  Returns TDS data
	 *
	 *  Applicable to:\n
	 *  DATA_TDS - returns the following elements of TDS in\n
	 *             tm     - tm\n
	 *             data   - dt\n
	 *             status - st\n
	 *
	 * @param tm     The pointer to store tm
	 * @param dt     The pointer to store data
	 * @param st     The pointer to store status
	 * @return       The TDS data
	**/
	int			get_tds           (time_t *tm, float *dt, u_char *st);

	/**
	 *  Returns TDS data from specified position
	 *
	 *  Applicable to:\n
	 *  DATA_TDS - returns the following fields of the element of TDS array in\n
	 *             tm     - tm\n
	 *             data   - dt\n
	 *             status - st\n
	 *
	 * @param index  The position of data
	 * @param tm     The pointer to store tm
	 * @param dt     The pointer to store data
	 * @param st     The pointer to store status
	 * @return       The TDS data
	**/
	int			get_tds      (time_t *tm, float *dt, u_char *st, int index); 

	/**
	 *  Returns XY data
	 *
	 *  Applicable to:\n
	 *  DATA_XY - returns the following elements of XY in\n
	 *            x_data - x\n
	 *            y_data - y\n
	 *
	 * @param x      The pointer to store x_data
	 * @param y      The pointer to store y_data
	 * @return       The XY data
	**/
	int			get_xy         (float *x, float *y);

	/**
	 *  Returns XY data from specified position
	 *
	 *  Applicable to:\n
	 *  DATA_XY - returns the following fields of the element of XY array in\n
	 *            x_data - x\n
	 *            y_data - y\n
	 *
	 * @param index  The position of data
	 * @param x      The pointer to store x_data
	 * @param y      The pointer to store y_data
	 * @return       The XY data
	**/
	int			get_xy         (float *x, float *y, int index);

	/**
	 *  Returns IIII data
	 *
	 *  Applicable to:\n
	 *  DATA_IIII    - returns the following elements of IIII in\n
	 *               \ti1_data - i1\n
	 *               \ti2_data - i2\n
	 *               \ti3_data - i3\n
	 *               \ti4_data - i4\n
	 *  DATA_A_INT   - returns the following elements of int array in\n
	 *               \tint array [0] - i1\n
	 *               \tint array [1] - i2\n
	 *               \tint array [2] - i3\n
	 *               \tint array [3] - i4\n
	 *  DATA_A_FLOAT - returns the following elements of float array in\n
	 *               \tfloat array [0] - i1\n
	 *               \tfloat array [1] - i2\n
	 *               \tfloat array [2] - i3\n
	 *               \tfloat array [3] - i4\n
	 *
	 * @param i1     The pointer to store int
	 * @param i2     The pointer to store int
	 * @param i3     The pointer to store int
	 * @param i4     The pointer to store int
	 * @return       The IIII data
	**/
	int			get_iiii        (int *i1, int *i2, int *i3, int *i4);

	/**
	 *  Returns IFFF data
	 *
	 *  Applicable to:\n
	 *  DATA_IFFF    - returns the following elements of IFFF in\n
	 *               \ti1_data - f1\n
	 *               \tf1_data - f2\n
	 *               \tf2_data - f3\n
	 *               \tf3_data - f4\n
	 *  DATA_A_INT   - returns the following elements of int array in\n
	 *               \tint array [0] - f1\n
	 *               \tint array [1] - f2\n
	 *               \tint array [2] - f3\n
	 *               \tint array [3] - f4\n
	 *  DATA_A_FLOAT - returns the following elements of float array in\n
	 *               \tfloat array [0] - f1\n
	 *               \tfloat array [1] - f2\n
	 *               \tfloat array [2] - f3\n
	 *               \tfloat array [3] - f4\n
	 *
	 * @param i1     The pointer to store int
	 * @param f1     The pointer to store float
	 * @param f2     The pointer to store float
	 * @param f3     The pointer to store float
	 * @return       The IFFF data
	**/
	int			get_ifff         (int *i1, float *f1, float *f2, float *f3);

	/**
	 *  Returns TTII data
	 *
	 *  Applicable to:
	 *  DATA_TTII    - returns the following elements of TTII in\n
	 *               \ttm1     - t1\n
	 *               \ttm2     - t2\n
	 *               \ti1_data - i1\n
	 *               \ti2_data - i2\n
	 *  DATA_A_INT   - returns the following elements of int array in\n
	 *               \tint array [0] - i1\n
	 *               \tint array [1] - i2\n
	 *               \tint array [2] - i3\n
	 *               \tint array [3] - i4\n
	 *
	 * @param t1     The pointer to store time_t
	 * @param t2     The pointer to store time_t
	 * @param i1     The pointer to store int
	 * @param i2     The pointer to store int
	 * @return       The TTII data
	**/
	int			get_ttii        (time_t *t1, time_t *t2, int *i1, int *i2);

	/**
	 *  Returns BYTE data from specified position
	 *
	 *  Applicable to:\n
	 *  DATA_A_BYTE - returns the element of BYTE array\n
	 *
	 * @param index  The position of data
	 * @return       The BYTE data
	**/
	u_char			get_byte       (int index);

	/**
	 *  Returns BYTE data descriptor
	 *
	 *  Applicable to:\n
	 *  DATA_A_BYTE - returns the descriptor of BYTE array in\n
	 *              \tx_dim                 - x\n
	 *              \ty_dim                 - y\n
	 *              \tx_offset              - xo\n
	 *              \ty_offset              - yo\n
	 *              \toption                - opt\n
	 *              \tpointer to BYTE array - cp\n
	 *
	 * @param x      The pointer to store x dimention
	 * @param y      The pointer to store y dimention
	 * @param xo     The pointer to store x offset
	 * @param yo     The pointer to store y offset
	 * @param opt    The pointer to store options
	 * @param cp     The pointer to store pointer to byte array
	 * @return       The descriptor of BYTE array
	**/
	int			get_byte         (int *x, int *y, int *xo, int *yo, int *opt, u_char **cp);


	/**
	 *  Returns pointer byte array and its length
	 *
	 *  Applicable to:\n
	 *  DATA_A_BYTE - returns the length and pointer of BYTE array in len and cp\n
	 *
	 * @param len    The pointer to store int
	 * @param cp     The pointer to store pointer to byte array
	 * @return       The pointer and length in len, cp
	**/
	int			get_byte         (int *len, u_char **cp);

	/**
	 *  Returns SPECTRUM, GSPECTRUM data descriptor
	 *
	 *  Applicable to:\n
	 *  DATA_SPECTRUM - returns the following elements of SPECTRUM descriptor in\n
	 *                \tpointer to comment string        - com\n
	 *                \tlength of comment string         - clen\n
	 *                \ttm                               - tm\n
	 *                \ts_start                          - str\n
	 *                \ts_inc                            - inc\n
	 *                \tstatus                           - st\n
	 *                \tpointer to SPECTRUM float array  - fa\n
	 *                \tlength of SPECTRUM float array   - flen\n
	 *  DATA_GSPECTRUM - returns the following elements of GSPECTRUM descriptor in\n
	 *                \tpointer to comment string        - com\n
	 *                \tlength of comment string         - clen\n
	 *                \ttm                               - tm\n
	 *                \ts_start                          - str\n
	 *                \ts_inc                            - inc\n
	 *                \tstatus                           - st\n
	 *                \tpointer to GSPECTRUM float array - fa\n
	 *                \tlength of GSPECTRUM float array  - flen\n
	 *
	 * @param com    The pointer to store pointer to comment string
	 * @param clen   The pointer to store length of comment string
	 * @param tm     The pointer to store time
	 * @param str    The pointer to store start
	 * @param inc    The pointer to store increment step
	 * @param st     The pointer to store status
	 * @param fa     The pointer to store pointer to spectrum array
	 * @param flen   The pointer to store length of spectrum array
	 * @return       1 for SPECTRUM data descriptor
	 *               2 for GSPECTRUM data descriptor
	 *               0 othewise
	 *
	**/
	int			get_spectrum     (char **com, int *clen, time_t *tm, float *str, float *inc, u_int *st, float **fa, int *flen);

	/**
	 *  Returns USTR data from specified position
	 *
	 *  Applicable to:\n
	 *  DATA_A_USTR - returns the following fields of the element of USTR array\n
	 *              \ti1_data                     - i1\n
	 *              \tf1_data                     - f1\n
	 *              \tf2_data                     - f2\n
	 *              \ttm                          - tm\n
	 *              \tpointer to string component - com\n
	 *
	 * @param index  The position of data
	 * @param i1     The pointer to store i1_data
	 * @param f1     The pointer to store f1_data
	 * @param f2     The pointer to store f2_data
	 * @param tm     The pointer to store tm
	 * @param com    The pointer to store comment string
	 * @return       The USTR data element
         *
	**/
	int			get_ustr         (int *i1, float *f1, float *f2, time_t *tm, char **com, int index);

	/**
	 *  Returns XYZS data from specified position
	 *
	 *  Applicable to:\n
	 *  DATA_A_XYZS - returns the following fields of the element of XYZS array\n
	 *              \tx                     - x\n
	 *              \ty                     - y\n
	 *              \tz                     - z\n
	 *              \tstatus                - st\n
	 *              \tpointer to string loc - com\n
	 *
	 * @param index  The position of data
	 * @param x      The pointer to store x
	 * @param y      The pointer to store y
	 * @param z      The pointer to store z
	 * @param st     The pointer to store status
	 * @param com    The pointer to store loc string
	 * @return       The XYZS data element
         *
	**/
	int			get_xyzs          (int *st, float *x, float *y, float *z, char **com, int index);

	/**
	 *  Returns MDA_FLOAT array
	 *
	 *  Applicable to:\n
	 *  DATA_MDA_FLOAT - returns pointer to MDA_FLOAT\n
	 *  others         - returns 0\n
	 *
	 * @return       The pointer to multidimentional MDA_FLOAT array
	**/
	MDA_FLOAT		*get_mda         (void);

	/**
	 *  Returns the configuration of a multidimentional array (number and length of each dimention)
	 *
	 *  Applicable to:\n
	 *  DATA_A_INT     - \n
	 *  DATA_A_FLOAT   - \n
	 *  DATA_A_TDS     - \n
	 *  DATA_A_XY      - \n
	 *  DATA_A_BYTE    - \n
	 *  DATA_SPECTRUM  - \n
	 *  DATA_GSPECTRUM - \n
	 *  DATA_MDA_FLOAT - \n
	 *
	 * @param dim    The indeces
	 * @param dims   Number of dimentions
	 * @param str    The buffer to store the text description of the array
	 * @param len    The length of the buffer to store the text description of the array
	 * @return       The 1 if OK, 0 if unsupported data type
	**/
	int			get_dims          (int *dim, int *dims, char *str, int len);

	/**
	 *  Returns the floating point data from the specified position
	 *
	 *  Applicable to:\n
	 *  DATA_A_INT     - \n
	 *  DATA_A_FLOAT   - \n
	 *  DATA_A_TDS     - \n
	 *  DATA_SPECTRUM  - \n
	 *  DATA_GSPECTRUM - \n
	 *  DATA_MDA_FLOAT - \n
	 *
	 * @param dim    The indeces
	 * @param dims   Number of dimentions
	 * @param valp   The pointer to store the value of the accessed element of multidimentional array
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			get               (int *dim, int dims, float *valp);

	/**
	 *  Returns the time stamp
	 *
	 *  Applicable to: DATA_IMAGE
	 *
	 * @param psec   The pointer to seconds
	 * @param pusec  The pointer to microseconds
	 * @param pstat  The pointer to status
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int				get_timestamp   (time_t *psec, time_t *pusec, int *pstat);

	/**
	 *  Returns the image data header
	 *
	 *  Applicable to: DATA_IMAGE
	 *
	 * @param valp   The pointer to store the value of the image data header
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_image_header   (IMH *valp);

	/**
	 *  Returns the image data and header
	 *
	 *  Applicable to: DATA_IMAGE
	 *
	 * @param valp   The pointer to store the value of the image data
	 * @param lenp   The pointer to store the length of the image
	 * @param hdp    The pointer to store the header of the image
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_image         (u_char **valp, int *lenp, IMH *hdp = 0);

	/**
	 *  Returns the image comment
	 *
	 *  Applicable to: DATA_IMAGE
	 *
	 * @param bufp   The pointer to store the comment of the image
	 * @param len    The length of the buffer to store the comment
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_image_comm    (char *bufp, size_t len);

	/**
	 *  Sets the length of data
	 *
	 *  Applicable to:\n
	 *  DATA_STRING         - sets length of the string\n
	 *  DATA_STRING16       - sets length of the short string\n
	 *  DATA_XML            - sets length of the XML string\n
	 *  DATA_SPECTRUM       - sets number of elements in the spectrum array\n
	 *  DATA_GSPECTRUM      - sets number of elements in the grouped spectrum array\n
	 *  DATA_IMAGE          - sets number of elements in the image array\n
	 *  DATA_A_SHORT        - sets number of elements in the short array\n
	 *  DATA_A_LONG         - sets number of elements in the long long array\n
	 *  DATA_A_BOOL         - sets number of elements in the bool array\n
	 *  DATA_A_INT          - sets number of elements in the int array\n
	 *  DATA_A_FLOAT        - sets number of elements in the float array\n
	 *  DATA_A_TDS          - sets number of elements in the TDS array\n
	 *  DATA_A_XY           - sets number of elements in the XY array\n
	 *  DATA_A_USTR         - sets number of elements in the USTR array\n
	 *  DATA_A_XYZS         - sets number of elements in the XYZS array\n
	 *  DATA_A_THUMBNAIL    - sets number of elements in THUMBNAIL array\n
	 *  DATA_A_BYTE         - sets number of elements in the BYTE array\n
	 *  DATA_A_TS_BOOL      - sets number of elements in the TS_BOOL array\n
	 *  DATA_A_TS_INT       - sets number of elements in the TS_INT array\n
	 *  DATA_A_TS_FLOAT     - sets number of elements in the TS_FLOAT array\n
	 *  DATA_A_TS_DOUBLE    - sets number of elements in the TS_DOUBLE array\n
	 *  DATA_A_TS_STRING16  - sets number of elements in the TS_STRING16 array\n
	 *  DATA_A_TS_STRING    - sets number of elements in the TS_STRING array\n
	 *  DATA_A_TS_USTR      - sets number of elements in the TS_USTR array\n
	 *  DATA_A_TS_XML       - sets number of elements in the TS_XML array\n
	 *  DATA_A_TS_XY        - sets number of elements in the TS_XY array\n
	 *  DATA_A_TS_IIII      - sets number of elements in the TS_IIII array\n
	 *  DATA_A_TS_IFFF      - sets number of elements in the TS_IFFF array\n
	 *  DATA_A_TS_TTII      - sets number of elements in the TS_TTII array\n
	 *  DATA_A_TS_XYZS      - sets number of elements in the TS_XYZS array\n
	 *  DATA_A_TS_SPECTRUM  - sets number of elements in the TS_SPECTRUM array\n
	 *  DATA_A_TS_GSPECTRUM - sets number of elements in the TS_GSPECTRUM array\n
	 *  DATA_KEYVAL         - sets number of elements in the KEYVAL table\n
	 *  DATA_IIII           - sets 4\n
	 *  DATA_IFFF           - sets 4\n
	 *  others              - sets 1\n
	 *
	 * @param len  The length of data
	 *
	**/
	int			length            (int len);

	/**
	 *  Stores time stamp of data
	 *
	 * @param tm  The time stamp
	 *
	**/
	int			time              (time_t tm);

	/**
	 *  Stores cycle mask of data
	 *
	 * @param cycle  The cycle mask
	 *
	**/
	int			cycle_mask        (int cycle);

	/**
	 *  Stores error code of data
	 *
	 * @param err  The error code
	 *
	**/
	int			error             (int err);

	/**
	 *  Stores error code of data and
	 *  sets DATA_STRING data type and stores error string
	 *
	 * @param err     The error code
	 * @param errstr  The error description string
	 *
	**/
	int			error             (int err, char *errstr);

	/**
	 *  Sets data type
	 *
	 * @param dtype  The data type
	 *
	**/
	int			set_type          (int dtype);

	/**
	 *  Sets DATA_BOOL data type and stores a value
	 *
	 * @param dt  The boolean value
	 *
	**/
	int			set_bool          (int dt);

	/**
	 *  Stores int field of data at specified position
	 *
	 *  Applicable to:\n
	 *       DATA_A_SHORT   - stores int at specified position in short array\n
	 *
	 * @param index  The position of data
	 * @param dt     The short value
	 *
	**/
	int			set               (short dt, int index);

	/**
	 *  Stores int field of data at specified position
	 *
	 *  Applicable to:\n
	 *       DATA_A_LONG    - stores int at specified position in long long array\n
	 *
	 * @param index  The position of data
	 * @param dt     The long long value
	 *
	**/
	int			set               (long long dt, int index);

	/**
	 *  Sets DATA_INT data type and stores a value
	 *
	 * @param dt  The int value
	 *
	**/
	int			set               (int dt);

	/**
	 *  Stores int field of data at specified position
	 *
	 *  Applicable to:\n
	 *       DATA_A_INT     - stores int at specified position in int_array\n
	 *       DATA_A_BOOL    - stores int at specified position in bool array\n
	 *       DATA_A_SHORT   - stores int at specified position in short array\n
	 *       DATA_A_LONG    - stores int at specified position in long long array\n
	 *       DATA_A_FLOAT   - stores int at specified position in float array\n
	 *       DATA_A_DOUBLE  - stores int at specified position in double array\n
	 *       DATA_A_USTR    - stores int in\n
	 *                        \ti1_data if index == 0\n
	 *                        \tf1_data if index == 1\n
	 *                        \tf2_data if index == 2\n
	 *                        \ttm      if index == 3\n
	 *       DATA_A_XYZS    - stores int in\n
	 *                        \tstatus if index == 0\n
	 *                        \tx      if index == 1\n
	 *                        \ty      if index == 2\n
	 *                        \tz      if index == 3\n
	 *       DATA_MDA_FLOAT - stores int in multidimentional array\n
	 *       DATA_IMAGE     - stores int in image\n
	 *
	 * @param index  The position of data
	 * @param dt     The int value
	 *
	**/
	int			set               (int dt, int index);

	/**
	 *  Sets DATA_FLOAT data type and stores a value
	 *
	 * @param dt     The float value
	 *
	**/
	int			set               (float dt);

	/**
	 *  Stores float value at specified position
	 *
	 *  Applicable to:\n
	 *       DATA_A_INT     - stores int at specified position in int array\n
	 *       DATA_A_FLOAT   - stores float at specified position in float array\n
	 *       DATA_A_DOUBLE  - stores double at specified position in double array\n
	 *       DATA_A_USTR    - stores int in\n
	 *                        \ti1_data if index == 0\n
	 *                        \tf1_data if index == 1\n
	 *                        \tf2_data if index == 2\n
	 *                        \ttm      if index == 3\n
	 *       DATA_A_XYZS    - stores int in\n
	 *                        \tstatus if index == 0\n
	 *                        \tx      if index == 1\n
	 *                        \ty      if index == 2\n
	 *                        \tz      if index == 3\n
	 *       DATA_SPECTRUM  - stores float at specified position in spectrum array\n
	 *       DATA_GSPECTRUM - stores float at specified position in gspectrum array\n
	 *       DATA_MDA_FLOAT - stores float in multidimentional array\n
	 *
	 * @param index  The position of data
	 * @param dt     The float value
	 *
	**/
	int			set               (float dt, int index);

	/**
	 *  Sets DATA_DOUBLE data type and stores a value
	 *
	 * @param dt     The double value
	 *
	**/
	int			set               (double dt);

	/**
	 *  Stores double value at specified position
	 *
	 *  Applicable to:\n
	 *       DATA_A_INT     - stores int at specified position in int array\n
	 *       DATA_A_FLOAT   - stores float at specified position in float array\n
	 *       DATA_A_DOUBLE  - stores double at specified position in double array\n
	 *       DATA_A_USTR    - stores int in\n
	 *                        \ti1_data if index == 0\n
	 *                        \tf1_data if index == 1\n
	 *                        \tf2_data if index == 2\n
	 *                        \ttm      if index == 3\n
	 *       DATA_A_XYZS    - stores int in\n
	 *                        \tstatus if index == 0\n
	 *                        \tx      if index == 1\n
	 *                        \ty      if index == 2\n
	 *                        \tz      if index == 3\n
	 *       DATA_SPECTRUM  - stores float at specified position in spectrum array\n
	 *       DATA_GSPECTRUM - stores float at specified position in gspectrum array\n
	 *       DATA_MDA_FLOAT - stores float in multidimentional array\n
	 *
	 * @param index  The position of data
	 * @param dt     The double value
	 *
	**/
	int			set               (double dt, int index);

	/**
	 *  Stores a string in data and sets DATA_STRING or\n
	 *  DATA_TEXT data type depending on the length of\n
	 *  the input string\n
	 *
	 * @param str  The string pointer
	 *
	**/
	int			set               (char *str);

	/**
	 *  Depending on the data type reads the input string and
	 *  fills the binary fields of data structure with converted
	 *  from string values. It does not change data type
	 *
	 * @param str  The string pointer
	 *
	**/
	int			set_from_string   (char *str);

	/**
	 *  Sets DATA_TDS data type and stores a value
	 *
	 * @param p      The pointer to TDS (if 0, null-TDS is stored)
	 *
	**/
	int			set               (TDS *p);

	/**
	 *  Sets DATA_TDS data type and stores values
	 *
	 * @param tm      The tm value
	 * @param dt      The data value
	 * @param st      The status value
	 *
	**/
	int			set               (time_t tm, float dt, u_char st);

	/**
	 *  Sets DATA_A_TDS data type and stores a value at specified position
	 *
	 * @param index  The position of data
	 * @param p      The pointer to TDS (if 0, null-TDS is stored)
	 *
	**/
	int			set               (TDS *p, int index);

	/**
	 *  Sets DATA_A_TDS data type and stores values at specified position
	 *
	 * @param index   The position of data
	 * @param tm      The tm value
	 * @param dt      The data value
	 * @param st      The status value
	 *
	**/
	int			set               (time_t tm, float dt, u_char st, int index);

	/**
	 *  Sets DATA_SPECTRUM data type and stores a value
	 *
	 * @param p      The pointer to SPECTRUM (if 0, null-SPECTRUM is stored)
	 *
	**/
	int			set               (SPECTRUM *p);

	/**
	 *  Sets DATA_GSPECTRUM data type and stores a value
	 *
	 * @param p      The pointer to GSPECTRUM (if 0, null-GSPECTRUM is stored)
	 *
	 **/
	int			set               (GSPECTRUM *p);
	
	/**
	 *  Sets the identification of spectrum data
	 *
	 *  Applicable to:\n
	 *  DATA_SPECTRUM  - sets the status filed of SPECTRUM\n
	 *  DATA_GSPECTRUM - sets the id field of GSPECTRUM\n
	 *  others         - returns 0\n
	 *
	 * @param id     The id
	 * @return       integer
	 **/
	int			set_spec_id        (u_int id);

	/**
	 *  Sets the statistics of grouped spectrum data
	 *
	 *  Applicable to:\n
	 *  DATA_GSPECTRUM - sets the statistics mask and block of GSPECTRUM\n
	 *  others         - returns 0\n
	 *
	 * @param stmsk  The statistics mask
	 * @param stp    The pointer to statistics data block to be stored
	 * @return       integer
	 **/
	int			set_spec_statistics  (int stmsk, float *stp);
	
	/**
	 *  Sets the number of groups and group increment in grouped spectrum data
	 *
	 *  Applicable to:\n
	 *  DATA_GSPECTRUM - sets the number of groupes and group increment of GSPECTRUM\n
	 *  others         - returns 0\n
	 *
	 * @param groups  The number of groups in spectrum
	 * @param grpinc  The group increment parameter
	 * @return        integer
	 **/
	int			set_spec_groups    (int groups, int grpsize, float grpinc);
	
	/**
	 *  Sets DATA_SPECTRUM data type and stores values
	 *
	 * @param com   The comment string (if 0 comment set to "")
	 * @param tm    The tm value
	 * @param str   The s_start value
	 * @param inc   The s_inc value
	 * @param st    The status value
	 * @param fp    The pointer to spectrum array copied to d_spect_array (if 0 all values of d_spect_array set to 0.0)
	 * @param flen  The length of spectrum array
	 *
	**/
	int			set       (char *com, time_t tm, float str, float inc, u_int st, float *fp, int flen);

	/**
	 *  Sets DATA_GSPECTRUM data type and stores values,
	 *  statistic mask and data block are zeroed
	 *
	 * @param id     The identification
	 * @param grp    The number of groups
	 * @param grpinc The group increment
	 * @param com    The comment string (if 0 comment set to "")
	 * @param tm     The tm value
	 * @param str    The s_start value
	 * @param inc    The s_inc value
	 * @param st     The status value
	 * @param fp     The pointer to spectrum array copied to d_spect_array (if 0 all values of d_spect_array set to 0.0)
	 * @param flen   The length of spectrum array
	 *
	 **/
	int			set        (u_int id, int grp, int grpsize, float grpinc, char *com, time_t tm, float str, float inc, u_int st, float *fp, int flen);
	
	/**
	 *  Sets DATA_A_USTR data type and stores a value at specified position
	 *
	 * @param index   The position of data
	 * @param p       The pointer to USTR (if 0, null-USTR is stored)
	 *
	**/
	int			set            (USTR *p, int index);

	/**
	 *  Sets DATA_A_USTR data type and stores values at specified position
	 *
	 * @param index   The position of data
	 * @param i1      The i1_data value
	 * @param f1      The f1_data value
	 * @param f2      The f2_data value
	 * @param tm      The tm value
	 * @param str     The pointer to string copied to str_data (if str 0 str_data set to "")
	 *
	**/
	int			set             (int i1, float f1, float f2, time_t tm, char *str, int index);

	/**
	 *  Sets DATA_XY data type and stores a value
	 *
	 * @param p       The pointer to XY (if 0, null-XY is stored)
	 *
	**/
	int			set               (XY *p);

	/**
	 *  Sets DATA_XY data type and stores values
	 *
	 * @param x       The x_data value
	 * @param y       The y_data value
	 *
	**/
	int			set               (float x, float y);

	/**
	 *  Sets DATA_A_XY data type and stores a value at specified position
	 *
	 * @param index   The position of data
	 * @param p       The pointer to XY (if 0, null-XY is stored)
	 *
	**/
	int			set               (XY *p, int index);

	/**
	 *  Sets DATA_A_XY data type and stores values at specified position
	 *
	 * @param index   The position of data
	 * @param x       The x_data value
	 * @param y       The y_data value
	 *
	**/
	int			set               (float x, float y, int index);

	/**
	 *  Sets DATA_IIII data type and stores a value
	 *
	 * @param p       The pointer to IIII (if 0, null-IIII is stored)
	 *
	**/
	int			set               (IIII *p);

	/**
	 *  Sets DATA_IIII data type and stores a value
	 *
	 * @param i1      The i1_data value
	 * @param i2      The i2_data value
	 * @param i3      The i3_data value
	 * @param i4      The i4_data value
	 *
	**/
	int			set               (int i1, int i2, int i3, int i4);

	/**
	 *  Sets DATA_IFFF data type and stores a value
	 *
	 * @param p       The pointer to IFFF (if 0, null-IFFF is stored)
	 *
	**/
	int			set               (IFFF *p);

	/**
	 *  Sets DATA_IFFF data type and stores a value
	 *
	 * @param i1      The i1_data value
	 * @param f1      The f1_data value
	 * @param f2      The f2_data value
	 * @param f3      The f3_data value
	 *
	**/
	int			set               (int i1, float f1, float f2, float f3);

	/**
	 *  Sets DATA_TTII data type and stores a value
	 *
	 * @param p       The pointer to TTII (if 0, null-TTII is stored)
	 *
	**/
	int			set               (TTII *p);

	/**
	 *  Sets DATA_TTII data type and stores a value
	 *
	 * @param tm1     The tm1 value
	 * @param tm2     The tm2 value
	 * @param i1      The i1_data value
	 * @param i2      The i2_data value
	 *
	**/
	int			set               (time_t tm1, time_t tm2, int i1, int i2);

	/**
	 *  Sets DATA_A_XYZS data type and stores a value at specified position
	 *
	 * @param index   The position of data
	 * @param p       The pointer to XYZS (if 0, null-XYZS is stored)
	 *
	**/
	int			set               (XYZS *p, int index);

	/**
	 *  Sets DATA_A_XYZS data type and stores values at specified position
	 *
	 * @param index   The position of data
	 * @param st      The status value
	 * @param x       The x value
	 * @param y       The y value
	 * @param z       The z value
	 * @param str     The pointer to string copied to loc (if str 0, loc is set to "")
	 *
	**/
	int			set_xyzs          (int st, float x, float y, float z, char *str, int index);

	/**
	 *  Sets DATA_A_BYTE data type and stores a value
	 *
	 * @param x       The x_dim value
	 * @param y       The y_dim value
	 * @param xo      The x_offset value
	 * @param yo      The y_offset value
	 * @param opt     The option
	 * @param ba      The pointer to byte array copied to d_byte_array (if ba 0, d_byte_array filled with 0)
	 *
	**/
	int			set_byte          (int x, int y, int xo, int yo, int opt, u_char *ba);

	/**
	 *  Sets DATA_A_BYTE data type and stores a byte array
	 *
	 * @param len     The length of byte array
	 * @param ba      The pointer to byte array (if ba 0, d_byte_array filled with 0)
	 *
	**/
	int			set_byte          (size_t len, u_char *ba);

	/**
	 *  Sets DATA_XML data type and stores a byte array
	 *
	 * @param str     The pointer to XML array
	 *
	**/
	int			set_xml	       (char *str);

	/**
	 *  Configure a miltidimentional array (number and length of each dimention)
	 *
	 *  Applicable to:\n
	 *  DATA_MDA_FLOAT  - \n
	 *
	 * @param dim    The indeces
	 * @param dims   Number of dimentions
	 * @param str    The buffer to store the text description of the array
	 * @param val    The initial value of the configured array
	 * @return       The 1 if OK, 0 if unsupported data type
	**/
	int			set_dims          (int *dim, int dims, char *str, float val = HUGE);

	/**
	 *  Store the value into the addressed position of the multidimentional array
	 *
	 *  Applicable to:\n
	 *  DATA_A_INT     - \n
	 *  DATA_A_FLOAT   - \n
	 *  DATA_A_TDS     - \n
	 *  DATA_SPECTRUM  - \n
	 *  DATA_GSPECTRUM  - \n
	 *  DATA_MDA_FLOAT - \n
	 *
	 * @param dim    The indeces
	 * @param dims   Number of dimentions
	 * @param val    The value to be assigned to the element of multidimentional array
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set               (int *dim, int dims, float val);

	/**
	 *  Extract the Time-Status-Data from the addressed position of the history array
	 *
	 *  Applicable to:\n
	 *  DATA_A_TS_BOOL      - \n
	 *  DATA_A_TS_INT       - \n
	 *  DATA_A_TS_FLOAT     - \n
	 *  DATA_A_TS_DOUBLE    - \n
	 *  DATA_A_TS_STRING16  - \n
	 *  DATA_A_TS_STRING    - \n
	 *  DATA_A_TS_USTR      - \n
	 *  DATA_A_TS_XML       - \n
	 *  DATA_A_TS_XY        - \n
	 *  DATA_A_TS_IIII      - \n
	 *  DATA_A_TS_IFFF      - \n
	 *  DATA_A_TS_TTII      - \n
	 *  DATA_A_TS_XYZS      - \n
	 *  DATA_A_TS_SPECTRUM  - \n
	 *  DATA_A_TS_GSPECTRUM - \n
	 *
	 * @param tsp   The pointer to time (sec) field
	 * @param tmsp  The pointer to time (msec) field
	 * @param sp    The pointer to status field
	 * @param vp    The pointer to the object where the extracted element of the history array stored
	 *              It should point to one of the types: BOOL, INT, FLOAT, DOUBLE, STRING16,
	 *                          STRING, USTR, XML, XY, IIII, IFFF, TTII, XYZS, SPECTRUM
	 * @param index The index (position)
	 * @return      The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_tsd           (u_int *tsp, u_short *tmsp, u_short *sp, void *vp, int index = 0);

	/**
	 *  Extract the Time-Status-Bool from the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_BOOL
	 *
	 * @param tsp   The pointer to time (sec) field
	 * @param tmsp  The pointer to time (msec) field
	 * @param sp    The pointer to status field
	 * @param vp    The pointer to the object where the extracted element of the history array stored
	 * @param index The index (position)
	 * @return      The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_ts_bool       (u_int *tsp, u_short *tmsp, u_short *sp, bool *vp, int index = 0);

	/**
	 *  Extract the Time-Status-Int from the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_INT
	 *
	 * @param tsp   The pointer to time (sec) field
	 * @param tmsp  The pointer to time (msec) field
	 * @param sp    The pointer to status field
	 * @param vp    The pointer to the object where the extracted element of the history array stored
	 * @param index The index (position)
	 * @return      The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_ts_int        (u_int *tsp, u_short *tmsp, u_short *sp, int *vp, int index = 0);

	/**
	 *  Extract the Time-Status-Float from the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_FLOAT
	 *
	 * @param tsp   The pointer to time (sec) field
	 * @param tmsp  The pointer to time (msec) field
	 * @param sp    The pointer to status field
	 * @param vp    The pointer to the object where the extracted element of the history array stored
	 * @param index The index (position)
	 * @return      The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_ts_float      (u_int *tsp, u_short *tmsp, u_short *sp, float *vp, int index = 0);

	/**
	 *  Extract the Time-Status-Double from the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_DOUBLE
	 *
	 * @param tsp   The pointer to time (sec) field
	 * @param tmsp  The pointer to time (msec) field
	 * @param sp    The pointer to status field
	 * @param vp    The pointer to the object where the extracted element of the history array stored
	 * @param index The index (position)
	 * @return      The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_ts_double     (u_int *tsp, u_short *tmsp, u_short *sp, double *vp, int index = 0);

	/**
	 *  Extract the Time-Status-String from the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_STRING
	 *
	 * @param tsp   The pointer to time (sec) field
	 * @param tmsp  The pointer to time (msec) field
	 * @param sp    The pointer to status field
	 * @param vp    The pointer to the object where the extracted element of the history array stored
	 * @param index The index (position)
	 * @return      The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_ts_str        (u_int *tsp, u_short *tmsp, u_short *sp, char *vp, int index = 0);

	/**
	 *  Extract the Time-Status-USTR from the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_USTR
	 *
	 * @param tsp   The pointer to time (sec) field
	 * @param tmsp  The pointer to time (msec) field
	 * @param sp    The pointer to status field
	 * @param vp    The pointer to the object where the extracted element of the history array stored
	 * @param index The index (position)
	 * @return      The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_ts_ustr       (u_int *tsp, u_short *tmsp, u_short *sp, USTR *vp, int index = 0);

	/**
	 *  Extract the Time-Status-XML from the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_XML
	 *
	 * @param tsp   The pointer to time (sec) field
	 * @param tmsp  The pointer to time (msec) field
	 * @param sp    The pointer to status field
	 * @param vp    The pointer to the object where the extracted element of the history array stored
	 * @param index The index (position)
	 * @return      The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_ts_xml        (u_int *tsp, u_short *tmsp, u_short *sp, XML *vp, int index = 0);

	/**
	 *  Extract the Time-Status-XY from the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_XY
	 *
	 * @param tsp   The pointer to time (sec) field
	 * @param tmsp  The pointer to time (msec) field
	 * @param sp    The pointer to status field
	 * @param vp    The pointer to the object where the extracted element of the history array stored
	 * @param index The index (position)
	 * @return      The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_ts_xy         (u_int *tsp, u_short *tmsp, u_short *sp, XY *vp, int index = 0);

	/**
	 *  Extract the Time-Status-IIII from the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_IIII
	 *
	 * @param tsp   The pointer to time (sec) field
	 * @param tmsp  The pointer to time (msec) field
	 * @param sp    The pointer to status field
	 * @param vp    The pointer to the object where the extracted element of the history array stored
	 * @param index The index (position)
	 * @return      The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_ts_iiii       (u_int *tsp, u_short *tmsp, u_short *sp, IIII *vp, int index = 0);

	/**
	 *  Extract the Time-Status-IFFF from the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_IFFF
	 *
	 * @param tsp   The pointer to time (sec) field
	 * @param tmsp  The pointer to time (msec) field
	 * @param sp    The pointer to status field
	 * @param vp    The pointer to the object where the extracted element of the history array stored
	 * @param index The index (position)
	 * @return      The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_ts_ifff       (u_int *tsp, u_short *tmsp, u_short *sp, IFFF *vp, int index = 0);

	/**
	 *  Extract the Time-Status-XYZS from the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_XYZS
	 *
	 * @param tsp   The pointer to time (sec) field
	 * @param tmsp  The pointer to time (msec) field
	 * @param sp    The pointer to status field
	 * @param vp    The pointer to the object where the extracted element of the history array stored
	 * @param index The index (position)
	 * @return      The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_ts_xyzs       (u_int *tsp, u_short *tmsp, u_short *sp, XYZS *vp, int index = 0);

	/**
	 *  Extract the Time-Status-SPECTRUM from the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_SPECTRUM
	 *
	 * @param tsp   The pointer to time (sec) field
	 * @param tmsp  The pointer to time (msec) field
	 * @param sp    The pointer to status field
	 * @param vp    The pointer to the object where the extracted element of the history array stored
	 * @param index The index (position)
	 * @return      The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_ts_spect      (u_int *tsp, u_short *tmsp, u_short *sp, SPECTRUM *vp, int index = 0);

	/**
	 *  Extract the Time-Status-GSPECTRUM from the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_GSPECTRUM
	 *
	 * @param tsp   The pointer to time (sec) field
	 * @param tmsp  The pointer to time (msec) field
	 * @param sp    The pointer to status field
	 * @param vp    The pointer to the object where the extracted element of the history array stored
	 * @param index The index (position)
	 * @return      The 1 if OK, 0 if unsupported data type or access
	 **/
	int			get_ts_gspect     (u_int *tsp, u_short *tmsp, u_short *sp, GSPECTRUM *vp, int index = 0);
    
	/**
	 *  Get the list of all registered keys in KEYVAL table
	 *
	 *  Applicable to: DATA_KEYVAL
	 *
	 * @param keyp  The pointer to table allocated for a list of keys (an entry per key)
	 * @param lenp  The pointer to the size of table allocated for a list of keys
	 * @return      The 1 if OK, 0 if unsupported data type or access violation
	 **/
         int             get_keys           (char ***keyp, int *lenp);

	/**
	 *  Get the value of the requested key in KEYVAL table
	 *
	 *  Applicable to: DATA_KEYVAL
	 *
	 * @param keyp  The pointer to the requested key
	 * @param valp  The pointer to buffer for the value of the requested key
	 * @param len   The length of the provided buffer for the requested value
	 * @return      The 1 if OK, 0 if unsupported data type or access
	 **/
         int             get_val            (char *keyp, char *valp, int len);

	/**
	 *  Store the Time-Status-Data into the addressed position of the history array
	 *
	 *  Applicable to:\n
	 *  DATA_A_TS_BOOL      - \n
	 *  DATA_A_TS_INT       - \n
	 *  DATA_A_TS_FLOAT     - \n
	 *  DATA_A_TS_DOUBLE    - \n
	 *  DATA_A_TS_STRING16  - \n
	 *  DATA_A_TS_STRING    - \n
	 *  DATA_A_TS_USTR      - \n
	 *  DATA_A_TS_XML       - \n
	 *  DATA_A_TS_XY        - \n
	 *  DATA_A_TS_IIII      - \n
	 *  DATA_A_TS_IFFF      - \n
	 *  DATA_A_TS_TTII      - \n
	 *  DATA_A_TS_XYZS      - \n
	 *  DATA_A_TS_SPECTRUM  - \n
	 *  DATA_A_TS_GSPECTRUM - \n
	 *
	 * @param tsec   Time (sec)
	 * @param tmsec  Time (msec)
	 * @param stat   Status
	 * @param vp     The address of the structure to be assigned to the element of the history array
	 *               It should point to one of the types: BOOL, INT, FLOAT, DOUBLE, STRING16,
	 *                           STRING, USTR, XML, XY, IIII, IFFF, TTII, XYZS, SPECTRUM
	 * @param index  The index (position)
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_tsd           (u_int tsec, u_short tmsec, u_short stat, void *vp, int index = 0);

	/**
	 *  Store the Time-Status-Bool into the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_BOOL
	 *
	 * @param tsec   Time (sec)
	 * @param tmsec  Time (msec)
	 * @param stat   Status
	 * @param vp     The value to be assigned to the element of the history array
	 * @param index  The index (position)
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_ts_bool       (u_int tsec, u_short tmsec, u_short stat, bool vp, int index = 0);

	/**
	 *  Store the Time-Status-Int into the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_INT
	 *
	 * @param tsec   Time (sec)
	 * @param tmsec  Time (msec)
	 * @param stat   Status
	 * @param vp     The value to be assigned to the element of the history array
	 * @param index  The index (position)
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_ts_int        (u_int tsec, u_short tmsec, u_short stat, int vp, int index = 0);

	/**
	 *  Store the Time-Status-Float into the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_FLOAT
	 *
	 * @param tsec   Time (sec)
	 * @param tmsec  Time (msec)
	 * @param stat   Status
	 * @param vp     The value to be assigned to the element of the history array
	 * @param index  The index (position)
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_ts_float      (u_int tsec, u_short tmsec, u_short stat, float vp, int index = 0);

	/**
	 *  Store the Time-Status-Double into the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_DOUBLE
	 *
	 * @param tsec   Time (sec)
	 * @param tmsec  Time (msec)
	 * @param stat   Status
	 * @param vp     The value to be assigned to the element of the history array
	 * @param index  The index (position)
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_ts_double     (u_int tsec, u_short tmsec, u_short stat, double vp, int index = 0);

	/**
	 *  Store the Time-Status-String into the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_STRING
	 *
	 * @param tsec   Time (sec)
	 * @param tmsec  Time (msec)
	 * @param stat   Status
	 * @param vp     The pointer to string to be assigned to the element of the history array
	 * @param index  The index (position)
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_ts_str        (u_int tsec, u_short tmsec, u_short stat, char *vp, int index = 0);

	/**
	 *  Store the Time-Status-USTR into the addressed position of the history array
	 *
	 *  Applicable to:DATA_A_TS_USTR
	 *
	 * @param tsec   Time (sec)
	 * @param tmsec  Time (msec)
	 * @param stat   Status
	 * @param vp     The pointer to USTR object to be assigned to the element of the history array
	 * @param index  The index (position)
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_ts_ustr       (u_int tsec, u_short tmsec, u_short stat, USTR *vp, int index = 0);

	/**
	 *  Store the Time-Status-XML into the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_XML
	 *
	 * @param tsec   Time (sec)
	 * @param tmsec  Time (msec)
	 * @param stat   Status
	 * @param vp     The pointer to XML object to be assigned to the element of the history array
	 * @param index  The index (position)
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_ts_xml        (u_int tsec, u_short tmsec, u_short stat, XML *vp, int index = 0);

	/**
	 *  Store the Time-Status-XML into the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_XY
	 *
	 * @param tsec   Time (sec)
	 * @param tmsec  Time (msec)
	 * @param stat   Status
	 * @param vp     The pointer to XY object to be assigned to the element of the history array
	 * @param index  The index (position)
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_ts_xy         (u_int tsec, u_short tmsec, u_short stat, XY *vp, int index = 0);

	/**
	 *  Store the Time-Status-IIII into the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_IIII
	 *
	 * @param tsec   Time (sec)
	 * @param tmsec  Time (msec)
	 * @param stat   Status
	 * @param vp     The pointer to IIII object to be assigned to the element of the history array
	 * @param index  The index (position)
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_ts_iiii       (u_int tsec, u_short tmsec, u_short stat, IIII *vp, int index = 0);

	/**
	 *  Store the Time-Status-IFFF into the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_IFFF
	 *
	 * @param tsec   Time (sec)
	 * @param tmsec  Time (msec)
	 * @param stat   Status
	 * @param vp     The pointer to IFFF object to be assigned to the element of the history array
	 * @param index  The index (position)
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_ts_ifff       (u_int tsec, u_short tmsec, u_short stat, IFFF *vp, int index = 0);

	/**
	 *  Store the Time-Status-XYZS into the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_XYZS
	 *
	 * @param tsec   Time (sec)
	 * @param tmsec  Time (msec)
	 * @param stat   Status
	 * @param vp     The pointer to XYZS object to be assigned to the element of the history array
	 * @param index  The index (position)
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_ts_xyzs       (u_int tsec, u_short tmsec, u_short stat, XYZS *vp, int index = 0);

	/**
	 *  Store the Time-Status-SPECTRUM into the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_SPECTRUM
	 *
	 * @param tsec   Time (sec)
	 * @param tmsec  Time (msec)
	 * @param stat   Status
	 * @param vp     The pointer to SPECTRUM object to be assigned to the element of the history array
	 * @param index  The index (position)
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_ts_spect      (u_int tsec, u_short tmsec, u_short stat, SPECTRUM *vp, int index = 0);

	/**
	 *  Store the Time-Status-GSPECTRUM into the addressed position of the history array
	 *
	 *  Applicable to: DATA_A_TS_GSPECTRUM
	 *
	 * @param tsec   Time (sec)
	 * @param tmsec  Time (msec)
	 * @param stat   Status
	 * @param vp     The pointer to GSPECTRUM object to be assigned to the element of the history array
	 * @param index  The index (position)
	 * @return       The 1 if OK, 0 if unsupported data type or access
	 **/
	int			set_ts_gspect     (u_int tsec, u_short tmsec, u_short stat, GSPECTRUM *vp, int index = 0);
	
	/**
	 *  Store the time stamp
	 *
	 *  Applicable to: DATA_IMAGE
	 *
	 * @param sec    The time (sec)
	 * @param usec   The time (usec)
	 * @param stat   The status
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_timestamp   (time_t sec, time_t usec, int stat);

	/**
	 *  Store the image data header
	 *
	 *  Applicable to: DATA_IMAGE
	 *
	 * @param valp   The pointer to the header to be assigned to the image data header
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_image_header   (IMH *valp);

	/**
	 *  Store the image data and header
	 *
	 *  Applicable to: DATA_IMAGE
	 *
	 * @param valp   The pointer to image data
	 * @param len    The length of the image data
	 * @param hdp    The pointer to the header to be assigned to the image data header
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_image          (u_char *valp, int len, IMH *hdp = 0);

	/**
	 *  Store the image comment
	 *
	 *  Applicable to: DATA_IMAGE
	 *
	 * @param bufp   The pointer to image comment
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			set_image_comm     (char *bufp);

	/**
	 *  Returns a thumbnail
	 *
	 *  Applicable to: DATA_A_THUMBNAIL
	 *
	 * @param pt     The pointer to thumbnail
	 * @param index  The position of thumbnail
	 * @return       The 1 if OK, 0 if unsupported data type or access
	**/
	int			get_thumbnail      (THUMBNAIL *pt, int index);

	/**
	 *  Store a thumbnail
	 *
	 *  Applicable to: DATA_A_THUMBNAIL
	 *
	 * @param pt     The pointer to thumbnail
	 * @param index  The position of thumbnail
	**/
	int			set_thumbnail      (THUMBNAIL *pt, int index);
    
	/**
	 *  Clear all keys of the table (clear table)
	 *
	 *  Applicable to: DATA_KEYVAL
	 *
	 * @return       The 1 if OK, 0 if unsupported data type
	**/
        int             clear_keys         (void);

	/**
	 *  Store a key-value pair
	 *
	 *  Applicable to: DATA_KEYVAL
	 *
	 * @param keyp   The pointer to the requested key
	 * @param valp   The pointer to the value of the requested key
	 * @return       The 1 if OK, 0 if unsupported data type or table full
	**/
        int             set_val            (char *keyp, char *valp);
};


// ========================================================================================
//. 
//:DESCRIPTION
//.The EqData class is the interface to the data structures transfered from the client
//.program to the device server and the returned data to the client. It contains an
//.error flag, a time stamp, the data type and the data length and a variable sized
//.amount of data of different type. The data is always send in is native form and
//.converted on the receiver. One may request a floating point value on a client in
//.one call and read the result as a float, int or string without redoing the transfer.
//.Some functions return an int as a return code: 1 is success, 0 is an error.
//!EqData ();
//!EqData (EqDataBlock *);
//.	Constructor. The second form is not used in client programs. The parameter is a pointer
//.	to the internal C structure which is transported over the network.
//.	Not recommended to use as it will be removed in the next versions.
//!~EqData ();
//.	Destructor.
//.
//!void  load_simple (EqData* from)
//.	Loads this EqData from the data of the from argument. The
//.	base data, the time stamp and error fields are copied.
//.	Base data means all data types with up to
//.	4 words (no arrays), i.e. short data types.
//!void	copy_from (EqData *from)
//.	Copies the data of from into this EqData.
//!void	add_from (EqData *from)
//.	Adds the data of from to this EqData structure. Works for all
//.	array data types. The first call to an initialized EqData defines
//.	the data type. In further add_from calls the data type of the
//.	from argument has to be identical, other wise nothing is added.
//.	If from has an error set the error code  is set in this EqData and
//.	no data is appended.
//.
//!time_t time ()	
//!int    time (time_t)
//.	Get and set the time stamp of the data block.
//!char *time_string (char *bufp, int len);
//.	Get the time stamp as an ASCII string in supplied external buffer.
//!int error ()
//!int error (int)
//!int error (int, char *)
//.	Get and set the error code of the data block. The second version set the corresponding 
//.	character string also. Setting of the error is mainly done on the server side.
//!int type ()
//!int set_type (int ty)
//.	Get and set the data type of the data block. The data type is set automaticaly in most of
//.	the set oparations. To destinguish between a fill of a float array or spectrum one has to
//.	specify the type before writing a float.
//.
//!Data Types:
//.
//. DATA_NULL		: no data
//. DATA_INT		: single integer
//. DATA_BOOL		: single boolean
//. DATA_FLOAT		: return "FLOAT";
//. DATA_STRING		: string of charaters
//. DATA_STRING16	: string of charaters with max 16 chars.
//. DATA_TDS		: time, data and status structure
//. DATA_XY		: two floats
//. DATA_IIII		: 4 integer
//. DATA_IFFF		: one integer and 3 floats
//. DATA_TTII		: 2 times and 2 integer
//. DATA_SPECTRUM	: spectrum structure
//. DATA_XML		: spectrum structure
//. DATA_A_FLOAT	: float array
//. DATA_A_INT		: integer array
//. DATA_A_TDS		: array of time, data and status structure
//. DATA_A_XY		: array of two floats
//. DATA_A_USTR		: array of "universal" strings, includes int and two floats
//. DATA_A_BYTE		: byte array
//. DATA_A_XYZS		: array of <int status:float x:float y:float z:string location name>
//. DATA_A_BOOL         : history array of type time, status and bool
//. DATA_A_INT          : history array of type time, status and int
//. DATA_A_FLOAT        : history array of type time, status and float
//. DATA_A_DOUBLE       : history array of type time, status and double
//. DATA_A_USTR         : history array of type time, status and universal string
//. DATA_A_XML          : history array of type time, status and xml string
//. DATA_A_XY           : history array of type time, status and xy structure
//. DATA_A_IIII         : history array of type time, status and iiii structure
//. DATA_A_IFFF         : history array of type time, status and ifff structure
//. DATA_A_TTII         : history array of type time, status and ttii structure
//. DATA_A_XYZS         : history array of type time, status and xyzs structure
//. DATA_A_SPECTRUM     : history array of type time, status and spectrum structure
//.
//!char *type_string ();
//.	Get the data type as ASCII name of current data block.
//!char *type_string (int)
//.	Get the data type as an ASCII name from int arg.
//!EqDataBlock* data_block ()
//.	Not recommended to use. Only for internal processing.
//!int  data_block (EqDataBlock *db)
//.	Get and set the data block pointer. The pointer points to a C structure.
//.	Not recommended to use. Only for internal processing.
//!int length ()
//!int length (int)
//.	Get and set the data length. If the data is an array of a struct the array length is returned.
//.	Setting of the length is normally done by writing into an array. The highest index defines
//.	the length.
//!int array_length ()
//.	Returns the length of arrays only. If the data type is not an
//.	array type 0 is returned.
//!int string_length (int index = 0)
//.	Returns the length of the string part. If no string is part of the
//.	data type 0 is returned. index is the array index for array data types.

//!int  init ()
//.	Reset content and set data type to DATA_NULL. Used on server only.

//!int  get_int ()
//!int  set (int)
//.	Get and set the data as int. If the native data is not int it is converted to int in the get call.

//!int  set_bool (int)
//.	Set the data block to a single boolean.
//!float get_float ()
//!int   set (float)
//.	Get and set the data as float. If the native data is not float it is converted to float in the
//.	get call.
//!float get_float (int index)
//!int   set (float data, int index)
//.	Get and set the data as float from/to index position. The data type must be a spectrum or 
//.	float array.
//!float *get_float_array ()
//.	Get the data as a pointer to the float array. The data type must be a spectrum or float array.
//!char *get_string (char *bufp, int len)
//.	Get the data as an ASCII string in supplied external buffer.
//!int   set (char *)
//.	Get and set the data as/from an ASCII string. If the native data is not float it is converted 
//.	to ASCII in the get call.
//!char *get_string (int index, char *bufp, int len)
//	Get the data from index position as an ASCII string in supplied external buffer.
//.	Data must be an array type.
//!char *get_string_extra (int index, char *bufp, int len)
//.	Like get_string but gets a second time point for plots if data type is A_TDS as an 
//.	ASCII string in supplied external buffer.
//!int  set_from_string (char *)
//.	Set the data block from a string. The data type must be defined before with set_type.
//!int  set_from_string (char *, int index)
//.	Set the data block at index position from a string. The data type must be defined before  
//.	with set_type.
//!char *get_time_string (char *bufp, int len)
//.	Get time as ASCII string in supplied external buffer.
//.	If the data has no time info: return the time stamp.
//!char *get_time_string (int index, char *bufp, int len)
//.	Get the time as an ASCII string from index position in supplied external buffer.
//.	If the data has no array with a time info: return the time stamp.
//!char *get_string_arg (char *bufp, int len)
//!char *get_string_arg (int index, char *bufp, int len)
//.	Get the string part of the data in supplied external buffer.
//.	STRING, USTR and SPECTRUM are the only data types with a string part.
//.	For other data types a zero length string is returned. With the
//.	index argument: get the string part from index position.
//!TDS  *get_tds ()
//!int  set (TDS *)
//.	Get a or set TimeDataStatus pointer. The TDS structure contains a time stamp, a float  
//.	value and a status word.
//!int  get_tds (time_t *, float *, u_char *)
//!int  set (time_t, float, u_char);
//.	Get or set the Time, Data and Status value of the TDS structure.
//!TDS  *get_tds (int index)
//!int  set (TDS *, int index);
//.	Get or set a TimeDataStatus pointer from/to index position.
//!int  get_tds (time_t *, float *, u_char *, int index)
//!int  set (time_t, float, u_char, int index);
//.	Get or set the Time, Data and Status value from/to index position.

//!SPECTRUM *get_spectrum ()
//!int  set (SPECTRUM *)
//.	Get or set a spectrum pointer.
//!int  get_spectrum (char **, int *, time_t *, float *, float *, u_char *, float **, int *)
//!int  set (char *, time_t, float, float, u_char, float *, int len);
//.	Get or set all values from/to a spectrum. A spectrum structure contains an ASCII  
//.	string as a comment and the length of this string, a time stamp, a float as the 
//.	start value and a float as the increment, a status info, a float array and the
//.	length of the float array. If a pointer is provided in a set call the data is 
//.	NOT copied. The second form does a copy.
//!int	get_spec_comm (char *comm);
//.	Get the comment of the spectrum.
	
//!USTR *get_ustr (int index)
//!int  set (USTR *, int index);/* set the data block to string array
//.	Get or set the pointer to the index position in a string array. The USTR data type 
//.	is a universal string which contains an int arg (type info e.g.), two float args 
//.	(x and y value e.g.), a  time stamp and the character string with its length. 
//.	The USTR data type is used in the  names call to retrieve property, location 
//.	and device info. The set call copies the data.
//!int  get_ustr (int *, float *, float *, time_t *, char **, int index);
//!int  set (int, float, float, time_t, char *, int index);
//.	Get and set the values at/to index position from/to a USTR array data type. 
//.	The set call copies the data.

//!XY   *get_xy ()
//!int  get_xy (float *x, float *y)
//!XY   *get_xy (int index)
//!int  get_xy (float *, float *, int index)
//.	Get a x/y value pair as a pointer to a XY structure or as single items. If the data 
//.	type is an array index is used to point to the entry.

//!IIII *get_iiii ()
//!int   set (IIII *)
//!int   get_iiii (int *, int *, int *, int *)
//!int   set (int, int, int, int)
//.	Get and set a pointer to the 4 integers or get/set the individual integers. 
//.	The IIII data type is used to hold device addresses (line, crete, module..) 
//.	for instance.

//!IFFF *get_ifff ()
//!int  set (IFFF *)
//!int  get_ifff (int *, float *, float *, float *)
//!int  set (int, float, float, float)
//.	Get or set a pointer to one integer and 3 floats or get/set the individual values. 
//.	The IFFF data type is used to set filter parameters or polynoms for instance.

//!TTII *get_ttii ()
//!int  set (TTII *)
//!int  get_ttii (time_t *, time_t *, int *, int *)
//!int  set (time_t, time_t, int, int)
//.	Get or set a pointer to the 2 times and 2 integers or get/set the individual values. 
//.	The TTII data type is used to specify a time range in historical readings for instance.

//!u_char get_byte (int index)
//!int  get_byte (int *ix, int *iy, int *ixo, int *iyo, int *iopt, u_char *cp)
//!int  set_byte (int  ix, int  iy, int  ixo, int  iyo, int  iopt, u_char *cp)
//.	Get or set a byte array. The byte array can by one or two dimensional. 
//.	The latter version is used to transfer images. The int args
//.	specify the length of the x and y dimension, the offset in x and y direction and
//	an optional parameter. No buffer is allocated
//.	by the set_byte call! The length is calculated from the x and y dimension
//.	args. and cp points to the start of the char buffer.

	
#endif
