#ifndef _EQ_TYPES_H_RPCGEN
#define	_EQ_TYPES_H_RPCGEN

#include <rpc/rpc.h>

#ifndef _KERNEL
#ifndef __GNUC__

#include <synch.h>
#include <thread.h>

#endif
#endif /* !_KERNEL */

#ifdef __cplusplus

extern "C" {

#endif

#ifndef eq_types_h
#define eq_types_h


#if defined (__APPLE__) || defined (__SunOS_5_6)

typedef   u_long    rpcprog_t;
typedef   u_long    rpcvers_t;
typedef   u_long    rpcproc_t;
    
#endif


/* protocol defined data field sizes */

#define	FACILITY_STRING_LENGTH  64
#define	DEVICE_STRING_LENGTH    64
#define	LOCATION_STRING_LENGTH  64
#define	PROPERTY_STRING_LENGTH  64

#define HOST_STRING_LENGTH      256
#define AUTH_STRING_LENGTH      128
#define SERVER_STRING_LENGTH    128
#define FILE_STRING_LENGTH      256


#define	DATA_NULL               0
#define	DATA_INT                1
#define	DATA_FLOAT              2
#define	DATA_STRING             3
#define	DATA_BOOL               4
#define	DATA_STRING16           5
#define	DATA_DOUBLE             6
#define	DATA_TEXT               7

#define	DATA_TDS                12
#define	DATA_XY                 13
#define	DATA_IIII               14
#define	DATA_IFFF               15
#define	DATA_USTR               16
#define	DATA_TTII               18
#define	DATA_SPECTRUM           19
#define	DATA_XML                20
#define	DATA_XYZS               21
#define	DATA_IMAGE              22
#define	DATA_GSPECTRUM          24

#define	DATA_A_FLOAT            100
#define	DATA_A_TDS              101
#define	DATA_A_XY               102
#define	DATA_A_USTR             103
#define	DATA_A_INT              105
#define	DATA_A_BYTE             106
#define	DATA_A_XYZS             108
#define	DATA_MDA_FLOAT          109
#define	DATA_A_DOUBLE           110
#define	DATA_A_BOOL             111
#define	DATA_A_STRING           112
#define DATA_A_SHORT            113
#define DATA_A_LONG             114

#define	DATA_A_THUMBNAIL        120

/*
  not yet implemented in C++ API

#define	DATA_A_NAME16FI	        111
#define	DATA_A_NAME16II	        112
#define	DATA_A_NAME16	        113
#define DATA_A_FI               114
#define DATA_A_SHORT            115
#define DATA_A_II               116
#define DATA_A_NAME64           117
#define DATA_A_NAME32           118
#define DATA_A_BOOL             119
*/

#define	DATA_A_TS_BOOL          1000
#define	DATA_A_TS_INT           1001
#define	DATA_A_TS_FLOAT         1002
#define	DATA_A_TS_DOUBLE        1003
#define	DATA_A_TS_STRING        1005
#define	DATA_A_TS_USTR          1006
#define	DATA_A_TS_XML           1007
#define	DATA_A_TS_XY            1008
#define	DATA_A_TS_IIII          1009
#define	DATA_A_TS_IFFF          1010
#define	DATA_A_TS_SPECTRUM      1013
#define	DATA_A_TS_XYZS          1014
#define DATA_A_TS_GSPECTRUM     1015

#define	DATA_KEYVAL             1016


/* configuration constants */

#define	SHORT_STRING_LENGTH     16
#define	STRING_LENGTH           80
#define	TEXT_LENGTH             1048576
#define	RB_LENGTH               2000
#define	TDS_LENGTH              2000
#define	BOOL_LENGTH             8192
#define	SHORT_LENGTH            262144  /* old value 8192 */
#define	LONG_LENGTH             262144  /* old value 8192 */
#define	INT_LENGTH              262144  /* old value 8192 */
#define	FLOAT_LENGTH            262144  /* old value 8192 */
#define	DOUBLE_LENGTH           131072  /* old value 8192 */
#define	IMAGE_LENGTH            8388608
#define	BYTE_ARRAY_LENGTH       8388608
#define	XML_ARRAY_LENGTH        4194304
#define	XY_LENGTH               8192
#define	SPECTRUM_LENGTH         8192
#define	SPECTRUM_LENGTH_MAX     1048576
#define	USTR_LENGTH             4000
#define	STRING_ARRAY_LENGTH     128
#define	XYZS_LENGTH             1000
#define	THUMBNAIL_LENGTH        1000
#define	THUMBNAIL_LENGTH_MAX    10000
#define	MDA_LENGTH              16384
#define	GSPECTRUM_LENGTH        4194304
#define KEYVAL_LENGTH           1000
#define KEYVAL_STRING_LENGTH    1000

#define	DIM_MAX                 6

#define	TS_BOOL_LENGTH          8192
#define	TS_INT_LENGTH           8192
#define	TS_FLOAT_LENGTH         8192
#define	TS_DOUBLE_LENGTH        8192
#define	TS_STRING_LENGTH        8192
#define	TS_USTR_LENGTH          1000
#define	TS_XML_LENGTH           128
#define	TS_XY_LENGTH            8192
#define	TS_IIII_LENGTH          8192
#define	TS_IFFF_LENGTH          8192
#define	TS_SPECTRUM_LENGTH      128
#define	TS_XYZS_LENGTH          8192
#define	TS_GSPECTRUM_LENGTH     64


#define CHANNEL_STRING_LENGTH   ((LOCATION_STRING_LENGTH) + \
                                 (PROPERTY_STRING_LENGTH))
	
#ifdef time_t
#undef time_t

typedef long time_t;

#endif


/* data structures */


enum {
    
    SPEC_STAT_MIN = 0,
    SPEC_STAT_MAX,
    SPEC_STAT_MEAN,
    SPEC_STAT_RMS,
    SPEC_STAT_SIGMA,
    SPEC_STAT_PP_ON,
    SPEC_STAT_PP_A0,
    SPEC_STAT_PP_A1,
    SPEC_STAT_PP_A2,
    SPEC_STAT_SPARE0,
    SPEC_STAT_SPARE1,
    SPEC_STAT_SPARE2,
    SPEC_STAT_SPARE3,
    SPEC_STAT_SPARE4,
    SPEC_STAT_SPARE5,
    SPEC_STAT_SPARE6,
    SPEC_STAT_SPARE7,
    SPEC_STAT_SPARE8,
    SPEC_STAT_SPARE9,
    SPEC_STAT_SPARE10,
    SPEC_STAT_SPARE11,
    SPEC_STAT_SPARE12,
    SPEC_STAT_SPARE13,
    SPEC_STAT_SPARE14,
    SPEC_STAT_SPARE15,
    SPEC_STAT_SPARE16,
    SPEC_STAT_SPARE17,
    SPEC_STAT_SPARE18,
    SPEC_STAT_SPARE19,
    SPEC_STAT_SPARE20,
    SPEC_STAT_SPARE21,
    SPEC_STAT_SPARE22
};

struct TDS {
	time_t tm;
	float data;
	u_char status;
};
typedef struct TDS TDS;

struct XY {
	float x_data;
	float y_data;
};
typedef struct XY XY;

struct IIII {
	int i1_data;
	int i2_data;
	int i3_data;
	int i4_data;
};
typedef struct IIII IIII;

struct IFFF {
	int i1_data;
	float f1_data;
	float f2_data;
	float f3_data;
};
typedef struct IFFF IFFF;

struct XML {
	struct {
		u_int d_xml_array_len;
		char *d_xml_array_val;
	} d_xml_array;
};
typedef struct XML XML;

struct USTR {
	int i1_data;
	float f1_data;
	float f2_data;
	time_t tm;
	struct {
		u_int str_data_len;
		char *str_data_val;
	} str_data;
};
typedef struct USTR USTR;

struct TTII {
	time_t tm1;
	time_t tm2;
	int i1_data;
	int i2_data;
};
typedef struct TTII TTII;

struct XYZS {
	int status;
	float x;
	float y;
	float z;
	struct {
		u_int loc_len;
		char *loc_val;
	} loc;
};
typedef struct XYZS XYZS;

struct SPECTRUM {
	struct {
		u_int comment_len;
		char *comment_val;
	} comment;
	time_t tm;
	float s_start;
	float s_inc;
	u_int status;
	struct {
		u_int d_spect_array_len;
		float *d_spect_array_val;
	} d_spect_array;
};
typedef struct SPECTRUM SPECTRUM;

struct MDA_FLOAT {
	struct {
		u_int comment_len;
		char *comment_val;
	} comment;
	int dimn;
	int dims[DIM_MAX];
	struct {
		u_int d_mdfa_len;
		float *d_mdfa_val;
	} d_mdfa;
};
typedef struct MDA_FLOAT MDA_FLOAT;

struct A_BYTE {
	int x_dim;
	int y_dim;
	int x_offset;
	int y_offset;
	int option;
	struct {
		u_int d_byte_array_len;
		u_char *d_byte_array_val;
	} d_byte_array;
};
typedef struct A_BYTE A_BYTE;

struct IMH {
	int width;
	int height;
	int aoi_width;
	int aoi_height;
	int x_start;
	int y_start;
	int bpp;
	int ebitpp;
	int hbin;
	int vbin;
	int source_format;
	int image_format;
	u_int frame;
	u_int event;
	float scale_x;
	float scale_y;
	float image_rotation;
	float fspare2;
	float fspare3;
	float fspare4;
	int image_flags;
	int ispare2;
	int ispare3;
	int ispare4;
	int length;
};
typedef struct IMH IMH;

struct IMAGE {
	IMH hdr;
	time_t sec;
	time_t usec;
	int status;
	struct {
		u_int comment_len;
		char *comment_val;
	} comment;
	struct {
		u_int val_len;
		u_char *val_val;
	} val;
};
typedef struct IMAGE IMAGE;

struct THUMBNAIL {
	double average;
	double deviation;
	u_int count;
	u_int time_sec_0;
	u_int time_sec_min;
	u_int time_sec_max;
	double data_min;
	double data_max;
};
typedef struct THUMBNAIL THUMBNAIL;

	
struct GSPECTRUM {
	struct {
		u_int comment_len;
		char *comment_val;
	} comment;
	u_int id;
	time_t tm;
	float s_start;
	float s_inc;
	u_int status;
	int groups;
        int group_size;
	float group_inc;
	int stat_mask;
	float statistics [32];
	struct {
		u_int d_gspect_array_len;
		float *d_gspect_array_val;
	} d_gspect_array;
};
typedef struct GSPECTRUM GSPECTRUM;

struct TSH {
	u_int   tsec;
	u_short tmsec;
	u_short status;
};
typedef struct TSH TSH;

struct TS_BOOL {
	struct TSH hdr;
	bool_t val;
};
typedef struct TS_BOOL TS_BOOL;

struct TS_INT {
	struct TSH hdr;
	int val;
};
typedef struct TS_INT TS_INT;

struct TS_FLOAT {
	struct TSH hdr;
	float val;
};
typedef struct TS_FLOAT TS_FLOAT;

struct TS_DOUBLE {
	struct TSH hdr;
	double val;
};
typedef struct TS_DOUBLE TS_DOUBLE;

struct TS_STRING {
	struct TSH hdr;
	struct {
		u_int val_len;
		char *val_val;
	} val;
};
typedef struct TS_STRING TS_STRING;

struct TS_USTR {
	struct TSH hdr;
	struct USTR val;
};
typedef struct TS_USTR TS_USTR;

struct TS_XML {
	struct TSH hdr;
	struct XML val;
};
typedef struct TS_XML TS_XML;

struct TS_XY {
	struct TSH hdr;
	struct XY val;
};
typedef struct TS_XY TS_XY;

struct TS_IIII {
	struct TSH hdr;
	struct IIII val;
};
typedef struct TS_IIII TS_IIII;

struct TS_IFFF {
	struct TSH hdr;
	struct IFFF val;
};
typedef struct TS_IFFF TS_IFFF;

struct TS_SPECTRUM {
	struct TSH hdr;
	struct SPECTRUM val;
};
typedef struct TS_SPECTRUM TS_SPECTRUM;

struct TS_XYZS {
	struct TSH hdr;
	struct XYZS val;
};
typedef struct TS_XYZS TS_XYZS;

struct TS_GSPECTRUM {
	struct TSH hdr;
	struct GSPECTRUM val;
};
typedef struct TS_GSPECTRUM TS_GSPECTRUM;

struct KEYVAL {
	struct {
		u_int val_len;
		char *val_val;
	} val;
};
typedef struct KEYVAL KEYVAL;



struct DataUnion {
	int data_sel;
	union {
		int       d_bool;
		int       d_int;
		float     d_float;
		double    d_double;
		TDS       d_tds;
		XY        d_xy;
		IIII      d_iiii;
		IFFF      d_ifff;
		TTII      d_ttii;
		SPECTRUM  d_spectrum;
		XML       d_xml;
		XYZS      d_xyzs;
		IMAGE     d_image;
		GSPECTRUM d_gspectrum;
		MDA_FLOAT d_mdfa;
		A_BYTE    d_byte_struct;
		char      d_char16 [SHORT_STRING_LENGTH];
		struct {
			u_int d_char_len;
			char *d_char_val;
		} d_char;
		struct {
			u_int d_short_array_len;
			short *d_short_array_val;
		} d_short_array;
		struct {
			u_int d_llong_array_len;
			long long *d_llong_array_val;
		} d_llong_array;
		struct {
			u_int d_int_array_len;
			int *d_int_array_val;
		} d_int_array;
                struct {
			u_int d_float_array_len;
			float *d_float_array_val;
		} d_float_array;
		struct {
			u_int d_double_array_len;
			double *d_double_array_val;
		} d_double_array;
		struct {
			u_int d_tds_array_len;
			TDS *d_tds_array_val;
		} d_tds_array;
		struct {
			u_int d_xy_array_len;
			XY *d_xy_array_val;
		} d_xy_array;
		struct {
			u_int d_ustr_array_len;
			USTR *d_ustr_array_val;
		} d_ustr_array;
		struct {
			u_int d_xyzs_array_len;
			XYZS *d_xyzs_array_val;
		} d_xyzs_array;
		struct {
			u_int d_thumbnail_array_len;
			THUMBNAIL *d_thumbnail_array_val;
		} d_thumbnail_array;
		struct {
			u_int d_ts_bool_array_len;
			TS_BOOL *d_ts_bool_array_val;
		} d_ts_bool_array;
		struct {
			u_int d_ts_int_array_len;
			TS_INT *d_ts_int_array_val;
		} d_ts_int_array;
		struct {
			u_int d_ts_float_array_len;
			TS_FLOAT *d_ts_float_array_val;
		} d_ts_float_array;
		struct {
			u_int d_ts_double_array_len;
			TS_DOUBLE *d_ts_double_array_val;
		} d_ts_double_array;
		struct {
			u_int d_ts_char_array_len;
			TS_STRING *d_ts_char_array_val;
		} d_ts_char_array;
		struct {
			u_int d_ts_ustr_array_len;
			TS_USTR *d_ts_ustr_array_val;
		} d_ts_ustr_array;
		struct {
			u_int d_ts_xml_array_len;
			TS_XML *d_ts_xml_array_val;
		} d_ts_xml_array;
		struct {
			u_int d_ts_xy_array_len;
			TS_XY *d_ts_xy_array_val;
		} d_ts_xy_array;
		struct {
			u_int d_ts_iiii_array_len;
			TS_IIII *d_ts_iiii_array_val;
		} d_ts_iiii_array;
		struct {
			u_int d_ts_ifff_array_len;
			TS_IFFF *d_ts_ifff_array_val;
		} d_ts_ifff_array;
		struct {
			u_int d_ts_spectrum_array_len;
			TS_SPECTRUM *d_ts_spectrum_array_val;
		} d_ts_spectrum_array;
		struct {
			u_int d_ts_xyzs_array_len;
			TS_XYZS *d_ts_xyzs_array_val;
		} d_ts_xyzs_array;
		struct {
			u_int d_ts_gspectrum_array_len;
			TS_GSPECTRUM *d_ts_gspectrum_array_val;
		} d_ts_gspectrum_array;
		struct {
			u_int d_keyval_array_len;
			KEYVAL *d_keyval_array_val;
		} d_keyval_array;
	} DataUnion_u;
};
typedef struct DataUnion DataUnion;

struct EqDataBlock {
	time_t tm;
	int error;
	DataUnion data_u;
};
typedef struct EqDataBlock EqDataBlock;

struct EqNameString {
	struct {
		u_int n_facility_len;
		char *n_facility_val;
	} n_facility;
	struct {
		u_int n_location_len;
		char *n_location_val;
	} n_location;
	struct {
		u_int n_object_len;
		char *n_object_val;
	} n_object;
	struct {
		u_int n_property_len;
		char *n_property_val;
	} n_property;
};
typedef struct EqNameString EqNameString;


#endif  /* eq_types_h */


/* the xdr functions */

extern  bool_t xdr_time_t(XDR *, time_t*);
extern  bool_t xdr_TDS(XDR *, TDS*);
extern  bool_t xdr_XY(XDR *, XY*);
extern  bool_t xdr_SPECTRUM(XDR *, SPECTRUM*);
extern  bool_t xdr_MDA_FLOAT(XDR *, MDA_FLOAT*);
extern  bool_t xdr_A_BYTE(XDR *, A_BYTE*);
extern  bool_t xdr_XML(XDR *, XML*);
extern  bool_t xdr_IIII(XDR *, IIII*);
extern  bool_t xdr_IFFF(XDR *, IFFF*);
extern  bool_t xdr_USTR(XDR *, USTR*);
extern  bool_t xdr_TTII(XDR *, TTII*);
extern  bool_t xdr_XYZS(XDR *, XYZS*);
extern  bool_t xdr_IMH(XDR *, IMH*);
extern  bool_t xdr_IMAGE(XDR *, IMAGE*);
extern  bool_t xdr_THUMBNAIL(XDR *, THUMBNAIL*);
extern  bool_t xdr_GSPECTRUM(XDR *, GSPECTRUM*);
extern  bool_t xdr_TSH(XDR *, TSH*);
extern  bool_t xdr_TS_BOOL(XDR *, TS_BOOL*);
extern  bool_t xdr_TS_INT(XDR *, TS_INT*);
extern  bool_t xdr_TS_FLOAT(XDR *, TS_FLOAT*);
extern  bool_t xdr_TS_DOUBLE(XDR *, TS_DOUBLE*);
extern  bool_t xdr_TS_STRING(XDR *, TS_STRING*);
extern  bool_t xdr_TS_USTR(XDR *, TS_USTR*);
extern  bool_t xdr_TS_XML(XDR *, TS_XML*);
extern  bool_t xdr_TS_XY(XDR *, TS_XY*);
extern  bool_t xdr_TS_IIII(XDR *, TS_IIII*);
extern  bool_t xdr_TS_IFFF(XDR *, TS_IFFF*);
extern  bool_t xdr_TS_SPECTRUM(XDR *, TS_SPECTRUM*);
extern  bool_t xdr_TS_XYZS(XDR *, TS_XYZS*);
extern  bool_t xdr_TS_GSPECTRUM(XDR *, TS_GSPECTRUM*);
extern  bool_t xdr_DataUnion(XDR *, DataUnion*);
extern  bool_t xdr_EqDataBlock(XDR *, EqDataBlock*);
extern  bool_t xdr_EqNameString(XDR *, EqNameString*);

#ifdef __cplusplus
}
#endif

#endif /* !_EQ_TYPES_H_RPCGEN */
