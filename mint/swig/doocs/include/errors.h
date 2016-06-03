/**< \file errors.h
 *
 * system error code definitions
 *
 */
#ifndef ERRORS_H_
#define ERRORS_H_

#ifndef NO_DOXYGEN
#include "project.def"
#if defined(__cplusplus) && !defined(FORCE_CPP)
extern "C" {
#endif

#ifndef _char
# define _char char
#endif

extern _char erlst[][32];

#define cc2str(c) erlst[(c)&0xff]
#endif /* NO_DOXYGEN */

#ifndef TINE_DECORATED_CONSTANTS

#define CE_REDIRECTED  0x8000
#define CE_SENDDATA    0x4000
#define CE_RENEGOTIATE 0x2000
#define CE_SENDCODE    0x1000

#else

#define TDC_CE_REDIRECTED  0x8000
#define TDC_CE_SENDDATA    0x4000
#define TDC_CE_RENEGOTIATE 0x2000

#endif

#ifndef tine_decorated_constants

/** \def Error Number Definitions */
#define operation_success          0  /**< successful completion */
#define illegal_line               1  /**< illegal line number */
#define illegal_format             2  /**< illegal format or data type */
#define illegal_arithmetic         3  /**< illegal arithmetic expression */
#define ambiguous                  4  /**< ambigous input information */
#define illegal_delimiter          5  /**< illegal delimiter in expression */
#define zero_divide                6  /**< divide by zero */
#define work_area_full             7  /**< working area buffer is full */
#define non_existent               8  /**< requested target does not exist */
#define invalid_transport          9  /**< transport medium is not valid */
#define data_size_mismatch         10 /**< incoming data payload larger than requested ? */

#define not_properly               11 /**< expression not properly terminated (unused) */
#define un_allocated               12 /**< requested target is not allocated in the given context */
#define no_such_line               13 /**< expected line not found in database */
#define illegal_data_size          14 /**< give data size is illegal or not valid */
#define io_error                   15 /**< bus or network io error */
#define illegal_context            16 /**< given context is illegal */
#define runtime_error              17 /**< runtime error */
#define system_error               18 /**< system internal error */
#define hardware_busy              19 /**< hardware is busy */
#define argument_list_error        20 /**< input argument or parameter is invalid (alias for invalid_parameter */
#define invalid_parameter          20 /**< input argument or parameter is invalid */

#define file_error                 21 /**< file access error */
#define use_stream_transport       22 /**< datagram message size exceeds the allowable mtu (reserved) */
#define dimension_error            23 /**< array dimension error */
#define square_root_negative       24 /**< square root of a negative number */
#define buffer_too_small           25 /**< supplied buffer is too small to handle the request */
#define string_too_long            26 /**< supplied string is too long to buffer and cannot be truncated */
#define ipx_socket_error           27 /**< ipx socket error (reserved) */
#define net_read_error             28 /**< network read error */
#define not_ready                  29 /**< call to equipment module dispatch is not ready (internal, reserved) */
#define invalid_transport_size     30 /**< total message size exceeds the configured allowable message size */

#define log_negative               31 /**< log of a negative or zero number */
#define device_not_connected       32 /**< requested device is not connected */
#define unauthorised_action        33 /**< requested operation is not authorized */
#define hardware_error             34 /**< hardware error */
#define illegal_equipment_number   35 /**< requested device does not have a table entry (alias for illegal_device_number) */
#define illegal_device_number      35 /**< requested device does not have a table entry */
#define illegal_device             35 /**< requested device is illegal */
#define illegal_property           36 /**< requested property is illegal */
#define out_of_range               37 /**< given value is out of range */
#define not_implemented            38 /**< requested operation has not been implemented */
#define no_such_computer           39 /**< requested host computer does not exist */
#define struct_sealed              40 /**< tagged structure has already been sealed (reserved) */

#define syntax_error               41 /**< expression syntax error (unused) */
#define no_such_file               42 /**< requested file does not exist */
#define file_already_exists        43 /**< requested file already exists */
#define no_file_space              44 /**< not enough capacity in pre-allocated random access file to complete the requested operation */
#define link_not_open              45 /**< presistent link has timed out (alias for link_timeout) */
#define link_timeout               45 /**< presistent link has timed out */
#define remitted_data_lost         46 /**< returned only a portion of the requested data set (reserved) */
#define end_of_file                47 /**< attempt to read past the end of a file (unused) */
#define archive_busy               48 /**< archive operation is in progress */
#define server_name_in_use         49 /**< device server name is already being used */
#define no_such_column             50 /**< expected column not found in database */

#define out_of_client_memory       51 /**< a client call could not allocate sufficient memory */
#define mcast_access_required      52 /**< access requires CA_NETWORK which was not supplied! */
#define address_unresolved         53 /**< requested target address could not be resolved */
#define invalid_property           54 /**< requested property is not properly formed */
#define address_unknown            55 /**< requested target address could not be found */
#define name_unknown               56 /**< requested name could not be found */
#define invalid_keyword            57 /**< input globals keyword is invalid */
#define invalid_link               58 /**< input link id is invalid */
#define string_expected            59 /**< database file read encountered an end-of-line while expecting a string */       
#define out_of_local_memory        60 /**< a local call to allocated memory was not successful */

#define out_of_shared_memory       61 /**< a call to obtain a shared-memory buffer was not successful */
#define invalid_structure_tag      62 /**< the given structure tag is not in the structure registry */
#define invalid_index              63 /**< the index given does not reference a table entry */
#define illegal_equipment_name     64 /**< the equipment name given is not valid */
#define link_error                 65 /**< internal error processing data link (reserved) */
#define code_failure               66 /**< internal failure due to unhandled exceptions, unavailable resource, OS failure, etc. */
#define non_existent_server        67 /**< requested server does not exist */
#define function_deprecated        68 /**< requested action or functionality has been deprecated */
#define not_supported              69 /**< requested action or functionality not supported */
#define illegal_character          70 /**< error parsing input */
#define parsing_error              70 /**< error parsing input */

#define illegal_operator           71 /**< input arithmetic or logical operator is not allowed */
#define not_allowed                72 /**< requested operation is not allowed */
#define illegal_read_write         73 /**< requested access is not allowed (typically CA_WRITE) */
#define out_of_server_memory       74 /**< a server call could not allocate sufficient memory */
#define database_not_loaded        75 /**< requested operation requires a database which has not been loaded */
#define illegal_command            76 /**< requested command is illegal */
#define resources_exhausted        77 /**< requested operation would exceed the configured capacity of an allocation table (clients, contracts, connections, etc.) */
#define file_not_open              78 /**< requested file io could not be made as the file is not open (unused) */
#define sedac_error                79 /**< sedac bus error */
#define reacquire_address          80 /**< signal to reacquire a server's address from ENS (reserved, internal) */

#define semaphore_error            81 /**< a semaphore or mutex could not be obtained or released */
#define driver_not_installed       82 /**< request requires a hardware driver which is not installed */
#define port_not_available         83 /**< request requires access to a port whichis not available */
#define async_access_required      84 /**< access to selected property should be asynchronous */
#define semaphore_busy             85 /**< a semahpre or mutex could not be obtained within the requested time */
#define non_existent_elem          86 /**< the requested equipment modules does not exist in the registry */
#define non_existent_fec           87 /**< the requested FEC does not exist */
#define non_existent_client        88 /**< (unused) */
#define cannot_lock                89 /**< not able to obtain an access lock or lock resources */
#define not_running                90 /**< the requested resource is not running */

#define not_posted                 91 /**< a dispatch routine did not post a reply */
#define not_accepted               92 /**< the requested operation was refused */
#define operation_timeout          93 /**< the requested operation did not complete within the time limit specified */
#define illegal_protocol           94 /**< the incoming tine protocol is not recognized (reserved, internal) */
#define gpib_error                 95 /**< gpib bus error */
#define rs232_error                96 /**< rs232 bus error */
#define operation_busy             97 /**< the requested operation is busy and needs more time to complete */
#define connection_timeout         98 /**< a synchronous link did not return within the time limit specified */
#define illegal_mode               99 /**< the requested access mode is not valid */
#define not_owner                 100 /**< caller is not the owner of the requested resource */

#define not_defined               101 /**< the requested operation has not yet been defined */
#define net_write_error           102 /**< unable to send data on the network */
#define invalid_data              103 /**< the input data are not valid */
#define software_error            104 /**< general software error */
#define access_denied             105 /**< call was refused by the security subsystem */
#define tcp_not_supported         106 /**< an attempt to use TCP/IP access was made and is not supported on this host (reserved) */
#define ipx_not_supported	        107 /**< an attempt to use IPX access was made and is not supported on this host (reserved) */
#define host_not_resolved	        108 /**< the requested host could not be resolved */
#define tcp_connect_error	        109 /**< unable to connect to a TCP/IP socket (reserved) */
#define tcp_socket_error 	        110 /**< general TCP/IP socket error (reserved) */

#define configuration_error       111 /**< resource configuration error (e.g. redirection information) */
#define invalid_return_code       112 /**< indicates a systematic error in processing the equipment dispatch handler (reserved) */
#define link_blacklisted 	        113 /**< requested link parameters cannot succeed (due to initial return of non_existent_elem, address_unknown, etc.) */
#define link_exists      	        114 /**< requested link already exists in the link table indicating a dependent link (reserved,internal) */
#define max_alarms_exceeded       115 /**< attempt to add an alarm to the local alarm system whose capacity has been reached */
#define not_exported              116 /**< (unused) */
#define is_switching              117 /**< switching transition is in progress (unused) */
#define at_limit                  118 /**< setting is at a maximum or minimum (unused) */
#define get_subscription_id       119 /**< used by access mode CM_NETWORK (reserved, internal) */
#define renegotiate_contract      120 /**< used to repackage output data when few bytes are returned than requested (reserved, internal) */

#define server_redirection        121 /**< used to supply redirection information (reserved, internal) */
#define value_too_high            122 /**< returned value is too high */
#define value_too_low             123 /**< returned value is too low  */
#define warn_too_high             124 /**< returned value is approaching too high  */
#define warn_too_low              125 /**< returned value is approaching too low  */
#define completion_waiting        126 /**< dependent link is waiting for parent to complete */
#define cannot_select             127 /**< select() operation failed (reserved) */
#define has_query_function        128 /**< query call to PROPERTIES or DEVICES has completed but indicated a hierachical preference (reserved) */
#define not_signalled             129 /**< contract is waiting to call the equipment module dispatch routine (reserved, internal) */
#define call_redirection          130 /**< nominal alias for server_redirection */

#define udp_socket_error 	        131 /**< general UDP socket error (reserved) */
#define mutex_error               132 /**< a mutex could not be obtained or released */
#define data_not_local            133 /**< indicates that a device withing a wildcard call needs to span multiple servers (reserved, internal) */
#define name_already_exists       134 /**< a structure or bitfield field name already exists */
#define already_assigned          135 /**< the requested resource has already been assigned */
#define invalid_field             136 /**< a structure or bitfield field is invalid */
#define mcast_init_error          137 /**< could not initialize a multicast socket (reserved) */
#define invalid_mcast_group       138 /**< multicast group address is invalid */
#define non_existent_function     139 /**< an equipment module's dispatch routine is not registered */
#define property_is_mca           140 /**< requested property is a multi-channel array (reserved, internal) */

#define invalid_name              141 /**< input name is invalid */
#define invalid_value             142 /**< input value is invalid */
#define device_error              143 /**< requested device is in an error state */
#define device_not_ready          144 /**< requested device is not ready */
#define device_busy               145 /**< requested device is busy */
#define invalid_interval          146 /**< requested polling interval is invalid */
#define notification_pending      147 /**< a link is pending notification (reserved, internal) */
#define thread_error              148 /**< a thread operation failed */
#define is_terminating            149 /**< the device server is terminating (reserved, internal) */
#define name_server_unknown       150 /**< the equipment name server address is unknown (reserved, internal) */

#define lock_required             151 /**< an access lock is required to process the operation */
#define not_initialized           152 /**< server is not yet initialized */
#define invalid_structure_size    153 /**< the input structure size does not match the size in the structure registry */
#define canbus_error              154 /**< can bus error */
#define profibus_error            155 /**< profibus error */
#define data_stale                156 /**< contract's data are stale (reserved) */
#define has_bitfield_tag          157 /**< reqested contract needs to request a bitfield (reserved, internal) */
#define not_locked                158 /**< the requested resource does not have an access lock (reserved) */
#define lock_expired              159 /**< the requested access lock has expired (reserved) */

#define not_registered            160 /**< the requested property (e.g. SystemScheduleProperty()) has not been registered */
#define non_existent_property     161 /**< the requested property does not exist */
#define memory_overwrite          162 /**< a heap error within a dispatch routine has been detected */
#define server_idle               163 /**< the server is in an idle state */
#define not_applicable            164 /**< the requested opertion is not applicable in the current context */
#define access_locked             165 /**< access temporarily suspended due to access lock */
#define invalid_reference         166 /**< invalid data reference handle */
#define has_structure_tag         167 /**< reqested contract needs to request a structure (reserved, internal) */
#define warn_disk_space           168 /**< warn that disk space is low (local alarm system) */
#define low_disk_space            169 /**< signal that disk space is very low (local alarm system)  */

#define information_static        170 /**< signal that static information is being 'polled'  */
#define reset_mca_property        171 /**< signal to re-establish all mca single element links */
#define record_too_long           172 /**< record length is too long to process within kernel */
#define unix_socket_error         173 /**< an error condition on a unix (pipe) socket */
#define pipe_error                173 /**< an error condition on a unix (pipe) socket */
#define shm_error                 174 /**< an error occured accessing shared memory */
#define tcp_connection_closed     175 /**< tcp peer has closed the connection */
#define service_unavailable       176 /**< requested service is unavailable */
#define device_offline            177 /**< requested device is offline */
#define out_of_sequence           178 /**< data value is out of sequence */

#else /* tine_decorated_constants defined */

/** \def decorated Error Number Definitions */
#define tdc_operation_success          0  /**< successful completion */
#define tdc_illegal_line               1  /**< illegal line number */
#define tdc_illegal_format             2  /**< illegal format or data type */
#define tdc_illegal_arithmetic         3  /**< illegal arithmetic expression */
#define tdc_ambiguous                  4  /**< ambigous input information */
#define tdc_illegal_delimiter          5  /**< illegal delimiter in expression */
#define tdc_zero_divide                6  /**< divide by zero */
#define tdc_work_area_full             7  /**< working area buffer is full */
#define tdc_non_existent               8  /**< requested target does not exist */
#define tdc_invalid_transport          9  /**< transport medium is not valid */
#define tdc_data_size_mismatch         10 /**< link exhausted (unused) */

#define tdc_not_properly               11 /**< expression not properly terminated (unused) */
#define tdc_un_allocated               12 /**< requested target is not allocated in the given context */
#define tdc_no_such_line               13 /**< expected line not found in database */
#define tdc_illegal_data_size          14 /**< give data size is illegal or not valid */
#define tdc_io_error                   15 /**< bus or network io error */
#define tdc_illegal_context            16 /**< given context is illegal */
#define tdc_runtime_error              17 /**< runtime error */
#define tdc_system_error               18 /**< system internal error */
#define tdc_hardware_busy              19 /**< hardware is busy */
#define tdc_argument_list_error        20 /**< input argument or parameter is invalid (alias for invalid_parameter */
#define tdc_invalid_parameter          20 /**< input argument or parameter is invalid */

#define tdc_file_error                 21 /**< file access error */
#define tdc_use_stream_transport       22 /**< datagram message size exceeds the allowable mtu (reserved) */
#define tdc_dimension_error            23 /**< array dimension error */
#define tdc_square_root_negative       24 /**< square root of a negative number */
#define tdc_buffer_too_small           25 /**< supplied buffer is too small to handle the request */
#define tdc_string_too_long            26 /**< supplied string is too long to buffer and cannot be truncated */
#define tdc_ipx_socket_error           27 /**< ipx socket error (reserved) */
#define tdc_net_read_error             28 /**< network read error */
#define tdc_not_ready                  29 /**< call to equipment module dispatch is not ready (internal, reserved) */
#define tdc_invalid_transport_size     30 /**< total message size exceeds the configured allowable message size */

#define tdc_log_negative               31 /**< log of a negative or zero number */
#define tdc_device_not_connected       32 /**< requested device is not connected */
#define tdc_unauthorised_action        33 /**< requested operation is not authorized */
#define tdc_hardware_error             34 /**< hardware error */
#define tdc_illegal_equipment_number   35 /**< requested device does not have a table entry (alias for illegal_device_number) */
#define tdc_illegal_device_number      35 /**< requested device does not have a table entry */
#define tdc_illegal_device             35 /**< requested device is illegal */
#define tdc_illegal_property           36 /**< requested property is illegal */
#define tdc_out_of_range               37 /**< given value is out of range */
#define tdc_not_implemented            38 /**< requested operation has not been implemented */
#define tdc_no_such_computer           39 /**< requested host computer does not exist */
#define tdc_struct_sealed              40 /**< tagged structure has already been sealed (reserved) */

#define tdc_syntax_error               41 /**< expression syntax error (unused) */
#define tdc_no_such_file               42 /**< requested file does not exist */
#define tdc_file_already_exists        43 /**< requested file already exists */
#define tdc_no_file_space              44 /**< not enough capacity in pre-allocated random access file to complete the requested operation */
#define tdc_link_not_open              45 /**< presistent link has timed out (alias for link_timeout) */
#define tdc_link_timeout               45 /**< presistent link has timed out */
#define tdc_remitted_data_lost         46 /**< returned only a portion of the requested data set (reserved) */
#define tdc_end_of_file                47 /**< attempt to read past the end of a file (unused) */
#define tdc_archive_busy               48 /**< archive operation is in progress */
#define tdc_server_name_in_use         49 /**< device server name is already being used */
#define tdc_no_such_column             50 /**< expected column not found in database */

#define tdc_out_of_client_memory       51 /**< a client call could not allocate sufficient memory */
#define tdc_mcast_access_required      52 /**< access requires CA_NETWORK which was not supplied! */
#define tdc_address_unresolved         53 /**< requested target address could not be resolved */
#define tdc_invalid_property           54 /**< requested property is not properly formed */
#define tdc_address_unknown            55 /**< requested target address could not be found */
#define tdc_name_unknown               56 /**< requested name could not be found */
#define tdc_invalid_keyword            57 /**< input globals keyword is invalid */
#define tdc_invalid_link               58 /**< input link id is invalid */
#define tdc_string_expected            59 /**< database file read encountered an end-of-line while expecting a string */       
#define tdc_out_of_local_memory        60 /**< a local call to allocated memory was not successful */

#define tdc_out_of_shared_memory       61 /**< a call to obtain a shared-memory buffer was not successful */
#define tdc_invalid_structure_tag      62 /**< the given structure tag is not in the structure registry */
#define tdc_invalid_index              63 /**< the index given does not reference a table entry */
#define tdc_illegal_equipment_name     64 /**< the equipment name given is not valid */
#define tdc_link_error                 65 /**< internal error processing data link (reserved) */
#define tdc_code_failure               66 /**< internal failure due to unhandled exceptions, unavailable resource, OS failure, etc. */
#define tdc_non_existent_server        67 /**< requested server does not exist */
#define tdc_function_deprecated        68 /**< requested action or functionality has been deprecated */
#define tdc_not_supported              69 /**< requested action or functionality not supported */
#define tdc_illegal_character          70 /**< input character is illegal (unused) */
#define tdc_parsing_error              70 /**< error parsing input */

#define tdc_illegal_operator           71 /**< input arithmetic or logical operator is not allowed */
#define tdc_not_allowed                72 /**< requested operation is not allowed */
#define tdc_illegal_read_write         73 /**< requested access is not allowed (typically CA_WRITE) */
#define tdc_out_of_server_memory       74 /**< a server call could not allocate sufficient memory */
#define tdc_database_not_loaded        75 /**< requested operation requires a database which has not been loaded */
#define tdc_illegal_command            76 /**< requested command is illegal */
#define tdc_resources_exhausted        77 /**< requested operation would exceed the configured capacity of an allocation table (clients, contracts, connections, etc.) */
#define tdc_file_not_open              78 /**< requested file io could not be made as the file is not open (unused) */
#define tdc_sedac_error                79 /**< sedac bus error */
#define tdc_reacquire_address          80 /**< signal to reacquire a server's address from ENS (reserved, internal) */

#define tdc_semaphore_error            81 /**< a semaphore or mutex could not be obtained or released */
#define tdc_driver_not_installed       82 /**< request requires a hardware driver which is not installed */
#define tdc_port_not_available         83 /**< request requires access to a port whichis not available */
#define tdc_scan_error                 84 /**< (unused) */
#define tdc_semaphore_busy             85 /**< a semahpre or mutex could not be obtained within the requested time */
#define tdc_non_existent_elem          86 /**< the requested equipment modules does not exist in the registry */
#define tdc_non_existent_fec           87 /**< the requested FEC does not exist */
#define tdc_non_existent_client        88 /**< (unused) */
#define tdc_cannot_lock                89 /**< not able to obtain an access lock or lock resources */
#define tdc_not_running                90 /**< the requested resource is not running */

#define tdc_not_posted                 91 /**< a dispatch routine did not post a reply */
#define tdc_not_accepted               92 /**< the requested operation was refused */
#define tdc_operation_timeout          93 /**< the requested operation did not complete within the time limit specified */
#define tdc_illegal_protocol           94 /**< the incoming tine protocol is not recognized (reserved, internal) */
#define tdc_gpib_error                 95 /**< gpib bus error */
#define tdc_rs232_error                96 /**< rs232 bus error */
#define tdc_operation_busy             97 /**< the requested operation is busy and needs more time to complete */
#define tdc_connection_timeout         98 /**< a synchronous link did not return within the time limit specified */
#define tdc_illegal_mode               99 /**< the requested access mode is not valid */
#define tdc_not_owner                 100 /**< (unused) */

#define tdc_not_defined               101 /**< the requested operation has not yet been defined */
#define tdc_net_write_error           102 /**< unable to send data on the network */
#define tdc_invalid_data              103 /**< the input data are not valid */
#define tdc_software_error            104 /**< general software error */
#define tdc_access_denied             105 /**< call was refused by the security subsystem */
#define tdc_tcp_not_supported         106 /**< an attempt to use TCP/IP access was made and is not supported on this host (reserved) */
#define tdc_ipx_not_supported	        107 /**< an attempt to use IPX access was made and is not supported on this host (reserved) */
#define tdc_host_not_resolved	        108 /**< the requested host could not be resolved */
#define tdc_tcp_connect_error	        109 /**< unable to connect to a TCP/IP socket (reserved) */
#define tdc_tcp_socket_error 	        110 /**< general TCP/IP socket error (reserved) */

#define tdc_configuration_error       111 /**< resource configuration error (e.g. redirection information) */
#define tdc_invalid_return_code       112 /**< indicates a systematic error in processing the equipment dispatch handler (reserved) */
#define tdc_link_blacklisted 	        113 /**< requested link parameters cannot succeed (due to initial return of non_existent_elem, address_unknown, etc.) */
#define tdc_link_exists      	        114 /**< requested link already exists in the link table indicating a dependent link (reserved,internal) */
#define tdc_max_alarms_exceeded       115 /**< attempt to add an alarm to the local alarm system whose capacity has been reached */
#define tdc_not_exported              116 /**< (unused) */
#define tdc_is_switching              117 /**< switching transition is in progress (unused) */
#define tdc_at_limit                  118 /**< setting is at a maximum or minimum (unused) */
#define tdc_get_subscription_id       119 /**< used by access mode CM_NETWORK (reserved, internal) */
#define tdc_renegotiate_contract      120 /**< used to repackage output data when few bytes are returned than requested (reserved, internal) */

#define tdc_server_redirection        121 /**< used to supply redirection information (reserved, internal) */
#define tdc_value_too_high            122 /**< returned value is too high */
#define tdc_value_too_low             123 /**< returned value is too low  */
#define tdc_warn_too_high             124 /**< returned value is approaching too high  */
#define tdc_warn_too_low              125 /**< returned value is approaching too low  */
#define tdc_completion_waiting        126 /**< dependent link is waiting for parent to complete */
#define tdc_cannot_select             127 /**< select() operation failed (reserved) */
#define tdc_has_query_function        128 /**< query call to PROPERTIES or DEVICES has completed but indicated a hierachical preference (reserved) */
#define tdc_not_signalled             129 /**< contract is waiting to call the equipment module dispatch routine (reserved, internal) */
#define tdc_call_redirection          130 /**< nominal alias for server_redirection */

#define tdc_udp_socket_error 	        131 /**< general UDP socket error (reserved) */
#define tdc_mutex_error               132 /**< a mutex could not be obtained or released */
#define tdc_data_not_local            133 /**< indicates that a device withing a wildcard call needs to span multiple servers (reserved, internal) */
#define tdc_name_already_exists       134 /**< a structure or bitfield field name already exists */
#define tdc_already_assigned          135 /**< the requested resource has already been assigned */
#define tdc_invalid_field             136 /**< a structure or bitfield field is invalid */
#define tdc_mcast_init_error          137 /**< could not initialize a multicast socket (reserved) */
#define tdc_invalid_mcast_group       138 /**< multicast group address is invalid */
#define tdc_non_existent_function     139 /**< an equipment module's dispatch routine is not registered */
#define tdc_property_is_mca           140 /**< requested property is a multi-channel array (reserved, internal) */

#define tdc_invalid_name              141 /**< input name is invalid */
#define tdc_invalid_value             142 /**< input value is invalid */
#define tdc_device_error              143 /**< requested device is in an error state */
#define tdc_device_not_ready          144 /**< requested device is not ready */
#define tdc_device_busy               145 /**< requested device is busy */
#define tdc_invalid_interval          146 /**< requested polling interval is invalid */
#define tdc_notification_pending      147 /**< a link is pending notification (reserved, internal) */
#define tdc_thread_error              148 /**< a thread operation failed */
#define tdc_is_terminating            149 /**< the device server is terminating (reserved, internal) */
#define tdc_name_server_unknown       150 /**< the equipment name server address is unknown (reserved, internal) */

#define tdc_lock_required             151 /**< an access lock is required to process the operation */
#define tdc_not_initialized           152 /**< server is not yet initialized */
#define tdc_invalid_structure_size    153 /**< the input structure size does not match the size in the structure registry */
#define tdc_canbus_error              154 /**< can bus error */
#define tdc_profibus_error            155 /**< profibus error */
#define tdc_data_stale                156 /**< contract's data are stale (reserved) */
#define tdc_has_bitfield_tag          157 /**< reqested contract needs to request a bitfield (reserved, internal) */
#define tdc_not_locked                158 /**< the requested resource does not have an access lock (reserved) */
#define tdc_lock_expired              159 /**< the requested access lock has expired (reserved) */

#define tdc_not_registered            160 /**< the requested property (e.g. SystemScheduleProperty()) has not been registered */
#define tdc_non_existent_property     161 /**< the requested property does not exist */
#define tdc_memory_overwrite          162 /**< a heap error within a dispatch routine has been detected */
#define tdc_server_idle               163 /**< the server is in an idle state */
#define tdc_not_applicable            164 /**< the requested opertion is not applicable in the current context */
#define tdc_access_locked             165 /**< access temporarily suspended due to exclusive read */
#define tdc_invalid_reference         166 /**< invalid data reference handle */
#define tdc_has_structure_tag         167 /**< reqested contract needs to request a structure (reserved, internal) */
#define tdc_warn_disk_space           168 /**< warn that disk space is low (local alarm system) */
#define tdc_low_disk_space            169 /**< signal that disk space is very low (local alarm system)  */

#define tdc_information_static        170 /**< signal that static information is being 'polled'  */
#define tdc_reset_mca_property        171 /**< signal to re-establish all mca single element links */
#define tdc_record_too_long           172 /**< record length is too long to process within kernel */
#define tdc_unix_socket_error         173 /**< an error condition on a unix (pipe) socket */
#define tdc_pipe_error                173 /**< an error condition on a unix (pipe) socket */
#define tdc_shm_error                 174 /**< an error occured accessing shared memory */
#define tdc_tcp_connection_closed     175 /**< tcp peer has closed the connection */
#define tdc_service_unavailable       176 /**< requested service is unavailable */
#define tdc_device_offline            177 /**< requested device is offline */
#define tdc_out_of_sequence           178 /**< data value is out of sequence */

#endif /* !defined tine_decorated_constants */

extern int numErr;

#if defined(__cplusplus) && !defined(FORCE_CPP)
}
#endif

#endif /* ERRORS_H_ */
