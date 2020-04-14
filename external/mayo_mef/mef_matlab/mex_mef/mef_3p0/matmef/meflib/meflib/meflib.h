
#ifndef MEFLIB_IN
#define MEFLIB_IN


/************************************************************************************/
/*******************************  MEF 3.0 C Library Header **************************/
/************************************************************************************/


// Specification for Multiscale Electrophysiology Format (MEF) version 3.0
// Copyright 2013, Mayo Foundation, Rochester MN. All rights reserved.
// Written by Matt Stead, Ben Brinkmann, and Dan Crepeau.

// Usage and modification of this source code is governed by the Apache 2.0 license.
// You may not use this file except in compliance with this License.
// A copy of the Apache 2.0 License may be obtained at http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under this License is distributed on an "as is" basis,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either expressed or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Thanks to all who acknowledge the Mayo Systems Electrophysiology Laboratory, Rochester, MN
// in academic publications of their work facilitated by this software.

// The encryption / decryption algorithm is the 128-bit AES standard ( http://www.csrc.nist.gov/publications/fips/fips197/fips-197.pdf ).
// AES routines (128 bit only) are included in the library, with attribution, for convenience.

// The hash algorithm is the SHA-256 standard ( http://csrc.nist.gov/publications/fips/fips180-4/fips-180-4.pdf ).
// Basic SHA-256 routines are included in the library, with attribution, for convenience.

// Strings in MEF 3.0 are encoded in the Universal Character Set standard, ISO/IEC 10646:2012 otherwise known as UTF-8.
// ( http://standards.iso.org/ittf/PubliclyAvailableStandards/c056921_ISO_IEC_10646_2012.zip )
// Basic UTF-8 manipulation routines are also included in the library, with attribution, for convenience.

// written with tab width = indent width = 8 spaces and a monospaced font
// set editor prefernces to these for intended alignment




/************************************************************************************/
/**********************************  Library Includes  ******************************/
/************************************************************************************/
#ifdef _WIN32
	#define _USE_MATH_DEFINES

	#include <stdlib.h>
	#include <stdio.h>
	#include <string.h>
	#include <time.h>
	#include <math.h>
	#include <float.h>
	#include <sys/stat.h>
	#include <stdarg.h>
	#include <errno.h>
	#include <fcntl.h>
	#include <limits.h>
	#include <malloc.h>  // for alloca()
	#include <stdint.h>
	#include <windows.h>

#else
	#include <stdlib.h>
	#include <unistd.h>
	#include <stdio.h>
	#include <string.h>
	#include <sys/time.h>
	#include <math.h>
	#include <float.h>
	#include <sys/stat.h>
	#include <stdarg.h>
	#include <errno.h>
	#include <fcntl.h>
	#include <limits.h>
	#include <dirent.h>
	//#include <pthread.h>
#endif





/************************************************************************************/
/********************************  Elemental Typedefs  ******************************/
/************************************************************************************/

#ifndef SIZE_TYPES_IN
#define SIZE_TYPES_IN

#ifdef _WIN32
	typedef char			si1;
	typedef unsigned char		ui1;
	typedef short			si2;
	typedef unsigned short		ui2;
	typedef int			si4;
	typedef unsigned int		ui4;
	typedef long long int		si8;
	typedef long long unsigned int	ui8;
	typedef float			sf4;
	typedef double			sf8;
	typedef long double		sf16;   // NOTE: it often requires an explicit compiler instruction to implement true long floating point math
						// in gcc and icc: "-Qoption,cpp,--extended_float_type"
	typedef uint8_t 		u_int8_t;
	typedef uint16_t 		u_int16_t;
	typedef uint32_t 		u_int32_t;
#else
	typedef char			si1;
	typedef unsigned char		ui1;
	typedef short			si2;
	typedef unsigned short		ui2;
	typedef int			si4;
	typedef unsigned int		ui4;
	typedef long int		si8;
	typedef long unsigned int	ui8;
	typedef float			sf4;
	typedef double			sf8;
	typedef long double		sf16;   // NOTE: it often requires an explicit compiler instruction to implement true long floating point math
						// in gcc and icc: "-Qoption,cpp,--extended_float_type"

#endif
#endif



/************************************************************************************/
/************************************  MEF Globals  *********************************/
/************************************************************************************/

typedef struct {
        // time constants
	si8	recording_time_offset;
	ui4	recording_time_offset_mode;
	si4	GMT_offset;
        si8	DST_start_time;
        si8	DST_end_time;
	// alignment fields
	si4	universal_header_aligned;
	si4	metadata_section_1_aligned;
	si4	time_series_metadata_section_2_aligned;
	si4	video_metadata_section_2_aligned;
	si4	metadata_section_3_aligned;
	si4	all_metadata_structures_aligned;
	si4	time_series_indices_aligned;
	si4	video_indices_aligned;
	si4	RED_block_header_aligned;
	si4	record_header_aligned;
	si4	record_indices_aligned;
	si4	all_record_structures_aligned;
	si4	all_structures_aligned;
	// RED
	sf8	*RED_normal_CDF_table;
	// CRC
	ui4	*CRC_table;
        ui4	CRC_mode;
	// AES tables
	si4	*AES_sbox_table;
	si4	*AES_rcon_table;
	si4	*AES_rsbox_table;
	// SHA256 tables
	ui4	*SHA256_h0_table;
	ui4	*SHA256_k_table;
	// UTF8 tables
	ui4	*UTF8_offsets_from_UTF8_table;
	si1	*UTF8_trailing_bytes_for_UTF8_table;
        // miscellaneous
        si4	verbose;
        ui4	behavior_on_fail;
        ui4	file_creation_umask;
} MEF_GLOBALS;



/************************************************************************************/
/************************  Error Checking Standard Functions  ***********************/
/************************************************************************************/

// Constants
#define USE_GLOBAL_BEHAVIOR	0
#define RESTORE_BEHAVIOR	1
#define EXIT_ON_FAIL		2
#define RETURN_ON_FAIL		4
#define SUPPRESS_ERROR_OUTPUT	8

// Function Prototypes
si4 	e_system(const char *command, const si1 *function, si4 line, ui4 behavior_on_fail);
void	*e_calloc(size_t n_members, size_t size, const si1 *function, si4 line, ui4 behavior_on_fail);
FILE	*e_fopen(si1 *path, si1 *mode, const si1 *function, si4 line, ui4 behavior_on_fail);
size_t	e_fread(void *ptr, size_t size, size_t n_members, FILE *stream, si1 *path, const si1 *function, si4 line, ui4 behavior_on_fail);
si4	e_fseek(FILE *stream, size_t offset, si4 whence, si1 *path, const si1 *function, si4 line, ui4 behavior_on_fail);
long	e_ftell(FILE *stream, const si1 *function, si4 line, ui4 behavior_on_fail);
size_t	e_fwrite(void *ptr, size_t size, size_t n_members, FILE *stream, si1 *path, const si1 *function, si4 line, ui4 behavior_on_fail);
void	*e_malloc(size_t n_bytes, const si1 *function, si4 line, ui4 behavior_on_fail);
void	*e_realloc(void *ptr, size_t n_bytes, const si1 *function, si4 line, ui4 behavior_on_fail);



/************************************************************************************/
/**********************************  MEF Constants  *********************************/
/************************************************************************************/

// Miscellaneous Constants
#define MEF_TRUE				1
#define MEF_UNKNOWN				0
#define MEF_FALSE				-1
#define MEF_BIG_ENDIAN				0
#define MEF_LITTLE_ENDIAN			1
#define MEF_VERSION_MAJOR			3
#define MEF_VERSION_MINOR			0
#define TYPE_BYTES				5
#define UUID_BYTES				16
#define NO_UUID					0
#define TIME_STRING_BYTES			32
#define MEF_BASE_FILE_NAME_BYTES		256	// utf8[63]
#define MEF_SEGMENT_BASE_FILE_NAME_BYTES	(MEF_BASE_FILE_NAME_BYTES + 8)
#define MEF_FULL_FILE_NAME_BYTES		1024	// utf8[255]
#define PAD_BYTE_VALUE				0x7e	// ascii tilde ("~")
#define FILE_NUMBERING_DIGITS			6
#define NO_TYPE_CODE				0
#define MAXIMUM_GMT_OFFSET			86400
#define MINIMUM_GMT_OFFSET			-86400
#define UNKNOWN_NUMBER_OF_ENTRIES		-1
#define UUTC_NO_ENTRY				0x8000000000000000
#define CRC_NO_ENTRY				0

// CRC Modes
#define CRC_IGNORE				0
#define CRC_VALIDATE				1
#define CRC_VALIDATE_ON_INPUT			2
#define CRC_VALIDATE_ON_OUTPUT			4
#define CRC_CALCULATE				8
#define CRC_CALCULATE_ON_INPUT			16
#define CRC_CALCULATE_ON_OUTPUT			32

// Encryption & Password Constants
#define ENCRYPTION_LEVEL_NO_ENTRY		-128
#define NO_ENCRYPTION				0
#define LEVEL_0_ENCRYPTION			NO_ENCRYPTION
#define LEVEL_1_ENCRYPTION			1
#define LEVEL_2_ENCRYPTION			2
#define LEVEL_0_ACCESS				LEVEL_0_ENCRYPTION
#define LEVEL_1_ACCESS				LEVEL_1_ENCRYPTION
#define LEVEL_2_ACCESS				LEVEL_2_ENCRYPTION
#define LEVEL_1_ENCRYPTION_DECRYPTED		-LEVEL_1_ENCRYPTION
#define LEVEL_2_ENCRYPTION_DECRYPTED		-LEVEL_2_ENCRYPTION
#define ENCRYPTION_BLOCK_BYTES			16		// AES-128
#define ENCRYPTION_KEY_BYTES			176		// AES-128   = ((AES_NR + 1) * AES_NK * AES_NB)
#define PASSWORD_BYTES				ENCRYPTION_BLOCK_BYTES
#define UTF8_PASSWORD_BYTES			(PASSWORD_BYTES * 4)
#define MAX_PASSWORD_CHARACTERS			(PASSWORD_BYTES - 1)
#define PASSWORD_VALIDATION_FIELD_BYTES		PASSWORD_BYTES

// Recording Time Offset Modes & Constants
#define RTO_USE_SYSTEM_TIME			-1
#define RTO_IGNORE				0
#define RTO_APPLY				1
#define RTO_REMOVE				2
#define RTO_APPLY_ON_OUTPUT			4
#define RTO_APPLY_ON_INPUT			8
#define RTO_REMOVE_ON_OUTPUT			16
#define RTO_REMOVE_ON_INPUT			32
#define RTO_INPUT_ACTION			1
#define RTO_OUTPUT_ACTION			2

// Global Defaults
#define MEF_GLOBALS_VERBOSE_DEFAULT			MEF_FALSE
#define MEF_GLOBALS_RECORDING_TIME_OFFSET_DEFAULT	0
#define MEF_GLOBALS_RECORDING_TIME_OFFSET_MODE_DEFAULT	(RTO_APPLY_ON_OUTPUT | RTO_REMOVE_ON_INPUT)
#define MEF_GLOBALS_GMT_OFFSET_DEFAULT			0
#define MEF_GLOBALS_DST_START_TIME_DEFAULT		UUTC_NO_ENTRY
#define MEF_GLOBALS_DST_END_TIME_DEFAULT		UUTC_NO_ENTRY
#define MEF_GLOBALS_FILE_CREATION_UMASK_DEFAULT		S_IWOTH  // defined in <sys/stat.h>
#define MEF_GLOBALS_BEHAVIOR_ON_FAIL_DEFAULT		EXIT_ON_FAIL
#define MEF_GLOBALS_CRC_MODE_DEFAULT			(CRC_CALCULATE_ON_OUTPUT)

// File Type Constants
#define NO_FILE_TYPE_STRING				""				// ascii[4]
#define NO_FILE_TYPE_CODE				(ui4) 0x00000000		// ui4 (big & little endian)
#define SESSION_DIRECTORY_TYPE_STRING			"mefd"				// ascii[4]
#define SESSION_DIRECTORY_TYPE_CODE			(ui4) 0x6466656d		// ui4 (little endian)
// #define SESSION_DIRECTORY_TYPE_CODE			(ui4) 0x6d656664		// ui4 (big endian)
#define SEGMENT_DIRECTORY_TYPE_STRING			"segd"				// ascii[4]
#define SEGMENT_DIRECTORY_TYPE_CODE			(ui4) 0x64676573		// ui4 (little endian)
// #define SEGMENT_DIRECTORY_TYPE_CODE			(ui4) 0x73656764		// ui4 (big endian)
#define RECORD_DATA_FILE_TYPE_STRING			"rdat"				// ascii[4]
#define RECORD_DATA_FILE_TYPE_CODE			(ui4) 0x74616472		// ui4 (little endian)
// #define RECORD_DATA_FILE_TYPE_CODE			(ui4) 0x72646174		// ui4 (big endian)
#define RECORD_INDICES_FILE_TYPE_STRING			"ridx"				// ascii[4]
#define RECORD_INDICES_FILE_TYPE_CODE			(ui4) 0x78646972		// ui4 (little endian)
// #define RECORD_INDICES_FILE_TYPE_CODE		(ui4) 0x72696478		// ui4 (big endian)
#define VIDEO_CHANNEL_DIRECTORY_TYPE_STRING		"vidd"				// ascii[4]
#define VIDEO_CHANNEL_DIRECTORY_TYPE_CODE		(ui4) 0x64646976		// ui4 (little endian)
// #define VIDEO_CHANNEL_DIRECTORY_TYPE_CODE		(ui4) 0x76696464		// ui4 (big endian)
#define VIDEO_METADATA_FILE_TYPE_STRING			"vmet"				// ascii[4]
#define VIDEO_METADATA_FILE_TYPE_CODE			(ui4) 0x74656d76		// ui4 (little endian)
// #define VIDEO_METADATA_FILE_TYPE_CODE		(ui4) 0x766d6574		// ui4 (big endian)
#define VIDEO_INDICES_FILE_TYPE_STRING			"vidx"				// ascii[4]
#define VIDEO_INDICES_FILE_TYPE_CODE			(ui4) 0x78646976		// ui4 (little endian)
// #define VIDEO_INDICES_FILE_TYPE_CODE			(ui4) 0x76696478		// ui4 (big endian)
#define TIME_SERIES_CHANNEL_DIRECTORY_TYPE_STRING	"timd"				// ascii[4]
#define TIME_SERIES_CHANNEL_DIRECTORY_TYPE_CODE		(ui4) 0x646d6974		// ui4 (little endian)
// #define TIME_SERIES_CHANNEL_DIRECTORY_TYPE_CODE	(ui4) 0x74696d64		// ui4 (big endian)
#define TIME_SERIES_METADATA_FILE_TYPE_STRING		"tmet"				// ascii[4]
#define TIME_SERIES_METADATA_FILE_TYPE_CODE		(ui4) 0x74656d74		// ui4 (little endian)
// #define TIME_SERIES_METADATA_FILE_TYPE_CODE		(ui4) 0x746d6574		// ui4 (big endian)
#define TIME_SERIES_DATA_FILE_TYPE_STRING		"tdat"				// ascii[4]
#define TIME_SERIES_DATA_FILE_TYPE_CODE			(ui4) 0x74616474		// ui4 (little endian)
// #define TIME_SERIES_DATA_FILE_TYPE_CODE		(ui4) 0x74646174		// ui4 (big endian)
#define TIME_SERIES_INDICES_FILE_TYPE_STRING		"tidx"				// ascii[4]
#define TIME_SERIES_INDICES_FILE_TYPE_CODE		(ui4) 0x78646974		// ui4 (little endian)
// #define TIME_SERIES_INDICES_FILE_TYPE_CODE		(ui4) 0x74696478		// ui4 (big endian)

// Channel Types
#define UNKNOWN_CHANNEL_TYPE		-1
#define TIME_SERIES_CHANNEL_TYPE	1
#define VIDEO_CHANNEL_TYPE		2

// File Processing constants
#define FPS_FILE_LENGTH_UNKNOWN					-1
#define FPS_FULL_FILE						-1
#define FPS_NO_LOCK_TYPE					~(F_RDLCK | F_WRLCK | F_UNLCK)  // from <fcntl.h>
#define FPS_NO_LOCK_MODE					0
#define FPS_READ_LOCK_ON_READ_OPEN				1
#define FPS_WRITE_LOCK_ON_READ_OPEN				2
#define FPS_WRITE_LOCK_ON_WRITE_OPEN				4
#define FPS_WRITE_LOCK_ON_READ_WRITE_OPEN			8
#define FPS_READ_LOCK_ON_READ					16
#define FPS_WRITE_LOCK_ON_WRITE					32
#define FPS_NO_OPEN_MODE					0
#define FPS_R_OPEN_MODE						1
#define FPS_R_PLUS_OPEN_MODE					2
#define FPS_W_OPEN_MODE						4
#define FPS_W_PLUS_OPEN_MODE					8
#define FPS_A_OPEN_MODE						16
#define FPS_A_PLUS_OPEN_MODE					32
#define FPS_GENERIC_READ_OPEN_MODE				(FPS_R_OPEN_MODE | FPS_R_PLUS_OPEN_MODE | FPS_W_PLUS_OPEN_MODE | FPS_A_PLUS_OPEN_MODE)
#define FPS_GENERIC_WRITE_OPEN_MODE				(FPS_R_PLUS_OPEN_MODE | FPS_W_OPEN_MODE | FPS_W_PLUS_OPEN_MODE | FPS_A_OPEN_MODE | FPS_A_PLUS_OPEN_MODE)

// File Processing Directives defaults
#define FPS_DIRECTIVE_CLOSE_FILE_DEFAULT			MEF_TRUE
#define FPS_DIRECTIVE_FREE_PASSWORD_DATA_DEFAULT		MEF_FALSE
#define FPS_DIRECTIVE_LOCK_MODE_DEFAULT				(FPS_READ_LOCK_ON_READ_OPEN | FPS_WRITE_LOCK_ON_WRITE_OPEN | FPS_WRITE_LOCK_ON_READ_WRITE_OPEN)
#define FPS_DIRECTIVE_OPEN_MODE_DEFAULT				FPS_NO_OPEN_MODE
#define FPS_DIRECTIVE_IO_BYTES_DEFAULT				FPS_FULL_FILE  // bytes to read or write

// Universal Header: File Format Constants
#define UNIVERSAL_HEADER_OFFSET						0
#define UNIVERSAL_HEADER_BYTES						1024			// 1 kb
#define UNIVERSAL_HEADER_HEADER_CRC_OFFSET				0			// ui4
#define UNIVERSAL_HEADER_HEADER_CRC_NO_ENTRY				CRC_NO_ENTRY
#define UNIVERSAL_HEADER_BODY_CRC_OFFSET				4			// ui4
#define UNIVERSAL_HEADER_BODY_CRC_NO_ENTRY				CRC_NO_ENTRY
#define UNIVERSAL_HEADER_FILE_TYPE_OFFSET				8			// ascii[4] or ui4
#define UNIVERSAL_HEADER_FILE_TYPE_NO_ENTRY				0			// zero as ui4 or zero-length string as ascii[4]
#define UNIVERSAL_HEADER_MEF_VERSION_MAJOR_OFFSET			13			// ui1
#define UNIVERSAL_HEADER_MEF_VERSION_MAJOR_NO_ENTRY			0xff
#define UNIVERSAL_HEADER_MEF_VERSION_MINOR_OFFSET			14			// ui1
#define UNIVERSAL_HEADER_MEF_VERSION_MINOR_NO_ENTRY			0xff
#define UNIVERSAL_HEADER_BYTE_ORDER_CODE_OFFSET				15			// ui1
#define UNIVERSAL_HEADER_BYTE_ORDER_CODE_NO_ENTRY			0xff
#define UNIVERSAL_HEADER_START_TIME_OFFSET				16			// si8
#define UNIVERSAL_HEADER_START_TIME_NO_ENTRY				UUTC_NO_ENTRY
#define UNIVERSAL_HEADER_END_TIME_OFFSET				24			// si8
#define UNIVERSAL_HEADER_END_TIME_NO_ENTRY				UUTC_NO_ENTRY
#define UNIVERSAL_HEADER_NUMBER_OF_ENTRIES_OFFSET			32			// si8
#define UNIVERSAL_HEADER_NUMBER_OF_ENTRIES_NO_ENTRY			-1
#define UNIVERSAL_HEADER_MAXIMUM_ENTRY_SIZE_OFFSET			40			// si8
#define UNIVERSAL_HEADER_MAXIMUM_ENTRY_SIZE_NO_ENTRY			-1
#define UNIVERSAL_HEADER_SEGMENT_NUMBER_OFFSET				48			// si4
#define UNIVERSAL_HEADER_SEGMENT_NUMBER_NO_ENTRY			-1
#define UNIVERSAL_HEADER_CHANNEL_LEVEL_CODE				-2
#define UNIVERSAL_HEADER_SESSION_LEVEL_CODE				-3
#define UNIVERSAL_HEADER_CHANNEL_NAME_OFFSET				52			// utf8[63]
#define UNIVERSAL_HEADER_SESSION_NAME_OFFSET				308			// utf8[63]
#define UNIVERSAL_HEADER_ANONYMIZED_NAME_OFFSET				564			// utf8[63]
#define UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES				256			// utf8[63]
#define UNIVERSAL_HEADER_LEVEL_UUID_OFFSET				820			// ui1
#define UNIVERSAL_HEADER_FILE_UUID_OFFSET				836			// ui1
#define UNIVERSAL_HEADER_PROVENANCE_UUID_OFFSET				852			// ui1
#define UNIVERSAL_HEADER_LEVEL_1_PASSWORD_VALIDATION_FIELD_OFFSET	868			// ui1
#define UNIVERSAL_HEADER_LEVEL_2_PASSWORD_VALIDATION_FIELD_OFFSET	884			// ui1
#define UNIVERSAL_HEADER_PROTECTED_REGION_OFFSET			900
#define UNIVERSAL_HEADER_PROTECTED_REGION_BYTES				60
#define UNIVERSAL_HEADER_DISCRETIONARY_REGION_OFFSET			960
#define UNIVERSAL_HEADER_DISCRETIONARY_REGION_BYTES			64

// Metadata: File Format Constants
#define METADATA_FILE_BYTES					16384	// 16 kb

// Metadata: File Format Constants - Section 1 Fields
#define METADATA_SECTION_1_BYTES				1536
#define METADATA_SECTION_2_ENCRYPTION_OFFSET			1024		// si1
#define METADATA_SECTION_2_ENCRYPTION_DEFAULT			LEVEL_1_ENCRYPTION
#define METADATA_SECTION_3_ENCRYPTION_OFFSET			1025		// si1
#define METADATA_SECTION_3_ENCRYPTION_DEFAULT			LEVEL_2_ENCRYPTION
#define METADATA_SECTION_1_PROTECTED_REGION_OFFSET		1026
#define METADATA_SECTION_1_PROTECTED_REGION_BYTES		766
#define METADATA_SECTION_1_DISCRETIONARY_REGION_OFFSET		1792
#define METADATA_SECTION_1_DISCRETIONARY_REGION_BYTES		768

// Metadata: File Format Constants - Section 2 Fields
#define METADATA_SECTION_2_OFFSET				2560
#define METADATA_SECTION_2_BYTES				10752

// Metadata: File Format Constants - Generic Section 2 Fields
#define METADATA_CHANNEL_DESCRIPTION_OFFSET			2560		// utf8[511]
#define METADATA_CHANNEL_DESCRIPTION_BYTES			2048
#define METADATA_SESSION_DESCRIPTION_OFFSET			4608		// utf8[511]
#define METADATA_SESSION_DESCRIPTION_BYTES			2048
#define METADATA_RECORDING_DURATION_OFFSET			6656		// si8
#define METADATA_RECORDING_DURATION_NO_ENTRY			-1

// Metadata: File Format Constants - Time Series Section 2 Fields
#define TIME_SERIES_METADATA_REFERENCE_DESCRIPTION_OFFSET		6664		// utf8[511]
#define TIME_SERIES_METADATA_REFERENCE_DESCRIPTION_BYTES		METADATA_CHANNEL_DESCRIPTION_BYTES
#define TIME_SERIES_METADATA_ACQUISITION_CHANNEL_NUMBER_OFFSET		8712		// si8
#define TIME_SERIES_METADATA_ACQUISITION_CHANNEL_NUMBER_NO_ENTRY	-1
#define TIME_SERIES_METADATA_SAMPLING_FREQUENCY_OFFSET			8720		// sf8
#define TIME_SERIES_METADATA_SAMPLING_FREQUENCY_NO_ENTRY		-1.0
#define TIME_SERIES_METADATA_LOW_FREQUENCY_FILTER_SETTING_OFFSET	8728		// sf8
#define TIME_SERIES_METADATA_LOW_FREQUENCY_FILTER_SETTING_NO_ENTRY	-1.0
#define TIME_SERIES_METADATA_HIGH_FREQUENCY_FILTER_SETTING_OFFSET	8736		// sf8
#define TIME_SERIES_METADATA_HIGH_FREQUENCY_FILTER_SETTING_NO_ENTRY	-1.0
#define TIME_SERIES_METADATA_NOTCH_FILTER_FREQUENCY_SETTING_OFFSET	8744		// sf8
#define TIME_SERIES_METADATA_NOTCH_FILTER_FREQUENCY_SETTING_NO_ENTRY	-1.0
#define TIME_SERIES_METADATA_AC_LINE_FREQUENCY_OFFSET			8752		// sf8
#define TIME_SERIES_METADATA_AC_LINE_FREQUENCY_NO_ENTRY			-1.0
#define TIME_SERIES_METADATA_UNITS_CONVERSION_FACTOR_OFFSET		8760		// sf8
#define TIME_SERIES_METADATA_UNITS_CONVERSION_FACTOR_NO_ENTRY		0.0
#define TIME_SERIES_METADATA_UNITS_DESCRIPTION_OFFSET			8768		// utf8[31]
#define TIME_SERIES_METADATA_UNITS_DESCRIPTION_BYTES			128		// utf8[31]
#define TIME_SERIES_METADATA_MAXIMUM_NATIVE_SAMPLE_VALUE_OFFSET		8896		// ssf8
#define TIME_SERIES_METADATA_MAXIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY	NAN // (nan(NULL))	// NOTE this value must be tested with isnan(), or an equivalent function, rather than ==
#define TIME_SERIES_METADATA_MINIMUM_NATIVE_SAMPLE_VALUE_OFFSET		8904		// sf8
#define TIME_SERIES_METADATA_MINIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY	NAN //(nan(NULL))	// NOTE this value must be tested with isnan(), or an equivalent function, rather than ==
#define TIME_SERIES_METADATA_START_SAMPLE_OFFSET			8912
#define TIME_SERIES_METADATA_START_SAMPLE_NO_ENTRY			-1
#define TIME_SERIES_METADATA_NUMBER_OF_SAMPLES_OFFSET			8920		// si8
#define TIME_SERIES_METADATA_NUMBER_OF_SAMPLES_NO_ENTRY			-1
#define TIME_SERIES_METADATA_NUMBER_OF_BLOCKS_OFFSET			8928		// si8
#define TIME_SERIES_METADATA_NUMBER_OF_BLOCKS_NO_ENTRY			-1
#define TIME_SERIES_METADATA_MAXIMUM_BLOCK_BYTES_OFFSET			8936		// si8
#define TIME_SERIES_METADATA_MAXIMUM_BLOCK_BYTES_NO_ENTRY		-1
#define TIME_SERIES_METADATA_MAXIMUM_BLOCK_SAMPLES_OFFSET		8944		// ui4
#define TIME_SERIES_METADATA_MAXIMUM_BLOCK_SAMPLES_NO_ENTRY		0xFFFFFFFF
#define TIME_SERIES_METADATA_MAXIMUM_DIFFERENCE_BYTES_OFFSET		8948		// ui4
#define TIME_SERIES_METADATA_MAXIMUM_DIFFERENCE_BYTES_NO_ENTRY		0xFFFFFFFF
#define TIME_SERIES_METADATA_BLOCK_INTERVAL_OFFSET			8952		// si8
#define TIME_SERIES_METADATA_BLOCK_INTERVAL_NO_ENTRY			-1
#define TIME_SERIES_METADATA_NUMBER_OF_DISCONTINUITIES_OFFSET		8960		// ui4
#define TIME_SERIES_METADATA_NUMBER_OF_DISCONTINUITIES_NO_ENTRY		-1
#define TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCKS_OFFSET		8968		// si8
#define TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCKS_NO_ENTRY		-1
#define TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCK_BYTES_OFFSET	8976		// si8
#define TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCK_BYTES_NO_ENTRY	-1
#define TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_SAMPLES_OFFSET		8984		// si8
#define TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_SAMPLES_NO_ENTRY	-1
#define TIME_SERIES_METADATA_SECTION_2_PROTECTED_REGION_OFFSET		8992
#define TIME_SERIES_METADATA_SECTION_2_PROTECTED_REGION_BYTES		2160
#define TIME_SERIES_METADATA_SECTION_2_DISCRETIONARY_REGION_OFFSET	11152
#define TIME_SERIES_METADATA_SECTION_2_DISCRETIONARY_REGION_BYTES	2160

// Metadata: File Format Constants - Video Section 2 Fields
#define VIDEO_METADATA_HORIZONTAL_RESOLUTION_OFFSET			6664		// si8
#define VIDEO_METADATA_HORIZONTAL_RESOLUTION_NO_ENTRY			-1
#define VIDEO_METADATA_VERTICAL_RESOLUTION_OFFSET			6672		// si8
#define VIDEO_METADATA_VERTICAL_RESOLUTION_NO_ENTRY			-1
#define VIDEO_METADATA_FRAME_RATE_OFFSET				6680		// sf8
#define VIDEO_METADATA_FRAME_RATE_NO_ENTRY				-1.0
#define VIDEO_METADATA_NUMBER_OF_CLIPS_OFFSET				6688
#define VIDEO_METADATA_NUMBER_OF_CLIPS_NO_ENTRY				-1
#define VIDEO_METADATA_MAXIMUM_CLIP_BYTES_OFFSET			6696
#define VIDEO_METADATA_MAXIMUM_CLIP_BYTES_NO_ENTRY			-1
#define VIDEO_METADATA_VIDEO_FORMAT_OFFSET				6704		// utf8[31]
#define VIDEO_METADATA_VIDEO_FORMAT_BYTES				128
#define VIDEO_METADATA_VIDEO_FILE_CRC_OFFSET				6832
#define VIDEO_METADATA_VIDEO_FILE_CRC_NO_ENTRY				CRC_NO_ENTRY
#define VIDEO_METADATA_SECTION_2_PROTECTED_REGION_OFFSET		6836
#define VIDEO_METADATA_SECTION_2_PROTECTED_REGION_BYTES			3236
#define VIDEO_METADATA_SECTION_2_DISCRETIONARY_REGION_OFFSET		10072
#define VIDEO_METADATA_SECTION_2_DISCRETIONARY_REGION_BYTES		3240

// Metadata: File Format Constants - Section 3 Fields
#define METADATA_SECTION_3_OFFSET				13312
#define METADATA_SECTION_3_BYTES				3072
#define METADATA_RECORDING_TIME_OFFSET_OFFSET			13312		//si8
#define METADATA_RECORDING_TIME_OFFSET_NO_ENTRY			UUTC_NO_ENTRY
#define METADATA_DST_START_TIME_OFFSET				13320		//si8
#define METADATA_DST_START_TIME_NO_ENTRY			UUTC_NO_ENTRY
#define METADATA_DST_END_TIME_OFFSET				13328		//si8
#define METADATA_DST_END_TIME_NO_ENTRY				UUTC_NO_ENTRY
#define METADATA_GMT_OFFSET_OFFSET				13336		//si4
#define GMT_OFFSET_NO_ENTRY					-86401
#define METADATA_SUBJECT_NAME_1_OFFSET				13340		// utf8[31]
#define METADATA_SUBJECT_NAME_BYTES				128		// utf8[31]
#define METADATA_SUBJECT_NAME_2_OFFSET				13468		// utf8[31]
#define METADATA_SUBJECT_ID_OFFSET				13596		// utf8[31]
#define METADATA_SUBJECT_ID_BYTES				128		// utf8[31]
#define METADATA_RECORDING_LOCATION_OFFSET			13724		// utf8[63]
#define METADATA_RECORDING_LOCATION_BYTES			512		// utf8[127]
#define METADATA_SECTION_3_PROTECTED_REGION_OFFSET		14236
#define METADATA_SECTION_3_PROTECTED_REGION_BYTES		1124
#define METADATA_SECTION_3_DISCRETIONARY_REGION_OFFSET		15360
#define METADATA_SECTION_3_DISCRETIONARY_REGION_BYTES		1024

// Records: Header Format Constants
#define RECORD_HEADER_BYTES				24
#define RECORD_HEADER_RECORD_CRC_OFFSET			0	// ui4
#define RECORD_HEADER_RECORD_CRC_NO_ENTRY		CRC_NO_ENTRY
#define RECORD_HEADER_TYPE_OFFSET			4	// ascii[4] or ui4
#define RECORD_HEADER_TYPE_NO_ENTRY			0	// as ui4
#define RECORD_HEADER_VERSION_MAJOR_OFFSET		9	// ui1
#define RECORD_HEADER_VERSION_MAJOR_NO_ENTRY		0xff
#define RECORD_HEADER_VERSION_MINOR_OFFSET		10	// ui1
#define RECORD_HEADER_VERSION_MINOR_NO_ENTRY		0xff
#define RECORD_HEADER_ENCRYPTION_OFFSET			11	// si1
#define RECORD_HEADER_BYTES_OFFSET			12	// ui4
#define RECORD_HEADER_BYTES_NO_ENTRY			0
#define RECORD_HEADER_TIME_OFFSET			16	// si8
#define RECORD_HEADER_TIME_NO_ENTRY			UUTC_NO_ENTRY	// si8

// Record Index: Format Constants
#define RECORD_INDEX_BYTES				24
#define RECORD_INDEX_TYPE_OFFSET			0	// ascii[4] or ui4
#define RECORD_INDEX_TYPE_NO_ENTRY			0	// as ui4
#define RECORD_INDEX_VERSION_MAJOR_OFFSET		5	// ui1
#define RECORD_INDEX_VERSION_MAJOR_NO_ENTRY		0xff
#define RECORD_INDEX_VERSION_MINOR_OFFSET		6	// ui1
#define RECORD_INDEX_VERSION_MINOR_NO_ENTRY		0xff
#define RECORD_INDEX_ENCRYPTION_OFFSET			7	// si1
#define RECORD_INDEX_FILE_OFFSET_OFFSET			8	// si8
#define RECORD_INDEX_FILE_OFFSET_NO_ENTRY		-1
#define RECORD_INDEX_TIME_OFFSET			16	// si8
#define RECORD_INDEX_TIME_NO_ENTRY			UUTC_NO_ENTRY

// Time Series Index: Format Constants
#define TIME_SERIES_INDEX_BYTES					56
#define TIME_SERIES_INDEX_FILE_OFFSET_OFFSET			0		// si8
#define TIME_SERIES_INDEX_FILE_OFFSET_NO_ENTRY			-1
#define TIME_SERIES_INDEX_START_TIME_OFFSET			8		// si8
#define TIME_SERIES_INDEX_START_TIME_NO_ENTRY			UUTC_NO_ENTRY
#define TIME_SERIES_INDEX_START_SAMPLE_OFFSET			16		// si8
#define TIME_SERIES_INDEX_START_SAMPLE_NO_ENTRY			-1
#define TIME_SERIES_INDEX_NUMBER_OF_SAMPLES_OFFSET		24		// ui4
#define TIME_SERIES_INDEX_NUMBER_OF_SAMPLES_NO_ENTRY		0xFFFFFFFF
#define TIME_SERIES_INDEX_BLOCK_BYTES_OFFSET			28		// ui4
#define TIME_SERIES_INDEX_BLOCK_BYTES_NO_ENTRY			0xFFFFFFFF
#define TIME_SERIES_INDEX_MAXIMUM_SAMPLE_VALUE_OFFSET		32		// si4
#define TIME_SERIES_INDEX_MAXIMUM_SAMPLE_VALUE_NO_ENTRY		RED_NAN
#define TIME_SERIES_INDEX_MINIMUM_SAMPLE_VALUE_OFFSET		36		// si4
#define TIME_SERIES_INDEX_MINIMUM_SAMPLE_VALUE_NO_ENTRY		RED_NAN
#define TIME_SERIES_INDEX_PROTECTED_REGION_OFFSET		40
#define TIME_SERIES_INDEX_PROTECTED_REGION_BYTES		4
#define TIME_SERIES_INDEX_RED_BLOCK_FLAGS_OFFSET		44		// ui1
#define RED_BLOCK_FLAGS_BYTES					1
#define TIME_SERIES_INDEX_RED_BLOCK_PROTECTED_REGION_OFFSET	45
#define RED_BLOCK_PROTECTED_REGION_BYTES			3
#define TIME_SERIES_INDEX_RED_BLOCK_DISCRETIONARY_REGION_OFFSET	48
#define RED_BLOCK_DISCRETIONARY_REGION_BYTES			8

// Video Index: Format Constants
#define VIDEO_INDEX_BYTES			64
#define VIDEO_INDEX_START_TIME_OFFSET		0	// si8
#define VIDEO_INDEX_START_TIME_NO_ENTRY		UUTC_NO_ENTRY
#define VIDEO_INDEX_END_TIME_OFFSET		8	// si8
#define VIDEO_INDEX_END_TIME_NO_ENTRY		UUTC_NO_ENTRY
#define VIDEO_INDEX_START_FRAME_OFFSET		16	// ui4
#define VIDEO_INDEX_START_FRAME_NO_ENTRY	0xFFFFFFFF
#define VIDEO_INDEX_END_FRAME_OFFSET		20	// ui4
#define VIDEO_INDEX_END_FRAME_NO_ENTRY		0xFFFFFFFF
#define VIDEO_INDEX_FILE_OFFSET_OFFSET		24	// si8
#define VIDEO_INDEX_FILE_OFFSET_NO_ENTRY	-1
#define VIDEO_INDEX_CLIP_BYTES_OFFSET		32	// si8
#define VIDEO_INDEX_CLIP_BYTES_NO_ENTRY		-1
#define VIDEO_INDEX_PROTECTED_REGION_OFFSET	40
#define VIDEO_INDEX_PROTECTED_REGION_BYTES	16
#define VIDEO_INDEX_DISCRETIONARY_REGION_OFFSET	56
#define VIDEO_INDEX_DISCRETIONARY_REGION_BYTES	8




/************************************************************************************/
/************************************  MEF Macros  **********************************/
/************************************************************************************/

#define ABS(x)			( ((x) >= 0) ? (x) : -(x) )
#define HEX_STRING_BYTES(x)	( ((x) + 1) * 3 )



/************************************************************************************/
/**********************************  MEF Structures  ********************************/
/************************************************************************************/

#pragma pack(1)

// Password Structures
typedef struct {
        ui1	level_1_encryption_key[ENCRYPTION_KEY_BYTES];
        ui1	level_2_encryption_key[ENCRYPTION_KEY_BYTES];
        ui1	access_level;
} PASSWORD_DATA;

// Universal Header Structures
typedef struct {
	ui4	header_CRC;
	ui4	body_CRC;
        si1	file_type_string[TYPE_BYTES];
        ui1	mef_version_major;
        ui1	mef_version_minor;
        ui1	byte_order_code;
	si8	start_time;
	si8	end_time;
	si8	number_of_entries;
	si8	maximum_entry_size;
        si4	segment_number;
        si1	channel_name[MEF_BASE_FILE_NAME_BYTES]; // utf8[63], base name only, no extension
        si1	session_name[MEF_BASE_FILE_NAME_BYTES]; // utf8[63], base name only, no extension
	si1	anonymized_name[UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES]; // utf8[63]
        ui1	level_UUID[UUID_BYTES];
	ui1	file_UUID[UUID_BYTES];
        ui1	provenance_UUID[UUID_BYTES];
	ui1	level_1_password_validation_field[PASSWORD_VALIDATION_FIELD_BYTES];
	ui1	level_2_password_validation_field[PASSWORD_VALIDATION_FIELD_BYTES];
        ui1	protected_region[UNIVERSAL_HEADER_PROTECTED_REGION_BYTES];
	ui1	discretionary_region[UNIVERSAL_HEADER_DISCRETIONARY_REGION_BYTES];
} UNIVERSAL_HEADER;


// Metadata Structures
typedef struct {
        si1			section_2_encryption;
        si1			section_3_encryption;
        ui1			protected_region[METADATA_SECTION_1_PROTECTED_REGION_BYTES];
        ui1			discretionary_region[METADATA_SECTION_1_DISCRETIONARY_REGION_BYTES];
} METADATA_SECTION_1;

typedef struct {
        // type-independent fields
        si1			channel_description[METADATA_CHANNEL_DESCRIPTION_BYTES]; // utf8[511];
        si1			session_description[METADATA_SESSION_DESCRIPTION_BYTES]; // utf8[511];
        si8			recording_duration;
        // type-specific fields
	si1			reference_description[METADATA_CHANNEL_DESCRIPTION_BYTES]; // utf8[511];
	si8			acquisition_channel_number;
        sf8			sampling_frequency;
        sf8			low_frequency_filter_setting;
        sf8			high_frequency_filter_setting;
	sf8			notch_filter_frequency_setting;
	sf8			AC_line_frequency;
        sf8			units_conversion_factor;
        si1			units_description[TIME_SERIES_METADATA_UNITS_DESCRIPTION_BYTES];	// utf8[31]
        sf8			maximum_native_sample_value;
        sf8			minimum_native_sample_value;
        si8			start_sample;
	si8			number_of_samples;
	si8			number_of_blocks;
	si8			maximum_block_bytes;
	ui4			maximum_block_samples;
        ui4			maximum_difference_bytes;
	si8			block_interval;
	si8			number_of_discontinuities;
	si8			maximum_contiguous_blocks;
	si8			maximum_contiguous_block_bytes;
        si8			maximum_contiguous_samples;
        ui1			protected_region[TIME_SERIES_METADATA_SECTION_2_PROTECTED_REGION_BYTES];
        ui1			discretionary_region[TIME_SERIES_METADATA_SECTION_2_DISCRETIONARY_REGION_BYTES];
} TIME_SERIES_METADATA_SECTION_2;

typedef struct {
        // type-independent fields
        si1			channel_description[METADATA_CHANNEL_DESCRIPTION_BYTES]; // utf8[511]
        si1			session_description[METADATA_SESSION_DESCRIPTION_BYTES];
        si8			recording_duration;
        // type-specific fields
        si8			horizontal_resolution;
        si8			vertical_resolution;
        sf8			frame_rate;
	si8			number_of_clips;
	si8			maximum_clip_bytes;
        si1			video_format[VIDEO_METADATA_VIDEO_FORMAT_BYTES]; // utf8[31]
	ui4			video_file_CRC;
        ui1			protected_region[VIDEO_METADATA_SECTION_2_PROTECTED_REGION_BYTES];
        ui1			discretionary_region[VIDEO_METADATA_SECTION_2_DISCRETIONARY_REGION_BYTES];
} VIDEO_METADATA_SECTION_2;

typedef struct {
	si8			recording_time_offset;
	si8			DST_start_time;
	si8			DST_end_time;
	si4			GMT_offset;
        si1			subject_name_1[METADATA_SUBJECT_NAME_BYTES]; // utf8[31]
        si1			subject_name_2[METADATA_SUBJECT_NAME_BYTES]; // utf8[31]
        si1			subject_ID[METADATA_SUBJECT_ID_BYTES]; // utf8[31]
        si1			recording_location[METADATA_RECORDING_LOCATION_BYTES]; // utf8[127]
        ui1			protected_region[METADATA_SECTION_3_PROTECTED_REGION_BYTES];
        ui1			discretionary_region[METADATA_SECTION_3_DISCRETIONARY_REGION_BYTES];
} METADATA_SECTION_3;

typedef struct {
        METADATA_SECTION_1		*section_1;
        TIME_SERIES_METADATA_SECTION_2	*time_series_section_2;
	VIDEO_METADATA_SECTION_2	*video_section_2;
        METADATA_SECTION_3		*section_3;
} METADATA;


// Record Structures
typedef struct {
	ui4	record_CRC;
	si1	type_string[TYPE_BYTES];
	ui1	version_major;
	ui1	version_minor;
	si1	encryption;
	ui4	bytes;
	si8	time;
} RECORD_HEADER;

typedef struct {
        si1	type_string[TYPE_BYTES];
        ui1	version_major;
        ui1	version_minor;
        si1	encryption;
        si8	file_offset;
        si8	time;
} RECORD_INDEX;

// Block Indices Structures
typedef struct {
	si8	file_offset;
	si8	start_time;
	si8	start_sample;
	ui4	number_of_samples;
	ui4	block_bytes;
	si4	maximum_sample_value;
	si4	minimum_sample_value;
	ui1	protected_region[TIME_SERIES_INDEX_PROTECTED_REGION_BYTES];
        ui1	RED_block_flags;
        ui1	RED_block_protected_region[RED_BLOCK_PROTECTED_REGION_BYTES];
	ui1	RED_block_discretionary_region[RED_BLOCK_DISCRETIONARY_REGION_BYTES];
} TIME_SERIES_INDEX;

// Frame Indices Structures
typedef struct {
	si8	start_time;
	si8	end_time;
	ui4	start_frame;
	ui4	end_frame;
	si8	file_offset;
	si8	clip_bytes;
	ui1	protected_region[VIDEO_INDEX_PROTECTED_REGION_BYTES];
	ui1	discretionary_region[VIDEO_INDEX_DISCRETIONARY_REGION_BYTES];
} VIDEO_INDEX;

// File Processing Structures
typedef struct {
	si1				close_file;
	si1				free_password_data;  // when freeing FPS
        si8				io_bytes;  // bytes to read or write
        ui4				lock_mode;
	ui4				open_mode;
} FILE_PROCESSING_DIRECTIVES;

typedef struct {
	si1				full_file_name[MEF_FULL_FILE_NAME_BYTES];  // full path including extension
	FILE				*fp;
        si4				fd;
	si8				file_length;
	ui4				file_type_code;
	UNIVERSAL_HEADER		*universal_header;
        FILE_PROCESSING_DIRECTIVES	directives;
	PASSWORD_DATA			*password_data;  // this will often be the same for all files
	METADATA			metadata;
	TIME_SERIES_INDEX		*time_series_indices;
	VIDEO_INDEX			*video_indices;
	ui1				*records;
        RECORD_INDEX			*record_indices;
        ui1				*RED_blocks;
	si8				raw_data_bytes;
	ui1				*raw_data;
} FILE_PROCESSING_STRUCT;

// Session, Channel, Segment Processing Structures
typedef struct {
        si4			channel_type;
	FILE_PROCESSING_STRUCT	*metadata_fps;
	FILE_PROCESSING_STRUCT	*time_series_data_fps;
	FILE_PROCESSING_STRUCT	*time_series_indices_fps;
        FILE_PROCESSING_STRUCT	*video_indices_fps;
	FILE_PROCESSING_STRUCT	*record_data_fps;
	FILE_PROCESSING_STRUCT	*record_indices_fps;
	si1			name[MEF_SEGMENT_BASE_FILE_NAME_BYTES];  // just base name, no extension
	si1			path[MEF_FULL_FILE_NAME_BYTES];  // full path to enclosing directory (channel directory)
	si1			channel_name[MEF_BASE_FILE_NAME_BYTES];  // just base name, no extension
	si1			session_name[MEF_BASE_FILE_NAME_BYTES];  // just base name, no extension
	ui1			level_UUID[UUID_BYTES];
} SEGMENT;

typedef struct {
        si4			channel_type;
	METADATA		metadata;
        FILE_PROCESSING_STRUCT	*record_data_fps;
	FILE_PROCESSING_STRUCT	*record_indices_fps;
        si8			number_of_segments;
	SEGMENT			*segments;
	si1			path[MEF_FULL_FILE_NAME_BYTES];  // full path to enclosing directory (session directory)
	si1			name[MEF_BASE_FILE_NAME_BYTES];  // just base name, no extension
	si1			extension[TYPE_BYTES];  // channel directory extension
	si1			session_name[MEF_BASE_FILE_NAME_BYTES];  // just base name, no extension
	ui1			level_UUID[UUID_BYTES];
	si1			anonymized_name[UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES];
	// below variables refer to segments
	si8			maximum_number_of_records;
	si8			maximum_record_bytes;
	si8			earliest_start_time;
	si8			latest_end_time;
} CHANNEL;

typedef struct {
	METADATA		time_series_metadata;
        si4			number_of_time_series_channels;
	CHANNEL			*time_series_channels;
        METADATA		video_metadata;
        si4			number_of_video_channels;
        CHANNEL			*video_channels;
        FILE_PROCESSING_STRUCT	*record_data_fps;
        FILE_PROCESSING_STRUCT	*record_indices_fps;
	si1			name[MEF_BASE_FILE_NAME_BYTES];  // just base name, no extension
	si1			path[MEF_FULL_FILE_NAME_BYTES];  // full path to enclosing directory
	si1			anonymized_name[UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES];
	ui1			level_UUID[UUID_BYTES];
	// below variables refer to channels
	si8			maximum_number_of_records;
	si8			maximum_record_bytes;
	si8			earliest_start_time;
	si8			latest_end_time;
} SESSION;

// Miscellaneous Structures
typedef struct NODE_STRUCT {
	sf8			val;
	si4			idx;
	struct NODE_STRUCT     *prev, *next;
} NODE;

#pragma pack()

/************************************************************************************/
/**********************************  MEF Prototypes  ********************************/
/************************************************************************************/

// Alignment Function Prototypes
si4			check_all_alignments(const si1 *function, si4 line);
si4			check_metadata_alignment(ui1 *bytes);
si4			check_metadata_section_1_alignment(ui1 *bytes);
si4			check_metadata_section_3_alignment(ui1 *bytes);
si4			check_record_header_alignment(ui1 *bytes);
si4			check_record_indices_alignment(ui1 *bytes);
si1			check_record_structure_alignments(ui1 *bytes);
si4			check_RED_block_header_alignment(ui1 *bytes);
si4			check_time_series_indices_alignment(ui1 *bytes);
si4			check_time_series_metadata_section_2_alignment(ui1 *bytes);
si4			check_universal_header_alignment(ui1 *bytes);
si4			check_video_indices_alignment(ui1 *bytes);
si4			check_video_metadata_section_2_alignment(ui1 *bytes);

// MEF Function Prototypes
si1			all_zeros(ui1 *bytes, si4 field_length);
FILE_PROCESSING_STRUCT	*allocate_file_processing_struct(si8 raw_data_bytes, ui4 file_type_code, FILE_PROCESSING_DIRECTIVES *directives, FILE_PROCESSING_STRUCT *proto_fps, si8 bytes_to_copy);
void			apply_recording_time_offset(si8 *time);
si4			check_password(si1 *password, const si1 *function, si4 line);
si4                     compare_sf8(const void *a, const void * b);
ui1			cpu_endianness(void);
si4			decrypt_metadata(FILE_PROCESSING_STRUCT *fps);
si4			decrypt_records(FILE_PROCESSING_STRUCT *fps);
si4			encrypt_metadata(FILE_PROCESSING_STRUCT *fps);
si4			encrypt_records(FILE_PROCESSING_STRUCT *fps);
si4			extract_path_parts(si1 *full_file_name, si1 *path, si1 *name, si1 *extension);
void			extract_terminal_password_bytes(si1 *password, si1 *password_bytes);
void			fill_empty_password_bytes(si1 *password_bytes);
si8			*find_discontinuity_indices(TIME_SERIES_INDEX *tsi, si8 num_disconts, si8 number_of_blocks);
si8			*find_discontinuity_samples(TIME_SERIES_INDEX *tsi, si8 num_disconts, si8 number_of_blocks, si1 add_tail);
void			force_behavior(ui4 behavior);
void			fps_close(FILE_PROCESSING_STRUCT *fps);
si4			fps_lock(FILE_PROCESSING_STRUCT *fps, si4 lock_type, const si1 *function, si4 line, ui4 behavior_on_fail);
si4			fps_open(FILE_PROCESSING_STRUCT *fps, const si1 *function, si4 line, ui4 behavior_on_fail);
si4			fps_read(FILE_PROCESSING_STRUCT *fps, const si1 *function, si4 line, ui4 behavior_on_fail);
si4			fps_unlock(FILE_PROCESSING_STRUCT *fps, const si1 *function, si4 line, ui4 behavior_on_fail);
si4			fps_write(FILE_PROCESSING_STRUCT *fps, const si1 *function, si4 line, ui4 behavior_on_fail);
void			free_channel(CHANNEL *channel, si4 free_channel_structure);
void			free_file_processing_struct(FILE_PROCESSING_STRUCT *fps);
void			free_segment(SEGMENT *segment, si4 free_segment_structure);
void			free_session(SESSION *session, si4 free_session_structure);
si1			**generate_file_list(si1 **file_list, si4 *num_files, si1 *enclosing_directory, si1 *extension);
si1			*generate_hex_string(ui1 *bytes, si4 num_bytes, si1 *string);
si8			generate_recording_time_offset(si8 recording_start_time_uutc, si4 GMT_offset);
si1			*generate_segment_name(FILE_PROCESSING_STRUCT *fps, si1 *segment_name);
ui1			*generate_UUID(ui1 *uuid);
FILE_PROCESSING_DIRECTIVES *initialize_file_processing_directives(FILE_PROCESSING_DIRECTIVES *directives);
void			initialize_MEF_globals(void);
si4			initialize_meflib(void);
si4			initialize_metadata(FILE_PROCESSING_STRUCT *fps);
si4			initialize_universal_header(FILE_PROCESSING_STRUCT *fps, si1 generate_level_UUID, si1 generate_file_UUID, si1 originating_file);
si1			*local_date_time_string(si8 uutc_time, si1 *time_str);
si8			MEF_pad(ui1 *buffer, si8 content_len, ui4 alignment);
si4			MEF_sprintf(si1 *target, si1 *format, ...);
void			MEF_snprintf(si1 *target, si4 target_field_bytes, si1 *format, ...);
si4			MEF_strcat(si1 *target_string, si1 *source_string);
si4			MEF_strcpy(si1 *target_string, si1 *source_string);
void			MEF_strncat(si1 *target_string, si1 *source_string, si4 target_field_bytes);
void			MEF_strncpy(si1 *target_string, si1 *source_string, si4 target_field_bytes);
si1			*numerical_fixed_width_string(si1 *string, si4 string_bytes, si4 number);
si4			offset_record_index_times(FILE_PROCESSING_STRUCT *fps, si4 action);
si4			offset_time_series_index_times(FILE_PROCESSING_STRUCT *fps, si4 action);
si4			offset_universal_header_times(FILE_PROCESSING_STRUCT *fps, si4 action);
si4			offset_video_index_times(FILE_PROCESSING_STRUCT *fps, si4 action);
PASSWORD_DATA		*process_password_data(si1 *unspecified_password, si1 *level_1_password, si1 *level_2_password, UNIVERSAL_HEADER *universal_header);
void			proportion_filt(sf8 *x, sf8 *px, si8 len, sf8 prop, si4 span);
ui1			random_byte(ui4 *m_w, ui4 *m_z);
CHANNEL			*read_MEF_channel(CHANNEL *channel, si1 *chan_path, si4 channel_type, si1 *password, PASSWORD_DATA *password_data, si1 read_time_series_data, si1 read_record_data);
FILE_PROCESSING_STRUCT	*read_MEF_file(FILE_PROCESSING_STRUCT *fps, si1 *file_name, si1 *password, PASSWORD_DATA *password_data, FILE_PROCESSING_DIRECTIVES *directives, ui4 behavior_on_fail);
SEGMENT			*read_MEF_segment(SEGMENT *segment, si1 *seg_path, si4 channel_type, si1 *password, PASSWORD_DATA *password_data, si1 read_time_series_data, si1 read_record_data);
SESSION			*read_MEF_session(SESSION *session, si1 *sess_path, si1 *password, PASSWORD_DATA *password_data, si1 read_time_series_data, si1 read_record_data);
si4			reallocate_file_processing_struct(FILE_PROCESSING_STRUCT *fps, si8 raw_data_bytes);
si4                     remove_line_noise(si4 *data, si8 n_samps, sf8 sampling_frequency, sf8 line_frequency, sf8 *template);
void			remove_line_noise_adaptive(si4 *data, si8 n_samps, sf8 sampling_frequency, sf8 line_frequency, si4 n_cycles);
void			remove_recording_time_offset(si8 *time);
void			show_file_processing_struct(FILE_PROCESSING_STRUCT *fps);
void			show_metadata(FILE_PROCESSING_STRUCT *fps);
void			show_password_data(FILE_PROCESSING_STRUCT *fps);
void			show_record(RECORD_HEADER *record_header, ui4 record_number, PASSWORD_DATA *pwd);
void			show_records(FILE_PROCESSING_STRUCT *fps);
void			show_universal_header(FILE_PROCESSING_STRUCT *fps);
si4			sort_by_idx(const void *n1, const void *n2);
si4			sort_by_val(const void *n1, const void *n2);
sf8			val_equals_prop(NODE *curr_node, NODE *prop_node);
si4			write_MEF_file(FILE_PROCESSING_STRUCT *fps);

#ifdef _WIN32
	void 		slash_to_backslash(si1*);
#endif

/************************************************************************************/
/**************************************  FILTER  ************************************/
/************************************************************************************/

// ATTRIBUTION
//
// Some of the filter code was adapted from Matlab functions.
// www.mathworks.com
//
// The c code was written entirely from scratch.
//
// NOTE: This code requres long double (sf16) math.
// It often requires an explicit compiler instruction to implement true long floating point math.
// in icc: "-Qoption,cpp,--extended_float_type"


// Constants
#define FILT_LOWPASS_TYPE				1
#define FILT_BANDPASS_TYPE				2
#define FILT_HIGHPASS_TYPE				3
#define FILT_BANDSTOP_TYPE				4
#define FILT_TYPE_DEFAULT				FILT_LOWPASS_TYPE
#define FILT_ORDER_DEFAULT				5
#define FILT_MAX_ORDER 					10
#define FILT_BAD_FILTER					-1
#define FILT_BAD_DATA					-2
#define FILT_EPS_SF8         				2.22045e-16
#define FILT_EPS_SF16        				1.0842e-19
#define FILT_RADIX           				2
#define FILT_ZERO            				((sf16) 0.0)
#define FILT_ONE					((sf16) 1.0)

// Macros
#define FILT_ABS(x)          				((x) >= FILT_ZERO ? (x) : (-x))
#define FILT_SIGN(x, y)      				((y) >= FILT_ZERO ? FILT_ABS(x) : -FILT_ABS(x))

// Typedefs & Structures
typedef struct {
        si4	order;
        si4	poles;
        si4	type;
        sf8	sampling_frequency;
        si8	data_length;
        sf8	cutoffs[2];
        sf8	*numerators;
        sf8	*denominators;
        sf8	*initial_conditions;
        si4	*orig_data;
	si4	*filt_data;
        sf8	*sf8_filt_data;
        sf8	*sf8_buffer;
} FILT_PROCESSING_STRUCT;

typedef struct {
	sf16	real;
	sf16	imag;
} FILT_LONG_COMPLEX;

// Prototypes
void			FILT_balance(sf16 **a, si4 poles);
si4			FILT_butter(FILT_PROCESSING_STRUCT *filtps);
void			FILT_complex_divl(FILT_LONG_COMPLEX *a, FILT_LONG_COMPLEX *b, FILT_LONG_COMPLEX *quotient);
void			FILT_complex_expl(FILT_LONG_COMPLEX *exponent, FILT_LONG_COMPLEX *ans);
void			FILT_complex_multl(FILT_LONG_COMPLEX *a, FILT_LONG_COMPLEX *b, FILT_LONG_COMPLEX *product);
void			FILT_elmhes(sf16 **a, si4 poles);
si4			FILT_filtfilt(FILT_PROCESSING_STRUCT *filtps);
void			FILT_free_processing_struct(FILT_PROCESSING_STRUCT *filtps, si1 free_orig_data, si1 free_filt_data);
FILT_PROCESSING_STRUCT	*FILT_initialize_processing_struct(si4 order, si4 type, sf8 samp_freq, si8 data_len, si1 alloc_orig_data, si1 alloc_filt_data, sf8 cutoff_1, ...);
void			FILT_generate_initial_conditions(FILT_PROCESSING_STRUCT *filtps);
void			FILT_hqr(sf16 **a, si4 poles, FILT_LONG_COMPLEX *eigs);
void			FILT_invert_matrix(sf16 **a, sf16 **inv_a, si4 order);
void			FILT_mat_multl(void *a, void *b, void *product, si4 outer_dim1, si4 inner_dim, si4 outer_dim2);
void			FILT_unsymmeig(sf16 **a, si4 poles, FILT_LONG_COMPLEX *eigs);



/************************************************************************************/
/****************************************  RED  *************************************/
/************************************************************************************/

// Constants
#define RED_BLOCK_HEADER_BYTES				304
#define RED_BLOCK_COMPRESSED_DATA_OFFSET		RED_BLOCK_HEADER_BYTES
#define RED_BLOCK_BLOCK_CRC_OFFSET			0		// ui4
#define RED_BLOCK_FLAGS_OFFSET				4		// ui1
#define RED_BLOCK_PROTECTED_REGION_OFFSET		5
// #define RED_BLOCK_PROTECTED_REGION_BYTES		3		// defined above because used in TIME_SERIES_INDEX structure
#define RED_BLOCK_DISCRETIONARY_REGION_OFFSET		8
// #define RED_BLOCK_DISCRETIONARY_REGION_BYTES		8		// defined above because used in TIME_SERIES_INDEX structure
#define RED_BLOCK_DETREND_SLOPE_OFFSET			16		// sf4
#define RED_BLOCK_DETREND_INTERCEPT_OFFSET		20		// sf4
#define RED_BLOCK_SCALE_FACTOR_OFFSET			24		// sf4
#define RED_BLOCK_DIFFERENCE_BYTES_OFFSET		28		// ui4
#define RED_BLOCK_NUMBER_OF_SAMPLES_OFFSET		32		// ui4
#define RED_BLOCK_BLOCK_BYTES_OFFSET			36		// ui4
#define RED_BLOCK_START_TIME_OFFSET			40		// si8
#define RED_BLOCK_STATISTICS_OFFSET			48		// ui1
#define RED_BLOCK_STATISTICS_BYTES			256

// RED Codec: Flag Masks
#define RED_DISCONTINUITY_MASK				((ui1) 1)	// Bit 0
#define RED_LEVEL_1_ENCRYPTION_MASK			((ui1) 2)	// Bit 1
#define RED_LEVEL_2_ENCRYPTION_MASK			((ui1) 4)	// Bit 2

// RED Codec: Reserved Values
#define RED_NAN						((si4) 0x80000000)
#define RED_NEGATIVE_INFINITY				((si4) 0x80000001)
#define RED_POSITIVE_INFINITY				((si4) 0x7FFFFFFF)
#define RED_MAXIMUM_SAMPLE_VALUE			((si4) 0x7FFFFFFE)
#define RED_MINIMUM_SAMPLE_VALUE			((si4) 0x80000002)

// RED Codec: Range Encoding Constants
#define TOP_VALUE					((ui4) 0x80000000)
#define TOP_VALUE_MINUS_1				((ui4) 0x7FFFFFFF)
#define CARRY_CHECK					((ui4) 0x7F800000)
#define SHIFT_BITS					23
#define EXTRA_BITS					7
#define BOTTOM_VALUE					((ui4) 0x800000)

// RED Codec: Compression Modes
#define RED_DECOMPRESSION				0  // any non-zero value is compression
#define	RED_LOSSLESS_COMPRESSION			1  // lossless (default)
#define RED_FIXED_SCALE_FACTOR				2  // apply this scale factor to the block, 1.0 results in lossless compression;
#define	RED_FIXED_COMPRESSION_RATIO			4  // e.g. 20% of original si4 size is 0.2 - if lossless satifisfies, no compression is done
#define RED_MEAN_RESIDUAL_RATIO				8  // sum(abs((scaled_data - original_data))) / sum(abs(original_data)), e.g. 5% difference is 0.05

// RED Codec: defaults
#define RED_COMPRESSION_MODE_DEFAULT				RED_LOSSLESS_COMPRESSION
#define RED_ENCRYPTION_LEVEL_DEFAULT				NO_ENCRYPTION
#define RED_SCALE_FACTOR_DEFAULT				1.0
#define RED_DETREND_SLOPE_DEFAULT				0.0
#define RED_DETREND_INTERCEPT_DEFAULT				0.0
#define RED_MAXIMUM_ROUNDS_PER_BLOCK_DEFAULT			20  // number of times the lossy compression routines will try to achieve the goal compression parameters
#define RED_RETURN_LOSSY_DATA_DEFAULT				MEF_FALSE
#define RED_DETREND_DATA_DEFAULT				MEF_FALSE
#define RED_VALIDATE_BLOCK_CRC_DEFAULT				MEF_FALSE
#define RED_DISCONTINUITY_DEFAULT				MEF_TRUE  // if true, the dicontinuity flag will be set to true in RED_allocate_processing_struct().
									  // Useful in combination with the "reset_discontinuity" directive as the first block in a
									  // segment is always a discontiuity.
#define RED_RESET_DISCONTINUITY_DEFAULT				MEF_TRUE
#define RED_GOAL_COMPRESSION_RATIO_DEFAULT			0.05
#define RED_GOAL_MEAN_RESIDUAL_RATIO_DEFAULT			0.05
#define RED_GOAL_TOLERANCE_DEFAULT				0.005
#define RED_REQUIRE_NORMALITY_DEFAULT				MEF_TRUE
#define RED_NORMAL_CORRELATION_DEFAULT				0.5  // range -1.0 to 1.0 with 1.0 being perfect

// RED Codec: Macros
#define RED_MAX_DIFFERENCE_BYTES(x)	(x * 5)	// full si4 plus 1 keysample flag byte per sample
#define RED_MAX_COMPRESSED_BYTES(x, y)	((RED_MAX_DIFFERENCE_BYTES(x) + RED_BLOCK_HEADER_BYTES + 7) * y) // no compression plus header plus maximum pad bytes, foy y blocks

// RED UPDATE RPS POINTERS FLAGS
#define	RED_UPDATE_ORIGINAL_PTR		1
#define	RED_UPDATE_BLOCK_HEADER_PTR	2
#define	RED_UPDATE_DECOMPRESSED_PTR	4

// Normal cumulative distribution fucntion values from -3 to +3 standard deviations in 0.1 sigma steps
#define RED_NORMAL_CDF_TABLE_ENTRIES	61
#define RED_NORMAL_CDF_TABLE	      {	0.00134989803163010, 0.00186581330038404, 0.00255513033042794, 0.00346697380304067, \
					0.00466118802371875, 0.00620966532577614, 0.00819753592459614, 0.01072411002167580, \
					0.01390344751349860, 0.01786442056281660, 0.02275013194817920, 0.02871655981600180, \
					0.03593031911292580, 0.04456546275854310, 0.05479929169955800, 0.06680720126885810, \
					0.08075665923377110, 0.09680048458561040, 0.11506967022170800, 0.13566606094638300, \
					0.15865525393145700, 0.18406012534676000, 0.21185539858339700, 0.24196365222307300, \
					0.27425311775007400, 0.30853753872598700, 0.34457825838967600, 0.38208857781104700, \
					0.42074029056089700, 0.46017216272297100, 0.50000000000000000, 0.53982783727702900, \
					0.57925970943910300, 0.61791142218895300, 0.65542174161032400, 0.69146246127401300, \
					0.72574688224992600, 0.75803634777692700, 0.78814460141660300, 0.81593987465324100, \
					0.84134474606854300, 0.86433393905361700, 0.88493032977829200, 0.90319951541439000, \
					0.91924334076622900, 0.93319279873114200, 0.94520070830044200, 0.95543453724145700, \
					0.96406968088707400, 0.97128344018399800, 0.97724986805182100, 0.98213557943718300, \
					0.98609655248650100, 0.98927588997832400, 0.99180246407540400, 0.99379033467422400, \
					0.99533881197628100, 0.99653302619695900, 0.99744486966957200, 0.99813418669961600, \
					0.99865010196837000 }

#define RED_SUM_NORMAL_CDF		30.5
#define RED_SUM_SQ_NORMAL_CDF		24.864467406647070

// Typedefs & Structures
typedef struct {
	ui4	block_CRC;
	ui1	flags;
        ui1	protected_region[RED_BLOCK_PROTECTED_REGION_BYTES];
        ui1	discretionary_region[RED_BLOCK_DISCRETIONARY_REGION_BYTES];
        sf4	detrend_slope;
        sf4	detrend_intercept;
	sf4	scale_factor;
	ui4	difference_bytes;
	ui4	number_of_samples;
	ui4	block_bytes;
	si8	start_time;
	ui1	statistics[RED_BLOCK_STATISTICS_BYTES];
} RED_BLOCK_HEADER;

typedef struct {
        si1			encryption_level;  // encryption level for data blocks, passed in compression, returned in decompression
        si1			discontinuity;  // set if block is first after a discontinuity, passed in compression, returned in decompression
        si1			detrend_data;  // set if block is to be detrended (somewhat useful in lossless, more useful in lossy compression)
	si1			return_lossy_data;  // if set, lossy data returned in decompressed_data during lossy compression
	si1			reset_discontinuity;  // if discontinuity directive == MEF_TRUE, reset to MEF_FALSE after compressing the block
        si1			require_normality;  // in lossy compression, lossless compression will be performed in blocks whose samples are not approximately normally distributed
        sf8			normal_correlation;  // if require_normality is set, the correlation of the sample distribution with a normal distribution must be >= this number (range -1.0 to 1.0)
} RED_PROCESSING_DIRECTIVES;

typedef struct {
	ui1			mode;  // compression mode
	sf8			goal_compression_ratio;  // goal value passed
	sf8			actual_compression_ratio;  // actual value returned in RED_FIXED_COMPRESSION_RATIO mode
	sf8			goal_mean_residual_ratio;  // goal value passed
	sf8			actual_mean_residual_ratio;  // actual value returned in RED_MEAN_RESIDUAL_RATIO mode
	sf8			goal_tolerance;  // tolerance for lossy compression mode goal, value of <= 0.0 uses default values, which are returned
	si4			maximum_rounds_per_block;  // maximum loops to attain goal compression
} RED_COMPRESSION_PARAMETERS;

typedef struct {
        ui4				counts[RED_BLOCK_STATISTICS_BYTES + 1];  // used by RED_encode & RED_decode
        PASSWORD_DATA			*password_data;  // passed in compression & decompression
	RED_COMPRESSION_PARAMETERS	compression;
	RED_PROCESSING_DIRECTIVES	directives;
        si1				*difference_buffer;  // passed in both compression & decompression
        ui1				*compressed_data;  // passed in decompression, returned in compression, should not be updated
	RED_BLOCK_HEADER		*block_header; // points to beginning of current block within compressed_data array, updatable
        si4				*decompressed_data;  // returned in decompression or if lossy data requested, used in some compression modes, should not be updated
        si4				*decompressed_ptr;  // points to beginning of current block within decompressed_data array, updatable
        si4				*original_data;  // passed in compression, should not be updated
        si4				*original_ptr;  // points to beginning of current block within original_data array, updatable
        si4				*detrended_buffer;  // used if needed in compression, size of decompressed block
        si4				*scaled_buffer;  // used if needed in compression, size of decompressed block
} RED_PROCESSING_STRUCT;

// Function Prototypes
RED_PROCESSING_STRUCT	*RED_allocate_processing_struct(si8 original_data_size, si8 compressed_data_size, si8 decompressed_data_size, si8 difference_buffer_size, si8 detrended_buffer_size, si8 scaled_buffer_size, PASSWORD_DATA *password_data);
sf8			RED_calculate_mean_residual_ratio(si4 *original_data, si4 *lossy_data, ui4 n_samps);
si1			RED_check_RPS_allocation(RED_PROCESSING_STRUCT *rps);
void			RED_decode(RED_PROCESSING_STRUCT *rps);
si4			*RED_detrend(RED_PROCESSING_STRUCT *rps, si4 *input_buffer, si4 *output_buffer);
void			RED_encode(RED_PROCESSING_STRUCT *rps);
void			RED_encode_exec(RED_PROCESSING_STRUCT *rps, si4 *input_buffer, si1 input_is_detrended);
void			RED_encode_lossy(RED_PROCESSING_STRUCT *rps);
void			RED_filter(FILT_PROCESSING_STRUCT *filtps);
void			RED_find_extrema(si4 *buffer, si8 number_of_samples, TIME_SERIES_INDEX *tsi);
void			RED_free_processing_struct(RED_PROCESSING_STRUCT *rps);
void			RED_generate_lossy_data(RED_PROCESSING_STRUCT *rps, si4 *input_buffer, si4 *output_buffer, si1 input_is_detrended);
sf8			*RED_initialize_normal_CDF_table(si4 global_flag);
si4			*RED_retrend(RED_PROCESSING_STRUCT *rps, si4 *input_buffer, si4 *output_buffer);
si4			RED_round(sf8 val);
si4			*RED_scale(RED_PROCESSING_STRUCT *rps, si4 *input_buffer, si4 *output_buffer);
void			RED_show_block_header(RED_BLOCK_HEADER *bh);
sf8			RED_test_normality(si4 *data, ui4 n_samps);
si4			*RED_unscale(RED_PROCESSING_STRUCT *rps, si4 *input_buffer, si4 *output_buffer);
RED_BLOCK_HEADER	*RED_update_RPS_pointers(RED_PROCESSING_STRUCT *rps, ui1 flags);


/************************************************************************************/
/****************************************  CRC  *************************************/
/************************************************************************************/

// Constants
#define CRC_BYTES				4
#define	KOOPMAN32				0xEB31D82E
#define CRC_TABLE_ENTRIES			256
#define CRC_START_VALUE				0xFFFFFFFF


#define CRC_KOOPMAN32_KEY     {	0x00000000, 0x9695C4CA, 0xFB4839C9, 0x6DDDFD03, \
				0x20F3C3CF, 0xB6660705, 0xDBBBFA06, 0x4D2E3ECC, \
				0x41E7879E, 0xD7724354, 0xBAAFBE57, 0x2C3A7A9D, \
				0x61144451, 0xF781809B, 0x9A5C7D98, 0x0CC9B952, \
				0x83CF0F3C, 0x155ACBF6, 0x788736F5, 0xEE12F23F, \
				0xA33CCCF3, 0x35A90839, 0x5874F53A, 0xCEE131F0, \
				0xC22888A2, 0x54BD4C68, 0x3960B16B, 0xAFF575A1, \
				0xE2DB4B6D, 0x744E8FA7,	0x199372A4, 0x8F06B66E, \
				0xD1FDAE25, 0x47686AEF, 0x2AB597EC, 0xBC205326, \
				0xF10E6DEA, 0x679BA920, 0x0A465423, 0x9CD390E9, \
				0x901A29BB, 0x068FED71, 0x6B521072, 0xFDC7D4B8, \
				0xB0E9EA74, 0x267C2EBE, 0x4BA1D3BD, 0xDD341777, \
				0x5232A119, 0xC4A765D3, 0xA97A98D0, 0x3FEF5C1A, \
				0x72C162D6, 0xE454A61C, 0x89895B1F, 0x1F1C9FD5, \
				0x13D52687, 0x8540E24D, 0xE89D1F4E, 0x7E08DB84, \
				0x3326E548, 0xA5B32182, 0xC86EDC81, 0x5EFB184B, \
				0x7598EC17, 0xE30D28DD, 0x8ED0D5DE, 0x18451114, \
				0x556B2FD8, 0xC3FEEB12, 0xAE231611, 0x38B6D2DB, \
				0x347F6B89, 0xA2EAAF43, 0xCF375240, 0x59A2968A, \
				0x148CA846, 0x82196C8C, 0xEFC4918F, 0x79515545, \
				0xF657E32B, 0x60C227E1, 0x0D1FDAE2, 0x9B8A1E28, \
				0xD6A420E4, 0x4031E42E, 0x2DEC192D, 0xBB79DDE7, \
				0xB7B064B5, 0x2125A07F, 0x4CF85D7C, 0xDA6D99B6, \
				0x9743A77A, 0x01D663B0, 0x6C0B9EB3, 0xFA9E5A79, \
				0xA4654232, 0x32F086F8, 0x5F2D7BFB, 0xC9B8BF31, \
				0x849681FD, 0x12034537, 0x7FDEB834, 0xE94B7CFE, \
				0xE582C5AC, 0x73170166, 0x1ECAFC65, 0x885F38AF, \
				0xC5710663, 0x53E4C2A9, 0x3E393FAA, 0xA8ACFB60, \
				0x27AA4D0E, 0xB13F89C4, 0xDCE274C7, 0x4A77B00D, \
				0x07598EC1, 0x91CC4A0B, 0xFC11B708, 0x6A8473C2, \
				0x664DCA90, 0xF0D80E5A, 0x9D05F359, 0x0B903793, \
				0x46BE095F, 0xD02BCD95, 0xBDF63096, 0x2B63F45C, \
				0xEB31D82E, 0x7DA41CE4, 0x1079E1E7, 0x86EC252D, \
				0xCBC21BE1, 0x5D57DF2B, 0x308A2228, 0xA61FE6E2, \
				0xAAD65FB0, 0x3C439B7A, 0x519E6679, 0xC70BA2B3, \
				0x8A259C7F, 0x1CB058B5, 0x716DA5B6, 0xE7F8617C, \
				0x68FED712, 0xFE6B13D8, 0x93B6EEDB, 0x05232A11, \
				0x480D14DD, 0xDE98D017, 0xB3452D14, 0x25D0E9DE, \
				0x2919508C, 0xBF8C9446, 0xD2516945, 0x44C4AD8F, \
				0x09EA9343, 0x9F7F5789, 0xF2A2AA8A, 0x64376E40, \
				0x3ACC760B, 0xAC59B2C1, 0xC1844FC2, 0x57118B08, \
				0x1A3FB5C4, 0x8CAA710E, 0xE1778C0D, 0x77E248C7, \
				0x7B2BF195, 0xEDBE355F, 0x8063C85C, 0x16F60C96, \
				0x5BD8325A, 0xCD4DF690, 0xA0900B93, 0x3605CF59, \
				0xB9037937, 0x2F96BDFD, 0x424B40FE, 0xD4DE8434, \
				0x99F0BAF8, 0x0F657E32, 0x62B88331, 0xF42D47FB, \
				0xF8E4FEA9, 0x6E713A63, 0x03ACC760, 0x953903AA, \
				0xD8173D66, 0x4E82F9AC, 0x235F04AF, 0xB5CAC065, \
				0x9EA93439, 0x083CF0F3, 0x65E10DF0, 0xF374C93A, \
				0xBE5AF7F6, 0x28CF333C, 0x4512CE3F, 0xD3870AF5, \
				0xDF4EB3A7, 0x49DB776D, 0x24068A6E, 0xB2934EA4, \
				0xFFBD7068, 0x6928B4A2, 0x04F549A1, 0x92608D6B, \
				0x1D663B05, 0x8BF3FFCF, 0xE62E02CC, 0x70BBC606, \
				0x3D95F8CA, 0xAB003C00, 0xC6DDC103, 0x504805C9, \
				0x5C81BC9B, 0xCA147851, 0xA7C98552, 0x315C4198, \
				0x7C727F54, 0xEAE7BB9E, 0x873A469D, 0x11AF8257, \
				0x4F549A1C, 0xD9C15ED6, 0xB41CA3D5, 0x2289671F, \
				0x6FA759D3, 0xF9329D19, 0x94EF601A, 0x027AA4D0, \
				0x0EB31D82, 0x9826D948, 0xF5FB244B, 0x636EE081, \
				0x2E40DE4D, 0xB8D51A87, 0xD508E784, 0x439D234E, \
				0xCC9B9520, 0x5A0E51EA, 0x37D3ACE9, 0xA1466823, \
				0xEC6856EF, 0x7AFD9225, 0x17206F26, 0x81B5ABEC, \
				0x8D7C12BE, 0x1BE9D674, 0x76342B77, 0xE0A1EFBD, \
				0xAD8FD171, 0x3B1A15BB, 0x56C7E8B8, 0xC0522C72 }

// Function Prototypes
ui4	CRC_calculate(ui1 *block_ptr, si8 block_bytes);
ui4	*CRC_initialize_table(si4 global_flag);
ui4	CRC_update(ui1 *block_ptr, si8 block_bytes, ui4 current_crc);
si4	CRC_validate(ui1 *block_ptr, si8 block_bytes, ui4 crc_to_validate);


/************************************************************************************/
/**************************************  UTF-8  *************************************/
/************************************************************************************/

// ATTRIBUTION
//
// Basic UTF-8 manipulation routines
// by Jeff Bezanson
// placed in the public domain Fall 2005
//
// "This code is designed to provide the utilities you need to manipulate
// UTF-8 as an internal string encoding. These functions do not perform the
// error checking normally needed when handling UTF-8 data, so if you happen
// to be from the Unicode Consortium you will want to flay me alive.
// I do this because error checking can be performed at the boundaries (I/O),
// with these routines reserved for higher performance on data known to be
// valid."
//
// downloaded from http://www.cprogramming.com
//
// Minor modifications for compatibility with the MEF Library.


// Macros
#define isutf(c) (((c) & 0xC0) != 0x80) // is c the start of a utf8 sequence?

#define OFFSETS_FROM_UTF8_TABLE_ENTRIES	6
#define OFFSETS_FROM_UTF8	{ 0x00000000UL, 0x00003080UL, 0x000E2080UL, 0x03C82080UL, 0xFA082080UL, 0x82082080UL }

#define TRAILING_BYTES_FOR_UTF8_TABLE_ENTRIES	256
#define TRAILING_BYTES_FOR_UTF8	{	0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, \
					0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, \
					0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, \
					0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, \
					0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, \
					0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, \
					1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, \
					2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2, 3,3,3,3,3,3,3,3,4,4,4,4,5,5,5,5 }



// Function Prototypes
si4	UTF8_charnum(si1 *s, si4 offset);  // byte offset to character number
void	UTF8_dec(si1 *s, si4 *i);  // move to previous character
si4	UTF8_escape(si1 *buf, si4 sz, si1 *src, si4 escape_quotes);  // convert UTF-8 "src" to ASCII with escape sequences.
si4	UTF8_escape_wchar(si1 *buf, si4 sz, ui4 ch);  // given a wide character, convert it to an ASCII escape sequence stored in buf, where buf is "sz" bytes. returns the number of characters output
si4	UTF8_fprintf(FILE *stream, si1 *fmt, ...);  // fprintf() where the format string and arguments may be in UTF-8. You can avoid this function and just use ordinary printf() if the current locale is UTF-8.
si4	UTF8_hex_digit(si1 c);  // utility predicates used by the above
void	UTF8_inc(si1 *s, si4 *i);  // move to next character
ui4	*UTF8_initialize_offsets_from_UTF8_table(si4 global_flag);
si1	*UTF8_initialize_trailing_bytes_for_UTF8_table(si4 global_flag);
si4	UTF8_is_locale_utf8(si1 *locale);  // boolean function returns if locale is UTF-8, 0 otherwise
si1	*UTF8_memchr(si1 *s, ui4 ch, size_t sz, si4 *charn);  // same as the above, but searches a buffer of a given size instead of a NUL-terminated string.
ui4	UTF8_nextchar(si1 *s, si4 *i);  // return next character, updating an index variable
si4	UTF8_octal_digit(si1 c);  // utility predicates used by the above
si4	UTF8_offset(si1 *str, si4 charnum);  // character number to byte offset
si4	UTF8_printf(si1 *fmt, ...);  // printf() where the format string and arguments may be in UTF-8. You can avoid this function and just use ordinary printf() if the current locale is UTF-8.
si4	UTF8_read_escape_sequence(si1 *str, ui4 *dest);  // assuming src points to the character after a backslash, read an escape sequence, storing the result in dest and returning the number of input characters processed
si4	UTF8_seqlen(si1 *s);  // returns length of next UTF-8 sequence
si1	*UTF8_strchr(si1 *s, ui4 ch, si4 *charn);  // return a pointer to the first occurrence of ch in s, or NULL if not found. character index of found character returned in *charn.
si4	UTF8_strlen(si1 *s);  // count the number of characters in a UTF-8 string
si4	UTF8_toucs(ui4 *dest, si4 sz, si1 *src, si4 srcsz);  // convert UTF-8 data to wide character
si4	UTF8_toutf8(si1 *dest, si4 sz, ui4 *src, si4 srcsz);  // convert wide character to UTF-8 data
si4	UTF8_unescape(si1 *buf, si4 sz, si1 *src);  // convert a string "src" containing escape sequences to UTF-8 if escape_quotes is nonzero, quote characters will be preceded by  backslashes as well.
si4	UTF8_vfprintf(FILE *stream, si1 *fmt, va_list ap);    // called by UTF8_fprintf()
si4	UTF8_vprintf(si1 *fmt, va_list ap);  // called by UTF8_printf()
si4	UTF8_wc_toutf8(si1 *dest, ui4 ch);  // single character to UTF-8



/************************************************************************************/
/***************************************  AES-128  **********************************/
/************************************************************************************/


// ATRIBUTION
//
// Advanced Encryption Standard implementation in C.
// By Niyaz PK
// E-mail: niyazlife@gmail.com
// Downloaded from Website: www.hoozi.com
//
// "This is the source code for encryption using the latest AES algorithm.
// AES algorithm is also called Rijndael algorithm. AES algorithm is
// recommended for non-classified use by the National Institute of Standards
// and Technology (NIST), USA. Now-a-days AES is being used for almost
// all encryption applications all around the world."
//
// For the complete description of the algorithm, see:
// http://www.csrc.nist.gov/publications/fips/fips197/fips-197.pdf
//
// THE CODE IN THIS FILE IS SET FOR 128-BIT ENCRYPTION / DECRYPTION ONLY
//
// Minor modifications for compatibility with the MEF Library.


#define AES_NR	10	// The number of rounds in AES Cipher
#define AES_NK	4	// The number of 32 bit words in the key
#define AES_NB	4	// The number of columns comprising a state in AES. This is a constant in AES.
#define AES_XTIME(x) ((x<<1) ^ (((x>>7) & 1) * 0x1b)) // AES_XTIME is a macro that finds the product of {02} and the argument to AES_XTIME modulo {1b}
#define AES_MULTIPLY(x,y) (((y & 1) * x) ^ ((y>>1 & 1) * AES_XTIME(x)) ^ ((y>>2 & 1) * AES_XTIME(AES_XTIME(x))) ^ ((y>>3 & 1) * AES_XTIME(AES_XTIME(AES_XTIME(x)))) ^ ((y>>4 & 1) * AES_XTIME(AES_XTIME(AES_XTIME(AES_XTIME(x)))))) // Multiplty is a macro used to multiply numbers in the field GF(2^8)

#define AES_SBOX_ENTRIES	256
#define AES_SBOX {	0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76, \
			0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0, \
			0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15, \
			0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75, \
			0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84, \
			0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf, \
			0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8, \
			0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2, \
			0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73, \
			0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb, \
			0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79, \
			0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08, \
			0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a, \
			0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e, \
			0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf, \
			0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 }

#define AES_RSBOX_ENTRIES	256
#define AES_RSBOX {	0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb, \
			0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb, \
			0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e, \
			0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25, \
			0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92, \
			0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84, \
			0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06, \
			0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b, \
			0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73, \
			0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e, \
			0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b, \
			0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4, \
			0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f, \
			0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef, \
			0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61, \
			0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d };

#define AES_RCON_ENTRIES	255
#define AES_RCON {	0x8d, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a, \
			0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91, 0x39, \
			0x72, 0xe4, 0xd3, 0xbd, 0x61, 0xc2, 0x9f, 0x25, 0x4a, 0x94, 0x33, 0x66, 0xcc, 0x83, 0x1d, 0x3a, \
			0x74, 0xe8, 0xcb, 0x8d, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, \
			0xab, 0x4d, 0x9a, 0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, \
			0xc5, 0x91, 0x39, 0x72, 0xe4, 0xd3, 0xbd, 0x61, 0xc2, 0x9f, 0x25, 0x4a, 0x94, 0x33, 0x66, 0xcc, \
			0x83, 0x1d, 0x3a, 0x74, 0xe8, 0xcb, 0x8d, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, \
			0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a, 0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, \
			0x7d, 0xfa, 0xef, 0xc5, 0x91, 0x39, 0x72, 0xe4, 0xd3, 0xbd, 0x61, 0xc2, 0x9f, 0x25, 0x4a, 0x94, \
			0x33, 0x66, 0xcc, 0x83, 0x1d, 0x3a, 0x74, 0xe8, 0xcb, 0x8d, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, \
			0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a, 0x2f, 0x5e, 0xbc, 0x63, 0xc6, 0x97, 0x35, \
			0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91, 0x39, 0x72, 0xe4, 0xd3, 0xbd, 0x61, 0xc2, 0x9f, \
			0x25, 0x4a, 0x94, 0x33, 0x66, 0xcc, 0x83, 0x1d, 0x3a, 0x74, 0xe8, 0xcb, 0x8d, 0x01, 0x02, 0x04, \
			0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36, 0x6c, 0xd8, 0xab, 0x4d, 0x9a, 0x2f, 0x5e, 0xbc, 0x63, \
			0xc6, 0x97, 0x35, 0x6a, 0xd4, 0xb3, 0x7d, 0xfa, 0xef, 0xc5, 0x91, 0x39, 0x72, 0xe4, 0xd3, 0xbd, \
			0x61, 0xc2, 0x9f, 0x25, 0x4a, 0x94, 0x33, 0x66, 0xcc, 0x83, 0x1d, 0x3a, 0x74, 0xe8, 0xcb };


// Function Prototypes
void	AES_add_round_key(si4 round, ui1 state[][4], ui1 *round_key);
void	AES_decrypt(ui1 *in, ui1 *out, si1 *password, ui1 *expanded_key);
void	AES_encrypt(ui1 *in, ui1 *out, si1 *password, ui1 *expanded_key);
void	AES_key_expansion(ui1 *round_key, si1 *key);
void	AES_cipher(ui1 *in, ui1 *out, ui1 state[][4], ui1 *round_key);
si4	AES_get_sbox_invert(si4 num);
si4	AES_get_sbox_value(si4 num);
si4	*AES_initialize_rcon_table(si4 global_flag);
si4	*AES_initialize_rsbox_table(si4 global_flag);
si4	*AES_initialize_sbox_table(si4 global_flag);
void	AES_inv_cipher(ui1 *in, ui1 *out, ui1 state[][4], ui1 *round_key);
void	AES_inv_mix_columns(ui1 state[][4]);
void	AES_inv_shift_rows(ui1 state[][4]);
void	AES_inv_sub_bytes(ui1 state[][4]);
void	AES_mix_columns(ui1 state[][4]);
void	AES_shift_rows(ui1 state[][4]);
void	AES_sub_bytes(ui1 state[][4]);


/************************************************************************************/
/**************************************  SHA-2  *************************************/
/************************************************************************************/

// ATTRIBUTION
//
// FIPS 180-2 SHA-224/256/384/512 implementation
// Last update: 02/02/2007
// Issue date:  04/30/2005
//
// Copyright (C) 2005, 2007 Olivier Gay <olivier.gay@a3.epfl.ch>
// All rights reserved.
//
// "Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the project nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE PROJECT AND CONTRIBUTORS ``AS IS'' AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED.  IN NO EVENT SHALL THE PROJECT OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
// OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE."
//
// ONLY SHA-256 FUNCTIONS ARE INCLUDED IN THE MEF 3.0 LIBRARY
//
// Minor modifications for compatibility with the MEF Library.


// Constants
#define SHA256_OUTPUT_SIZE	256
#define SHA256_DIGEST_SIZE	(256 / 8)
#define SHA256_BLOCK_SIZE	(512 / 8)

#define SHA256_H0_ENTRIES	8
#define SHA256_H0 {	0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, \
			0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 }

#define SHA256_K_ENTRIES	64
#define SHA256_K {	0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, \
			0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, \
			0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, \
			0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174, \
			0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, \
			0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da, \
			0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, \
			0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, \
			0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, \
			0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, \
			0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, \
			0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, \
			0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, \
			0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3, \
			0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, \
			0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2 }

// Macros
#define SHFR(x, n)	(x >> n)
#define ROTR(x, n)	((x >> n) | (x << ((sizeof(x) << 3) - n)))
#define ROTL(x, n)	((x << n) | (x >> ((sizeof(x) << 3) - n)))
#define CH(x, y, z)	((x & y) ^ (~x & z))
#define MAJ(x, y, z)	((x & y) ^ (x & z) ^ (y & z))

#define SHA256_F1(x)	(ROTR(x,  2) ^ ROTR(x, 13) ^ ROTR(x, 22))
#define SHA256_F2(x)	(ROTR(x,  6) ^ ROTR(x, 11) ^ ROTR(x, 25))
#define SHA256_F3(x)	(ROTR(x,  7) ^ ROTR(x, 18) ^ SHFR(x,  3))
#define SHA256_F4(x)	(ROTR(x, 17) ^ ROTR(x, 19) ^ SHFR(x, 10))

#define UNPACK32(x, str)	{ *((str) + 3) = (ui1) (x); *((str) + 2) = (ui1) ((x) >>  8); *((str) + 1) = (ui1) ((x) >> 16); *((str) + 0) = (ui1) ((x) >> 24); }

#define PACK32(str, x)		{ *(x) = ((ui4) *((str) + 3)) | ((ui4) *((str) + 2) <<  8) | ((ui4) *((str) + 1) << 16) | ((ui4) *((str) + 0) << 24); }

#define SHA256_SCR(i)		{  w[i] =  SHA256_F4(w[i -  2]) + w[i -  7] + SHA256_F3(w[i - 15]) + w[i - 16]; }

#define SHA256_EXP(a, b, c, d, e, f, g, h, j)	{ t1 = wv[h] + SHA256_F2(wv[e]) + CH(wv[e], wv[f], wv[g]) + MEF_globals->SHA256_k_table[j] + w[j];	t2 = SHA256_F1(wv[a]) + MAJ(wv[a], wv[b], wv[c]); wv[d] += t1; wv[h] = t1 + t2; }

// Typedefs & Structures
typedef struct {
	ui4	tot_len;
	ui4	len;
	ui1	block[2 * SHA256_BLOCK_SIZE];
	ui4	h[8];
} SHA256_ctx;

// Function Prototypes
void	sha256(const ui1 *message, ui4 len, ui1 *digest);
void	SHA256_final(SHA256_ctx *ctx, ui1 *digest);
void	SHA256_init(SHA256_ctx *ctx);
ui4	*SHA256_initialize_h0_table(si4 global_flag);
ui4	*SHA256_initialize_k_table(si4 global_flag);
void	SHA256_transf(SHA256_ctx *ctx, const ui1 *message, ui4 block_nb);
void	SHA256_update(SHA256_ctx *ctx, const ui1 *message, ui4 len);



/************************************************************************************/
/******************  Library Includes (that depend on meflib.h)   *******************/
/************************************************************************************/

#include "mefrec.h"



#endif   // MEFLIB_IN







