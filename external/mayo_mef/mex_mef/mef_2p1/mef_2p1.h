/*									mef_header_2_1.h
*
* Specification for Mayo EEG Format (MEF) version 2.1, 
 # Copyright 2012, Mayo Foundation, Rochester MN. All rights reserved
 # Written by Ben Brinkmann, Matt Stead, and Dan Crepeau
 # usage and modification of this source code is governed by the Apache 2.0 license
 # Licensed under the Apache License, Version 2.0 (the "License");
 # you may not use this file except in compliance with the License.
 # You may obtain a copy of the License at
 # 
 #   http://www.apache.org/licenses/LICENSE-2.0
 # 
 # Unless required by applicable law or agreed to in writing, software
 # distributed under the License is distributed on an "AS IS" BASIS,
 # WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 # See the License for the specific language governing permissions and
 # limitations under the License.
 
 
 Thanks to all who acknowledge the Mayo Systems Electrophysiology Laboratory, Rochester, MN
 in academic publications of their work facilitated by this software.

 This file specifies the offsets to and sizes of (in bytes) all header values needed for generation of .mef files, as
 well as a structure for the header. Data types are specified in comments where applicable, shorthand notation as follows 
 (see size_types.h):
  signed char		si1
  unsigned char		ui1
  signed short		si2
  unsigned short	ui2
  signed int		si4
  unsigned int		ui4
  float			sf4
  long signed int	si8
  long unsigned int	ui8
  double		sf8
  long double		sf16
  n-char string		$(n)  -allow 1 space for termination character

 Header Encryption:
 The header begins with 176 unencrypted bytes, including two text fields and a series of numeric values defining the file’s 
 format and characteristics. The remainder of the header is encrypted with either a “subject” or “file” password. The passwords 
 are zero-terminated strings with a 15 character limit at this time. The subject password is used to access the file password 
 for decryption. The file password decrypts all but subject identifying header fields. The encryption / decryption algorithm 
 is the 128-bit AES standard (http://www.csrc.nist.gov/publications/fips/fips197/fips-197.pdf).

 Header Alignment:
 Fields in the header have required byte alignments relative to its start. 16-byte alignment facilitates encryption/decryption 
 beginning at that offset. Other alignment requirements are determined by the data-types: e.g. 8-byte alignment facilitates 
 reading si8, ui8, and sf8 data types
 
 Time Data:
 Each mef file ends in a block of recording time data. The offset to this data and the number of data entries are given in the file header. 
 This block contains triplets of times, file offsets, and sample indices of EEG data. Triplets are ui8 values containing the elapsed microseconds 
 since January 1, 1970 at 00:00:00 in the GMT (Greenwich, England) time zone. These values are easily converted to UTC time format (seconds since 
 1/1/1970 at 00:00:00 GMT), referred to hereafter as uUTC for "micro-UTC"
 
*/

#ifndef _MEF_H
#define _MEF_H

#ifndef SIZE_TYPES_IN
#define SIZE_TYPES_IN
	/* 64-bit typedefs - clearer nomenclature */
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
	typedef long double		sf16;
#endif

#define MEF_FALSE 0
#define MEF_TRUE 1


/************* header version & constants *******************/

#define HEADER_MAJOR_VERSION		2
#define HEADER_MINOR_VERSION		1
#define MEF_HEADER_LENGTH		1024
#define DATA_START_OFFSET		MEF_HEADER_LENGTH
#define UNENCRYPTED_REGION_OFFSET	0
#define UNENCRYPTED_REGION_LENGTH	176
#define SUBJECT_ENCRYPTION_OFFSET	176
#define SUBJECT_ENCRYPTION_LENGTH	160
#define SESSION_ENCRYPTION_OFFSET	352
#define SESSION_ENCRYPTION_LENGTH	512 //maintain multiple of 16
#define ENCRYPTION_BLOCK_BITS		128
#define ENCRYPTION_BLOCK_BYTES		(ENCRYPTION_BLOCK_BITS / 8)

#define SAMPLE_VALUE_NAN		(si4) 0xFF800000
#define SAMPLE_VALUE_NEG_INFINITY	(si4) 0xFF800001
#define SAMPLE_VALUE_POS_INFINITY	(si4) 0x007FFFFF

#define DMA_HIGH_FREQUENCY_FILTER	9000.0


/******************** header fields *************************/

// Begin Unencrypted Block
#define INSTITUTION_OFFSET			0
#define INSTITUTION_LENGTH			64		// $(63)
#define UNENCRYPTED_TEXT_FIELD_OFFSET		64
#define UNENCRYPTED_TEXT_FIELD_LENGTH		64		// $(63)
#define ENCRYPTION_ALGORITHM_OFFSET		128
#define ENCRYPTION_ALGORITHM_LENGTH		32		// $(29)
#define SUBJECT_ENCRYPTION_USED_OFFSET		160
#define SUBJECT_ENCRYPTION_USED_LENGTH		1		// ui1
#define SESSION_ENCRYPTION_USED_OFFSET		161
#define SESSION_ENCRYPTION_USED_LENGTH		1		// ui1
#define DATA_ENCRYPTION_USED_OFFSET		162
#define DATA_ENCRYPTION_USED_LENGTH		1		// ui1
#define BYTE_ORDER_CODE_OFFSET			163
#define BYTE_ORDER_CODE_LENGTH			1		// ui1
#define HEADER_MAJOR_VERSION_OFFSET		164
#define HEADER_MAJOR_VERSION_LENGTH		1		// ui1
#define HEADER_MINOR_VERSION_OFFSET		165
#define HEADER_MINOR_VERSION_LENGTH		1		// ui1
#define HEADER_LENGTH_OFFSET			166
#define HEADER_LENGTH_LENGTH			2		// ui2
#define SESSION_UNIQUE_ID_OFFSET		168
#define SESSION_UNIQUE_ID_LENGTH		8		// ui1
// End Unencrypted Block

// Begin Subject Encrypted Block
#define SUBJECT_FIRST_NAME_OFFSET		176
#define SUBJECT_FIRST_NAME_LENGTH		32		// $(31)
#define SUBJECT_SECOND_NAME_OFFSET		208
#define SUBJECT_SECOND_NAME_LENGTH		32		// $(31)
#define SUBJECT_THIRD_NAME_OFFSET		240
#define SUBJECT_THIRD_NAME_LENGTH		32		// $(31)
#define SUBJECT_ID_OFFSET			272
#define SUBJECT_ID_LENGTH			32		// $(31)
#define SESSION_PASSWORD_OFFSET			304
#define SESSION_PASSWORD_LENGTH			ENCRYPTION_BLOCK_BYTES		// $(15)
#define SUBJECT_VALIDATION_FIELD_OFFSET		320
#define SUBJECT_VALIDATION_FIELD_LENGTH		16
// End Subject Encrypted Block

// Begin Protected Block
#define PROTECTED_REGION_OFFSET			336
#define PROTECTED_REGION_LENGTH			16
// End Protected Block

// Begin Session Encrypted Block
#define SESSION_VALIDATION_FIELD_OFFSET		352
#define SESSION_VALIDATION_FIELD_LENGTH		16		// ui1
#define NUMBER_OF_SAMPLES_OFFSET		368
#define NUMBER_OF_SAMPLES_LENGTH		8		// ui8
#define CHANNEL_NAME_OFFSET			376
#define CHANNEL_NAME_LENGTH			32		// $(31)	
#define RECORDING_START_TIME_OFFSET		408
#define RECORDING_START_TIME_LENGTH		8		// ui8
#define RECORDING_END_TIME_OFFSET		416
#define RECORDING_END_TIME_LENGTH		8		// ui8
#define SAMPLING_FREQUENCY_OFFSET		424
#define SAMPLING_FREQUENCY_LENGTH		8		// sf8
#define LOW_FREQUENCY_FILTER_SETTING_OFFSET	432
#define LOW_FREQUENCY_FILTER_SETTING_LENGTH	8		// sf8
#define HIGH_FREQUENCY_FILTER_SETTING_OFFSET	440
#define HIGH_FREQUENCY_FILTER_SETTING_LENGTH	8		// sf8
#define NOTCH_FILTER_FREQUENCY_OFFSET		448
#define NOTCH_FILTER_FREQUENCY_LENGTH		8		// sf8
#define VOLTAGE_CONVERSION_FACTOR_OFFSET	456
#define VOLTAGE_CONVERSION_FACTOR_LENGTH	8		// sf8
#define ACQUISITION_SYSTEM_OFFSET		464
#define ACQUISITION_SYSTEM_LENGTH		32		// $(31)
#define CHANNEL_COMMENTS_OFFSET			496
#define CHANNEL_COMMENTS_LENGTH			128		// $(127)
#define STUDY_COMMENTS_OFFSET			624
#define STUDY_COMMENTS_LENGTH			128		// $(127)
#define PHYSICAL_CHANNEL_NUMBER_OFFSET		752
#define PHYSICAL_CHANNEL_NUMBER_LENGTH		4		// si4
#define COMPRESSION_ALGORITHM_OFFSET		756
#define COMPRESSION_ALGORITHM_LENGTH		32		// $(31)
#define MAXIMUM_COMPRESSED_BLOCK_SIZE_OFFSET	788
#define MAXIMUM_COMPRESSED_BLOCK_SIZE_LENGTH	4		// ui4
#define MAXIMUM_BLOCK_LENGTH_OFFSET		792
#define MAXIMUM_BLOCK_LENGTH_LENGTH		8		// ui8
#define BLOCK_INTERVAL_OFFSET			800
#define BLOCK_INTERVAL_LENGTH			8		// sf8
#define MAXIMUM_DATA_VALUE_OFFSET		808
#define MAXIMUM_DATA_VALUE_LENGTH		4		// si4
#define MINIMUM_DATA_VALUE_OFFSET		812
#define MINIMUM_DATA_VALUE_LENGTH		4		// si4
#define INDEX_DATA_OFFSET_OFFSET		816
#define	INDEX_DATA_OFFSET_LENGTH		8		// ui8
#define NUMBER_OF_INDEX_ENTRIES_OFFSET		824
#define NUMBER_OF_INDEX_ENTRIES_LENGTH		8		// ui8
#define BLOCK_HEADER_LENGTH_OFFSET		832
#define BLOCK_HEADER_LENGTH_LENGTH		2		// ui2
#define GMT_OFFSET_OFFSET			836
#define GMT_OFFSET_LENGTH			4		//sf4
#define DISCONTINUITY_DATA_OFFSET_OFFSET	840
#define DISCONTINUITY_DATA_OFFSET_LENGTH	8		//ui8
#define NUMBER_OF_DISCONTINUITY_ENTRIES_OFFSET	848
#define NUMBER_OF_DISCONTINUITY_ENTRIES_LENGTH	8		//ui8
// End Session Encrypted Block

// Begin Unencrypted Block
#define FILE_UNIQUE_ID_OFFSET			948
#define FILE_UNIQUE_ID_LENGTH			8
#define ANONYMIZED_SUBJECT_NAME_OFFSET		956
#define ANONYMIZED_SUBJECT_NAME_LENGTH		64		//$(63)
#define HEADER_CRC_OFFSET			1020
#define HEADER_CRC_LENGTH			4		//ui4
// End Unencrypted Block

/******************** structure & type definitions *****************/


typedef struct {
	ui8	time;
	ui8	file_offset;
	ui8	sample_number;
} INDEX_DATA;


typedef struct {
	si1		institution[INSTITUTION_LENGTH];
	si1		unencrypted_text_field[UNENCRYPTED_TEXT_FIELD_LENGTH];
	si1		encryption_algorithm[ENCRYPTION_ALGORITHM_LENGTH];
	ui1		subject_encryption_used;
	ui1		session_encryption_used;
	ui1		data_encryption_used;
	ui1		byte_order_code;
	ui1		header_version_major;
	ui1		header_version_minor;
	ui1		session_unique_ID[SESSION_UNIQUE_ID_LENGTH];
	ui2		header_length;
	si1		subject_first_name[SUBJECT_FIRST_NAME_LENGTH];
	si1		subject_second_name[SUBJECT_SECOND_NAME_LENGTH];
	si1		subject_third_name[SUBJECT_THIRD_NAME_LENGTH];
	si1		subject_id[SUBJECT_ID_LENGTH];
	si1		session_password[SESSION_PASSWORD_LENGTH];
	si1		subject_validation_field[SUBJECT_VALIDATION_FIELD_LENGTH];
	si1		session_validation_field[SESSION_VALIDATION_FIELD_LENGTH];
	si1		protected_region[PROTECTED_REGION_LENGTH];
	ui8		number_of_samples;
	si1		channel_name[CHANNEL_NAME_LENGTH];
	ui8		recording_start_time;
	ui8		recording_end_time;
	sf8		sampling_frequency;
	sf8		low_frequency_filter_setting;
	sf8		high_frequency_filter_setting;
	sf8		notch_filter_frequency;
	sf8		voltage_conversion_factor;
	si1		acquisition_system[ACQUISITION_SYSTEM_LENGTH];
	si1		channel_comments[CHANNEL_COMMENTS_LENGTH];
	si1		study_comments[STUDY_COMMENTS_LENGTH];
	si4		physical_channel_number;
	si1		compression_algorithm[COMPRESSION_ALGORITHM_LENGTH];
	ui4		maximum_compressed_block_size;
	ui8		maximum_block_length; 
	ui8		block_interval;
	si4		maximum_data_value;
	si4		minimum_data_value;
	ui8		index_data_offset;
	ui8		number_of_index_entries;
	ui2		block_header_length;
	sf4		GMT_offset;
	ui8		discontinuity_data_offset;
	ui8		number_of_discontinuity_entries;
	ui1		file_unique_ID[FILE_UNIQUE_ID_LENGTH];
	si1		anonymized_subject_name[ANONYMIZED_SUBJECT_NAME_LENGTH];
	ui4		header_crc;
	INDEX_DATA	*file_index;
	ui8		*discontinuity_data;
} MEF_HEADER_INFO;

//RED Codec
#define TOP_VALUE		(ui4) 0x80000000
#define TOP_VALUE_M_1		(ui4) 0x7FFFFFFF
#define CARRY_CHECK		(ui4) 0x7F800000
#define SHIFT_BITS		23
#define EXTRA_BITS		7
#define BOTTOM_VALUE		(ui4) 0x800000
#define BOTTOM_VALUE_M_1	(ui4) 0x7FFFFF
#define FILLER_BYTE		(ui1) 0x55 

//
/* 4 byte checksum, 4 byte compressed byte count, 8 byte time value, 4 byte difference count,  */
/* 4 byte sample count, 3 byte data maximum, 3 byte data minimum, 1 byte discontinuity flag, 256 byte model counts */
#define BLOCK_HEADER_BYTES			287
#define RED_CHECKSUM_OFFSET			0
#define RED_CHECKSUM_LENGTH			4
#define RED_COMPRESSED_BYTE_COUNT_OFFSET	4
#define RED_UUTC_TIME_OFFSET			8
#define RED_DIFFERENCE_COUNT_OFFSET		16
#define RED_SAMPLE_COUNT_OFFSET			20
#define RED_DATA_MAX_OFFSET			24
#define RED_DATA_MIN_OFFSET			27
#define RED_DISCONTINUITY_OFFSET		30
#define RED_STAT_MODEL_OFFSET			31
#define RED_DATA_OFFSET				BLOCK_HEADER_BYTES

/****************************************************************************************************/
/***  block size defines desired packet spacing - do not exceed 2^23 = 8388608 samples per block  ***/
/****************************************************************************************************/

typedef struct {
	ui4 CRC_32;
	ui1 CRC_validated;
	si4 compressed_bytes;
	ui8 block_start_time;
	si4 difference_count;
	si4 sample_count;
	si4 max_value; //NOTE: max and min are stored in block header as si3's
	si4 min_value;
	ui1 discontinuity;
} RED_BLOCK_HDR_INFO; 

// Globals
//static int crc_tab32_init = 0;//initialize to FALSE

//CRC polynomial and expanded key
#define	Koopman32	0xEB31D82E
#define CRC_KOOPMAN32_KEY {0x0, 0x09695c4ca, 0xfb4839c9, 0x6dddfd03, 0x20f3c3cf, 0xb6660705, 0xdbbbfa06, 0x4d2e3ecc, 0x41e7879e, 0xd7724354, 0xbaafbe57, 0x2c3a7a9d, 0x61144451, 0xf781809b, 0x9a5c7d98, 0xcc9b952, 0x83cf0f3c, 0x155acbf6, 0x788736f5, 0xee12f23f, 0xa33cccf3, 0x35a90839, 0x5874f53a, 0xcee131f0, 0xc22888a2, 0x54bd4c68, 0x3960b16b, 0xaff575a1, 0xe2db4b6d, 0x744e8fa7, 0x199372a4, 0x8f06b66e, 0xd1fdae25, 0x47686aef, 0x2ab597ec, 0xbc205326, 0xf10e6dea, 0x679ba920, 0xa465423, 0x9cd390e9, 0x901a29bb, 0x68fed71, 0x6b521072, 0xfdc7d4b8, 0xb0e9ea74, 0x267c2ebe, 0x4ba1d3bd, 0xdd341777, 0x5232a119, 0xc4a765d3, 0xa97a98d0, 0x3fef5c1a, 0x72c162d6, 0xe454a61c, 0x89895b1f, 0x1f1c9fd5, 0x13d52687, 0x8540e24d, 0xe89d1f4e, 0x7e08db84, 0x3326e548, 0xa5b32182, 0xc86edc81, 0x5efb184b, 0x7598ec17, 0xe30d28dd, 0x8ed0d5de, 0x18451114, 0x556b2fd8, 0xc3feeb12, 0xae231611, 0x38b6d2db, 0x347f6b89, 0xa2eaaf43, 0xcf375240, 0x59a2968a, 0x148ca846, 0x82196c8c, 0xefc4918f, 0x79515545, 0xf657e32b, 0x60c227e1, 0xd1fdae2, 0x9b8a1e28, 0xd6a420e4, 0x4031e42e, 0x2dec192d, 0xbb79dde7, 0xb7b064b5, 0x2125a07f, 0x4cf85d7c, 0xda6d99b6, 0x9743a77a, 0x1d663b0, 0x6c0b9eb3, 0xfa9e5a79, 0xa4654232, 0x32f086f8, 0x5f2d7bfb, 0xc9b8bf31, 0x849681fd, 0x12034537, 0x7fdeb834, 0xe94b7cfe, 0xe582c5ac, 0x73170166, 0x1ecafc65, 0x885f38af, 0xc5710663, 0x53e4c2a9, 0x3e393faa, 0xa8acfb60, 0x27aa4d0e, 0xb13f89c4, 0xdce274c7, 0x4a77b00d, 0x7598ec1, 0x91cc4a0b, 0xfc11b708, 0x6a8473c2, 0x664dca90, 0xf0d80e5a, 0x9d05f359, 0xb903793, 0x46be095f, 0xd02bcd95, 0xbdf63096, 0x2b63f45c, 0xeb31d82e, 0x7da41ce4, 0x1079e1e7, 0x86ec252d, 0xcbc21be1, 0x5d57df2b, 0x308a2228, 0xa61fe6e2, 0xaad65fb0, 0x3c439b7a, 0x519e6679, 0xc70ba2b3, 0x8a259c7f, 0x1cb058b5, 0x716da5b6, 0xe7f8617c, 0x68fed712, 0xfe6b13d8, 0x93b6eedb, 0x5232a11, 0x480d14dd, 0xde98d017, 0xb3452d14, 0x25d0e9de, 0x2919508c, 0xbf8c9446, 0xd2516945, 0x44c4ad8f, 0x9ea9343, 0x9f7f5789, 0xf2a2aa8a, 0x64376e40, 0x3acc760b, 0xac59b2c1, 0xc1844fc2, 0x57118b08, 0x1a3fb5c4, 0x8caa710e, 0xe1778c0d, 0x77e248c7, 0x7b2bf195, 0xedbe355f, 0x8063c85c, 0x16f60c96, 0x5bd8325a, 0xcd4df690, 0xa0900b93, 0x3605cf59, 0xb9037937, 0x2f96bdfd, 0x424b40fe, 0xd4de8434, 0x99f0baf8, 0xf657e32, 0x62b88331, 0xf42d47fb, 0xf8e4fea9, 0x6e713a63, 0x3acc760, 0x953903aa, 0xd8173d66, 0x4e82f9ac, 0x235f04af, 0xb5cac065, 0x9ea93439, 0x83cf0f3, 0x65e10df0, 0xf374c93a, 0xbe5af7f6, 0x28cf333c, 0x4512ce3f, 0xd3870af5, 0xdf4eb3a7, 0x49db776d, 0x24068a6e, 0xb2934ea4, 0xffbd7068, 0x6928b4a2, 0x4f549a1, 0x92608d6b, 0x1d663b05, 0x8bf3ffcf, 0xe62e02cc, 0x70bbc606, 0x3d95f8ca, 0xab003c00, 0xc6ddc103, 0x504805c9, 0x5c81bc9b, 0xca147851, 0xa7c98552, 0x315c4198, 0x7c727f54, 0xeae7bb9e, 0x873a469d, 0x11af8257, 0x4f549a1c, 0xd9c15ed6, 0xb41ca3d5, 0x2289671f, 0x6fa759d3, 0xf9329d19, 0x94ef601a, 0x27aa4d0, 0xeb31d82, 0x9826d948, 0xf5fb244b, 0x636ee081, 0x2e40de4d, 0xb8d51a87, 0xd508e784, 0x439d234e, 0xcc9b9520, 0x5a0e51ea, 0x37d3ace9, 0xa1466823, 0xec6856ef, 0x7afd9225, 0x17206f26, 0x81b5abec, 0x8d7c12be, 0x1be9d674, 0x76342b77, 0xe0a1efbd, 0xad8fd171, 0x3b1a15bb, 0x56c7e8b8, 0xc0522c72}
static ui4 crc_tab32[256] = CRC_KOOPMAN32_KEY;

// Range Encoding
typedef struct {
	ui4	low_bound;
	ui4	range;
	ui1	out_byte;
	ui4	underflow_bytes;
	ui1	*ob_p;
} RANGE_STATS;

// AES
#define AES_ENCRYPTION_KEY_LENGTH	240

// mef_lib function prototypes
si4		build_mef_header_block(ui1 *, MEF_HEADER_INFO *, si1 *);
si4		read_mef_header_block(ui1 *, MEF_HEADER_INFO *, si1 *);
ui4		calculate_header_CRC();
si4		validate_password(ui1 *, si1 *);
void		showHeader(MEF_HEADER_INFO *);
ui8		generate_unique_ID(ui1 *);
void		set_hdr_unique_ID(MEF_HEADER_INFO *, ui1 *);
void		set_block_hdr_unique_ID(ui1 *, ui1 *);
ui8		set_session_unique_ID(char *, ui1 *);
si4		check_header_block_alignment(ui1 *, si4);
void		strncpy2(si1 *, si1 *, si4);
void		init_hdr_struct(MEF_HEADER_INFO *);
si4		write_mef(si4 *, MEF_HEADER_INFO *, ui8, si1 *, si1 *);
si4		build_RED_block_header(ui1 *, RED_BLOCK_HDR_INFO *);
si4		read_RED_block_header(ui1 *, RED_BLOCK_HDR_INFO *);
ui4		calculate_compressed_block_CRC(ui1 *);
ui4		update_crc_32(ui4, si1);
void		init_crc32_tab(void);
ui8		RED_decompress_block(ui1 *, si4 *, si1 *, ui1 *, ui1, ui1,  RED_BLOCK_HDR_INFO *);
inline void	dec_normalize(ui4 *, ui4 *, ui1 *, ui1 **);
ui8		RED_compress_block(si4 *, ui1 *, ui4, ui8, ui1, ui1 *, ui1, RED_BLOCK_HDR_INFO *);
void		done_encoding(RANGE_STATS *);
inline void	encode_symbol(ui1, ui4, ui4, ui4, RANGE_STATS *);
inline void	enc_normalize(RANGE_STATS *);
ui1		cpu_endianness();
void		reverse_in_place(void *, si4);
si2		rev_si2(si2);
ui2		rev_ui2(ui2);
ui8		rev_ui8(ui8);
sf8		rev_sf8(sf8);
si4		rev_si4(si4);
ui4		rev_ui4(ui4);
sf4		rev_sf4(sf4);
void		AES_decryptWithKey(ui1 *, ui1 *, ui1 *);
void		AES_decrypt(ui1 *, ui1 *, si1 *);
void		AES_encryptWithKey(ui1 *, ui1 *, ui1 *);
void		AES_encrypt(ui1 *, ui1 *, si1 *);
void		InvCipher(si4, ui1 *, ui1 *, ui1 [][4], ui1 *);
void		Cipher(si4, ui1 *, ui1 *, ui1 [][4], ui1 *);
void		InvMixColumns(ui1 [][4]);
void		MixColumns(ui1 [][4]);
void		InvShiftRows(ui1 [][4]);
void		ShiftRows(ui1 [][4]);
void		SubBytes(ui1 [][4]);
void		AddRoundKey(si4, ui1 [][4], ui1*);
void		AES_KeyExpansion(si4, si4, ui1 *, si1 *);
si4		getSBoxInvert(si4);
si4		getSBoxValue(si4);

#endif



