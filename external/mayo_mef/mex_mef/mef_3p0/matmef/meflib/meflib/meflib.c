

/************************************************************************************/
/**********************************  MEF 3.0 C Library  *****************************/
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
// set editor preferences to these for intended alignment


#include "meflib.h"

// global
MEF_GLOBALS	*MEF_globals = NULL;

#ifdef _WIN32
	void bzero(void *dest, size_t num)
	{
		memset(dest, 0, num);
	}

	int random()
	{
		return rand();
	}

	void srandom(int num)
	{
		srand(num);
	}

	#if !HAVE_TIME_R
	struct tm *gmtime_r(time_t *_clock, struct tm *_result)
	{
		struct tm *p;

	    p = gmtime(_clock);

		if (p)
			*(_result) = *p;

		return p;
	}

	struct tm *localtime_r(time_t *_clock, struct tm *_result)
	{
		struct tm *p;

	    p = localtime(_clock);

		if (p)
			*(_result) = *p;

		return p;
	}
	#endif // !HAVE_TIME_R

	char *asctime_r(struct tm *timeptr, char *buf)
	{
		strcpy(buf, asctime(timeptr));

		return buf;
	}

	time_t timegm(struct tm * a_tm)
	{
	    int offset;
		time_t ltime;
	    time_t utc;
	    struct tm tm_val;

	    ltime = mktime(a_tm);
		gmtime_s(&tm_val, &ltime);
		offset = (tm_val.tm_hour - a_tm->tm_hour);
		if (offset > 12)
		{
			offset = 24 - offset;
		}
		utc = mktime(a_tm) - offset * 3600;
		return utc;
	}
#endif


/*************************************************************************/
/****************************  AES-128 FUNCTIONS  ************************/
/*************************************************************************/


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
// Minor modifications have been made for compatibility with the MEF Library.


// This function adds the round key to state.
// The round key is added to the state by an XOR function.
void	AES_add_round_key(si4 round, ui1 state[][4], ui1 *round_key)
{
	si4	i, j;
	
	
	for (i = 0; i < 4; i++) {
		for (j = 0;j < 4; j++) {
			state[j][i] ^= round_key[round * AES_NB * 4 + i * AES_NB + j];
		}
	}

	
	return;
}


// In is encrypted buffer (16 bytes / 128 bits).
// Out is decrypted buffer (16 bytes / 128 bits).
// In can equal out, i.e. can be done in place.
// Pass in expanded key externally - this is more efficient than passing the password
// if encrypting multiple times with the same encryption key
void	AES_decrypt(ui1 *in, ui1 *out, si1 *password, ui1 *expanded_key)
{
	si1	key[16] = {0};
	ui1	state[4][4]; // the array that holds the intermediate results during encryption
	ui1	round_key[240]; // The array that stores the round keys
	
	
	if (expanded_key != NULL) {
		AES_inv_cipher(in, out, state, expanded_key);
	} else if (password != NULL) {
		// password becomes the key (16 bytes, zero-padded if shorter, truncated if longer)
		strncpy((si1 *) key, password, 16);
		
		//The Key-Expansion routine must be called before the decryption routine.
		AES_key_expansion(round_key, key);
		
		// The next function call decrypts the CipherText with the Key using AES algorithm.
		AES_inv_cipher(in, out, state, round_key);
	} else {
		fprintf(stderr, "Error: No password or expanded key => exiting [function \"%s\", line %d]\n", __FUNCTION__, __LINE__);
		exit(-1);
	}
	
        
	return;
}


// in is buffer to be encrypted (16 bytes)
// out is encrypted buffer (16 bytes)
// in can equal out, i.e. can be done in place
// Pass in expanded key externally - this is more efficient tahn passing the password
// if encrypting multiple times with the same encryption key
void	AES_encrypt(ui1 *in, ui1 *out, si1 *password, ui1 *expanded_key)
{
	si1	key[16] = {0};
	ui1	state[4][4]; // the array that holds the intermediate results during encryption
	ui1	round_key[240]; // The array that stores the round keys
	
	
	if (expanded_key != NULL) {
		AES_cipher(in, out, state, expanded_key);
	} else if (password != NULL) {
		// password becomes the key (16 bytes, zero-padded if shorter, truncated if longer)
		strncpy((si1 *) key, password, 16);
		
		// The KeyExpansion routine must be called before encryption.
		AES_key_expansion(round_key, key);
		
		// The next function call encrypts the PlainText with the Key using AES algorithm.
		AES_cipher(in, out, state, round_key);
	} else {
		fprintf(stderr, "Error: No password or expanded key => exiting [function \"%s\", line %d]\n", __FUNCTION__, __LINE__);
		exit(-1);
	}
	
        
	return;
}


// This function produces AES_NB * (AES_NR + 1) round keys. The round keys are used in each round to encrypt the states.
// NOTE: make sure any terminal unused bytes in key array (password) are zeroed
void	AES_key_expansion(ui1 *expanded_key, si1 *key)
{
	// The round constant word array, Rcon[i], contains the values given by
	// x to the power (i - 1) being powers of x (x is denoted as {02}) in the field GF(28)
	// Note that i starts at 1, not 0).
	si4	i, j;
	ui1	temp[4], k;
        
	
	if (MEF_globals->AES_rcon_table == NULL)
		(void) AES_initialize_rcon_table(MEF_TRUE);
	
	// The first round key is the key itself.
	for (i = 0; i < AES_NK; i++) {
		expanded_key[i * 4] = key[i * 4];
		expanded_key[i * 4 + 1] = key[i * 4 + 1];
		expanded_key[i * 4 + 2] = key[i * 4 + 2];
		expanded_key[i * 4 + 3] = key[i * 4 + 3];
	}
	
	// All other round keys are found from the previous round keys.
	while (i < (AES_NB * (AES_NR + 1))) {
		
		for (j = 0; j < 4; j++) {
			temp[j] = expanded_key[(i - 1) * 4 + j];
		}
		
		if (i % AES_NK == 0) {
			// This rotates the 4 bytes in a word to the left once.
			// [a0,a1,a2,a3] becomes [a1,a2,a3,a0]
			k = temp[0];
			temp[0] = temp[1];
			temp[1] = temp[2];
			temp[2] = temp[3];
			temp[3] = k;
			
			// This takes a four-byte input word and applies the S-box
			// to each of the four bytes to produce an output word.
			temp[0] = AES_get_sbox_value(temp[0]);
			temp[1] = AES_get_sbox_value(temp[1]);
			temp[2] = AES_get_sbox_value(temp[2]);
			temp[3] = AES_get_sbox_value(temp[3]);
			
			temp[0] = temp[0] ^ MEF_globals->AES_rcon_table[i / AES_NK];
		} else if (AES_NK > 6 && i % AES_NK == 4) {
			// This takes a four-byte input word and applies the S-box
			// to each of the four bytes to produce an output word.
			temp[0] = AES_get_sbox_value(temp[0]);
			temp[1] = AES_get_sbox_value(temp[1]);
			temp[2] = AES_get_sbox_value(temp[2]);
			temp[3] = AES_get_sbox_value(temp[3]);
		}
		
		expanded_key[i * 4] = expanded_key[(i - AES_NK) * 4] ^ temp[0];
		expanded_key[i * 4 + 1] = expanded_key[(i - AES_NK) * 4 + 1] ^ temp[1];
		expanded_key[i * 4 + 2] = expanded_key[(i - AES_NK) * 4 + 2] ^ temp[2];
		expanded_key[i * 4 + 3] = expanded_key[(i - AES_NK) * 4 + 3] ^ temp[3];
		
		i++;
	}
	
        
	return;
}


// Cipher is the main function that encrypts the PlainText.
void	AES_cipher(ui1 *in, ui1 *out, ui1 state[][4], ui1 *round_key)
{
	si4	i, j, round = 0;
	
	
	//Copy the input PlainText to state array.
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			state[j][i] = in[i * 4 + j];
		}
	}
	
	// Add the First round key to the state before starting the rounds.
	AES_add_round_key(0, state, round_key);
	
	// There will be AES_NR rounds.
	// The first AES_NR - 1 rounds are identical.
	// These AES_NR - 1 rounds are executed in the loop below.
	for (round = 1; round < AES_NR; round++) {
		AES_sub_bytes(state);
		AES_shift_rows(state);
		AES_mix_columns(state);
		AES_add_round_key(round, state, round_key);
	}
	
	// The last round is given below.
	// The MixColumns function is not here in the last round.
	AES_sub_bytes(state);
	AES_shift_rows(state);
	AES_add_round_key(AES_NR, state, round_key);
	
	// The encryption process is over.
	// Copy the state array to output array.
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			out[i * 4 + j] = state[j][i];
		}
	}
	
        
	return;
}


inline si4	AES_get_sbox_invert(si4 num)
{
	if (MEF_globals->AES_rsbox_table == NULL)
		(void) AES_initialize_rsbox_table(MEF_TRUE);
	
        
	return(MEF_globals->AES_rsbox_table[num]);
}


inline si4	AES_get_sbox_value(si4 num)
{
	if (MEF_globals->AES_sbox_table == NULL)
		(void) AES_initialize_sbox_table(MEF_TRUE);
	
        
	return(MEF_globals->AES_sbox_table[num]);
}


si4	*AES_initialize_rcon_table(si4 global_flag)
{
	si4	*rcon_table;
	
	
	rcon_table = (si4 *) e_calloc((size_t) AES_RCON_ENTRIES, sizeof(si4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	{
		si4 temp[AES_RCON_ENTRIES] = AES_RCON;
		memcpy(rcon_table, temp, AES_RCON_ENTRIES * sizeof(si4));
	}
	
	if (global_flag == MEF_TRUE) {
		MEF_globals->AES_rcon_table = rcon_table;
		return(NULL);
	}
	
        
	return(rcon_table);
}


si4	*AES_initialize_rsbox_table(si4 global_flag)
{
	si4	*rsbox_table;
	
	
	rsbox_table = (si4 *) e_calloc((size_t) AES_RSBOX_ENTRIES, sizeof(si4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	{
		si4 temp[AES_RSBOX_ENTRIES] = AES_RSBOX;
		memcpy(rsbox_table, temp, AES_RSBOX_ENTRIES * sizeof(si4));
	}
	
	if (global_flag == MEF_TRUE) {
		MEF_globals->AES_rsbox_table = rsbox_table;
		return(NULL);
	}
	
        
	return(rsbox_table);
}


si4	*AES_initialize_sbox_table(si4 global_flag)
{
	si4	*sbox_table;
	
	
	sbox_table = (si4 *) e_calloc((size_t) AES_SBOX_ENTRIES, sizeof(si4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	{
		si4 temp[AES_SBOX_ENTRIES] = AES_SBOX;
		memcpy(sbox_table, temp, AES_SBOX_ENTRIES * sizeof(si4));
	}
	
	if (global_flag == MEF_TRUE) {
		MEF_globals->AES_sbox_table = sbox_table;
		return(NULL);
	}
	
        
	return(sbox_table);
}


// AES_inv_cipher is the main decryption function
void	AES_inv_cipher(ui1 *in, ui1 *out, ui1 state[][4], ui1 *round_key)
{
	si4	i, j, round = 0;
	
	
	// Copy the input encrypted text to state array.
	for (i = 0; i < 4; i++) {
		for (j = 0;j < 4; j++) {
			state[j][i] = in[i * 4 + j];
		}
	}
	
	// Add the First round key to the state before starting the rounds.
	AES_add_round_key(AES_NR, state, round_key);
	
	// There will be AES_NR rounds.
	// The first AES_NR - 1 rounds are identical.
	// These AES_NR - 1 rounds are executed in the loop below.
	for (round = AES_NR - 1; round > 0; round--) {
		AES_inv_shift_rows(state);
		AES_inv_sub_bytes(state);
		AES_add_round_key(round, state, round_key);
		AES_inv_mix_columns(state);
	}
	
	// The last round is given below.
	// The MixColumns function is not here in the last round.
	AES_inv_shift_rows(state);
	AES_inv_sub_bytes(state);
	AES_add_round_key(0, state, round_key);
	
	// The decryption process is over.
	// Copy the state array to output array.
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			out[i * 4 + j] = state[j][i];
		}
	}
	
        
	return;
}


// The method used to multiply may be difficult to understand.
// Please use the references to gain more information.
void	AES_inv_mix_columns(ui1 state[][4])
{
	si4	i;
	ui1	a, b, c, d;
	
	
	for (i = 0; i < 4; i++) {
		a = state[0][i];
		b = state[1][i];
		c = state[2][i];
		d = state[3][i];
		state[0][i] = AES_MULTIPLY(a, 0x0e) ^ AES_MULTIPLY(b, 0x0b) ^ AES_MULTIPLY(c, 0x0d) ^ AES_MULTIPLY(d, 0x09);
		state[1][i] = AES_MULTIPLY(a, 0x09) ^ AES_MULTIPLY(b, 0x0e) ^ AES_MULTIPLY(c, 0x0b) ^ AES_MULTIPLY(d, 0x0d);
		state[2][i] = AES_MULTIPLY(a, 0x0d) ^ AES_MULTIPLY(b, 0x09) ^ AES_MULTIPLY(c, 0x0e) ^ AES_MULTIPLY(d, 0x0b);
		state[3][i] = AES_MULTIPLY(a, 0x0b) ^ AES_MULTIPLY(b, 0x0d) ^ AES_MULTIPLY(c, 0x09) ^ AES_MULTIPLY(d, 0x0e);
	}
	
        
	return;
}


void	AES_inv_shift_rows(ui1 state[][4])
{
	ui1	temp;
	
	
	// Rotate first row 1 columns to right
	temp = state[1][3];
	state[1][3] = state[1][2];
	state[1][2] = state[1][1];
	state[1][1] = state[1][0];
	state[1][0] = temp;
	
	// Rotate second row 2 columns to right
	temp = state[2][0];
	state[2][0] = state[2][2];
	state[2][2] = temp;
	
	temp = state[2][1];
	state[2][1] = state[2][3];
	state[2][3] = temp;
	
	// Rotate third row 3 columns to right
	temp = state[3][0];
	state[3][0] = state[3][1];
	state[3][1] = state[3][2];
	state[3][2] = state[3][3];
	state[3][3] = temp;
	
        
	return;
}


void	AES_inv_sub_bytes(ui1 state[][4])
{
	si4	i, j;
	
	
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			state[i][j] = AES_get_sbox_invert(state[i][j]);
		}
	}
	
        
	return;
}


// MixColumns function mixes the columns of the state matrix
// The method used may look complicated, but it is easy if you know the underlying theory.
// Refer the documents specified above.
void	AES_mix_columns(ui1 state[][4])
{
	si4	i;
	ui1	Tmp, Tm, t;
	
	
	for (i = 0; i < 4; i++) {
		t = state[0][i];
		Tmp = state[0][i] ^ state[1][i] ^ state[2][i] ^ state[3][i];
		Tm = state[0][i] ^ state[1][i];
		Tm = AES_XTIME(Tm);
		state[0][i] ^= Tm ^ Tmp;
		Tm = state[1][i] ^ state[2][i];
		Tm = AES_XTIME(Tm);
		state[1][i] ^= Tm ^ Tmp;
		Tm = state[2][i] ^ state[3][i];
		Tm = AES_XTIME(Tm);
		state[2][i] ^= Tm ^ Tmp;
		Tm = state[3][i] ^ t;
		Tm = AES_XTIME(Tm);
		state[3][i] ^= Tm ^ Tmp;
	}
	
        
	return;
}


// The ShiftRows() function shifts the rows in the state to the left.
// Each row is shifted with different offset.
// Offset = Row number. So the first row is not shifted.
void	AES_shift_rows(ui1 state[][4])
{
	ui1	temp;
	
	
	// Rotate first row 1 columns to left
	temp = state[1][0];
	state[1][0] = state[1][1];
	state[1][1] = state[1][2];
	state[1][2] = state[1][3];
	state[1][3] = temp;
	
	// Rotate second row 2 columns to left
	temp = state[2][0];
	state[2][0] = state[2][2];
	state[2][2] = temp;
	
	temp = state[2][1];
	state[2][1] = state[2][3];
	state[2][3] = temp;
	
	// Rotate third row 3 columns to left
	temp = state[3][0];
	state[3][0] = state[3][3];
	state[3][3] = state[3][2];
	state[3][2] = state[3][1];
	state[3][1] = temp;
	
        
	return;
}


// The SubBytes Function Substitutes the values in the
// state matrix with values in an S-box.
void	AES_sub_bytes(ui1 state[][4])
{
	si4	i, j;
	
	
	for (i = 0; i < 4; i++) {
		for (j = 0; j < 4; j++) {
			state[i][j] = AES_get_sbox_value(state[i][j]);
		}
	}
	
        
	return;
}


/*************************************************************************/
/**************************  End AES-128 FUNCTIONS  **********************/
/*************************************************************************/


si1	all_zeros(ui1 *bytes, si4 field_length)
{
	while (field_length--)
		if (*bytes++)
			return(MEF_FALSE);
	
        
	return(MEF_TRUE);
}


FILE_PROCESSING_STRUCT	*allocate_file_processing_struct(si8 raw_data_bytes, ui4 file_type_code, FILE_PROCESSING_DIRECTIVES *directives, FILE_PROCESSING_STRUCT *proto_fps, si8 bytes_to_copy)
{
        FILE_PROCESSING_STRUCT	*fps;
        void			*data_ptr;
	
	// allocate
        fps = (FILE_PROCESSING_STRUCT *) e_calloc((size_t) 1, sizeof(FILE_PROCESSING_STRUCT), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	// zero metadata structure because calloc() does not zero substructues and arrays
	fps->metadata.section_1 = NULL; fps->metadata.time_series_section_2 = NULL; fps->metadata.video_section_2 = NULL; fps->metadata.section_3 = NULL;
	
	if (raw_data_bytes > 0) {
        	fps->raw_data = (ui1 *) e_calloc((size_t) raw_data_bytes, sizeof(ui1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		if (raw_data_bytes >= UNIVERSAL_HEADER_BYTES)
			fps->universal_header = (UNIVERSAL_HEADER *) fps->raw_data; // all files start with universal header
        }
        fps->raw_data_bytes = raw_data_bytes;
        fps->file_type_code = file_type_code;
	fps->file_length = FPS_FILE_LENGTH_UNKNOWN;
	
        // if a prototype FILE_PROCESSING_STRUCT was passed - copy its directives, password data, and raw data
        if (proto_fps != NULL) {
		if (directives == NULL)
			fps->directives = proto_fps->directives;
                fps->password_data = proto_fps->password_data;
                if ((bytes_to_copy > proto_fps->raw_data_bytes) || (bytes_to_copy > fps->raw_data_bytes)) {
                        fprintf(stderr, "Error: copy request size exceeds avaiable data => no copying done [function \"%s\", line %d]\n", __FUNCTION__, __LINE__);
                        if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL) {
                                (void) fprintf(stderr, "\t=> exiting program\n\n");
        			exit(1);
                        }
		} else
                	memcpy(fps->raw_data, proto_fps->raw_data, bytes_to_copy);
        }
	
	// explicitly passed directives supersede the prototype's directives
	if (directives != NULL)
		fps->directives = *directives;
	else if (proto_fps == NULL)
		(void) initialize_file_processing_directives(&fps->directives);  // set directives to defaults
	
	// set up the universal header
	if (fps->universal_header != NULL) {
		if (proto_fps == NULL)
			initialize_universal_header(fps, MEF_FALSE, MEF_FALSE, MEF_FALSE);
		// passed file type code supersedes copied on
		MEF_strncpy(fps->universal_header->file_type_string, (si1 *) &fps->file_type_code, TYPE_BYTES);
		// set CRCs to no entry values
		fps->universal_header->header_CRC = fps->universal_header->body_CRC = CRC_NO_ENTRY;
	}
        
        // set appropriate pointers
        data_ptr = NULL;
		if (fps->raw_data != NULL)
			data_ptr = (ui1 *)fps->raw_data + UNIVERSAL_HEADER_BYTES;
        switch (file_type_code) {
                case NO_TYPE_CODE:
                        break;
                case TIME_SERIES_INDICES_FILE_TYPE_CODE:
                        fps->time_series_indices = (TIME_SERIES_INDEX *) data_ptr;
                        break;
                case TIME_SERIES_DATA_FILE_TYPE_CODE:
                        fps->RED_blocks = (ui1 *) data_ptr;
                        break;
                case TIME_SERIES_METADATA_FILE_TYPE_CODE:
                        fps->metadata.section_1 = (METADATA_SECTION_1 *) data_ptr;
                        fps->metadata.time_series_section_2 = (TIME_SERIES_METADATA_SECTION_2 *) (fps->raw_data + METADATA_SECTION_2_OFFSET);
                        fps->metadata.section_3 = (METADATA_SECTION_3 *) (fps->raw_data + METADATA_SECTION_3_OFFSET);
			if (raw_data_bytes == METADATA_FILE_BYTES && bytes_to_copy != METADATA_FILE_BYTES)
				initialize_metadata(fps);
                        break;
                case VIDEO_METADATA_FILE_TYPE_CODE:
                        fps->metadata.section_1 = (METADATA_SECTION_1 *) data_ptr;
                        fps->metadata.video_section_2 = (VIDEO_METADATA_SECTION_2 *) (fps->raw_data + METADATA_SECTION_2_OFFSET);
                        fps->metadata.section_3 = (METADATA_SECTION_3 *) (fps->raw_data + METADATA_SECTION_3_OFFSET);
			if (raw_data_bytes == METADATA_FILE_BYTES && bytes_to_copy != METADATA_FILE_BYTES)
				initialize_metadata(fps);
                        break;
                case VIDEO_INDICES_FILE_TYPE_CODE:
                        fps->video_indices = (VIDEO_INDEX *) data_ptr;
                	break;
                case RECORD_DATA_FILE_TYPE_CODE:
                        fps->records = (ui1 *) data_ptr;
                        break;
                case RECORD_INDICES_FILE_TYPE_CODE:
                        fps->record_indices = (RECORD_INDEX *) data_ptr;
                        break;
                default:
                        free(fps->raw_data);
                        free(fps);
                        fprintf(stderr, "Error: unrecognized type code \"0x%x\" [function \"%s\", line %d]\n", file_type_code, __FUNCTION__, __LINE__);
                        if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL) {
                                (void) fprintf(stderr, "\t=> exiting program\n\n");
        			exit(1);
                        }
                        return(NULL);
        }
        
        
        return(fps);
}


inline void	apply_recording_time_offset(si8 *time)
{
        if (*time == UUTC_NO_ENTRY)
                return;
        
	if (*time < 0)  // negative times indicate recording time offset already applied
		return;
	
	// apply recording time offset & make negative to indicate application
	*time = -(*time - MEF_globals->recording_time_offset);
	
        
	return;
}


si4	channel_type_from_path(si1 *path)
{
	si1	*c, temp_path[MEF_FULL_FILE_NAME_BYTES], name[MEF_SEGMENT_BASE_FILE_NAME_BYTES], extension[TYPE_BYTES];
	
	
	// move pointer to end of string
	c = path + strlen(path) - 1;
	
	// ignore terminal "/" if present
	if (*c == '/')
		c--;
	
	// copy extension
	c -= 4;
	if (*c != '.')
		return(UNKNOWN_CHANNEL_TYPE);
	MEF_strncpy(extension, ++c, TYPE_BYTES);
	
	// compare extension: record types => get extension of next level up the hierarchy
	if (!(strcmp(extension, RECORD_DATA_FILE_TYPE_STRING)) || !(strcmp(extension, RECORD_INDICES_FILE_TYPE_STRING))) {
		extract_path_parts(path, temp_path, NULL, NULL);
		extract_path_parts(temp_path, temp_path, name, extension);
	}
	
	// compare extension: segment directory => get extension of next level up the hierarchy
	else if (!(strcmp(extension, SEGMENT_DIRECTORY_TYPE_STRING))) {
		extract_path_parts(path, temp_path, NULL, NULL);
		extract_path_parts(temp_path, temp_path, name, extension);
	}
	
	// compare extension: TIMES_SERIES_CHANNEL_TYPE
	if (!(strcmp(extension, TIME_SERIES_CHANNEL_DIRECTORY_TYPE_STRING)))
		return(TIME_SERIES_CHANNEL_TYPE);
	else if (!(strcmp(extension, TIME_SERIES_METADATA_FILE_TYPE_STRING)))
		return(TIME_SERIES_CHANNEL_TYPE);
	else if (!(strcmp(extension, TIME_SERIES_DATA_FILE_TYPE_STRING)))
		return(TIME_SERIES_CHANNEL_TYPE);
	else if (!(strcmp(extension, TIME_SERIES_INDICES_FILE_TYPE_STRING)))
		return(TIME_SERIES_CHANNEL_TYPE);
	
	// compare extension: VIDEO_CHANNEL_TYPE
	else if (!(strcmp(extension, VIDEO_CHANNEL_DIRECTORY_TYPE_STRING)))
		return(VIDEO_CHANNEL_TYPE);
	else if (!(strcmp(extension, VIDEO_METADATA_FILE_TYPE_STRING)))
		return(VIDEO_CHANNEL_TYPE);
	else if (!(strcmp(extension, VIDEO_INDICES_FILE_TYPE_STRING)))
		return(VIDEO_CHANNEL_TYPE);

	// unknown channel type
	return(UNKNOWN_CHANNEL_TYPE);
}


/*************************************************************************/
/**********************  ALIGNMENT CHECKING FUNCTIONS  *******************/
/*************************************************************************/


si4	check_all_alignments(const si1 *function, si4 line)
{
	si4		return_value;
	ui1		*bytes;
	
        
	// see if already checked
	if (MEF_globals->all_structures_aligned != MEF_UNKNOWN)
		return(MEF_globals->all_structures_aligned);
        
	return_value = MEF_TRUE;
	bytes = (ui1 *) e_malloc(METADATA_FILE_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);  // METADATA is largest file structure
	
	// check all structures
	if ((check_universal_header_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_metadata_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_RED_block_header_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_time_series_indices_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_video_indices_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_record_indices_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_record_header_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_record_structure_alignments(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	
	free(bytes);
	
	if (return_value == MEF_TRUE) {
		MEF_globals->all_structures_aligned = MEF_TRUE;
		if (MEF_globals->verbose == MEF_TRUE)
			(void) printf("%s(): All MEF Library structures are aligned\n", __FUNCTION__);
	} else {
		if (!(MEF_globals->behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
			(void) fprintf(stderr, "%c\n%s(): unaligned MEF structures (code update required)\n", 7, __FUNCTION__);
			(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
			if (MEF_globals->behavior_on_fail & RETURN_ON_FAIL)
				(void) fprintf(stderr, "\t=> returning MEF_FALSE\n\n");
			else if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL)
				(void) fprintf(stderr, "\t=> exiting program\n\n");
		}
		if (MEF_globals->behavior_on_fail & RETURN_ON_FAIL)
			return(MEF_FALSE);
		else if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL)
			exit(1);
	}
	
        
	return(return_value);
}


si4	check_metadata_alignment(ui1 *bytes)
{
	si4		return_value;
	si4		free_flag = MEF_FALSE;
	
	
	// see if already checked
	if (MEF_globals->all_metadata_structures_aligned != MEF_UNKNOWN)
		return(MEF_globals->all_metadata_structures_aligned);
	
	return_value = MEF_TRUE;
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(METADATA_FILE_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	
	// check substructures
	if ((check_metadata_section_1_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_time_series_metadata_section_2_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_video_metadata_section_2_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_metadata_section_3_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (return_value == MEF_TRUE)
		MEF_globals->all_metadata_structures_aligned = MEF_TRUE;
	
        
	return(return_value);
}


si4	check_metadata_section_1_alignment(ui1 *bytes)
{
	METADATA_SECTION_1	*md1;
	si4			free_flag = MEF_FALSE;
	
	
	// see if already checked
	if (MEF_globals->metadata_section_1_aligned == MEF_UNKNOWN)
		MEF_globals->metadata_section_1_aligned = MEF_FALSE;
	else
		return(MEF_globals->metadata_section_1_aligned);
	
	// check overall size
	if (sizeof(METADATA_SECTION_1) != METADATA_SECTION_1_BYTES)
		goto METADATA_SECTION_1_NOT_ALIGNED;
	
	// check fields
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(METADATA_FILE_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	md1 = (METADATA_SECTION_1 *) (bytes + UNIVERSAL_HEADER_BYTES);
	if (&md1->section_2_encryption != (si1 *) (bytes + METADATA_SECTION_2_ENCRYPTION_OFFSET))
		goto METADATA_SECTION_1_NOT_ALIGNED;
	if (&md1->section_3_encryption != (si1 *) (bytes + METADATA_SECTION_3_ENCRYPTION_OFFSET))
		goto METADATA_SECTION_1_NOT_ALIGNED;
	if (md1->protected_region != (ui1 *) (bytes + METADATA_SECTION_1_PROTECTED_REGION_OFFSET))
		goto METADATA_SECTION_1_NOT_ALIGNED;
	if (md1->discretionary_region != (ui1 *) (bytes + METADATA_SECTION_1_DISCRETIONARY_REGION_OFFSET))
		goto METADATA_SECTION_1_NOT_ALIGNED;
	
	// aligned
        MEF_globals->metadata_section_1_aligned = MEF_TRUE;
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): METADATA_SECTION_1 structure is aligned\n", __FUNCTION__);
	
        return(MEF_TRUE);
	
	// not aligned
METADATA_SECTION_1_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) fprintf(stderr, "%c%s(): METADATA_SECTION_1 structure is not aligned\n", 7, __FUNCTION__);
	
        
        return(MEF_FALSE);
}


si4	check_metadata_section_3_alignment(ui1 *bytes)
{
	METADATA_SECTION_3	*md3;
	si4			free_flag = MEF_FALSE;
	
	
	// see if already checked
	if (MEF_globals->metadata_section_3_aligned == MEF_UNKNOWN)
		MEF_globals->metadata_section_3_aligned = MEF_FALSE;
	else
		return(MEF_globals->metadata_section_3_aligned);
	
	// check overall size
	if (sizeof(METADATA_SECTION_3) != METADATA_SECTION_3_BYTES)
		goto METADATA_SECTION_3_NOT_ALIGNED;
	
	// check fields
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(METADATA_FILE_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	md3 = (METADATA_SECTION_3 *) (bytes + METADATA_SECTION_3_OFFSET);
	if (&md3->recording_time_offset != (si8 *) (bytes + METADATA_RECORDING_TIME_OFFSET_OFFSET))
		goto METADATA_SECTION_3_NOT_ALIGNED;
	if (&md3->DST_start_time != (si8 *) (bytes + METADATA_DST_START_TIME_OFFSET))
		goto METADATA_SECTION_3_NOT_ALIGNED;
	if (&md3->DST_end_time != (si8 *) (bytes + METADATA_DST_END_TIME_OFFSET))
		goto METADATA_SECTION_3_NOT_ALIGNED;
	if (&md3->GMT_offset != (si4 *) (bytes + METADATA_GMT_OFFSET_OFFSET))
		goto METADATA_SECTION_3_NOT_ALIGNED;
	if (md3->subject_name_1 != (si1 *) (bytes + METADATA_SUBJECT_NAME_1_OFFSET))
		goto METADATA_SECTION_3_NOT_ALIGNED;
	if (md3->subject_name_2 != (si1 *) (bytes + METADATA_SUBJECT_NAME_2_OFFSET))
		goto METADATA_SECTION_3_NOT_ALIGNED;
	if (md3->subject_ID != (si1 *) (bytes + METADATA_SUBJECT_ID_OFFSET))
		goto METADATA_SECTION_3_NOT_ALIGNED;
	if (md3->recording_location != (si1 *) (bytes + METADATA_RECORDING_LOCATION_OFFSET))
		goto METADATA_SECTION_3_NOT_ALIGNED;
	if (md3->protected_region != (ui1 *) (bytes + METADATA_SECTION_3_PROTECTED_REGION_OFFSET))
		goto METADATA_SECTION_3_NOT_ALIGNED;
	if (md3->discretionary_region != (ui1 *) (bytes + METADATA_SECTION_3_DISCRETIONARY_REGION_OFFSET))
		goto METADATA_SECTION_3_NOT_ALIGNED;
	
	// aligned
        MEF_globals->metadata_section_3_aligned = MEF_TRUE;
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): METADATA_SECTION_3 structure is aligned\n", __FUNCTION__);
	
        return(MEF_TRUE);
	
	// not aligned
METADATA_SECTION_3_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	(void) fprintf(stderr, "%c%s(): METADATA_SECTION_3 structure is not aligned\n", 7, __FUNCTION__);
        
        
	return(MEF_FALSE);
}


si4	check_password(si1 *password, const si1 *function, si4 line)
{
        si4	pw_len;
        
        
        // check pointer
        if (password == NULL) {
		if (MEF_globals->verbose == MEF_TRUE)
			printf("%s(): password field points to NULL [called from function \"%s\", line %d]\n", __FUNCTION__, function, line);
                return(1);
        }
        
        // check password length
        pw_len = UTF8_strlen(password);
        if (pw_len == 0) {
		if (MEF_globals->verbose == MEF_TRUE)
			fprintf(stderr, "%s(): password has no characters [called from function \"%s\", line %d]\n", __FUNCTION__, function, line);
                return(1);
        }
        if (pw_len > MAX_PASSWORD_CHARACTERS) {
		if (MEF_globals->verbose == MEF_TRUE)
			fprintf(stderr, "%s() Error: password too long [called from function \"%s\", line %d]\n", __FUNCTION__, function, line);
                return(1);
        }
        
	if (MEF_globals->verbose == MEF_TRUE)
		fprintf(stderr, "%s(): password is of valid length [called from function \"%s\", line %d]\n", __FUNCTION__, function, line);
	
        
        return(0);
}


si4	check_record_header_alignment(ui1 *bytes)
{
	RECORD_HEADER	*rh;
	si4		free_flag = MEF_FALSE;
	
	
	// see if already checked
	if (MEF_globals->record_header_aligned == MEF_UNKNOWN)
		MEF_globals->record_header_aligned = MEF_FALSE;
	else
		return(MEF_globals->record_header_aligned);
	
	// check overall size
	if (sizeof(RECORD_HEADER) != RECORD_HEADER_BYTES)
		goto RECORD_HEADER_NOT_ALIGNED;
	
	// check fields
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(RECORD_HEADER_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	rh = (RECORD_HEADER *) bytes;
	if (&rh->record_CRC != (ui4 *) (bytes + RECORD_HEADER_RECORD_CRC_OFFSET))
		goto RECORD_HEADER_NOT_ALIGNED;
	if (rh->type_string != (si1 *) (bytes + RECORD_HEADER_TYPE_OFFSET))
		goto RECORD_HEADER_NOT_ALIGNED;
	if (&rh->version_major != (ui1 *) (bytes + RECORD_HEADER_VERSION_MAJOR_OFFSET))
		goto RECORD_HEADER_NOT_ALIGNED;
	if (&rh->version_minor != (ui1 *) (bytes + RECORD_HEADER_VERSION_MINOR_OFFSET))
		goto RECORD_HEADER_NOT_ALIGNED;
	if (&rh->encryption != (si1 *) (bytes + RECORD_HEADER_ENCRYPTION_OFFSET))
		goto RECORD_HEADER_NOT_ALIGNED;
	if (&rh->bytes != (ui4 *) (bytes + RECORD_HEADER_BYTES_OFFSET))
		goto RECORD_HEADER_NOT_ALIGNED;
	if (&rh->time != (si8 *) (bytes + RECORD_HEADER_TIME_OFFSET))
		goto RECORD_HEADER_NOT_ALIGNED;
	
	// aligned
        MEF_globals->record_header_aligned = MEF_TRUE;
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): RECORD_HEADER structure is aligned\n", __FUNCTION__);
	
        return(MEF_TRUE);
	
	// not aligned
RECORD_HEADER_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	(void) fprintf(stderr, "%c%s(): RECORD_HEADER structure is not aligned\n", 7, __FUNCTION__);
        
        
        return(MEF_FALSE);
}


si4	check_record_indices_alignment(ui1 *bytes)
{
	RECORD_INDEX	*ri;
	si4		free_flag = MEF_FALSE;
	
	
	// see if already checked
	if (MEF_globals->record_indices_aligned == MEF_UNKNOWN)
		MEF_globals->record_indices_aligned = MEF_FALSE;
	else
		return(MEF_globals->record_indices_aligned);
	
	// check overall size
	if (sizeof(RECORD_INDEX) != RECORD_INDEX_BYTES)
		goto RECORD_INDICES_NOT_ALIGNED;
	
	// check fields
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(RECORD_INDEX_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	ri = (RECORD_INDEX *) bytes;
	if (ri->type_string != (si1 *) (bytes + RECORD_INDEX_TYPE_OFFSET))
		goto RECORD_INDICES_NOT_ALIGNED;
	if (&ri->version_major != (ui1 *) (bytes + RECORD_INDEX_VERSION_MAJOR_OFFSET))
		goto RECORD_INDICES_NOT_ALIGNED;
	if (&ri->version_minor != (ui1 *) (bytes + RECORD_INDEX_VERSION_MINOR_OFFSET))
		goto RECORD_INDICES_NOT_ALIGNED;
	if (&ri->encryption != (si1 *) (bytes + RECORD_INDEX_ENCRYPTION_OFFSET))
		goto RECORD_INDICES_NOT_ALIGNED;
	if (&ri->file_offset != (si8 *) (bytes + RECORD_INDEX_FILE_OFFSET_OFFSET))
		goto RECORD_INDICES_NOT_ALIGNED;
	if (&ri->time != (si8 *) (bytes + RECORD_INDEX_TIME_OFFSET))
		goto RECORD_INDICES_NOT_ALIGNED;
	
	// aligned
        MEF_globals->record_indices_aligned = MEF_TRUE;
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): RECORD_INDEX structure is aligned\n", __FUNCTION__);
        
        return(MEF_TRUE);
	
	// not aligned
RECORD_INDICES_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	(void) fprintf(stderr, "%c%s(): RECORD_INDEX structure is not aligned\n", 7, __FUNCTION__);
        
        
        return(MEF_FALSE);
}


si4	check_RED_block_header_alignment(ui1 *bytes)
{
	RED_BLOCK_HEADER	*rbh;
	si4			free_flag = MEF_FALSE;
	
	
	// see if already checked
	if (MEF_globals->RED_block_header_aligned == MEF_UNKNOWN)
		MEF_globals->RED_block_header_aligned = MEF_FALSE;
	else
		return(MEF_globals->RED_block_header_aligned);
	
        // check overall size
	if (sizeof(RED_BLOCK_HEADER) != RED_BLOCK_HEADER_BYTES)
		goto RED_BLOCK_HEADER_NOT_ALIGNED;
	
	// check fields
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(RED_BLOCK_HEADER_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	rbh = (RED_BLOCK_HEADER *) bytes;
	if (&rbh->block_CRC != ((ui4 *) bytes + RED_BLOCK_BLOCK_CRC_OFFSET))
		goto RED_BLOCK_HEADER_NOT_ALIGNED;
	if (&rbh->flags != (ui1 *) (bytes + RED_BLOCK_FLAGS_OFFSET))
		goto RED_BLOCK_HEADER_NOT_ALIGNED;
	if (rbh->protected_region != (ui1 *) (bytes + RED_BLOCK_PROTECTED_REGION_OFFSET))
		goto RED_BLOCK_HEADER_NOT_ALIGNED;
	if (rbh->discretionary_region != (ui1 *) (bytes + RED_BLOCK_DISCRETIONARY_REGION_OFFSET))
		goto RED_BLOCK_HEADER_NOT_ALIGNED;
	if (&rbh->detrend_slope != (sf4 *) (bytes + RED_BLOCK_DETREND_SLOPE_OFFSET))
		goto RED_BLOCK_HEADER_NOT_ALIGNED;
	if (&rbh->detrend_intercept != (sf4 *) (bytes + RED_BLOCK_DETREND_INTERCEPT_OFFSET))
		goto RED_BLOCK_HEADER_NOT_ALIGNED;
	if (&rbh->scale_factor != (sf4 *) (bytes + RED_BLOCK_SCALE_FACTOR_OFFSET))
		goto RED_BLOCK_HEADER_NOT_ALIGNED;
	if (&rbh->difference_bytes != (ui4 *) (bytes + RED_BLOCK_DIFFERENCE_BYTES_OFFSET))
		goto RED_BLOCK_HEADER_NOT_ALIGNED;
	if (&rbh->number_of_samples != (ui4 *) (bytes + RED_BLOCK_NUMBER_OF_SAMPLES_OFFSET))
		goto RED_BLOCK_HEADER_NOT_ALIGNED;
	if (&rbh->block_bytes != (ui4 *) (bytes + RED_BLOCK_BLOCK_BYTES_OFFSET))
		goto RED_BLOCK_HEADER_NOT_ALIGNED;
	if (&rbh->start_time != (si8 *) (bytes + RED_BLOCK_START_TIME_OFFSET))
		goto RED_BLOCK_HEADER_NOT_ALIGNED;
	if (rbh->statistics != (ui1 *) (bytes + RED_BLOCK_STATISTICS_OFFSET))
		goto RED_BLOCK_HEADER_NOT_ALIGNED;
	
	// aligned
        MEF_globals->RED_block_header_aligned = MEF_TRUE;
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): RED_BLOCK_HEADER structure is aligned\n", __FUNCTION__);
	
        return(MEF_TRUE);
	
	// not aligned
RED_BLOCK_HEADER_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
        
	(void) fprintf(stderr, "%c%s(): RED_BLOCK_HEADER structure is not aligned\n", 7, __FUNCTION__);
	
        
	return(MEF_FALSE);
}


si4	check_time_series_indices_alignment(ui1 *bytes)
{
	TIME_SERIES_INDEX	*tsi;
	si4			free_flag = MEF_FALSE;
	
	
	// see if already checked
	if (MEF_globals->time_series_indices_aligned == MEF_UNKNOWN)
		MEF_globals->time_series_indices_aligned = MEF_FALSE;
	else
		return(MEF_globals->time_series_indices_aligned);
	
	// check overall size
	if (sizeof(TIME_SERIES_INDEX) != TIME_SERIES_INDEX_BYTES)
		goto TIME_SERIES_INDICES_NOT_ALIGNED;
	
	// check fields
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(TIME_SERIES_INDEX_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	tsi = (TIME_SERIES_INDEX *) bytes;
	if (&tsi->file_offset != (si8 *) (bytes + TIME_SERIES_INDEX_FILE_OFFSET_OFFSET))
		goto TIME_SERIES_INDICES_NOT_ALIGNED;
	if (&tsi->start_time != (si8 *) (bytes + TIME_SERIES_INDEX_START_TIME_OFFSET))
		goto TIME_SERIES_INDICES_NOT_ALIGNED;
	if (&tsi->start_sample != (si8 *) (bytes + TIME_SERIES_INDEX_START_SAMPLE_OFFSET))
		goto TIME_SERIES_INDICES_NOT_ALIGNED;
	if (&tsi->number_of_samples != (ui4 *) (bytes + TIME_SERIES_INDEX_NUMBER_OF_SAMPLES_OFFSET))
		goto TIME_SERIES_INDICES_NOT_ALIGNED;
	if (&tsi->block_bytes != (ui4 *) (bytes + TIME_SERIES_INDEX_BLOCK_BYTES_OFFSET))
		goto TIME_SERIES_INDICES_NOT_ALIGNED;
	if (&tsi->maximum_sample_value != (si4 *) (bytes + TIME_SERIES_INDEX_MAXIMUM_SAMPLE_VALUE_OFFSET))
		goto TIME_SERIES_INDICES_NOT_ALIGNED;
	if (&tsi->minimum_sample_value != (si4 *) (bytes + TIME_SERIES_INDEX_MINIMUM_SAMPLE_VALUE_OFFSET))
		goto TIME_SERIES_INDICES_NOT_ALIGNED;
	if (tsi->protected_region != (ui1 *) (bytes + TIME_SERIES_INDEX_PROTECTED_REGION_OFFSET))
		goto TIME_SERIES_INDICES_NOT_ALIGNED;
	if (&tsi->RED_block_flags != (ui1 *) (bytes + TIME_SERIES_INDEX_RED_BLOCK_FLAGS_OFFSET))
		goto TIME_SERIES_INDICES_NOT_ALIGNED;
	if (tsi->RED_block_protected_region != (ui1 *) (bytes + TIME_SERIES_INDEX_RED_BLOCK_PROTECTED_REGION_OFFSET))
		goto TIME_SERIES_INDICES_NOT_ALIGNED;
	if (tsi->RED_block_discretionary_region != (ui1 *) (bytes + TIME_SERIES_INDEX_RED_BLOCK_DISCRETIONARY_REGION_OFFSET))
		goto TIME_SERIES_INDICES_NOT_ALIGNED;
	
	// aligned
	MEF_globals->time_series_indices_aligned = MEF_TRUE;
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): TIME_SERIES_INDEX structure is aligned\n", __FUNCTION__);
	
	return(MEF_TRUE);
	
	// not aligned
TIME_SERIES_INDICES_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	(void) fprintf(stderr, "%c%s(): TIME_SERIES_INDEX structure is not aligned\n", 7, __FUNCTION__);
	
	
	return(MEF_FALSE);
}


si4	check_time_series_metadata_section_2_alignment(ui1 *bytes)
{
	TIME_SERIES_METADATA_SECTION_2	*md2;
	si4				free_flag = MEF_FALSE;
	
	
	// see if already checked
	if (MEF_globals->time_series_metadata_section_2_aligned == MEF_UNKNOWN)
		MEF_globals->time_series_metadata_section_2_aligned = MEF_FALSE;
	else
		return(MEF_globals->time_series_metadata_section_2_aligned);
	
	// check overall size
	if (sizeof(TIME_SERIES_METADATA_SECTION_2) != METADATA_SECTION_2_BYTES)
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	
	// check fields
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(METADATA_FILE_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	md2 = (TIME_SERIES_METADATA_SECTION_2 *) (bytes + METADATA_SECTION_2_OFFSET);
        // type-independent fields
        if (md2->channel_description != (si1 *) (bytes + METADATA_CHANNEL_DESCRIPTION_OFFSET))
                goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
        if (md2->session_description != (si1 *) (bytes + METADATA_SESSION_DESCRIPTION_OFFSET))
                goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->recording_duration != (si8 *) (bytes + METADATA_RECORDING_DURATION_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
        // type-specific fields
	if (md2->reference_description != (si1 *) (bytes + TIME_SERIES_METADATA_REFERENCE_DESCRIPTION_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->acquisition_channel_number != (si8 *) (bytes + TIME_SERIES_METADATA_ACQUISITION_CHANNEL_NUMBER_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->sampling_frequency != (sf8 *) (bytes + TIME_SERIES_METADATA_SAMPLING_FREQUENCY_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->low_frequency_filter_setting != (sf8 *) (bytes + TIME_SERIES_METADATA_LOW_FREQUENCY_FILTER_SETTING_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->high_frequency_filter_setting != (sf8 *) (bytes + TIME_SERIES_METADATA_HIGH_FREQUENCY_FILTER_SETTING_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->notch_filter_frequency_setting != (sf8 *) (bytes + TIME_SERIES_METADATA_NOTCH_FILTER_FREQUENCY_SETTING_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->AC_line_frequency != (sf8 *) (bytes + TIME_SERIES_METADATA_AC_LINE_FREQUENCY_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->units_conversion_factor != (sf8 *) (bytes + TIME_SERIES_METADATA_UNITS_CONVERSION_FACTOR_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (md2->units_description != (si1 *) (bytes + TIME_SERIES_METADATA_UNITS_DESCRIPTION_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->maximum_native_sample_value != (sf8 *) (bytes + TIME_SERIES_METADATA_MAXIMUM_NATIVE_SAMPLE_VALUE_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->minimum_native_sample_value != (sf8 *) (bytes + TIME_SERIES_METADATA_MINIMUM_NATIVE_SAMPLE_VALUE_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
        if (&md2->start_sample != (si8 *) (bytes + TIME_SERIES_METADATA_START_SAMPLE_OFFSET))
                goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->number_of_samples != (si8 *) (bytes + TIME_SERIES_METADATA_NUMBER_OF_SAMPLES_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->number_of_blocks != (si8 *) (bytes + TIME_SERIES_METADATA_NUMBER_OF_BLOCKS_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->maximum_block_samples != (ui4 *) (bytes + TIME_SERIES_METADATA_MAXIMUM_BLOCK_SAMPLES_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->maximum_difference_bytes != (ui4 *) (bytes + TIME_SERIES_METADATA_MAXIMUM_DIFFERENCE_BYTES_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->block_interval != (si8 *) (bytes + TIME_SERIES_METADATA_BLOCK_INTERVAL_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->number_of_discontinuities != (si8 *) (bytes + TIME_SERIES_METADATA_NUMBER_OF_DISCONTINUITIES_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->maximum_contiguous_blocks != (si8 *) (bytes + TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCKS_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->maximum_contiguous_block_bytes != (si8 *) (bytes + TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCK_BYTES_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (&md2->maximum_contiguous_samples != (si8 *) (bytes + TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_SAMPLES_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (md2->protected_region != (ui1 *) (bytes + TIME_SERIES_METADATA_SECTION_2_PROTECTED_REGION_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	if (md2->discretionary_region != (ui1 *) (bytes + TIME_SERIES_METADATA_SECTION_2_DISCRETIONARY_REGION_OFFSET))
		goto TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED;
	
	// aligned
	MEF_globals->time_series_metadata_section_2_aligned = MEF_TRUE;
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): TIME_SERIES_METADATA_SECTION_2 structure is aligned\n", __FUNCTION__);
	
	return(MEF_TRUE);
	
	// not aligned
TIME_SERIES_METADATA_SECTION_2_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) fprintf(stderr, "%c%s(): TIME_SERIES_METADATA_SECTION_2 structure is not aligned\n", 7, __FUNCTION__);
	
	
	return(MEF_FALSE);
}


si4	check_universal_header_alignment(ui1 *bytes)
{
	UNIVERSAL_HEADER	*uh;
	si4			free_flag = MEF_FALSE;
	
	
	// see if already checked
	if (MEF_globals->universal_header_aligned == MEF_UNKNOWN)
		MEF_globals->universal_header_aligned = MEF_FALSE;
	else
		return(MEF_globals->universal_header_aligned);
	
	// check overall size
	if (sizeof(UNIVERSAL_HEADER) != UNIVERSAL_HEADER_BYTES)
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	
	// check fields
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(UNIVERSAL_HEADER_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	uh = (UNIVERSAL_HEADER *) bytes;
	if (&uh->header_CRC != (ui4 *) (bytes + UNIVERSAL_HEADER_HEADER_CRC_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (&uh->body_CRC != (ui4 *) (bytes + UNIVERSAL_HEADER_BODY_CRC_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (uh->file_type_string != (si1 *) (bytes + UNIVERSAL_HEADER_FILE_TYPE_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (&uh->mef_version_major != (ui1 *) (bytes + UNIVERSAL_HEADER_MEF_VERSION_MAJOR_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (&uh->mef_version_minor != (ui1 *) (bytes + UNIVERSAL_HEADER_MEF_VERSION_MINOR_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (&uh->byte_order_code != (ui1 *) (bytes + UNIVERSAL_HEADER_BYTE_ORDER_CODE_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (&uh->start_time != (si8 *) (bytes + UNIVERSAL_HEADER_START_TIME_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (&uh->end_time != (si8 *) (bytes + UNIVERSAL_HEADER_END_TIME_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (&uh->number_of_entries != (si8 *) (bytes + UNIVERSAL_HEADER_NUMBER_OF_ENTRIES_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (&uh->maximum_entry_size != (si8 *) (bytes + UNIVERSAL_HEADER_MAXIMUM_ENTRY_SIZE_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (&uh->segment_number != (si4 *) (bytes + UNIVERSAL_HEADER_SEGMENT_NUMBER_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (uh->channel_name != (si1 *) (bytes + UNIVERSAL_HEADER_CHANNEL_NAME_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (uh->session_name != (si1 *) (bytes + UNIVERSAL_HEADER_SESSION_NAME_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (uh->anonymized_name != (si1 *) (bytes + UNIVERSAL_HEADER_ANONYMIZED_NAME_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (uh->level_UUID != (ui1 *) (bytes + UNIVERSAL_HEADER_LEVEL_UUID_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (uh->file_UUID != (ui1 *) (bytes + UNIVERSAL_HEADER_FILE_UUID_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (uh->provenance_UUID != (ui1 *) (bytes + UNIVERSAL_HEADER_PROVENANCE_UUID_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (uh->level_1_password_validation_field != (ui1 *) (bytes + UNIVERSAL_HEADER_LEVEL_1_PASSWORD_VALIDATION_FIELD_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (uh->level_2_password_validation_field != (ui1 *) (bytes + UNIVERSAL_HEADER_LEVEL_2_PASSWORD_VALIDATION_FIELD_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (uh->protected_region != (ui1 *) (bytes + UNIVERSAL_HEADER_PROTECTED_REGION_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	if (uh->discretionary_region != (ui1 *) (bytes + UNIVERSAL_HEADER_DISCRETIONARY_REGION_OFFSET))
		goto UNIVERSAL_HEADER_NOT_ALIGNED;
	
	// aligned
        MEF_globals->universal_header_aligned = MEF_TRUE;
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): UNIVERSAL_HEADER structure is aligned\n", __FUNCTION__);
	
        return(MEF_TRUE);
	
	// not aligned
UNIVERSAL_HEADER_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) fprintf(stderr, "%c%s(): UNIVERSAL_HEADER structure is not aligned\n", 7, __FUNCTION__);
	
        
        return(MEF_FALSE);
}


si4	check_video_indices_alignment(ui1 *bytes)
{
	VIDEO_INDEX	*vi;
	si4		free_flag = MEF_FALSE;
	
	
	// see if already checked
	if (MEF_globals->video_indices_aligned == MEF_UNKNOWN)
		MEF_globals->video_indices_aligned = MEF_FALSE;
	else
		return(MEF_globals->video_indices_aligned);
	
	// check overall size
	if (sizeof(VIDEO_INDEX) != VIDEO_INDEX_BYTES)
		goto VIDEO_INDICES_NOT_ALIGNED;
	
	// check fields
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(VIDEO_INDEX_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	vi = (VIDEO_INDEX *) bytes;
	if (&vi->start_time != (si8 *) (bytes + VIDEO_INDEX_START_TIME_OFFSET))
		goto VIDEO_INDICES_NOT_ALIGNED;
	if (&vi->end_time != (si8 *) (bytes + VIDEO_INDEX_END_TIME_OFFSET))
		goto VIDEO_INDICES_NOT_ALIGNED;
	if (&vi->start_frame != (ui4 *) (bytes + VIDEO_INDEX_START_FRAME_OFFSET))
		goto VIDEO_INDICES_NOT_ALIGNED;
	if (&vi->end_frame != (ui4 *) (bytes + VIDEO_INDEX_END_FRAME_OFFSET))
		goto VIDEO_INDICES_NOT_ALIGNED;
	if (&vi->file_offset != (si8 *) (bytes + VIDEO_INDEX_FILE_OFFSET_OFFSET))
		goto VIDEO_INDICES_NOT_ALIGNED;
	if (&vi->clip_bytes != (si8 *) (bytes + VIDEO_INDEX_CLIP_BYTES_OFFSET))
		goto VIDEO_INDICES_NOT_ALIGNED;
	if (vi->protected_region != (ui1 *) (bytes + VIDEO_INDEX_PROTECTED_REGION_OFFSET))
		goto VIDEO_INDICES_NOT_ALIGNED;
	if (vi->discretionary_region != (ui1 *) (bytes + VIDEO_INDEX_DISCRETIONARY_REGION_OFFSET))
		goto VIDEO_INDICES_NOT_ALIGNED;

	// aligned
	MEF_globals->video_indices_aligned = MEF_TRUE;
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): VIDEO_INDEX structure is aligned\n", __FUNCTION__);
	
	return(MEF_TRUE);
	
	// not aligned
VIDEO_INDICES_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	(void) fprintf(stderr, "%c%s(): VIDEO_INDEX structure is not aligned\n", 7, __FUNCTION__);
	
	
	return(MEF_FALSE);
}


si4	check_video_metadata_section_2_alignment(ui1 *bytes)
{
	VIDEO_METADATA_SECTION_2	*vmd2;
	si4				free_flag = MEF_FALSE;
	
	
	// see if already checked
	if (MEF_globals->video_metadata_section_2_aligned == MEF_UNKNOWN)
		MEF_globals->video_metadata_section_2_aligned = MEF_FALSE;
	else
		return(MEF_globals->video_metadata_section_2_aligned);
	
	// check overall size
	if (sizeof(VIDEO_METADATA_SECTION_2) != METADATA_SECTION_2_BYTES)
		goto VIDEO_METADATA_SECTION_2_NOT_ALIGNED;
	
	// check fields
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(METADATA_FILE_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	vmd2 = (VIDEO_METADATA_SECTION_2 *) (bytes + METADATA_SECTION_2_OFFSET);
        // type-independent fields
        if (vmd2->channel_description != (si1 *) (bytes + METADATA_CHANNEL_DESCRIPTION_OFFSET))
                goto VIDEO_METADATA_SECTION_2_NOT_ALIGNED;
        if (vmd2->session_description != (si1 *) (bytes + METADATA_SESSION_DESCRIPTION_OFFSET))
                goto VIDEO_METADATA_SECTION_2_NOT_ALIGNED;
        if (&vmd2->recording_duration != (si8 *) (bytes + METADATA_RECORDING_DURATION_OFFSET))
                goto VIDEO_METADATA_SECTION_2_NOT_ALIGNED;
        // type-specific fields
	if (&vmd2->horizontal_resolution != (si8 *) (bytes + VIDEO_METADATA_HORIZONTAL_RESOLUTION_OFFSET))
		goto VIDEO_METADATA_SECTION_2_NOT_ALIGNED;
	if (&vmd2->vertical_resolution != (si8 *) (bytes + VIDEO_METADATA_VERTICAL_RESOLUTION_OFFSET))
		goto VIDEO_METADATA_SECTION_2_NOT_ALIGNED;
	if (&vmd2->frame_rate != (sf8 *) (bytes + VIDEO_METADATA_FRAME_RATE_OFFSET))
		goto VIDEO_METADATA_SECTION_2_NOT_ALIGNED;
	if (&vmd2->number_of_clips != (si8 *) (bytes + VIDEO_METADATA_NUMBER_OF_CLIPS_OFFSET))
		goto VIDEO_METADATA_SECTION_2_NOT_ALIGNED;
	if (&vmd2->maximum_clip_bytes != (si8 *) (bytes + VIDEO_METADATA_MAXIMUM_CLIP_BYTES_OFFSET))
		goto VIDEO_METADATA_SECTION_2_NOT_ALIGNED;
	if (vmd2->video_format != (si1 *) (bytes + VIDEO_METADATA_VIDEO_FORMAT_OFFSET))
		goto VIDEO_METADATA_SECTION_2_NOT_ALIGNED;
	if (&vmd2->video_file_CRC != (ui4 *) (bytes + VIDEO_METADATA_VIDEO_FILE_CRC_OFFSET))
		goto VIDEO_METADATA_SECTION_2_NOT_ALIGNED;
	if (vmd2->protected_region != (ui1 *) (bytes + VIDEO_METADATA_SECTION_2_PROTECTED_REGION_OFFSET))
		goto VIDEO_METADATA_SECTION_2_NOT_ALIGNED;
	if (vmd2->discretionary_region != (ui1 *) (bytes + VIDEO_METADATA_SECTION_2_DISCRETIONARY_REGION_OFFSET))
		goto VIDEO_METADATA_SECTION_2_NOT_ALIGNED;
	
	// aligned
	MEF_globals->video_metadata_section_2_aligned = MEF_TRUE;
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): VIDEO_METADATA_SECTION_2 structure is aligned\n", __FUNCTION__);
	
	return(MEF_TRUE);
	
	// not aligned
VIDEO_METADATA_SECTION_2_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) fprintf(stderr, "%c%s(): VIDEO_METADATA_SECTION_2 structure is not aligned\n", 7, __FUNCTION__);
	
	
	return(MEF_FALSE);
}


/*************************************************************************/
/********************  END ALIGNMENT CHECKING FUNCTIONS  *****************/
/*************************************************************************/


inline si4     compare_sf8(const void *a, const void * b)
{
        if (*((sf8 *) a) > *((sf8 *) b))
                return(1);
        else if (*((sf8 *) a) < *((sf8 *) b))
                return(-1);
        return(0);
}


ui1	cpu_endianness()
{
	ui2	x = 1;
	
        
	return(*((ui1 *) &x));
}


/*************************************************************************/
/*****************************  CRC FUNCTIONS  ***************************/
/*************************************************************************/


inline ui4	CRC_calculate(ui1 *block_ptr, si8 block_bytes)
{
	ui4	crc;
	
	
	crc = CRC_update(block_ptr, block_bytes, CRC_START_VALUE);
	
        
	return(crc);
}


ui4	*CRC_initialize_table(si4 global_flag)
{
	ui4	*crc_table;
	
	
	crc_table = (ui4 *) e_calloc((size_t) CRC_TABLE_ENTRIES, sizeof(ui4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	{
		ui4 temp[CRC_TABLE_ENTRIES] = CRC_KOOPMAN32_KEY;
		memcpy(crc_table, temp, CRC_TABLE_ENTRIES * sizeof(ui4));
	}
	
	if (global_flag == MEF_TRUE) {
		MEF_globals->CRC_table = crc_table;
		return(NULL);
	}
	
        
	return(crc_table);
}


inline ui4	CRC_update(ui1 *block_ptr, si8 block_bytes, ui4 current_crc)
{
	si8	i;
	ui4	tmp;
	
	
	if (MEF_globals->CRC_table == NULL)
		(void) CRC_initialize_table(MEF_TRUE);
	
	for (i = block_bytes; i--;) {
		tmp = current_crc ^ (ui4) *block_ptr++;
		current_crc = (current_crc >> 8) ^ MEF_globals->CRC_table[tmp & 0xff];
	}
	
        
	return(current_crc);
}


inline si4	CRC_validate(ui1 *block_ptr, si8 block_bytes, ui4 crc_to_validate)
{
	ui4	crc;
	
	
	crc = CRC_calculate(block_ptr, block_bytes);
	
	if (crc == crc_to_validate)
		return(MEF_TRUE);
	
        
	return(MEF_FALSE);
}


/*************************************************************************/
/***************************  END CRC FUNCTIONS  *************************/
/*************************************************************************/


si4	decrypt_metadata(FILE_PROCESSING_STRUCT *fps)
{
	ui1		*ui1_p, *decryption_key;
	si4		i, decryption_blocks;
        PASSWORD_DATA	*pwd;
	
	
        pwd = fps->password_data;
        
        // section 2 decryption
        if (fps->metadata.section_1->section_2_encryption > NO_ENCRYPTION) {  // natively encrypted and currently encrypted
                if (pwd->access_level >= fps->metadata.section_1->section_2_encryption) {
                        if (fps->metadata.section_1->section_2_encryption == LEVEL_1_ENCRYPTION)
                                decryption_key = pwd->level_1_encryption_key;
                        else
                                decryption_key = pwd->level_2_encryption_key;
                        decryption_blocks = METADATA_SECTION_2_BYTES / ENCRYPTION_BLOCK_BYTES;
                        ui1_p = fps->raw_data + METADATA_SECTION_2_OFFSET;
                        for (i = 0; i < decryption_blocks; ++i) {
                                AES_decrypt(ui1_p, ui1_p, NULL, decryption_key);
                                ui1_p += ENCRYPTION_BLOCK_BYTES;
                        }
                        fps->metadata.section_1->section_2_encryption = -fps->metadata.section_1->section_2_encryption;  // mark as currently decrypted
                }
        }

        // section 3 decryption
        if (fps->metadata.section_1->section_3_encryption > NO_ENCRYPTION) {  // natively encrypted and currently encrypted
                if (pwd->access_level >= fps->metadata.section_1->section_3_encryption) {
                        if (fps->metadata.section_1->section_3_encryption == LEVEL_1_ENCRYPTION)
                                decryption_key = pwd->level_1_encryption_key;
                        else
                                decryption_key = pwd->level_2_encryption_key;
                        decryption_blocks = METADATA_SECTION_3_BYTES / ENCRYPTION_BLOCK_BYTES;
                        ui1_p = fps->raw_data + METADATA_SECTION_3_OFFSET;
                        for (i = 0; i < decryption_blocks; ++i) {
                                AES_decrypt(ui1_p, ui1_p, NULL, decryption_key);
                                ui1_p += ENCRYPTION_BLOCK_BYTES;
                        }
                        fps->metadata.section_1->section_3_encryption = -fps->metadata.section_1->section_3_encryption;  // mark as currently decrypted
		}
        }
	
	// set global RTOs
        if (fps->metadata.section_1->section_3_encryption <= NO_ENCRYPTION) {
		MEF_globals->recording_time_offset = fps->metadata.section_3->recording_time_offset;
		MEF_globals->DST_start_time = fps->metadata.section_3->DST_start_time;
		MEF_globals->DST_end_time = fps->metadata.section_3->DST_end_time;
		MEF_globals->GMT_offset = fps->metadata.section_3->GMT_offset;
	}
	
	
	return(0);
}


si4	decrypt_records(FILE_PROCESSING_STRUCT *fps)
{
        si1		CRC_validity;
        ui4		i, j, decryption_blocks, r_cnt;
        ui1		*ui1_p, *end_p, *decryption_key;
	si8		number_of_records;
        RECORD_HEADER	*record_header;
        PASSWORD_DATA	*pwd;
        
        
        ui1_p = fps->records;
        pwd = fps->password_data;
	number_of_records = fps->universal_header->number_of_entries;
	
	if (number_of_records == UNKNOWN_NUMBER_OF_ENTRIES) {  // can still process if not passed, but will fail on incomplete final record
		end_p = fps->raw_data + fps->raw_data_bytes;
                r_cnt = 0;
		while (ui1_p < end_p) {
			record_header = (RECORD_HEADER *) ui1_p;
			// Validate record CRC
                        if (MEF_globals->CRC_mode & (CRC_VALIDATE | CRC_VALIDATE_ON_INPUT)) {
                                if (record_header->encryption >= NO_ENCRYPTION) { // CRCs calculated on encrypted records if encryption is specified
                                	CRC_validity = CRC_validate(ui1_p + CRC_BYTES, RECORD_HEADER_BYTES + record_header->bytes - CRC_BYTES, record_header->record_CRC);
                                        if (CRC_validity == MEF_FALSE)
                                                fprintf(stderr, "Invalid record CRC detected in record %d\n", r_cnt);
                                } else
                                        fprintf(stderr, "Can't validate CRC on decrypted record %d\n", r_cnt);
                        }
			// apply or remove recording time offset if requested
                        if (MEF_globals->recording_time_offset_mode & (RTO_APPLY | RTO_APPLY_ON_INPUT))
                                apply_recording_time_offset(&record_header->time);
                        else if (MEF_globals->recording_time_offset_mode & (RTO_REMOVE | RTO_REMOVE_ON_INPUT))
                                remove_recording_time_offset(&record_header->time);
			// decrypt
			if ((record_header->encryption > NO_ENCRYPTION) && (pwd->access_level >= record_header->encryption)) {
				if (record_header->encryption == LEVEL_1_ENCRYPTION)
					decryption_key = pwd->level_1_encryption_key;
				else
					decryption_key = pwd->level_2_encryption_key;
				decryption_blocks = record_header->bytes / ENCRYPTION_BLOCK_BYTES;
				ui1_p += RECORD_HEADER_BYTES;
				for (i = 0; i < decryption_blocks; ++i) {
					AES_decrypt(ui1_p, ui1_p, NULL, decryption_key);
					ui1_p += ENCRYPTION_BLOCK_BYTES;
				}
				record_header->encryption = -record_header->encryption;  // mark as currently decrypted
			} else {
				ui1_p += (RECORD_HEADER_BYTES + record_header->bytes);
			}
                        ++r_cnt;
		}
	} else {   // use number_of_records if known
                for (i = 0; i < number_of_records; ++i) {
			record_header = (RECORD_HEADER *) ui1_p;
                        if (MEF_globals->CRC_mode & (CRC_VALIDATE | CRC_VALIDATE_ON_INPUT)) {
                                if (record_header->encryption >= NO_ENCRYPTION) { // CRCs calculated on encrypted records if encryption is specified
                                	CRC_validity = CRC_validate(ui1_p + CRC_BYTES, RECORD_HEADER_BYTES + record_header->bytes - CRC_BYTES, record_header->record_CRC);
                                        if (CRC_validity == MEF_FALSE)
                                                fprintf(stderr, "Invalid record CRC detected in record %d\n", i);
                                } else
                                        fprintf(stderr, "Can't validate CRC on decrypted record %d\n", i);
                        }
			// apply or remove recording time offset if requested
                        if (MEF_globals->recording_time_offset_mode & (RTO_APPLY | RTO_APPLY_ON_INPUT))
                                apply_recording_time_offset(&record_header->time);
                        else if (MEF_globals->recording_time_offset_mode & (RTO_REMOVE | RTO_REMOVE_ON_INPUT))
                                remove_recording_time_offset(&record_header->time);
			// decrypt
			if ((record_header->encryption > NO_ENCRYPTION) && (pwd->access_level >= record_header->encryption)) {
				if (record_header->encryption == LEVEL_1_ENCRYPTION)
					decryption_key = pwd->level_1_encryption_key;
				else
					decryption_key = pwd->level_2_encryption_key;
				decryption_blocks = record_header->bytes / ENCRYPTION_BLOCK_BYTES;
				ui1_p += RECORD_HEADER_BYTES;
				for (j = 0; j < decryption_blocks; ++j) {
					AES_decrypt(ui1_p, ui1_p, NULL, decryption_key);
					ui1_p += ENCRYPTION_BLOCK_BYTES;
				}
				record_header->encryption = -record_header->encryption;  // mark as currently decrypted
			} else {
				ui1_p += (RECORD_HEADER_BYTES + record_header->bytes);
			}
		}
	}
	
        
        return(0);
}


/*************************************************************************/
/******************  ERROR CHECKING STANDARD FUNCTIONS  ******************/
/*************************************************************************/

si4 	e_system(const char *command, const si1 *function, si4 line, ui4 behavior_on_fail)
{
	si4 sys_res;

	if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
		behavior_on_fail = MEF_globals->behavior_on_fail;

	if ((sys_res = system(command)) != 0) {
		if (!(behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
			(void) fprintf(stderr, "\tsystem error number %d (%s)\n", errno, strerror(errno));
			if (function != NULL)
				(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
			if (behavior_on_fail & RETURN_ON_FAIL)
				(void) fprintf(stderr, "\t=> returning -1\n\n");
			else if (behavior_on_fail & EXIT_ON_FAIL)
				(void) fprintf(stderr, "\t=> exiting program\n\n");
		}
		if (behavior_on_fail & RETURN_ON_FAIL)
			return(-1);
		else if (behavior_on_fail & EXIT_ON_FAIL)
			exit(1);
	}

	return(sys_res);
}

void	*e_calloc(size_t n_members, size_t size, const si1 *function, si4 line, ui4 behavior_on_fail)
{
	void	*ptr;
	
	
	if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
		behavior_on_fail = MEF_globals->behavior_on_fail;
	
	if ((ptr = calloc(n_members, size)) == NULL) {
		if (!(behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
			#ifdef _WIN32
				(void) fprintf(stderr, "%c\n\t%s() failed to allocate the requested array (%lld members of size %lld)\n", 7, __FUNCTION__, n_members, size);
			#else
				(void) fprintf(stderr, "%c\n\t%s() failed to allocate the requested array (%ld members of size %ld)\n", 7, __FUNCTION__, n_members, size);
			#endif
			(void) fprintf(stderr, "\tsystem error number %d (%s)\n", errno, strerror(errno));
			if (function != NULL)
				(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
			if (behavior_on_fail & RETURN_ON_FAIL)
				(void) fprintf(stderr, "\t=> returning NULL\n\n");
			else if (behavior_on_fail & EXIT_ON_FAIL)
				(void) fprintf(stderr, "\t=> exiting program\n\n");
		}
		if (behavior_on_fail & RETURN_ON_FAIL)
			return(NULL);
		else if (behavior_on_fail & EXIT_ON_FAIL)
			exit(1);
	}
	
        
	return(ptr);
}


FILE	*e_fopen(si1 *path, si1 *mode, const si1 *function, si4 line, ui4 behavior_on_fail)
{
	FILE	 *fp;
	
	
	if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
		behavior_on_fail = MEF_globals->behavior_on_fail;
	
	if ((fp = fopen(path, mode)) == NULL) {
		if (!(behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
			(void) UTF8_fprintf(stderr, "%c\n\t%s() failed to open file \"%s\"\n", 7, __FUNCTION__, path);
			(void) fprintf(stderr, "\tsystem error number %d (%s)\n", errno, strerror(errno));
			if (function != NULL)
				(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
			if (behavior_on_fail & RETURN_ON_FAIL)
				(void) fprintf(stderr, "\t=> returning NULL\n\n");
			else if (behavior_on_fail & EXIT_ON_FAIL)
				(void) fprintf(stderr, "\t=> exiting program\n\n");
		}
		if (behavior_on_fail & RETURN_ON_FAIL)
			return(NULL);
		else if (behavior_on_fail & EXIT_ON_FAIL)
			exit(1);
	}
	
        
	return(fp);
}


size_t	e_fread(void *ptr, size_t size, size_t n_members, FILE *stream, si1 *path, const si1 *function, si4 line, ui4 behavior_on_fail)
{
	size_t	nr;
	
	
	if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
		behavior_on_fail = MEF_globals->behavior_on_fail;
	
	if ((nr = fread(ptr, size, n_members, stream)) != n_members) {
		if (!(behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
			(void) UTF8_fprintf(stderr, "%c\n\t%s() failed to read file \"%s\"\n", 7, __FUNCTION__, path);
			(void) fprintf(stderr, "\tsystem error number %d (%s)\n", errno, strerror(errno));
			if (function != NULL)
				(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
			if (behavior_on_fail & RETURN_ON_FAIL)
				(void) fprintf(stderr, "\t=> returning number of items read\n\n");
			else if (behavior_on_fail & EXIT_ON_FAIL)
				(void) fprintf(stderr, "\t=> exiting program\n\n");
		}
		if (behavior_on_fail & RETURN_ON_FAIL)
			return(nr);
		else if (behavior_on_fail & EXIT_ON_FAIL)
			exit(1);
	}
	
        
	return(nr);
}


si4	e_fseek(FILE *stream, size_t offset, si4 whence, si1 *path, const si1 *function, si4 line, ui4 behavior_on_fail)
{
	if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
		behavior_on_fail = MEF_globals->behavior_on_fail;
	
	if ((fseek(stream, offset, whence)) == -1) {
		if (!(behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
			#ifdef _WIN32
				(void) fprintf(stderr, "%c\n\t%s() failed to move the file pointer to requested location (offset %lld, whence %d)\n", 7, __FUNCTION__, offset, whence);
			#else
				(void) fprintf(stderr, "%c\n\t%s() failed to move the file pointer to requested location (offset %ld, whence %d)\n", 7, __FUNCTION__, offset, whence);
			#endif
			(void) UTF8_fprintf(stderr, "%\tin file \"%s\"\n", path);
			(void) fprintf(stderr, "\tsystem error number %d (%s)\n", errno, strerror(errno));
			if (function != NULL)
				(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
			if (behavior_on_fail & RETURN_ON_FAIL)
				(void) fprintf(stderr, "\t=> returning -1\n\n");
			else if (behavior_on_fail & EXIT_ON_FAIL)
				(void) fprintf(stderr, "\t=> exiting program\n\n");
		}
		if (behavior_on_fail & RETURN_ON_FAIL)
			return(-1);
		else if (behavior_on_fail & EXIT_ON_FAIL)
			exit(1);
	}
	
        
	return(0);
}


long	e_ftell(FILE *stream, const si1 *function, si4 line, ui4 behavior_on_fail)
{
	long	pos;
	
	
	if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
		behavior_on_fail = MEF_globals->behavior_on_fail;
	
	if ((pos = ftell(stream)) == -1) {
		if (!(behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
			(void) fprintf(stderr, "%c\n\t%s() failed obtain the current location\n", 7, __FUNCTION__);
			(void) fprintf(stderr, "\tsystem error number %d (%s)\n", errno, strerror(errno));
			if (function != NULL)
				(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
			if (behavior_on_fail & RETURN_ON_FAIL)
				(void) fprintf(stderr, "\t=> returning -1\n\n");
			else if (behavior_on_fail & EXIT_ON_FAIL)
				(void) fprintf(stderr, "\t=> exiting program\n\n");
		}
		if (behavior_on_fail & RETURN_ON_FAIL)
			return(-1);
		else if (behavior_on_fail & EXIT_ON_FAIL)
			exit(1);
	}
	
        
	return(pos);
}


size_t	e_fwrite(void *ptr, size_t size, size_t n_members, FILE *stream, si1 *path, const si1 *function, si4 line, ui4 behavior_on_fail)
{
	size_t	nw;
	
	
	if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
		behavior_on_fail = MEF_globals->behavior_on_fail;
	
	if ((nw = fwrite(ptr, size, n_members, stream)) != n_members) {
		if (!(behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
			(void) UTF8_fprintf(stderr, "%c\n\t%s() failed to write file \"%s\"\n", 7, __FUNCTION__, path);
			(void) fprintf(stderr, "\tsystem error number %d (%s)\n", errno, strerror(errno));
			if (function != NULL)
				(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
			if (behavior_on_fail & RETURN_ON_FAIL)
				(void) fprintf(stderr, "\t=> returning number of items written\n\n");
			else if (behavior_on_fail & EXIT_ON_FAIL)
				(void) fprintf(stderr, "\t=> exiting program\n\n");
		}
		if (behavior_on_fail & RETURN_ON_FAIL)
			return(nw);
		else if (behavior_on_fail & EXIT_ON_FAIL)
			exit(1);
	}
	
        
	return(nw);
}


void	*e_malloc(size_t n_bytes, const si1 *function, si4 line, ui4 behavior_on_fail)
{
	void	 *ptr;
	
	
	if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
		behavior_on_fail = MEF_globals->behavior_on_fail;
	
	if ((ptr = malloc(n_bytes)) == NULL) {
		if (!(behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
			#ifdef _WIN32
				(void) fprintf(stderr, "%c\n\t%s() failed to allocate the requested array (%lld bytes)\n", 7, __FUNCTION__, n_bytes);
			#else
				(void) fprintf(stderr, "%c\n\t%s() failed to allocate the requested array (%ld bytes)\n", 7, __FUNCTION__, n_bytes);
			#endif
			(void) fprintf(stderr, "\tsystem error number %d (%s)\n", errno, strerror(errno));
			if (function != NULL)
				(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
			if (behavior_on_fail & RETURN_ON_FAIL)
				(void) fprintf(stderr, "\t=> returning NULL\n\n");
			else if (behavior_on_fail & EXIT_ON_FAIL)
				(void) fprintf(stderr, "\t=> exiting program\n\n");
		}
		if (behavior_on_fail & RETURN_ON_FAIL)
			return(NULL);
		else if (behavior_on_fail & EXIT_ON_FAIL)
			exit(1);
	}
	
	
	return(ptr);
}


void	*e_realloc(void *ptr, size_t n_bytes, const si1 *function, si4 line, ui4 behavior_on_fail)
{
	if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
		behavior_on_fail = MEF_globals->behavior_on_fail;
	
	if ((ptr = realloc(ptr, n_bytes)) == NULL) {
		if (!(behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
			#ifdef _WIN32
				(void) fprintf(stderr, "%c\n\t%s() failed to reallocate the requested array (%lld bytes)\n", 7, __FUNCTION__, n_bytes);
			#else
				(void) fprintf(stderr, "%c\n\t%s() failed to reallocate the requested array (%ld bytes)\n", 7, __FUNCTION__, n_bytes);
			#endif
			(void) fprintf(stderr, "\tsystem error number %d (%s)\n", errno, strerror(errno));
			if (function != NULL)
				(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
			if (behavior_on_fail & RETURN_ON_FAIL)
				(void) fprintf(stderr, "\t=> returning NULL\n\n");
			else if (behavior_on_fail & EXIT_ON_FAIL)
				(void) fprintf(stderr, "\t=> exiting program\n\n");
		}
		if (behavior_on_fail & RETURN_ON_FAIL)
			return(NULL);
		else if (behavior_on_fail & EXIT_ON_FAIL)
			exit(1);
	}
	
	
	return(ptr);
}


/*************************************************************************/
/****************  END ERROR CHECKING STANDARD FUNCTIONS  ****************/
/*************************************************************************/


si4	encrypt_metadata(FILE_PROCESSING_STRUCT *fps)
{
	ui1		*ui1_p, *encryption_key;
	si4		i, encryption_blocks;
        PASSWORD_DATA	*pwd;
        
        
	pwd = fps->password_data;
	
	// section 2 encrypt
	if (fps->metadata.section_1->section_2_encryption < NO_ENCRYPTION) {  // natively encrypted and currently decrypted
		if (pwd->access_level >= -fps->metadata.section_1->section_2_encryption) {
			fps->metadata.section_1->section_2_encryption = -fps->metadata.section_1->section_2_encryption;  // mark as currently encrypted
			if (fps->metadata.section_1->section_2_encryption == LEVEL_1_ENCRYPTION)
				encryption_key = pwd->level_1_encryption_key;
			else
				encryption_key = pwd->level_2_encryption_key;
			encryption_blocks = METADATA_SECTION_2_BYTES / ENCRYPTION_BLOCK_BYTES;
			ui1_p = fps->raw_data + METADATA_SECTION_2_OFFSET;
			for (i = 0; i < encryption_blocks; ++i) {
				AES_encrypt(ui1_p, ui1_p, NULL, encryption_key);
				ui1_p += ENCRYPTION_BLOCK_BYTES;
			}
		}
	}
        
	// section 3 encrypt
	if (fps->metadata.section_1->section_3_encryption < NO_ENCRYPTION) {  // natively encrypted and currently decrypted
		if (pwd->access_level >= -fps->metadata.section_1->section_3_encryption) {
			fps->metadata.section_1->section_3_encryption = -fps->metadata.section_1->section_3_encryption;  // mark as currently encrypted
			if (fps->metadata.section_1->section_3_encryption == LEVEL_1_ENCRYPTION)
				encryption_key = pwd->level_1_encryption_key;
			else
				encryption_key = pwd->level_2_encryption_key;
			encryption_blocks = METADATA_SECTION_3_BYTES / ENCRYPTION_BLOCK_BYTES;
			ui1_p = fps->raw_data + METADATA_SECTION_3_OFFSET;
			for (i = 0; i < encryption_blocks; ++i) {
				AES_encrypt(ui1_p, ui1_p, NULL, encryption_key);
				ui1_p += ENCRYPTION_BLOCK_BYTES;
			}
		}
	}
	
        
	return(0);
}


si4	encrypt_records(FILE_PROCESSING_STRUCT *fps)
{
        ui4		i, j, encryption_blocks;
        ui1		*ui1_p, *end_p, *encryption_key;
	si8		number_of_records;
        RECORD_HEADER	*record_header;
        PASSWORD_DATA	*pwd;
        
        
	ui1_p = fps->records;
        pwd = fps->password_data;
	number_of_records = fps->universal_header->number_of_entries;
        
	if (number_of_records == UNKNOWN_NUMBER_OF_ENTRIES) {  // can still process if not passed, but will fail on incomplete final record
		end_p = fps->raw_data + fps->raw_data_bytes;
		while (ui1_p < end_p) {
			record_header = (RECORD_HEADER *) ui1_p;
			if (record_header->encryption < NO_ENCRYPTION) {
				record_header->encryption = -record_header->encryption;  // mark as currently encrypted
				if (record_header->encryption == LEVEL_1_ENCRYPTION)
					encryption_key = pwd->level_1_encryption_key;
				else
					encryption_key = pwd->level_2_encryption_key;
				encryption_blocks = record_header->bytes / ENCRYPTION_BLOCK_BYTES;
				ui1_p += RECORD_HEADER_BYTES;
				for (i = 0; i < encryption_blocks; ++i) {
					AES_encrypt(ui1_p, ui1_p, NULL, encryption_key);
					ui1_p += ENCRYPTION_BLOCK_BYTES;
				}
			} else {
				ui1_p += (RECORD_HEADER_BYTES + record_header->bytes);
			}
			// apply or remove recording time offset if requested
                        if (record_header->time != UUTC_NO_ENTRY) {
                        	if (MEF_globals->recording_time_offset_mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT))
					apply_recording_time_offset(&record_header->time);
				else if (MEF_globals->recording_time_offset_mode & (RTO_REMOVE | RTO_REMOVE_ON_OUTPUT))
					remove_recording_time_offset(&record_header->time);
                        }
			// Calculate record CRC
                        if (MEF_globals->CRC_mode & (CRC_CALCULATE | CRC_CALCULATE_ON_OUTPUT))
                        	record_header->record_CRC = CRC_calculate((ui1 *) record_header + CRC_BYTES, RECORD_HEADER_BYTES + record_header->bytes - CRC_BYTES);
		}
	} else {   // use number_of_records if known
                for (i = 0; i < number_of_records; ++i) {
			record_header = (RECORD_HEADER *) ui1_p;
			if (record_header->encryption < NO_ENCRYPTION) {
				record_header->encryption = -record_header->encryption;  // mark as currently encrypted
				if (record_header->encryption == LEVEL_1_ENCRYPTION)
					encryption_key = pwd->level_1_encryption_key;
				else
					encryption_key = pwd->level_2_encryption_key;
				encryption_blocks = record_header->bytes / ENCRYPTION_BLOCK_BYTES;
				ui1_p += RECORD_HEADER_BYTES;
				for (j = 0; j < encryption_blocks; ++j) {
					AES_encrypt(ui1_p, ui1_p, NULL, encryption_key);
					ui1_p += ENCRYPTION_BLOCK_BYTES;
				}
			} else {
				ui1_p += (RECORD_HEADER_BYTES + record_header->bytes);
			}
			// apply or remove recording time offset if requested
                        if (record_header->time != UUTC_NO_ENTRY) {
                        	if (MEF_globals->recording_time_offset_mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT))
					apply_recording_time_offset(&record_header->time);
				else if (MEF_globals->recording_time_offset_mode & (RTO_REMOVE | RTO_REMOVE_ON_OUTPUT))
					remove_recording_time_offset(&record_header->time);
                        }
			// Calculate record CRC
                        if (MEF_globals->CRC_mode & (CRC_CALCULATE | CRC_CALCULATE_ON_OUTPUT))
                        	record_header->record_CRC = CRC_calculate((ui1 *) record_header + CRC_BYTES, RECORD_HEADER_BYTES + record_header->bytes - CRC_BYTES);
		}
        }
        
        
        return(0);
}

#ifdef _WIN32
	si4	extract_path_parts(si1 *full_file_name, si1 *path, si1 *name, si1 *extension)
	{
		si1	*c, *cc, *cwd, temp_full_file_name[MEF_FULL_FILE_NAME_BYTES];
		
	    slash_to_backslash(full_file_name);
		MEF_strncpy(temp_full_file_name, full_file_name, MEF_FULL_FILE_NAME_BYTES);  // do non-destructively
		
		// move pointer to end of string
		c = temp_full_file_name + strlen(temp_full_file_name) - 1;
	        
		// remove terminal "\\" if present
		if (*c == '\\')
			*c-- = 0;
		
		// step back to first extension
		cc = c;
		while (*--c != '.') {
			if (*c == '\\') {
				c = cc;
				break;
			}
		}
		
		// copy extension if allocated
		if (extension != NULL) {
			if (*c == '.')
				MEF_strncpy(extension, c + 1, TYPE_BYTES);
			else
				bzero(extension, TYPE_BYTES);
		}
		if (*c == '.')
			*c-- = 0;
			
		// step back to next directory break
		while (*--c != '\\');
		
		// copy name if allocated
		if (name != NULL)
			MEF_strncpy(name, c + 1, MEF_BASE_FILE_NAME_BYTES);
		*c = 0;
		
		// copy path if allocated
		if (path != NULL)
			MEF_strncpy(path, temp_full_file_name, MEF_FULL_FILE_NAME_BYTES);
		
		return(0);
	}
#else
	si4	extract_path_parts(si1 *full_file_name, si1 *path, si1 *name, si1 *extension)
	{
		si1	*c, *cc, *cwd, temp_full_file_name[MEF_FULL_FILE_NAME_BYTES];
		
		// check that path starts from root
	    if (*full_file_name == '/') {
			MEF_strncpy(temp_full_file_name, full_file_name, MEF_FULL_FILE_NAME_BYTES);  // do non-destructively
		} else {
			if (!(MEF_globals->behavior_on_fail & SUPPRESS_ERROR_OUTPUT))
				(void) fprintf(stderr, "%s() Warning: path \"%s\" does not start from root => prepending current working directory\n", __FUNCTION__, full_file_name);
			cwd = getenv("PWD");
			c = full_file_name;
			if (*c == '.' && *(c + 1) != '.')
				++c;
			if (*c == '/')
				++c;
			MEF_snprintf(temp_full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s", cwd, c);
		}
		
		// move pointer to end of string
		c = temp_full_file_name + strlen(temp_full_file_name) - 1;
	        
		// remove terminal "/" if present
		if (*c == '/')
			*c-- = 0;
		
		// step back to first extension
		cc = c;
		while (*--c != '.') {
			if (*c == '/') {
				c = cc;
				break;
			}
		}
		
		// copy extension if allocated
		if (extension != NULL) {
			if (*c == '.')
				MEF_strncpy(extension, c + 1, TYPE_BYTES);
			else
				bzero(extension, TYPE_BYTES);
		}
		if (*c == '.')
			*c-- = 0;
			
		// step back to next directory break
		while (*--c != '/');
		
		// copy name if allocated
		if (name != NULL)
			MEF_strncpy(name, c + 1, MEF_BASE_FILE_NAME_BYTES);
		*c = 0;
		
		// copy path if allocated
		if (path != NULL)
			MEF_strncpy(path, temp_full_file_name, MEF_FULL_FILE_NAME_BYTES);

		
		return(0);
	}
#endif

void	extract_terminal_password_bytes(si1 *password, si1 *password_bytes)
{
        si1			*s;		// terminal (most unique) bytes of UTF-8 password
	si4			i, j;
	ui4			ch;
	
        
        s = password;
        i = j = 0;
        do {
                ch = UTF8_nextchar(s, &i);
                password_bytes[j++] = (ui1) (ch & 0xFF);
        } while (ch);
        for (; j < PASSWORD_BYTES; ++j)
                password_bytes[j] = 0;
        
	
        return;
}


void	fill_empty_password_bytes(si1 *password_bytes)
{
        ui4	m_w, m_z;
        si4	i;
        
        
        // initialize random number generator
        m_w = m_z = 0;
        for (i = 0; i < PASSWORD_BYTES; ++i) {
                if (password_bytes[i] == 0)
                        break;
                m_w += password_bytes[i];
                m_z -= password_bytes[i];
        }
        if (m_w == 0 || m_w == 0x464FFFFF)
                m_w += 1;
        if (m_z == 0 || m_z == 0x9068FFFF)
                m_z -= 1;
        
        // fill in bytes
        for (; i < PASSWORD_BYTES; ++i)
                password_bytes[i] = random_byte(&m_w, &m_z);
  
        
        return;
}


/*************************************************************************/
/****************************  FILTER FUNCTIONS  *************************/
/*************************************************************************/


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


void	FILT_balance(sf16 **a, si4 poles)
{
        sf16    radix, sqrdx, c, r, g, f, s;
        si4     i, j, done;
        
        
        radix = (sf16) FILT_RADIX;
        sqrdx = radix * radix;
        done = 0;
        while (!done) {
                done = 1;
                for (i = 0; i < poles; i++) {
                        r = c = FILT_ZERO;
                        for (j = 0; j < poles; j++)
                                if (j != i) {
                                        c += FILT_ABS(a[j][i]);
                                        r += FILT_ABS(a[i][j]);
                                }
                        if (c != FILT_ZERO && r != FILT_ZERO) {
                                g = r / radix;
                                f = FILT_ONE;
                                s = c + r;
                                while (c < g) {
                                        f *= radix;
                                        c *= sqrdx;
                                }
                                g = r * radix;
                                while (c > g) {
                                        f /= radix;
                                        c /= sqrdx;
                                }
                                if (((c + r) / f) < ((sf16) 0.95 * s)) {
                                        done = 0;
                                        g = 1.0 / f;
                                        for (j = 0; j < poles; j++)
                                                a[i][j] *= g;
                                        for (j = 0; j < poles; j++)
                                                a[j][i] *= f;
                                }
                        }
                }
        }
      
        
        return;
}

si4	FILT_butter(FILT_PROCESSING_STRUCT *filtps)
{
	si4			i, j, n_fcs, offset, idx, order, poles, is_odd;
        sf8 			*d_num, *d_den;
	sf16			samp_freq, fcs[2], *den, sum_num, sum_den;
	sf16			u[2], pi, half_pi, bw, wn, w, *r, *num, ratio;
        sf16            	**a, **inv_a, **ta1, **ta2, *b, *bt, *c, d, t;
	FILT_LONG_COMPLEX	csum_num, csum_den, cratio, *ckern;
	FILT_LONG_COMPLEX	*p, tc, *eigs, *cden, *rc, *cnum;

	w = wn = bw = 0;

        
	// check input
        switch(filtps->type) {
        	case FILT_LOWPASS_TYPE:
        	case FILT_BANDPASS_TYPE:
        	case FILT_HIGHPASS_TYPE:
        	case FILT_BANDSTOP_TYPE:
                        break;
                default:
                        if (!(MEF_globals->behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
                                fprintf(stderr, "Unrecognized filter type: %d [function %s, line %d]\n", filtps->type, __FUNCTION__, __LINE__);
                                if (MEF_globals->behavior_on_fail & RETURN_ON_FAIL)
                                        (void) fprintf(stderr, "\t=> returning NULL\n\n");
                                else if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL)
                                        (void) fprintf(stderr, "\t=> exiting program\n\n");
                        }
                        if (MEF_globals->behavior_on_fail & RETURN_ON_FAIL)
                                return(-1);
                        else if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL)
                                exit(1);
	}
	samp_freq = (sf16) filtps->sampling_frequency;
	fcs[0] = (sf16) filtps->cutoffs[0];
	n_fcs = ((filtps->type == FILT_LOWPASS_TYPE) || (filtps->type == FILT_HIGHPASS_TYPE)) ? 1 : 2;
	if (n_fcs == 2)
		fcs[1] = (sf16) filtps->cutoffs[1];
        order = filtps->order;
	filtps->poles = poles = n_fcs * order;
	is_odd = order % 2;
	
	// step 1: get analog, pre-warped frequencies
	pi = (sf16) M_PI;
        half_pi = pi / (sf16) 2.0;
        for (i = 0; i < n_fcs; ++i)
                u[i] = (sf16) 4.0 * tanl((pi * fcs[i]) / samp_freq);
	
	// step 2: convert to low-pass prototype estimate
	switch (filtps->type) {
		case FILT_LOWPASS_TYPE:
			wn = u[0];
			break;
		case FILT_BANDPASS_TYPE:
			bw = u[1] - u[0];
			wn = sqrtl(u[0] * u[1]);
			break;
		case FILT_HIGHPASS_TYPE:
			wn = u[0];
			break;
		case FILT_BANDSTOP_TYPE:
			bw = u[1] - u[0];
			wn = sqrtl(u[0] * u[1]);
			break;
	}
	
	// step 3: Get N-th order Butterworth analog lowpass prototype
	p = (FILT_LONG_COMPLEX *) e_calloc((size_t) order, sizeof(FILT_LONG_COMPLEX), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	for (i = 1; i < order; i += 2) {
                p[i - 1].imag = ((pi * (sf16) i) / (sf16) (2 * order)) + half_pi;
                FILT_complex_expl(p + i - 1, p + i - 1);
	}
        for (i = 1; i < order; i += 2) {
                p[i].real = p[i - 1].real;
                p[i].imag = -p[i - 1].imag;
        }
        if (is_odd)
                p[order - 1].real = (sf16) -1.0;
        
        j = order - 1;  // sort into ascending order, by real values
        if (is_odd) --j;
        for (i = 0; j > i; ++i, --j) {
                tc = p[i];
                p[i] = p[j];
                p[j] = tc;
        }
        
        // Transform to state-space
        a = (sf16 **) e_calloc((size_t) poles, sizeof(sf16 *), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        inv_a = (sf16 **) e_calloc((size_t) poles, sizeof(sf16 *), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        ta1 = (sf16 **) e_calloc((size_t) poles, sizeof(sf16 *), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        ta2 = (sf16 **) e_calloc((size_t) poles, sizeof(sf16 *), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        for (i = 0; i < poles; ++i) {
                a[i] = (sf16 *) e_calloc((size_t) poles, sizeof(sf16), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		inv_a[i] = (sf16 *) e_calloc((size_t) poles, sizeof(sf16), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
                ta1[i] = (sf16 *) e_calloc((size_t) poles, sizeof(sf16), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		ta2[i] = (sf16 *) e_calloc((size_t) poles, sizeof(sf16), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	}
        b = (sf16 *) e_calloc((size_t) poles, sizeof(sf16), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        bt = (sf16 *) e_calloc((size_t) poles, sizeof(sf16), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        c = (sf16 *) e_calloc((size_t) poles, sizeof(sf16), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	if ((offset = is_odd))
		a[0][0] = (sf16) -1.0;
	for (i = 0; i < order - 1; i += 2) {
		if ((idx = i + offset))
			a[i + offset][idx - 1] = FILT_ONE;
		a[i + offset][i + offset] = p[i].real + p[i + 1].real;
		a[i + offset][i + offset + 1] = (sf16) -1.0;
		a[i + offset + 1][i + offset] = FILT_ONE;
	}
	b[order - 1] = FILT_ONE;
	c[order - 1] = FILT_ONE;
	d = FILT_ZERO;
	
	// step 4: Transform to lowpass, bandpass, highpass, or bandstop of desired Wn
	switch (filtps->type) {
		case FILT_LOWPASS_TYPE:
			for (i = 0; i < order; ++i) {
				for (j = 0; j < order; ++j)
					a[i][j] *= wn;
				b[i] *= wn;
			}
			break;
		case FILT_BANDPASS_TYPE:
			for (i = 0; i < order; ++i) {
				for (j = 0; j < order; ++j) {
					a[i][j] *= bw;
				}
				a[i][i + order] = wn;
				a[i + order][i] = -wn;
			}
			b[0] = bw;
			break;
		case FILT_HIGHPASS_TYPE:
			for (i = 0; i < order; ++i) {
				c[i] = (sf16) -1.0;
				b[i] = wn;
			}
			for (i = is_odd; i < order; i += 2) {
				c[i + 1] = a[i][i];
				b[i] = FILT_ZERO;
			}
			d = FILT_ONE;
			FILT_invert_matrix(a, inv_a, order);
			
			for (i = 0; i < order; ++i)
				for (j = 0; j < order; ++j)
					a[i][j] = wn * inv_a[i][j];
			break;
		case FILT_BANDSTOP_TYPE:
			for (i = 0; i < order; ++i) {
				c[i] = (sf16) -1.0;
				b[i] = bw;
			}
			for (i = is_odd; i < order; i += 2) {
				c[i + 1] = a[i][i];
				b[i] = FILT_ZERO;
			}
			FILT_invert_matrix(a, inv_a, order);
			for (i = 0; i < order; ++i) {
				for (j = 0; j < order; ++j) {
					a[i][j] = bw * inv_a[i][j];
				}
				a[i][i + order] = wn;
				a[i + order][i] = -wn;
			}
			d = FILT_ONE;
			break;
	}
	
	// step 5: Use bilinear transformation to find discrete equivalent
	t = (sf16) 0.25;
	for (i = 0; i < poles; ++i) {
		for (j = 0; j < poles; ++j) {
			ta1[i][j] = t * a[i][j];
			ta2[i][j] = -ta1[i][j];
		}
		ta1[i][i] += FILT_ONE;
		ta2[i][i] += FILT_ONE;
	}
        
	FILT_invert_matrix(ta2, inv_a, poles);
	
	FILT_mat_multl((void *) inv_a, (void *) ta1, (void *) a, poles, poles, poles);
        
	FILT_mat_multl((void *) c, (void *) inv_a, (void *) bt, 1, poles, poles);
	t = sqrtl((sf16) 0.5);
	for (i = 0; i < poles; ++i)
		c[i] = bt[i] * t;
	
	FILT_mat_multl((void *) bt, (void *) b, (void *) &t, 1, poles, 1);
	d += (t * (sf16) 0.25);
	
	FILT_mat_multl((void *) inv_a, (void *) b, (void *) bt, poles, poles, 1);
	t = FILT_ONE / sqrtl((sf16) 2.0);
	for (i = 0; i < poles; ++i)
		b[i] = bt[i] * t;
        
	// Transform to zero-pole-gain and polynomial forms
	eigs = (FILT_LONG_COMPLEX *) e_calloc((size_t) poles, sizeof(FILT_LONG_COMPLEX), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	FILT_unsymmeig(a, poles, eigs);
	
	den = (sf16 *) e_calloc((size_t) (poles + 1), sizeof(sf16), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	cden = (FILT_LONG_COMPLEX *) e_calloc((size_t) (poles + 1), sizeof(FILT_LONG_COMPLEX), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	cden[0].real = FILT_ONE;
	for (i = 0; i < poles; ++i) {
		for (j = i + 1; j--;) {
			FILT_complex_multl(eigs + i, cden + j, &tc);
			cden[j + 1].real -= tc.real;
			cden[j + 1].imag -= tc.imag;
		}
	}
	for (i = 0; i <= poles; ++i)
		den[i] = cden[i].real;
	
	// generate numerator
	r = (sf16 *) e_calloc((size_t) (poles + 1), sizeof(sf16), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	rc = (FILT_LONG_COMPLEX *) e_calloc((size_t) (poles + 1), sizeof(FILT_LONG_COMPLEX), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	wn = (sf16) 2.0 * atan2l(wn, 4.0);
	
	switch (filtps->type) {
		case FILT_LOWPASS_TYPE:
			for (i = 0; i < poles; ++i)
				r[i] = (sf16) -1.0;
			break;
		case FILT_BANDPASS_TYPE:
			for (i = 0; i < order; ++i) {
				r[i] = FILT_ONE;
				r[i + order] = (sf16) -1.0;
			}
			w = -wn;
			break;
		case FILT_HIGHPASS_TYPE:
			for (i = 0; i < poles; ++i)
				r[i] = FILT_ONE;
			w = -pi;
			break;
		case FILT_BANDSTOP_TYPE:
			tc.real = FILT_ZERO;
			tc.imag = wn;
			FILT_complex_expl(&tc, &tc);
			for (i = 0; i < poles; i += 2) {
				rc[i].real = tc.real;
				rc[i].imag = tc.imag;
				rc[i + 1].real = tc.real;
				rc[i + 1].imag = -tc.imag;
			}
			break;
	}
        
	num = (sf16 *) e_calloc((size_t) (poles + 1), sizeof(sf16), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	cnum = (FILT_LONG_COMPLEX *) e_calloc((size_t) (poles + 1), sizeof(FILT_LONG_COMPLEX), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	if (filtps->type == FILT_BANDSTOP_TYPE) {
		cnum[0].real = FILT_ONE;
		for (i = 0; i < poles; ++i) {
			for (j = i + 1; j--;) {
				FILT_complex_multl(rc + i, cnum + j, &tc);
				cnum[j + 1].real -= tc.real;
				cnum[j + 1].imag -= tc.imag;
			}
		}
		for (i = 0; i <= poles; ++i)
			num[i] = cnum[i].real;
	} else {
		num[0] = FILT_ONE;
		for (i = 0; i < poles; ++i)
			for (j = i + 1; j--;)
				num[j + 1] -= r[i] * num[j];
	}
        
	// normalize
	ckern = (FILT_LONG_COMPLEX *) e_calloc((size_t) (poles + 1), sizeof(FILT_LONG_COMPLEX), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	if ((filtps->type == FILT_LOWPASS_TYPE) || (filtps->type == FILT_BANDSTOP_TYPE)) {
		sum_num = sum_den = FILT_ZERO;
		for (i = 0; i <= poles; ++i) {
			sum_num += num[i];
			sum_den += den[i];
		}
		ratio = sum_den / sum_num;
		for (i = 0; i <= poles; ++i)
			num[i] *= ratio;
	} else {
		tc.real = FILT_ZERO;
		for (i = 0; i <= poles; ++i) {
			tc.imag = w * (sf16) i;
			FILT_complex_expl(&tc, ckern + i);
			cnum[i].real = num[i];
			cnum[i].imag = FILT_ZERO;
			cden[i].real = den[i];
			cden[i].imag = FILT_ZERO;
		}
		csum_num.real = csum_den.real = csum_num.imag = csum_den.imag = FILT_ZERO;
		for (i = 0; i <= poles; ++i) {
			FILT_complex_multl(ckern + i, cnum + i, &tc);
			csum_num.real += tc.real;
			csum_num.imag += tc.imag;
			FILT_complex_multl(ckern + i, cden + i, &tc);
			csum_den.real += tc.real;
			csum_den.imag += tc.imag;
		}
		FILT_complex_divl(&csum_den, &csum_num, &cratio);
		for (i = 0; i <= poles; ++i) {
			FILT_complex_multl(cnum + i, &cratio, &tc);
			num[i] = tc.real;
		}
	}
	
	// set & check output
	d_num = (sf8 *) e_calloc((size_t) (poles + 1), sizeof(sf8), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	d_den = (sf8 *) e_calloc((size_t) (poles + 1), sizeof(sf8), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	for (i = 0; i <= poles; ++i) {
		d_den[i] = (sf8) den[i];
		d_num[i] = (sf8) num[i];
		if (isnan(d_num[i]) || isinf(d_num[i]) || isnan(d_den[i]) || isinf(d_den[i])) {
                        if (!(MEF_globals->behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
                                fprintf(stderr, "Bad filter: [function %s, line %d]\n", __FUNCTION__, __LINE__);
                                if (MEF_globals->behavior_on_fail & RETURN_ON_FAIL)
                                        (void) fprintf(stderr, "\t=> returning FILT_BAD_FILTER\n\n");
                                else if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL)
                                        (void) fprintf(stderr, "\t=> exiting program\n\n");
                        }
                        if (MEF_globals->behavior_on_fail & RETURN_ON_FAIL)
                                return(FILT_BAD_FILTER);
                        else if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL)
                                exit(FILT_BAD_FILTER);

			exit(FILT_BAD_FILTER);
                }
	}
        filtps->numerators = d_num;
        filtps->denominators = d_den;
        
        // clean up
        for (i = 0; i < poles; ++i) {
                free(a[i]);
		free(inv_a[i]);
		free(ta1[i]);
		free(ta2[i]);
	}
        free(a);
	free(inv_a);
	free(ta1);
	free(ta2);
	free(p);
        free(b);
	free(bt);
        free(c);
        free(den);
	free(cden);
	free(num);
	free(cnum);
	free(eigs);
	free(r);
	free(rc);
	free(ckern);
        
        
        return(0);
}


void	FILT_complex_divl(FILT_LONG_COMPLEX *a, FILT_LONG_COMPLEX *b, FILT_LONG_COMPLEX *quotient)  //  returns a / b
{
        FILT_LONG_COMPLEX	ta, tb;
	sf16			den;
        
        
        ta = *a;  // copy in case in place
        tb = *b;
	den = (tb.real * tb.real) + (tb.imag * tb.imag);
        quotient->real = ((ta.real * tb.real) + (ta.imag * tb.imag)) / den;
        quotient->imag = ((ta.imag * tb.real) - (ta.real * tb.imag)) / den;
	
        
        return;
}


void	FILT_complex_expl(FILT_LONG_COMPLEX *exponent, FILT_LONG_COMPLEX *ans)
{
        FILT_LONG_COMPLEX    t;
        sf16            	c;
        
        
        t = *exponent;  // copy in case in place
        c = expl(t.real);
        ans->real = c * cosl(t.imag);
        ans->imag = c * sinl(t.imag);
        
        
        return;
}


void	FILT_complex_multl(FILT_LONG_COMPLEX *a, FILT_LONG_COMPLEX *b, FILT_LONG_COMPLEX *product)
{
        FILT_LONG_COMPLEX    ta, tb;
        
        
        ta = *a;  // copy in case in place
        tb = *b;
        product->real = (ta.real * tb.real) - (ta.imag * tb.imag);
        product->imag = (ta.real * tb.imag) + (ta.imag * tb.real);
        
        
        return;
}


void	FILT_elmhes(sf16 **a, si4 poles)
{
        si4     i, j, m;
        sf16    x, y, t1;
        
        
        for (m = 1; m < (poles - 1); m++) {
                x = FILT_ZERO;
                i = m;
                for (j = m; j < poles; j++) {
                        if (FILT_ABS(a[j][m-1]) > FILT_ABS(x)) {
                                x = a[j][m - 1];
                                i = j;
                        }
                }
                if (i != m) {
                        for (j = m - 1; j < poles; j++) {
                                t1 = a[i][j];
                                a[i][j] = a[m][j];
                                a[m][j] = t1;
                        }
                        for (j = 0; j < poles; j++) {
                                t1 = a[j][i];
                                a[j][i] = a[j][m];
                                a[j][m] = t1;
                        }
                }
                if (x != FILT_ZERO) {
                        for (i = m + 1; i < poles; i++) {
                                y = a[i][m - 1];
                                if (y != FILT_ZERO) {
                                        y /= x;
                                        a[i][m - 1] = y;
                                        for (j = m; j < poles; j++)
                                                a[i][j] -= (y * a[m][j]);
                                        for (j = 0; j < poles; j++)
                                                a[j][m] += (y * a[j][i]);
                                }
                        }
                }
        }
        
        
        return;
}


si4	FILT_filtfilt(FILT_PROCESSING_STRUCT *filtps)
{
        si4	free_z_flag, free_filt_buf_flag, *data;
        si4	padded_data_len, pad_len, pad_lenx2, poles;
        si8	i, j, k, m, data_len;
	sf8	dx2, zc[FILT_MAX_ORDER * 2], t1, t2;
        sf8 	*num, *den, *filt_data, *z, *filt_buf;
        
  
        poles = filtps->poles;
        pad_len = poles * 3;
        pad_lenx2 = pad_len * 2;
	data_len = filtps->data_length;
	if (data_len < pad_lenx2) {
		if (!(MEF_globals->behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
			fprintf(stderr, "At least %d data points for a filter of order %d: [function %s, line %d]\n", pad_lenx2, poles, __FUNCTION__, __LINE__);
			if (MEF_globals->behavior_on_fail & RETURN_ON_FAIL)
				(void) fprintf(stderr, "\t=> returning without filtering\n\n");
			else if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL)
				(void) fprintf(stderr, "\t=> exiting program\n\n");
		}
		if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL)
			exit(FILT_BAD_DATA);
		return(FILT_BAD_DATA);
	}
	num = filtps->numerators;
	den = filtps->denominators;
	data = filtps->orig_data;
        filt_data = filtps->sf8_filt_data;
        if (filtps->initial_conditions == NULL) {
                FILT_generate_initial_conditions(filtps);
                free_z_flag = 1;
        } else
                free_z_flag = 0;
        z = filtps->initial_conditions;
        if (filtps->sf8_buffer == NULL) {
                free_filt_buf_flag = 1;
                filtps->sf8_buffer = (sf8 *) e_calloc((size_t) data_len + pad_lenx2, sizeof(sf8), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        } else
                free_filt_buf_flag = 0;
        filt_buf = filtps->sf8_buffer;
        
        // front pad
        dx2 = (sf8) data[0] * (sf8) 2.0;
        for (i = 0, j = pad_len; j; ++i, --j)
                filt_data[i] = dx2 - (sf8) data[j];
        // copy data
        for (i = 0, j = pad_len; i < data_len; ++i, ++j)
                filt_data[j] = (sf8) data[i];
        // back pad
        padded_data_len = data_len + pad_lenx2;
        dx2 = (sf8) data[data_len - 1] * (sf8) 2.0;
        for (i = data_len + (pad_len), j = data_len - 2; i < padded_data_len; ++i, --j)
                filt_data[i] = dx2 - (sf8) data[j];
        
        // copy and initialize initial conditions
        for (i = 0; i < poles; ++i)
                zc[i] = z[i] * filt_data[0];
	
        // forward filter from filt_data to filt_buf
        for (i = 0; i < padded_data_len; ++i) {
                t1 = filt_data[i];
                t2 = (num[0] * t1) + zc[0];
                for (j = 1; j < poles; ++j) {
                        zc[j - 1] = (num[j] * t1) - (den[j] * t2) + zc[j];
                }
                zc[poles - 1] = (num[poles] * t1) - (den[poles] * t2);
                filt_buf[i] = t2;
        }
	
        // copy and initialize initial conditions
        for (i = 0; i < poles; ++i)
                zc[i] = z[i] * filt_buf[padded_data_len - 1];
	
        // reverse filter from filt_buf to filt_data
	for (i = padded_data_len - 1, k = pad_len; k--;) {
                t1 = filt_buf[i--];
                t2 = (num[0] * t1) + zc[0];
                for (j = 1; j < poles; ++j) {
                        zc[j - 1] = (num[j] * t1) - (den[j] * t2) + zc[j];
                }
                zc[poles - 1] = (num[poles] * t1) - (den[poles] * t2);
        }
        for (m = i - pad_len, k = data_len; k--;) {
                t1 = filt_buf[i--];
                t2 = (num[0] * t1) + zc[0];
                for (j = 1; j < poles; ++j) {
                        zc[j - 1] = (num[j] * t1) - (den[j] * t2) + zc[j];
                }
                zc[poles - 1] = (num[poles] * t1) - (den[poles] * t2);
                filt_data[m--] = t2;
        }
        
	// free as needed
        if (free_z_flag)
                free(z);
        if (free_filt_buf_flag)
                free(filt_buf);
        
        
	return(0);
}


void	FILT_free_processing_struct(FILT_PROCESSING_STRUCT *filtps, si1 free_orig_data, si1 free_filt_data)
{
	if (filtps->numerators != NULL)
		free(filtps->numerators);
	if (filtps->denominators != NULL)
		free(filtps->denominators);
	if (filtps->initial_conditions != NULL)
		free(filtps->initial_conditions);
	if (filtps->orig_data != NULL && free_orig_data == MEF_TRUE)
		free(filtps->orig_data);
	if (filtps->filt_data != NULL && free_filt_data == MEF_TRUE)
		free(filtps->filt_data);
	if (filtps->sf8_filt_data != NULL)
		free(filtps->sf8_filt_data);
	if (filtps->sf8_buffer != NULL)
		free(filtps->sf8_buffer);
        
        free(filtps);
        
        
        return;
}


void	FILT_generate_initial_conditions(FILT_PROCESSING_STRUCT *filtps)
{
        si4     i, j, poles;
        sf16    **q, *rhs, *z;
        sf8	*num, *den;


        poles = filtps->poles;
        num = filtps->numerators;
        den = filtps->denominators;
        q = (sf16 **) e_calloc((size_t) poles, sizeof(sf16 *), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        for (i = 0; i < poles; ++i)
                q[i] = (sf16 *) e_calloc((size_t) poles, sizeof(sf16), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        
        rhs = (sf16 *) e_calloc((size_t) poles, sizeof(sf16), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        z = (sf16 *) e_calloc((size_t) poles, sizeof(sf16), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        filtps->initial_conditions = (sf8 *) e_calloc((size_t) poles, sizeof(sf8), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        
        q[0][0] = (sf16) 1.0 + (sf16) den[1];
        for (i = 1, j = 2; i < poles; ++i, ++j)
                q[i][0] = (sf16) den[j];
        for (i = 1; i < poles; ++i) {
                q[i - 1][i] = (sf16) -1.0;
                q[i][i] = (sf16) 1.0;
        }
        for (i = 0, j = 1; i < poles; ++i, ++j)
                rhs[i] = (sf16) num[j] - ((sf16) num[0] * (sf16) den[j]);
        
        FILT_invert_matrix(q, q, poles);
        FILT_mat_multl(q, rhs, z, poles, poles, 1);
        
        for (i = 0; i < poles; ++i)
                filtps->initial_conditions[i] = (sf8) z[i];
        
        for (i = 0; i < poles; ++i)
                free(q[i]);
	
        free(q);
        free(rhs);
        free(z);
        
        
        return;
}


void	FILT_hqr(sf16 **a, si4 poles, FILT_LONG_COMPLEX *eigs)
{
        si4     nn, m, l, k, j, its, i, mmin, max;
        sf16    z, y, x, w, v, u, t, s, r, q, p, anorm, eps, t1, t2, t3, t4;
	
        
        anorm = FILT_ZERO;
        eps = (sf16) FILT_EPS_SF16;

        p = q = r = 0;
        
        for (i = 0; i < poles; i++) {
                max = ((i - 1) > 0) ? (i - 1) : 0;
                for (j = max; j < poles; j++) {
                        anorm += FILT_ABS(a[i][j]);
                }
        }
        
        nn = poles - 1;
        t = FILT_ZERO;
        while (nn >= 0) {
                its = 0;
                do {
                        for (l = nn; l > 0; l--) {
                                t1 = FILT_ABS(a[l - 1][l - 1]);
                                t2 = FILT_ABS(a[l][l]);
                                s = t1 + t2;
                                if (s == FILT_ZERO)
                                        s = anorm;
                                t1 = FILT_ABS(a[l][l - 1]);
                                if (t1 <= (eps * s)) {
                                        a[l][l - 1] = FILT_ZERO;
                                        break;
                                }
                        }
                        x = a[nn][nn];
                        if (l == nn) {
                                eigs[nn].real = x + t;
                                eigs[nn--].imag = FILT_ZERO;
                        } else {
                                y = a[nn - 1][nn - 1];
                                w = a[nn][nn - 1] * a[nn - 1][nn];
                                if (l == (nn - 1)) {
                                        p = (sf16) 0.5 * (y - x);
                                        q = (p * p) + w;
                                        z = sqrtl(FILT_ABS(q));
                                        x += t;
                                        if (q >= FILT_ZERO) {
                                                t1 = FILT_SIGN(z, p);
                                                z = p + t1;
                                                eigs[nn - 1].real = eigs[nn].real = x + z;
                                                if (z != FILT_ZERO)
                                                        eigs[nn].real = x - w / z;
                                        } else {
                                                eigs[nn].real = x + p;
                                                eigs[nn].imag = -z;
                                                eigs[nn - 1].real = eigs[nn].real;
                                                eigs[nn - 1].imag = -eigs[nn].imag;
                                        }
                                        nn -= 2;
                                } else {
                                        if (its == 30) {
                                                fprintf(stderr, "Too many iterations in hqr\n");
                                                exit(1);
                                        }
                                        if (its == 10 || its == 20) {
                                                t += x;
                                                for (i = 0; i < nn + 1; i++)
                                                        a[i][i] -= x;
                                                t1 = FILT_ABS(a[nn][nn - 1]);
                                                t2 = FILT_ABS(a[nn - 1][nn - 2]);
                                                s = t1 + t2;
                                                y = x = (sf16) 0.75 * s;
                                                w = (sf16) -0.4375 * s * s;
                                        }
                                        ++its;
                                        for (m = nn - 2; m >= l; m--) {
                                                z = a[m][m];
                                                r = x - z;
                                                s = y - z;
                                                p = ((r * s - w) / a[m + 1][m]) + a[m][m + 1];
                                                q = a[m + 1][m + 1] - z - r - s;
                                                r = a[m + 2][m + 1];
                                                t1 = FILT_ABS(p);
                                                t2 = FILT_ABS(q);
                                                t3 = FILT_ABS(r);
                                                s = t1 + t2 + t3;
                                                p /= s;
                                                q /= s;
                                                r /= s;
                                                if (m == l)
                                                        break;
                                                t1 = FILT_ABS(a[m][m - 1]);
                                                t2 = FILT_ABS(q);
                                                t3 = FILT_ABS(r);
                                                u = t1 * (t2 + t3);
                                                t1 = FILT_ABS(p);
                                                t2 = FILT_ABS(a[m - 1][m - 1]);
                                                t3 = FILT_ABS(z);
                                                t4 = FILT_ABS(a[m + 1][m + 1]);
                                                v = t1 * (t2 + t3 + t4);
                                                if (u <= (eps * v))
                                                        break;
                                        }
                                        for (i = m; i < (nn - 1); i++) {
                                                a[i + 2][i] = FILT_ZERO;
                                                if (i != m)
                                                        a[i + 2][i - 1] = FILT_ZERO;
                                        }
                                        for (k = m; k < nn; k++) {
                                                if (k != m) {
                                                        p = a[k][k - 1];
                                                        q = a[k + 1][k - 1];
                                                        r = FILT_ZERO;
                                                        if (k + 1 != nn)
                                                                r = a[k + 2][k - 1];
                                                        t1 = FILT_ABS(p);
                                                        t2 = FILT_ABS(q);
                                                        t3 = FILT_ABS(r);
                                                        if ((x = t1 + t2 + t3) != FILT_ZERO) {
                                                                p /= x;
                                                                q /= x;
                                                                r /= x;
                                                        }
                                                }
                                                t1 = sqrtl((p * p) + (q * q) + (r * r));
                                                s = FILT_SIGN(t1, p);
                                                if (s != 0.0) {
                                                        if (k == m) {
                                                                if (l != m)
                                                                        a[k][k - 1] = -a[k][k - 1];
                                                        } else
                                                                a[k][k - 1] = -s * x;
                                                        p += s;
                                                        x = p / s;
                                                        y = q / s;
                                                        z = r / s;
                                                        q /= p;
                                                        r /= p;
                                                        for (j = k; j < (nn + 1); j++) {
                                                                p = a[k][j] + (q * a[k + 1][j]);
                                                                if ((k + 1) != nn) {
                                                                        p += r * a[k + 2][j];
                                                                        a[k + 2][j] -= p * z;
                                                                }
                                                                a[k + 1][j] -= p * y;
                                                                a[k][j] -= p * x;
                                                        }
                                                        mmin = nn < (k + 3) ? nn : k + 3;
                                                        for (i = l; i < (mmin + 1); i++) {
                                                                p = (x * a[i][k]) + (y * a[i][k + 1]);
                                                                if ((k + 1) != nn) {
                                                                        p += z * a[i][k + 2];
                                                                        a[i][k + 2] -= p * r;
                                                                }
                                                                a[i][k + 1] -= p * q;
                                                                a[i][k] -= p;
                                                        }
                                                }
                                        }
                                }
                        }
                } while ((l + 1) < nn);
        }
        
        
        return;
}


FILT_PROCESSING_STRUCT	*FILT_initialize_processing_struct(si4 order, si4 type, sf8 samp_freq, si8 data_len, si1 alloc_orig_data, si1 alloc_filt_data, sf8 cutoff_1, ...)
{
	si8			buff_len;
	FILT_PROCESSING_STRUCT	*filtps;
	va_list			argp;
	
	
	// NOTE: THIS FUNCTION DOES NOT ALLOCATE THE "data" ARRAY OF THE FILT_PROCESSING_STRUCT
	
	// allocate
	filtps = (FILT_PROCESSING_STRUCT *) e_calloc((size_t) 1, sizeof(FILT_PROCESSING_STRUCT), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	// populate
	filtps->order = order;
	filtps->type = type;
	filtps->sampling_frequency = samp_freq;
	filtps->data_length = data_len;
	filtps->cutoffs[0] = cutoff_1;
	if (type == FILT_BANDPASS_TYPE || type == FILT_BANDSTOP_TYPE) {
		va_start(argp, cutoff_1);
		filtps->cutoffs[1] = va_arg(argp, sf8);
		va_end(argp);
	}
	
	// geneerate coefficients
	FILT_butter(filtps);
	FILT_generate_initial_conditions(filtps);
	
	// allocate
	if (alloc_orig_data == MEF_TRUE)
		filtps->orig_data = (si4 *) e_calloc((size_t) data_len, sizeof(si4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	if (alloc_filt_data == MEF_TRUE)
		filtps->filt_data = (si4 *) e_calloc((size_t) data_len, sizeof(si4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	buff_len = data_len + (6 * filtps->poles);
	filtps->sf8_filt_data = (sf8 *) e_calloc((size_t) buff_len, sizeof(sf8), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	filtps->sf8_buffer = (sf8 *) e_calloc((size_t) buff_len, sizeof(sf8), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	
	return(filtps);
}


void	FILT_invert_matrix(sf16 **a, sf16 **inv_a, si4 order)  // done in place if a == inv_a
{
        si4	*indxc, *indxr, *ipiv;
	si4	i, icol, irow, j, k, l, ll;
	sf16	big, dum, pivinv, temp;
        
    icol = irow = 0;


	indxc = (si4 *) e_calloc((size_t) order, sizeof(si4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	indxr = (si4 *) e_calloc((size_t) order, sizeof(si4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	ipiv = (si4 *) e_calloc((size_t) order, sizeof(si4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        
        if (inv_a != a) {
                for (i = 0; i < order; i++)
                        for (j = 0; j < order; j++)
                                inv_a[i][j] = a[i][j];
        }
	
	for (i = 0; i < order; i++) {
		big = FILT_ZERO;
		for (j = 0; j < order; j++)
			if (ipiv[j] != 1)
				for (k = 0; k < order; k++) {
					if (ipiv[k] == 0) {
						if (FILT_ABS(inv_a[j][k]) >= big) {
							big = FILT_ABS(inv_a[j][k]);
							irow = j;
							icol = k;
						}
					}
				}
		++ipiv[icol];
		if (irow != icol) {
			for (l = 0; l < order; l++) {
                                temp = inv_a[irow][l];
                                inv_a[irow][l] = inv_a[icol][l];
                                inv_a[icol][l] = temp;
                        }
                }
		indxr[i] = irow;
		indxc[i] = icol;
		if (inv_a[icol][icol] == FILT_ZERO) {
                        fprintf(stderr, "invert_matrix: Singular Matrix\n");
                        exit(1);
                }
		pivinv = FILT_ONE / inv_a[icol][icol];
		inv_a[icol][icol] = FILT_ONE;
		for (l = 0; l < order; l++)
                        inv_a[icol][l] *= pivinv;
		for (ll = 0; ll < order; ll++) {
			if (ll != icol) {
				dum = inv_a[ll][icol];
				inv_a[ll][icol] = FILT_ZERO;
				for (l = 0; l < order;l++)
                                        inv_a[ll][l] -= inv_a[icol][l] * dum;
			}
                }
	}
        
	for (l = order - 1;l >= 0; l--) {
		if (indxr[l] != indxc[l]) {
			for (k = 0; k < order; k++) {
                                temp = inv_a[k][indxr[l]];
                                inv_a[k][indxr[l]] = inv_a[k][indxc[l]];
                                inv_a[k][indxc[l]] = temp;
                        }
		}
	}
        
	free(ipiv);
	free(indxr);
	free(indxc);
        
        
        return;
}


void	FILT_mat_multl(void *a, void *b, void *product, si4 outer_dim1, si4 inner_dim, si4 outer_dim2)
{
        si4	i, j, k, v1, v2, vp;
	sf16	sum, t1, t2, *av, **am, *bv, **bm, *pv, **pm;
        
	
	if ((outer_dim1 == 1) || (inner_dim == 1)) {
		av = (sf16 *) a;
		am = NULL;
		v1 = 1;
	} else {
		am = (sf16 **) a;
		av = NULL;
		v1 = 0;
	}
	if ((outer_dim2 == 1) || (inner_dim == 1)) {
		bv = (sf16 *) b;
		bm = NULL;
		v2 = 1;
	} else {
		bm = (sf16 **) b;
		bv = NULL;
		v2 = 0;
	}
	if ((outer_dim1 == 1) || (outer_dim2 == 1)) {
		pv = (sf16 *) product;
		pm = NULL;
		vp = 1;
	} else {
		pm = (sf16 **) product;
		pv = NULL;
		vp = 0;
	}
        
	for (i = 0; i < outer_dim1; ++i) {
		for (j = 0; j < outer_dim2; ++j) {
			sum = 0.0;
			for (k = 0; k < inner_dim; ++k) {
				t1 = (v1) ? av[k] : am[i][k];
				t2 = (v2) ? bv[k] : bm[k][j];
				sum += t1 * t2;
			}
			if (vp) {
				if (outer_dim1 == 1)
					pv[j] = sum;
				else
					pv[i] = sum;
			} else {
				pm[i][j] = sum;
			}
		}
	}
	
        
        return;
}


void	FILT_unsymmeig(sf16 **a, si4 poles, FILT_LONG_COMPLEX *eigs)
{
        FILT_balance(a, poles);
        FILT_elmhes(a, poles);
        FILT_hqr(a, poles, eigs);
	
        
        return;
}


/*************************************************************************/
/**************************  END FILTER FUNCTIONS  ***********************/
/*************************************************************************/


si8	*find_discontinuity_indices(TIME_SERIES_INDEX *tsi, si8 num_disconts, si8 number_of_blocks)
{
	si8	i, j, *disconts;
	
	
	disconts = (si8 *) malloc(num_disconts * sizeof(si8));
	for (i = j = 0; i < number_of_blocks; ++i) {
	    if (tsi[i].RED_block_flags & RED_DISCONTINUITY_MASK) {
		    disconts[j++] = i;
	    }
	}
	
	
	return(disconts);
}


si8	*find_discontinuity_samples(TIME_SERIES_INDEX *tsi, si8 num_disconts, si8 number_of_blocks, si1 add_tail)
{
	si8				i, *discont_inds, *discont_samps;
	
	
	// get indices of discontinuity blocks
	discont_inds = find_discontinuity_indices(tsi, num_disconts, number_of_blocks);
	
	// get sample numbers
	discont_samps = (si8 *) e_calloc((size_t) (num_disconts + 1), sizeof(si8), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	for (i = 0; i < num_disconts; ++i)
		discont_samps[i] = tsi[discont_inds[i]].start_sample;
	
	// clean up
	free(discont_inds);
	
	// add a tail
	if (add_tail == MEF_TRUE)
		discont_samps[num_disconts] = discont_samps[num_disconts - 1] + tsi[num_disconts - 1].number_of_samples;
	
	
	return(discont_samps);
}


/*** THIS ROUTINE IS NOT THREAD SAFE - USE ONLY IN UNTHREADED APPLICATIONS ***/
void	force_behavior(ui4 behavior)
{
	static ui4	saved_behavior = MEF_GLOBALS_BEHAVIOR_ON_FAIL_DEFAULT;
	
	
	
	if (behavior == RESTORE_BEHAVIOR) {
		MEF_globals->behavior_on_fail = saved_behavior;
		return;
	}
	
	saved_behavior = MEF_globals->behavior_on_fail;
	MEF_globals->behavior_on_fail = behavior;
	
	
	return;
}


/*************************************************************************/
/**************  FILE PROCESSING STRUCT STANDARD FUNCTIONS  **************/
/*************************************************************************/


void	fps_close(FILE_PROCESSING_STRUCT *fps) {
	
	fclose(fps->fp);
	fps->fp = NULL;
	fps->fd = -1;
	
	return;
}

#ifndef _WIN32
	si4	fps_lock(FILE_PROCESSING_STRUCT *fps, si4 lock_type, const si1 *function, si4 line, ui4 behavior_on_fail)
	{
		struct flock	fl;
		
		
		if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
			behavior_on_fail = MEF_globals->behavior_on_fail;
		
		fl.l_type = lock_type;
		fl.l_whence = SEEK_SET;
		fl.l_start = 0;
		fl.l_len = 0;
		fl.l_pid = getpid();
		if (fcntl(fps->fd, F_SETLKW, &fl) == -1) {
			if (!(behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
				(void) fprintf(stderr, "%c\n\tfcntl() failed to lock file\n", 7);
				(void) fprintf(stderr, "\tsystem error number %d (%s)\n", errno, strerror(errno));
				if (function != NULL)
					(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
				if (behavior_on_fail & RETURN_ON_FAIL)
					(void) fprintf(stderr, "\t=> returning -1\n\n");
				else if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL)
					(void) fprintf(stderr, "\t=> exiting program\n\n");
			}
			if (behavior_on_fail & RETURN_ON_FAIL)
				return(-1);
			else if (behavior_on_fail & EXIT_ON_FAIL)
				exit(1);
		}
		
		
		return(0);
	}
#endif

si4	fps_open(FILE_PROCESSING_STRUCT *fps, const si1 *function, si4 line, ui4 behavior_on_fail)
{
	si1		*mode, path[MEF_FULL_FILE_NAME_BYTES], command[MEF_FULL_FILE_NAME_BYTES + 16];
	si1		name[MEF_BASE_FILE_NAME_BYTES], extension[TYPE_BYTES];
	#ifndef _WIN32
		si4		lock_type;
	#endif
	struct stat	sb;
	
	
	if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
		behavior_on_fail = MEF_globals->behavior_on_fail;
	
	// open
	mode = NULL;
	switch(fps->directives.open_mode) {
		case FPS_R_OPEN_MODE:
			mode = "rb";
			break;
		case FPS_R_PLUS_OPEN_MODE:
			mode = "rb+";
			break;
		case FPS_W_OPEN_MODE:
			mode = "wb";
			break;
		case FPS_W_PLUS_OPEN_MODE:
			mode = "wb+";
			break;
		case FPS_A_OPEN_MODE:
			mode = "ab";
			break;
		case FPS_A_PLUS_OPEN_MODE:
			mode = "ab+";
			break;
		case FPS_NO_OPEN_MODE:
		default:
			if (!(behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
				(void) fprintf(stderr, "%c\n\t%s(): invalid open mode (%u)\n", 7, __FUNCTION__, fps->directives.open_mode);
				if (function != NULL)
					(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
				if (behavior_on_fail & RETURN_ON_FAIL)
					(void) fprintf(stderr, "\t=> returning -1\n\n");
				else if (behavior_on_fail & EXIT_ON_FAIL)
					(void) fprintf(stderr, "\t=> exiting program\n\n");
			}
			if (behavior_on_fail & RETURN_ON_FAIL)
				return(-1);
			else if (behavior_on_fail & EXIT_ON_FAIL)
				exit(1);
	}
	
	fps->fp = e_fopen(fps->full_file_name, mode, function, line, RETURN_ON_FAIL | SUPPRESS_ERROR_OUTPUT);
	if (fps->fp == NULL && errno == ENOENT) {
		// A component of the required directory tree does not exist - build it & try again
		extract_path_parts(fps->full_file_name, path, name, extension);
		#ifdef _WIN32
			sprintf(command, "if not exist \"%s\" mkdir \"%s\"", path, path);
		#else
			sprintf(command, "mkdir -p \"%s\"", path);
		#endif
		(void) e_system(command, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		// system(command);
		fps->fp = e_fopen(fps->full_file_name, mode, function, line, behavior_on_fail);
	}
	if (fps->fp == NULL) {
		if (!(behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
			(void) UTF8_fprintf(stderr, "%c\n\t%s() failed to open file \"%s\"\n", 7, __FUNCTION__, fps->full_file_name);
			if (function != NULL)
				(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
			if (behavior_on_fail & RETURN_ON_FAIL)
				(void) fprintf(stderr, "\t=> returning -1\n\n");
			else if (behavior_on_fail & EXIT_ON_FAIL)
				(void) fprintf(stderr, "\t=> exiting program\n\n");
		}
		if (behavior_on_fail & RETURN_ON_FAIL)
			return(-1);
		else if (behavior_on_fail & EXIT_ON_FAIL)
			exit(1);
	}
	
	// get file descriptor
	fps->fd = fileno(fps->fp);
	
	#ifndef _WIN32
		// lock
		if (fps->directives.lock_mode != FPS_NO_LOCK_MODE) {
			lock_type = FPS_NO_LOCK_TYPE;
			if (fps->directives.open_mode == FPS_R_OPEN_MODE) {
				if (fps->directives.lock_mode & FPS_READ_LOCK_ON_READ_OPEN)
					lock_type = F_RDLCK;
				else if (fps->directives.lock_mode & FPS_WRITE_LOCK_ON_READ_OPEN)
					lock_type = F_WRLCK;
			} else if (fps->directives.lock_mode & (FPS_WRITE_LOCK_ON_WRITE_OPEN | FPS_WRITE_LOCK_ON_READ_WRITE_OPEN)) {
				lock_type = F_WRLCK;
			} else {
				if (!(MEF_globals->behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
					(void) fprintf(stderr, "%c\n\t%s(): incompatible lock (%u) and open (%u) modes\n", 7, __FUNCTION__, fps->directives.lock_mode, fps->directives.open_mode);
					if (function != NULL)
						(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
					if (behavior_on_fail & RETURN_ON_FAIL)
						(void) fprintf(stderr, "\t=> returning -1\n\n");
					else if (behavior_on_fail & EXIT_ON_FAIL)
						(void) fprintf(stderr, "\t=> exiting program\n\n");
				}
				if (behavior_on_fail & RETURN_ON_FAIL)
					return(-1);
				else if (behavior_on_fail & EXIT_ON_FAIL)
					exit(1);
			}
			fps_lock(fps, lock_type, function, line, behavior_on_fail);
		}
	#endif
	
	// get file length
	fstat(fps->fd, &sb);
	fps->file_length = sb.st_size;
	
	
	return(0);
}


si4	fps_read(FILE_PROCESSING_STRUCT *fps, const si1 *function, si4 line, ui4 behavior_on_fail)
{
	si8		i_bytes;
	
	
	if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
		behavior_on_fail = MEF_globals->behavior_on_fail;
	
	#ifndef _WIN32
		// lock
		if (fps->directives.lock_mode & FPS_READ_LOCK_ON_READ)
			fps_lock(fps, F_RDLCK, function, line, behavior_on_fail);
	#endif
	// read
	if (fps->directives.io_bytes == FPS_FULL_FILE)
		i_bytes = fps->file_length;
	else
		i_bytes = fps->directives.io_bytes;
	(void) e_fread(fps->raw_data, sizeof(ui1), (size_t) i_bytes, fps->fp, fps->full_file_name, __FUNCTION__, __LINE__, behavior_on_fail);
	
	#ifndef _WIN32
		// unlock
		if (fps->directives.lock_mode & FPS_READ_LOCK_ON_READ)
			fps_unlock(fps, function, line, behavior_on_fail);
	#endif
	
	return(0);
}

#ifndef _WIN32
	si4	fps_unlock(FILE_PROCESSING_STRUCT *fps, const si1 *function, si4 line, ui4 behavior_on_fail)
	{
		struct flock	fl;
		
		
		if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
			behavior_on_fail = MEF_globals->behavior_on_fail;
		
		fl.l_type = F_UNLCK;
		fl.l_whence = SEEK_SET;
		fl.l_start = 0;
		fl.l_len = 0;
		fl.l_pid = getpid();
		if (fcntl(fps->fd, F_SETLKW, &fl) == -1) {
			if (!(behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
				(void) fprintf(stderr, "%c\n\tfcntl() failed to unlock file\n", 7);
				(void) fprintf(stderr, "\tsystem error number %d (%s)\n", errno, strerror(errno));
				if (function != NULL)
					(void) fprintf(stderr, "\tcalled from function \"%s\", line %d\n", function, line);
				if (behavior_on_fail & RETURN_ON_FAIL)
					(void) fprintf(stderr, "\t=> returning -1\n\n");
				else if (behavior_on_fail & EXIT_ON_FAIL)
					(void) fprintf(stderr, "\t=> exiting program\n\n");
			}
			if (behavior_on_fail & RETURN_ON_FAIL)
				return(-1);
			else if (behavior_on_fail & EXIT_ON_FAIL)
				exit(1);
		}
		
		
		return(0);
	}
#endif

si4	fps_write(FILE_PROCESSING_STRUCT *fps, const si1 *function, si4 line, ui4 behavior_on_fail)
{
	si8		o_bytes;
	struct stat	sb;
	
        
	if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
		behavior_on_fail = MEF_globals->behavior_on_fail;
	
	#ifndef _WIN32
		// lock
		if (fps->directives.lock_mode & FPS_WRITE_LOCK_ON_WRITE)
			fps_lock(fps, F_WRLCK, function, line, behavior_on_fail);
	#endif
	
	// write
	if (fps->directives.io_bytes == FPS_FULL_FILE)
		o_bytes = fps->raw_data_bytes;
	else
		o_bytes = fps->directives.io_bytes;
	(void) e_fwrite(fps->raw_data, sizeof(ui1), (size_t) o_bytes, fps->fp, fps->full_file_name, __FUNCTION__, __LINE__, behavior_on_fail);
	
	#ifndef _WIN32
		// unlock
		if (fps->directives.lock_mode & FPS_WRITE_LOCK_ON_WRITE)
			fps_unlock(fps, function, line, behavior_on_fail);
	#endif
        
	// get file length
	fflush(fps->fp);  // have to flush to update stat structure after write
	fstat(fps->fd, &sb);
	fps->file_length = sb.st_size;
	
	
	return(0);
}


/*************************************************************************/
/************  END FILE PROCESSING STRUCT STANDARD FUNCTIONS  ************/
/*************************************************************************/


void	free_channel(CHANNEL *channel, si4 free_channel_structure)
{
        si4	i;
        
        
        for (i = 0; i < channel->number_of_segments; ++i)
                free_segment(channel->segments + i, MEF_FALSE);
        free(channel->segments);
	
        free(channel->metadata.section_1);
	if (channel->metadata.time_series_section_2 != NULL)
		free(channel->metadata.time_series_section_2);
	if (channel->metadata.video_section_2 != NULL)
		free(channel->metadata.video_section_2);
	free(channel->metadata.section_3);
	
	if (channel->record_data_fps != NULL)
		free_file_processing_struct(channel->record_data_fps);
	if (channel->record_indices_fps != NULL)
		free_file_processing_struct(channel->record_indices_fps);
	
        if (free_channel_structure == MEF_TRUE)
		free(channel);
        
        
        return;
}


void	free_file_processing_struct(FILE_PROCESSING_STRUCT *fps)
{
        if (fps == NULL) {
                if (!(MEF_globals->behavior_on_fail & SUPPRESS_ERROR_OUTPUT))
                	fprintf(stderr, "Warning: trying to free a NULL FILE_PROCESSING_STRUCT => returning with no action\n");
                return;
        }
        
	if (fps->password_data != NULL && fps->directives.free_password_data == MEF_TRUE)
                free(fps->password_data);
	
        if (fps->raw_data != NULL && fps->raw_data_bytes > 0)
                free(fps->raw_data);
        
	if (fps->fp != NULL && fps->directives.close_file == MEF_TRUE)
		(void) fclose(fps->fp);
        
        free(fps);
        
        
        return;
}


void	free_segment(SEGMENT *segment, si4 free_segment_structure)
{
	free_file_processing_struct(segment->metadata_fps);
	if (segment->time_series_data_fps != NULL) {
		segment->time_series_data_fps->directives.close_file = MEF_TRUE;
		free_file_processing_struct(segment->time_series_data_fps);
	}
	if (segment->time_series_indices_fps != NULL)
		free_file_processing_struct(segment->time_series_indices_fps);
	if (segment->video_indices_fps != NULL)
		free_file_processing_struct(segment->video_indices_fps);
	if (segment->record_data_fps != NULL)
		free_file_processing_struct(segment->record_data_fps);
	if (segment->record_indices_fps != NULL)
		free_file_processing_struct(segment->record_indices_fps);
	
        if (free_segment_structure == MEF_TRUE)
                free(segment);
        
        
        return;
}


void	free_session(SESSION *session, si4 free_session_structure)
{
        si4	i;
        

	if (session->number_of_time_series_channels > 0) {
		free(session->time_series_metadata.section_1);
		free(session->time_series_metadata.time_series_section_2);
		free(session->time_series_metadata.section_3);
		for (i = 0; i < session->number_of_time_series_channels; ++i)
			free_channel(session->time_series_channels + i, MEF_FALSE);
		free(session->time_series_channels);
	}

	if (session->number_of_video_channels > 0) {
		free(session->video_metadata.section_1);
		free(session->video_metadata.video_section_2);
		free(session->video_metadata.section_3);
		for (i = 0; i < session->number_of_video_channels; ++i)
			free_channel(session->video_channels + i, MEF_FALSE);
		free(session->video_channels);
	}
	
	if (session->record_data_fps != NULL)
		free_file_processing_struct(session->record_data_fps);
	if (session->record_indices_fps != NULL)
		free_file_processing_struct(session->record_indices_fps);
	
        if (free_session_structure == MEF_TRUE)
		free(session);
        
        
        return;
}

#ifdef _WIN32
	si1	**generate_file_list(si1 **file_list, si4 *num_files, si1 *enclosing_directory, si1 *extension)  // can be used to get a directory list also
	{
		si4 i, j;
		si1 temp_str[MEF_FULL_FILE_NAME_BYTES];
		si1 *ext;
        struct _stat sb;
        si4 skip_segment;
        si1 temp_path[MEF_FULL_FILE_NAME_BYTES], temp_name[MEF_SEGMENT_BASE_FILE_NAME_BYTES], temp_extension[TYPE_BYTES];

		// Windows structures
        WIN32_FIND_DATA fdFile; 
	    HANDLE hFind = NULL; 

	    // free previous file list
		if (file_list != NULL) {
			for (i = 0; i < *num_files; ++i)
				free(file_list[i]);
			free(file_list);
		}

	    // get the files / directoris with required extension and count by building a mask
	    *num_files = 0;
	    sprintf(temp_path, "%s\\*.*", enclosing_directory); 
	    if((hFind = FindFirstFile(temp_path, &fdFile)) == INVALID_HANDLE_VALUE) 
	    { 
	        (void) UTF8_fprintf(stderr, "%c\n\t%s() failed to open directory \"%s\"\n", 7, __FUNCTION__, enclosing_directory);
			return 0;
	    } 

	    while (FindNextFile(hFind, &fdFile)){
	    	//Skip initial "." and ".." directories
	        if((strcmp((si1 *) fdFile.cFileName, ".") != 0) && (strcmp((si1 *)fdFile.cFileName, "..") != 0))
	        { 
	        	sprintf(temp_path, "%s\\%s", enclosing_directory, (si1 *) fdFile.cFileName); 

	        	// Get extension
	        	ext = strrchr((si1 *) fdFile.cFileName, '.');
				if (ext != NULL && strlen(ext) != 1)
					ext++;
				
	        	if (!((ext == NULL) || (ext == (si1 *) fdFile.cFileName)) &&  (!strcmp(ext, extension)))
					++(*num_files);
	        }
	    }

		// clean up
	    FindClose(hFind);
		
	    // now read again and allocate and build
	    sprintf(temp_path, "%s\\*.*", enclosing_directory);
		hFind = FindFirstFile(temp_path, &fdFile);
		if ( file_list == NULL ){
			file_list = (si1 **) e_calloc((size_t) *num_files, sizeof(si1 *), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
			i = 0;
			while (FindNextFile(hFind, &fdFile)) {
				// Get extension

				if((strcmp((si1 *) fdFile.cFileName, ".") != 0) && (strcmp((si1 *)fdFile.cFileName, "..") != 0)){

					ext = strrchr((si1 *) fdFile.cFileName, '.');
					if (ext != NULL && strlen(ext) != 1)
						ext++;
					
					if (!((ext == NULL) || (ext == (si1 *) fdFile.cFileName)) &&  (!strcmp(ext, extension))){
						
						file_list[i] = (si1 *) e_malloc((size_t) MEF_FULL_FILE_NAME_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
						
						MEF_strcpy(temp_str, enclosing_directory);
						
                        MEF_strcat(temp_str, "/");

						MEF_strcat(temp_str, (si1 *) fdFile.cFileName);
						MEF_strncpy(file_list[i], temp_str, MEF_FULL_FILE_NAME_BYTES);
						memset(temp_str, 0, MEF_FULL_FILE_NAME_BYTES);

						// check for empty segment situation
			            skip_segment = 0;
			            if (!strcmp(extension, SEGMENT_DIRECTORY_TYPE_STRING))
			            {

			                extract_path_parts(file_list[i], temp_path, temp_name, temp_extension);
			                sprintf(temp_str, "%s/%s.tdat", file_list[i], temp_name);

			                // Chekf if the channel is time series
			                extract_path_parts(temp_path, NULL, NULL, temp_extension);
			                if (!strcmp(temp_extension, TIME_SERIES_CHANNEL_DIRECTORY_TYPE_STRING)){

				                // get file length
				                _stat(temp_str, &sb);

				                if (sb.st_size <= UNIVERSAL_HEADER_BYTES)
				                {
				                    skip_segment = 1;
				                    (*num_files)--;

				                    // normally calling function will free memeory, but here we have to do
				                    // it since calling function won't know about this entry
				                    free (file_list[*num_files]);
				                }
				            }
			            }
			            if (skip_segment == 0)
							++i;
					}
				}
			}
		}

		// clean up
	    FindClose(hFind);

	    // Finally sort file_list alphabetically
		for(i = 0; i <= *num_files; ++i){
			for(j = i+1; j <= *num_files - 1; ++j){
			  	if(strcmp(file_list[i], file_list[j])>0){
			        strcpy(temp_str, file_list[i]);
			        strcpy(file_list[i], file_list[j]);
			        strcpy(file_list[j], temp_str);
			    }
			}
		}

	    return(file_list);
	}
#else
	si1	**generate_file_list(si1 **file_list, si4 *num_files, si1 *enclosing_directory, si1 *extension)  // can be used to get a directory list also
	{
		si4	i, n_entries, n;
		si1	temp_str[MEF_FULL_FILE_NAME_BYTES];
		struct dirent **contents_list;
		si1 *ext;
        struct stat sb;
        si4 skip_segment;
        si1 temp_path[MEF_FULL_FILE_NAME_BYTES], temp_name[MEF_SEGMENT_BASE_FILE_NAME_BYTES], temp_extension[TYPE_BYTES];


		// free previous file list
		if (file_list != NULL) {
			for (i = 0; i < *num_files; ++i)
				free(file_list[i]);
			free(file_list);
		}
		
		// get the files / directoris with required extension and count
		*num_files = 0;
		n_entries = scandir(enclosing_directory, &contents_list, NULL, alphasort);
		n = n_entries;

		if (n < 0){
       		(void) UTF8_fprintf(stderr, "%c\n\t%s() failed to open directory \"%s\"\n", 7, __FUNCTION__, enclosing_directory);
			return 0;}
		n = 0;
		while (n < n_entries) {
			
			// Get extension
			ext = strrchr(contents_list[n]->d_name, '.');
			if (ext != NULL && strlen(ext) != 1)
				ext++;
			
			if (!((ext == NULL) || (ext == contents_list[n]->d_name)) &&  (!strcmp(ext, extension)))
				++(*num_files);
			++n;
	    }

		// now read again and allocate and build
		n = 0;
		if ( file_list == NULL ){
			file_list = (si1 **) e_calloc((size_t) *num_files, sizeof(si1 *), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
			i = 0;

			while (n < n_entries) {
				ext = strrchr(contents_list[n]->d_name, '.');
				if (ext != NULL && strlen(ext) != 1)
					ext++;
				
				if (!((ext == NULL) || (ext == contents_list[n]->d_name)) &&  (!strcmp(ext, extension))){
					file_list[i] = (si1 *) e_malloc((size_t) MEF_FULL_FILE_NAME_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
					MEF_strcpy(temp_str, enclosing_directory);
					
                    MEF_strcat(temp_str, "/");

					MEF_strcat(temp_str, contents_list[n]->d_name);
					MEF_strncpy(file_list[i], temp_str, MEF_FULL_FILE_NAME_BYTES);
					memset(temp_str, 0, MEF_FULL_FILE_NAME_BYTES);

					// check for empty segment situation
		            skip_segment = 0;
		            if (!strcmp(extension, SEGMENT_DIRECTORY_TYPE_STRING))
		            {
		                extract_path_parts(file_list[i], temp_path, temp_name, temp_extension);
		                sprintf(temp_str, "%s/%s.tdat", file_list[i], temp_name);

		                // Chekf if the channel is time series
		                extract_path_parts(temp_path, NULL, NULL, temp_extension);
		                if (!strcmp(temp_extension, TIME_SERIES_CHANNEL_DIRECTORY_TYPE_STRING)){

			                // get file length
			                stat(temp_str, &sb);

			                if (sb.st_size <= UNIVERSAL_HEADER_BYTES)
			                {
			                    skip_segment = 1;
			                    (*num_files)--;

			                    // normally calling function will free memeory, but here we have to do
			                    // it since calling function won't know about this entry
			                    free (file_list[*num_files]);
			                }
			            }
		            }
		            if (skip_segment == 0)
						++i;
				}
				free(contents_list[n]);
				++n;
			}
			free(contents_list);
		}
		
		return(file_list);
	}
#endif

si1	*generate_hex_string(ui1 *bytes, si4 num_bytes, si1 *string)
{
	si4	i;
	si1	*s;
	
	
	if (string == NULL)  // allocate if NULL is passed
		string = (si1 *) e_calloc((size_t) ((num_bytes + 1) * 3), sizeof(si1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	s = string;
	*s++ = '0';
	*s++ = 'x';
        
	for (i = 0; i < num_bytes; ++i) {
		sprintf(s, " %2x", bytes[i]);
		if (bytes[i] < 0x10)
			*(s + 1) = '0';
		s += 3;
	}
	*s = 0;
	
	
	return(string);
}


si8	generate_recording_time_offset(si8 recording_start_time_uutc, si4 GMT_offset)
{
	time_t		recording_start_time_utc;
	struct tm	time_info;
	si8		recording_time_offset_utc, recording_time_offset_uutc;
	
	
	if (recording_start_time_uutc == RTO_USE_SYSTEM_TIME) { // use current system time and time zone
		recording_start_time_utc = time(NULL);
		(void) localtime_r(&recording_start_time_utc, &time_info);
		#ifdef _WIN32
			GMT_offset = -_timezone;
		#else
			GMT_offset = time_info.tm_gmtoff;
		#endif
	} else {
                recording_start_time_utc = recording_start_time_uutc / (si8) 1e6;
                if (GMT_offset > MAXIMUM_GMT_OFFSET || GMT_offset < MINIMUM_GMT_OFFSET) {
                        fprintf(stderr, "%c\n%s(): invalid GMT offset => using 0\n", 7, __FUNCTION__);
                        GMT_offset = 0;
                }
	}
        
        (void) gmtime_r(&recording_start_time_utc, &time_info);
	time_info.tm_sec = 0;
	time_info.tm_min = 0;
	time_info.tm_hour = 0;
        
        recording_time_offset_utc = timegm(&time_info) + (86400 - GMT_offset);  // 86400 = 24 hr in secs
        while (recording_time_offset_utc > recording_start_time_utc)
                recording_time_offset_utc -= 86400;  // 24 hr
        
        recording_time_offset_uutc = recording_time_offset_utc * (si8) 1e6;
	
        MEF_globals->recording_time_offset = recording_time_offset_uutc;
        MEF_globals->GMT_offset = GMT_offset;
        
        if (MEF_globals->verbose == MEF_TRUE) {
		#ifdef _WIN32
			printf("Recording Time Offset = %lld\n", recording_time_offset_uutc);
		#else
			printf("Recording Time Offset = %ld\n", recording_time_offset_uutc);
		#endif
		printf("GMT Offset = %d\n", GMT_offset);
	}

	
	return(recording_time_offset_uutc);
}


si1	*generate_segment_name(FILE_PROCESSING_STRUCT *fps, si1 *segment_name)
{
	si1	segment_number_str[FILE_NUMBERING_DIGITS + 1];
	
	
	if (segment_name == NULL)  // if NULL is passed, this will be allocated, but the calling function has the responsibility to free it.
		segment_name = (si1 *) e_malloc((size_t) MEF_SEGMENT_BASE_FILE_NAME_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	numerical_fixed_width_string(segment_number_str, FILE_NUMBERING_DIGITS, fps->universal_header->segment_number);
	
	MEF_snprintf(segment_name, MEF_SEGMENT_BASE_FILE_NAME_BYTES, "%s-%s", fps->universal_header->channel_name, segment_number_str);

	
	return(segment_name);
}


ui1	*generate_UUID(ui1 *uuid)
{
	si4	i, zero_count;
	
	
	// allocate here if not passed, but caller must free
	if (uuid == NULL)
		uuid = (ui1 *) e_malloc((size_t) UUID_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	do {
		for (i = zero_count = 0; i < UUID_BYTES; ++i) {
			uuid[i] = (ui1) (random() % 0x100);
			if (!uuid[i])
				++zero_count;
		}
	} while (zero_count == UUID_BYTES);
	
        
	return(uuid);
}


/*************************************************************************/
/******************************  MEF GLOBALS  ****************************/
/*************************************************************************/


void	initialize_MEF_globals()
{
	if (MEF_globals == NULL)
		MEF_globals = (MEF_GLOBALS *) e_calloc((size_t) 1, sizeof(MEF_GLOBALS), __FUNCTION__, __LINE__, EXIT_ON_FAIL);
	
	// time constants
	MEF_globals->recording_time_offset = MEF_GLOBALS_RECORDING_TIME_OFFSET_DEFAULT;
	MEF_globals->recording_time_offset_mode = MEF_GLOBALS_RECORDING_TIME_OFFSET_MODE_DEFAULT;
	MEF_globals->GMT_offset = MEF_GLOBALS_GMT_OFFSET_DEFAULT;
        MEF_globals->DST_start_time = MEF_GLOBALS_DST_START_TIME_DEFAULT;
        MEF_globals->DST_end_time = MEF_GLOBALS_DST_END_TIME_DEFAULT;
	// alignment fields
	MEF_globals->universal_header_aligned = MEF_UNKNOWN;
	MEF_globals->metadata_section_1_aligned = MEF_UNKNOWN;
	MEF_globals->time_series_metadata_section_2_aligned = MEF_UNKNOWN;
	MEF_globals->video_metadata_section_2_aligned = MEF_UNKNOWN;
	MEF_globals->metadata_section_3_aligned = MEF_UNKNOWN;
	MEF_globals->all_metadata_structures_aligned = MEF_UNKNOWN;
	MEF_globals->time_series_indices_aligned = MEF_UNKNOWN;
	MEF_globals->video_indices_aligned = MEF_UNKNOWN;
	MEF_globals->RED_block_header_aligned = MEF_UNKNOWN;
	MEF_globals->record_header_aligned = MEF_UNKNOWN;
	MEF_globals->record_indices_aligned = MEF_UNKNOWN;
	MEF_globals->all_record_structures_aligned = MEF_UNKNOWN;
	MEF_globals->all_structures_aligned = MEF_UNKNOWN;
	// RED
	MEF_globals->RED_normal_CDF_table = NULL;
	// CRC
	MEF_globals->CRC_table = NULL;
	MEF_globals->CRC_mode = MEF_GLOBALS_CRC_MODE_DEFAULT;
	// AES
	MEF_globals->AES_sbox_table = NULL;
	MEF_globals->AES_rsbox_table = NULL;
	MEF_globals->AES_rcon_table = NULL;
	// SHA256
	MEF_globals->SHA256_h0_table = NULL;
	MEF_globals->SHA256_k_table = NULL;
	 // UTF8
	MEF_globals->UTF8_offsets_from_UTF8_table = NULL;
	MEF_globals->UTF8_trailing_bytes_for_UTF8_table = NULL;
	// miscellaneous
	MEF_globals->verbose = MEF_GLOBALS_VERBOSE_DEFAULT;
        MEF_globals->behavior_on_fail = MEF_GLOBALS_BEHAVIOR_ON_FAIL_DEFAULT;
        #ifndef _WIN32
		MEF_globals->file_creation_umask = MEF_GLOBALS_FILE_CREATION_UMASK_DEFAULT;
	#endif
	
	
	return;
}


/*************************************************************************/
/****************************  END MEF GLOBALS  **************************/
/*************************************************************************/


FILE_PROCESSING_DIRECTIVES *initialize_file_processing_directives(FILE_PROCESSING_DIRECTIVES *directives)
{
	if (directives == NULL)
		directives = (FILE_PROCESSING_DIRECTIVES *) e_calloc((size_t) 1, sizeof(FILE_PROCESSING_DIRECTIVES), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	// set directives to defaults
	directives->close_file = FPS_DIRECTIVE_CLOSE_FILE_DEFAULT;
	directives->free_password_data = FPS_DIRECTIVE_FREE_PASSWORD_DATA_DEFAULT;
        directives->io_bytes = FPS_DIRECTIVE_IO_BYTES_DEFAULT;
	directives->lock_mode = FPS_DIRECTIVE_LOCK_MODE_DEFAULT;
	directives->open_mode = FPS_DIRECTIVE_OPEN_MODE_DEFAULT;
	
	
	return(directives);
}


si4	initialize_meflib()
{
	si4	return_value;
	
	
	// set up globals
	initialize_MEF_globals();
	
	// check endianess
	if (cpu_endianness() != MEF_LITTLE_ENDIAN) {
		fprintf(stderr, "Error: Library only coded for little-endian machines currently => exiting [function \"%s\", line %d]\n", __FUNCTION__, __LINE__);
		exit(-1);
	}
	
	// check structure alignments
	return_value = check_all_alignments(__FUNCTION__, __LINE__);
	
	// seed random number generator
	srandom((ui4) time(NULL));
	
	// set file creation umask
	umask(MEF_globals->file_creation_umask);
	
	// make RED table global
	(void) RED_initialize_normal_CDF_table(MEF_TRUE);
	
	// make CRC table global
	(void) CRC_initialize_table(MEF_TRUE);
	
	// make UTF-8 tables global
	(void) UTF8_initialize_offsets_from_UTF8_table(MEF_TRUE);
	(void) UTF8_initialize_trailing_bytes_for_UTF8_table(MEF_TRUE);
	
	// make AES-128 tables global
	(void) AES_initialize_sbox_table(MEF_TRUE);
	(void) AES_initialize_rsbox_table(MEF_TRUE);
	(void) AES_initialize_rcon_table(MEF_TRUE);
	
	// make SHA-256 tables global
	(void) SHA256_initialize_h0_table(MEF_TRUE);
	(void) SHA256_initialize_k_table(MEF_TRUE);
	
	
	return(return_value);
}


si4	initialize_metadata(FILE_PROCESSING_STRUCT *fps)
{
	METADATA_SECTION_1		*md1;
	TIME_SERIES_METADATA_SECTION_2	*tmd2;
        VIDEO_METADATA_SECTION_2	*vmd2;
	METADATA_SECTION_3		*md3;
	
        
	// section 1 fields
	md1 = fps->metadata.section_1;
	md1->section_2_encryption = METADATA_SECTION_2_ENCRYPTION_DEFAULT;
	md1->section_3_encryption = METADATA_SECTION_3_ENCRYPTION_DEFAULT;
	
	// section 2 fields
	switch (fps->file_type_code) {
		case TIME_SERIES_METADATA_FILE_TYPE_CODE:
			tmd2 = fps->metadata.time_series_section_2;
			tmd2->recording_duration = METADATA_RECORDING_DURATION_NO_ENTRY;
			tmd2->acquisition_channel_number = TIME_SERIES_METADATA_ACQUISITION_CHANNEL_NUMBER_NO_ENTRY;
			tmd2->number_of_samples = TIME_SERIES_METADATA_NUMBER_OF_SAMPLES_NO_ENTRY;
			tmd2->sampling_frequency = TIME_SERIES_METADATA_SAMPLING_FREQUENCY_NO_ENTRY;
			tmd2->low_frequency_filter_setting = TIME_SERIES_METADATA_LOW_FREQUENCY_FILTER_SETTING_NO_ENTRY;
			tmd2->high_frequency_filter_setting = TIME_SERIES_METADATA_HIGH_FREQUENCY_FILTER_SETTING_NO_ENTRY;
			tmd2->notch_filter_frequency_setting = TIME_SERIES_METADATA_NOTCH_FILTER_FREQUENCY_SETTING_NO_ENTRY;
			tmd2->AC_line_frequency = TIME_SERIES_METADATA_AC_LINE_FREQUENCY_NO_ENTRY;
			tmd2->units_conversion_factor = TIME_SERIES_METADATA_UNITS_CONVERSION_FACTOR_NO_ENTRY;
			tmd2->maximum_native_sample_value = TIME_SERIES_METADATA_MAXIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY;
			tmd2->minimum_native_sample_value = TIME_SERIES_METADATA_MINIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY;
                        tmd2->start_sample = TIME_SERIES_METADATA_START_SAMPLE_NO_ENTRY;
			tmd2->maximum_difference_bytes = TIME_SERIES_METADATA_MAXIMUM_DIFFERENCE_BYTES_NO_ENTRY;
			tmd2->maximum_block_samples = TIME_SERIES_METADATA_MAXIMUM_BLOCK_SAMPLES_NO_ENTRY;
			tmd2->block_interval = TIME_SERIES_METADATA_BLOCK_INTERVAL_NO_ENTRY;
			tmd2->maximum_contiguous_samples = TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_SAMPLES_NO_ENTRY;
			tmd2->maximum_contiguous_block_bytes = TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCK_BYTES_NO_ENTRY;
			tmd2->number_of_discontinuities = TIME_SERIES_METADATA_NUMBER_OF_DISCONTINUITIES_NO_ENTRY;
			break;
		case VIDEO_METADATA_FILE_TYPE_CODE:
			vmd2 = fps->metadata.video_section_2;
			vmd2->horizontal_resolution = VIDEO_METADATA_HORIZONTAL_RESOLUTION_NO_ENTRY;
			vmd2->vertical_resolution = VIDEO_METADATA_VERTICAL_RESOLUTION_NO_ENTRY;
			vmd2->frame_rate = VIDEO_METADATA_FRAME_RATE_NO_ENTRY;
			vmd2->video_file_CRC = VIDEO_METADATA_VIDEO_FILE_CRC_NO_ENTRY;
			break;
		default:
			if (!(MEF_globals->behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
				fprintf(stderr, "Unrecognized METADATA SECTION 2 type \"%s\" [function \"%s\", line %d]\n", fps->full_file_name, __FUNCTION__, __LINE__);
				if (MEF_globals->behavior_on_fail & RETURN_ON_FAIL)
					(void) fprintf(stderr, "\t=> returning without initializing section 2\n\n");
				else if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL)
					(void) fprintf(stderr, "\t=> exiting program\n\n");
			}
			if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL)
				exit(1);
			break;
	}
        
	// section 3 fields
	md3 = fps->metadata.section_3;
	md3->recording_time_offset = MEF_globals->recording_time_offset;
	md3->DST_start_time = MEF_globals->DST_start_time;
	md3->DST_end_time = MEF_globals->DST_end_time;
        md3->GMT_offset = MEF_globals->GMT_offset;
	
        
	return(0);
}


si4	initialize_universal_header(FILE_PROCESSING_STRUCT *fps, si1 generate_level_UUID, si1 generate_file_UUID, si1 originating_file)
{
	UNIVERSAL_HEADER	*uh;
	
	
	uh = fps->universal_header;
	
	uh->header_CRC = CRC_NO_ENTRY;
	uh->body_CRC = CRC_NO_ENTRY;
	MEF_strncpy(uh->file_type_string, (si1 *) &fps->file_type_code, TYPE_BYTES);
	uh->mef_version_major = MEF_VERSION_MAJOR;
        uh->mef_version_minor = MEF_VERSION_MINOR;
        uh->byte_order_code = MEF_LITTLE_ENDIAN;
	uh->start_time = UUTC_NO_ENTRY;
	uh->end_time = UUTC_NO_ENTRY;
	uh->number_of_entries = UNIVERSAL_HEADER_NUMBER_OF_ENTRIES_NO_ENTRY;
	uh->maximum_entry_size = UNIVERSAL_HEADER_MAXIMUM_ENTRY_SIZE_NO_ENTRY;
	uh->segment_number = UNIVERSAL_HEADER_SEGMENT_NUMBER_NO_ENTRY;
	
	if (generate_level_UUID == MEF_TRUE)
		generate_UUID(uh->level_UUID);
	if (generate_file_UUID == MEF_TRUE)
		generate_UUID(uh->file_UUID);
	if (originating_file == MEF_TRUE)
		memcpy(uh->provenance_UUID, uh->file_UUID, UUID_BYTES);
	
        
	return(0);
}


si1	*local_date_time_string(si8 uutc_time, si1 *time_str)  // time_str buffer sould be of length 32
{
        si8		utc_time;
        si4		microseconds;
        si1		microseconds_str[7], year[5];
        struct tm	time_info;
        
        
	if (time_str == NULL)
		time_str = (si1 *) e_calloc((size_t) 32, sizeof(si1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	remove_recording_time_offset(&uutc_time);
        uutc_time += (si8) (MEF_globals->GMT_offset * 1e6);
        
        utc_time = (si8) uutc_time / 1000000;
        microseconds = (si4) (uutc_time % 1000000);
        numerical_fixed_width_string(microseconds_str, 6, microseconds);
        
        gmtime_r(&utc_time, &time_info);
        asctime_r(&time_info, time_str);
        
        time_str[24] = 0;
        strcpy(year, time_str + 20);
        time_str[19] = 0;
        sprintf(time_str, "%s.%s %s", time_str, microseconds_str, year);
        
        
        return(time_str);
}


si8	MEF_pad(ui1 *buffer, si8 content_len, ui4 alignment)
{
        si8	i, pad_bytes;
        
        
        pad_bytes = content_len % (si8) alignment;
        if (pad_bytes) {
                i = pad_bytes = (alignment - pad_bytes);
                buffer += content_len;
                while (i--)
                        *buffer++ = PAD_BYTE_VALUE;
        }
        
        
        return(content_len + pad_bytes);
}


void	MEF_snprintf(si1 *target, si4 target_field_bytes, si1 *format, ...)
{
	va_list	args;
	si1	*c;
	si4	bytes_to_zero;
        
	
	va_start(args, format);
        vsnprintf(target, target_field_bytes, format, args);
	va_end(args);
        
	c = target;
	while (*c++);
	
	bytes_to_zero = target_field_bytes - (c - target);
	if (bytes_to_zero > 0) {
		while (bytes_to_zero--) {
			*c++ = 0;
		}
	} else {
		target[target_field_bytes - 1] = 0;
	}
	
	
	return;
}


si4	MEF_sprintf(si1 *target, si1 *format, ...)
{
	va_list	args;
	si1	*c;
	
	
	va_start(args, format);
	vsprintf(target, format, args);
	va_end(args);
	
	c = target;
	while (*c++);
	

	return((si4) (c - target));
}


si4	MEF_strcat(si1 *target_string, si1 *source_string)
{
	si1 	*target_start;
	
	
	target_start = target_string;
	while ((*target_string++));
	--target_string;
	while ((*target_string++ = *source_string++));
	
	
	return((si4) (target_string - target_start));
}


si4	MEF_strcpy(si1 *target_string, si1 *source_string)
{
        si1 	*source_start;
        
        
        source_start = source_string;
        while ((*target_string++ = *source_string++));
        
        
        return((si4) (source_string - source_start));
}


void	MEF_strncat(si1 *target_string, si1 *source_string, si4 target_field_bytes)
{
        while (--target_field_bytes)
                if (*target_string++ == '\0')
                        break;
	
	--target_string;
	++target_field_bytes;
	       
	while (--target_field_bytes)
		if ((*target_string++ = *source_string++) == '\0')
			break;
	
	if (target_field_bytes)
		while (--target_field_bytes)
			*target_string++ = '\0';
        
        *target_string = '\0';
        
        
        return;
}


void	MEF_strncpy(si1 *target_string, si1 *source_string, si4 target_field_bytes)
{
        while (--target_field_bytes)
                if ((*target_string++ = *source_string++) == '\0')
                        break;
        
        if (target_field_bytes)
		while (--target_field_bytes)
			*target_string++ = '\0';
        
        *target_string = '\0';
        
        
        return;
}


si1	*numerical_fixed_width_string(si1 *string, si4 string_bytes, si4 number)
{
	si4	native_numerical_length, temp;
	si1	*c;
	
	// string bytes does not include terminal zero
	
	if (string == NULL)
		string = (si1 *) e_calloc((size_t) (string_bytes + 1), sizeof(si1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);

	native_numerical_length = 0;
	temp = number;
	while(temp) {
		temp /= 10;
		++native_numerical_length;
	}
	if (number <= 0)
		++native_numerical_length;
	
	c = string;
	temp = string_bytes - native_numerical_length;
	while (temp--)
		*c++ = '0';
        
	(void) sprintf(c, "%d", number);
	
        
	return(string);
}


si4	offset_record_index_times(FILE_PROCESSING_STRUCT *fps, si4 action)
{
	si1			apply, remove;
	ui4			mode;
	si8			i;
	RECORD_INDEX		*ri;
	
	
	mode = MEF_globals->recording_time_offset_mode;
	
	if (mode == RTO_IGNORE)
		return(0);
	
	apply = remove = MEF_UNKNOWN;
	if (action == RTO_INPUT_ACTION) {
		if (mode & (RTO_APPLY | RTO_APPLY_ON_INPUT))
			apply = MEF_TRUE;
		else
			apply = MEF_FALSE;
		if (mode & (RTO_REMOVE | RTO_REMOVE_ON_INPUT))
			remove = MEF_TRUE;
		else
			remove = MEF_FALSE;
	}
	else if (action == RTO_OUTPUT_ACTION) {
		if (mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT))
			apply = MEF_TRUE;
		else
			apply = MEF_FALSE;
		if (mode & (RTO_REMOVE | RTO_REMOVE_ON_OUTPUT))
			remove = MEF_TRUE;
		else
			remove = MEF_FALSE;
	}
	
	if (apply == MEF_TRUE && remove == MEF_TRUE) {
		fprintf(stderr, "%s(), line %d: both apply and remove specified for recording time offset mode => returning without doing either\n", __FUNCTION__, __LINE__);
		return(-1);
	}
	
	ri = fps->record_indices;
	if (apply == MEF_TRUE) {
		for (i = fps->universal_header->number_of_entries; i--;)
			apply_recording_time_offset(&(ri++)->time);
	} else if (remove == MEF_TRUE) {
		for (i = fps->universal_header->number_of_entries; i--;)
			remove_recording_time_offset(&(ri++)->time);
	}
	
	
	return(0);
}


si4	offset_time_series_index_times(FILE_PROCESSING_STRUCT *fps, si4 action)
{
	si1			apply, remove;
	ui4			mode;
	si8			i;
	TIME_SERIES_INDEX	*ti;
	
	
	mode = MEF_globals->recording_time_offset_mode;
	
	if (mode == RTO_IGNORE)
		return(0);
	
	apply = remove = MEF_UNKNOWN;
	if (action == RTO_INPUT_ACTION) {
		if (mode & (RTO_APPLY | RTO_APPLY_ON_INPUT))
			apply = MEF_TRUE;
		else
			apply = MEF_FALSE;
		if (mode & (RTO_REMOVE | RTO_REMOVE_ON_INPUT))
			remove = MEF_TRUE;
		else
			remove = MEF_FALSE;
	}
	else if (action == RTO_OUTPUT_ACTION) {
		if (mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT))
			apply = MEF_TRUE;
		else
			apply = MEF_FALSE;
		if (mode & (RTO_REMOVE | RTO_REMOVE_ON_OUTPUT))
			remove = MEF_TRUE;
		else
			remove = MEF_FALSE;
	}
	
	if (apply == MEF_TRUE && remove == MEF_TRUE) {
		fprintf(stderr, "%s(), line %d: both apply and remove specified for recording time offset mode => returning without doing either\n", __FUNCTION__, __LINE__);
		return(-1);
	}
	
	ti = fps->time_series_indices;
	if (apply == MEF_TRUE) {
		for (i = fps->universal_header->number_of_entries; i--;)
			apply_recording_time_offset(&(ti++)->start_time);
	} else if (remove == MEF_TRUE) {
		for (i = fps->universal_header->number_of_entries; i--;)
			remove_recording_time_offset(&(ti++)->start_time);
	}
	
	
	return(0);
}


si4	offset_universal_header_times(FILE_PROCESSING_STRUCT *fps, si4 action)
{
	si1			apply, remove;
	ui4			mode;
	UNIVERSAL_HEADER	*uh;
	
	
	mode = MEF_globals->recording_time_offset_mode;
	
	if (mode == RTO_IGNORE)
		return(0);
	
	apply = remove = MEF_UNKNOWN;
	if (action == RTO_INPUT_ACTION) {
		if (mode & (RTO_APPLY | RTO_APPLY_ON_INPUT))
			apply = MEF_TRUE;
		else
			apply = MEF_FALSE;
		if (mode & (RTO_REMOVE | RTO_REMOVE_ON_INPUT))
			remove = MEF_TRUE;
		else
			remove = MEF_FALSE;
	}
	else if (action == RTO_OUTPUT_ACTION) {
		if (mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT))
			apply = MEF_TRUE;
		else
			apply = MEF_FALSE;
		if (mode & (RTO_REMOVE | RTO_REMOVE_ON_OUTPUT))
			remove = MEF_TRUE;
		else
			remove = MEF_FALSE;
	}
	
	if (apply == MEF_TRUE && remove == MEF_TRUE) {
		fprintf(stderr, "%s(), line %d: both apply and remove specified for recording time offset mode => returning without doing either\n", __FUNCTION__, __LINE__);
		return(-1);
	}
	
	uh = fps->universal_header;
	if (apply == MEF_TRUE) {
		apply_recording_time_offset(&uh->start_time);
		apply_recording_time_offset(&uh->end_time);
	} else if (remove == MEF_TRUE) {
		remove_recording_time_offset(&uh->start_time);
		remove_recording_time_offset(&uh->end_time);
	}
	
	
	return(0);
}


si4	offset_video_index_times(FILE_PROCESSING_STRUCT *fps, si4 action)
{
	si1		apply, remove;
	ui4		mode;
	si8		i;
	VIDEO_INDEX	*vi;
	
	
	mode = MEF_globals->recording_time_offset_mode;
	
	if (mode == RTO_IGNORE)
		return(0);
	
	apply = remove = MEF_UNKNOWN;
	if (action == RTO_INPUT_ACTION) {
		if (mode & (RTO_APPLY | RTO_APPLY_ON_INPUT))
			apply = MEF_TRUE;
		else
			apply = MEF_FALSE;
		if (mode & (RTO_REMOVE | RTO_REMOVE_ON_INPUT))
			remove = MEF_TRUE;
		else
			remove = MEF_FALSE;
	}
	else if (action == RTO_OUTPUT_ACTION) {
		if (mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT))
			apply = MEF_TRUE;
		else
			apply = MEF_FALSE;
		if (mode & (RTO_REMOVE | RTO_REMOVE_ON_OUTPUT))
			remove = MEF_TRUE;
		else
			remove = MEF_FALSE;
	}
	
	if (apply == MEF_TRUE && remove == MEF_TRUE) {
		fprintf(stderr, "%s(), line %d: both apply and remove specified for recording time offset mode => returning without doing either\n", __FUNCTION__, __LINE__);
		return(-1);
	}
	
	vi = fps->video_indices;
	if (apply == MEF_TRUE) {
		for (i = fps->universal_header->number_of_entries; i--;) {
			apply_recording_time_offset(&vi->start_time);
			apply_recording_time_offset(&(vi++)->end_time);
		}
	} else if (remove == MEF_TRUE) {
		for (i = fps->universal_header->number_of_entries; i--;) {
			remove_recording_time_offset(&vi->start_time);
			remove_recording_time_offset(&(vi++)->end_time);
		}
	}
	
	
	return(0);
}


PASSWORD_DATA	*process_password_data(si1 *unspecified_password, si1 *level_1_password, si1 *level_2_password, UNIVERSAL_HEADER *universal_header)
{
	PASSWORD_DATA		*pwd;
	ui1			sha[SHA256_OUTPUT_SIZE];
	si1			password_bytes[PASSWORD_BYTES], l2_password_bytes[PASSWORD_BYTES];
	si1			putative_level_1_password_bytes[PASSWORD_BYTES];
	si4			i;
	
        
        // allocate
        pwd = (PASSWORD_DATA *) e_calloc((size_t) 1, sizeof(PASSWORD_DATA), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        pwd->access_level = LEVEL_0_ACCESS; // default access level
        
        // user passed single password for reading: validate against validation fields and generate encryption keys
	if (unspecified_password != NULL) {
		if (check_password(unspecified_password, __FUNCTION__, __LINE__) == 0) {
			
			// get terminal bytes
			extract_terminal_password_bytes(unspecified_password, password_bytes);
			
			// check for level 1 access
			sha256((ui1 *) password_bytes, PASSWORD_BYTES, sha);  // generate SHA-256 hash of password bytes
			
			for (i = 0; i < PASSWORD_VALIDATION_FIELD_BYTES; ++i)  // compare with stored level 1 hash
				if (sha[i] != universal_header->level_1_password_validation_field[i])
					break;
			if (i == PASSWORD_BYTES) {  // Level 1 password valid - cannot be level 2 password
				pwd->access_level = LEVEL_1_ACCESS;
				AES_key_expansion(pwd->level_2_encryption_key, password_bytes);  // generate key
				if (MEF_globals->verbose == MEF_TRUE)
					printf("Unspecified password is valid for Level 1 access\n");
				return(pwd);
			}
			
			// invalid level 1 => check if level 2 password
			for (i = 0; i < PASSWORD_BYTES; ++i)  // xor with level 2 password validation field
				putative_level_1_password_bytes[i] = sha[i] ^ universal_header->level_2_password_validation_field[i];
			
			sha256((ui1 *) putative_level_1_password_bytes, PASSWORD_BYTES, sha); // generate SHA-256 hash of putative level 1 password
			
			for (i = 0; i < PASSWORD_VALIDATION_FIELD_BYTES; ++i)  // compare with stored level 1 hash
				if (sha[i] != universal_header->level_1_password_validation_field[i])
					break;
			if (i == PASSWORD_VALIDATION_FIELD_BYTES) {  // Level 2 password valid
				pwd->access_level = LEVEL_2_ACCESS;
				AES_key_expansion(pwd->level_1_encryption_key, putative_level_1_password_bytes);  // generate key
				AES_key_expansion(pwd->level_2_encryption_key, password_bytes);  // generate key
				if (MEF_globals->verbose == MEF_TRUE)
					printf("Unspecified password is valid for Level 1 and Level 2 access\n");
			} else {
				fprintf(stderr, "%s(), line %d: unspecified password is not valid for Level 1 or Level 2 access\n", __FUNCTION__, __LINE__);
			}
		} else {
			fprintf(stderr, "%s(), line %d: unspecified password is not of valid form\n", __FUNCTION__, __LINE__);
			if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL) {
				(void) fprintf(stderr, "\t=> exiting program\n\n");
				exit(1);
			}
		}
		return(pwd);
	}
	
    // user passed level 1 password for writing: generate validation field and encryption key
    if (level_1_password != NULL){
		if (check_password(level_1_password, __FUNCTION__, __LINE__) == 0) {
	                
	                // passed a level 1 password - at least level 1 access
	                pwd->access_level = LEVEL_1_ACCESS;
	                
	                // get terminal bytes
	                extract_terminal_password_bytes(level_1_password, password_bytes);
			
	                // generate Level 1 password validation field
	                sha256((ui1 *) password_bytes, PASSWORD_BYTES, sha);
	                memcpy(universal_header->level_1_password_validation_field, sha, PASSWORD_VALIDATION_FIELD_BYTES);
	                if (MEF_globals->verbose == MEF_TRUE)
	                        printf("Level 1 password validation field generated\n");
	                
	                // generate encryption key
	                AES_key_expansion(pwd->level_1_encryption_key, password_bytes);
			if (MEF_globals->verbose == MEF_TRUE)
                printf("Level 1 encryption key generated\n");
                
                // user also passed level 2 password for writing: generate validation field and encryption key
                // Level 2 encryption requires a level 1 password, even if level 1 encryption is not used
            	if (level_2_password != NULL){
	                if (check_password(level_2_password, __FUNCTION__, __LINE__) == 0) {
	                        
	                        // passed a level 2 password - level 2 access
	                        pwd->access_level = LEVEL_2_ACCESS;
	                        
	                        // get terminal bytes
	                        extract_terminal_password_bytes(level_2_password, l2_password_bytes);
				
	                        // generate Level 2 password validation field
	                        sha256((ui1 *) l2_password_bytes, PASSWORD_BYTES, sha);
	                        memcpy(universal_header->level_2_password_validation_field, sha, PASSWORD_VALIDATION_FIELD_BYTES);
	                        for (i = 0; i < PASSWORD_VALIDATION_FIELD_BYTES; ++i) // exclusive or with level 1 password bytes
	                                universal_header->level_2_password_validation_field[i] ^= password_bytes[i];
	                        if (MEF_globals->verbose == MEF_TRUE)
	                                printf("Level 2 password validation field generated\n");
	                        
	                        // generate encryption key
	                        AES_key_expansion(pwd->level_2_encryption_key, l2_password_bytes);
	                        if (MEF_globals->verbose == MEF_TRUE)
	                                printf("Level 2 encryption key generated\n");
				
	                        
	                } else {
						fprintf(stderr, "%s(), line %d: Level 2 password is not of valid form\n", __FUNCTION__, __LINE__);
						if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL) {
							(void) fprintf(stderr, "\t=> exiting program\n\n");
							exit(1);
						}
					}
				}
		} else {
			fprintf(stderr, "%s(), line %d: Level 1 password is not of valid form\n", __FUNCTION__, __LINE__);
			if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL) {
				(void) fprintf(stderr, "\t=> exiting program\n\n");
				exit(1);
			}
		}
    }  
        
	return(pwd);
}


void    proportion_filt(sf8 *x, sf8 *px, si8 len, sf8 prop, si4 span)
{
	NODE    *nodes, *prop_node, *first_node, *last_node, *empty_node;
	NODE	*curr_node, *next_node, *prev_node, *head, *tail;
	si4     first_idx, last_idx, empty_idx, half_span, prop_span, prop_idx, span_plus_one;
	si8     i, j;
	sf8     prop_val, new_val, first_val, prop_shift, temp_shift;
	
	
	/* setup */
	if ((span % 2) == 0) // require odd-numbered window
		++span;
	if (px == NULL) // caller responsible for freeing px
		px = (sf8 *) calloc((size_t) len, sizeof(sf8));
	nodes = (NODE *) calloc((size_t) (span + 3), sizeof(NODE));
	first_idx = 0;
	first_node = nodes + first_idx;
	last_idx = span - 1;
	last_node = nodes + last_idx;
	empty_idx = span;
	empty_node = nodes + empty_idx;
	head = nodes + (span + 1);
	tail = nodes + (span + 2);
	half_span = span / 2;
	prop_span = (si4) (((sf8) (span - 1) * prop) + 0.5);
	span_plus_one = span + 1;
	
	/* build linked list */
	head->val = -DBL_MAX;
	tail->val = DBL_MAX;
	for (i = 0; i < span; ++i) {
		nodes[i].val = x[i];
		nodes[i].idx = i;
	}
	qsort(nodes, span, sizeof(NODE), sort_by_val);
	for (i = 1; i < span; ++i)
		nodes[i - 1].next = nodes + nodes[i].idx;
	nodes[span - 1].next = tail;
	head->next = nodes + nodes[0].idx;
	for (i = 1; i < span; ++i)
		nodes[i].prev = nodes + nodes[i - 1].idx;
	nodes[0].prev = head;
	tail->prev = nodes + nodes[span - 1].idx;
	prop_idx = nodes[prop_span].idx;
	qsort(nodes, span, sizeof(NODE), sort_by_idx);
	prop_node = nodes + prop_idx;
	
	/* fill in initial segment */
	prop_val = prop_node->val;
	for (i = 0; i <= half_span; ++i)
		px[i] = prop_val;
	
	/* slide window */
	for (i = span, j = half_span + 1; i < len; ++i, ++j) {
		
		/* insert new value into list using empty node: Note routine doesn't handle NaNs */
		
		if (isnan(empty_node->val = new_val = x[i])) {
			if (!(MEF_globals->behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
				fprintf(stderr, "Proportion_filt() does not currently handle NaN values [function \"%s\", line %d]\n", __FUNCTION__, __LINE__);
				if (MEF_globals->behavior_on_fail & RETURN_ON_FAIL)
					(void) fprintf(stderr, "\t=> returning without filtering\n\n");
				else if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL)
					(void) fprintf(stderr, "\t=> exiting program\n\n");
			}
			if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL)
				exit(1);
			return;
		}
		if (new_val > prop_val)
			prop_shift = 0.5;
		else if (new_val < prop_val)
			prop_shift = -0.5;
		else
			prop_shift = 0.0;
		curr_node = last_node;
		if (new_val > last_node->val) {  // search forward
			while (1) {
				next_node = curr_node->next;
				if (new_val <= next_node->val)
					break;
				curr_node = next_node;
			}
			if (prop_shift == 0.0) {
				prop_shift = val_equals_prop(curr_node, prop_node);
				if (prop_shift == 0.0)
					prop_shift = 0.5;
			}
			// insert after curr_node
			empty_node->prev = curr_node;
			empty_node->next = next_node;
			next_node->prev = curr_node->next = empty_node;
		} else {  // search backward
			while (1) {
				prev_node = curr_node->prev;
				if (new_val >= prev_node->val)
					break;
				curr_node = prev_node;
			}
			if (prop_shift == 0.0) {
				prop_shift = val_equals_prop(curr_node, prop_node);
				if (prop_shift == 0.0)
					prop_shift = -0.5;
			}
			// insert before curr_node
			empty_node->next = curr_node;
			empty_node->prev = prev_node;
			prev_node->next = curr_node->prev = empty_node;
		}
		
		/* update proportion value */
		first_val = first_node->val;
		if (first_val > prop_val)
			prop_shift -= 0.5;
		else if (first_val < prop_val)
			prop_shift += 0.5;
		else {
			temp_shift = val_equals_prop(first_node, prop_node);
			if (temp_shift == 0.0) {
				if (prop_shift > 0.0)
					prop_shift = 1.0;
				else
					prop_shift = -1.0;
			} else {
				prop_shift -= temp_shift;
			}
		}
		if (prop_shift == 1.0) {
			prop_node = prop_node->next;
			prop_val = prop_node->val;
		} else if (prop_shift == -1.0) {
			prop_node = prop_node->prev;
			prop_val = prop_node->val;
		}
		px[j] = prop_val;
		
		/* excise first (oldest) node */
		prev_node = first_node->prev;
		next_node = first_node->next;
		prev_node->next = next_node;
		next_node->prev = prev_node;
		
		/* update window */
		++first_idx; first_idx %= span_plus_one;
		first_node = nodes + first_idx;
		++last_idx; last_idx %= span_plus_one;
		last_node = nodes + last_idx;
		++empty_idx; empty_idx %= span_plus_one;
		empty_node = nodes + empty_idx;
	}
	
	/* fill in terminal segment */
	for (; j < len; ++j)
		px[j] = prop_val;
	
	/* clean up */
	free(nodes);
	
	return;
}


inline ui1	random_byte(ui4 *m_w, ui4 *m_z)
{
	ui1	rb;
        
        
        *m_z = 0x00009069 * (*m_z & 0x0000FFFF) + (*m_z >> 0x10);
        *m_w = 0x00004650 * (*m_w & 0x0000FFFF) + (*m_w >> 0x10);
	rb = (ui1) (((*m_z << 0x10) + *m_w) % 0x00000100);
        
        
        return(rb);
}


CHANNEL	*read_MEF_channel(CHANNEL *channel, si1 *chan_path, si4 channel_type, si1 *password, PASSWORD_DATA *password_data, si1 read_time_series_data, si1 read_record_data)
{
	si4				i, n_segments;
	si1				full_file_name[MEF_FULL_FILE_NAME_BYTES], **segment_names;
	METADATA_SECTION_1		*smd1, *cmd1;
        TIME_SERIES_METADATA_SECTION_2	*ctmd, *stmd;
        VIDEO_METADATA_SECTION_2	*cvmd, *svmd;
	METADATA_SECTION_3		*smd3, *cmd3;
        SEGMENT				*seg;
	FILE_PROCESSING_STRUCT		*temp_fps;
	
	
	// allocate channel if not passed
	if (channel == NULL)
		channel = (CHANNEL *) e_calloc((size_t) 1, sizeof(CHANNEL), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        
	// get channel path & name
	extract_path_parts(chan_path, channel->path, channel->name, channel->extension);
        
        // get session name - if channel records present, we could get it from one of the universal headers
        extract_path_parts(channel->path, NULL, channel->session_name, NULL);
        
        // get_channel_type
        if (channel_type == UNKNOWN_CHANNEL_TYPE)
                channel_type = channel_type_from_path(chan_path);
        channel->channel_type = channel_type;
	
	// initialize
	channel->maximum_number_of_records = 0;
	channel->maximum_record_bytes = 0;
	bzero(channel->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);

	#ifdef _WIN32
		channel->earliest_start_time = LLONG_MAX;
		channel->latest_end_time = LLONG_MIN;
	#else
		channel->earliest_start_time = LONG_MAX;
		channel->latest_end_time = LONG_MIN;
	#endif
	
	// loop over segments
	segment_names = generate_file_list(NULL, &n_segments, chan_path, SEGMENT_DIRECTORY_TYPE_STRING);
	channel->segments = (SEGMENT *) e_calloc((size_t) n_segments, sizeof(SEGMENT), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        channel->number_of_segments = n_segments;
	for (i = 0; i < n_segments; ++i) {
		(void) read_MEF_segment(channel->segments + i, segment_names[i], channel_type, password, password_data, read_time_series_data, read_record_data);
		if (password_data == NULL)
			password_data = channel->segments[i].metadata_fps->password_data;
		free(segment_names[i]);
	}
	free(segment_names);
        
        // fill in channel metadata
        if (channel->metadata.section_1 == NULL)
                channel->metadata.section_1 = (METADATA_SECTION_1 *) e_calloc((size_t) 1, sizeof(METADATA_SECTION_1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        if (channel->metadata.section_3 == NULL)
                channel->metadata.section_3 = (METADATA_SECTION_3 *) e_calloc((size_t) 1, sizeof(METADATA_SECTION_3), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	if (channel->channel_type == TIME_SERIES_CHANNEL_TYPE) {
		if (channel->metadata.time_series_section_2 == NULL)
			channel->metadata.time_series_section_2 = (TIME_SERIES_METADATA_SECTION_2 *) e_calloc((size_t) 1, sizeof(TIME_SERIES_METADATA_SECTION_2), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	} else if (channel->channel_type == VIDEO_CHANNEL_TYPE) {
		if (channel->metadata.video_section_2 == NULL)
			channel->metadata.video_section_2 = (VIDEO_METADATA_SECTION_2 *) e_calloc((size_t) 1, sizeof(VIDEO_METADATA_SECTION_2), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	}

        // times series channel
        if (channel->channel_type == TIME_SERIES_CHANNEL_TYPE) {
		
                for (i = 0; i < n_segments; ++i) {
                        seg = channel->segments + i;
			smd1 = seg->metadata_fps->metadata.section_1;
			cmd1 = channel->metadata.section_1;
			stmd = seg->metadata_fps->metadata.time_series_section_2;
                        ctmd = channel->metadata.time_series_section_2;
			smd3 = seg->metadata_fps->metadata.section_3;
			cmd3 = channel->metadata.section_3;
			if (i == 0) {
                                memcpy(cmd1, smd1, METADATA_SECTION_1_BYTES);
                                memcpy(ctmd, stmd, METADATA_SECTION_2_BYTES);
                                memcpy(cmd3, smd3, METADATA_SECTION_3_BYTES);
                                if (ABS(channel->earliest_start_time) > ABS(seg->metadata_fps->universal_header->start_time))
					channel->earliest_start_time = seg->metadata_fps->universal_header->start_time;
                                if (ABS(channel->latest_end_time) < ABS(seg->metadata_fps->universal_header->end_time))
					channel->latest_end_time = seg->metadata_fps->universal_header->end_time;
                                if (seg->record_data_fps != NULL) {
                                	if (channel->maximum_number_of_records < seg->record_data_fps->universal_header->number_of_entries)
						channel->maximum_number_of_records = seg->record_data_fps->universal_header->number_of_entries;
                                        if (channel->maximum_record_bytes < seg->record_data_fps->universal_header->maximum_entry_size)
						channel->maximum_record_bytes = seg->record_data_fps->universal_header->maximum_entry_size;
                                }
				if (strlen(channel->anonymized_name) == 0)
					MEF_strncpy(channel->anonymized_name, seg->metadata_fps->universal_header->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
				else if (strcmp(channel->anonymized_name, seg->metadata_fps->universal_header->anonymized_name))
					bzero(channel->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
                                continue;
                        }
			// universal header
			if (ABS(seg->metadata_fps->universal_header->start_time) < ABS(channel->earliest_start_time))
				channel->earliest_start_time = seg->metadata_fps->universal_header->start_time;
			if (ABS(seg->metadata_fps->universal_header->end_time) > ABS(channel->latest_end_time))
				channel->latest_end_time = seg->metadata_fps->universal_header->end_time;
			if (strcmp(channel->anonymized_name, seg->metadata_fps->universal_header->anonymized_name))
				bzero(channel->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
			if (seg->record_data_fps != NULL) {
				if (seg->record_data_fps->universal_header->number_of_entries > channel->maximum_number_of_records)
					channel->maximum_number_of_records = seg->record_data_fps->universal_header->number_of_entries;
				if (seg->record_data_fps->universal_header->maximum_entry_size > channel->maximum_record_bytes)
					channel->maximum_record_bytes = seg->record_data_fps->universal_header->maximum_entry_size;
			}
			// section 1
			if (smd1->section_2_encryption != cmd1->section_2_encryption)
				cmd1->section_2_encryption = ENCRYPTION_LEVEL_NO_ENTRY;
			if (smd1->section_3_encryption != cmd1->section_3_encryption)
				cmd1->section_3_encryption = ENCRYPTION_LEVEL_NO_ENTRY;
			if (memcmp(smd1->protected_region, cmd1->protected_region, METADATA_SECTION_1_PROTECTED_REGION_BYTES))
				bzero(cmd1->protected_region, METADATA_SECTION_1_PROTECTED_REGION_BYTES);
			if (memcmp(smd1->discretionary_region, cmd1->discretionary_region, METADATA_SECTION_1_DISCRETIONARY_REGION_BYTES))
				bzero(cmd1->discretionary_region, METADATA_SECTION_1_DISCRETIONARY_REGION_BYTES);
			// section 2
			if (strcmp(stmd->channel_description, ctmd->channel_description))
				bzero(ctmd->channel_description, METADATA_CHANNEL_DESCRIPTION_BYTES);
			if (strcmp(stmd->session_description, ctmd->session_description))
				bzero(ctmd->session_description, METADATA_SESSION_DESCRIPTION_BYTES);
			if (channel->latest_end_time == UUTC_NO_ENTRY || channel->earliest_start_time == UUTC_NO_ENTRY)
				ctmd->recording_duration = UUTC_NO_ENTRY;
			else
				ctmd->recording_duration = ABS(channel->latest_end_time) - ABS(channel->earliest_start_time) + 1;
			if (strcmp(stmd->reference_description, ctmd->reference_description))
				bzero(ctmd->reference_description, METADATA_CHANNEL_DESCRIPTION_BYTES);
			if (stmd->acquisition_channel_number != ctmd->acquisition_channel_number)
				ctmd->acquisition_channel_number = TIME_SERIES_METADATA_ACQUISITION_CHANNEL_NUMBER_NO_ENTRY;
                        if (stmd->sampling_frequency != ctmd->sampling_frequency)
                                ctmd->sampling_frequency = TIME_SERIES_METADATA_SAMPLING_FREQUENCY_NO_ENTRY;
                        if (stmd->low_frequency_filter_setting != ctmd->low_frequency_filter_setting)
                                ctmd->low_frequency_filter_setting = TIME_SERIES_METADATA_LOW_FREQUENCY_FILTER_SETTING_NO_ENTRY;
                        if (stmd->high_frequency_filter_setting != ctmd->high_frequency_filter_setting)
                                ctmd->high_frequency_filter_setting = TIME_SERIES_METADATA_HIGH_FREQUENCY_FILTER_SETTING_NO_ENTRY;
                        if (stmd->notch_filter_frequency_setting != ctmd->notch_filter_frequency_setting)
                                ctmd->notch_filter_frequency_setting = TIME_SERIES_METADATA_NOTCH_FILTER_FREQUENCY_SETTING_NO_ENTRY;
			if (stmd->AC_line_frequency != ctmd->AC_line_frequency)
				ctmd->AC_line_frequency = TIME_SERIES_METADATA_AC_LINE_FREQUENCY_NO_ENTRY;
			if (stmd->units_conversion_factor != ctmd->units_conversion_factor)
                                ctmd->units_conversion_factor = TIME_SERIES_METADATA_UNITS_CONVERSION_FACTOR_NO_ENTRY;
                        if (strcmp(stmd->units_description, ctmd->units_description))
                                bzero(ctmd->units_description, TIME_SERIES_METADATA_UNITS_DESCRIPTION_BYTES);
                        if (stmd->maximum_native_sample_value > ctmd->maximum_native_sample_value)
                                ctmd->maximum_native_sample_value = stmd->maximum_native_sample_value;
                        if (stmd->minimum_native_sample_value < ctmd->minimum_native_sample_value)
                                ctmd->minimum_native_sample_value = stmd->minimum_native_sample_value;
                        if (stmd->start_sample < ctmd->start_sample)
                        	ctmd->start_sample = stmd->start_sample;
			if (i > 0)
				ctmd->number_of_samples += stmd->number_of_samples;
			if (i > 0)
				ctmd->number_of_blocks += stmd->number_of_blocks;
			if (stmd->maximum_block_bytes > ctmd->maximum_block_bytes)
				ctmd->maximum_block_bytes = stmd->maximum_block_bytes;
			if (stmd->maximum_block_samples > ctmd->maximum_block_samples)
				ctmd->maximum_block_samples = stmd->maximum_block_samples;
			if (stmd->maximum_difference_bytes > ctmd->maximum_difference_bytes)
                                ctmd->maximum_difference_bytes = stmd->maximum_difference_bytes;
                        if (stmd->block_interval != ctmd->block_interval)
                                ctmd->block_interval = TIME_SERIES_METADATA_BLOCK_INTERVAL_NO_ENTRY;
			if (stmd->number_of_discontinuities > ctmd->number_of_discontinuities)
				ctmd->number_of_discontinuities = stmd->number_of_discontinuities;
			if (stmd->maximum_contiguous_blocks > ctmd->maximum_contiguous_blocks)
				ctmd->maximum_contiguous_blocks = stmd->maximum_contiguous_blocks;
			if (stmd->maximum_contiguous_block_bytes > ctmd->maximum_contiguous_block_bytes)
				ctmd->maximum_contiguous_block_bytes = stmd->maximum_contiguous_block_bytes;
			if (stmd->maximum_contiguous_samples > ctmd->maximum_contiguous_samples)
                                ctmd->maximum_contiguous_samples = stmd->maximum_contiguous_samples;
                        if (memcmp(stmd->protected_region, ctmd->protected_region, TIME_SERIES_METADATA_SECTION_2_PROTECTED_REGION_BYTES))
                                bzero(ctmd->protected_region, TIME_SERIES_METADATA_SECTION_2_PROTECTED_REGION_BYTES);
                        if (memcmp(stmd->discretionary_region, ctmd->discretionary_region, TIME_SERIES_METADATA_SECTION_2_DISCRETIONARY_REGION_BYTES))
                                bzero(ctmd->discretionary_region, TIME_SERIES_METADATA_SECTION_2_DISCRETIONARY_REGION_BYTES);
			// section 3
			if (smd3->recording_time_offset != cmd3->recording_time_offset)
				cmd3->recording_time_offset = UUTC_NO_ENTRY;
			if (smd3->DST_start_time != cmd3->DST_start_time)
				cmd3->DST_start_time = UUTC_NO_ENTRY;
			if (smd3->DST_end_time != cmd3->DST_end_time)
				cmd3->DST_end_time = UUTC_NO_ENTRY;
			if (smd3->GMT_offset != cmd3->GMT_offset)
				cmd3->GMT_offset = GMT_OFFSET_NO_ENTRY;
			if (strcmp(smd3->subject_name_1, cmd3->subject_name_1))
				bzero(cmd3->subject_name_1, METADATA_SUBJECT_NAME_BYTES);
			if (strcmp(smd3->subject_name_2, cmd3->subject_name_2))
				bzero(cmd3->subject_name_2, METADATA_SUBJECT_NAME_BYTES);
			if (strcmp(smd3->subject_ID, cmd3->subject_ID))
				bzero(cmd3->subject_ID, METADATA_SUBJECT_ID_BYTES);
			if (strcmp(smd3->recording_location, cmd3->recording_location))
				bzero(cmd3->recording_location, METADATA_RECORDING_LOCATION_BYTES);
			if (memcmp(smd3->protected_region, cmd3->protected_region, METADATA_SECTION_3_PROTECTED_REGION_BYTES))
				bzero(cmd3->protected_region, METADATA_SECTION_3_PROTECTED_REGION_BYTES);
			if (memcmp(smd3->discretionary_region, cmd3->discretionary_region, METADATA_SECTION_3_DISCRETIONARY_REGION_BYTES))
				bzero(cmd3->discretionary_region, METADATA_SECTION_3_DISCRETIONARY_REGION_BYTES);
                }
        }
        
        // video channel
        if (channel->channel_type == VIDEO_CHANNEL_TYPE) {

                for (i = 0; i < n_segments; ++i) {
			seg = channel->segments + i;
			smd1 = seg->metadata_fps->metadata.section_1;
			cmd1 = channel->metadata.section_1;
			svmd = seg->metadata_fps->metadata.video_section_2;
			cvmd = channel->metadata.video_section_2;
			smd3 = seg->metadata_fps->metadata.section_3;
			cmd3 = channel->metadata.section_3;
			if (i == 0) {
				memcpy(cmd1, smd1, METADATA_SECTION_1_BYTES);
				memcpy(cvmd, svmd, METADATA_SECTION_2_BYTES);
				memcpy(cmd3, smd3, METADATA_SECTION_3_BYTES);
				if (ABS(channel->earliest_start_time) > ABS(seg->metadata_fps->universal_header->start_time))
					channel->earliest_start_time = seg->metadata_fps->universal_header->start_time;
				if (ABS(channel->latest_end_time) < ABS(seg->metadata_fps->universal_header->end_time))
					channel->latest_end_time = seg->metadata_fps->universal_header->end_time;
				if (seg->record_data_fps != NULL) {
					if (channel->maximum_number_of_records < seg->record_data_fps->universal_header->number_of_entries)
						channel->maximum_number_of_records = seg->record_data_fps->universal_header->number_of_entries;
					if (channel->maximum_record_bytes < seg->record_data_fps->universal_header->maximum_entry_size)
						channel->maximum_record_bytes = seg->record_data_fps->universal_header->maximum_entry_size;
				}
				if (strlen(channel->anonymized_name) == 0)
					MEF_strncpy(channel->anonymized_name, seg->metadata_fps->universal_header->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
				else if (strcmp(channel->anonymized_name, seg->metadata_fps->universal_header->anonymized_name))
					bzero(channel->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
				continue;
			}
			// universal header
			if (ABS(seg->metadata_fps->universal_header->start_time) < ABS(channel->earliest_start_time))
				channel->earliest_start_time = seg->metadata_fps->universal_header->start_time;
			if (ABS(seg->metadata_fps->universal_header->end_time) > ABS(channel->latest_end_time))
				channel->latest_end_time = seg->metadata_fps->universal_header->end_time;
			if (strcmp(channel->anonymized_name, seg->metadata_fps->universal_header->anonymized_name))
				bzero(channel->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
			if (seg->record_data_fps != NULL) {
				if (seg->record_data_fps->universal_header->number_of_entries > channel->maximum_number_of_records)
					channel->maximum_number_of_records = seg->record_data_fps->universal_header->number_of_entries;
				if (seg->record_data_fps->universal_header->maximum_entry_size > channel->maximum_record_bytes)
					channel->maximum_record_bytes = seg->record_data_fps->universal_header->maximum_entry_size;
			}
			// section 1
			if (smd1->section_2_encryption != cmd1->section_2_encryption)
				cmd1->section_2_encryption = ENCRYPTION_LEVEL_NO_ENTRY;
			if (smd1->section_3_encryption != cmd1->section_3_encryption)
				cmd1->section_3_encryption = ENCRYPTION_LEVEL_NO_ENTRY;
			if (memcmp(smd1->protected_region, cmd1->protected_region, METADATA_SECTION_1_PROTECTED_REGION_BYTES))
				bzero(cmd1->protected_region, METADATA_SECTION_1_PROTECTED_REGION_BYTES);
			if (memcmp(smd1->discretionary_region, cmd1->discretionary_region, METADATA_SECTION_1_DISCRETIONARY_REGION_BYTES))
				bzero(cmd1->discretionary_region, METADATA_SECTION_1_DISCRETIONARY_REGION_BYTES);
			// section 2
			if (strcmp(svmd->channel_description, cvmd->channel_description))
				bzero(cvmd->channel_description, METADATA_CHANNEL_DESCRIPTION_BYTES);
			if (strcmp(svmd->session_description, cvmd->session_description))
				bzero(cvmd->session_description, METADATA_SESSION_DESCRIPTION_BYTES);
			if (channel->latest_end_time == UUTC_NO_ENTRY || channel->earliest_start_time == UUTC_NO_ENTRY)
				cvmd->recording_duration = UUTC_NO_ENTRY;
			else
				cvmd->recording_duration = ABS(channel->latest_end_time) - ABS(channel->earliest_start_time) + 1;
                        if (svmd->horizontal_resolution != cvmd->horizontal_resolution)
                                cvmd->horizontal_resolution = VIDEO_METADATA_HORIZONTAL_RESOLUTION_NO_ENTRY;
                        if (svmd->vertical_resolution != cvmd->vertical_resolution)
                                cvmd->vertical_resolution = VIDEO_METADATA_VERTICAL_RESOLUTION_NO_ENTRY;
                        if (svmd->frame_rate != cvmd->frame_rate)
                                cvmd->frame_rate = VIDEO_METADATA_FRAME_RATE_NO_ENTRY;
			if (i > 0)
				cvmd->number_of_clips += svmd->number_of_clips;
			if (svmd->maximum_clip_bytes > cvmd->maximum_clip_bytes)
				cvmd->maximum_clip_bytes = svmd->maximum_clip_bytes;
			if (strcmp(svmd->video_format, cvmd->video_format))
				bzero(cvmd->video_format, VIDEO_METADATA_VIDEO_FORMAT_BYTES);
                        if (svmd->video_file_CRC != cvmd->video_file_CRC)
                                cvmd->video_file_CRC = CRC_NO_ENTRY;
                        if (memcmp(svmd->protected_region, cvmd->protected_region, VIDEO_METADATA_SECTION_2_PROTECTED_REGION_BYTES))
                                bzero(cvmd->protected_region, VIDEO_METADATA_SECTION_2_PROTECTED_REGION_BYTES);
                        if (memcmp(svmd->discretionary_region, cvmd->discretionary_region, VIDEO_METADATA_SECTION_2_DISCRETIONARY_REGION_BYTES))
                                bzero(cvmd->discretionary_region, VIDEO_METADATA_SECTION_2_DISCRETIONARY_REGION_BYTES);
			// section 3
			if (smd3->recording_time_offset != cmd3->recording_time_offset)
				cmd3->recording_time_offset = UUTC_NO_ENTRY;
			if (smd3->DST_start_time != cmd3->DST_start_time)
				cmd3->DST_start_time = UUTC_NO_ENTRY;
			if (smd3->DST_end_time != cmd3->DST_end_time)
				cmd3->DST_end_time = UUTC_NO_ENTRY;
			if (smd3->GMT_offset != cmd3->GMT_offset)
				cmd3->GMT_offset = GMT_OFFSET_NO_ENTRY;
			if (strcmp(smd3->subject_name_1, cmd3->subject_name_1))
				bzero(cmd3->subject_name_1, METADATA_SUBJECT_NAME_BYTES);
			if (strcmp(smd3->subject_name_2, cmd3->subject_name_2))
				bzero(cmd3->subject_name_2, METADATA_SUBJECT_NAME_BYTES);
			if (strcmp(smd3->subject_ID, cmd3->subject_ID))
				bzero(cmd3->subject_ID, METADATA_SUBJECT_ID_BYTES);
			if (strcmp(smd3->recording_location, cmd3->recording_location))
				bzero(cmd3->recording_location, METADATA_RECORDING_LOCATION_BYTES);
			if (memcmp(smd3->protected_region, smd3->protected_region, METADATA_SECTION_3_PROTECTED_REGION_BYTES))
				bzero(cmd3->protected_region, METADATA_SECTION_3_PROTECTED_REGION_BYTES);
			if (memcmp(smd3->discretionary_region, cmd3->discretionary_region, METADATA_SECTION_3_DISCRETIONARY_REGION_BYTES))
				bzero(cmd3->discretionary_region, METADATA_SECTION_3_DISCRETIONARY_REGION_BYTES);
                }
        }
	
        // read channel record indices if present
		MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s/%s.%s", channel->path, channel->name, channel->extension, channel->name, RECORD_INDICES_FILE_TYPE_STRING);
		channel->record_indices_fps = read_MEF_file(NULL, full_file_name, password, password_data, NULL, RETURN_ON_FAIL | SUPPRESS_ERROR_OUTPUT);
	    if (channel->record_indices_fps != NULL) {
			if (password_data == NULL)
				password_data = channel->record_indices_fps->password_data;
			// copy level UUID
			memcpy(channel->level_UUID, channel->record_indices_fps->universal_header->level_UUID, UUID_BYTES);
            // read channel record data
            MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s/%s.%s", channel->path, channel->name, channel->extension, channel->name, RECORD_DATA_FILE_TYPE_STRING);
            channel->record_data_fps = allocate_file_processing_struct(0, RECORD_DATA_FILE_TYPE_CODE, NULL, NULL, 0);
            if (read_record_data == MEF_FALSE) {
                    channel->record_data_fps->directives.io_bytes = UNIVERSAL_HEADER_BYTES;
            }
            (void) read_MEF_file(channel->record_data_fps, full_file_name, password, password_data, NULL, RETURN_ON_FAIL | SUPPRESS_ERROR_OUTPUT);
            if (channel->record_data_fps == NULL)
                    UTF8_fprintf(stderr, "%s() Warning: Channel record indices file, but no channel record data file (\"%s\") in channel directory\n\n", __FUNCTION__, full_file_name);
			
			if (channel->maximum_number_of_records < channel->record_data_fps->universal_header->number_of_entries)
				channel->maximum_number_of_records = channel->record_data_fps->universal_header->number_of_entries;
			if (channel->maximum_record_bytes < channel->record_data_fps->universal_header->maximum_entry_size)
				channel->maximum_record_bytes = channel->record_data_fps->universal_header->maximum_entry_size;
			if (ABS(channel->record_data_fps->universal_header->start_time) < ABS(channel->earliest_start_time))
				channel->earliest_start_time = channel->record_data_fps->universal_header->start_time;
			if (ABS(channel->record_data_fps->universal_header->end_time) > ABS(channel->latest_end_time))
				channel->latest_end_time = channel->record_data_fps->universal_header->end_time;
			MEF_strncpy(channel->anonymized_name, channel->record_data_fps->universal_header->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
		}

	if (MEF_globals->verbose == MEF_TRUE) {
		if (channel_type == TIME_SERIES_CHANNEL_TYPE) {
			printf("------------ Time Series Channel Metadata --------------\n");
			temp_fps = allocate_file_processing_struct(0, TIME_SERIES_METADATA_FILE_TYPE_CODE, NULL, NULL, 0);
		} else if (channel_type == VIDEO_CHANNEL_TYPE) {
			printf("--------------- Video Channel Metadata -----------------\n");
			temp_fps = allocate_file_processing_struct(0, VIDEO_METADATA_FILE_TYPE_CODE, NULL, NULL, 0);
		} else {
			return(channel);
		}
		temp_fps->metadata = channel->metadata;
		temp_fps->password_data = password_data;
		show_metadata(temp_fps);
		free_file_processing_struct(temp_fps);
	}


	
	
	return(channel);
}

#ifdef _WIN32
	void slash_to_backslash(si1* orig_string)
	{
	    si1     *c;

	    c = orig_string;

	    while(*c != '\0')
	    {
	        if (*c == '/'){
	            *c = '\\';
	        }
	        *c++;
	    }
	}
#endif

FILE_PROCESSING_STRUCT	*read_MEF_file(FILE_PROCESSING_STRUCT *fps, si1 *file_name, si1 *password, PASSWORD_DATA *password_data, FILE_PROCESSING_DIRECTIVES *directives, ui4 behavior_on_fail)
{
	si8	i_bytes;
	ui4 *file_type_string_int;
	si4	allocated_fps, CRC_result;
    void	*data_ptr;
	
    if (access(file_name, 0) == -1)
	{
	   // file doesn't exist
	   return (NULL);
	}
	
	#ifdef _WIN32
        	if (fps->full_file_name != NULL)
        		slash_to_backslash(fps->full_file_name);
	#endif
        
	if (behavior_on_fail == USE_GLOBAL_BEHAVIOR)
		behavior_on_fail = MEF_globals->behavior_on_fail;
	
	// allocate FILE_PROCESSING_STRUCT if required
	if (fps == NULL) {
		allocated_fps = MEF_TRUE;
		fps = allocate_file_processing_struct(0, 0, directives, NULL, 0);
	} else
		allocated_fps = MEF_FALSE;
        
	// copy file name if passed
	if (file_name != NULL)
		MEF_strncpy(fps->full_file_name, file_name, MEF_FULL_FILE_NAME_BYTES);
	
	// open file if not already open
	if (fps->fp == NULL) {
		if (!(fps->directives.open_mode & FPS_GENERIC_READ_OPEN_MODE))
			fps->directives.open_mode = FPS_R_OPEN_MODE;
		fps_open(fps, __FUNCTION__, __LINE__, behavior_on_fail);
		if (fps->fp == NULL) {
			if (allocated_fps == MEF_TRUE)
				free_file_processing_struct(fps);
			return(NULL);
		}
	} else
		e_fseek(fps->fp, 0, SEEK_SET, fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	// check file not empty
	if (fps->file_length == 0) {
		if (!(fps->directives.open_mode & FPS_GENERIC_READ_OPEN_MODE))
			fps->directives.open_mode = FPS_R_OPEN_MODE;
		fps_close(fps);
		if (allocated_fps == MEF_TRUE)
			free_file_processing_struct(fps);
		return(NULL);
	}
	
	// get read size
	if (fps->directives.io_bytes == FPS_FULL_FILE)
		i_bytes = fps->file_length;
	else
		i_bytes = fps->directives.io_bytes;
        
	// allocate raw data
	if (fps->raw_data == NULL) {
		fps->raw_data = (ui1 *) e_calloc((size_t) i_bytes, sizeof(ui1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		fps->raw_data_bytes = i_bytes;
	}
        
	// read in raw data
	fps_read(fps, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        
	// close
	if (fps->directives.close_file == MEF_TRUE)
		fps_close(fps);
	
	// can't go any further if read was too small
	if (i_bytes < UNIVERSAL_HEADER_BYTES) {
		fps->universal_header = NULL;
		if (!(behavior_on_fail & SUPPRESS_ERROR_OUTPUT)) {
			UTF8_fprintf(stderr, "Warning: fewer than UNIVERSAL_HEADER_BYTES read from file \"%s\" [function \"%s\", line %d]\n", fps->full_file_name, __FUNCTION__, __LINE__);
			if (behavior_on_fail & RETURN_ON_FAIL)
				(void) fprintf(stderr, "\t=> returning FILE_PROCESSING_STRUCT *\n\n");
			else if (behavior_on_fail & EXIT_ON_FAIL)
				(void) fprintf(stderr, "\t=> exiting program\n\n");
		}
		if (behavior_on_fail & RETURN_ON_FAIL)
			return(fps);
		else if (behavior_on_fail & EXIT_ON_FAIL)
			exit(1);
	}
	
	// cast universal header
	fps->universal_header = (UNIVERSAL_HEADER *) fps->raw_data;
	
	// set file type code
	file_type_string_int = (ui4 *) fps->universal_header->file_type_string;
	fps->file_type_code = *file_type_string_int;
	
	// process password data
       if (fps->password_data == NULL) {
		if (password_data == NULL)
			fps->password_data = process_password_data(password, NULL, NULL, fps->universal_header);
		else
			fps->password_data = password_data;
	}
        
    // CRCs
    if (MEF_globals->CRC_mode & (CRC_VALIDATE | CRC_VALIDATE_ON_INPUT)) {
        if (fps->directives.io_bytes == FPS_FULL_FILE) {
            CRC_result = CRC_validate(fps->raw_data + UNIVERSAL_HEADER_BYTES, fps->raw_data_bytes - UNIVERSAL_HEADER_BYTES, fps->universal_header->body_CRC);
            if (CRC_result == MEF_TRUE)
            {
                if (MEF_globals->verbose == MEF_TRUE)
                    UTF8_printf("Body CRC is valid in file \"%s\".\n", fps->full_file_name);
            }
            else
                UTF8_fprintf(stderr, "Warning: body CRC is invalid in file \"%s\".\n", fps->full_file_name);
            CRC_result = CRC_validate(fps->raw_data + CRC_BYTES, UNIVERSAL_HEADER_BYTES - CRC_BYTES, fps->universal_header->header_CRC);
            if (CRC_result == MEF_TRUE)
            {
                if (MEF_globals->verbose == MEF_TRUE)
                    UTF8_printf("Header CRC is valid in file \"%s\".\n", fps->full_file_name);
            }
            else
                UTF8_fprintf(stderr, "Warning: header CRC is invalid in file \"%s\".\n", fps->full_file_name);
        }
    }
	
	// if just reading UNIVERSAL HEADER (e.g. large files like records or segment data files)
        if (fps->directives.io_bytes == UNIVERSAL_HEADER_BYTES) {
        	if (MEF_globals->verbose == MEF_TRUE)
			show_file_processing_struct(fps);
		offset_universal_header_times(fps, RTO_INPUT_ACTION);
		return(fps);
	}

	// cast rest of data
    data_ptr = NULL;
	if (fps->raw_data != NULL)
		data_ptr = (ui1 *)fps->raw_data + UNIVERSAL_HEADER_BYTES;
	switch (fps->file_type_code) {
		case TIME_SERIES_INDICES_FILE_TYPE_CODE:
			fps->time_series_indices = (TIME_SERIES_INDEX *) data_ptr;
			break;
		case TIME_SERIES_DATA_FILE_TYPE_CODE:
			fps->RED_blocks = (ui1 *) data_ptr;
			break;
		case TIME_SERIES_METADATA_FILE_TYPE_CODE:
			fps->metadata.section_1 = (METADATA_SECTION_1 *) data_ptr;
			fps->metadata.time_series_section_2 = (TIME_SERIES_METADATA_SECTION_2 *) (fps->raw_data + METADATA_SECTION_2_OFFSET);
			fps->metadata.section_3 = (METADATA_SECTION_3 *) (fps->raw_data + METADATA_SECTION_3_OFFSET);
			break;
		case VIDEO_METADATA_FILE_TYPE_CODE:
			fps->metadata.section_1 = (METADATA_SECTION_1 *) data_ptr;
			fps->metadata.video_section_2 = (VIDEO_METADATA_SECTION_2 *) (fps->raw_data + METADATA_SECTION_2_OFFSET);
			fps->metadata.section_3 = (METADATA_SECTION_3 *) (fps->raw_data + METADATA_SECTION_3_OFFSET);
			break;
		case VIDEO_INDICES_FILE_TYPE_CODE:
			fps->video_indices = (VIDEO_INDEX *) data_ptr;
			break;
		case RECORD_DATA_FILE_TYPE_CODE:
			fps->records = (ui1 *) data_ptr;
			break;
		case RECORD_INDICES_FILE_TYPE_CODE:
			fps->record_indices = (RECORD_INDEX *) data_ptr;
			break;
		default:
			UTF8_fprintf(stderr, "Error: unrecognized type code in file \"%s\" [function \"%s\", line %d]\n", fps->full_file_name, __FUNCTION__, __LINE__);
                        if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL) {
                                (void) fprintf(stderr, "\t=> exiting program\n\n");
                                exit(1);
                        }
			return(NULL);
	}
	
	// decrypt encrypted data
	switch (fps->file_type_code) {
		case TIME_SERIES_METADATA_FILE_TYPE_CODE:
			decrypt_metadata(fps);  // also sets global recording time offsets
			break;
		case VIDEO_METADATA_FILE_TYPE_CODE:
			decrypt_metadata(fps);  // also sets global recording time offsets
			break;
		case RECORD_DATA_FILE_TYPE_CODE:
			decrypt_records(fps);    // also does time offsets, for efficiency
			break;
		default:
			break;
	}
	
	// offset times
	offset_universal_header_times(fps, RTO_INPUT_ACTION);
	switch (fps->file_type_code) {
		case TIME_SERIES_INDICES_FILE_TYPE_CODE:
			offset_time_series_index_times(fps, RTO_INPUT_ACTION);
			break;
		case VIDEO_INDICES_FILE_TYPE_CODE:
			offset_video_index_times(fps, RTO_INPUT_ACTION);
			break;
		case RECORD_INDICES_FILE_TYPE_CODE:
			offset_record_index_times(fps, RTO_INPUT_ACTION);
			break;
		default:
			break;
	}
	
	// show
	if (MEF_globals->verbose == MEF_TRUE)
		show_file_processing_struct(fps);
        
        
	return(fps);
}


SEGMENT	*read_MEF_segment(SEGMENT *segment, si1 *seg_path, si4 channel_type, si1 *password, PASSWORD_DATA *password_data, si1 read_time_series_data, si1 read_record_data)
{
	si1		full_file_name[MEF_FULL_FILE_NAME_BYTES];
	
	
	// allocate segment if not passed
	if (segment == NULL)
		segment = (SEGMENT *) e_calloc((size_t) 1, sizeof(SEGMENT), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        
	// get segment path & name
	extract_path_parts(seg_path, segment->path, segment->name, NULL);
	
	// get channel type
        if (channel_type == UNKNOWN_CHANNEL_TYPE)
		channel_type = channel_type_from_path(seg_path);
        segment->channel_type = channel_type;
	
	// read segment metadata
	switch (channel_type) {
		case TIME_SERIES_CHANNEL_TYPE:
			MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s/%s.%s", segment->path, segment->name, SEGMENT_DIRECTORY_TYPE_STRING, segment->name, TIME_SERIES_METADATA_FILE_TYPE_STRING);
			break;
		case VIDEO_CHANNEL_TYPE:
			MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s/%s.%s", segment->path, segment->name, SEGMENT_DIRECTORY_TYPE_STRING, segment->name, VIDEO_METADATA_FILE_TYPE_STRING);
			break;
		default:
			UTF8_fprintf(stderr, "Error: unrecognized type code in file \"%s\" [function \"%s\", line %d]\n", full_file_name, __FUNCTION__, __LINE__);
			if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL) {
				(void) fprintf(stderr, "\t=> exiting program\n\n");
				exit(1);
			}
			break;
	}
	segment->metadata_fps = read_MEF_file(NULL, full_file_name, password, password_data, NULL, USE_GLOBAL_BEHAVIOR);
	password_data = segment->metadata_fps->password_data;  // if password was passed we should have PASSWORD_DATA now
	
	// copy level UUID
	memcpy(segment->level_UUID, segment->metadata_fps->universal_header->level_UUID, UUID_BYTES);
        
        // copy session and channel names
        strcpy(segment->channel_name, segment->metadata_fps->universal_header->channel_name);
        strcpy(segment->session_name, segment->metadata_fps->universal_header->session_name);
	
	// read segment data
	switch (channel_type) {
		case TIME_SERIES_CHANNEL_TYPE:
			MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s/%s.%s", segment->path, segment->name, SEGMENT_DIRECTORY_TYPE_STRING, segment->name, TIME_SERIES_DATA_FILE_TYPE_STRING);
			segment->time_series_data_fps = allocate_file_processing_struct(0, TIME_SERIES_DATA_FILE_TYPE_CODE, NULL, NULL, 0);
			if (read_time_series_data == MEF_FALSE) {
				segment->time_series_data_fps->directives.io_bytes = UNIVERSAL_HEADER_BYTES;
				//segment->time_series_data_fps->directives.close_file = MEF_FALSE;
			}
			(void) read_MEF_file(segment->time_series_data_fps, full_file_name, password, password_data, NULL, USE_GLOBAL_BEHAVIOR);
			// update metadata if metadata conflicts with actual data
			if (segment->metadata_fps->metadata.time_series_section_2->number_of_blocks > segment->time_series_data_fps->universal_header->number_of_entries)
                                segment->metadata_fps->metadata.time_series_section_2->number_of_blocks = segment->time_series_data_fps->universal_header->number_of_entries;
			break;
		case VIDEO_CHANNEL_TYPE:
			// video channel data is video file
			break;
		default:
			UTF8_fprintf(stderr, "Error: unrecognized type code in file \"%s\" [function \"%s\", line %d]\n", full_file_name, __FUNCTION__, __LINE__);
			if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL) {
				(void) fprintf(stderr, "\t=> exiting program\n\n");
				exit(1);
			}
			break;
	}
	
	// read segment indices
	switch (channel_type) {
		case TIME_SERIES_CHANNEL_TYPE:
			MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s/%s.%s", segment->path, segment->name, SEGMENT_DIRECTORY_TYPE_STRING, segment->name, TIME_SERIES_INDICES_FILE_TYPE_STRING);
			segment->time_series_indices_fps = read_MEF_file(NULL, full_file_name, password, password_data, NULL, USE_GLOBAL_BEHAVIOR);
			// update metadata if metadata conflicts with actual data
			if (segment->metadata_fps->metadata.time_series_section_2->number_of_blocks > segment->time_series_indices_fps->universal_header->number_of_entries)
                                segment->metadata_fps->metadata.time_series_section_2->number_of_blocks = segment->time_series_indices_fps->universal_header->number_of_entries;
			break;
		case VIDEO_CHANNEL_TYPE:
			MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s/%s.%s", segment->path, segment->name, SEGMENT_DIRECTORY_TYPE_STRING, segment->name, VIDEO_INDICES_FILE_TYPE_STRING);
			segment->video_indices_fps = read_MEF_file(NULL, full_file_name, password, password_data, NULL, USE_GLOBAL_BEHAVIOR);
			break;
		default:
			UTF8_fprintf(stderr, "Error: unrecognized type code in file \"%s\" [function \"%s\", line %d]\n", full_file_name, __FUNCTION__, __LINE__);
			if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL) {
				(void) fprintf(stderr, "\t=> exiting program\n\n");
				exit(1);
			}
			break;
	}
	
	// read segment records
	MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s/%s.%s", segment->path, segment->name, SEGMENT_DIRECTORY_TYPE_STRING, segment->name, RECORD_INDICES_FILE_TYPE_STRING);
	segment->record_indices_fps = read_MEF_file(NULL, full_file_name, password, password_data, NULL, RETURN_ON_FAIL | SUPPRESS_ERROR_OUTPUT);
	if (segment->record_indices_fps != NULL) {
		// read segment record data
		MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s/%s.%s", segment->path, segment->name, SEGMENT_DIRECTORY_TYPE_STRING, segment->name, RECORD_DATA_FILE_TYPE_STRING);
		segment->record_data_fps = allocate_file_processing_struct(0, RECORD_DATA_FILE_TYPE_CODE, NULL, NULL, 0);
		if (read_record_data == MEF_FALSE) {
			segment->record_data_fps->directives.io_bytes = UNIVERSAL_HEADER_BYTES;
		}
		(void) read_MEF_file(segment->record_data_fps, full_file_name, password, password_data, NULL, RETURN_ON_FAIL | SUPPRESS_ERROR_OUTPUT);
		if (segment->record_data_fps == NULL)
			UTF8_fprintf(stderr, "%s() Warning: Segment record indices file, but no segment record data file (\"%s\") in segment directory\n\n", __FUNCTION__, full_file_name);
	}
	
	
	return(segment);
}


SESSION	*read_MEF_session(SESSION *session, si1 *sess_path, si1 *password, PASSWORD_DATA *password_data, si1 read_time_series_data, si1 read_record_data)
{
	si4				i, n_channels;
	si1				full_file_name[MEF_FULL_FILE_NAME_BYTES], **channel_names;
	CHANNEL				*chan;
	METADATA_SECTION_1		*smd1, *cmd1;
	TIME_SERIES_METADATA_SECTION_2	*ctmd, *stmd;
	VIDEO_METADATA_SECTION_2	*cvmd, *svmd;
	METADATA_SECTION_3		*smd3, *cmd3;
	FILE_PROCESSING_STRUCT		*temp_fps;
	
	
	// allocate session if not passed
	if (session == NULL)
		session = (SESSION *) e_calloc((size_t) 1, sizeof(SESSION), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        
	// get session path & name
	extract_path_parts(sess_path, session->path, session->name, NULL);
	
	// initialize
	session->maximum_number_of_records = 0;
	session->maximum_record_bytes = 0;
	bzero(session->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
	#ifdef _WIN32
		session->earliest_start_time = LLONG_MAX;
		session->latest_end_time = LLONG_MIN;
	#else
		session->earliest_start_time = LONG_MAX;
		session->latest_end_time = LONG_MIN;
	#endif
	
	// loop over time series channels
	channel_names = generate_file_list(NULL, &n_channels, sess_path, TIME_SERIES_CHANNEL_DIRECTORY_TYPE_STRING);
	session->time_series_channels = (CHANNEL *) e_calloc((size_t) n_channels, sizeof(CHANNEL), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	for (i = 0; i < n_channels; ++i) {
		(void) read_MEF_channel(session->time_series_channels + i, channel_names[i], TIME_SERIES_CHANNEL_TYPE, password, password_data, read_time_series_data, read_record_data);
		if ((password_data == NULL) && (session->time_series_channels[i].number_of_segments > 0))
            password_data = session->time_series_channels[i].segments[0].metadata_fps->password_data;
		free(channel_names[i]);
	}
	session->number_of_time_series_channels = n_channels;
	free(channel_names);

	// loop over video channels
	channel_names = generate_file_list(NULL, &n_channels, sess_path, VIDEO_CHANNEL_DIRECTORY_TYPE_STRING);
	session->video_channels = (CHANNEL *) e_calloc((size_t) n_channels, sizeof(CHANNEL), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	for (i = 0; i < n_channels; ++i) {
		(void) read_MEF_channel(session->video_channels + i, channel_names[i], VIDEO_CHANNEL_TYPE, password, password_data, read_time_series_data, read_record_data);
		if (password_data == NULL)
			password_data = session->video_channels[i].segments[0].metadata_fps->password_data;
		free(channel_names[i]);
	}
	session->number_of_video_channels = n_channels;
	free(channel_names);
	
	// fill in session metadata: times series channels
	if (session->number_of_time_series_channels > 0) {
		if (session->time_series_metadata.section_1 == NULL)
			session->time_series_metadata.section_1 = (METADATA_SECTION_1 *) e_calloc((size_t) 1, sizeof(METADATA_SECTION_1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		if (session->time_series_metadata.time_series_section_2 == NULL)
			session->time_series_metadata.time_series_section_2 = (TIME_SERIES_METADATA_SECTION_2 *) e_calloc((size_t) 1, sizeof(TIME_SERIES_METADATA_SECTION_2), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		if (session->time_series_metadata.section_3 == NULL)
			session->time_series_metadata.section_3 = (METADATA_SECTION_3 *) e_calloc((size_t) 1, sizeof(METADATA_SECTION_3), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	}

	for (i = 0; i < session->number_of_time_series_channels; ++i) {
		chan = session->time_series_channels + i;
		cmd1 = chan->metadata.section_1;
		smd1 = session->time_series_metadata.section_1;
		ctmd = chan->metadata.time_series_section_2;
		stmd = session->time_series_metadata.time_series_section_2;
		cmd3 = chan->metadata.section_3;
		smd3 = session->time_series_metadata.section_3;
		if (i == 0) {
			memcpy(smd1, cmd1, METADATA_SECTION_1_BYTES);
			memcpy(stmd, ctmd, METADATA_SECTION_2_BYTES);
			memcpy(smd3, cmd3, METADATA_SECTION_3_BYTES);
			if (ABS(session->earliest_start_time) > ABS(chan->earliest_start_time))
				session->earliest_start_time = chan->earliest_start_time;
			if (ABS(session->latest_end_time) < ABS(chan->latest_end_time))
				session->latest_end_time = chan->latest_end_time;
			if (chan->record_data_fps != NULL) {
				if (session->maximum_number_of_records < chan->maximum_number_of_records)
					session->maximum_number_of_records = chan->maximum_number_of_records;
				if (session->maximum_record_bytes < chan->maximum_record_bytes)
					session->maximum_record_bytes = chan->maximum_record_bytes;
			}
			if (strlen(session->anonymized_name) == 0)
				MEF_strncpy(session->anonymized_name, chan->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
			else if (strcmp(session->anonymized_name, chan->anonymized_name))
				bzero(session->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
			continue;
		}
		// universal header
		if (ABS(chan->earliest_start_time) < ABS(session->earliest_start_time))
			session->earliest_start_time = chan->earliest_start_time;
		if (ABS(chan->latest_end_time) > ABS(session->latest_end_time))
			session->latest_end_time = chan->latest_end_time;
		if (strcmp(session->anonymized_name, chan->anonymized_name))
			bzero(session->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
		if (chan->record_data_fps != NULL) {
			if (chan->record_data_fps->universal_header->number_of_entries > session->maximum_number_of_records)
				session->maximum_number_of_records = chan->record_data_fps->universal_header->number_of_entries;
			if (chan->record_data_fps->universal_header->maximum_entry_size > session->maximum_record_bytes)
				session->maximum_record_bytes = chan->record_data_fps->universal_header->maximum_entry_size;
		}
		// section 1
		if (smd1->section_2_encryption != cmd1->section_2_encryption)
			smd1->section_2_encryption = ENCRYPTION_LEVEL_NO_ENTRY;
		if (smd1->section_3_encryption != cmd1->section_3_encryption)
			smd1->section_3_encryption = ENCRYPTION_LEVEL_NO_ENTRY;
		if (memcmp(smd1->protected_region, cmd1->protected_region, METADATA_SECTION_1_PROTECTED_REGION_BYTES))
			bzero(smd1->protected_region, METADATA_SECTION_1_PROTECTED_REGION_BYTES);
		if (memcmp(smd1->discretionary_region, cmd1->discretionary_region, METADATA_SECTION_1_DISCRETIONARY_REGION_BYTES))
			bzero(smd1->discretionary_region, METADATA_SECTION_1_DISCRETIONARY_REGION_BYTES);
		// section 2
		if (strcmp(ctmd->channel_description, stmd->channel_description))
			bzero(stmd->channel_description, METADATA_CHANNEL_DESCRIPTION_BYTES);
		if (strcmp(ctmd->session_description, stmd->session_description))
			bzero(stmd->session_description, METADATA_SESSION_DESCRIPTION_BYTES);
		if (session->latest_end_time == UUTC_NO_ENTRY || session->earliest_start_time == UUTC_NO_ENTRY)
			stmd->recording_duration = UUTC_NO_ENTRY;
		else
			stmd->recording_duration = ABS(session->latest_end_time) - ABS(session->earliest_start_time) + 1;
		if (strcmp(stmd->reference_description, ctmd->reference_description))
			bzero(stmd->reference_description, METADATA_CHANNEL_DESCRIPTION_BYTES);
		if (ctmd->acquisition_channel_number != stmd->acquisition_channel_number)
			stmd->acquisition_channel_number = TIME_SERIES_METADATA_ACQUISITION_CHANNEL_NUMBER_NO_ENTRY;
		if (ctmd->sampling_frequency != stmd->sampling_frequency)
			stmd->sampling_frequency = TIME_SERIES_METADATA_SAMPLING_FREQUENCY_NO_ENTRY;
		if (ctmd->low_frequency_filter_setting != stmd->low_frequency_filter_setting)
			stmd->low_frequency_filter_setting = TIME_SERIES_METADATA_LOW_FREQUENCY_FILTER_SETTING_NO_ENTRY;
		if (ctmd->high_frequency_filter_setting != stmd->high_frequency_filter_setting)
			stmd->high_frequency_filter_setting = TIME_SERIES_METADATA_HIGH_FREQUENCY_FILTER_SETTING_NO_ENTRY;
		if (ctmd->notch_filter_frequency_setting != stmd->notch_filter_frequency_setting)
			stmd->notch_filter_frequency_setting = TIME_SERIES_METADATA_NOTCH_FILTER_FREQUENCY_SETTING_NO_ENTRY;
		if (ctmd->AC_line_frequency != stmd->AC_line_frequency)
			stmd->AC_line_frequency = TIME_SERIES_METADATA_AC_LINE_FREQUENCY_NO_ENTRY;
		if (ctmd->units_conversion_factor != stmd->units_conversion_factor)
			stmd->units_conversion_factor = TIME_SERIES_METADATA_UNITS_CONVERSION_FACTOR_NO_ENTRY;
		if (strcmp(ctmd->units_description, stmd->units_description))
                	bzero(stmd->units_description, TIME_SERIES_METADATA_UNITS_DESCRIPTION_BYTES);
		if (ctmd->maximum_native_sample_value > stmd->maximum_native_sample_value)
			stmd->maximum_native_sample_value = ctmd->maximum_native_sample_value;
		if (ctmd->minimum_native_sample_value < stmd->minimum_native_sample_value)
			stmd->minimum_native_sample_value = ctmd->minimum_native_sample_value;
                if (ctmd->start_sample < stmd->start_sample)
                        stmd->start_sample = ctmd->start_sample;
		if (ctmd->number_of_samples > stmd->number_of_samples)
			stmd->number_of_samples = ctmd->number_of_samples;
		if (ctmd->number_of_blocks > stmd->number_of_blocks)
			stmd->number_of_blocks = ctmd->number_of_blocks;
		if (ctmd->maximum_block_samples > stmd->maximum_block_samples)
			stmd->maximum_block_samples = ctmd->maximum_block_samples;
		if (ctmd->maximum_difference_bytes > stmd->maximum_difference_bytes)
			stmd->maximum_difference_bytes = ctmd->maximum_difference_bytes;
		if (ctmd->block_interval != stmd->block_interval)
			stmd->block_interval = TIME_SERIES_METADATA_BLOCK_INTERVAL_NO_ENTRY;
		if (ctmd->number_of_discontinuities > stmd->number_of_discontinuities)
			stmd->number_of_discontinuities = ctmd->number_of_discontinuities;
		if (ctmd->maximum_contiguous_blocks > stmd->maximum_contiguous_blocks)
			stmd->maximum_contiguous_blocks = ctmd->maximum_contiguous_blocks;
		if (ctmd->maximum_contiguous_block_bytes > stmd->maximum_contiguous_block_bytes)
			stmd->maximum_contiguous_block_bytes = ctmd->maximum_contiguous_block_bytes;
		if (ctmd->maximum_contiguous_samples > stmd->maximum_contiguous_samples)
			stmd->maximum_contiguous_samples = ctmd->maximum_contiguous_samples;
		if (memcmp(ctmd->protected_region, stmd->protected_region, TIME_SERIES_METADATA_SECTION_2_PROTECTED_REGION_BYTES))
			bzero(stmd->protected_region, TIME_SERIES_METADATA_SECTION_2_PROTECTED_REGION_BYTES);
		if (memcmp(ctmd->discretionary_region, stmd->discretionary_region, TIME_SERIES_METADATA_SECTION_2_DISCRETIONARY_REGION_BYTES))
			bzero(stmd->discretionary_region, TIME_SERIES_METADATA_SECTION_2_DISCRETIONARY_REGION_BYTES);
		// section 3
		if (smd3->recording_time_offset != cmd3->recording_time_offset)
			smd3->recording_time_offset = UUTC_NO_ENTRY;
		if (smd3->DST_start_time != cmd3->DST_start_time)
			smd3->DST_start_time = UUTC_NO_ENTRY;
		if (smd3->DST_end_time != cmd3->DST_end_time)
			smd3->DST_end_time = UUTC_NO_ENTRY;
		if (smd3->GMT_offset != cmd3->GMT_offset)
			smd3->GMT_offset = GMT_OFFSET_NO_ENTRY;
		if (strcmp(smd3->subject_name_1, cmd3->subject_name_1))
			bzero(smd3->subject_name_1, METADATA_SUBJECT_NAME_BYTES);
		if (strcmp(smd3->subject_name_2, cmd3->subject_name_2))
			bzero(smd3->subject_name_2, METADATA_SUBJECT_NAME_BYTES);
		if (strcmp(smd3->subject_ID, cmd3->subject_ID))
			bzero(smd3->subject_ID, METADATA_SUBJECT_ID_BYTES);
		if (strcmp(smd3->recording_location, cmd3->recording_location))
			bzero(smd3->recording_location, METADATA_RECORDING_LOCATION_BYTES);
		if (memcmp(smd3->protected_region, cmd3->protected_region, METADATA_SECTION_3_PROTECTED_REGION_BYTES))
			bzero(smd3->protected_region, METADATA_SECTION_3_PROTECTED_REGION_BYTES);
		if (memcmp(smd3->discretionary_region, cmd3->discretionary_region, METADATA_SECTION_3_DISCRETIONARY_REGION_BYTES))
			bzero(smd3->discretionary_region, METADATA_SECTION_3_DISCRETIONARY_REGION_BYTES);
		
	}

	// fill in session metadata: video channels
	if (session->number_of_video_channels > 0) {
		if (session->video_metadata.section_1 == NULL)
			session->video_metadata.section_1 = (METADATA_SECTION_1 *) e_calloc((size_t) 1, sizeof(METADATA_SECTION_1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		if (session->video_metadata.video_section_2 == NULL)
			session->video_metadata.video_section_2 = (VIDEO_METADATA_SECTION_2 *) e_calloc((size_t) 1, sizeof(VIDEO_METADATA_SECTION_2), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		if (session->video_metadata.section_3 == NULL)
			session->video_metadata.section_3 = (METADATA_SECTION_3 *) e_calloc((size_t) 1, sizeof(METADATA_SECTION_3), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	}
	
	for (i = 0; i < session->number_of_video_channels; ++i) {
		chan = session->video_channels + i;
		cmd1 = chan->metadata.section_1;
		smd1 = session->time_series_metadata.section_1;
		cvmd = chan->metadata.video_section_2;
		svmd = session->video_metadata.video_section_2;
		cmd3 = chan->metadata.section_3;
		smd3 = session->time_series_metadata.section_3;
		if (i == 0) {
			memcpy(smd1, cmd1, METADATA_SECTION_1_BYTES);
			memcpy(svmd, cvmd, METADATA_SECTION_2_BYTES);
			memcpy(smd3, cmd3, METADATA_SECTION_3_BYTES);
			if (ABS(session->earliest_start_time) > ABS(chan->earliest_start_time))
				session->earliest_start_time = chan->earliest_start_time;
			if (ABS(session->latest_end_time) < ABS(chan->latest_end_time))
				session->latest_end_time = chan->latest_end_time;
			if (chan->record_data_fps != NULL) {
				if (session->maximum_number_of_records < chan->maximum_number_of_records)
					session->maximum_number_of_records = chan->maximum_number_of_records;
				if (session->maximum_record_bytes < chan->maximum_record_bytes)
					session->maximum_record_bytes = chan->maximum_record_bytes;
			}
			if (strlen(session->anonymized_name) == 0)
				MEF_strncpy(session->anonymized_name, chan->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
			else if (strcmp(session->anonymized_name, chan->anonymized_name))
				bzero(session->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
			continue;
		}
		// universal header
		if (ABS(chan->earliest_start_time) < ABS(session->earliest_start_time))
			session->earliest_start_time = chan->earliest_start_time;
		if (ABS(chan->latest_end_time) > ABS(session->latest_end_time))
			session->latest_end_time = chan->latest_end_time;
		if (strcmp(session->anonymized_name, chan->anonymized_name))
			bzero(session->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
		if (chan->record_data_fps != NULL) {
			if (chan->record_data_fps->universal_header->number_of_entries > session->maximum_number_of_records)
				session->maximum_number_of_records = chan->record_data_fps->universal_header->number_of_entries;
			if (chan->record_data_fps->universal_header->maximum_entry_size > session->maximum_record_bytes)
				session->maximum_record_bytes = chan->record_data_fps->universal_header->maximum_entry_size;
		}
		// section 1
		if (smd1->section_2_encryption != cmd1->section_2_encryption)
			smd1->section_2_encryption = ENCRYPTION_LEVEL_NO_ENTRY;
		if (smd1->section_3_encryption != cmd1->section_3_encryption)
			smd1->section_3_encryption = ENCRYPTION_LEVEL_NO_ENTRY;
		if (memcmp(smd1->protected_region, cmd1->protected_region, METADATA_SECTION_1_PROTECTED_REGION_BYTES))
			bzero(smd1->protected_region, METADATA_SECTION_1_PROTECTED_REGION_BYTES);
		if (memcmp(smd1->discretionary_region, cmd1->discretionary_region, METADATA_SECTION_1_DISCRETIONARY_REGION_BYTES))
			bzero(smd1->discretionary_region, METADATA_SECTION_1_DISCRETIONARY_REGION_BYTES);
		// section 2
		if (strcmp(cvmd->channel_description, svmd->channel_description))
			bzero(svmd->channel_description, METADATA_CHANNEL_DESCRIPTION_BYTES);
		if (strcmp(cvmd->session_description, svmd->session_description))
			bzero(svmd->session_description, METADATA_SESSION_DESCRIPTION_BYTES);
		if (session->latest_end_time == UUTC_NO_ENTRY || session->earliest_start_time == UUTC_NO_ENTRY)
			svmd->recording_duration = UUTC_NO_ENTRY;
		else
			svmd->recording_duration = ABS(session->latest_end_time) - ABS(session->earliest_start_time) + 1;
		if (cvmd->horizontal_resolution != svmd->horizontal_resolution)
			svmd->horizontal_resolution = VIDEO_METADATA_HORIZONTAL_RESOLUTION_NO_ENTRY;
		if (cvmd->vertical_resolution != svmd->vertical_resolution)
			svmd->vertical_resolution = VIDEO_METADATA_VERTICAL_RESOLUTION_NO_ENTRY;
		if (cvmd->frame_rate != svmd->frame_rate)
			svmd->vertical_resolution = VIDEO_METADATA_FRAME_RATE_NO_ENTRY;
		if (cvmd->number_of_clips > svmd->number_of_clips)
			svmd->number_of_clips = cvmd->number_of_clips;
		if (cvmd->maximum_clip_bytes > svmd->maximum_clip_bytes)
			svmd->maximum_clip_bytes = cvmd->maximum_clip_bytes;
		if (strcmp(cvmd->video_format, svmd->video_format))
                	bzero(svmd->video_format, VIDEO_METADATA_VIDEO_FORMAT_BYTES);
		if (svmd->video_file_CRC != cvmd->video_file_CRC)
			svmd->video_file_CRC = CRC_NO_ENTRY;
		if (memcmp(cvmd->protected_region, svmd->protected_region, VIDEO_METADATA_SECTION_2_PROTECTED_REGION_BYTES))
			bzero(svmd->protected_region, VIDEO_METADATA_SECTION_2_PROTECTED_REGION_BYTES);
		if (memcmp(cvmd->discretionary_region, svmd->discretionary_region, VIDEO_METADATA_SECTION_2_DISCRETIONARY_REGION_BYTES))
			bzero(svmd->discretionary_region, VIDEO_METADATA_SECTION_2_DISCRETIONARY_REGION_BYTES);
		// section 3
		if (smd3->recording_time_offset != cmd3->recording_time_offset)
			smd3->recording_time_offset = UUTC_NO_ENTRY;
		if (smd3->DST_start_time != cmd3->DST_start_time)
			smd3->DST_start_time = UUTC_NO_ENTRY;
		if (smd3->DST_end_time != cmd3->DST_end_time)
			smd3->DST_end_time = UUTC_NO_ENTRY;
		if (smd3->GMT_offset != cmd3->GMT_offset)
			smd3->GMT_offset = GMT_OFFSET_NO_ENTRY;
		if (strcmp(smd3->subject_name_1, cmd3->subject_name_1))
			bzero(smd3->subject_name_1, METADATA_SUBJECT_NAME_BYTES);
		if (strcmp(smd3->subject_name_2, cmd3->subject_name_2))
			bzero(smd3->subject_name_2, METADATA_SUBJECT_NAME_BYTES);
		if (strcmp(smd3->subject_ID, cmd3->subject_ID))
			bzero(smd3->subject_ID, METADATA_SUBJECT_ID_BYTES);
		if (strcmp(smd3->recording_location, cmd3->recording_location))
			bzero(smd3->recording_location, METADATA_RECORDING_LOCATION_BYTES);
		if (memcmp(smd3->protected_region, smd3->protected_region, METADATA_SECTION_3_PROTECTED_REGION_BYTES))
			bzero(smd3->protected_region, METADATA_SECTION_3_PROTECTED_REGION_BYTES);
		if (memcmp(smd3->discretionary_region, cmd3->discretionary_region, METADATA_SECTION_3_DISCRETIONARY_REGION_BYTES))
			bzero(smd3->discretionary_region, METADATA_SECTION_3_DISCRETIONARY_REGION_BYTES);
	}

	// read session record indices if present
	MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s/%s.%s", session->path, session->name, SESSION_DIRECTORY_TYPE_STRING, session->name, RECORD_INDICES_FILE_TYPE_STRING);
	session->record_indices_fps = read_MEF_file(NULL, full_file_name, password, password_data, NULL, RETURN_ON_FAIL | SUPPRESS_ERROR_OUTPUT);
    if (session->record_indices_fps != NULL) {
		if (password_data == NULL)
			password_data = session->record_indices_fps->password_data;
		// copy level UUID
		memcpy(session->level_UUID, session->record_indices_fps->universal_header->level_UUID, UUID_BYTES);
        // read session records data
        MEF_snprintf(full_file_name, MEF_FULL_FILE_NAME_BYTES, "%s/%s.%s/%s.%s", session->path, session->name, SESSION_DIRECTORY_TYPE_STRING, session->name, RECORD_DATA_FILE_TYPE_STRING);
        session->record_data_fps = allocate_file_processing_struct(0, RECORD_DATA_FILE_TYPE_CODE, NULL, NULL, 0);
        if (read_record_data == MEF_FALSE) {
                session->record_data_fps->directives.io_bytes = UNIVERSAL_HEADER_BYTES;
        }
        (void) read_MEF_file(session->record_data_fps, full_file_name, password, password_data, NULL, RETURN_ON_FAIL | SUPPRESS_ERROR_OUTPUT);
        if (session->record_data_fps == NULL)
            UTF8_fprintf(stderr, "%s() Warning: Session record indices file, but no session records data file (\"%s\") in session directory\n\n", __FUNCTION__, full_file_name);
		
		if (session->record_data_fps->universal_header->number_of_entries > session->maximum_number_of_records)
			session->maximum_number_of_records = session->record_data_fps->universal_header->number_of_entries;
		if (session->maximum_record_bytes < session->record_data_fps->universal_header->maximum_entry_size)
			session->maximum_record_bytes = session->record_data_fps->universal_header->maximum_entry_size;
		if (ABS(session->record_data_fps->universal_header->start_time) < ABS(session->earliest_start_time))
			session->earliest_start_time = session->record_data_fps->universal_header->start_time;
		if (ABS(session->latest_end_time) < session->record_data_fps->universal_header->end_time)
			session->latest_end_time = session->record_data_fps->universal_header->end_time;
		MEF_strncpy(session->anonymized_name, session->record_data_fps->universal_header->anonymized_name, UNIVERSAL_HEADER_ANONYMIZED_NAME_BYTES);
	}
	
	if (MEF_globals->verbose == MEF_TRUE) {
		if (session->number_of_time_series_channels > 0) {
			printf("------------ Session Time Series Metadata --------------\n");
			temp_fps = allocate_file_processing_struct(0, TIME_SERIES_METADATA_FILE_TYPE_CODE, NULL, NULL, 0);
			temp_fps->metadata = session->time_series_metadata;
			temp_fps->password_data = password_data;
			show_metadata(temp_fps);
			free_file_processing_struct(temp_fps);
		}
		if (session->number_of_video_channels > 0) {
			printf("--------------- Session Video Metadata -----------------\n");
			temp_fps = allocate_file_processing_struct(0, VIDEO_METADATA_FILE_TYPE_CODE, NULL, NULL, 0);
			temp_fps->metadata = session->video_metadata;
			temp_fps->password_data = password_data;
			show_metadata(temp_fps);
			free_file_processing_struct(temp_fps);
		}
	}

	
	return(session);
}


si4	reallocate_file_processing_struct(FILE_PROCESSING_STRUCT *fps, si8 raw_data_bytes)
{
	void	*data_ptr;
	
	
	// reallocate
	fps->raw_data = (ui1 *) e_realloc((void *) fps->raw_data, (size_t) raw_data_bytes, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	// zero additional memory
	if (raw_data_bytes > fps->raw_data_bytes)
		bzero(fps->raw_data + fps->raw_data_bytes, raw_data_bytes - fps->raw_data_bytes);
        fps->raw_data_bytes = raw_data_bytes;
	
	// reset universal header pointer
	if (raw_data_bytes >= UNIVERSAL_HEADER_BYTES)
		fps->universal_header = (UNIVERSAL_HEADER *) fps->raw_data; // all files start with universal header
	else
		return(0);
        
	// set appropriate pointers
	data_ptr = NULL;
	if (fps->raw_data != NULL)
		data_ptr = (ui1 *)fps->raw_data + UNIVERSAL_HEADER_BYTES;

        switch (fps->file_type_code) {
                case NO_TYPE_CODE:
                        break;
                case TIME_SERIES_INDICES_FILE_TYPE_CODE:
                        fps->time_series_indices = (TIME_SERIES_INDEX *) data_ptr;
                        break;
		case TIME_SERIES_METADATA_FILE_TYPE_CODE:
			fps->metadata.section_1 = (METADATA_SECTION_1 *) data_ptr;
			fps->metadata.time_series_section_2 = (TIME_SERIES_METADATA_SECTION_2 *) (fps->raw_data + METADATA_SECTION_2_OFFSET);
			fps->metadata.section_3 = (METADATA_SECTION_3 *) (fps->raw_data + METADATA_SECTION_3_OFFSET);
			break;
		case TIME_SERIES_DATA_FILE_TYPE_CODE:
			fps->RED_blocks = (ui1 *) data_ptr;
			break;
                case VIDEO_METADATA_FILE_TYPE_CODE:
                        fps->metadata.section_1 = (METADATA_SECTION_1 *) data_ptr;
                        fps->metadata.video_section_2 = (VIDEO_METADATA_SECTION_2 *) (fps->raw_data + METADATA_SECTION_2_OFFSET);
                        fps->metadata.section_3 = (METADATA_SECTION_3 *) (fps->raw_data + METADATA_SECTION_3_OFFSET);
                        break;
                case VIDEO_INDICES_FILE_TYPE_CODE:
                        fps->video_indices = (VIDEO_INDEX *) data_ptr;
                	break;
                case RECORD_DATA_FILE_TYPE_CODE:
                        fps->records = (ui1 *) data_ptr;
                        break;
                case RECORD_INDICES_FILE_TYPE_CODE:
                        fps->record_indices = (RECORD_INDEX *) data_ptr;
                        break;
                default:
                        fprintf(stderr, "Error: unrecognized type code \"0x%x\" [function \"%s\", line %d]\n", fps->file_type_code, __FUNCTION__, __LINE__);
                        if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL) {
                                (void) fprintf(stderr, "\t=> exiting program\n\n");
        			exit(1);
                        }
                        return(-1);
        }
	
        
	return(0);
}


/*************************************************************************/
/******************************  RED FUNCTIONS  **************************/
/*************************************************************************/


RED_PROCESSING_STRUCT	*RED_allocate_processing_struct(si8 original_data_size, si8 compressed_data_size, si8 decompressed_data_size, si8 difference_buffer_size, si8 detrended_buffer_size, si8 scaled_buffer_size, PASSWORD_DATA *password_data)
{
        RED_PROCESSING_STRUCT	*rps;
        
	
	// allocate
        rps = (RED_PROCESSING_STRUCT *) e_calloc((size_t) 1, sizeof(RED_PROCESSING_STRUCT), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        
        if (original_data_size)
		rps->original_ptr = rps->original_data = (si4 *) e_calloc((size_t) original_data_size, sizeof(si4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
        if (compressed_data_size) {
		rps->compressed_data = (ui1 *) e_calloc((size_t) compressed_data_size, sizeof(ui1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
        	rps->block_header = (RED_BLOCK_HEADER *) rps->compressed_data;
	}
		
	if (decompressed_data_size)
		rps->decompressed_ptr = rps->decompressed_data = (si4 *) e_calloc((size_t) decompressed_data_size, sizeof(si4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	if (difference_buffer_size)
		rps->difference_buffer = (si1 *) e_calloc((size_t) difference_buffer_size, sizeof(ui1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	if (detrended_buffer_size)
		rps->detrended_buffer = (si4 *) e_calloc((size_t) detrended_buffer_size, sizeof(si4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	if (scaled_buffer_size)
		rps->scaled_buffer = (si4 *) e_calloc((size_t) scaled_buffer_size, sizeof(si4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	// assign password data
	rps->password_data = password_data;
        
	// set directives defaults
	rps->directives.return_lossy_data = RED_RETURN_LOSSY_DATA_DEFAULT;
	rps->directives.encryption_level = RED_ENCRYPTION_LEVEL_DEFAULT;
	rps->directives.detrend_data = RED_DETREND_DATA_DEFAULT;
	rps->directives.discontinuity = RED_DISCONTINUITY_DEFAULT;
	rps->directives.reset_discontinuity = RED_RESET_DISCONTINUITY_DEFAULT;
        rps->directives.require_normality = RED_REQUIRE_NORMALITY_DEFAULT;
        rps->directives.normal_correlation = RED_NORMAL_CORRELATION_DEFAULT;
	
	// set compression defaults
	rps->compression.mode = RED_COMPRESSION_MODE_DEFAULT;
	rps->compression.goal_compression_ratio = RED_GOAL_COMPRESSION_RATIO_DEFAULT;
	rps->compression.goal_mean_residual_ratio = RED_GOAL_MEAN_RESIDUAL_RATIO_DEFAULT;
	rps->compression.goal_tolerance = RED_GOAL_TOLERANCE_DEFAULT;
	rps->compression.maximum_rounds_per_block = RED_MAXIMUM_ROUNDS_PER_BLOCK_DEFAULT;
	
        
        return(rps);
}


inline sf8	RED_calculate_mean_residual_ratio(si4 *original_data, si4 *lossy_data, ui4 n_samps)
{
	sf8	sum, mrr, diff, r;
	si8	i;
	
	
	sum = (sf8) 0.0;
	for (i = n_samps; i--;) {
		if (*original_data) {
			diff = (sf8) (*original_data - *lossy_data++);
                        r = diff / (sf8) *original_data++;
			sum += ABS(r);
		} else {
			--n_samps;
			++original_data;
			++lossy_data;
		}
	}
	
	if (sum == (sf8) 0.0)
		mrr = (sf8) 0.0;
	else
		mrr = sum / (sf8) n_samps;
	
	
	return(mrr);
}


si1	RED_check_RPS_allocation(RED_PROCESSING_STRUCT *rps)
{
        si1	error = MEF_FALSE;
	si1	need_compressed_data = MEF_FALSE;
	si1	need_decompressed_data = MEF_FALSE;
        si1	need_original_data = MEF_FALSE;
	si1	need_detrended_buffer = MEF_FALSE;
	si1	need_scaled_buffer = MEF_FALSE;
	si1	need_difference_buffer = MEF_FALSE;
	
        
        // decompression
        if (rps->compression.mode == RED_DECOMPRESSION) {
                need_compressed_data = MEF_TRUE;
                need_decompressed_data = MEF_TRUE;
                need_difference_buffer = MEF_TRUE;
        }
        
        // compression
        else {
                need_compressed_data = MEF_TRUE;
                need_original_data = MEF_TRUE;
                need_difference_buffer = MEF_TRUE;
                
                // lossy compression
                if (rps->compression.mode > RED_LOSSLESS_COMPRESSION) {
                        need_scaled_buffer = MEF_TRUE;
                        if (rps->directives.detrend_data == MEF_TRUE)
                                need_detrended_buffer = MEF_TRUE;
                        if (rps->compression.mode == RED_MEAN_RESIDUAL_RATIO)
                                need_decompressed_data = MEF_TRUE;
                        if (rps->directives.return_lossy_data == MEF_TRUE)
                                need_decompressed_data = MEF_TRUE;
                }
        }
        
        // check compressed_data
        if (need_compressed_data == MEF_TRUE && rps->compressed_data == NULL) {
                fprintf(stderr, "\"compressed_data\" is not allocated in the RED_PROCESSING_STRUCT [function %s, line %d]\n", __FUNCTION__, __LINE__);
		error = MEF_TRUE;
        }
	
        // check difference_buffer
        if (need_difference_buffer == MEF_TRUE && rps->difference_buffer == NULL) {
                fprintf(stderr, "\"difference_buffer\" is not allocated in the RED_PROCESSING_STRUCT [function %s, line %d]\n", __FUNCTION__, __LINE__);
		error = MEF_TRUE;
        }
	
        // check original_data
        if (need_original_data == MEF_TRUE && rps->original_data == NULL) {
                fprintf(stderr, "\"original_data\" is not allocated in the RED_PROCESSING_STRUCT [function %s, line %d]\n", __FUNCTION__, __LINE__);
		error = MEF_TRUE;
        }
	if (need_original_data == MEF_FALSE && rps->original_data != NULL) {
                fprintf(stderr, "\"original_data\" is needlessly allocated in the RED_PROCESSING_STRUCT => freeing [function %s, line %d]\n", __FUNCTION__, __LINE__);
		free(rps->original_data);
		rps->original_ptr = rps->original_data = NULL;
	}
        
        // check decompressed_data
        if (need_decompressed_data == MEF_TRUE && rps->decompressed_data == NULL) {
                fprintf(stderr, "\"decompressed_data\" is not allocated in the RED_PROCESSING_STRUCT [function %s, line %d]\n", __FUNCTION__, __LINE__);
		error = MEF_TRUE;
        }
	if (need_decompressed_data == MEF_FALSE && rps->decompressed_data != NULL) {
		fprintf(stderr, "\"decompressed_data\" is needlessly allocated in the RED_PROCESSING_STRUCT => freeing [function %s, line %d]\n", __FUNCTION__, __LINE__);
		free(rps->decompressed_data);
		rps->decompressed_ptr = rps->decompressed_data = NULL;
	}

        // check detrended_buffer
        if (need_detrended_buffer == MEF_TRUE && rps->detrended_buffer == NULL) {
                fprintf(stderr, "\"detrended_buffer\" is not allocated in the RED_PROCESSING_STRUCT [function %s, line %d]\n", __FUNCTION__, __LINE__);
		error = MEF_TRUE;
        }
	if (need_detrended_buffer == MEF_FALSE && rps->detrended_buffer != NULL) {
                fprintf(stderr, "\"detrended_buffer\" is needlessly allocated in the RED_PROCESSING_STRUCT => freeing [function %s, line %d]\n", __FUNCTION__, __LINE__);
		free(rps->detrended_buffer);
		rps->detrended_buffer = NULL;
	}
        
        // check scaled_buffer
        if (need_scaled_buffer == MEF_TRUE && rps->scaled_buffer == NULL) {
                fprintf(stderr, "\"scaled_buffer\" is not allocated in the RED_PROCESSING_STRUCT [function %s, line %d]\n", __FUNCTION__, __LINE__);
		error = MEF_TRUE;
        }
	if (need_scaled_buffer == MEF_FALSE && rps->scaled_buffer != NULL) {
                fprintf(stderr, "\"scaled_buffer\" is needlessly allocated in the RED_PROCESSING_STRUCT => freeing [function %s, line %d]\n", __FUNCTION__, __LINE__);
		free(rps->scaled_buffer);
		rps->scaled_buffer = NULL;
	}
	
	// error
	if (error == MEF_TRUE) {
		if (MEF_globals->behavior_on_fail & EXIT_ON_FAIL) {
			(void) fprintf(stderr, "\t=> exiting program\n\n");
			exit(1);
		}
	}
        

        return(error);
}


void 	RED_decode(RED_PROCESSING_STRUCT *rps)
{
        si1			*si1_p1, *si1_p2, *diff_buffer_p, CRC_valid;
        ui1			*ui1_p, *ib_p, in_byte, *scaled_counts, *key;
        si4			*si4_p, current_val;
	si8			i;
        ui4			cc, *cumulative_counts, low_bound, range, symbol;
        ui4			scaled_total_counts, temp_ui4, range_per_count, *ui4_p1, *ui4_p2;
        RED_BLOCK_HEADER	*block_header;
        
        
        // RED decompress from compressed_ptr to decompressed_ptr
	block_header = rps->block_header;
	
        // check CRC
	if (MEF_globals->CRC_mode & (CRC_VALIDATE | CRC_VALIDATE_ON_INPUT)) {
                CRC_valid = CRC_validate((ui1 *) block_header + CRC_BYTES, block_header->block_bytes - CRC_BYTES, block_header->block_CRC);
                if (CRC_valid == MEF_FALSE) {
                        (void) fprintf(stderr, "%c\n%s(): invalid RED block CRC => returning without decoding\n", 7, __FUNCTION__);
                        return;
                }
        }
        
        // decrypt
	if (block_header->flags & RED_LEVEL_1_ENCRYPTION_MASK) {
		rps->directives.encryption_level = LEVEL_1_ENCRYPTION;
		key = rps->password_data->level_1_encryption_key;
	} else if (block_header->flags & RED_LEVEL_2_ENCRYPTION_MASK) {
		rps->directives.encryption_level = LEVEL_2_ENCRYPTION;
		key = rps->password_data->level_2_encryption_key;
	} else
		rps->directives.encryption_level = NO_ENCRYPTION;
	if (rps->directives.encryption_level > NO_ENCRYPTION) {
		if (rps->password_data->access_level >= rps->directives.encryption_level) {
			AES_decrypt(block_header->statistics, block_header->statistics, NULL, key);
			block_header->flags &= ~RED_LEVEL_1_ENCRYPTION_MASK;
			block_header->flags &= ~RED_LEVEL_2_ENCRYPTION_MASK;
			rps->directives.encryption_level = -rps->directives.encryption_level;   // mark as decrypted
		} else {
			(void) fprintf(stderr, "%c\n%s(): No access to encrypted data => returning without decoding\n", 7, __FUNCTION__);
			return;
		}
	}
        
        // offset recording time
        if (MEF_globals->recording_time_offset_mode & (RTO_APPLY | RTO_APPLY_ON_INPUT))
                apply_recording_time_offset(&block_header->start_time);
        else if (MEF_globals->recording_time_offset_mode & (RTO_REMOVE | RTO_REMOVE_ON_INPUT))
                remove_recording_time_offset(&block_header->start_time);
	
	// discontinuity
	if (block_header->flags & RED_DISCONTINUITY_MASK)
		rps->directives.discontinuity = MEF_TRUE;
	else
		rps->directives.discontinuity = MEF_FALSE;
        
	// if no samples, just return
	if (block_header->number_of_samples == 0)
		return;
	
        // range decode difference data
        ui1_p = scaled_counts = block_header->statistics;
        *(ui4_p1 = cumulative_counts = rps->counts) = 0;
        ui4_p2 = ui4_p1 + 1;
	for (i = RED_BLOCK_STATISTICS_BYTES; i--;)
		*ui4_p2++ = *ui4_p1++ + (ui4) *ui1_p++;
	scaled_total_counts = *ui4_p1;
	
        diff_buffer_p = rps->difference_buffer;
        *diff_buffer_p++ = -128; // initial keysample flag not coded in encode (should be a low frequency symbol)
	
	ib_p = (ui1 *) rps->block_header + RED_BLOCK_HEADER_BYTES;
	in_byte = *ib_p++;
	low_bound = in_byte >> (8 - EXTRA_BITS);
	range = (ui4) 1 << EXTRA_BITS;
	ui4_p2 = cumulative_counts + 256;
	for (i = block_header->difference_bytes; i--;) {
                while (range <= BOTTOM_VALUE) {
			low_bound = (low_bound << 8) | ((in_byte << EXTRA_BITS) & 0xff);
                        in_byte = *ib_p++;
                        low_bound |= in_byte >> (8 - EXTRA_BITS);
                        range <<= 8;
                }
		temp_ui4 = low_bound / (range_per_count = range / scaled_total_counts);
		cc = (temp_ui4 >= scaled_total_counts ? (scaled_total_counts - 1) : temp_ui4);
		if (cc > cumulative_counts[128]) {
			for (ui4_p1 = ui4_p2; *--ui4_p1 > cc;);
			symbol = ui4_p1 - cumulative_counts;
		} else {
			for (ui4_p1 = cumulative_counts; *++ui4_p1 <= cc;);
			symbol = ui4_p1 - cumulative_counts - 1;
		}
		low_bound -= (temp_ui4 = range_per_count * cumulative_counts[symbol]);
		if (symbol < 255)
			range = range_per_count * scaled_counts[symbol];
		else
			range -= temp_ui4;
		*diff_buffer_p++ = symbol;
	}
	
	// generate output from difference data
	si1_p1 = (si1 *) rps->difference_buffer;
	si4_p = rps->decompressed_ptr;
	for (i = block_header->number_of_samples; i--;) {
		if (*si1_p1 == -128) {
			++si1_p1;
			si1_p2 = (si1 *) &current_val;
			*si1_p2++ = *si1_p1++; *si1_p2++ = *si1_p1++; *si1_p2++ = *si1_p1++; *si1_p2 = *si1_p1++;
		} else
			current_val += (si4) *si1_p1++;
		*si4_p++ = current_val;
	}
	
        // unscale decompressed_data if scaled (in place)
	if (block_header->scale_factor > (sf4) 1.0)
                RED_unscale(rps, rps->decompressed_ptr, rps->decompressed_ptr);
        
        // add trend to decompressed_data if detrended (in place)
	if ((block_header->detrend_slope != (sf4) 0.0) || (block_header->detrend_intercept != (sf4) 0.0))
                RED_retrend(rps, rps->decompressed_ptr, rps->decompressed_ptr);
	
	
        return;
}


si4	*RED_detrend(RED_PROCESSING_STRUCT *rps, si4 *input_buffer, si4 *output_buffer)
{
        si4			*si4_p1, *si4_p2;
        sf8			sx, sy, sxx, syy, sxy, n, mx, my, c, sf8_m, sf8_b, val;
        sf4			m, b;
        si8			i;
        RED_BLOCK_HEADER	*block_header;
	
	
        // detrend from input_buffer to output_buffer
	// slope and intercept values entered into block_header
	// if input_buffer == output_buffer detrending will be done in place
        
        block_header = rps->block_header;
	
        // calculate slope & intercept
        n = block_header->number_of_samples;
        sx = (n * (n + (sf8) 1.0)) / (sf8) 2.0;
        sxx = (n * (n + (sf8) 1.0) * ((n * (sf8) 2.0) + (sf8) 1.0)) / (sf8) 6.0;
        
        sy = syy = sxy = 0.0;
        c = (sf8) 1.0;
        si4_p1 = input_buffer;
        for (i = block_header->number_of_samples; i--;) {
                val = (sf8) *si4_p1++;
                sy += val;
                syy += val * val;
                sxy += val * c;
                c += (sf8) 1.0;
        }
        
        mx = sx / n;
        my = sy / n;
        sf8_m = (((sx * my) - sxy) / ((sx * mx) - sxx));
        sf8_b = (my - (sf8_m * mx));
        
        m = block_header->detrend_slope = (sf4) sf8_m;
        b = block_header->detrend_intercept = (sf4) sf8_b;
        
        // subtract trend from input_buffer to output_buffer
        sf8_m = (sf8) m;
        sf8_b = (sf8) b;
        c = (sf8) 0.0;
        si4_p1 = input_buffer;
	si4_p2 = output_buffer;
        for (i = block_header->number_of_samples; i--;) {
                c += (sf8) 1.0;
                *si4_p2++ = RED_round((sf8) *si4_p1++ - (sf8_m * c) - sf8_b);
        }
        
        
        return(output_buffer);
}


inline void	RED_encode(RED_PROCESSING_STRUCT *rps)
{
	// RED compress from original_ptr to block_header pointer (compressed data array)
	
	rps->block_header->detrend_slope = rps->block_header->detrend_intercept = (sf4) 0.0;
	if (rps->compression.mode == RED_LOSSLESS_COMPRESSION) {
		rps->block_header->scale_factor = (sf4) 1.0;
		RED_encode_exec(rps, rps->original_ptr, MEF_FALSE);
	} else if (rps->compression.mode > RED_LOSSLESS_COMPRESSION) {
		RED_encode_lossy(rps);
	} else {
		(void) fprintf(stderr, "%c\n%s(): Invalid compression mode => returning without encoding\n", 7, __FUNCTION__);
	}
	
	
	return;
}
	       
	       
void	RED_encode_exec(RED_PROCESSING_STRUCT *rps, si4 *input_buffer, si1 input_is_detrended)
{
	ui4			max_count, extra_bytes, range, r, underflow_bytes, low_bound, temp_ui4;
	ui4			*counts, *cumulative_counts, scaled_total_counts, *ui4_p1, *ui4_p2;
        ui1			*compressed_buffer_p, *diff_buffer_p, out_byte, *key, last_byte_val, *last_byte_ptr;
	sf8			stats_scale;
	si4			*si4_p1, *si4_p2, diff;
	si1			*si1_p1, *si1_p2;
	ui1			*ui1_p, *scaled_counts;
	si8			i;
	RED_BLOCK_HEADER	*block_header;
	
 
	// RED compress from input_buffer to to block_header pointer (compressed data array)
	block_header = rps->block_header;
	
        // apply recording time offset time
        if (MEF_globals->recording_time_offset_mode & (RTO_APPLY | RTO_APPLY_ON_OUTPUT))
                apply_recording_time_offset(&block_header->start_time);
        else if (MEF_globals->recording_time_offset_mode & (RTO_REMOVE | RTO_REMOVE_ON_OUTPUT))
                remove_recording_time_offset(&block_header->start_time);
	
	// if no samples: fill in an empty block header & return;
	if (block_header->number_of_samples <= 0) {
		block_header->flags = 0;
		block_header->detrend_slope = block_header->detrend_intercept = block_header->scale_factor = (sf4) 0.0;
		block_header->difference_bytes = block_header->number_of_samples = (ui4) 0;
		block_header->block_bytes = RED_BLOCK_HEADER_BYTES;
		bzero(block_header + RED_BLOCK_STATISTICS_OFFSET, RED_BLOCK_STATISTICS_BYTES);
		block_header->block_CRC = CRC_calculate((ui1 *) block_header + CRC_BYTES, block_header->block_bytes - CRC_BYTES);
		(void) fprintf(stderr, "%s(): No samples in block => returning without encoding\n", __FUNCTION__);
		return;
	}
        
        // detrend from input_buffer to detrended_buffer (lossless)
	if (rps->directives.detrend_data == MEF_TRUE && input_is_detrended == MEF_FALSE)
		input_buffer = RED_detrend(rps, input_buffer, rps->detrended_buffer);
	
        // scale if scale factor > 1 in block_header (lossy)
        if (block_header->scale_factor > (sf4) 1.0)
                input_buffer = RED_scale(rps, input_buffer, rps->scaled_buffer);
        
	// generate differences
        si4_p1 = input_buffer;
	si4_p2 = si4_p1 + 1;
	si1_p1 = rps->difference_buffer;  // first 4 bytes are keysample without keysample flag (-128)
	si1_p2 = (si1 *) si4_p1;
	*si1_p1++ = *si1_p2++; *si1_p1++ = *si1_p2++; *si1_p1++ = *si1_p2++; *si1_p1++ = *si1_p2;
	for (i = block_header->number_of_samples; --i;) {
		diff = *si4_p2++ - *si4_p1++;
		if (diff > 127 || diff < -127) {
			si1_p2 = (si1 *) si4_p1;
			*si1_p1++ = -128;
			*si1_p1++ = *si1_p2++; *si1_p1++ = *si1_p2++; *si1_p1++ = *si1_p2++; *si1_p1++ = *si1_p2;
		} else
			*si1_p1++ = (si1) diff;
	}
	block_header->difference_bytes = (ui4) (si1_p1 - rps->difference_buffer);
        
	// generate statistics
	counts = rps->counts;
	bzero(counts, RED_BLOCK_STATISTICS_BYTES * sizeof(ui4));
	ui1_p = (ui1 *) rps->difference_buffer;
	for (i = block_header->difference_bytes; i--;)
		++counts[*ui1_p++];
	ui4_p1 = counts;
	max_count = *ui4_p1;
	for (i = RED_BLOCK_STATISTICS_BYTES; --i;)
		if (*++ui4_p1 > max_count)
			max_count = *ui4_p1;
	ui1_p = scaled_counts = block_header->statistics;
	ui4_p1 = counts;
	if (max_count > 255) {
		stats_scale = (sf8) 254.999999999 / (sf8) max_count;
		for (i = RED_BLOCK_STATISTICS_BYTES; i--; ++ui4_p1) {
			if (*ui4_p1)
				*ui1_p++ = (ui1) ceil((sf8) *ui4_p1 * stats_scale);
			else
				*ui1_p++ = 0;
		}
	} else {
		for (i = RED_BLOCK_STATISTICS_BYTES; i--;)
			*ui1_p++ = (ui1) *ui4_p1++;
	}
	*(ui4_p1 = cumulative_counts = rps->counts) = 0;
	ui4_p2 = ui4_p1 + 1;
	ui1_p = scaled_counts;
	for (i = RED_BLOCK_STATISTICS_BYTES; i--;)
		*ui4_p2++ = (ui4) *ui1_p++ + *ui4_p1++;
	scaled_total_counts = *ui4_p1;
	
	// range encode
	low_bound = out_byte = underflow_bytes = 0;
	range = TOP_VALUE;
	last_byte_val = *(last_byte_ptr = compressed_buffer_p = ((ui1 *) rps->block_header + RED_BLOCK_HEADER_BYTES - 1)); // first byte is junk - save & overwrite last statistics byte
	diff_buffer_p = (ui1 *) rps->difference_buffer;
	for(i = block_header->difference_bytes; i--; ++diff_buffer_p) {
                while (range <= BOTTOM_VALUE) {
                        if (low_bound < CARRY_CHECK) {
                                *compressed_buffer_p++ = out_byte;
                                for (; underflow_bytes; --underflow_bytes)
                                        *compressed_buffer_p++ = 0xff;
                                out_byte = (ui1) (low_bound >> SHIFT_BITS);
                        } else if (low_bound & TOP_VALUE) {
                                *compressed_buffer_p++ = out_byte + 1;
                                for (; underflow_bytes; --underflow_bytes)
                                        *compressed_buffer_p++ = 0;
                                out_byte = (ui1) (low_bound >> SHIFT_BITS);
                        } else
                                ++underflow_bytes;
                        range <<= 8;
                        low_bound = (low_bound << 8) & TOP_VALUE_MINUS_1;
                }
                low_bound += (temp_ui4 = (r = range / scaled_total_counts) * cumulative_counts[*diff_buffer_p]);
                if (*diff_buffer_p < 0xff)
                        range = r * (ui4) scaled_counts[*diff_buffer_p];
                else
                        range -= temp_ui4;
        }
	while (range <= BOTTOM_VALUE) {
		if (low_bound < CARRY_CHECK) {
			*compressed_buffer_p++ = out_byte;
			for (; underflow_bytes; --underflow_bytes)
				*compressed_buffer_p++ = 0xff;
			out_byte = (ui1) (low_bound >> SHIFT_BITS);
		} else if (low_bound & TOP_VALUE) {
			*compressed_buffer_p++ = out_byte + 1;
			for (; underflow_bytes; --underflow_bytes)
				*compressed_buffer_p++ = 0;
			out_byte = (ui1) (low_bound >> SHIFT_BITS);
		} else
			++underflow_bytes;
		range <<= 8;
		low_bound = (low_bound << 8) & TOP_VALUE_MINUS_1;
	}
	temp_ui4 = (low_bound >> SHIFT_BITS) + 1;
	if (temp_ui4 > 0xff) {
		*compressed_buffer_p++ = out_byte + 1;
		for (; underflow_bytes; --underflow_bytes)
			*compressed_buffer_p++ = 0;
	} else {
		*compressed_buffer_p++ = out_byte;
		for (; underflow_bytes; --underflow_bytes)
			*compressed_buffer_p++ = 0xff;
	}
	*compressed_buffer_p++ = temp_ui4 & 0xff;
        *compressed_buffer_p++ = 0;
	
	// replace last statistics byte
	*last_byte_ptr = last_byte_val;
        
        // add one for uncoded initial keysample flag - need that byte allocated in decode
        ++block_header->difference_bytes;
        
        // align to 8-byte boundary
	block_header->block_bytes = (ui4) (compressed_buffer_p - (ui1 *) block_header);
	extra_bytes = block_header->block_bytes % 8;
	if (extra_bytes) {
		extra_bytes = 8 - extra_bytes;
		for (i = extra_bytes; i--;)
			*compressed_buffer_p++ = PAD_BYTE_VALUE;
                block_header->block_bytes += extra_bytes;
	}

        // flags
        block_header->flags = 0;
        
        // encrypt
        if (rps->directives.encryption_level > NO_ENCRYPTION) {
                if (rps->password_data->access_level >= rps->directives.encryption_level) {
                        if (rps->directives.encryption_level == LEVEL_1_ENCRYPTION) {
                                key = rps->password_data->level_1_encryption_key;
				block_header->flags |= RED_LEVEL_1_ENCRYPTION_MASK;
				block_header->flags &= ~RED_LEVEL_2_ENCRYPTION_MASK;
			} else {
                                key = rps->password_data->level_2_encryption_key;
				block_header->flags &= ~RED_LEVEL_1_ENCRYPTION_MASK;
				block_header->flags |= RED_LEVEL_2_ENCRYPTION_MASK;
			}
                        AES_encrypt(block_header->statistics, block_header->statistics, NULL, key);
                } else {
                        (void) fprintf(stderr, "%c\n%s(): Cannot encrypt data => returning without encrypting\n", 7, __FUNCTION__);
                }
	} else {
		block_header->flags &= ~RED_LEVEL_1_ENCRYPTION_MASK;
		block_header->flags &= ~RED_LEVEL_2_ENCRYPTION_MASK;
		rps->directives.encryption_level = NO_ENCRYPTION;
	}
        
        // discontinuity
	if (rps->directives.discontinuity == MEF_TRUE) {
                block_header->flags |= RED_DISCONTINUITY_MASK;
		if (rps->directives.reset_discontinuity == MEF_TRUE)
			rps->directives.discontinuity = MEF_FALSE;
	}
	
        // calculate CRC
        block_header->block_CRC = CRC_calculate((ui1 *) block_header + CRC_BYTES, block_header->block_bytes - CRC_BYTES);
	
	
	return;
}


void 	RED_encode_lossy(RED_PROCESSING_STRUCT *rps)
{
	si1			input_is_detrended, compression_mode;
	si8			i;
        si4			*input_buffer;
	sf8			original_size, goal_compression_ratio, r;
	sf8			low_sf, high_sf, mrr, mrr2, mrr5, sf_per_mrr;
	sf8			goal_low_bound, goal_high_bound, goal_mrr, goal_tol;
	sf4			new_scale_factor;
	RED_BLOCK_HEADER	*block_header;
	
	mrr = 0;


	// RED compress from original_ptr to block_header pointer (compressed data array)
	input_buffer = rps->original_ptr;
	block_header = rps->block_header;
	
	// detrend data from original_ptr to detrended_buffer
	if (rps->directives.detrend_data == MEF_TRUE) {
		input_buffer = RED_detrend(rps, input_buffer, rps->detrended_buffer);
		input_is_detrended = MEF_TRUE;
	} else {
		input_is_detrended = MEF_FALSE;
	}
        
        // check block normality
        compression_mode = rps->compression.mode;
        if (rps->directives.require_normality == MEF_TRUE) {
                r = RED_test_normality(input_buffer, block_header->number_of_samples);
		if (r < rps->directives.normal_correlation)
                        compression_mode = RED_LOSSLESS_COMPRESSION;
        }
	
	// RED encode
	switch(compression_mode) {
		case RED_LOSSLESS_COMPRESSION:
			block_header->scale_factor = (sf4) 1.0; // fall through to RED_FIXED_SCALE_FACTOR
		case RED_FIXED_SCALE_FACTOR:
			RED_encode_exec(rps, input_buffer, input_is_detrended);
			break;
		case RED_FIXED_COMPRESSION_RATIO:
			goal_compression_ratio = rps->compression.goal_compression_ratio;
			goal_low_bound = goal_compression_ratio - rps->compression.goal_tolerance;
			goal_high_bound = goal_compression_ratio + rps->compression.goal_tolerance;
			block_header->scale_factor = (sf4) 1.0;
			RED_encode_exec(rps, input_buffer, input_is_detrended);
			original_size = (sf8) block_header->number_of_samples * (sf8) sizeof(si4);
			rps->compression.actual_compression_ratio = (sf8) block_header->block_bytes / original_size;
			if (rps->compression.actual_compression_ratio <= goal_high_bound)
				break;
			// loop until acceptable scale factor found
			for (i = rps->compression.maximum_rounds_per_block; i--;) {
				new_scale_factor = block_header->scale_factor * (sf4) (rps->compression.actual_compression_ratio / goal_compression_ratio);
				if ((ABS(new_scale_factor - block_header->scale_factor) <= (sf4) 0.000001) || (new_scale_factor <= (sf4) 1.0))
					break;
				block_header->scale_factor = new_scale_factor;
				RED_encode_exec(rps, input_buffer, input_is_detrended);
				rps->compression.actual_compression_ratio = (sf8) block_header->block_bytes / original_size;
				if ((rps->compression.actual_compression_ratio <= goal_high_bound) && (rps->compression.actual_compression_ratio >= goal_low_bound))
					break;
			}
			break;
		case RED_MEAN_RESIDUAL_RATIO:
                        // get residual ratio at sf 2 & 5 (roughly linear relationship: reasonable sample points)
                        block_header->scale_factor = (sf4) 2.0;
                        RED_generate_lossy_data(rps, input_buffer, rps->decompressed_ptr, input_is_detrended);
                        mrr2 = RED_calculate_mean_residual_ratio(rps->original_ptr, rps->decompressed_ptr, block_header->number_of_samples);
                        if (mrr2 == (sf8) 0.0) {  // all zeros in block
                                block_header->scale_factor = (sf4) 1.0;
                                rps->compression.actual_mean_residual_ratio = (sf8) 0.0;
                                RED_encode_exec(rps, input_buffer, input_is_detrended);
                                break;
                        }
                        block_header->scale_factor = (sf4) 5.0;
                        RED_generate_lossy_data(rps, input_buffer, rps->decompressed_ptr, input_is_detrended);
                        mrr5 = RED_calculate_mean_residual_ratio(rps->original_ptr, rps->decompressed_ptr, block_header->number_of_samples);
                        sf_per_mrr = (sf8) 3.0 / (mrr5 - mrr2);
                        // estimate starting points
                        goal_mrr = rps->compression.goal_mean_residual_ratio;
                        goal_tol = rps->compression.goal_tolerance;
                        goal_low_bound = goal_mrr - goal_tol;
                        goal_high_bound = goal_mrr + goal_tol;
                        block_header->scale_factor = (sf4) (((goal_mrr - mrr2) * sf_per_mrr) + (sf8) 2.0);
                        high_sf = ((goal_high_bound - mrr2) * sf_per_mrr) + (sf8) 2.0;
                        high_sf *= (sf8) 2.0;  // empirically reasonable
                        low_sf = (sf8) 1.0;
			for (i = rps->compression.maximum_rounds_per_block; i--;) {
				RED_generate_lossy_data(rps, input_buffer, rps->decompressed_ptr, input_is_detrended);
				mrr = RED_calculate_mean_residual_ratio(rps->original_ptr, rps->decompressed_ptr, block_header->number_of_samples);
				if (mrr < goal_low_bound)
					low_sf = (sf8) block_header->scale_factor;
				else if (mrr > goal_high_bound)
					high_sf = (sf8) block_header->scale_factor;
                                else
                                        break;
                                new_scale_factor = (sf4) ((low_sf + high_sf) / (sf8) 2.0);
                                if (new_scale_factor <= (sf4) 1.0)
                                        break;
                                if ((high_sf - low_sf) < (sf8) 0.005) {
                                        block_header->scale_factor = new_scale_factor;
                                        break;
                                }
                                block_header->scale_factor = new_scale_factor;
			}
			rps->compression.actual_mean_residual_ratio = mrr;
			RED_encode_exec(rps, input_buffer, input_is_detrended);
			break;
		default:
			(void) fprintf(stderr, "%c\n%s(): unrecognized compression mode (%d) => exiting\n", 7, __FUNCTION__, rps->compression.mode);
			exit(1);
			break;
	}
	
	// generate lossy data from input_buffer to decompressed_ptr
	if (rps->directives.return_lossy_data == MEF_TRUE)
		RED_generate_lossy_data(rps, input_buffer, rps->decompressed_ptr, input_is_detrended);
	
	
	return;
}


void	RED_filter(FILT_PROCESSING_STRUCT *filtps)
{
	si4	*si4_p;
	si8	i;
	sf8	*sf8_p;
	
	
	// filter
	FILT_filtfilt(filtps);
	
	// make si4 version of filtered data
	si4_p = filtps->filt_data;
	sf8_p = filtps->sf8_filt_data;
	for (i = filtps->data_length; i--;)
		*si4_p++ = RED_round(*sf8_p++);
	
	
	return;
}


void	RED_find_extrema(si4 *buffer, si8 number_of_samples, TIME_SERIES_INDEX *tsi)
{
	si4	min, max;
	si8	i;
	
	
	min = max = *buffer;
	for (i = number_of_samples - 1; i--;) {
		if (*++buffer > max)
			max = *buffer;
		else if (*buffer < min)
			min = *buffer;
	}
	
	tsi->maximum_sample_value = max;
	tsi->minimum_sample_value = min;

	
	return;
}


void	RED_free_processing_struct(RED_PROCESSING_STRUCT *rps)
{
	if (rps->original_data != NULL)
		free(rps->original_data);
	
	if (rps->decompressed_data != NULL)
		free(rps->decompressed_data);
	
	if (rps->compressed_data != NULL)
		free(rps->compressed_data);
	
	if (rps->difference_buffer != NULL)
		free(rps->difference_buffer);
	
	if (rps->detrended_buffer != NULL)
		free(rps->detrended_buffer);
	
	if (rps->scaled_buffer != NULL)
		free(rps->scaled_buffer);
	
	free(rps);
	
        
	return;
}


void	RED_generate_lossy_data(RED_PROCESSING_STRUCT *rps, si4 *input_buffer, si4 *output_buffer, si1 input_is_detrended)
{	
	
        // generates lossy data from input_buffer to output_buffer
	// if input_buffer == output_buffer lossy data will be made in place
        	
	// detrend from input_buffer to output_buffer (lossless)
	if (rps->directives.detrend_data == MEF_TRUE && input_is_detrended == MEF_FALSE) {
		input_buffer = RED_detrend(rps, input_buffer, output_buffer);
		input_is_detrended = MEF_TRUE;
	}

        // scale from input_buffer to output_buffer (lossy)
	input_buffer = RED_scale(rps, input_buffer, output_buffer);
                
	// unscale from input_buffer to output_buffer
	input_buffer = RED_unscale(rps, input_buffer, output_buffer);
	
        // retrend from input_buffer to output_buffer
	if (input_is_detrended == MEF_TRUE)
		RED_retrend(rps, input_buffer, output_buffer);

	
        return;
}

	       
sf8	*RED_initialize_normal_CDF_table(si4 global_flag)
{
	sf8	*cdf_table;
	
	
	cdf_table = (sf8 *) e_calloc((size_t) RED_NORMAL_CDF_TABLE_ENTRIES, sizeof(sf8), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	{
		sf8 temp[RED_NORMAL_CDF_TABLE_ENTRIES] = RED_NORMAL_CDF_TABLE;
		memcpy(cdf_table, temp, RED_NORMAL_CDF_TABLE_ENTRIES * sizeof(sf8));
	}
	
	if (global_flag == MEF_TRUE) {
		MEF_globals->RED_normal_CDF_table = cdf_table;
		return(NULL);
	}
	
	
	return(cdf_table);
}


si4	*RED_retrend(RED_PROCESSING_STRUCT *rps, si4 *input_buffer, si4 *output_buffer)
{
	si4			*si4_p1, *si4_p2;
	si8			i;
	sf8			m, b, c;
	RED_BLOCK_HEADER	*block_header;
	
	
	// retrend data from input_buffer to output_buffer
	// if input_buffer == output_buffer retrending data will be done in place
	
	block_header = rps->block_header;

	m = (sf8) block_header->detrend_slope;
	b = (sf8) (sf8) block_header->detrend_intercept;
	c = (sf8) 0.0;
	si4_p1 = input_buffer;
	si4_p2 = output_buffer;
	for (i = block_header->number_of_samples; i--;) {
                c += (sf8) 1.0;
		*si4_p2++ = RED_round( (sf8) *si4_p1++ + (m * c) + b);
	}
	

	return(output_buffer);
}


inline si4 RED_round(sf8 val)
{
	if (val >= 0.0) {
		if ((val += 0.5) >= (sf8) RED_POSITIVE_INFINITY)
			return(RED_POSITIVE_INFINITY);
	} else {
		if ((val -= 0.5) <= (sf8) RED_NEGATIVE_INFINITY)
			return(RED_NEGATIVE_INFINITY);
	}
	
	
	return((si4) val);
}


si4	*RED_scale(RED_PROCESSING_STRUCT *rps, si4 *input_buffer, si4 *output_buffer)
{
        si4	*si4_p1, *si4_p2;
        si8	i;
        sf8	sf;
	
	
	// scale from input_buffer to output_buffer
	// if input_buffer == output_buffer scaling will be done in place
	
        si4_p1 = input_buffer;
        si4_p2 = output_buffer;
        sf = (sf8) rps->block_header->scale_factor;
        for (i = rps->block_header->number_of_samples; i--;)
                *si4_p2++ = RED_round((sf8) *si4_p1++ / sf);
	
        
        return(output_buffer);
}


void	RED_show_block_header(RED_BLOCK_HEADER *bh)
{
	si1			hex_str[HEX_STRING_BYTES(RED_BLOCK_STATISTICS_BYTES)], time_str[TIME_STRING_BYTES];
	
	
	printf("--------------- RED Block Header - START ---------------\n");
	if (bh->block_CRC == CRC_NO_ENTRY)
		printf("Block CRC: no entry\n");
	else {
		generate_hex_string((ui1 *) &bh->block_CRC, CRC_BYTES, hex_str);
		printf("Block CRC: %s\n", hex_str);
	}
	generate_hex_string((ui1 *) &bh->flags, 1, hex_str);
	printf("Flags: %s\n", hex_str);
	printf("Detrend Slope: %f\n", bh->detrend_slope);
	printf("Detrend Intercept: %f\n", bh->detrend_intercept);
	printf("Scale Factor: %f\n", bh->scale_factor);
	printf("Difference Bytes: %u\n", bh->difference_bytes);
	printf("Number of Samples: %u\n", bh->number_of_samples);
	printf("Block Bytes: %u\n", bh->block_bytes);
	if (bh->start_time == UUTC_NO_ENTRY)
		printf("Start Time: no entry\n");
	else {
		local_date_time_string(bh->start_time, time_str);
		if (bh->start_time < 0)
			printf("Offset ");
		#ifdef _WIN32
			printf("Start Time: %lld (uUTC), %s (ascii, local)\n", ABS(bh->start_time), time_str);
		#else
			printf("Start Time: %ld (uUTC), %s (ascii, local)\n", ABS(bh->start_time), time_str);
		#endif
	}
	generate_hex_string((ui1 *) &bh->statistics, RED_BLOCK_STATISTICS_BYTES, hex_str);
	printf("Statistics: %s\n", hex_str);
	printf("---------------- RED Block Header - END ----------------\n\n");

	
	return;
}


sf8	RED_test_normality(si4 *data, ui4 n_samps)
{
	sf8	sx, sx2, sy, sy2, sxy, mx, mx2, sd, val, z, r, n, *norm_cdf;
	sf8	num, den1, den2, cdf[RED_NORMAL_CDF_TABLE_ENTRIES];
	si8	i, counts[RED_NORMAL_CDF_TABLE_ENTRIES] = {0};
        si4	*si4_p, bin;
        
        
        // returns the correlation of the distribution in the data to that expected from a normal distribution
        
        // calculate mean & standard deviation
	n = (sf8) n_samps;
        si4_p = data;
        sx = sx2 = (sf8) 0.0;
        for (i = n_samps; i--;) {
                val = (sf8) *si4_p++;
                sx += val;
                sx2 += val * val;
        }
        mx = sx / n;
        mx2 = sx2 / n;
        sd = sqrt(mx2 - (mx * mx));
	
	// bin the samples
	si4_p = data;
	for (i = n_samps; i--;) {
		val = (sf8) *si4_p++;
		z = (val - mx) / sd;
		bin = (si4) ((z + (sf8) 3.1) * (sf8) 10.0);
		if (bin < 0)
			bin = 0;
		if (bin >= RED_NORMAL_CDF_TABLE_ENTRIES)
			continue;
		++counts[bin];
	}
	
	// generate data CDF
	cdf[0] = (sf8) counts[0] / (sf8) n_samps;
	for (i = 1; i < RED_NORMAL_CDF_TABLE_ENTRIES; ++i)
		cdf[i] = (sf8) (counts[i] + counts[i - 1]) / n;
		
        // calculate correlation between data CDF and normal CDF
	sx = sx2 = sxy = (sf8) 0.0;
	sy = (sf8) RED_SUM_NORMAL_CDF;
	sy2 = (sf8) RED_SUM_SQ_NORMAL_CDF;
	norm_cdf = MEF_globals->RED_normal_CDF_table;
	for (i = 0; i < RED_NORMAL_CDF_TABLE_ENTRIES; ++i) {
		sx += cdf[i];
		sx2 += cdf[i] * cdf[i];
		sxy += cdf[i] * norm_cdf[i];
	}
	
	num = (n * sxy) - (sx * sy);
	den1 = (n * sx2) - (sx * sx);
	den2 = (n * sy2) - (sy * sy);
	
	if (den1 <= (sf8) 0.0)
		r = (sf8) 0.0;
	else
		r = num / (sqrt(den1) * sqrt(den2));
	
        return(r);
}


si4	*RED_unscale(RED_PROCESSING_STRUCT *rps, si4 *input_buffer, si4 *output_buffer)
{
	si4			*si4_p1, *si4_p2;
	si8			i;
	sf8			sf;
	
	
	// unscale from input_buffer to output_buffer
	// if input_buffer == output_buffer unscaling will be done in place
	
	si4_p1 = input_buffer;
	si4_p2 = output_buffer;
	sf = (sf8) rps->block_header->scale_factor;
	for (i = rps->block_header->number_of_samples; i--;)
		*si4_p2++ = RED_round((sf8) *si4_p1++ * sf);
	
	
	return(output_buffer);
}


inline RED_BLOCK_HEADER	*RED_update_RPS_pointers(RED_PROCESSING_STRUCT *rps, ui1 flags)
{
        RED_BLOCK_HEADER	*block_header;

        
        block_header = rps->block_header;
	if (flags & RED_UPDATE_ORIGINAL_PTR)
                rps->original_ptr += block_header->number_of_samples;
	if (flags & RED_UPDATE_BLOCK_HEADER_PTR)
                rps->block_header = (RED_BLOCK_HEADER *) ((ui1 *) rps->block_header + block_header->block_bytes);
        if (flags & RED_UPDATE_DECOMPRESSED_PTR)
        	rps->decompressed_ptr += block_header->number_of_samples;
        
        
        return(rps->block_header);
}


/*************************************************************************/
/****************************  END RED FUNCTIONS  ************************/
/*************************************************************************/


si4	remove_line_noise(si4 *data, si8 n_samps, sf8 sampling_frequency, sf8 line_frequency, sf8 *template)
{
        FILT_PROCESSING_STRUCT	*filtps;
        si8			i, j, si8_curr_samp, leftovers, median_pt;
        sf8			*filt_data, sf8_curr_samp, *point_arrays, *pa, sf8_template_len;
        si4			template_len, n_waveforms;
        si1			free_template;
        
        
        // filter
	filtps = FILT_initialize_processing_struct(5, FILT_BANDPASS_TYPE, sampling_frequency, n_samps, MEF_FALSE, MEF_FALSE, line_frequency - 10.0, (line_frequency * 5.0) + 10.0);
        filtps->orig_data = data;
        filtps->data_length = n_samps;
        FILT_filtfilt(filtps);
        
        // allocate
        filt_data = filtps->sf8_filt_data;
        sf8_template_len = sampling_frequency / line_frequency;
        template_len = (si4) (sf8_template_len + 0.5);
        n_waveforms = (si4) ((sf8) n_samps / sf8_template_len);
        point_arrays = (sf8 *) e_calloc((size_t) template_len * n_waveforms, sizeof(sf8), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	if (template == NULL) {
                template = (sf8 *) e_calloc((size_t) template_len, sizeof(sf8), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_template = MEF_TRUE;
	} else {
		free_template = MEF_FALSE;
	}
	
        // reorder points
        sf8_curr_samp = 0.0;
        for  (i = 0; i < n_waveforms; ++i) {
                si8_curr_samp = (si8) (sf8_curr_samp + 0.5);
                pa = point_arrays + i;
                for (j = 0; j < template_len; ++j) {
                        *pa = filt_data[si8_curr_samp++];
                        pa += n_waveforms;
                }
                sf8_curr_samp += sf8_template_len;
        }
        
        // sort
        pa = point_arrays;
        for (i = 0; i < template_len; ++i) {
                qsort(pa, n_waveforms, sizeof(sf8), compare_sf8);
                pa += n_waveforms;
        }
        
        // build template from medians
        median_pt = n_waveforms / 2;
        pa = point_arrays + median_pt;
        for (i = 0; i < template_len; ++i) {
                template[i] = *pa;
                pa += n_waveforms;
        }
	
        // subtract template
        sf8_curr_samp = 0.0;
	si8_curr_samp = 0;
        for  (i = 0; i < n_waveforms; ++i) {
                for (j = 0; j < template_len; ++j, ++si8_curr_samp)
                        data[si8_curr_samp] = RED_round((sf8) data[si8_curr_samp] - template[j]);
                sf8_curr_samp += sf8_template_len;
		if (si8_curr_samp < (si8) (sf8_curr_samp + 0.5)) {
                        data[si8_curr_samp] -= template[0];
			++si8_curr_samp;
		}
        }
        leftovers = n_samps - ++si8_curr_samp;
        for (i = 0; i < leftovers; ++i, si8_curr_samp++)
                data[si8_curr_samp] = RED_round((sf8) data[si8_curr_samp] - template[i]);
	
        // clean up
        free(point_arrays);
        FILT_free_processing_struct(filtps, MEF_FALSE, MEF_FALSE);
        if (free_template == MEF_TRUE)
                free(template);
        
        
        return(template_len);
}


void	remove_line_noise_adaptive(si4 *data, si8 n_samps, sf8 sampling_frequency, sf8 line_frequency, si4 n_cycles)
{
	FILT_PROCESSING_STRUCT	*filtps;
	si8			i, j, si8_curr_samp, leftovers;
	sf8			*filt_data, sf8_curr_samp, *point_arrays, *pa, *ma, *tma, sf8_template_len;
	si4			template_len, n_waveforms;
	
	
	// filter
	filtps = FILT_initialize_processing_struct(5, FILT_BANDPASS_TYPE, sampling_frequency, n_samps, MEF_FALSE, MEF_FALSE, line_frequency - 10.0, (line_frequency * 5.0) + 10.0);
	filtps->orig_data = data;
	filtps->data_length = n_samps;
	FILT_filtfilt(filtps);
	
	// allocate
	filt_data = filtps->sf8_filt_data;
	sf8_template_len = sampling_frequency / line_frequency;
	template_len = (si4) (sf8_template_len + 0.5);
	n_waveforms = (si4) ((sf8) n_samps / sf8_template_len);
	point_arrays = (sf8 *) e_calloc((size_t) template_len * n_waveforms, sizeof(sf8), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	// reorder points
	sf8_curr_samp = 0.0;
	for  (i = 0; i < n_waveforms; ++i) {
		si8_curr_samp = (si8) (sf8_curr_samp + 0.5);
		pa = point_arrays + i;
		for (j = 0; j < template_len; ++j) {
			*pa = filt_data[si8_curr_samp++];
			pa += n_waveforms;
		}
		sf8_curr_samp += sf8_template_len;
	}
	
	// sort
	pa = point_arrays;
	ma = filt_data;
	for (i = 0; i < template_len; ++i) {
		proportion_filt(pa, ma, n_waveforms, 0.5, n_cycles);
		pa += n_waveforms;
		ma += n_waveforms;
	}
	
	// subtract template
	sf8_curr_samp = 0.0;
	si8_curr_samp = 0;
	ma = filt_data - 1;
	for  (i = 0; i < n_waveforms; ++i, ++ma) {
		tma = ma - n_waveforms;
		for (j = 0; j < template_len; ++j, ++si8_curr_samp)
			data[si8_curr_samp] = RED_round((sf8) data[si8_curr_samp] - *(tma += n_waveforms));
		sf8_curr_samp += sf8_template_len;
		if (si8_curr_samp < (si8) (sf8_curr_samp + 0.5)) {
			data[si8_curr_samp] -= *ma;
			++si8_curr_samp;
		}
	}
	leftovers = n_samps - ++si8_curr_samp;
	tma = ma - n_waveforms;
	for (i = 0; i < leftovers; ++i, si8_curr_samp++)
		data[si8_curr_samp] = RED_round((sf8) data[si8_curr_samp] - *(tma += n_waveforms));
	
	// clean up
	free(point_arrays);
	FILT_free_processing_struct(filtps, MEF_FALSE, MEF_FALSE);
	
	
	return;
}


inline void	remove_recording_time_offset(si8 *time)
{
        if (*time == UUTC_NO_ENTRY)
                return;
	
	if (*time > 0)  // positive times indicate recording time offset not applied, 0 - unlikely any records started at 1970.
		return;
	
	// remove recording time offset & make positive to indicate removal
	*time = (-*time) + MEF_globals->recording_time_offset;
	
	
	return;
}


/*************************************************************************/
/****************************  SHA-256 FUNCTIONS  ************************/
/*************************************************************************/

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


void sha256(const ui1 *message, ui4 len, ui1 *digest)
{
	SHA256_ctx	ctx;
	
	SHA256_init(&ctx);
	SHA256_update(&ctx, message, len);
	SHA256_final(&ctx, digest);
	
        
	return;
}


void SHA256_final(SHA256_ctx *ctx, ui1 *digest)
{
	ui4	block_nb;
	ui4	pm_len;
	ui4	len_b;
	
	block_nb = (1 + ((SHA256_BLOCK_SIZE - 9) < (ctx->len % SHA256_BLOCK_SIZE)));
	
	len_b = (ctx->tot_len + ctx->len) << 3;
	pm_len = block_nb << 6;
	
	memset(ctx->block + ctx->len, 0, pm_len - ctx->len);
	ctx->block[ctx->len] = 0x80;
	UNPACK32(len_b, ctx->block + pm_len - 4);
	
	SHA256_transf(ctx, ctx->block, block_nb);
	
	UNPACK32(ctx->h[0], &digest[ 0]);
	UNPACK32(ctx->h[1], &digest[ 4]);
	UNPACK32(ctx->h[2], &digest[ 8]);
	UNPACK32(ctx->h[3], &digest[12]);
	UNPACK32(ctx->h[4], &digest[16]);
	UNPACK32(ctx->h[5], &digest[20]);
	UNPACK32(ctx->h[6], &digest[24]);
	UNPACK32(ctx->h[7], &digest[28]);
	
        
	return;
}


void SHA256_init(SHA256_ctx *ctx)
{
	if (MEF_globals->SHA256_h0_table == NULL)
		(void) SHA256_initialize_h0_table(MEF_TRUE);
	
	ctx->h[0] = MEF_globals->SHA256_h0_table[0]; ctx->h[1] = MEF_globals->SHA256_h0_table[1];
	ctx->h[2] = MEF_globals->SHA256_h0_table[2]; ctx->h[3] = MEF_globals->SHA256_h0_table[3];
	ctx->h[4] = MEF_globals->SHA256_h0_table[4]; ctx->h[5] = MEF_globals->SHA256_h0_table[5];
	ctx->h[6] = MEF_globals->SHA256_h0_table[6]; ctx->h[7] = MEF_globals->SHA256_h0_table[7];
	
	ctx->len = 0;
	ctx->tot_len = 0;
	
        
	return;
}


ui4	*SHA256_initialize_h0_table(si4 global_flag)
{
	ui4	*SHA256_h0_table;
	
	
	SHA256_h0_table = (ui4 *) e_calloc((size_t) SHA256_H0_ENTRIES, sizeof(ui4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	{
		ui4 temp[SHA256_H0_ENTRIES] = SHA256_H0;
		memcpy(SHA256_h0_table, temp, SHA256_H0_ENTRIES * sizeof(ui4));
	}
	
	if (global_flag == MEF_TRUE) {
		MEF_globals->SHA256_h0_table = SHA256_h0_table;
		return(NULL);
	}
	
        
	return(SHA256_h0_table);
}


ui4	*SHA256_initialize_k_table(si4 global_flag)
{
	ui4	*SHA256_k_table;
	
	
	SHA256_k_table = (ui4 *) e_calloc((size_t) SHA256_K_ENTRIES, sizeof(ui4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	{
		ui4 temp[SHA256_K_ENTRIES] = SHA256_K;
		memcpy(SHA256_k_table, temp, SHA256_K_ENTRIES * sizeof(ui4));
	}
	
	if (global_flag == MEF_TRUE) {
		MEF_globals->SHA256_k_table = SHA256_k_table;
		return(NULL);
	}
	
        
	return(SHA256_k_table);
}


void SHA256_transf(SHA256_ctx *ctx, const ui1 *message, ui4 block_nb)
{
	ui4		w[64];
	ui4		wv[8];
	ui4		t1, t2;
	const ui1	*sub_block;
	si4		i;
        
	
	if (MEF_globals->SHA256_k_table == NULL)
		(void) SHA256_initialize_k_table(MEF_TRUE);
	
	for (i = 0; i < (si4) block_nb; i++) {
		sub_block = message + (i << 6);
		
		PACK32(&sub_block[ 0], &w[ 0]); PACK32(&sub_block[ 4], &w[ 1]);
		PACK32(&sub_block[ 8], &w[ 2]); PACK32(&sub_block[12], &w[ 3]);
		PACK32(&sub_block[16], &w[ 4]); PACK32(&sub_block[20], &w[ 5]);
		PACK32(&sub_block[24], &w[ 6]); PACK32(&sub_block[28], &w[ 7]);
		PACK32(&sub_block[32], &w[ 8]); PACK32(&sub_block[36], &w[ 9]);
		PACK32(&sub_block[40], &w[10]); PACK32(&sub_block[44], &w[11]);
		PACK32(&sub_block[48], &w[12]); PACK32(&sub_block[52], &w[13]);
		PACK32(&sub_block[56], &w[14]); PACK32(&sub_block[60], &w[15]);
		
		SHA256_SCR(16); SHA256_SCR(17); SHA256_SCR(18); SHA256_SCR(19);
		SHA256_SCR(20); SHA256_SCR(21); SHA256_SCR(22); SHA256_SCR(23);
		SHA256_SCR(24); SHA256_SCR(25); SHA256_SCR(26); SHA256_SCR(27);
		SHA256_SCR(28); SHA256_SCR(29); SHA256_SCR(30); SHA256_SCR(31);
		SHA256_SCR(32); SHA256_SCR(33); SHA256_SCR(34); SHA256_SCR(35);
		SHA256_SCR(36); SHA256_SCR(37); SHA256_SCR(38); SHA256_SCR(39);
		SHA256_SCR(40); SHA256_SCR(41); SHA256_SCR(42); SHA256_SCR(43);
		SHA256_SCR(44); SHA256_SCR(45); SHA256_SCR(46); SHA256_SCR(47);
		SHA256_SCR(48); SHA256_SCR(49); SHA256_SCR(50); SHA256_SCR(51);
		SHA256_SCR(52); SHA256_SCR(53); SHA256_SCR(54); SHA256_SCR(55);
		SHA256_SCR(56); SHA256_SCR(57); SHA256_SCR(58); SHA256_SCR(59);
		SHA256_SCR(60); SHA256_SCR(61); SHA256_SCR(62); SHA256_SCR(63);
		
		wv[0] = ctx->h[0]; wv[1] = ctx->h[1];
		wv[2] = ctx->h[2]; wv[3] = ctx->h[3];
		wv[4] = ctx->h[4]; wv[5] = ctx->h[5];
		wv[6] = ctx->h[6]; wv[7] = ctx->h[7];
		
		SHA256_EXP(0,1,2,3,4,5,6,7, 0); SHA256_EXP(7,0,1,2,3,4,5,6, 1);
		SHA256_EXP(6,7,0,1,2,3,4,5, 2); SHA256_EXP(5,6,7,0,1,2,3,4, 3);
		SHA256_EXP(4,5,6,7,0,1,2,3, 4); SHA256_EXP(3,4,5,6,7,0,1,2, 5);
		SHA256_EXP(2,3,4,5,6,7,0,1, 6); SHA256_EXP(1,2,3,4,5,6,7,0, 7);
		SHA256_EXP(0,1,2,3,4,5,6,7, 8); SHA256_EXP(7,0,1,2,3,4,5,6, 9);
		SHA256_EXP(6,7,0,1,2,3,4,5,10); SHA256_EXP(5,6,7,0,1,2,3,4,11);
		SHA256_EXP(4,5,6,7,0,1,2,3,12); SHA256_EXP(3,4,5,6,7,0,1,2,13);
		SHA256_EXP(2,3,4,5,6,7,0,1,14); SHA256_EXP(1,2,3,4,5,6,7,0,15);
		SHA256_EXP(0,1,2,3,4,5,6,7,16); SHA256_EXP(7,0,1,2,3,4,5,6,17);
		SHA256_EXP(6,7,0,1,2,3,4,5,18); SHA256_EXP(5,6,7,0,1,2,3,4,19);
		SHA256_EXP(4,5,6,7,0,1,2,3,20); SHA256_EXP(3,4,5,6,7,0,1,2,21);
		SHA256_EXP(2,3,4,5,6,7,0,1,22); SHA256_EXP(1,2,3,4,5,6,7,0,23);
		SHA256_EXP(0,1,2,3,4,5,6,7,24); SHA256_EXP(7,0,1,2,3,4,5,6,25);
		SHA256_EXP(6,7,0,1,2,3,4,5,26); SHA256_EXP(5,6,7,0,1,2,3,4,27);
		SHA256_EXP(4,5,6,7,0,1,2,3,28); SHA256_EXP(3,4,5,6,7,0,1,2,29);
		SHA256_EXP(2,3,4,5,6,7,0,1,30); SHA256_EXP(1,2,3,4,5,6,7,0,31);
		SHA256_EXP(0,1,2,3,4,5,6,7,32); SHA256_EXP(7,0,1,2,3,4,5,6,33);
		SHA256_EXP(6,7,0,1,2,3,4,5,34); SHA256_EXP(5,6,7,0,1,2,3,4,35);
		SHA256_EXP(4,5,6,7,0,1,2,3,36); SHA256_EXP(3,4,5,6,7,0,1,2,37);
		SHA256_EXP(2,3,4,5,6,7,0,1,38); SHA256_EXP(1,2,3,4,5,6,7,0,39);
		SHA256_EXP(0,1,2,3,4,5,6,7,40); SHA256_EXP(7,0,1,2,3,4,5,6,41);
		SHA256_EXP(6,7,0,1,2,3,4,5,42); SHA256_EXP(5,6,7,0,1,2,3,4,43);
		SHA256_EXP(4,5,6,7,0,1,2,3,44); SHA256_EXP(3,4,5,6,7,0,1,2,45);
		SHA256_EXP(2,3,4,5,6,7,0,1,46); SHA256_EXP(1,2,3,4,5,6,7,0,47);
		SHA256_EXP(0,1,2,3,4,5,6,7,48); SHA256_EXP(7,0,1,2,3,4,5,6,49);
		SHA256_EXP(6,7,0,1,2,3,4,5,50); SHA256_EXP(5,6,7,0,1,2,3,4,51);
		SHA256_EXP(4,5,6,7,0,1,2,3,52); SHA256_EXP(3,4,5,6,7,0,1,2,53);
		SHA256_EXP(2,3,4,5,6,7,0,1,54); SHA256_EXP(1,2,3,4,5,6,7,0,55);
		SHA256_EXP(0,1,2,3,4,5,6,7,56); SHA256_EXP(7,0,1,2,3,4,5,6,57);
		SHA256_EXP(6,7,0,1,2,3,4,5,58); SHA256_EXP(5,6,7,0,1,2,3,4,59);
		SHA256_EXP(4,5,6,7,0,1,2,3,60); SHA256_EXP(3,4,5,6,7,0,1,2,61);
		SHA256_EXP(2,3,4,5,6,7,0,1,62); SHA256_EXP(1,2,3,4,5,6,7,0,63);
		
		ctx->h[0] += wv[0]; ctx->h[1] += wv[1];
		ctx->h[2] += wv[2]; ctx->h[3] += wv[3];
		ctx->h[4] += wv[4]; ctx->h[5] += wv[5];
		ctx->h[6] += wv[6]; ctx->h[7] += wv[7];
	}
	
        
	return;
}


void SHA256_update(SHA256_ctx *ctx, const ui1 *message, unsigned int len)
{
	ui4		block_nb;
	ui4		new_len, rem_len, tmp_len;
	const ui1	*shifted_message;
	
	tmp_len = SHA256_BLOCK_SIZE - ctx->len;
	rem_len = len < tmp_len ? len : tmp_len;
	
	memcpy(&ctx->block[ctx->len], message, rem_len);
	
	if (ctx->len + len < SHA256_BLOCK_SIZE) {
		ctx->len += len;
		return;
	}
	
	new_len = len - rem_len;
	block_nb = new_len / SHA256_BLOCK_SIZE;
	
	shifted_message = message + rem_len;
	
	SHA256_transf(ctx, ctx->block, 1);
	SHA256_transf(ctx, shifted_message, block_nb);
	
	rem_len = new_len % SHA256_BLOCK_SIZE;
	
	memcpy(ctx->block, &shifted_message[block_nb << 6], rem_len);
	
	ctx->len = rem_len;
	ctx->tot_len += (block_nb + 1) << 6;
	
        
	return;
}


/*************************************************************************/
/**************************  END SHA-256 FUNCTIONS  **********************/
/*************************************************************************/


void	show_file_processing_struct(FILE_PROCESSING_STRUCT *fps)
{
        si1	hex_str[HEX_STRING_BYTES(4)], *s;
	si4	i;
        
        
	printf("----------- File Processing Structure - START ----------\n");
	UTF8_printf("Full File Name: %s\n", fps->full_file_name);
	if (fps->fd == -1)
		printf("File Descriptor: %d (closed)\n", fps->fd);
	else
		printf("File Descriptor: %d (open)\n", fps->fd);
	printf("File Length: ");
	if (fps->file_length == FPS_FILE_LENGTH_UNKNOWN)
		printf("unknown\n");
	else
		#ifdef _WIN32
			printf("%lld\n", fps->file_length);
		#else
			printf("%ld\n", fps->file_length);
		#endif
	s = (si1 *) &fps->file_type_code;
	generate_hex_string((ui1 *) s, 4, hex_str);
        printf("File Type Code: %s    (", hex_str);
	for (i = 0; i < 4; ++i)
		printf(" %c ", *s++);
	printf(")\n");
	#ifdef _WIN32
		printf("Raw Data Bytes: %lld\n", fps->raw_data_bytes);
	#else
		printf("Raw Data Bytes: %ld\n", fps->raw_data_bytes);
	#endif
	show_password_data(fps);
	show_universal_header(fps);
	if (fps->raw_data_bytes > UNIVERSAL_HEADER_BYTES) {
		switch (fps->file_type_code) {
			case TIME_SERIES_METADATA_FILE_TYPE_CODE:
			case VIDEO_METADATA_FILE_TYPE_CODE:
				show_metadata(fps);
				break;
			case RECORD_DATA_FILE_TYPE_CODE:
				show_records(fps);
				break;
			default:
				break;
		}
	}
	printf("------------ File Processing Structure - END -----------\n\n");
        
        
	return;
}


void	show_metadata(FILE_PROCESSING_STRUCT *fps)
{
	METADATA_SECTION_1		*md1;
	TIME_SERIES_METADATA_SECTION_2	*tmd2;
	VIDEO_METADATA_SECTION_2	*vmd2;
	METADATA_SECTION_3		*md3;
	si1				hex_str[HEX_STRING_BYTES(CRC_BYTES)];
	
	
	// assign
	md1 = fps->metadata.section_1;
	tmd2 = fps->metadata.time_series_section_2;
	vmd2 = fps->metadata.video_section_2;
	md3 = fps->metadata.section_3;
	
	// decrypt if needed
	if (md1->section_2_encryption > NO_ENCRYPTION || md1->section_3_encryption > NO_ENCRYPTION)
		decrypt_metadata(fps);
        
	// show
	printf("------------------- Metadata - START -------------------\n");
	printf("------------------ Section 1 - START -------------------\n");
	printf("Section 2 Encryption: %d ", md1->section_2_encryption);
	if (md1->section_2_encryption == NO_ENCRYPTION)
		printf("(none)\n");
	else if (md1->section_2_encryption == LEVEL_1_ENCRYPTION)
		printf("(level 1, currently encrypted)\n");
	else if (md1->section_2_encryption == LEVEL_2_ENCRYPTION)
		printf("(level 2, currently encrypted)\n");
	else if (md1->section_2_encryption == -LEVEL_1_ENCRYPTION)
		printf("(level 1, currently decrypted)\n");
	else if (md1->section_2_encryption == -LEVEL_2_ENCRYPTION)
		printf("(level 2, currently decrypted)\n");
	else
	printf("(unrecognized code)\n");
	printf("Section 3 Encryption: %d ", md1->section_3_encryption);
	if (md1->section_3_encryption == NO_ENCRYPTION)
		printf("(none)\n");
	else if (md1->section_3_encryption == LEVEL_1_ENCRYPTION)
		printf("(level 1, currently encrypted)\n");
	else if (md1->section_3_encryption == LEVEL_2_ENCRYPTION)
		printf("(level 2, currently encrypted)\n");
	else if (md1->section_3_encryption == -LEVEL_1_ENCRYPTION)
		printf("(level 1, currently decrypted)\n");
	else if (md1->section_3_encryption == -LEVEL_2_ENCRYPTION)
		printf("(level 2, currently decrypted)\n");
	else
		printf("(unrecognized code)\n");
	printf("------------------- Section 1 - END --------------------\n\n");
	printf("------------------ Section 2 - START -------------------\n");
	if (fps->password_data->access_level >= ABS(md1->section_2_encryption)) {
		if (fps->file_type_code == TIME_SERIES_METADATA_FILE_TYPE_CODE) {
                        // type-independent fields
                        if (strlen(tmd2->channel_description))
                                UTF8_printf("Channel Description: %s\n", tmd2->channel_description);
                        else
                                printf("Channel Description: no entry\n");
                        if (strlen(tmd2->session_description))
                                UTF8_printf("Session Description: %s\n", tmd2->session_description);
                        else
                                printf("Session Description: no entry\n");
			if (tmd2->recording_duration == METADATA_RECORDING_DURATION_NO_ENTRY)
				printf("Recording Duration: no entry\n");
			else
				#ifdef _WIN32
					printf("Recording Duration: %lld (microseconds)\n", tmd2->recording_duration);
				#else
					printf("Recording Duration: %ld (microseconds)\n", tmd2->recording_duration);
				#endif
                        // type-specific fields
			if (strlen(tmd2->reference_description))
				UTF8_printf("Reference Description: %s\n", tmd2->reference_description);
			else
				printf("Reference Description: no entry\n");
			if (tmd2->acquisition_channel_number == TIME_SERIES_METADATA_ACQUISITION_CHANNEL_NUMBER_NO_ENTRY)
				printf("Acquisition Channel Number: no entry\n");
			else
				#ifdef _WIN32
					printf("Acquisition Channel Number: %lld\n", tmd2->acquisition_channel_number);
				#else
					printf("Acquisition Channel Number: %ld\n", tmd2->acquisition_channel_number);
				#endif
			if (tmd2->sampling_frequency == TIME_SERIES_METADATA_SAMPLING_FREQUENCY_NO_ENTRY)
				printf("Sampling Frequency: no entry\n");
			else
				printf("Sampling Frequency: %lf\n", tmd2->sampling_frequency);
                        if (tmd2->low_frequency_filter_setting == TIME_SERIES_METADATA_LOW_FREQUENCY_FILTER_SETTING_NO_ENTRY)
                                printf("Low Frequency Filter Setting: no entry\n");
                        else
				printf("Low Frequency Filter Setting: %lf\n", tmd2->low_frequency_filter_setting);
                        if (tmd2->high_frequency_filter_setting == TIME_SERIES_METADATA_HIGH_FREQUENCY_FILTER_SETTING_NO_ENTRY)
                                printf("High Frequency Filter Setting: no entry\n");
                        else
				printf("High Frequency Filter Setting: %lf\n", tmd2->high_frequency_filter_setting);
                        if (tmd2->notch_filter_frequency_setting == TIME_SERIES_METADATA_NOTCH_FILTER_FREQUENCY_SETTING_NO_ENTRY)
                                printf("Notch Filter Frequency Setting: no entry\n");
                        else
				printf("Notch Filter Frequency Setting: %lf\n", tmd2->notch_filter_frequency_setting);
			if (tmd2->AC_line_frequency == TIME_SERIES_METADATA_AC_LINE_FREQUENCY_NO_ENTRY)
				printf("AC Line Frequency: no entry\n");
			else
				printf("AC Line Frequency: %lf\n", tmd2->AC_line_frequency);
			if (tmd2->units_conversion_factor == TIME_SERIES_METADATA_UNITS_CONVERSION_FACTOR_NO_ENTRY)
                                printf("Units Conversion Factor: no entry\n");
                        else
				printf("Units Conversion Factor: %lf\n", tmd2->units_conversion_factor);
                        if (strlen(tmd2->units_description))
				UTF8_printf("Units Description: %s\n", tmd2->units_description);
                        else
                                printf("Units Description: no entry\n");
                        if (tmd2->maximum_native_sample_value == TIME_SERIES_METADATA_MAXIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY)
                                printf("Maximum Native Sample Value: no entry\n");
                        else
				printf("Maximum Native Sample Value: %lf\n", tmd2->maximum_native_sample_value);
                        if (tmd2->minimum_native_sample_value == TIME_SERIES_METADATA_MINIMUM_NATIVE_SAMPLE_VALUE_NO_ENTRY)
                                printf("Minimum Native Sample Value: no entry\n");
                        else
				printf("Minimum Native Sample Value: %lf\n", tmd2->minimum_native_sample_value);
                        if (tmd2->start_sample == TIME_SERIES_METADATA_START_SAMPLE_NO_ENTRY)
                                printf("Start Sample: no entry\n");
                        else
                        	#ifdef _WIN32
					printf("Start Sample: %lld\n", tmd2->start_sample);
				#else
					printf("Start Sample: %ld\n", tmd2->start_sample);
				#endif
                        if (tmd2->number_of_samples == TIME_SERIES_METADATA_NUMBER_OF_SAMPLES_NO_ENTRY)
				printf("Number of Samples: no entry\n");
			else
				#ifdef _WIN32
					printf("Number of Samples: %lld\n", tmd2->number_of_samples);
				#else
					printf("Number of Samples: %ld\n", tmd2->number_of_samples);
				#endif
			if (tmd2->number_of_blocks == TIME_SERIES_METADATA_NUMBER_OF_BLOCKS_NO_ENTRY)
				printf("Number of Blocks: no entry\n");
			else
				#ifdef _WIN32
					printf("Number of Blocks: %lld\n", tmd2->number_of_blocks);
				#else
					printf("Number of Blocks: %ld\n", tmd2->number_of_blocks);
				#endif
			if (tmd2->maximum_block_bytes == TIME_SERIES_METADATA_MAXIMUM_BLOCK_BYTES_NO_ENTRY)
				printf("Maximum Block Bytes: no entry\n");
			else
				#ifdef _WIN32
					printf("Maximum Block Bytes: %lld\n", tmd2->maximum_block_bytes);
				#else
					printf("Maximum Block Bytes: %ld\n", tmd2->maximum_block_bytes);
				#endif
			if (tmd2->maximum_block_samples == TIME_SERIES_METADATA_MAXIMUM_BLOCK_SAMPLES_NO_ENTRY)
				printf("Maximum Block Samples: no entry\n");
			else
				printf("Maximum Block Samples: %u\n", tmd2->maximum_block_samples);
                        if (tmd2->maximum_difference_bytes == TIME_SERIES_METADATA_MAXIMUM_DIFFERENCE_BYTES_NO_ENTRY)
                                printf("Maximum Difference Bytes: no entry\n");
                        else
				printf("Maximum Difference Bytes: %u\n", tmd2->maximum_difference_bytes);
                       if (tmd2->block_interval == TIME_SERIES_METADATA_BLOCK_INTERVAL_NO_ENTRY)
                                printf("Block Interval: no entry\n");
                        else
				#ifdef _WIN32
					printf("Block Interval: %lld (microseconds)\n", tmd2->block_interval);
				#else
					printf("Block Interval: %ld (microseconds)\n", tmd2->block_interval);
				#endif
			if (tmd2->number_of_discontinuities == TIME_SERIES_METADATA_NUMBER_OF_DISCONTINUITIES_NO_ENTRY)
				printf("Number of Discontinuities: no entry\n");
			else
				#ifdef _WIN32
					printf("Number of Discontinuities: %lld\n", tmd2->number_of_discontinuities);
				#else
					printf("Number of Discontinuities: %ld\n", tmd2->number_of_discontinuities);
				#endif
			if (tmd2->maximum_contiguous_blocks == TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCKS_NO_ENTRY)
				printf("Maximum Contiguous Blocks: no entry\n");
			else
				#ifdef _WIN32
					printf("Maximum Contiguous Blocks: %lld\n", tmd2->maximum_contiguous_blocks);
				#else
					printf("Maximum Contiguous Blocks: %ld\n", tmd2->maximum_contiguous_blocks);
				#endif
			if (tmd2->maximum_contiguous_block_bytes == TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_BLOCK_BYTES_NO_ENTRY)
				printf("Maximum Contiguous Block Bytes: no entry\n");
			else
				#ifdef _WIN32
					printf("Maximum Contiguous Block Bytes: %lld\n", tmd2->maximum_contiguous_block_bytes);
				#else
					printf("Maximum Contiguous Block Bytes: %ld\n", tmd2->maximum_contiguous_block_bytes);
				#endif
			if (tmd2->maximum_contiguous_samples == TIME_SERIES_METADATA_MAXIMUM_CONTIGUOUS_SAMPLES_NO_ENTRY)
                                printf("Maximum Contiguous Samples: no entry\n");
                        else
				#ifdef _WIN32
					printf("Maximum Contiguous Samples: %lld\n", tmd2->maximum_contiguous_samples);
				#else
					printf("Maximum Contiguous Samples: %ld\n", tmd2->maximum_contiguous_samples);
				#endif
		}
		else if (fps->file_type_code == VIDEO_METADATA_FILE_TYPE_CODE) {
                        // type-independent fields
                        if (strlen(vmd2->channel_description))
                                UTF8_printf("Channel Description: %s\n", vmd2->channel_description);
                        else
                                printf("Channel Description: no entry\n");
                        if (strlen(vmd2->session_description))
                                UTF8_printf("Session Description: %s\n", vmd2->session_description);
                        else
                                printf("Session Description: no entry\n");
                       if (vmd2->recording_duration == METADATA_RECORDING_DURATION_NO_ENTRY)
                                printf("Recording Duration: no entry\n");
                        else
                        	#ifdef _WIN32
					printf("Recording Duration: %lld (microseconds)\n", vmd2->recording_duration);
				#else
					printf("Recording Duration: %ld (microseconds)\n", vmd2->recording_duration);
				#endif
                        // type-specific fields
			if (vmd2->horizontal_resolution == VIDEO_METADATA_HORIZONTAL_RESOLUTION_NO_ENTRY)
				printf("Horizontal Resolution: no entry\n");
			else
				#ifdef _WIN32
					printf("Horizontal Resolution: %lld\n", vmd2->horizontal_resolution);
				#else
					printf("Horizontal Resolution: %ld\n", vmd2->horizontal_resolution);
				#endif
			if (vmd2->vertical_resolution == VIDEO_METADATA_VERTICAL_RESOLUTION_NO_ENTRY)
				printf("Vertical Resolution: no entry\n");
			else
				#ifdef _WIN32
					printf("Vertical Resolution: %lld\n", vmd2->vertical_resolution);
				#else
					printf("Vertical Resolution: %ld\n", vmd2->vertical_resolution);
				#endif
                        if (vmd2->frame_rate == VIDEO_METADATA_FRAME_RATE_NO_ENTRY)
				printf("Frame Rate: no entry\n");
			else
				printf("Frame Rate: %lf (frames per second)\n", vmd2->frame_rate);
			if (vmd2->number_of_clips == VIDEO_METADATA_NUMBER_OF_CLIPS_NO_ENTRY)
				printf("Number of Clips: no entry\n");
			else
				#ifdef _WIN32
					printf("Number of Clips: %lld (= number of video indices)\n", vmd2->number_of_clips);
				#else
					printf("Number of Clips: %ld (= number of video indices)\n", vmd2->number_of_clips);
				#endif
			if (vmd2->maximum_clip_bytes == VIDEO_METADATA_MAXIMUM_CLIP_BYTES_NO_ENTRY)
				printf("Maximum Clip Bytes: no entry\n");
			else
				#ifdef _WIN32
					printf("Maximum Clip Bytes: %lld\n", vmd2->maximum_clip_bytes);
				#else
					printf("Maximum Clip Bytes: %ld\n", vmd2->maximum_clip_bytes);
				#endif
			if (strlen(vmd2->video_format))
				UTF8_printf("Video Format: %s\n", vmd2->video_format);
			else
				printf("Video Format: no entry\n");
			if (vmd2->video_file_CRC == CRC_NO_ENTRY)
				printf("Video File CRC: no entry\n");
			else {
				generate_hex_string((ui1 *) &vmd2->video_file_CRC, TYPE_BYTES, hex_str);
				printf("Video File CRC: %s\n", hex_str);
			}
		}
		else {
			printf("(unrecognized metadata section 2 type)\n");
		}
	} else {
		printf("No access to section 2\n");
	}
	printf("------------------- Section 2 - END --------------------\n\n");
	printf("------------------ Section 3 - START -------------------\n");
	if (fps->password_data->access_level >= ABS(md1->section_3_encryption)) {
		if (md3->recording_time_offset == UUTC_NO_ENTRY)
			printf("Recording Time Offset: no entry\n");
		else
			#ifdef _WIN32
				printf("Recording Time Offset: %lld\n", md3->recording_time_offset);
			#else
				printf("Recording Time Offset: %ld\n", md3->recording_time_offset);
			#endif
		if (md3->DST_start_time == UUTC_NO_ENTRY)
			printf("DST Start Time: no entry\n");
		else
			#ifdef _WIN32
				printf("DST Start Time: %lld\n", md3->DST_start_time);
			#else
				printf("DST Start Time: %ld\n", md3->DST_start_time);
			#endif
		if (md3->DST_end_time == UUTC_NO_ENTRY)
			printf("DST End Time: no entry\n");
		else
			#ifdef _WIN32
				printf("DST End Time: %lld\n", md3->DST_end_time);
			#else
				printf("DST End Time: %ld\n", md3->DST_end_time);
			#endif
		if (md3->GMT_offset == GMT_OFFSET_NO_ENTRY)
			printf("GMT: no entry\n");
		else
			printf("GMT Offset: %d\n", md3->GMT_offset);
		if (strlen(md3->subject_name_1))
			UTF8_printf("Subject Name 1: %s\n", md3->subject_name_1);
		else
			printf("Subject Name 1: no entry\n");
		if (strlen(md3->subject_name_2))
			UTF8_printf("Subject Name 2: %s\n", md3->subject_name_2);
		else
			printf("Subject Name 2: no entry\n");
		if (strlen(md3->subject_ID))
			UTF8_printf("Subject ID: %s\n", md3->subject_ID);
		else
			printf("Subject ID: no entry\n");
		if (strlen(md3->recording_location))
			UTF8_printf("Recording Location: %s\n", md3->recording_location);
		else
			printf("Recording Location: no entry\n");
	} else {
		printf("No access to section 3\n");
	}
	printf("------------------- Section 3 - END --------------------\n\n");
	printf("-------------------- Metadata - END --------------------\n\n");
        
	
	return;
}


void	show_password_data(FILE_PROCESSING_STRUCT *fps)
{
	PASSWORD_DATA	*pwd;
        si1		hex_str[HEX_STRING_BYTES(ENCRYPTION_KEY_BYTES)];
	
	
	pwd = fps->password_data;
	
	printf("------------------ Password Data - START -----------------\n");
	if (pwd->access_level >= LEVEL_1_ACCESS) {
                generate_hex_string(pwd->level_1_encryption_key, ENCRYPTION_KEY_BYTES, hex_str);
		printf("Level 1 Encryption Key: %s\n", hex_str);
	}
	if (pwd->access_level == LEVEL_2_ACCESS) {
                generate_hex_string(pwd->level_2_encryption_key, ENCRYPTION_KEY_BYTES, hex_str);
		printf("Level 2 Encryption Key: %s\n", hex_str);
	}
	printf("Access Level: %d\n", pwd->access_level);
	printf("------------------- Password Data - END ------------------\n\n");
	
        
	return;
}


void	show_records(FILE_PROCESSING_STRUCT *fps)
{
	si8		number_of_records;
        ui4		i, r_cnt;
        ui1		*ui1_p, *end_p;
        RECORD_HEADER	*record_header;
        PASSWORD_DATA	*pwd;
        
	
	number_of_records = fps->universal_header->number_of_entries;
        ui1_p = fps->raw_data + UNIVERSAL_HEADER_BYTES;
        pwd = fps->password_data;
	
	// number_of_records obtained from session metadata - use if passed
	if (number_of_records == UNKNOWN_NUMBER_OF_ENTRIES) {   // can still process if not passed, but will fail on incomplete final record
		end_p = fps->raw_data + fps->file_length;
		r_cnt = 0;
		while (ui1_p < end_p) {
			record_header = (RECORD_HEADER *) ui1_p;
			show_record(record_header, r_cnt, pwd);
			ui1_p += (RECORD_HEADER_BYTES + record_header->bytes);
			++r_cnt;
		}
		fps->universal_header->number_of_entries = r_cnt;
	} else {
                for (i = 0; i < number_of_records; ++i) {
			record_header = (RECORD_HEADER *) ui1_p;
			show_record(record_header, i, pwd);
			ui1_p += (RECORD_HEADER_BYTES + record_header->bytes);
		}
	}
        
        
        return;
}


void	show_universal_header(FILE_PROCESSING_STRUCT *fps)
{
	UNIVERSAL_HEADER	*uh;
        si1			hex_str[HEX_STRING_BYTES(PASSWORD_VALIDATION_FIELD_BYTES)], time_str[TIME_STRING_BYTES];
	ui4			type_code, *file_type_string_int;
	
	
	uh = fps->universal_header;
	file_type_string_int = (ui4 *) uh->file_type_string;
	type_code = *file_type_string_int;
	
	printf("---------------- Universal Header - START ----------------\n");
	if (uh->header_CRC == CRC_NO_ENTRY)
		printf("Header CRC: no entry\n");
	else {
                generate_hex_string((ui1 *) &uh->header_CRC, CRC_BYTES, hex_str);
                printf("Header CRC: %s\n", hex_str);
        }
	if (uh->body_CRC == CRC_NO_ENTRY)
		printf("Body CRC: no entry\n");
	else {
		generate_hex_string((ui1 *) &uh->body_CRC, CRC_BYTES, hex_str);
		printf("Body CRC: %s\n", hex_str);
	}
	if (strlen(uh->file_type_string))
		printf("File Type String: %s\n", uh->file_type_string);
	else
		printf("File Type String: no entry\n");
	if (uh->mef_version_major == UNIVERSAL_HEADER_MEF_VERSION_MAJOR_NO_ENTRY || uh->mef_version_minor == UNIVERSAL_HEADER_MEF_VERSION_MINOR_NO_ENTRY) {
		if (uh->mef_version_major == UNIVERSAL_HEADER_MEF_VERSION_MAJOR_NO_ENTRY)
			printf("MEF Version Major: no entry\n");
		else
			printf("MEF Version Major: %u\n", uh->mef_version_major);
		if (uh->mef_version_minor == UNIVERSAL_HEADER_MEF_VERSION_MINOR_NO_ENTRY)
			printf("MEF Version Minor: no entry\n");
		else
			printf("MEF Version Minor: %u\n", uh->mef_version_minor);
	} else
		printf("MEF Version: %u.%u\n", uh->mef_version_major, uh->mef_version_minor);
	if (uh->byte_order_code == UNIVERSAL_HEADER_BYTE_ORDER_CODE_NO_ENTRY)
		printf("Byte Order Code: no entry ");
	else {
		printf("Byte Order Code: %u ", uh->byte_order_code);
		if (uh->byte_order_code == MEF_LITTLE_ENDIAN)
			printf("(little endian)\n");
		else if (uh->byte_order_code == MEF_BIG_ENDIAN)
			printf("(big endian)\n");
		else
			printf("(unrecognized code)\n");
	}
	if (uh->start_time == UUTC_NO_ENTRY)
		printf("Start Time: no entry\n");
	else {
		local_date_time_string(uh->start_time, time_str);
		if (uh->start_time < 0)
			printf("Offset ");
		#ifdef _WIN32
			printf("Start Time: %lld (uUTC), %s (ascii, local)\n", ABS(uh->start_time), time_str);
		#else
			printf("Start Time: %ld (uUTC), %s (ascii, local)\n", ABS(uh->start_time), time_str);
		#endif
	}
	if (uh->end_time == UUTC_NO_ENTRY)
		printf("End Time: no entry\n");
	else {
		local_date_time_string(uh->end_time, time_str);
		if (uh->start_time < 0)
			printf("Offset ");
		#ifdef _WIN32
			printf("End Time: %lld (uUTC), %s (ascii, local)\n", ABS(uh->end_time), time_str);
		#else
			printf("End Time: %ld (uUTC), %s (ascii, local)\n", ABS(uh->end_time), time_str);
		#endif
	}
	if (uh->number_of_entries == UNIVERSAL_HEADER_NUMBER_OF_ENTRIES_NO_ENTRY)
		printf("Number of Entries: no entry\n");
	else {
		#ifdef _WIN32
			printf("Number of Entries: %lld  ", uh->number_of_entries);
		#else
			printf("Number of Entries: %ld  ", uh->number_of_entries);
		#endif
		switch (type_code) {
			case RECORD_DATA_FILE_TYPE_CODE:
				printf("(number of records in the file)\n");
				break;
			case RECORD_INDICES_FILE_TYPE_CODE:
				printf("(number of records indices in the file)\n");
				break;
			case TIME_SERIES_METADATA_FILE_TYPE_CODE:
			case VIDEO_METADATA_FILE_TYPE_CODE:
				printf("(one metadata entry per metadata file)\n");
				break;
			case VIDEO_INDICES_FILE_TYPE_CODE:
				printf("(number of video indices in the file)\n");
				break;
			case TIME_SERIES_DATA_FILE_TYPE_CODE:
				printf("(number of RED blocks in the file)\n");
				break;
			case TIME_SERIES_INDICES_FILE_TYPE_CODE:
				printf("(number of time series indices in the file)\n");
				break;
			default:
				printf("\n");
				break;
		}
	}
	if (uh->maximum_entry_size == UNIVERSAL_HEADER_MAXIMUM_ENTRY_SIZE_NO_ENTRY)
		printf("Maximum Entry Size: no entry\n");
	else {
		#ifdef _WIN32
			printf("Maximum Entry Size: %lld  ", uh->maximum_entry_size);
		#else
			printf("Maximum Entry Size: %ld  ", uh->maximum_entry_size);
		#endif
		switch (type_code) {
			case RECORD_DATA_FILE_TYPE_CODE:
				printf("(number of bytes in the largest record in the file)\n");
				break;
			case RECORD_INDICES_FILE_TYPE_CODE:
				printf("(number of bytes in a record index)\n");
				break;
			case TIME_SERIES_METADATA_FILE_TYPE_CODE:
			case VIDEO_METADATA_FILE_TYPE_CODE:
				printf("(number of bytes in a metadata file)\n");
				break;
			case VIDEO_INDICES_FILE_TYPE_CODE:
				printf("(number of bytes in the largest clip in the video data file)\n");
				break;
			case TIME_SERIES_DATA_FILE_TYPE_CODE:
				printf("(number of samples in the largest RED block in the file)\n");
				break;
			case TIME_SERIES_INDICES_FILE_TYPE_CODE:
				printf("(number of bytes in a time series index)\n");
				break;
			default:
				printf("\n");
				break;
				
		}
	}
	if (uh->segment_number == UNIVERSAL_HEADER_SEGMENT_NUMBER_NO_ENTRY)
		printf("Segment Number: no entry\n");
        else if (uh->segment_number == UNIVERSAL_HEADER_CHANNEL_LEVEL_CODE)
                printf("Segment Number: channel level\n");
        else if (uh->segment_number == UNIVERSAL_HEADER_SESSION_LEVEL_CODE)
                printf("Segment Number: session level\n");
	else
		printf("Segment Number: %d\n", uh->segment_number);
	if (strlen(uh->channel_name))
		UTF8_printf("Channel Name: %s\n", uh->channel_name);
	else
		printf("Channel Name: no entry\n");
	if (strlen(uh->session_name))
		UTF8_printf("Session Name: %s\n", uh->session_name);
	else
		printf("Session Name: no entry\n");
	if (strlen(uh->anonymized_name))
		UTF8_printf("Anonymized Name: %s\n", uh->anonymized_name);
	else
		printf("Anonymized Name: no entry\n");
	if (all_zeros(uh->level_UUID, UUID_BYTES) == MEF_TRUE)
		printf("Level UUID: no entry\n");
	else {
		generate_hex_string(uh->level_UUID, UUID_BYTES, hex_str);
		printf("Level UUID: %s\n", hex_str);
	}
	if (all_zeros(uh->file_UUID, UUID_BYTES) == MEF_TRUE)
		printf("File UUID: no entry\n");
	else {
		generate_hex_string(uh->file_UUID, UUID_BYTES, hex_str);
		printf("File UUID: %s\n", hex_str);
	}
	if (all_zeros(uh->provenance_UUID, UUID_BYTES) == MEF_TRUE)
		printf("Provenance UUID: no entry\n");
	else {
		generate_hex_string(uh->provenance_UUID, UUID_BYTES, hex_str);
		printf("Provenance UUID: %s  ", hex_str);
		if (memcmp(uh->provenance_UUID, uh->file_UUID, UUID_BYTES))
			printf("(derived data\n)");
		else
			printf("(original data\n)");
	}
	if (all_zeros(uh->level_1_password_validation_field, PASSWORD_VALIDATION_FIELD_BYTES) == MEF_TRUE)
		printf("Level 1 Password Validation_Field: no entry\n");
	else {
                generate_hex_string(uh->level_1_password_validation_field, PASSWORD_VALIDATION_FIELD_BYTES, hex_str);
		printf("Level 1 Password Validation_Field: %s\n", hex_str);
	}
	if (all_zeros(uh->level_2_password_validation_field, PASSWORD_VALIDATION_FIELD_BYTES) == MEF_TRUE)
		printf("Level 2 Password Validation_Field: no entry\n");
	else {
                generate_hex_string(uh->level_2_password_validation_field, PASSWORD_VALIDATION_FIELD_BYTES, hex_str);
		printf("Level 2 Password Validation_Field: %s\n", hex_str);
	}
	printf("---------------- Universal Header - END ----------------\n\n");
	
        
	return;
}


si4     sort_by_idx(const void *n1, const void *n2)
{
	si4     i1, i2;
	
	
	i1 = ((NODE *) n1)->idx;
	i2 = ((NODE *) n2)->idx;
	
	return(i1 - i2);
}


si4     sort_by_val(const void *n1, const void *n2)
{
	sf8     v1, v2;
	
	
	v1 = ((NODE *) n1)->val;
	v2 = ((NODE *) n2)->val;
	
	if (v1 > v2)
		return(1);
	else if (v1 < v2)
		return(-1);
	return(0);
}


/*************************************************************************/
/********************************  UTF-8 FUNCTIONS  **********************/
/*************************************************************************/

// ATTRIBUTION
//
// Basic UTF-8 manipulation routines
// by Jeff Bezanson
// placed in the public domain Fall 2005

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


// byte offset => charnum
si4 UTF8_charnum(si1 *s, si4 offset)
{
	si4	charnum = 0, offs = 0;
	
	while (offs < offset && s[offs]) {
		(void)(isutf(s[++offs]) || isutf(s[++offs]) || isutf(s[++offs]) || ++offs);
		charnum++;
	}
	
        
	return(charnum);
}


inline void UTF8_dec(si1 *s, si4 *i)
{
	(void) (isutf(s[--(*i)]) || isutf(s[--(*i)]) || isutf(s[--(*i)]) || --(*i));
        
        
        return;
}


si4 UTF8_escape(si1 *buf, si4 sz, si1 *src, si4 escape_quotes)
{
	si4	c = 0, i = 0, amt;
	
	while (src[i] && c < sz) {
		if (escape_quotes && src[i] == '"') {
			amt = snprintf(buf, sz - c, "\\\"");
			i++;
		}
		else {
			amt = UTF8_escape_wchar(buf, sz - c, UTF8_nextchar(src, &i));
		}
		c += amt;
		buf += amt;
	}
	if (c < sz)
		*buf = '\0';
	
        
	return(c);
}


si4 UTF8_escape_wchar(si1 *buf, si4 sz, ui4 ch)
{
	if (ch == L'\n')
		return(snprintf(buf, sz, "\\n"));
	else if (ch == L'\t')
		return(snprintf(buf, sz, "\\t"));
	else if (ch == L'\r')
		return(snprintf(buf, sz, "\\r"));
	else if (ch == L'\b')
		return(snprintf(buf, sz, "\\b"));
	else if (ch == L'\f')
		return(snprintf(buf, sz, "\\f"));
	else if (ch == L'\v')
		return(snprintf(buf, sz, "\\v"));
	else if (ch == L'\a')
		return(snprintf(buf, sz, "\\a"));
	else if (ch == L'\\')
		return(snprintf(buf, sz, "\\\\"));
	else if (ch < 32 || ch == 0x7f)
		return(snprintf(buf, sz, "\\x%hhX", (ui1) ch));
	else if (ch > 0xFFFF)
		return(snprintf(buf, sz, "\\U%.8X", (ui4) ch));
	else if (ch >= 0x80 && ch <= 0xFFFF)
		return(snprintf(buf, sz, "\\u%.4hX", (ui2) ch));
	
        
	return(snprintf(buf, sz, "%c", (si1) ch));
}


si4 UTF8_fprintf(FILE *stream, si1 *fmt, ...)
{
	si4	cnt;
	va_list args;
	
	va_start(args, fmt);
	
	cnt = UTF8_vfprintf(stream, fmt, args);
	
	va_end(args);
	
        
	return(cnt);
}


inline si4 UTF8_hex_digit(si1 c)
{
	return((c >= '0' && c <= '9') || (c >= 'A' && c <= 'F') || (c >= 'a' && c <= 'f'));
}


inline void UTF8_inc(si1 *s, si4 *i)
{
	(void) (isutf(s[++(*i)]) || isutf(s[++(*i)]) || isutf(s[++(*i)]) || ++(*i));
}


ui4	*UTF8_initialize_offsets_from_UTF8_table(si4 global_flag)
{
	ui4	*offsetsFromUTF8_table;
	
	
	offsetsFromUTF8_table = (ui4 *) e_calloc((size_t) OFFSETS_FROM_UTF8_TABLE_ENTRIES, sizeof(ui4), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	{
		ui4 temp[OFFSETS_FROM_UTF8_TABLE_ENTRIES] = OFFSETS_FROM_UTF8;
		memcpy(offsetsFromUTF8_table, temp, OFFSETS_FROM_UTF8_TABLE_ENTRIES * sizeof(ui4));
	}
	
	if (global_flag == MEF_TRUE) {
		MEF_globals->UTF8_offsets_from_UTF8_table = offsetsFromUTF8_table;
		return(NULL);
	}
	
        
	return(offsetsFromUTF8_table);
}


si1	*UTF8_initialize_trailing_bytes_for_UTF8_table(si4 global_flag)
{
	si1	*trailingBytesForUTF8_table;
	
	
	trailingBytesForUTF8_table = (si1 *) e_calloc((size_t) TRAILING_BYTES_FOR_UTF8_TABLE_ENTRIES, sizeof(si1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	{
		ui4 temp[TRAILING_BYTES_FOR_UTF8_TABLE_ENTRIES] = TRAILING_BYTES_FOR_UTF8;
		memcpy(trailingBytesForUTF8_table, temp, TRAILING_BYTES_FOR_UTF8_TABLE_ENTRIES * sizeof(si1));
	}
	
	if (global_flag == MEF_TRUE) {
		MEF_globals->UTF8_trailing_bytes_for_UTF8_table = trailingBytesForUTF8_table;
		return(NULL);
	}
	
        
	return(trailingBytesForUTF8_table);
}


si4 UTF8_is_locale_utf8(si1 *locale)
{
	// this code based on libutf8
	const si1	*cp = locale;
	
        
	for (; *cp != '\0' && *cp != '@' && *cp != '+' && *cp != ','; cp++) {
		if (*cp == '.') {
			const si1 *encoding = ++cp;
			for (; *cp != '\0' && *cp != '@' && *cp != '+' && *cp != ','; cp++)
				;
			if ((cp - encoding == 5 && !strncmp(encoding, "UTF-8", 5)) || (cp - encoding == 4 && !strncmp(encoding, "utf8", 4)))
				return(1); // it's UTF-8
			break;
		}
	}
	
        
	return(0);
}


si1 *UTF8_memchr(si1 *s, ui4 ch, size_t sz, si4 *charn)
{
	si4	i = 0, lasti = 0;
	ui4	c;
	si4	csz;
	
        
	if (MEF_globals->UTF8_offsets_from_UTF8_table == NULL)
		(void) UTF8_initialize_offsets_from_UTF8_table(MEF_TRUE);
	
	*charn = 0;
	while (i < sz) {
		c = csz = 0;
		do {
			c <<= 6;
			c += (ui1) s[i++];
			csz++;
		} while (i < sz && !isutf(s[i]));
		c -= MEF_globals->UTF8_offsets_from_UTF8_table[csz - 1];
		
		if (c == ch) {
			return(&s[lasti]);
		}
		lasti = i;
		(*charn)++;
	}
	
        
	return(NULL);
}


// reads the next utf-8 sequence out of a string, updating an index
ui4 UTF8_nextchar(si1 *s, si4 *i)
{
	ui4	ch = 0;
	si4	sz = 0;
	
        
	if (MEF_globals->UTF8_offsets_from_UTF8_table == NULL)
		(void) UTF8_initialize_offsets_from_UTF8_table(MEF_TRUE);
	
	do {
		ch <<= 6;
		ch += (ui1) s[(*i)++];
		sz++;
	} while (s[*i] && !isutf(s[*i]));
	
	ch -= MEF_globals->UTF8_offsets_from_UTF8_table[sz - 1];
	
        
	return(ch);
}


inline si4 UTF8_octal_digit(si1 c)
{
	return(c >= '0' && c <= '7');
}


// charnum => byte offset
si4 UTF8_offset(si1 *str, si4 charnum)
{
	si4	offs = 0;
	
        
	while (charnum > 0 && str[offs]) {
		(void)(isutf(str[++offs]) || isutf(str[++offs]) || isutf(str[++offs]) || ++offs);
		charnum--;
	}
	
        
	return(offs);
}


si4 UTF8_printf(si1 *fmt, ...)
{
	si4	cnt;
	va_list args;
	
        
	va_start(args, fmt);
	
	cnt = UTF8_vprintf(fmt, args);
	
	va_end(args);
	
        
	return(cnt);
}


// assumes that src points to the character after a backslash
// returns number of input characters processed
si4 UTF8_read_escape_sequence(si1 *str, ui4 *dest)
{
	ui4	ch;
	si1	digs[9] = "\0\0\0\0\0\0\0\0";
	si4	dno = 0, i = 1;
	
        
	ch = (ui4) str[0];    // take literal character
	if (str[0] == 'n')
		ch = L'\n';
	else if (str[0] == 't')
		ch = L'\t';
	else if (str[0] == 'r')
		ch = L'\r';
	else if (str[0] == 'b')
		ch = L'\b';
	else if (str[0] == 'f')
		ch = L'\f';
	else if (str[0] == 'v')
		ch = L'\v';
	else if (str[0] == 'a')
		ch = L'\a';
	else if (UTF8_octal_digit(str[0])) {
		i = 0;
		do {
			digs[dno++] = str[i++];
		} while (UTF8_octal_digit(str[i]) && dno < 3);
		ch = strtol(digs, NULL, 8);
	}
	else if (str[0] == 'x') {
		while (UTF8_hex_digit(str[i]) && dno < 2) {
			digs[dno++] = str[i++];
		}
		if (dno > 0)
			ch = strtol(digs, NULL, 16);
	}
	else if (str[0] == 'u') {
		while (UTF8_hex_digit(str[i]) && dno < 4) {
			digs[dno++] = str[i++];
		}
		if (dno > 0)
			ch = strtol(digs, NULL, 16);
	}
	else if (str[0] == 'U') {
		while (UTF8_hex_digit(str[i]) && dno < 8) {
			digs[dno++] = str[i++];
		}
		if (dno > 0)
			ch = strtol(digs, NULL, 16);
	}
	*dest = ch;
	
        
	return(i);
}


// returns length of next utf-8 sequence
inline si4	UTF8_seqlen(si1 *s)
{
	if (MEF_globals->UTF8_trailing_bytes_for_UTF8_table == NULL)
		(void) UTF8_initialize_trailing_bytes_for_UTF8_table(MEF_TRUE);
	
        
	return(MEF_globals->UTF8_trailing_bytes_for_UTF8_table[(si4) (ui1) s[0]] + 1);
}


si1 *UTF8_strchr(si1 *s, ui4 ch, si4 *charn)
{
	si4	i = 0, lasti = 0;
	ui4	c;
	
        
	*charn = 0;
	while (s[i]) {
		c = UTF8_nextchar(s, &i);
		if (c == ch) {
			return(&s[lasti]);
		}
		lasti = i;
		(*charn)++;
	}
	
        
	return(NULL);
}


// number of characters
si4 UTF8_strlen(si1 *s)
{
	si4	count = 0;
	si4	i = 0;
	
        
	while (UTF8_nextchar(s, &i) != 0)
		count++;
	
        
	return(count);
}


// conversions without error checking
// only works for valid UTF-8, i.e. no 5- or 6-byte sequences
// srcsz = source size in bytes, or -1 if 0-terminated
// sz = dest size in # of wide characters

// returns # characters converted
// dest will always be L'\0'-terminated, even if there isn't enough room
// for all the characters.
// if sz = srcsz+1 (i.e. 4*srcsz+4 bytes), there will always be enough space
si4 UTF8_toucs(ui4 *dest, si4 sz, si1 *src, si4 srcsz)
{
	ui4	ch;
	si1	*src_end = src + srcsz;
	si4	nb;
	si4	i = 0;
	
        
	if (MEF_globals->UTF8_offsets_from_UTF8_table == NULL)
		(void) UTF8_initialize_offsets_from_UTF8_table(MEF_TRUE);
	
	if (MEF_globals->UTF8_trailing_bytes_for_UTF8_table == NULL)
		(void) UTF8_initialize_trailing_bytes_for_UTF8_table(MEF_TRUE);
	
	while (i < sz - 1) {
		nb = MEF_globals->UTF8_trailing_bytes_for_UTF8_table[(ui1) *src];
		if (srcsz == -1) {
			if (*src == 0)
				goto UTF8_DONE_TOUCS;
		}
		else {
			if (src + nb >= src_end)
				goto UTF8_DONE_TOUCS;
		}
		ch = 0;
		switch (nb) {
				// these fall through deliberately
			case 3: ch += (ui1) *src++; ch <<= 6;
			case 2: ch += (ui1) *src++; ch <<= 6;
			case 1: ch += (ui1) *src++; ch <<= 6;
			case 0: ch += (ui1) *src++;
		}
		ch -= MEF_globals->UTF8_offsets_from_UTF8_table[nb];
		dest[i++] = ch;
	}
	
	UTF8_DONE_TOUCS:
	
	dest[i] = 0;
	
        
	return(i);
}


// srcsz = number of source characters, or -1 if 0-terminated
// sz = size of dest buffer in bytes

// returns # characters converted
// dest will only be '\0'-terminated if there is enough space. this is
// for consistency; imagine there are 2 bytes of space left, but the next
// character requires 3 bytes. in this case we could NUL-terminate, but in
// general we can't when there's insufficient space. therefore this function
// only NUL-terminates if all the characters fit, and there's space for
// the NUL as well.
// the destination string will never be bigger than the source string
si4 UTF8_toutf8(si1 *dest, si4 sz, ui4 *src, si4 srcsz)
{
	ui4	ch;
	si4	i = 0;
	si1	*dest_end = dest + sz;
	
        
	while (srcsz < 0 ? src[i] != 0 : i < srcsz) {
		ch = src[i];
		if (ch < 0x80) {
			if (dest >= dest_end)
				return(i);
			*dest++ = (si1) ch;
		}
		else if (ch < 0x800) {
			if (dest >= dest_end - 1)
				return(i);
			*dest++ = (ch >> 6) | 0xC0;
			*dest++ = (ch & 0x3F) | 0x80;
		}
		else if (ch < 0x10000) {
			if (dest >= dest_end - 2)
				return(i);
			*dest++ = (ch >> 12) | 0xE0;
			*dest++ = ((ch >> 6) & 0x3F) | 0x80;
			*dest++ = (ch & 0x3F) | 0x80;
		}
		else if (ch < 0x110000) {
			if (dest >= dest_end - 3)
				return(i);
			*dest++ = (ch >> 18) | 0xF0;
			*dest++ = ((ch >> 12) & 0x3F) | 0x80;
			*dest++ = ((ch >> 6) & 0x3F) | 0x80;
			*dest++ = (ch & 0x3F) | 0x80;
		}
		i++;
	}
	if (dest < dest_end)
		*dest = '\0';
	
        
	return(i);
}


// convert a string with literal \uxxxx or \Uxxxxxxxx characters to UTF-8
// example: UTF8_unescape(mybuf, 256, "hello\\u220e")
// note the double backslash is needed if called on a C string literal
si4 UTF8_unescape(si1 *buf, si4 sz, si1 *src)
{
	si4	c = 0, amt;
	ui4	ch;
	si1	temp[4];
	
        
	while (*src && c < sz) {
		if (*src == '\\') {
			src++;
			amt = UTF8_read_escape_sequence(src, &ch);
		}
		else {
			ch = (u_int32_t)*src;
			amt = 1;
		}
		src += amt;
		amt = UTF8_wc_toutf8(temp, ch);
		if (amt > sz - c)
			break;
		memcpy(&buf[c], temp, amt);
		c += amt;
	}
	if (c < sz)
		buf[c] = '\0';
	
        
	return(c);
}


si4 UTF8_vfprintf(FILE *stream, si1 *fmt, va_list ap)
{
	si4	cnt, sz = 512;
	si1	*buf;
	ui4	*wcs;
	
	
	buf = (si1 *) alloca(sz);
	
	UTF8_TRY_FPRINT:
        
	cnt = vsnprintf(buf, sz, fmt, ap);
	if (cnt >= sz) {
		buf = (si1 *) alloca(cnt - sz + 1);
		sz = cnt + 1;
		goto UTF8_TRY_FPRINT;
	}
	wcs = (ui4 *) alloca((cnt + 1) * sizeof(ui4));
	cnt = UTF8_toucs(wcs, cnt + 1, buf, cnt);
	fprintf(stream, "%ls", (wchar_t *) wcs);
	
        
	return(cnt);
}


si4 UTF8_vprintf(si1 *fmt, va_list ap)
{
	si4	cnt, sz = 512;
	si1	*buf;
	ui4	*wcs;
	
        
	buf = (si1 *) alloca(sz);
	
	UTF8_TRY_PRINT:
        
	cnt = vsnprintf(buf, sz, fmt, ap);
	if (cnt >= sz) {
		buf = (si1 *) alloca(cnt - sz + 1);
		sz = cnt + 1;
		goto UTF8_TRY_PRINT;
	}
	wcs = (ui4 *) alloca((cnt + 1) * sizeof(ui4));
	cnt = UTF8_toucs(wcs, cnt + 1, buf, cnt);
	printf("%ls", (wchar_t *) wcs);
	
        
	return(cnt);
}


si4 UTF8_wc_toutf8(si1 *dest, ui4 ch)
{
	if (ch < 0x80) {
		dest[0] = (char)ch;
		return(1);
	}
	if (ch < 0x800) {
		dest[0] = (ch >> 6) | 0xC0;
		dest[1] = (ch & 0x3F) | 0x80;
		return(2);
	}
	if (ch < 0x10000) {
		dest[0] = (ch >> 12) | 0xE0;
		dest[1] = ((ch >> 6) & 0x3F) | 0x80;
		dest[2] = (ch & 0x3F) | 0x80;
		return(3);
	}
	if (ch < 0x110000) {
		dest[0] = (ch >> 18) | 0xF0;
		dest[1] = ((ch >> 12) & 0x3F) | 0x80;
		dest[2] = ((ch >> 6) & 0x3F) | 0x80;
		dest[3] = (ch & 0x3F) | 0x80;
		return(4);
	}
	
        
	return(0);
}


/*************************************************************************/
/******************************  END UTF-8 FUNCTIONS  ********************/
/*************************************************************************/


sf8     val_equals_prop(NODE *curr_node, NODE *prop_node)
{
	sf8     prop_val;
	
	
	if (curr_node == prop_node)
		return(0.0);
	
	prop_val = prop_node->val;
	while (1) {
		curr_node = curr_node->next;
		if (curr_node == prop_node)
			return(-0.5);
		if (curr_node->val != prop_val)
			return(0.5);
	}
}


si4	write_MEF_file(FILE_PROCESSING_STRUCT *fps)
{
	si4	CRC_result;
	
	#ifdef _WIN32
		if (fps->full_file_name != NULL)
	        	slash_to_backslash(fps->full_file_name);
	#endif

	// clobber file if exists and is closed, create if non-existent
	if (fps->fp == NULL) {
		if (!(fps->directives.open_mode & FPS_GENERIC_WRITE_OPEN_MODE))
			fps->directives.open_mode = FPS_W_OPEN_MODE;
		fps_open(fps, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	}
	
	// offset times
	offset_universal_header_times(fps, RTO_OUTPUT_ACTION);
	switch (fps->file_type_code) {
		case TIME_SERIES_INDICES_FILE_TYPE_CODE:
			offset_time_series_index_times(fps, RTO_OUTPUT_ACTION);
			break;
		case VIDEO_INDICES_FILE_TYPE_CODE:
			offset_video_index_times(fps, RTO_OUTPUT_ACTION);
			break;
		case RECORD_INDICES_FILE_TYPE_CODE:
			offset_record_index_times(fps, RTO_OUTPUT_ACTION);
			break;
		default:
			break;
	}
	
	// encrypt
	switch (fps->file_type_code) {
		case TIME_SERIES_METADATA_FILE_TYPE_CODE:
		case VIDEO_METADATA_FILE_TYPE_CODE:
			encrypt_metadata(fps);
			break;
		case RECORD_DATA_FILE_TYPE_CODE:
			encrypt_records(fps);   // also does time offsets, for efficiency
			break;
		default:
			break;
	}
	
	// CRCs
	if (MEF_globals->CRC_mode & (CRC_CALCULATE | CRC_CALCULATE_ON_OUTPUT)) {
		if (fps->directives.io_bytes == FPS_FULL_FILE)  // if doing piecemeal writes, body CRC calculation should be done explicitly in the code
			fps->universal_header->body_CRC = CRC_calculate(fps->raw_data + UNIVERSAL_HEADER_BYTES, fps->raw_data_bytes - UNIVERSAL_HEADER_BYTES);
		fps->universal_header->header_CRC = CRC_calculate(fps->raw_data + CRC_BYTES, UNIVERSAL_HEADER_BYTES - CRC_BYTES);

	}
	if (MEF_globals->CRC_mode & (CRC_VALIDATE | CRC_VALIDATE_ON_OUTPUT)) {
		if (fps->directives.io_bytes == FPS_FULL_FILE) {
			CRC_result = CRC_validate(fps->raw_data + UNIVERSAL_HEADER_BYTES, fps->raw_data_bytes - UNIVERSAL_HEADER_BYTES, fps->universal_header->body_CRC);
			if (CRC_result == MEF_TRUE && MEF_globals->verbose == MEF_TRUE)
				UTF8_printf("Body CRC is valid in file \"%s\".\n", fps->full_file_name);
			else
				UTF8_fprintf(stderr, "Warning: body CRC is invalid in file \"%s\".\n", fps->full_file_name);
			CRC_result = CRC_validate(fps->raw_data + CRC_BYTES, UNIVERSAL_HEADER_BYTES - CRC_BYTES, fps->universal_header->header_CRC);
			if (CRC_result == MEF_TRUE && MEF_globals->verbose == MEF_TRUE)
				UTF8_printf("Header CRC is valid in file \"%s\".\n", fps->full_file_name);
			else
				UTF8_fprintf(stderr, "Warning: header CRC is invalid in file \"%s\".\n", fps->full_file_name);
		}
	}
	
	// write
	fps_write(fps, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
	
	// close
	if (fps->directives.close_file == MEF_TRUE) {
		if (fps->directives.io_bytes != FPS_FULL_FILE) {
			rewind(fps->fp);
			e_fwrite(fps->universal_header, sizeof(UNIVERSAL_HEADER), (size_t) 1, fps->fp, fps->full_file_name, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		}
		fps_close(fps);
	}
	
	// show
	if (MEF_globals->verbose == MEF_TRUE)
		show_file_processing_struct(fps);
	
	
	return(0);	
}


