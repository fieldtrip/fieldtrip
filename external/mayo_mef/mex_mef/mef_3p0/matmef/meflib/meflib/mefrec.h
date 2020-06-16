
#ifndef MEFREC_IN
#define MEFREC_IN


#include "meflib.h"


/*************************************************************************************/
/***************************   General Record Constants   ****************************/
/*************************************************************************************/

// This needs to be kept up to date.  If it exceeds METADATA_FILE_BYTES then check_all_alignments() should be updated.
#define LARGEST_RECORD_BYTES	(RECORD_HEADER_BYTES + MEFREC_Seiz_1_0_BYTES + MEFREC_Seiz_1_0_CHANNEL_BYTES)  // Seiz Record type



/*************************************************************************************/
/********************************   Record Types   ***********************************/
/*************************************************************************************/


// Use CamelCase for 4-letter type string (case sensitive types)
// Translate to type code one byte at a time with ascii table
// (note that bytes, and therefore letters, are reversed in little-endian versions)
// #defines should follow the format in examples below
// Function names should follow the format in the prototypes below


/******************************   EDFA: EDF Annotation   *****************************/

// Constants
#define MEFREC_EDFA_TYPE_STRING		"EDFA"			// ascii[4]
#define MEFREC_EDFA_TYPE_CODE		(ui4) 0x41464445	// ui4 (little endian)
// #define MEFREC_EDFA_TYPE_CODE	(ui4) 0x45444641	// ui4 (big endian)

// Version 1.0
#define MEFREC_EDFA_1_0_OFFSET			(RECORD_HEADER_BYTES + 0)
#define MEFREC_EDFA_1_0_BYTES			8
#define MEFREC_EDFA_1_0_DURATION_OFFSET		(RECORD_HEADER_BYTES + 0)	// si8
#define MEFREC_EDFA_1_0_ANNOTATION_OFFSET	(RECORD_HEADER_BYTES + 8)	// si1

// Structures
typedef struct {
	si8	duration;  // microseconds (single element structure not really necessary, but done for consistency).  Annotation follows structure - aribitrary length array of si1s
} MEFREC_EDFA_1_0;

// Prototypes
void	show_mefrec_EDFA_type(RECORD_HEADER *record_header);
si4	check_mefrec_EDFA_type_alignment(ui1 *bytes);

/*************************************************************************************/



/****************************   LNTP: Line Noise Template   **************************/

// Constants
#define MEFREC_LNTP_TYPE_STRING		"LNTP"			// ascii[4]
#define MEFREC_LNTP_TYPE_CODE		(ui4) 0x50544e4c	// ui4 (little endian)
// #define MEFREC_LNTP_TYPE_CODE	(ui4) 0x4c4e5450	// ui4 (big endian)

// Version 1.0
#define MEFREC_LNTP_1_0_OFFSET			(RECORD_HEADER_BYTES + 0)
#define MEFREC_LNTP_1_0_BYTES			8
#define MEFREC_LNTP_1_0_LENGTH_OFFSET		(RECORD_HEADER_BYTES + 0)	// si8
#define MEFREC_LNTP_1_0_TEMPLATE_OFFSET		(RECORD_HEADER_BYTES + 8)	// si4

// Structures
typedef struct {
	si8	length;  // (single element structure not really necessary, but done for consistency).  Template follows structure - array of si4s.
} MEFREC_LNTP_1_0;

// Prototypes
void	show_mefrec_LNTP_type(RECORD_HEADER *record_header);
si4	check_mefrec_LNTP_type_alignment(ui1 *bytes);

/*************************************************************************************/



/********************************   Note: Note Record   ******************************/

// Constants
#define MEFREC_Note_TYPE_STRING		"Note"			// ascii[4]
#define MEFREC_Note_TYPE_CODE		(ui4) 0x65746f4e	// ui4 (little endian)
// #define MEFREC_Note_TYPE_CODE	(ui4) 0x4e6f7465	// ui4 (big endian)

// Version 1.0
#define MEFREC_Note_1_0_TEXT_OFFSET	RECORD_HEADER_BYTES

// Structures
// (none)

// Prototypes
void	show_mefrec_Note_type(RECORD_HEADER *record_header);
si4	check_mefrec_Note_type_alignment(ui1 *bytes);

/*************************************************************************************/



/*******************************   Seiz: Seizure Record   ****************************/

// Constants
#define MEFREC_Seiz_TYPE_STRING				"Seiz"							// ascii[4]
#define MEFREC_Seiz_TYPE_CODE				(ui4) 0x7a696553					// ui4 (little endian)
// #define MEFREC_Seiz_TYPE_CODE			(ui4) 0x5365697a					// ui4 (big endian)

// Version 1.0
// MEFREC_Seiz_1_0 offsets below apply to base address of record
#define MEFREC_Seiz_1_0_OFFSET				RECORD_HEADER_BYTES					// version 1.0
#define MEFREC_Seiz_1_0_BYTES				1312
#define MEFREC_Seiz_1_0_EARLIEST_ONSET_OFFSET		(RECORD_HEADER_BYTES + 0)				// si8
#define MEFREC_Seiz_1_0_LATEST_OFFSET_OFFSET		(RECORD_HEADER_BYTES + 8)				// si8
#define MEFREC_Seiz_1_0_DURATION_OFFSET			(RECORD_HEADER_BYTES + 16)				// si8
#define MEFREC_Seiz_1_0_NUMBER_OF_CHANNELS_OFFSET	(RECORD_HEADER_BYTES + 24)				// si4
#define MEFREC_Seiz_1_0_ONSET_CODE_OFFSET		(RECORD_HEADER_BYTES + 28)				// si4
#define MEFREC_Seiz_1_0_MARKER_NAME_1_OFFSET		(RECORD_HEADER_BYTES + 32)				// utf8[31]
#define MEFREC_Seiz_1_0_MARKER_NAME_BYTES		128
#define MEFREC_Seiz_1_0_MARKER_NAME_2_OFFSET		(RECORD_HEADER_BYTES + 160)				// utf8[31]
#define MEFREC_Seiz_1_0_ANNOTATION_OFFSET		(RECORD_HEADER_BYTES + 288)				// utf8[255]
#define MEFREC_Seiz_1_0_ANNOTATION_BYTES		1024
#define MEFREC_Seiz_1_0_CHANNELS_OFFSET			(MEFREC_Seiz_1_0_OFFSET + MEFREC_Seiz_1_0_BYTES)
// MEFREC_Seiz_1_0_CHANNEL offsets below apply to base address of channel
#define MEFREC_Seiz_1_0_CHANNEL_BYTES			272
#define MEFREC_Seiz_1_0_CHANNEL_NAME_OFFSET		0	// utf8[31]
#define MEFREC_Seiz_1_0_CHANNEL_ONSET_OFFSET		256	// si8
#define MEFREC_Seiz_1_0_CHANNEL_ONSET_NO_ENTRY		-1
#define MEFREC_Seiz_1_0_CHANNEL_OFFSET_OFFSET		264	// si8
#define MEFREC_Seiz_1_0_CHANNEL_OFFSET_NO_ENTRY		-1
// Onset Codes
#define MEFREC_Seiz_1_0_ONSET_NO_ENTRY		-1
#define MEFREC_Seiz_1_0_ONSET_UNKNOWN		0
#define MEFREC_Seiz_1_0_ONSET_FOCAL		1
#define MEFREC_Seiz_1_0_ONSET_GENERALIZED	2
#define MEFREC_Seiz_1_0_ONSET_PROPAGATED	3
#define MEFREC_Seiz_1_0_ONSET_MIXED		4

// Structures
typedef struct {
        si8	earliest_onset;						// uutc
        si8	latest_offset;						// uutc
        si8	duration;						// microseconds
        si4	number_of_channels;
        si4	onset_code;
        si1	marker_name_1[MEFREC_Seiz_1_0_MARKER_NAME_BYTES];	// utf8[31]
        si1	marker_name_2[MEFREC_Seiz_1_0_MARKER_NAME_BYTES];	// utf8[31]
        si1	annotation[MEFREC_Seiz_1_0_ANNOTATION_BYTES];		// utf8[255]
} MEFREC_Seiz_1_0;

typedef struct {
 	si1	name[MEF_BASE_FILE_NAME_BYTES];	// no extension
       	si8	onset;				// uutc
        si8	offset;				// uutc
} MEFREC_Seiz_1_0_CHANNEL;

// Prototypes
void	show_mefrec_Seiz_type(RECORD_HEADER *record_header);
si4	check_mefrec_Seiz_type_alignment(ui1 *bytes);

/*************************************************************************************/



/*****************************   SyLg: System Log Record   ***************************/

// Constants
#define MEFREC_SyLg_TYPE_STRING		"SyLg"			// ascii[4]
#define MEFREC_SyLg_TYPE_CODE		(ui4) 0x674c7953	// ui4 (little endian)
// #define MEFREC_SyLg_TYPE_CODE	(ui4) 0x53794c67	// ui4 (big endian)

// Version 1.0
#define MEFREC_SyLg_1_0_TEXT_OFFSET	RECORD_HEADER_BYTES

// Structures
// (none)

// Prototype
void	show_mefrec_SyLg_type(RECORD_HEADER *record_header);
si4	check_mefrec_SyLg_type_alignment(ui1 *bytes);

/*************************************************************************************/


/********************************   Csti: Cognitive stimulation   ******************************/

// Constants
#define MEFREC_CSti_TYPE_STRING     "CSti"              // ascii[4]
#define MEFREC_CSti_TYPE_CODE       (ui4) 0x69745343    // ui4 (little endian)
// #define MEFREC_ESti_TYPE_CODE    (ui4) 0x43537469    // ui4 (big endian)

// Version 1.0
#define MEFREC_CSti_1_0_OFFSET          RECORD_HEADER_BYTES
#define MEFREC_CSti_1_0_BYTES           200 // Watch out here - need to pad when writing.
#define MEFREC_CSti_1_0_TASK_TYPE_OFFSET              (RECORD_HEADER_BYTES + 0)        // utf8[15]
#define MEFREC_CSti_1_0_TASK_TYPE_BYTES             64  
#define MEFREC_CSti_1_0_STIMULUS_DURATION_OFFSET   (RECORD_HEADER_BYTES + 64)   // si8
#define MEFREC_CSti_1_0_STIMULUS_TYPE_OFFSET       (RECORD_HEADER_BYTES + 72)     // utf8[15]
#define MEFREC_CSti_1_0_STIMULUS_TYPE_BYTES      64
#define MEFREC_CSti_1_0_PATIENT_RESPONSE_OFFSET       (RECORD_HEADER_BYTES + 136)     // utf8[15]
#define MEFREC_CSti_1_0_PATIENT_RESPONSE_BYTES      64

// Structures
typedef struct {
    si1 task_type[MEFREC_CSti_1_0_TASK_TYPE_BYTES];
    si8 stimulus_duration;
    si1 stimulus_type[MEFREC_CSti_1_0_STIMULUS_TYPE_BYTES];
    si1 patient_response[MEFREC_CSti_1_0_PATIENT_RESPONSE_BYTES];
} MEFREC_CSti_1_0;

// Prototypes
void    show_mefrec_CSti_type(RECORD_HEADER *record_header);
si4 check_mefrec_CSti_type_alignment(ui1 *bytes);

/*************************************************************************************/


/********************************   Esti: Electric stimulation   ******************************/

// Constants
#define MEFREC_ESti_TYPE_STRING     "ESti"              // ascii[4]
#define MEFREC_ESti_TYPE_CODE       (ui4) 0x69745345    // ui4 (little endian)
// #define MEFREC_ESti_TYPE_CODE    (ui4) 0x45537469    // ui4 (big endian)

// Version 1.0
#define MEFREC_ESti_1_0_OFFSET          RECORD_HEADER_BYTES
#define MEFREC_ESti_1_0_BYTES           416
#define MEFREC_ESti_1_0_AMPLITUDE_OFFSET    (RECORD_HEADER_BYTES + 0)   //sf8
#define MEFREC_ESti_1_0_FREQUENCY_OFFSET    (RECORD_HEADER_BYTES + 8)   //sf8
#define MEFREC_ESti_1_0_PULSE_WIDTH_OFFSET  (RECORD_HEADER_BYTES + 16)  //si8
#define MEFREC_ESti_1_0_AMPUNIT_CODE_OFFSET (RECORD_HEADER_BYTES + 24)  //si4
#define MEFREC_ESti_1_0_MODE_CODE_OFFSET    (RECORD_HEADER_BYTES + 28)  //si4
#define MEFREC_ESti_1_0_WAVEFORM_OFFSET     (RECORD_HEADER_BYTES + 32)  //utf8[31]
#define MEFREC_ESti_1_0_WAVEFORM_BYTES      128
#define MEFREC_ESti_1_0_ANODE_OFFSET        (RECORD_HEADER_BYTES + 160) //utf8[31]
#define MEFREC_ESti_1_0_ANODE_BYTES         128
#define MEFREC_ESti_1_0_CATODE_OFFSET       (RECORD_HEADER_BYTES + 288) //utf8[31]
#define MEFREC_ESti_1_0_CATODE_BYTES        128

// Unit codes
#define MEFREC_ESti_1_0_AMPUNIT_NO_ENTRY      -1
#define MEFREC_ESti_1_0_AMPUNIT_UNKNOWN       0
#define MEFREC_ESti_1_0_AMPUNIT_MA     1
#define MEFREC_ESti_1_0_AMPUNIT_V   2
// Mode codes
#define MEFREC_ESti_1_0_MODE_NO_ENTRY      -1
#define MEFREC_ESti_1_0_MODE_UNKNOWN       0
#define MEFREC_ESti_1_0_MODE_CURRENT     1
#define MEFREC_ESti_1_0_MODE_VOLTAGE   2

// Structures
typedef struct {
    sf8 amplitude;
    sf8 frequency;
    si8 pulse_width;
    si4 ampunit_code;
    si4 mode_code;
    si1 waveform[MEFREC_ESti_1_0_WAVEFORM_BYTES];
    si1 anode[MEFREC_ESti_1_0_ANODE_BYTES];
    si1 catode[MEFREC_ESti_1_0_CATODE_BYTES];
} MEFREC_ESti_1_0;

// Prototypes
void    show_mefrec_ESti_type(RECORD_HEADER *record_header);
si4 check_mefrec_ESti_type_alignment(ui1 *bytes);

/*************************************************************************************/



/***************************   UnRc: Unrecognized Record   ***************************/

// Constants
#define MEFREC_UnRc_TYPE_STRING		"UnRc"			// ascii[4]
#define MEFREC_UnRc_TYPE_CODE		(ui4) 0x63526e55	// ui4 (little endian)
// #define MEFREC_UnRc_TYPE_CODE	(ui4) 0x556e5263	// ui4 (big endian)

/*************************************************************************************/


#endif // MEFREC_IN

