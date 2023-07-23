
#include "mefrec.h"



/**********************************   show_record()   ********************************/

void	show_record(RECORD_HEADER *record_header, ui4 record_number, PASSWORD_DATA *pwd)
{
        si4	i, decryption_blocks;
	ui4	type_code, *type_string_int;
        ui1	*ui1_p, *decryption_key;
	si1	time_str[TIME_STRING_BYTES], hex_str[HEX_STRING_BYTES(CRC_BYTES)];
        
	
	// display record header fields
        printf("Record Number: %u\n", record_number);
	printf("---------------- Record Header - START ----------------\n");
	if (record_header->record_CRC == RECORD_HEADER_RECORD_CRC_NO_ENTRY) {
		printf("Record CRC: no entry\n");
	} else {
                generate_hex_string((ui1 *) &record_header->record_CRC, CRC_BYTES, hex_str);
		printf("Record CRC: %s\n", hex_str);
        }
	if (strlen(record_header->type_string)) {
		type_string_int = (ui4 *) record_header->type_string;
		type_code = *type_string_int;
                generate_hex_string((ui1 *) record_header->type_string, CRC_BYTES, hex_str);
		printf("Record Type String: %s (%s)\n", record_header->type_string, hex_str);
	} else {
		type_code = MEFREC_UnRc_TYPE_CODE;
		printf("Record Type String: no entry\n");
	}
	if (record_header->version_major == RECORD_HEADER_VERSION_MAJOR_NO_ENTRY || record_header->version_minor == RECORD_HEADER_VERSION_MINOR_NO_ENTRY) {
		if (record_header->version_major == RECORD_HEADER_VERSION_MAJOR_NO_ENTRY)
			printf("Record Version Major: no entry\n");
		else
			printf("Record Version Major: %u\n", record_header->version_major);
		if (record_header->version_minor == RECORD_HEADER_VERSION_MINOR_NO_ENTRY)
			printf("Record Version Minor: no entry\n");
		else
			printf("Record Version Minor: %u\n", record_header->version_minor);
	} else {
		printf("Record Version: %u.%u\n", record_header->version_major, record_header->version_minor);
		printf("Record Encryption: %d ", record_header->encryption);
	}
	if (record_header->encryption == NO_ENCRYPTION)
		printf("(none)\n");
	else if (record_header->encryption == LEVEL_1_ENCRYPTION)
		printf("(level 1, currently encrypted)\n");
	else if (record_header->encryption == LEVEL_2_ENCRYPTION)
		printf("(level 2, currently encrypted)\n");
	else if (record_header->encryption == LEVEL_1_ENCRYPTION_DECRYPTED)
		printf("(level 1, currently decrypted)\n");
	else if (record_header->encryption == LEVEL_2_ENCRYPTION_DECRYPTED)
		printf("(level 2, currently decrypted)\n");
	else
		printf("(unrecognized code)\n");
	if (record_header->bytes == RECORD_HEADER_BYTES_NO_ENTRY)
		printf("Record Bytes: no entry\n");
	else
		printf("Record Bytes: %u\n", record_header->bytes);

	if (record_header->time == RECORD_HEADER_TIME_NO_ENTRY)
		printf("Record Time: no entry\n");
	else {
		local_date_time_string(record_header->time, time_str);
		#ifdef _WIN32
			printf("Record Time: %lld (uUTC), %s (ascii, local)\n", ABS(record_header->time), time_str);
		#else
			printf("Record Time: %ld (uUTC), %s (ascii, local)\n", ABS(record_header->time), time_str);
		#endif
	}
	printf("----------------- Record Header - END -----------------\n\n");

        
	// decrypt record body if necesary & access level is sufficient
	printf("----------------- Record Body - START -----------------\n");
	if (record_header->encryption > NO_ENCRYPTION) {
                if (pwd->access_level >= record_header->encryption) {
                        if (record_header->encryption == LEVEL_1_ENCRYPTION)
                                decryption_key = pwd->level_1_encryption_key;
                        else
                                decryption_key = pwd->level_2_encryption_key;
                        decryption_blocks = record_header->bytes / ENCRYPTION_BLOCK_BYTES;
                        ui1_p = (ui1 *) record_header + RECORD_HEADER_BYTES;
                        for (i = decryption_blocks; i--;) {
                                AES_decrypt(ui1_p, ui1_p, NULL, decryption_key);
                                ui1_p += ENCRYPTION_BLOCK_BYTES;
                        }
                        record_header->encryption = -record_header->encryption;  // mark as currently decrypted
                        printf("                (record now decrypted)\n");
                } else {
                        printf("No access to this record\n");
                        printf("------------------ Record Body - END ------------------\n\n");
                        return;
                }
        }

	
	// pass the display off to custom functions - new records types should be added here (maintain alphabetical order of record types)
	switch (type_code) {
		case MEFREC_Note_TYPE_CODE:
			show_mefrec_Note_type(record_header);
			break;
		case MEFREC_EDFA_TYPE_CODE:
			show_mefrec_EDFA_type(record_header);
			break;
		case MEFREC_Seiz_TYPE_CODE:
			show_mefrec_Seiz_type(record_header);
			break;
		case MEFREC_SyLg_TYPE_CODE:
			show_mefrec_SyLg_type(record_header);
			break;
		case MEFREC_CSti_TYPE_CODE:
			show_mefrec_CSti_type(record_header);
			break;
		case MEFREC_UnRc_TYPE_CODE:
			printf("\"%s\" (0x%x) is an unrecognized record type\n", record_header->type_string, type_code);
			break;
		default:
			printf("\"%s\" (0x%x) is an unrecognized record type\n", record_header->type_string, type_code);
			break;
	}
	printf("------------------ Record Body - END ------------------\n\n");
        
        return;
}

/*************************************************************************************/



/**************************   check_record_alignments()   *****************************/

si1	check_record_structure_alignments(ui1 *bytes)
{
 	si4			return_value;
	si4			free_flag = MEF_FALSE;
        extern MEF_GLOBALS	*MEF_globals;

        
	// see if already checked
	if (MEF_globals->all_record_structures_aligned != MEF_UNKNOWN)
		return(MEF_globals->all_record_structures_aligned);
        
	return_value = MEF_TRUE;
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(LARGEST_RECORD_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	
	// check all structures - add new functions here
	if ((check_mefrec_Note_type_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_mefrec_EDFA_type_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_mefrec_LNTP_type_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_mefrec_Seiz_type_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_mefrec_SyLg_type_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	if ((check_mefrec_CSti_type_alignment(bytes)) == MEF_FALSE)
		return_value = MEF_FALSE;
	
        if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (return_value == MEF_TRUE) {
		MEF_globals->all_record_structures_aligned = MEF_TRUE;
		if (MEF_globals->verbose == MEF_TRUE)
			(void) printf("%s(): All Record structures are aligned\n", __FUNCTION__);
	} else {
		MEF_globals->all_record_structures_aligned = MEF_FALSE;
		if (MEF_globals->verbose == MEF_TRUE)
			(void) printf("%s(): One or more Record structures are not aligned\n", __FUNCTION__);
	}
	
	return(return_value);
}

/*************************************************************************************/



/********************************   Note: Note Record   ******************************/

void	show_mefrec_Note_type(RECORD_HEADER *record_header)
{
        si1	*Note;
  
        
        // Version 1.0
        if (record_header->version_major == 1 && record_header->version_minor == 0) {
                Note = (si1 *) record_header + MEFREC_Note_1_0_TEXT_OFFSET;
                UTF8_printf("Note text: %s\n", Note);
        }
        // Unrecognized record version
        else {
                printf("Unrecognized Note version\n");
        }
        
        return;
}

si4	check_mefrec_Note_type_alignment(ui1 *bytes)
{
        // no structures to check
        return(MEF_TRUE);
}

/*************************************************************************************/



/********************************   EDFA: EDF Annotation Record   ******************************/

void	show_mefrec_EDFA_type(RECORD_HEADER *record_header)
{
	MEFREC_EDFA_1_0	*edfa;
	si1		*annotation;
	
	
	// Version 1.0
	if (record_header->version_major == 1 && record_header->version_minor == 0) {
		edfa = (MEFREC_EDFA_1_0 *) ((ui1 *) record_header + MEFREC_EDFA_1_0_OFFSET);
		annotation = (si1 *) record_header + MEFREC_EDFA_1_0_ANNOTATION_OFFSET;
		#ifdef _WIN32
			UTF8_printf("Annotation: %s\nDuration %lld microseconds\n", annotation, edfa->duration);
		#else
			UTF8_printf("Annotation: %s\nDuration %ld microseconds\n", annotation, edfa->duration);
		#endif
	}
	// Unrecognized record version
	else {
		printf("Unrecognized Note version\n");
	}
	
	return;
}

si4	check_mefrec_EDFA_type_alignment(ui1 *bytes)
{
	MEFREC_EDFA_1_0		*edfa;
	si4			free_flag = MEF_FALSE;
	extern MEF_GLOBALS	*MEF_globals;


	// check overall size
	if (sizeof(MEFREC_EDFA_1_0) != MEFREC_EDFA_1_0_BYTES)
		goto MEFREC_EDFA_1_0_NOT_ALIGNED;
	
	// check fields
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(LARGEST_RECORD_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	edfa = (MEFREC_EDFA_1_0 *) (bytes + MEFREC_EDFA_1_0_OFFSET);
	if (&edfa->duration != (si8 *) (bytes + MEFREC_EDFA_1_0_DURATION_OFFSET))
		goto MEFREC_EDFA_1_0_NOT_ALIGNED;
	
	// aligned
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): MEFREC_EDFA_1_0 structure is aligned\n", __FUNCTION__);
	
	return(MEF_TRUE);
	
	// not aligned
	MEFREC_EDFA_1_0_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	(void) fprintf(stderr, "%c%s(): MEFREC_EDFA_1_0 structure is not aligned\n", 7, __FUNCTION__);
	
	return(MEF_FALSE);

}

/*************************************************************************************/



/*************************   LNTP: Line Noise Template Record   **********************/

void	show_mefrec_LNTP_type(RECORD_HEADER *record_header)
{
	MEFREC_LNTP_1_0	*lntp;
	si4		*template;
	si8		i;
	
	
	// Version 1.0
	if (record_header->version_major == 1 && record_header->version_minor == 0) {
		lntp = (MEFREC_LNTP_1_0 *) ((ui1 *) record_header + MEFREC_LNTP_1_0_OFFSET);
		template = (si4 *) record_header + MEFREC_LNTP_1_0_TEMPLATE_OFFSET;
		printf("Line Noise Template:\n");
		for (i = 0; i < lntp->length; ++i)
			printf("%d\n", template[i]);
	}
	// Unrecognized record version
	else {
		printf("Unrecognized LNTP version\n");
	}
	
	return;
}

si4	check_mefrec_LNTP_type_alignment(ui1 *bytes)
{
	MEFREC_LNTP_1_0		*lntp;
	si4			free_flag = MEF_FALSE;
	extern MEF_GLOBALS	*MEF_globals;
	
	
	// check overall size
	if (sizeof(MEFREC_LNTP_1_0) != MEFREC_LNTP_1_0_BYTES)
		goto MEFREC_LNTP_1_0_NOT_ALIGNED;
	
	// check fields
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(LARGEST_RECORD_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	lntp = (MEFREC_LNTP_1_0 *) (bytes + MEFREC_LNTP_1_0_OFFSET);
	if (&lntp->length != (si8 *) (bytes + MEFREC_LNTP_1_0_LENGTH_OFFSET))
		goto MEFREC_LNTP_1_0_NOT_ALIGNED;
	
	// aligned
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): MEFREC_LNTP_1_0 structure is aligned\n", __FUNCTION__);
	
	return(MEF_TRUE);
	
	// not aligned
	MEFREC_LNTP_1_0_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	(void) fprintf(stderr, "%c%s(): MEFREC_LNTP_1_0 structure is not aligned\n", 7, __FUNCTION__);
	
	return(MEF_FALSE);
	
}

/*************************************************************************************/



/*******************************   Seiz: Seizure Record   ****************************/

void	show_mefrec_Seiz_type(RECORD_HEADER *record_header)
{
        si4			i, mn1 = MEF_FALSE, mn2 = MEF_FALSE;
        MEFREC_Seiz_1_0		*seizure;
        MEFREC_Seiz_1_0_CHANNEL	*channels;
	si1			time_str[32];
        
        
        // Version 1.0
        if (record_header->version_major == 1 && record_header->version_minor == 0) {
                seizure = (MEFREC_Seiz_1_0 *) ((ui1 *) record_header + MEFREC_Seiz_1_0_OFFSET);
		local_date_time_string(seizure->earliest_onset, time_str);
		#ifdef _WIN32
			printf("Earliest Onset: %lld (uUTC), %s (ascii, local)\n", ABS(seizure->earliest_onset), time_str);
		#else
			printf("Earliest Onset: %ld (uUTC), %s (ascii, local)\n", ABS(seizure->earliest_onset), time_str);
		#endif
		local_date_time_string(seizure->latest_offset, time_str);
		#ifdef _WIN32
			printf("Latest Offset: %lld (uUTC), %s (ascii, local)\n", ABS(seizure->latest_offset), time_str);
			printf("Duration: %lld (microseconds)\n", seizure->duration);
		#else
			printf("Latest Offset: %ld (uUTC), %s (ascii, local)\n", ABS(seizure->latest_offset), time_str);
			printf("Duration: %ld (microseconds)\n", seizure->duration);
		#endif
		printf("Number of Channels: %d\n", seizure->number_of_channels);
		printf("Onset Code: %d ", seizure->onset_code);
                switch (seizure->onset_code) {
                        case MEFREC_Seiz_1_0_ONSET_NO_ENTRY:
                                printf("(no entry)\n");
                                break;
                        case MEFREC_Seiz_1_0_ONSET_UNKNOWN:
                                printf("(unknown)\n");
                                break;
                        case MEFREC_Seiz_1_0_ONSET_FOCAL:
                                printf("(focal)\n");
                                break;
                        case MEFREC_Seiz_1_0_ONSET_GENERALIZED:
                                printf("(generalized)\n");
                                break;
                        case MEFREC_Seiz_1_0_ONSET_PROPAGATED:
                                printf("(propagated)\n");
                                break;
                        case MEFREC_Seiz_1_0_ONSET_MIXED:
                                printf("(mixed)\n");
                                break;
                        default:
                                printf("(unrecognized code)\n");
                                break;
                }
                if (strlen(seizure->marker_name_1))
                        mn1 = MEF_TRUE;
                if (strlen(seizure->marker_name_2))
                        mn2 = MEF_TRUE;
                if (mn1 == MEF_TRUE && mn2 == MEF_TRUE)
                	UTF8_printf("Marker Name: %s %s\n", seizure->marker_name_1, seizure->marker_name_2);
                else if (mn1 == MEF_TRUE)
			UTF8_printf("Marker Name 1: %s\nMarker Name 2: no entry\n", seizure->marker_name_1);
                else if (mn2 == MEF_TRUE)
			UTF8_printf("Marker Name 1: no_entry\nMarker Name 2: %s\n", seizure->marker_name_2);
                else
                        printf("Marker Name: no_entry\n");
                if (strlen(seizure->annotation))
                        UTF8_printf("Annotation: %s\n", seizure->annotation);
                else
                        printf("Annotation: no entry\n");
                channels = (MEFREC_Seiz_1_0_CHANNEL *) ((ui1 *) record_header + MEFREC_Seiz_1_0_CHANNELS_OFFSET);
                for (i = 0; i < seizure->number_of_channels; ++i) {
                        if (strlen(channels[i].name))
                                UTF8_printf("Channel Name: %s\n", channels[i].name);
                        else
                                printf("Channel Name: no entry\n");
                        local_date_time_string(channels[i].onset, time_str);
                        #ifdef _WIN32
							printf("\tOnset: %lld (uUTC), %s (ascii, local)\n", ABS(channels[i].onset), time_str);
						#else
							printf("\tOnset: %ld (uUTC), %s (ascii, local)\n", ABS(channels[i].onset), time_str);
						#endif
                        local_date_time_string(channels[i].offset, time_str);
                        #ifdef _WIN32
							printf("\tOffset: %lld (uUTC), %s (ascii, local)\n", ABS(channels[i].offset), time_str);
						#else
							printf("\tOffset: %ld (uUTC), %s (ascii, local)\n", ABS(channels[i].offset), time_str);
						#endif
                }
        }
        // Unrecognized record version
        else {
                printf("Unrecognized Seiz version\n");
        }
         
        return;
}

si4	check_mefrec_Seiz_type_alignment(ui1 *bytes)
{
	MEFREC_Seiz_1_0		*seiz;
        MEFREC_Seiz_1_0_CHANNEL	*chan;
	si4			free_flag = MEF_FALSE;
        ui1			*chan_bytes;
        extern MEF_GLOBALS	*MEF_globals;
	
	
	// check overall sizes
	if (sizeof(MEFREC_Seiz_1_0) != MEFREC_Seiz_1_0_BYTES)
		goto MEFREC_Seiz_1_0_NOT_ALIGNED;
	if (sizeof(MEFREC_Seiz_1_0_CHANNEL) != MEFREC_Seiz_1_0_CHANNEL_BYTES)
		goto MEFREC_Seiz_1_0_NOT_ALIGNED;
	
	// check fields - base structure
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(LARGEST_RECORD_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	seiz = (MEFREC_Seiz_1_0 *) (bytes + MEFREC_Seiz_1_0_OFFSET);
	if (&seiz->earliest_onset != (si8 *) (bytes + MEFREC_Seiz_1_0_EARLIEST_ONSET_OFFSET))
		goto MEFREC_Seiz_1_0_NOT_ALIGNED;
	if (&seiz->latest_offset != (si8 *) (bytes + MEFREC_Seiz_1_0_LATEST_OFFSET_OFFSET))
		goto MEFREC_Seiz_1_0_NOT_ALIGNED;
	if (&seiz->duration != (si8 *) (bytes + MEFREC_Seiz_1_0_DURATION_OFFSET))
		goto MEFREC_Seiz_1_0_NOT_ALIGNED;
	if (&seiz->number_of_channels != (si4 *) (bytes + MEFREC_Seiz_1_0_NUMBER_OF_CHANNELS_OFFSET))
		goto MEFREC_Seiz_1_0_NOT_ALIGNED;
	if (&seiz->onset_code != (si4 *) (bytes + MEFREC_Seiz_1_0_ONSET_CODE_OFFSET))
		goto MEFREC_Seiz_1_0_NOT_ALIGNED;
	if (seiz->marker_name_1 != (si1 *) (bytes + MEFREC_Seiz_1_0_MARKER_NAME_1_OFFSET))
		goto MEFREC_Seiz_1_0_NOT_ALIGNED;
	if (seiz->marker_name_2 != (si1 *) (bytes + MEFREC_Seiz_1_0_MARKER_NAME_2_OFFSET))
		goto MEFREC_Seiz_1_0_NOT_ALIGNED;
	if (seiz->annotation != (si1 *) (bytes + MEFREC_Seiz_1_0_ANNOTATION_OFFSET))
		goto MEFREC_Seiz_1_0_NOT_ALIGNED;
	// check fields - channel structures
        chan_bytes = bytes + MEFREC_Seiz_1_0_CHANNELS_OFFSET;
	chan = (MEFREC_Seiz_1_0_CHANNEL *) chan_bytes;
	if (chan->name != (si1 *) (chan_bytes + MEFREC_Seiz_1_0_CHANNEL_NAME_OFFSET))
		goto MEFREC_Seiz_1_0_NOT_ALIGNED;
	if (&chan->onset != (si8 *) (chan_bytes + MEFREC_Seiz_1_0_CHANNEL_ONSET_OFFSET))
		goto MEFREC_Seiz_1_0_NOT_ALIGNED;
	if (&chan->offset != (si8 *) (chan_bytes + MEFREC_Seiz_1_0_CHANNEL_OFFSET_OFFSET))
		goto MEFREC_Seiz_1_0_NOT_ALIGNED;
        
	// aligned
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): MEFREC_Seiz_1_0 structure is aligned\n", __FUNCTION__);
	
        return(MEF_TRUE);
	
	// not aligned
	MEFREC_Seiz_1_0_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	(void) fprintf(stderr, "%c%s(): MEFREC_Seiz_1_0 structure is not aligned\n", 7, __FUNCTION__);
        
        return(MEF_FALSE);
}


/*************************************************************************************/



/*****************************   SyLg: System Log Record   ***************************/

void	show_mefrec_SyLg_type(RECORD_HEADER *record_header)
{
        si1	*log_entry;
        
        
        // Version 1.0
        if (record_header->version_major == 1 && record_header->version_minor == 0) {
                log_entry = (si1 *) record_header + MEFREC_SyLg_1_0_TEXT_OFFSET;
                UTF8_printf("System Log text: %s\n", log_entry);
        }
        // Unrecognized record version
        else {
                printf("Unrecognized SyLg version\n");
        }
        
        return;
}

si4	check_mefrec_SyLg_type_alignment(ui1 *bytes)
{
        // no structures to check	
        return(MEF_TRUE);
}

/*************************************************************************************/


//show_mefrec_CSti_type not tested!!!
/********************************   CSti: Cognitive stimulation   ******************************/

void	show_mefrec_CSti_type(RECORD_HEADER *record_header)
{
        MEFREC_CSti_1_0	*cog_stim;
        
        // Version 1.0
        if (record_header->version_major == 1 && record_header->version_minor == 0) {
                cog_stim = (MEFREC_CSti_1_0 *) ((ui1 *) record_header + MEFREC_CSti_1_0_OFFSET);
  
                if (strlen(cog_stim->task_type))
						UTF8_printf("Task type: %s\n", cog_stim->task_type);
                else
                    	printf("Task type: no_entry\n");

				#ifdef _WIN32
					printf("Duration: %lld (microseconds)\n", cog_stim->stimulus_duration);
				#else
					printf("Duration: %ld (microseconds)\n", cog_stim->stimulus_duration);
				#endif
				
                if (strlen(cog_stim->stimulus_type))
						UTF8_printf("Stimulation type: %s\n", cog_stim->stimulus_type);
                else
                    	printf("Stimulation type: no_entry\n");
                if (strlen(cog_stim->patient_response))
                    	UTF8_printf("Patient response: %s\n", cog_stim->patient_response);
                else
                    	printf("Patient response: no entry\n");
        }
        // Unrecognized record version
        else {
                printf("Unrecognized Note version\n");
        }
        
        return;
}

si4	check_mefrec_CSti_type_alignment(ui1 *bytes)
{
	MEFREC_CSti_1_0		*cog_stim;
	si4			free_flag = MEF_FALSE;
    extern MEF_GLOBALS	*MEF_globals;
	
	
	// check overall sizes
	if (sizeof(MEFREC_CSti_1_0) != MEFREC_CSti_1_0_BYTES)
		goto MEFREC_CSti_1_0_NOT_ALIGNED;
	
	// check fields - base structure
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(LARGEST_RECORD_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	cog_stim = (MEFREC_CSti_1_0 *) (bytes + MEFREC_CSti_1_0_OFFSET);
	if (cog_stim->task_type != (si1 *) (bytes + MEFREC_CSti_1_0_TASK_TYPE_OFFSET))
		goto MEFREC_CSti_1_0_NOT_ALIGNED;
	if (&cog_stim->stimulus_duration != (si8 *) (bytes + MEFREC_CSti_1_0_STIMULUS_DURATION_OFFSET))
		goto MEFREC_CSti_1_0_NOT_ALIGNED;
	if (cog_stim->stimulus_type != (si1 *) (bytes + MEFREC_CSti_1_0_STIMULUS_TYPE_OFFSET))
		goto MEFREC_CSti_1_0_NOT_ALIGNED;
	if (cog_stim->patient_response != (si1 *) (bytes + MEFREC_CSti_1_0_PATIENT_RESPONSE_OFFSET))
		goto MEFREC_CSti_1_0_NOT_ALIGNED;

	// aligned
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): MEFREC_CSti_1_0 structure is aligned\n", __FUNCTION__);
	
        return(MEF_TRUE);
	
	// not aligned
	MEFREC_CSti_1_0_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	(void) fprintf(stderr, "%c%s(): MEFREC_Csti_1_0 structure is not aligned\n", 7, __FUNCTION__);
        
        return(MEF_FALSE);
}


/*************************************************************************************/


//show_mefrec_ESti_type not tested!!!
/********************************   ESti: Electrical stimulation   ******************************/

void	show_mefrec_ESti_type(RECORD_HEADER *record_header)
{
        MEFREC_ESti_1_0	*el_stim;
        
        // Version 1.0
        if (record_header->version_major == 1 && record_header->version_minor == 0) {
                el_stim = (MEFREC_ESti_1_0 *) ((ui1 *) record_header + MEFREC_ESti_1_0_OFFSET);
  				

                printf("Amplitude: %f ", el_stim->amplitude);
                switch (el_stim->ampunit_code) {
                        case MEFREC_ESti_1_0_AMPUNIT_NO_ENTRY:
                                printf("(no entry)\n");
                                break;
                        case MEFREC_ESti_1_0_AMPUNIT_UNKNOWN:
                                printf("(unknown)\n");
                                break;
                        case MEFREC_ESti_1_0_AMPUNIT_MA:
                                printf("(mA)\n");
                                break;
                        case MEFREC_ESti_1_0_AMPUNIT_V:
                                printf("(V)\n");
                                break;
                        default:
                                printf("(unrecognized code)\n");
                                break;
                }
                printf("Frequency: %f (Hz)\n", el_stim->frequency);
				#ifdef _WIN32
					printf("Pulse width: %lld (microseconds)\n", el_stim->pulse_width);
				#else
					printf("Pulse width: %ld (microseconds)\n", el_stim->pulse_width);
				#endif
				printf("Operating mode: ");
				switch(el_stim->mode_code){
						case MEFREC_ESti_1_0_MODE_NO_ENTRY:
                                printf("(no entry)\n");
                                break;
                        case MEFREC_ESti_1_0_MODE_UNKNOWN:
                                printf("(unknown)\n");
                                break;
                        case MEFREC_ESti_1_0_MODE_CURRENT:
                                printf("(current)\n");
                                break;
                        case MEFREC_ESti_1_0_MODE_VOLTAGE:
                                printf("(voltage)\n");
                                break;
                        default:
                                printf("(unrecognized code)\n");
                                break;
				}
				
                if (strlen(el_stim->waveform))
						UTF8_printf("Waveform: %s\n", el_stim->waveform);
                else
                    	printf("Waveform: no_entry\n");
                if (strlen(el_stim->anode))
                    	UTF8_printf("Anode: %s\n", el_stim->anode);
                else
                    	printf("Anode: no entry\n");
                if (strlen(el_stim->catode))
                    	UTF8_printf("Catode: %s\n", el_stim->catode);
                else
                    	printf("Catode: no entry\n");
        }
        // Unrecognized record version
        else {
                printf("Unrecognized Note version\n");
        }
        
        return;
}


si4	check_mefrec_ESti_type_alignment(ui1 *bytes)
{
	MEFREC_ESti_1_0		*el_stim;
	si4			free_flag = MEF_FALSE;
    extern MEF_GLOBALS	*MEF_globals;
	
	
	// check overall sizes
	if (sizeof(MEFREC_ESti_1_0) != MEFREC_ESti_1_0_BYTES)
		goto MEFREC_ESti_1_0_NOT_ALIGNED;
	
	// check fields - base structure
	if (bytes == NULL) {
		bytes = (ui1 *) e_malloc(LARGEST_RECORD_BYTES, __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
		free_flag = MEF_TRUE;
	}
	el_stim = (MEFREC_ESti_1_0 *) (bytes + MEFREC_ESti_1_0_OFFSET);
	if (&el_stim->amplitude != (sf8 *) (bytes + MEFREC_ESti_1_0_AMPLITUDE_OFFSET))
		goto MEFREC_ESti_1_0_NOT_ALIGNED;
	if (&el_stim->frequency != (sf8 *) (bytes + MEFREC_ESti_1_0_FREQUENCY_OFFSET))
		goto MEFREC_ESti_1_0_NOT_ALIGNED;
	if (&el_stim->pulse_width != (si8 *) (bytes + MEFREC_ESti_1_0_PULSE_WIDTH_OFFSET))
		goto MEFREC_ESti_1_0_NOT_ALIGNED;
	if (&el_stim->ampunit_code != (si4 *) (bytes + MEFREC_ESti_1_0_AMPUNIT_CODE_OFFSET))
		goto MEFREC_ESti_1_0_NOT_ALIGNED;
	if (&el_stim->mode_code != (si4 *) (bytes + MEFREC_ESti_1_0_MODE_CODE_OFFSET))
		goto MEFREC_ESti_1_0_NOT_ALIGNED;
	if (el_stim->waveform != (si1 *) (bytes + MEFREC_ESti_1_0_WAVEFORM_OFFSET))
		goto MEFREC_ESti_1_0_NOT_ALIGNED;
	if (el_stim->anode != (si1 *) (bytes + MEFREC_ESti_1_0_ANODE_OFFSET))
		goto MEFREC_ESti_1_0_NOT_ALIGNED;
	if (el_stim->catode != (si1 *) (bytes + MEFREC_ESti_1_0_CATODE_OFFSET))
		goto MEFREC_ESti_1_0_NOT_ALIGNED;

	// aligned
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	if (MEF_globals->verbose == MEF_TRUE)
		(void) printf("%s(): MEFREC_ESti_1_0 structure is aligned\n", __FUNCTION__);
	
        return(MEF_TRUE);
	
	// not aligned
	MEFREC_ESti_1_0_NOT_ALIGNED:
	
	if (free_flag == MEF_TRUE)
		free(bytes);
	
	(void) fprintf(stderr, "%c%s(): MEFREC_Esti_1_0 structure is not aligned\n", 7, __FUNCTION__);
        
        return(MEF_FALSE);
}

/*************************************************************************************/