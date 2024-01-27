/**
 * 	@file
 * 	MEF 3.0 Library Matlab Wrapper
 * 	Functions to load data from MEF3 datafiles
 *	
 *  Copyright 2020, Max van den Boom (Multimodal Neuroimaging Lab, Mayo Clinic, Rochester MN)
 *	Adapted from PyMef (by Jan Cimbalnik, Matt Stead, Ben Brinkmann, and Dan Crepeau)
 *  
 *	
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 *  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "matmef_data.h"
#include "mex.h"

#include "meflib/meflib/meflib.c"
#include "meflib/meflib/mefrec.c"


/**
 * 	Read the channel data from a channel filepath, given a range of data to read.
 *  The range is defined as a type (RANGE_BY_SAMPLES or RANGE_BY_TIME), a startpoint and an endpoint.
 * 	
 *
 * 	@param channel_path		The channel filepath
 * 	@param password			Password for the MEF3 datafiles (no password = NULL or empty string)
 *	@param range_type		Modality that is used to define the data-range to read [either 'time' or 'samples']
 *	@param range_start		Start-point for the reading of data (either as an epoch/unix timestamp or samplenumber; -1 for first)
 *	@param range_end		End-point to stop the of reading data (either as an epoch/unix timestamp or samplenumber; -1 for last)
 * 	@return					Pointer to a matlab double matrix object (mxArray) containing the data, or NULL on failure
 */
mxArray *read_channel_data_from_path(si1 *channel_path, si1 *password, bool range_type, si8 range_start, si8 range_end) {

	// check if the password is empty, correct to NULL if it is
	if (password != NULL && password[0] == '\0') {
		password = NULL;
	}
	
	// initialize MEF library
	(void) initialize_meflib();
	
	// read the channel metadata
	MEF_globals->behavior_on_fail = SUPPRESS_ERROR_OUTPUT;
	CHANNEL *channel = read_MEF_channel(NULL, channel_path, TIME_SERIES_CHANNEL_TYPE, password, NULL, MEF_FALSE, MEF_FALSE);
	
	// check the number of segments
	if (channel->number_of_segments == 0) {
		mexPrintf("Error: no segments in channel, most likely due to an invalid channel folder, exiting...\n"); 
        return NULL;
	}
	
	// check if the data is encrypted and/or the correctness of password
	if (channel->metadata.section_1->section_2_encryption > 0 || channel->metadata.section_1->section_2_encryption > 0) {
		if (password == NULL)
			mexPrintf("Error: data is encrypted, but no password is given, exiting...\n"); 
		else
			mexPrintf("Error: wrong password for encrypted data, exiting...\n"); 
		return NULL;
	}
	
	// check if the channel is indeed of a time-series channel
	if (channel->channel_type != TIME_SERIES_CHANNEL_TYPE) {
		mexPrintf("Error: not a time series channel, exiting...\n"); 
		return NULL;
	}
	
	// read the data by the channel object
	mxArray *samples_read = read_channel_data_from_object(channel, range_type, range_start, range_end);
			
	// free the channel object memory
	if (channel->number_of_segments > 0)	channel->segments[0].metadata_fps->directives.free_password_data = MEF_TRUE;
	free_channel(channel, MEF_TRUE);

	// return the number of samples that were read
	return samples_read;
	
}

/**
 * 	Read the channel data based on a channel object (pointer) and a range of data to read.
 *  The range is defined as a type (RANGE_BY_SAMPLES or RANGE_BY_TIME), a startpoint and an endpoint.
 * 	
 *	Note: this function does not free the memory of the given channel object (that is up to the function's caller)
 *
 * 	@param channel			Pointer to the MEF channel object
 *	@param range_type		Modality that is used to define the data-range to read [either 'time' or 'samples']
 *	@param range_start		Start-point for the reading of data (either as an epoch/unix timestamp or samplenumber; -1 for first)
 *	@param range_end		End-point to stop the of reading data (either as an epoch/unix timestamp or samplenumber; -1 for last)	
 * 	@return					Pointer to a matlab double matrix object (mxArray) containing the data, or NULL on failure
 */
mxArray *read_channel_data_from_object(CHANNEL *channel, bool range_type, si8 range_start, si8 range_end) {
	ui8     i, j;
	ui8		num_blocks;
	ui8		num_block_in_segment;
	
	// check if the channel is indeed of a time-series channel
	if (channel->channel_type != TIME_SERIES_CHANNEL_TYPE) {
		mexPrintf("Error: not a time series channel, exiting...\n"); 
        return NULL;
    }
	
	// check the number of segments
	if (channel->number_of_segments == 0) {
		mexPrintf("Error: no segments in channel, exiting...\n"); 
        return NULL;
	}
	
	// set the default ranges for the samples and time to all
	si8 start_samp = 0;
    si8 start_time = channel->earliest_start_time;
    si8 end_samp = channel->metadata.time_series_section_2->number_of_samples;
    si8 end_time = channel->latest_end_time;
	
	// update the ranges if available (> -1)
	if (range_start > -1)	start_samp 	= start_time 	= range_start;
	if (range_end > -1)		end_samp 	= end_time 		= range_end;
	
	// check if valid data range
    if (range_type == RANGE_BY_TIME && start_time >= end_time) {
		mexPrintf("Error: start time later than end time, exiting...\n");
        return NULL;
    }
    if (range_type == RANGE_BY_SAMPLES && start_samp >= end_samp) {
        mexPrintf("Error: start sample larger than end sample, exiting...\n");
        return NULL;
    }
	
    // fire warnings if start or stop or both are out of file
    if (range_type == RANGE_BY_TIME) {
		
        if (((start_time < channel->earliest_start_time) & (end_time < channel->earliest_start_time)) |
            ((start_time > channel->latest_end_time) & (end_time > channel->latest_end_time))) {
            mexPrintf("Error: start and stop times are out of file.\n");
            return NULL;
        }
		
        if (end_time > channel->latest_end_time)			mexPrintf("Warning: stop uutc later than latest end time. Will insert NaNs\n");
        if (start_time < channel->earliest_start_time)		mexPrintf("Warning: start uutc earlier than earliest start time. Will insert NaNs\n");
		
    } else {
		
        if (((start_samp < 0) & (end_samp < 0)) |
            ((start_samp > channel->metadata.time_series_section_2->number_of_samples) & (end_samp > channel->metadata.time_series_section_2->number_of_samples))) {
            mexPrintf("Error: start and stop samples are out of file\n");
            return NULL;
        }
        if (end_samp > channel->metadata.time_series_section_2->number_of_samples) {
            mexPrintf("Error: stop sample larger than number of samples. Setting end sample to number of samples in channel\n");
			return NULL;
        }
        if (start_samp < 0) {
            mexPrintf("Error: start sample smaller than 0. Setting start sample to 0\n");
			return NULL;
        }
		
    }
	
    // determine the number of samples
    ui8 num_samps = 0;
    if (range_type == RANGE_BY_TIME)
        num_samps = (ui4)((((end_time - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency) + 0.5);
    else
        num_samps = (ui4) (end_samp - start_samp);
	
	// check if the range has no samples
	if (num_samps == 0) {
		
		// message
		mexPrintf("Warning: a range of 0 samples was given, returning empty array\n");
		
		// return an empty array
		return mxCreateDoubleMatrix(1, 1, mxREAL);
		
	}
	
    // iterate through segments, looking for data that matches our criteria
    ui4 n_segments = (ui4) channel->number_of_segments;
    ui4 start_segment = -1;
	ui4 end_segment = -1;
    
	// convert the range in samples to time or vise versa
    if (range_type == RANGE_BY_TIME) {
        start_samp = sample_for_uutc_c(start_time, channel);
        end_samp = sample_for_uutc_c(end_time, channel);
    } else {
        start_time = uutc_for_sample_c(start_samp, channel);
        end_time = uutc_for_sample_c(end_samp, channel);
    }
	

    // find start and stop segments by uutc time
    // find start segment by finding first segment whose ending is past the start time.
    // then find stop segment by using the previous segment of the (first segment whose start is past the end time)
    for (i = 0; i < n_segments; ++i) {
        
        if (range_type == RANGE_BY_TIME) {
			
            si8 segment_start_time = channel->segments[i].time_series_data_fps->universal_header->start_time;
            si8 segment_end_time   = channel->segments[i].time_series_data_fps->universal_header->end_time;
            remove_recording_time_offset( &segment_start_time);
            remove_recording_time_offset( &segment_end_time);
			
            if ((segment_end_time >= start_time) && (start_segment == -1)){
                start_segment = i;
                end_segment = i;
            }
            if ((end_segment != -1) && (segment_start_time <= end_time))
                end_segment = i;
			
        } else {
				
            si8 segment_start_sample = channel->segments[i].metadata_fps->metadata.time_series_section_2->start_sample;
            si8 segment_end_sample   = channel->segments[i].metadata_fps->metadata.time_series_section_2->start_sample +
										channel->segments[i].metadata_fps->metadata.time_series_section_2->number_of_samples;
			
            if ((start_samp >= segment_start_sample) && (start_samp <= segment_end_sample))
                start_segment = i;
            if ((end_samp >= segment_start_sample) && (end_samp <= segment_end_sample))
                end_segment = i;
            
        }
		
    }

	// check if both the start- and endsegment were found
	if (start_segment == -1 || end_segment == -1) {

		// message
		mexPrintf("Error: unable to find the start segment (%i) or end segment (%i), existing...\n", start_segment, end_segment);
		return NULL;
		
	}

    // find start block in start segment
    ui8 start_idx = 0;
	ui8 end_idx = 0;
    for (j = 1; j < channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks; j++) {
		
        si8 block_start_time = channel->segments[start_segment].time_series_indices_fps->time_series_indices[j].start_time;
        remove_recording_time_offset( &block_start_time);
        
        if (block_start_time > start_time) {
            start_idx = j - 1;
            break;
        }
        // starting point is in last block in segment
        start_idx = j;
    }
    
    // find stop block in stop segment
    for (j = 1; j < channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks; j++) {
        
        si8 block_start_time = channel->segments[end_segment].time_series_indices_fps->time_series_indices[j].start_time;
        remove_recording_time_offset( &block_start_time);
        
        if (block_start_time > end_time) {
            end_idx = j - 1;
            break;
        }
        // ending point is in last block in segment
        end_idx = j;
    }
    
    // find total_samps and total_data_bytes, so we can allocate buffers
    si8 total_samps = 0;
    ui8 total_data_bytes = 0;
    
    // check if the data is in one segment or multiple
    if (start_segment == end_segment) {
		// normal case - everything is in one segment
		
        if (end_idx < (ui8) (channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks - 1)) {
            total_samps += channel->segments[start_segment].time_series_indices_fps->time_series_indices[end_idx+1].start_sample -
            channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample;
            total_data_bytes += channel->segments[start_segment].time_series_indices_fps->time_series_indices[end_idx+1].file_offset -
            channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        } else {
            // case where end_idx is last block in segment
            total_samps += channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_samples -
            channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample;
            total_data_bytes += channel->segments[start_segment].time_series_data_fps->file_length -
            channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        }
        num_blocks = end_idx - start_idx + 1;
        
    } else {
		// spans across segments
		
		// start with first segment
        num_block_in_segment = (ui8) channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks;
        total_samps += channel->segments[start_segment].metadata_fps->metadata.time_series_section_2->number_of_samples -
        channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample;
        total_data_bytes +=  channel->segments[start_segment].time_series_data_fps->file_length -
        channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        num_blocks = num_block_in_segment - start_idx;

        if (channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset < 1024){
            mexPrintf("Error: Invalid index file offset, exiting....\n");
            return NULL;
        }
        
        // this loop will only run if there are segments in between the start and stop segments
        for (i = (start_segment + 1); i <= (end_segment - 1); i++) {
            num_block_in_segment = (ui8) channel->segments[i].metadata_fps->metadata.time_series_section_2->number_of_blocks;
            total_samps += channel->segments[i].metadata_fps->metadata.time_series_section_2->number_of_samples;
            total_data_bytes += channel->segments[i].time_series_data_fps->file_length -
            channel->segments[i].time_series_indices_fps->time_series_indices[0].file_offset;
            num_blocks += num_block_in_segment;

            if (channel->segments[i].time_series_indices_fps->time_series_indices[0].file_offset < 1024){
                mexPrintf("Error: Invalid index file offset, exiting....\n");
                return NULL;
            }
        }
        
        // then last segment
        num_block_in_segment = (ui8) channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks;
        if (end_idx < (ui8) (channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks - 1)) {
            total_samps += channel->segments[end_segment].time_series_indices_fps->time_series_indices[end_idx+1].start_sample -
            channel->segments[end_segment].time_series_indices_fps->time_series_indices[0].start_sample;
            total_data_bytes += channel->segments[end_segment].time_series_indices_fps->time_series_indices[end_idx+1].file_offset -
            channel->segments[end_segment].time_series_indices_fps->time_series_indices[0].file_offset;
            num_blocks += end_idx + 1;
        } else {
            // case where end_idx is last block in segment
            total_samps += channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_samples -
            channel->segments[end_segment].time_series_indices_fps->time_series_indices[0].start_sample;
            total_data_bytes += channel->segments[end_segment].time_series_data_fps->file_length -
            channel->segments[end_segment].time_series_indices_fps->time_series_indices[0].file_offset;
            num_blocks += end_idx + 1;
        }

        if (channel->segments[end_segment].time_series_indices_fps->time_series_indices[end_idx].file_offset < 1024){
			mexPrintf("Error: Invalid index file offset, exiting....\n");
            return NULL;
        }
		
	}
	
    // allocate a buffer for the compressed data
    ui1 *compressed_data_buffer = (ui1*) malloc((size_t) total_data_bytes);
	if (compressed_data_buffer == NULL) {
		mexPrintf("Error: could not allocated enough memory for the compressed data, exiting....\n");
		return NULL;
	}
    ui1 *cdp = compressed_data_buffer;
	
	// allocate the samples buffer
    si4 *decomp_data = (si4*) malloc((size_t) (num_samps * sizeof(si4)));
	if (decomp_data == NULL) {
		free (compressed_data_buffer);
		mexPrintf("Error: could not allocated enough memory for the sample buffer, exiting....\n");
		return NULL;
	}
	
	// initialize the entire sample buffer to nan
	memset_int(decomp_data, RED_NAN, num_samps);
	
    // read in RED data
    if (start_segment == end_segment) {
		// normal case - everything is in one segment
		
        if (channel->segments[start_segment].time_series_data_fps->fp == NULL){
            channel->segments[start_segment].time_series_data_fps->fp = fopen(channel->segments[start_segment].time_series_data_fps->full_file_name, "rb");
			#ifdef _WIN32
				channel->segments[start_segment].time_series_data_fps->fd = _fileno(channel->segments[start_segment].time_series_data_fps->fp);
			#else
				channel->segments[start_segment].time_series_data_fps->fd = fileno(channel->segments[start_segment].time_series_data_fps->fp);
			#endif
        }
        FILE *fp = channel->segments[start_segment].time_series_data_fps->fp;
        
        #ifdef _WIN32
            _fseeki64(fp, channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset, SEEK_SET);
        #else
            fseek(fp, channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset, SEEK_SET);
        #endif
        ui8 n_read = fread(cdp, sizeof(si1), (size_t) total_data_bytes, fp);
        if (n_read != total_data_bytes) {
			mexPrintf("Warning: read in fewer than expected bytes from data file in segment %d.\n", start_segment);
		}
        if (channel->segments[start_segment].time_series_data_fps->directives.close_file == MEF_TRUE)
            fps_close(channel->segments[start_segment].time_series_data_fps);
		
    } else {
		// spans across segments
		
		// start with first segment
        if (channel->segments[start_segment].time_series_data_fps->fp == NULL){
            channel->segments[start_segment].time_series_data_fps->fp = fopen(channel->segments[start_segment].time_series_data_fps->full_file_name, "rb");
			#ifdef _WIN32
				channel->segments[start_segment].time_series_data_fps->fd = _fileno(channel->segments[start_segment].time_series_data_fps->fp);
			#else
				channel->segments[start_segment].time_series_data_fps->fd = fileno(channel->segments[start_segment].time_series_data_fps->fp);
			#endif
        }
        FILE *fp = channel->segments[start_segment].time_series_data_fps->fp;
        #ifdef _WIN32
             _fseeki64(fp, channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset, SEEK_SET);
        #else
            fseek(fp, channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset, SEEK_SET);
        #endif
        ui8 bytes_to_read = channel->segments[start_segment].time_series_data_fps->file_length -
        channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].file_offset;
        ui8 n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
        if (n_read != bytes_to_read) {
			mexPrintf("Warning: read in fewer than expected bytes from data file in segment %d.\n", start_segment);
        }
        cdp += n_read;
        if (channel->segments[start_segment].time_series_data_fps->directives.close_file == MEF_TRUE)
            fps_close(channel->segments[start_segment].time_series_data_fps);
        
        // this loop will only run if there are segments in between the start and stop segments
        for (i = (start_segment + 1); i <= (end_segment - 1); i++) {
            if (channel->segments[i].time_series_data_fps->fp == NULL){
                channel->segments[i].time_series_data_fps->fp = fopen(channel->segments[i].time_series_data_fps->full_file_name, "rb");
				#ifdef _WIN32
					channel->segments[i].time_series_data_fps->fd = _fileno(channel->segments[i].time_series_data_fps->fp);
				#else
					channel->segments[i].time_series_data_fps->fd = fileno(channel->segments[i].time_series_data_fps->fp);
				#endif
            }
            fp = channel->segments[i].time_series_data_fps->fp;
            fseek(fp, UNIVERSAL_HEADER_BYTES, SEEK_SET);
            bytes_to_read = channel->segments[i].time_series_data_fps->file_length - 
            channel->segments[i].time_series_indices_fps->time_series_indices[0].file_offset;
            n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
            if (n_read != bytes_to_read) {
				mexPrintf("Warning: read in fewer than expected bytes from data file in segment %d.\n", i);
            }
            cdp += n_read;
            if (channel->segments[i].time_series_data_fps->directives.close_file == MEF_TRUE)
                fps_close(channel->segments[i].time_series_data_fps);
        }
        
        // then last segment
        if (channel->segments[end_segment].time_series_data_fps->fp == NULL){
            channel->segments[end_segment].time_series_data_fps->fp = fopen(channel->segments[end_segment].time_series_data_fps->full_file_name, "rb");
			#ifdef _WIN32
				channel->segments[end_segment].time_series_data_fps->fd = _fileno(channel->segments[end_segment].time_series_data_fps->fp);
			#else
				channel->segments[end_segment].time_series_data_fps->fd = fileno(channel->segments[end_segment].time_series_data_fps->fp);
			#endif
        }
        num_block_in_segment = channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks;
        if (end_idx < (ui8) (channel->segments[end_segment].metadata_fps->metadata.time_series_section_2->number_of_blocks - 1)) {
            fp = channel->segments[end_segment].time_series_data_fps->fp;
            fseek(fp, UNIVERSAL_HEADER_BYTES, SEEK_SET);
            bytes_to_read = channel->segments[end_segment].time_series_indices_fps->time_series_indices[end_idx+1].file_offset -
            channel->segments[end_segment].time_series_indices_fps->time_series_indices[0].file_offset;
            n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
            if (n_read != bytes_to_read) {
				mexPrintf("Warning: read in fewer than expected bytes from data file in segment %d.\n", end_segment);
            }
            cdp += n_read;
        } else {
            // case where end_idx is last block in segment
            fp = channel->segments[end_segment].time_series_data_fps->fp;
            fseek(fp, UNIVERSAL_HEADER_BYTES, SEEK_SET);
            bytes_to_read = channel->segments[end_segment].time_series_data_fps->file_length -
            channel->segments[end_segment].time_series_indices_fps->time_series_indices[0].file_offset;
            n_read = fread(cdp, sizeof(si1), (size_t) bytes_to_read, fp);
            if (n_read != bytes_to_read) {
				mexPrintf("Warning: read in fewer than expected bytes from data file in segment %d.\n", end_segment);
            }
            cdp += n_read;
        }

        if (channel->segments[end_segment].time_series_data_fps->directives.close_file == MEF_TRUE)
            fps_close(channel->segments[end_segment].time_series_data_fps);


	}
	
    // set up RED processing struct
    cdp = compressed_data_buffer;
	ui4 max_samps = channel->metadata.time_series_section_2->maximum_block_samples;
	
    // create RED processing struct
    RED_PROCESSING_STRUCT *rps = (RED_PROCESSING_STRUCT *) calloc((size_t) 1, sizeof(RED_PROCESSING_STRUCT));
    rps->compression.mode = RED_DECOMPRESSION;
    rps->decompressed_ptr = rps->decompressed_data = decomp_data;
    rps->difference_buffer = (si1 *) e_calloc((size_t) RED_MAX_DIFFERENCE_BYTES(max_samps), sizeof(ui1), __FUNCTION__, __LINE__, USE_GLOBAL_BEHAVIOR);
    
    // reset the pointer back to the start of the array
    cdp = compressed_data_buffer;
    
	// 
	si8 sample_counter = 0;
	si8 offset_into_output_buffer;
	si8 block_start_time_offset;

	//
	// decode the first block
	// 
    si4 *temp_data_buf = (int *) malloc((max_samps * 1.1) * sizeof(si4));
    if (temp_data_buf == NULL) {
        free (compressed_data_buffer);
        free (decomp_data);
        mexPrintf("Error: could not allocated enough memory for the block buffer, exiting....\n");
        return NULL;
    }
    rps->decompressed_ptr = rps->decompressed_data = temp_data_buf;
    rps->compressed_data = cdp;
    rps->block_header = (RED_BLOCK_HEADER *) rps->compressed_data;
    if (!check_block_crc((ui1 *)(rps->block_header), max_samps, compressed_data_buffer, total_data_bytes)) {
		// incorrect crc
		
		// message
		mexPrintf("Error: RED block %lu has 0 bytes, or CRC failed, data likely corrupt...", start_idx);

		//
        free (compressed_data_buffer);
        free (decomp_data);
        free (temp_data_buf);
        return NULL;
		
    }

	// 
	RED_decode(rps);
	cdp += rps->block_header->block_bytes;
	
	// 
	if (range_type == RANGE_BY_TIME)
		offset_into_output_buffer = (si4) ((((rps->block_header->start_time - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency) + 0.5);
	else
		offset_into_output_buffer = (si4) channel->segments[start_segment].time_series_indices_fps->time_series_indices[start_idx].start_sample - start_samp;
	
	// copy requested samples from first block to output buffer
	// TODO: this loop could be optimized
	for (i=0;i<rps->block_header->number_of_samples;i++) {
		if (offset_into_output_buffer < 0) {
			offset_into_output_buffer++;
			continue;
		}
		
		if ((ui4) offset_into_output_buffer >= num_samps)
			break;
		
		*(decomp_data + offset_into_output_buffer) = temp_data_buf[i];
		offset_into_output_buffer++;
	}

	// 
	sample_counter = offset_into_output_buffer;
	
	
	//
    // decode blocks in between the first and the last
	//
	for (i=1;i<num_blocks-1;i++) {
		
		// 
        rps->compressed_data = cdp;
        rps->block_header = (RED_BLOCK_HEADER *) rps->compressed_data;
        // check that block fits fully within output array
        // this should be true, but it's possible a stray block exists out-of-order, or with a bad timestamp
        
        // we need to manually remove offset, since we are using the time value of the block before decoding the block
        // (normally the offset is removed during the decoding process)
        if ((rps->block_header->block_bytes == 0) || !check_block_crc((ui1*)(rps->block_header), max_samps, compressed_data_buffer, total_data_bytes)) {
			// incorrect crc
						
			// message
			mexPrintf("Error: RED block %lu has 0 bytes, or CRC failed, data likely corrupt...", start_idx + i);

			//
			free (compressed_data_buffer);
			free (decomp_data);
			free (temp_data_buf);
			return NULL;
			
        }
		
		if (range_type == RANGE_BY_TIME) {
			block_start_time_offset = rps->block_header->start_time;
			remove_recording_time_offset(&block_start_time_offset);
			
			// The next two checks see if the block contains out-of-bounds samples.
			// In that case, skip the block and move on
			if (block_start_time_offset < start_time) {
				cdp += rps->block_header->block_bytes;
				continue;
			}
			if (block_start_time_offset + ((rps->block_header->number_of_samples / channel->metadata.time_series_section_2->sampling_frequency) * 1e6) >= end_time) {
				// Comment this out for now, it creates a strange boundary condition
				// cdp += rps->block_header->block_bytes;
				continue;
			}
			
			rps->decompressed_ptr = rps->decompressed_data = decomp_data + (int)((((block_start_time_offset - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency) + 0.5);
			
		} else {
			
			// buffer overflow check
			if ((sample_counter + rps->block_header->number_of_samples) > num_samps) {
	
				// message
				// TODO: better fix for buffer overflow, should not happen
				mexPrintf("Error: buffer overflow prevented, this should be fixed in the code");

				//
				free (compressed_data_buffer);
				free (decomp_data);
				free (temp_data_buf);
				return NULL;				
			
			}
			
			// 
			rps->decompressed_ptr = rps->decompressed_data = decomp_data + sample_counter;
		}
		
		// 
		RED_decode(rps);
		sample_counter += rps->block_header->number_of_samples;

		//
        cdp += rps->block_header->block_bytes;
		
    }	
	
	// 
	// decode last block to temp array
	// 	
    if (num_blocks > 1) {
		
		//
        rps->compressed_data = cdp;
        rps->block_header = (RED_BLOCK_HEADER *) rps->compressed_data;
        rps->decompressed_ptr = rps->decompressed_data = temp_data_buf;
        if (!check_block_crc((ui1*)(rps->block_header), max_samps, compressed_data_buffer, total_data_bytes)) {
			// incorrect crc
			
			// message
			mexPrintf("Error: RED block %lu has 0 bytes, or CRC failed, data likely corrupt...", start_idx + i);

			//
			free (compressed_data_buffer);
			free (decomp_data);
			free (temp_data_buf);
			return NULL;
			
        }
		
		// 
        RED_decode(rps);
        
		// 
        if (range_type == RANGE_BY_TIME)
            offset_into_output_buffer = (int)((((rps->block_header->start_time - start_time) / 1000000.0) * channel->metadata.time_series_section_2->sampling_frequency) + 0.5);
        else
            offset_into_output_buffer = sample_counter;
        
        // copy requested samples from last block to output buffer
        for (i=0;i<rps->block_header->number_of_samples;i++) {
            if (offset_into_output_buffer < 0) {
                offset_into_output_buffer++;
                continue;
            }
            
            if ((ui4) offset_into_output_buffer >= num_samps)
                break;
            
            *(decomp_data + offset_into_output_buffer) = temp_data_buf[i];
            offset_into_output_buffer++;
			
        }
		
    }
    
	// 
    // Decompressed data are integers (si4) and represent the "real" data; Integers technically do not have a NaN value as it exists for float datatypes.
    // Internally a si4 emulated NaN value ('RED_NAN') is used, however this value is not standard for Matlab (or Python)
    //
    // When range is indicated in time, then gaps/discontinuities in the data need to be filled with NaNs. Therefore, we return
    // the data as doubles. Integer to float cannot by byte-copied, therefore it needs to be cast per element, which is slow.
    //
    // TODO: could be speed up by directly casting or copying the values in the above loops to the matlab array; or
    //       when the range is indicated in samples, no nans exist, so we could allocate a matlab int array at the
    //       start and read data directly into that (above loops)
    //
    
	// allocate matlab double array
    // using doubles so we can use NaN values for discontinuities (when range is indicated in time)
	mxArray *mat_array = mxCreateDoubleMatrix(1, num_samps, mxREAL);
	mxDouble *ptr_mat_array = mxGetPr(mat_array);
    
    // copy/cast the data to the matlab array
	mxDouble mxNaN = mxGetNaN();
    for (i = 0; i < num_samps; i++) {
        if (*(decomp_data + i) == RED_NAN)
			*(ptr_mat_array + i) = mxNaN;
        else
            *(ptr_mat_array + i) = (sf8) *(decomp_data + i);
    }
    
    // free the memory holding the compressed data
    free (decomp_data);
    free (temp_data_buf);
    free (compressed_data_buffer);
    free (rps->difference_buffer);
    free (rps);
	
	// return the data
	return mat_array;
	
}


si8 sample_for_uutc_c(si8 uutc, CHANNEL *channel) {
    ui8 i, j, sample;
    sf8 native_samp_freq;
    ui8 prev_sample_number;
    si8 prev_time, seg_start_sample;
    
    native_samp_freq = channel->metadata.time_series_section_2->sampling_frequency;
    prev_sample_number = channel->segments[0].metadata_fps->metadata.time_series_section_2->start_sample;
    prev_time = channel->segments[0].time_series_indices_fps->time_series_indices[0].start_time;
    
    for (j = 0; j < channel->number_of_segments; j++) {
        seg_start_sample = channel->segments[j].metadata_fps->metadata.time_series_section_2->start_sample;
        for (i = 0; i < channel->segments[j].metadata_fps->metadata.time_series_section_2->number_of_blocks; ++i) {
            if (channel->segments[j].time_series_indices_fps->time_series_indices[i].start_time > uutc)
                goto done;
            prev_sample_number = channel->segments[j].time_series_indices_fps->time_series_indices[i].start_sample + seg_start_sample;
            prev_time = channel->segments[j].time_series_indices_fps->time_series_indices[i].start_time;
        }
    }
    
done:
    sample = prev_sample_number + (ui8) (((((sf8) (uutc - prev_time)) / 1000000.0) * native_samp_freq) + 0.5);
    
    return(sample);
}
			   

si8 uutc_for_sample_c(si8 sample, CHANNEL *channel) {
    ui8 i, j, uutc;
    sf8 native_samp_freq;
    ui8 prev_sample_number;
    si8 prev_time, seg_start_sample;
    
    
    native_samp_freq = channel->metadata.time_series_section_2->sampling_frequency;
    prev_sample_number = channel->segments[0].metadata_fps->metadata.time_series_section_2->start_sample;
    prev_time = channel->segments[0].time_series_indices_fps->time_series_indices[0].start_time;
    
    for (j = 0; j < channel->number_of_segments; j++) {
        seg_start_sample = channel->segments[j].metadata_fps->metadata.time_series_section_2->start_sample;
        for (i = 0; i < channel->segments[j].metadata_fps->metadata.time_series_section_2->number_of_blocks; ++i){
            if (channel->segments[j].time_series_indices_fps->time_series_indices[i].start_sample + seg_start_sample > sample)
                goto done;
            prev_sample_number = channel->segments[j].time_series_indices_fps->time_series_indices[i].start_sample + seg_start_sample;
            prev_time = channel->segments[j].time_series_indices_fps->time_series_indices[i].start_time;
        }
    }
    
done:
    uutc = prev_time + (ui8) ((((sf8) (sample - prev_sample_number) / native_samp_freq) * 1000000.0) + 0.5);
    
    return(uutc);
}

void memset_int(si4 *ptr, si4 value, size_t num) {
    si4 *temp_ptr;
    
    if (num < 1)
        return;
    
	si4 *limit = ptr + num;
    for (temp_ptr = ptr; temp_ptr < limit; ++temp_ptr)
		*temp_ptr = value;
	
}

si4 check_block_crc(ui1 *block_hdr_ptr, ui4 max_samps, ui1 *total_data_ptr, ui8 total_data_bytes) {
    ui8 offset_into_data, remaining_buf_size;
    si1 CRC_valid;
    RED_BLOCK_HEADER *block_header;
    
    offset_into_data = block_hdr_ptr - total_data_ptr;
    remaining_buf_size = total_data_bytes - offset_into_data;
    
    // check if remaining buffer at least contains the RED block header
    if (remaining_buf_size < RED_BLOCK_HEADER_BYTES)
        return 0;
    
    block_header = (RED_BLOCK_HEADER*) block_hdr_ptr;
    
    // check if entire block, based on size specified in header, can possibly fit in the remaining buffer
    if (block_header->block_bytes > remaining_buf_size)
        return 0;
    
    // check if size specified in header is absurdly large
    if (block_header->block_bytes > RED_MAX_COMPRESSED_BYTES(max_samps, 1))
        return 0;
    
    // at this point we know we have enough data to actually run the CRC calculation, so do it
    CRC_valid = CRC_validate((ui1*) block_header + CRC_BYTES, block_header->block_bytes - CRC_BYTES, block_header->block_CRC);
    
    // return output of CRC heck
    if (CRC_valid == MEF_TRUE)
        return 1;
    else
        return 0;
	
}


