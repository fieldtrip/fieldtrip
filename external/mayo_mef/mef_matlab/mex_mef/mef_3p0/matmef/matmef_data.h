#ifndef MATMEF_DATA_
#define MATMEF_DATA_
/**
 * 	@file - headers
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
#include "mex.h"
#include "meflib/meflib/meflib.h"


// Range Types
#define RANGE_BY_SAMPLES	0
#define RANGE_BY_TIME		1

// 
// Functions
//

mxArray *read_channel_data_from_path(si1 *channel_path, si1 *password, bool range_type, si8 range_start, si8 range_end);
mxArray *read_channel_data_from_object(CHANNEL *channel, bool range_type, si8 range_start, si8 range_end);

si8 sample_for_uutc_c(si8 uutc, CHANNEL *channel);
si8 uutc_for_sample_c(si8 sample, CHANNEL *channel);
void memset_int(si4 *ptr, si4 value, size_t num);
si4 check_block_crc(ui1 *block_hdr_ptr, ui4 max_samps, ui1 *total_data_ptr, ui8 total_data_bytes);




#endif   // MATMEF_DATA_