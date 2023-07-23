#ifndef MATMEF_MAPPING_
#define MATMEF_MAPPING_
/**
 * 	@file - headers
 * 	MEF 3.0 Library Matlab Wrapper
 * 	Functions to map MEF3 objects and structures to Matlab structures
 *	
 *  Copyright 2020, Max van den Boom (Multimodal Neuroimaging Lab, Mayo Clinic, Rochester MN)
 *	Adapted from PyMef (by Jan Cimbalnik, Matt Stead, Ben Brinkmann, and Dan Crepeau)
 *  Includes updates from Richard J. Cui - richard.cui@utoronto.ca (4 apr 2020)
 *	
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 *  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "mex.h"
#include "meflib/meflib/meflib.h"


void map_mef3_segment_tostruct(SEGMENT *segment, si1 map_indices_flag, mxArray *mat_segment, int mat_index);
mxArray *map_mef3_segment(SEGMENT *segment, si1 map_indices_flag);
void map_mef3_channel_tostruct(CHANNEL *channel, si1 map_indices_flag, mxArray *mat_channel, int mat_index);
mxArray *map_mef3_channel(CHANNEL *channel, si1 map_indices_flag);
mxArray *map_mef3_session(SESSION *session, si1 map_indices_flag);

mxArray *map_mef3_md1(METADATA_SECTION_1 *md1);
mxArray *map_mef3_tmd2(TIME_SERIES_METADATA_SECTION_2 *tmd2);
mxArray *map_mef3_vmd2(VIDEO_METADATA_SECTION_2 *vmd2);
mxArray *map_mef3_md3(METADATA_SECTION_3 *md3);

mxArray *map_mef3_ti(TIME_SERIES_INDEX *ti, si8 number_of_entries);
mxArray *map_mef3_vi(VIDEO_INDEX *vi, si8 number_of_entries);
mxArray *map_mef3_records(FILE_PROCESSING_STRUCT *ri_fps, FILE_PROCESSING_STRUCT *rd_fps);

mxArray *map_mef3_csti(RECORD_HEADER *rh);
	
mxArray *map_mef3_uh(UNIVERSAL_HEADER *universal_header);

#endif   // MATMEF_MAPPING_