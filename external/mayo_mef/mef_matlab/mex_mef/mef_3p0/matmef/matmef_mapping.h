#ifndef MATMEF_MAPPING_
#define MATMEF_MAPPING_
/**
 * 	@file - headers
 * 	MEF 3.0 Library Matlab Wrapper
 * 	Functions to map MEF3 objects and structures to Matlab structures
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

// Copyright 2020 Richard J. Cui. Created: Sun 02/16/2020 10:34:49.777 PM
// $Revision: 0.2 $  $Date: Mon 02/17/2020  7:38:26.233 PM $
//
// 1026 Rocky Creek Dr NE
// Rochester, MN 55906, USA
//
// Email: richard.cui@utoronto.ca

#include "mex.h"
#include "meflib.h"

// structures

// function declaration
void map_mef3_segment_universal_header_tostruct(
        UNIVERSAL_HEADER *universal_header, // universal header of 1st segment
        mxArray *mat_universal_header,
        int mat_index // index of structure matrix
        );
mxArray *map_mef3_segment_universal_header(UNIVERSAL_HEADER *universal_header);
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


#endif   // MATMEF_MAPPING_