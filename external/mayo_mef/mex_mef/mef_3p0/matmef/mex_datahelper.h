#ifndef MEX_DATAHELPER_
#define MEX_DATAHELPER_
/**
 * 	@file - headers
 * 	Functions to convert primitive c-datatypes to matlab primitive (1x1) arrays/matrices
 *	
 *  Copyright 2020, Max van den Boom
 *  Includes updates from Richard J. Cui - richard.cui@utoronto.ca (4 apr 2020)
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

mxArray *mxUint8ArrayByValue(ui1 *array, int num_bytes);

mxArray *mxUint8ByValue(ui1 value);
mxArray *mxInt8ByValue(si1 value);
mxArray *mxUint32ByValue(ui4 value);
mxArray *mxInt32ByValue(si4 value);
mxArray *mxUint64ByValue(ui8 value);
mxArray *mxInt64ByValue(si8 value);

mxArray *mxDoubleByValue(sf8 value);
mxArray *mxStringByUTF8Value(char *str);


#endif   // MEX_DATAHELPER_