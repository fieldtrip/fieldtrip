/**
 * 	@file
 * 	Functions to convert primitive c-datatypes to matlab primitive (1x1) arrays/matrices
 *	
 *  Copyright 2020, Max van den Boom
 *
 *  
 *  This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 *  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

// Copyright 2020 Richard J. Cui. Created: Sun 02/16/2020 10:34:49.777 PM
// $Revision: 0.3 $  $Date: Mon 02/17/2020  8:15:44.891 PM $
//
// 1026 Rocky Creek Dr NE
// Rochester, MN 55906, USA
//
// Email: richard.cui@utoronto.ca

#ifndef MEX_DATAHELPER_
#define MEX_DATAHELPER_
/**
 *
 * 	@file - header
 *
 *
 */
#include "mex.h"
#include "meflib.h"

mxArray *mxFillNumericArray(ui1 *array, int num_bytes);
mxArray *mxStringByUTF8Value(char *str);

mxArray *mxUint8ByValue(ui1 value);
mxArray *mxInt8ByValue(si1 value);
mxArray *mxUint32ByValue(ui4 value);
mxArray *mxInt32ByValue(si4 value);
mxArray *mxUint64ByValue(ui8 value);
mxArray *mxInt64ByValue(si8 value);
mxArray *mxDoubleByValue(sf8 value);

#endif   // MEX_DATAHELPER_
