/**
 * 	@file 
 * 	MEF 3.0 Library Matlab Wrapper
 * 	Read the MEF3 data from a time-series channel
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
#include <ctype.h>
#include "mex.h"
#include "matmef_data.h"


/**
 * Main entry point for 'read_mef_ts_data'
 *
 * @param channelPath	path (absolute or relative) to the MEF3 channel folder
 * @param password		Password to the MEF3 data; Pass empty string/variable if not encrypted
 * @param rangeType		Modality that is used to define the data-range to read [either 'time' or 'samples']
 * @param rangeStart	Start-point for the reading of data (either as an epoch/unix timestamp or samplenumber; -1 for first)
 * @param rangeEnd		End-point to stop the of reading data (either as an epoch/unix timestamp or samplenumber; -1 for last)
 * @return				A vector of doubles holding the channel data
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	
	//
	// channel path
	// 
	
	si1 channel_path[MEF_FULL_FILE_NAME_BYTES];
	
	// check the channel path input argument
    if(nrhs < 1) {
        mexErrMsgIdAndTxt( "MATLAB:read_mef_ts_data:noChannelPathArg", "channelPath input argument not set");
	} else {
		if(!mxIsChar(prhs[0])) {
			mexErrMsgIdAndTxt( "MATLAB:read_mef_ts_data:invalidChannelPathArg", "channelPath input argument invalid, should string (array of characters)");
		}
		if(mxIsEmpty(prhs[0])) {
			mexErrMsgIdAndTxt( "MATLAB:read_mef_ts_data:invalidChannelPathArg", "channelPath input argument invalid, argument is empty");
		}
	}
	
	// set the channel path
	char *mat_channel_path = mxArrayToString(prhs[0]);
	MEF_strncpy(channel_path, mat_channel_path, MEF_FULL_FILE_NAME_BYTES);
	

	// 
	// password (optional)
	// 
	
	si1 *password = NULL;
	si1 password_arr[PASSWORD_BYTES] = {0};
	
	// check if a password input argument is given
    if (nrhs > 1) {
		
		// note: if the password passed to any of the meflib read function is an empty string, than 
		//		 the 'process_password_data' function in 'meflib.c' will crash everything, so make
		// 		 sure it is either NULL or a string with at least one character
		
		// check if the password input argument is not empty
		if (!mxIsEmpty(prhs[1])) {
		
			// check the password input argument
			if (!mxIsChar(prhs[1])) {
				mexErrMsgIdAndTxt( "MATLAB:read_mef_session_metadata:invalidPasswordArg", "password input argument invalid, should string (array of characters)");
			}
			
			// TODO: really need a MEF3 dataset (which cannot be read without a password) to check
			// set the password
			//char *mat_password = mxArrayToUTF8String(prhs[1]);
			char *mat_password = mxArrayToString(prhs[1]);
			password = strcpy(password_arr, mat_password);
	
		}
		
	}
	
	
	//
	// range
	//
	
	bool range_type = RANGE_BY_SAMPLES;
	si8 range_start = -1;
	si8 range_end = -1;
	
	// check if a range=type input argument is given
    if (nrhs > 2) {
		
		// check valid range type
		if (!mxIsChar(prhs[2])) {
			mexErrMsgIdAndTxt( "MATLAB:read_mef_session_metadata:invalidRangeTypeArg", "range-type input argument invalid, should string (array of characters)");
		}
		char *mat_range_type = mxArrayToString(prhs[2]);
		for(int i = 0; mat_range_type[i]; i++)	mat_range_type[i] = tolower(mat_range_type[i]);
		if (strcmp(mat_range_type, "time") != 0 && strcmp(mat_range_type, "samples") != 0) {
			mexErrMsgIdAndTxt( "MATLAB:read_mef_session_metadata:invalidRangeTypeArg", "range-type input argument invalid, allowed values are 'time' or 'samples'");
		}
		
		// set the range type
		if (strcmp(mat_range_type, "time") == 0)
			range_type = RANGE_BY_TIME;
		
		
		// check if a range-start input argument is given
		if (nrhs > 3) {
			
			// check if single numeric
			if (!mxIsNumeric(prhs[3]) || mxGetNumberOfElements(prhs[3]) > 1) {
				mexErrMsgIdAndTxt( "MATLAB:read_mef_session_metadata:invalidRangeStartArg", "range-start input argument invalid, should be a single value numeric (either -1 or >=0)");
			}
			
			// set the range-start value
			range_start = mxGetScalar(prhs[3]);
            
            // check if -1 or positive value
            if (range_start != -1 && range_start < 0) {
				mexErrMsgIdAndTxt( "MATLAB:read_mef_session_metadata:invalidRangeStartArg", "range-start input argument invalid, should be a single value numeric (either -1 or >=0)");
			}
            
		}
		
		// check if a range-end input argument is given
		if (nrhs > 4) {
			
			// check if single numeric
			if (!mxIsNumeric(prhs[4]) || mxGetNumberOfElements(prhs[4]) > 1) {
				mexErrMsgIdAndTxt( "MATLAB:read_mef_session_metadata:invalidRangeEndArg", "range-end input argument invalid, should be a single value numeric (either -1 or >=0)");
			}
			
			// set the range-end value
			range_end = mxGetScalar(prhs[4]);
            
            // check if -1 or positive value
            if (range_end != -1 && range_end < 0) {
				mexErrMsgIdAndTxt( "MATLAB:read_mef_session_metadata:invalidRangeEndArg", "range-end input argument invalid, should be a single value numeric (either -1 or >=0)");
			}
            
		}
		
	}
	
	// 
	// read the data
	// 
	mxArray *data = read_channel_data_from_path(channel_path, password, range_type, range_start, range_end);
	
	// check for error
	if (data == NULL)	mexErrMsgTxt("Error while reading channel data");

	// check if output is expected
	if (nlhs > 0) {
		
		// set the data as output
		plhs[0] = data;
		
	}
	
	// succesfull return from call
	return;
	
}
