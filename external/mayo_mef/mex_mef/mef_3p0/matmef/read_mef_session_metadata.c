/**
 * 	@file 
 * 	MEF 3.0 Library Matlab Wrapper
 * 	Read a MEF3 folder and retrieve the session, channel(s), segment(s) and record(s) metadata
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
#include "matmef_mapping.h"

#include "meflib/meflib/meflib.c"
#include "meflib/meflib/mefrec.c"


/**
 * Main entry point for 'read_mef_session_metadata'
 *
 * @param sessionPath	Path (absolute or relative) to the MEF3 session folder
 * @param password		Password to the MEF3 data; Pass empty string/variable if not encrypted
 * @param mapIndices	Flag whether indices should be mapped [0 or 1; default is 0]
 * @return				Structure containing session metadata, channels metadata, segments metadata and records
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	
	//
	// session path
	// 
	
	// check the session path input argument
    if (nrhs < 1) {
        mexErrMsgIdAndTxt( "MATLAB:read_mef_session_metadata:noSessionPathArg", "sessionPath input argument not set");
	} else {
		if (!mxIsChar(prhs[0])) {
			mexErrMsgIdAndTxt( "MATLAB:read_mef_session_metadata:invalidSessionPathArg", "sessionPath input argument invalid, should string (array of characters)");
		}
		if (mxIsEmpty(prhs[0])) {
			mexErrMsgIdAndTxt( "MATLAB:read_mef_session_metadata:invalidSessionPathArg", "sessionPath input argument invalid, argument is empty");
		}
	}
	
	// set the session path
	si1 session_path[MEF_FULL_FILE_NAME_BYTES];
	char *mat_session_path = mxArrayToString(prhs[0]);
	MEF_strncpy(session_path, mat_session_path, MEF_FULL_FILE_NAME_BYTES);
	
	
	// 
	// password (optional)
	// 
	
	si1 *password = NULL;
	si1 password_arr[PASSWORD_BYTES] = {0};
	
	// check if a password input argument is given
    if (nrhs > 1) {

		// note: if the password passed to any of the meflib read function is an empty string, than 
		//		 the 'process_password_data' function in 'meflib.c' will crash everything, so make
		// 		 sure it is either NULL or a string

		// check if the password input argument is not empty
		if (!mxIsEmpty(prhs[1])) {
				
			// check the password input argument data type
			if (!mxIsChar(prhs[1])) {
				mexErrMsgIdAndTxt( "MATLAB:read_mef_session_metadata:invalidPasswordArg", "password input argument invalid, should string (array of characters)");
			}			
			
			// TODO: really need a MEF3 dataset (which cannot be read without a password) to check
			//char *mat_password = mxArrayToUTF8String(prhs[1]);
			char *mat_password = mxArrayToString(prhs[1]);
			password = strcpy(password_arr, mat_password);
	
		}
		
	}
	
	
	// 
	// map indices (optional)
	// 
	
	// Read indices flag
    si1 map_indices_flag = 0;	
	
	// check if a map indices input argument is given
    if (nrhs > 2) {
		
		// check if not empty
		if (!mxIsEmpty(prhs[2])) {

			// check if single numeric or logical
			if ((!mxIsNumeric(prhs[2]) && !mxIsLogical(prhs[2])) || mxGetNumberOfElements(prhs[2]) > 1) {
				mexErrMsgIdAndTxt( "MATLAB:read_mef_session_metadata:invalidMapIndicesArg", "mapIndices input argument invalid, should be a single value logical or numeric");
			}
			
			// retrieve the map indices flag value
			int mat_map_indices_flag = mxGetScalar(prhs[2]);
			
			// check the value
			if (mat_map_indices_flag != 0 && mat_map_indices_flag != 1) {
				mexErrMsgIdAndTxt( "MATLAB:read_mef_session_metadata:invalidMapIndicesArg", "mapIndices input argument invalid, allowed values are 0, false, 1 or true");
			}
			
			// set the flag
			map_indices_flag = mat_map_indices_flag;
			
		}
		
	}
	
	
	//
	// read session metadata
	//
	
    // initialize MEF library
	initialize_meflib();

	// read the session metadata
	MEF_globals->behavior_on_fail = SUPPRESS_ERROR_OUTPUT;
	SESSION *session = read_MEF_session(	NULL, 					// allocate new session object
											session_path, 			// session filepath
											password, 				// password
											NULL, 					// empty password
											MEF_FALSE, 				// do not read time series data
											MEF_TRUE				// read record data
										);
	MEF_globals->behavior_on_fail = EXIT_ON_FAIL;
	
	// check for error
	if (session == NULL)	mexErrMsgTxt("Error while reading session metadata");
	
	// check if the data is encrypted and/or the correctness of password
	if (session->time_series_metadata.section_1 != NULL) {
		if (session->time_series_metadata.section_1->section_2_encryption > 0 || session->time_series_metadata.section_1->section_2_encryption > 0) {
			free_session(session, MEF_TRUE);
			if (password == NULL)
				mexErrMsgTxt("Error: data is encrypted, but no password is given, exiting...\n");
			else
				mexErrMsgTxt("Error: wrong password for encrypted data, exiting...\n");
			
		}
	}
	if (session->video_metadata.section_1 != NULL) {
		if (session->video_metadata.section_1->section_2_encryption > 0 || session->video_metadata.section_1->section_2_encryption > 0) {
			free_session(session, MEF_TRUE);
			if (password == NULL)
				mexErrMsgTxt("Error: data is encrypted, but no password is given, exiting...\n");
			else
				mexErrMsgTxt("Error: wrong password for encrypted data, exiting...\n");
		}
	}
	
	// check if a session-struct should be returned
	if (nlhs > 0) {
		
		// map session object to matlab output struct
		// and assign to the first return argument
		plhs[0] = map_mef3_session(session, map_indices_flag);
		
	}		
	
	// free the session memory
	free_session(session, MEF_TRUE);
	
	// 
	return;
	
}

