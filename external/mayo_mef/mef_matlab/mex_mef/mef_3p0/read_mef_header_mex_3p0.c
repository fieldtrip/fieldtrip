// READ_MEF_HEADER_MEX_3P0.C mex gatefun to read universal header of data series segment
//
// See read_mef_header_mex_3p0.m for details of usage.

// Copyright 2020 Richard J. Cui. Created: Sun 02/02/2020  5:18:29.851 PM
// $Revision: 0.5 $  $Date: Wed 04/01/2020 10:09:33.460 PM$
//
// Multimodel Neuroimaging Lab (Dr. Dora Hermes)
// Mayo Clinic St. Mary Campus
// Rochester, MN 55905, USA
//
// Email: richard.cui@utoronto.ca

#include <ctype.h>
#include "mex.h"
#include "matmef/meflib.h"
#include "matmef/matmef_mapping.h"
#include "matmef/mex_datahelper.h"

/**************************************************************************
 * subroutines
 *************************************************************************/

/**************************************************************************
 * main entrence
 *************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // vars declare
    si1 channel_path[MEF_FULL_FILE_NAME_BYTES];
    char *mat_channel_path, *mat_password;
    si1 *password = NULL;
    si1 password_arr[PASSWORD_BYTES] = {0};
    si1 map_indices_flag = 1;
    si1 mat_map_indices_flag;
    CHANNEL *channel;
    
    // ***** channel_path *****
    // check input prhs[0]: channel_path (required)
    if (nrhs < 1) {
        mexErrMsgIdAndTxt("read_mef_header_mex_3p0:noChannelPathArg",
                "channel_path input argument is required");
    }
    else {
        if (mxIsEmpty(prhs[0])) {
            mexErrMsgIdAndTxt("read_mef_header_mex_3p0:emptyChannelPathArg",
                    "channel_path is empty");
        }
        if (!mxIsChar(prhs[0])) { // if channel_path not char string
            mexErrMsgIdAndTxt("read_mef_header_mex_3p0:invalidChannelPathArg",
                    "channel_path should be character string");
        }
    }
    
    // set the channel path
    mat_channel_path = mxArrayToString(prhs[0]);
    MEF_strncpy(channel_path, mat_channel_path, MEF_FULL_FILE_NAME_BYTES);
    
    // ***** password *****
    // check input prhs[1]: password (optional)
    if (nrhs > 1) {
		// note: if the password passed to any of the meflib read function is an empty string, than 
		//		 the 'process_password_data' function in 'meflib.c' will crash everything, so make
		// 		 sure it is either NULL or a string with at least one character
		
		// check if the password input argument is not empty
        if (!mxIsEmpty(prhs[1])) { // if password is not empty
            if (!mxIsChar(prhs[1])) { // if password not char string
                mexErrMsgIdAndTxt("read_mef_header_mex_3p0:invalidPasswordArg",
                        "password should be character string");
            }
            
            // set the password
            mat_password = mxArrayToString(prhs[1]);
            password = strcpy(password_arr, mat_password);
        }
    }
    
    // ***** map indices *****
    // check input prhs[2]: map_indices_flag (optional)
    if (nrhs > 2) {
        if (!mxIsEmpty(prhs[2])) { // if not empty, process
            // check if single numeric or logical
            if ((!mxIsNumeric(prhs[2])
            && !mxIsLogical(prhs[2]))
            || mxGetNumberOfElements(prhs[2]) > 1) {
                mexErrMsgIdAndTxt("read_mef_header_mex_3p0:invalidMapIndicesArg",
                        "map_indices_flag should be a single value, logical or numeric");
            }
            
            // retrieve the map indices flag value
            mat_map_indices_flag = (si1) mxGetScalar(prhs[2]);
            
            // check the value
            if (mat_map_indices_flag != 0 && mat_map_indices_flag != 1) {
                mexErrMsgIdAndTxt("read_mef_header_mex_3p0:invalidMapIndicesArg",
                        "map_indices_flag should be 0 (false) or 1 (true)");
            }
            
            // set the flag
            map_indices_flag = mat_map_indices_flag;
        }
    }
    
    // ***** get channel metadata *****
    initialize_meflib(); // initialize MEF library
    
    // read the channel object
    channel = read_MEF_channel( NULL, // allocate new channel object
            channel_path, // channel file path
            TIME_SERIES_CHANNEL_TYPE, // channel type
            password, // password
            NULL, // empty password data
            MEF_FALSE, // do not read time series data
            MEF_FALSE // do not read record data
            );
    
    // error check
    if (channel == NULL) { // nothing read
        mexErrMsgIdAndTxt("read_mef_header_mex_3p0:invalidChannel",
                "error while reading channel information");
    }
    
    if (channel->channel_type != TIME_SERIES_CHANNEL_TYPE) {
        mexErrMsgIdAndTxt("read_mef_header_mex_3p0:invalidChannelType",
                "not a time series channel");
    }
    
    // check encryption requirement
    if (channel->metadata.section_1 != NULL) {
        if (channel->metadata.section_1->section_2_encryption > 0
                && channel->metadata.section_1->section_3_encryption > 0) {
            // if the data is still encrypted
            if (password == NULL)
                mexErrMsgIdAndTxt("read_mef_header_mex_3p0:noPassword",
                        "data are encrypted, but no password is provided");
            else
                mexErrMsgIdAndTxt("read_mef_header_mex_3p0:invalidPassword",
                        "data are encrypted, but the password is invalid");
        }
    }
    
    // check number of segments
    if (channel->number_of_segments < 1)
        mexErrMsgIdAndTxt("read_mef_header_mex_3p0:noSegment",
                "no data segment in channel");
    
    // ***** output results *****
    // output plhs[0]: segment universal header
    if (nlhs > 0) plhs[0] = map_mef3_segment_universal_header(
            channel->segments[0].time_series_data_fps->universal_header);
    // output plhs[1] channel metadata
    // TODO: if map_indices_flag = 0, map_mef3_channel will crash
    if (nlhs > 1) plhs[1] = map_mef3_channel(channel, map_indices_flag);
    
    return;
}

// [EOF]