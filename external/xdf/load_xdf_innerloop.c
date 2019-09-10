#include "mex.h"
/*
 * [Values,Timestamps] = load_xdf_innerloop(Data, NumChannels, FormatString, SamplingInterval, LastTimestamp);
 * MEX kernel that implements the inner loop of the load_xdf function.
 * 
 * In:
 *   Data : byte array that contains the data of the chunk
 *
 *   NumChannels : number of channels per sample (double)
 *
 *   FormatString : format string of the data (determines the format of the data values); can be one of the following:
 *                  - '*int8'
 *                  - '*int16'
 *                  - '*int32'
 *                  - '*int64'
 *                  - '*float32'
 *                  - '*double64'
 *                  - '*string'
 *
 *   SamplingInterval : sampling interval in seconds (can be 0, if irregular)
 *
 *   LastTimestamp : time stamp of the last sample preceding the data, in seconds
 *
 * Out:
 *   Values : [#Channels x #Samples] array containing the chunk data; either numeric or cell array
 *
 *   Timestamps: #Samples array of the time stamps for every sample, in seconds
 *
 *                               Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
 *                               2012-08-06
 */

/* Value format of a channel. */
typedef enum {
	cft_undefined = 0,
	cft_float32 = 1,
	cft_double64 = 2,
	cft_string = 3,
	cft_int32 = 4,
	cft_int16 = 5,
	cft_int8 = 6,
	cft_int64 = 7,
} value_format_t;

/* Number of bytes occupied by a single value (excluding strings, which are variable-length). */
unsigned value_bytes[] = {0,4,8,0,4,2,1,8};

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] ) 
{
    /* input arguments */
    unsigned char *data;
    unsigned num_channels;
    char format_string[1024];
    double sampling_interval;
    double last_timestamp;

    /* output arguments */
    unsigned char *values;
    double *timestamps;
    
    /* variables */
    unsigned num_samples, s, c, sample_bytes;
    value_format_t format_code;
    mwSize string_dims[2];
    
    /* temporaries */
    unsigned bytes, num_chars;
    char *tmp;
    
    /* check input/output argument formats */   
    if (nrhs != 5)
        mexErrMsgTxt("5 input arguments required."); 
    if (nlhs != 2)
        mexErrMsgTxt("2 output arguments required."); 
   
    /* get inputs */
    data = (unsigned char*)mxGetData(prhs[0]);
    num_channels = (unsigned)(*((double*)mxGetData(prhs[1])));
    mxGetNChars_700(prhs[2], format_string, mxGetNumberOfElements(prhs[2])+1);
    sampling_interval = *((double*)mxGetData(prhs[3]));
    last_timestamp = *((double*)mxGetData(prhs[4]));

    /* set the format code */
    format_code = cft_undefined;
    if (strcmp(format_string,"*float32") == 0)
        format_code = cft_float32;
    if (strcmp(format_string,"*double64") == 0)
        format_code = cft_double64;
    if (strcmp(format_string,"*string") == 0)
        format_code = cft_string;
    if (strcmp(format_string,"*int32") == 0)
        format_code = cft_int32;
    if (strcmp(format_string,"*int16") == 0)
        format_code = cft_int16;
    if (strcmp(format_string,"*int8") == 0)
        format_code = cft_int8;
    if (strcmp(format_string,"*int64") == 0)
        format_code = cft_int64;
    if (format_code == cft_undefined)
        mexErrMsgTxt("The FormatString does not contain a recognized format.");     
    
    /* read number of samples (as a varlen int) */
    bytes = *data++;
    if (bytes == 1)
        num_samples = *data++;
    if (bytes == 4) {
        num_samples = *(unsigned*)data;
        data += 4;
    }
    if (bytes == 8)
        mexErrMsgTxt("This importer cannot yet handle chunks with more than 4 billion samples.");
    
    /* allocate memory for the time stamps*/
    plhs[1] = mxCreateNumericMatrix(1,num_samples,mxDOUBLE_CLASS,mxREAL);    
    timestamps = (double*)mxGetData(plhs[1]);    
    if (format_code != cft_string) {
        /* numeric data case: allocate output array */
        switch(format_code) {
            case cft_float32:
                plhs[0] = mxCreateNumericMatrix(num_channels,num_samples,mxSINGLE_CLASS,mxREAL);
                break;
            case cft_double64:
                plhs[0] = mxCreateNumericMatrix(num_channels,num_samples,mxDOUBLE_CLASS,mxREAL);
                break;
            case cft_int32:
                plhs[0] = mxCreateNumericMatrix(num_channels,num_samples,mxINT32_CLASS,mxREAL);
                break;
            case cft_int16:
                plhs[0] = mxCreateNumericMatrix(num_channels,num_samples,mxINT16_CLASS,mxREAL);
                break;
            case cft_int8:
                plhs[0] = mxCreateNumericMatrix(num_channels,num_samples,mxINT8_CLASS,mxREAL);
                break;
            case cft_int64:
                plhs[0] = mxCreateNumericMatrix(num_channels,num_samples,mxINT64_CLASS,mxREAL);
                break;
            default:
                mexErrMsgTxt("Unrecognized data type.");
        }
        sample_bytes = value_bytes[format_code]*num_channels;
        values = (unsigned char*)mxGetData(plhs[0]);
        /* for each sample... */
        for (s=0;s<num_samples;s++) {
            /* deduce timestamp */
            if (bytes = *data++) {
                if (bytes != 8)
                    mexErrMsgTxt("Parsing error in the file format.");
                last_timestamp = *(double*)data;
                data += bytes;
            } else
                last_timestamp += sampling_interval;
            *timestamps++ = last_timestamp;
            /* copy data (assumes little endian...) */
            memcpy(values,data,sample_bytes);
            values += sample_bytes;
            data += sample_bytes;
        }
    } else {
        /* string-formatted data case: allocate output array */
        plhs[0] = mxCreateCellMatrix(num_channels,num_samples);
        string_dims[0] = 1;
        /* for each sample... */
        for (s=0;s<num_samples;s++) {
            /* deduce timestamp */
            if (bytes = *data++) {
                last_timestamp = *(double*)data;
                data += 8;
            } else
                last_timestamp += sampling_interval;
            *timestamps++ = last_timestamp;
            
            /* for each channel... */
            for (c=0;c<num_channels;c++) {
                /* read string length */
                bytes = *data++;
                if (bytes == 1)
                    num_chars = *data++;
                if (bytes == 4) {
                    num_chars = *(unsigned*)data;
                    data += bytes;
                }
                if (bytes == 8)
                    mexErrMsgTxt("This importer cannot yet handle strings with more than 4 billion characters.");
                string_dims[1] = num_chars;
                /* copy data */
                tmp = malloc(num_chars+1); 
                memcpy(tmp, data, num_chars);
                tmp[num_chars] = 0;
                mxSetCell(plhs[0],c + s*num_channels,mxCreateString(tmp));
                free(tmp);
                data += num_chars;
            }
        }
    }
}
