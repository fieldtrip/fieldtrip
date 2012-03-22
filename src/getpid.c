#include <sys/types.h> 

/*
 * This mex file fails to compile using the LCC compiler that is included 
 * in 32-bit version of MATLAB. With MSVC 2008 it does compile.
 *
 * See also http://bugzilla.fcdonders.nl/show_bug.cgi?id=1384
 */

#if defined(_WIN32) || defined(_WIN64)
    #include <process.h>
#else
    #include <unistd.h> 
#endif

#include "mex.h" 

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) 
{ 
  unsigned int *pr; 
  plhs[0] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL); 
  pr = mxGetData(plhs[0]); 
  *pr = (unsigned int)getpid(); 
} 
