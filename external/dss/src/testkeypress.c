#include <unistd.h>
#include <fcntl.h>
#include <stdio.h>
#include "mex.h"


/* Input Arguments */

#define C_IN    prhs[0]

/* Output Arguments */

#define C_OUT  plhs[0]

void mexFunction( int nlhs, mxArray *plhs[], 
                  int nrhs, const mxArray*prhs[] )
     
{
  char str[16];
  int flags;
  int fd;
  int mychar = 'q';
  int temp = 0;
    
  /* Check for proper number of arguments */
  if (nrhs != 1) { 
    mexErrMsgTxt("One input argument required.");
  } else if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments."); 
  } 

  if(! mxGetString(C_IN, str, 15)) {
    mychar = str[0];
  }

  fd = fileno(stdin);
  flags = fcntl(fd, F_GETFL);
  fcntl(fd, F_SETFL, flags|O_NONBLOCK);
  while (temp != EOF && temp != mychar) {
    temp = getchar();
  }
  fcntl(fd, F_SETFL, flags);

  C_OUT = mxCreateScalarDouble(temp == mychar);

  return;
}
