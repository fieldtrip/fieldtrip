
/*
 * Copyright (C) 2006, Robert Oostenveld, F.C. Donders Ccentre for Cognitive Neuroimaging
 * Copyright (C) 2007, Peter Desain, Nijmegen Institute for Cognition and Information
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *
 * $Log: read_24bit.c,v $
 * Revision 1.4  2009/06/08 16:33:08  roboos
 * include sys/types.h for all platforms, not only WIN32/64, since also required for OSX 10.4
 *
 * Revision 1.3  2009/03/30 07:23:03  roboos
 * changed ifdef syntax
 *
 * Revision 1.2  2009/03/28 05:32:46  roboos
 * changed ifdef handling to accomodate matlab71, which did not like the WIN32 || WIN64 combination
 *
 * Revision 1.1  2009/01/14 09:43:37  roboos
 * moved source code for mex file from fileio/mex to file/private
 * compiling the source code from within Matlab now ensures that the mex file will be located immediately at the right position
 *
 * Revision 1.7  2008/12/11 17:32:33  roboos
 * added conditional include and define for win32 & 64, thanks to Guillaume
 *
 * Revision 1.6  2008/10/07 07:53:05  roboos
 * replaced fseek with fseeko, which I had forgotten to fix yesterday
 *
 * Revision 1.5  2008/10/06 09:25:24  roboos
 * changed from long int to off_t and included stdint for uint32_t and uint64_t
 * this should allow the mex file also to work on large filed (>2GB) on supported platforms
 *
 * Revision 1.4  2008/04/09 10:10:20  roboos
 * corrected the text of an error message
 *
 * Revision 1.3  2007/10/02 08:49:37  roboos
 * added brackets around the bitshift operations
 *
 * Revision 1.2  2007/09/12 12:43:39  roboos
 * added peter to copyright, merged suggestions from nici (see email 20070821), changed to using bitshift, changed the handing of the sign bit
 *
 *
 */

/*
 * _LARGEFILE_SOURCE takes care that the offset is represented as 64 bit on supported platforms
 * it should be defined prior to including the system header files
 */

#define _LARGEFILE_SOURCE

#include <math.h>
#include <sys/types.h>
#include "mex.h"
#include "matrix.h"

#if defined(_WIN32) || defined(_WIN64)
#define int32_t INT32_T
#define int64_t INT64_T
#define fseeko fseek
#else
#include <stdint.h>
#endif

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  FILE *fp;
  mxArray *dat;
  double *dat_p;
  char *filename, *buf;
  int32_t fnlen, b1, b2, b3;
  
  /* these have to be able to handle large files (>2GB) on supported platforms */
  int64_t count, indx;
  size_t  num, numwords;
  off_t   offset;
  
  if (nrhs != 3)
    mexErrMsgTxt ("Invalid number of input arguments");
  
  if (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=1)
    mexErrMsgTxt ("Invalid dimension for input argument 2");
  if (mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
    mexErrMsgTxt ("Invalid dimension for input argument 3");
  
  /* get the filename */
  fnlen    = mxGetN(prhs[0]);
  filename = mxCalloc(fnlen+1, sizeof(char));
  mxGetString(prhs[0], filename, fnlen+1);
  
  /* get values and typecast to integer */
  offset   = (off_t)(mxGetScalar (prhs[1]));
  numwords = (off_t)(mxGetScalar (prhs[2]));
  
  /*
  printf("filename = %s\n", filename);
  printf("fnlen    = %d\n", fnlen);
  printf("offset   = %d\n", offset);
  printf("numwords = %d\n", numwords);
   */
  
  /* open the file and read the desired bytes */
  fp = fopen(filename, "rb");
  if (fp==NULL)
  {
    printf("error opening file: %s\n", filename);
    return;
  }
  fseeko(fp, offset, SEEK_SET);
  
  buf   = mxCalloc(3*numwords, sizeof(char));
  count = fread(buf, sizeof(char), 3*numwords, fp);
  if (count<3*numwords)
  {
    printf("error reading from %s\n", filename);
    fclose(fp);
    return;
  }
  else
  {
    fclose(fp);
  }
  
  /* convert every thee bytes into one word */
  dat   = mxCreateDoubleMatrix (1, numwords, mxREAL);
  dat_p = mxGetData (dat);
  
  for (count=0; count<numwords; count++)
  {
    indx = count*3;
    b1 = (0x000000FF & ((int32_t)buf[indx  ]));
    b2 = (0x000000FF & ((int32_t)buf[indx+1]));
    b3 = (0x000000FF & ((int32_t)buf[indx+2]));
    dat_p[count] = ((int32_t) ((b3 << 24) | (b2 << 16) | (b1 << 8)))/256;
  }
  
  /* assign the output parameters */
  plhs[0] = dat;
  return;
}

