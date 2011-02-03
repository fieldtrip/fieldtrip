#include <math.h>
#include <string.h>
#include "mex.h"
#include "matrix.h"

/*
 * Copyright (C) 2008, Robert Oostenveld, F.C. Donders Ccentre for Cognitive Neuroimaging
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
 * $Log: read_16bit.c,v $
 * Revision 1.1  2009/01/14 09:43:37  roboos
 * moved source code for mex file from fileio/mex to file/private
 * compiling the source code from within Matlab now ensures that the mex file will be located immediately at the right position
 *
 * Revision 1.1  2008/04/09 10:04:37  roboos
 * new code based in 24bit, not tested
 *
 *
 *
 */

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  FILE *fp;
  mxArray *dat;
  double *dat_p;
  char *filename;
  short int *buf;
  long int offset, numwords, fnlen, i, j, num, b1, b2, b3;
  
  if (nrhs != 3)
    mexErrMsgTxt ("Invalid number of input arguments");
  
  if (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=1)
    mexErrMsgTxt ("Invalid dimension for input argument 2");
  if (mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=1)
    mexErrMsgTxt ("Invalid dimension for input argument 3");
  
  /* get the filename */
  fnlen    = mxGetN      (prhs[0]);
  filename = mxCalloc(fnlen+1, sizeof(char));
  mxGetString(prhs[0], filename, fnlen+1);
  
  /* get values and typecast to integer */
  offset   = (long int)(mxGetScalar (prhs[1]));
  numwords = (long int)(mxGetScalar (prhs[2]));
  
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
  fseek(fp, offset, SEEK_SET);
  
  buf = mxCalloc(numwords, sizeof(short));
  num = fread(buf, sizeof(short), numwords, fp);
  if (num<numwords)
  {
    printf("error reading from %s\n", filename);
    fclose(fp);
    return;
  }
  else
  {
    fclose(fp);
  }
  
  /* convert into matlab array */
  dat   = mxCreateDoubleMatrix (1, numwords, mxREAL);
  dat_p = mxGetData (dat);
  
  for (i=0; i<numwords; i++)
  {
    dat_p[i] = (double)(buf[i]);
  }
  
  /* assign the output parameters */
  plhs[0] = dat;
  return;
}

