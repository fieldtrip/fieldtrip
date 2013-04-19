/* DIST_EUCLIDEAN     Evaluate the distance matrix.
 *
 *         Description
 *         C = DIST_EUCLIDEAN(METRIC, TX) takes in metric structure METRIC and 
 *         matrix TX that contains input vectors to GP. Returns 
 *         covariance matrix C. Every element ij of C contains covariance 
 *         between inputs i and j in TX.
 *
 *
 * Last modified: 2011-09-07 10:21:14 EEST
 *
 */

/* Copyright (C) 1998-2001 Aki Vehtari
 * Copyright (C) 2008-2011 Jarno Vanhatalo
 * 
 *This software is distributed under the GNU General Public 
 *License (version 3 or later); please refer to the file 
 *License.txt, included with the software, for details.
 *
 */

#include<string.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"
#define max(a,b) (((a) > (b)) ? (a) : (b))

void mexFunction(const int nlhs, mxArray *plhs[],
		 const int nrhs, const mxArray *prhs[])
{

  if (nlhs!=1)
    mexErrMsgTxt( "Wrong number of output arguments.");
  
  if (nrhs!=2)
    mexErrMsgTxt( "Wrong number of input arguments." );
  
  {
    double *x, *l, *C, rr, d, *cc;
    const mwSize *dims;
    mxLogical *deltad;
    char *type;
    mwIndex i, j, k, i2;
    mwSize m, n, ncomp, lr;
    mxArray *field;
    const mxArray *components, *components_element;
    
    dims = mxGetDimensions(prhs[1]);
    x = mxGetPr(prhs[1]);
    m = dims[0];
    n = dims[1];
   
    if((field=mxGetField(*prhs, 0, "type"))==NULL)
      mexErrMsgTxt("Could not get metric.type");
    if (mxIsChar(field) !=1)
      mexErrMsgTxt( "metric.type must be a string." );
    type = mxArrayToString(field);

    
    /*
     * squared exponential covariance
     */
    if( strcmp( type, "metric_euclidean" ) == 0 ) {
      if((field=mxGetField(*prhs, 0, "lengthScale"))==NULL)
	mexErrMsgTxt("Could not get metric.lengthScale");
      dims = mxGetDimensions(field);
      if (dims[0]!=1 && dims[1]!=1)
	mexErrMsgTxt( "metric.lengthScale must be a scalar or a vector." );
      lr=max(dims[0],dims[1]);
      l = mxGetPr(field);

      if((components=mxGetField(*prhs, 0, "components"))==NULL)
	mexErrMsgTxt("Could not get metric.components");
      dims = mxGetDimensions(components);
      if (dims[0]!=1 && dims[1]!=1)
	mexErrMsgTxt( "metric.components must be a one dimensional cell array." );
      ncomp = max(dims[0],dims[1]);
      if (lr!=1 && lr!=ncomp)
	mexErrMsgTxt( "metric.lengthScale must be scalar or its lenght must be the same as the number of components." );

      if((field=mxGetField(*prhs, 0, "deltadist"))==NULL)
	mexErrMsgTxt("Could not get metric.deltadist");
      dims = mxGetDimensions(field);
      if (dims[0]!=1 && dims[1]!=1)
	mexErrMsgTxt( "metric.deltadist must be a logical vector." );
      if (max(dims[0],dims[1])!=ncomp)
	mexErrMsgTxt( "metric.logical must be of the same length as the number of components." );
      deltad = (mxLogical *)mxGetData(field);      

      /* evaluate the distance */
      plhs[0]=mxCreateDoubleMatrix(m, m, mxREAL);
      C = mxGetPr(plhs[0]);
      for (i=0;i<ncomp;i++) {
	rr=(lr>1)?(l[i]*l[i]):(l[0]*l[0]);
	components_element = mxGetCell(components, i);
	dims = mxGetDimensions(components_element);
	cc = mxGetPr(components_element);
	/*printf("hep l[%d]= %.2f \n", i, rr);*/
	for (i2=0;i2<max(dims[0],dims[1]);i2++) {
	  /*printf("cc[%d] = %d \n", i2, (mwIndex)cc[i2]);*/	  
	  for (j=0;j<m;j++) {
	    for (k=0;k<j;k++) {
	      d=x[j+m*((mwIndex)cc[i2]-1)]-x[k+m*((mwIndex)cc[i2]-1)];
	      if (deltad[i]==true){
		C[j*m+k]+=(d!=0)?(1/rr):(0); 
	      }else{ 
		C[j*m+k]+=d*d/rr; 
	      }
	    }
	  }
	}
      }
      for (j=0;j<m;j++) {
	for (k=0;k<j;k++) {
	  d=sqrt(C[j*m+k]);
	  C[j*m+k]=d;
	  C[j+k*m]=d;
	}
	C[j*(m+1)]=0;
      }
    }
    else{
      mexErrMsgTxt( "Wrong type of metric." );
    }
  }
  return;
}
