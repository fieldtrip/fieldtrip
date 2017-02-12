#include <math.h>
#include "mex.h"
#include "geometry.h"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  char str[256];
  mxArray *sa;
  double *v1, *v2, *v3, p1[3], p2[3], p3[3], *pnt, *tri, *sa_p;
  int i1, i2, i3, npnt, ntri, i, on_triangle;

  if (nrhs==2)
  {
    /* the input consists of a vertex array and a triangle array */
    if (mxGetN(prhs[0])!=3)
      mexErrMsgTxt ("Invalid dimension for input argument 1");
    if (mxGetN(prhs[1])!=3)
      mexErrMsgTxt ("Invalid dimension for input argument 2");
    npnt = mxGetM(prhs[0]);
    ntri = mxGetM(prhs[1]);
    pnt = mxGetData (prhs[0]);
    tri = mxGetData (prhs[1]);
    sa = mxCreateDoubleMatrix (ntri, 1, mxREAL);
    sa_p = mxGetData (sa);
    /* compute the solid angle for all triangles */
    for (i=0; i<ntri; i++)
    {
      i1 = (int)(tri[i]) - 1;
      i2 = (int)(tri[i+ntri]) - 1;
      i3 = (int)(tri[i+ntri+ntri]) - 1;
      p1[0] = pnt[i1]; p1[1] = pnt[i1+npnt]; p1[2] = pnt[i1+npnt+npnt];
      p2[0] = pnt[i2]; p2[1] = pnt[i2+npnt]; p2[2] = pnt[i2+npnt+npnt];
      p3[0] = pnt[i3]; p3[1] = pnt[i3+npnt]; p3[2] = pnt[i3+npnt+npnt];
      *(sa_p+i) = solang(p1, p2, p3, &on_triangle);
      if (on_triangle)
        *(sa_p+i) = mxGetNaN();
    }
    /* assign the output parameter */
    plhs[0] = sa;
  }
  else if (nrhs==3)
  {
    /* the input consists of three vertices */
    if (mxGetM(prhs[0])!=1 || mxGetN(prhs[0])!=3)
      mexErrMsgTxt ("Invalid dimension for input argument 1");
    if (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=3)
      mexErrMsgTxt ("Invalid dimension for input argument 2");
    if (mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=3)
      mexErrMsgTxt ("Invalid dimension for input argument 3");
    v1 = mxGetData (prhs[0]);
    v2 = mxGetData (prhs[1]);
    v3 = mxGetData (prhs[2]);
    sa = mxCreateDoubleMatrix (1, 1, mxREAL);
    sa_p = mxGetData (sa);
    /* compute the solid angle for a single triangle */
    *(sa_p) = solang(v1, v2, v3, &on_triangle);
    if (on_triangle)
      *(sa_p) = mxGetNaN();
    /* assign the output parameter */
    plhs[0] = sa;
  }
  else
    mexErrMsgTxt ("Invalid number of input arguments");

  return;
}

