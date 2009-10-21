#include <string.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  mxArray *lf;
  double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, eps, alpha, beta, A, B, C;
  double r[3], u[3], d[3];
  double *lf_p, *rm_p, *um_p, *R;
  int nchan, i;
  char str[256];

  if (nrhs<3)
    mexErrMsgTxt("Not enough input arguments");
  if (nrhs>3)
    mexErrMsgTxt("Too many input arguments");

  if (mxGetM(prhs[0])!=1 || mxGetN(prhs[0])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 1");
  if (mxGetN(prhs[1])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 2");
  if (mxGetN(prhs[2])!=3)
    mexErrMsgTxt ("Invalid dimension for input argument 3");
  nchan = mxGetM(prhs[1]);
  if (mxGetM(prhs[2])!=nchan)
    mexErrMsgTxt ("Number of channels does not match between argument 2 and 3");

  lf = mxCreateDoubleMatrix(nchan, 3, mxREAL);
  lf_p = mxGetData(lf);
  R    = mxGetData(prhs[0]);
  rm_p = mxGetData(prhs[1]);
  um_p = mxGetData(prhs[2]);
  eps  = mxGetEps();

  tmp2 = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);	/* norm(R) */

  for (i=0; i<nchan; i++)
  {
    /* get the position of this single channel */
    r[0] = rm_p[i];
    r[1] = rm_p[i+nchan];
    r[2] = rm_p[i+nchan+nchan];
    /* get the orientation of this single channel */
    u[0] = um_p[i];
    u[1] = um_p[i+nchan];
    u[2] = um_p[i+nchan+nchan];

    /* compute the difference between this channel and the dipole position */
    d[0] = r[0] - R[0];
    d[1] = r[1] - R[1];
    d[2] = r[2] - R[2];

    tmp1 = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);	/* norm(r) */
    /* tmp2 = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]); */
    tmp3 = sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);	/* norm(d) */
    tmp4 = r[0]*R[0] + r[1]*R[1] + r[2]*R[2];		/* dot(r, R) */
    tmp5 = r[0]*d[0] + r[1]*d[1] + r[2]*d[2];		/* dot(r, d) */
    tmp6 = R[0]*d[0] + R[1]*d[1] + R[2]*d[2];		/* dot(R, d) */
    tmp7 = tmp1*tmp2*tmp1*tmp2 - tmp4*tmp4;		/* norm(cross(r,R))^2 */

    alpha = 1 / (-tmp3 * (tmp1*tmp3+tmp5));
    A = 1/tmp3 - 2*alpha*tmp2*tmp2 - 1/tmp1;
    B = 2*alpha*tmp4;
    C = -tmp6/(tmp3*tmp3*tmp3);
    
    if (tmp7<eps)
      beta = 0;
    else
      beta = ((A*r[0] + B*R[0] + C*d[0])*u[0] +
              (A*r[1] + B*R[1] + C*d[1])*u[1] +
              (A*r[2] + B*R[2] + C*d[2])*u[2]) / tmp7;

    /* re-use the temporary array d for something else */
    d[0] = alpha*u[0] + beta*r[0];
    d[1] = alpha*u[1] + beta*r[1];
    d[2] = alpha*u[2] + beta*r[2];

    /* lf(i,:) = 1e-7 * cross(alpha*u  + beta*r, R) */
    lf_p[i            ] = (d[1]*R[2] - d[2]*R[1])*1e-7; 
    lf_p[i+nchan      ] = (d[2]*R[0] - d[0]*R[2])*1e-7; 
    lf_p[i+nchan+nchan] = (d[0]*R[1] - d[1]*R[0])*1e-7; 
  }

  /* assign the output parameter */
  plhs[0] = lf;

  return;
}
