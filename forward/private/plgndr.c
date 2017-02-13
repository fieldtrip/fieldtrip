/*
 *  The original plgndr function was implemented using some code from
 *  Numerical Recipes in C. However, it was not allowed to release that code
 *  as open source. The new implementation below is using some code from
 *  the GNU Scientific Library (http://www.gnu.org/software/gsl).
 *
 *  Copyright (C) 2002-2006 Robert Oostenveld
 *  Copyright (C) 2006, Thomas Hartmann
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include <math.h>
#include "mex.h"

double legendre_Pmm(double m, double x)
{
  if(m == 0)
  {
    return 1.0;
  }
  else
  {
    double p_mm = 1.0;
    double root_factor = sqrt(1.0-x)*sqrt(1.0+x);
    double fact_coeff = 1.0;
    int i;
    for(i=1; i<=m; i++)
    {
      p_mm *= -fact_coeff * root_factor;
      fact_coeff += 2.0;
    }
    return p_mm;
  }
}

double plgndr(int l, int m, double x)
{
  /* these are constant, but can only be assigned after checking the input arguments */
  double dif, sum, t_d, t_s, exp_check, err_amp;  

  double p_mm, p_mmp1;
  double result;
  
  /* determine whether we have correct input arguments */
  if (m < 0 || m > l || fabs(x) > 1.0)
    mexErrMsgTxt ("Bad arguments in routine plgndr");
  
  dif = l-m;
  sum = l+m;
  t_d = ( dif == 0.0 ? 0.0 : 0.5 * dif * (log(dif)-1.0) );
  t_s = ( dif == 0.0 ? 0.0 : 0.5 * sum * (log(sum)-1.0) );
  exp_check = 0.5 * log(2.0*l+1.0) + t_d - t_s;
  err_amp = 1.0 / (0.00000000001 + fabs(1.0-fabs(x)));
  p_mm   = legendre_Pmm(m, x);
  p_mmp1 = x * (2*m + 1) * p_mm;
  
  /* P_m^m(x) and P_{m+1}^m(x) */
  
  if(l == m)
  {
    result = p_mm;
  }
  else if(l == m + 1)
  {
    result = p_mmp1;
  }
  else
  {
    /*
     * upward recurrence: (l-m) P(l,m) = (2l-1) z P(l-1,m) - (l+m-1) P(l-2,m)
     * start at P(m,m), P(m+1,m)
     */
    
    double p_ellm2 = p_mm;
    double p_ellm1 = p_mmp1;
    double p_ell = 0.0;
    int ell;
    
    for(ell=(int)m+2; ell <= l; ell++)
    {
      p_ell = (x*(2*ell-1)*p_ellm1 - (ell+m-1)*p_ellm2) / (ell-m);
      p_ellm2 = p_ellm1;
      p_ellm1 = p_ell;
    }
    result = p_ell;
  }
  return result;
}

void
mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[])
{
  int l, m;
  double x;
  double *pd;
  
  if (nrhs != 3)
    mexErrMsgTxt ("invalid number of arguments for PLGNDR");
  
  l = mxGetScalar (prhs[0]);
  m = mxGetScalar (prhs[1]);
  x = mxGetScalar (prhs[2]);
  
  plhs[0] = mxCreateDoubleMatrix (1, 1, mxREAL);
  pd = mxGetData (plhs[0]);
  pd[0] = plgndr(l, m, x);
  
  return;
}
