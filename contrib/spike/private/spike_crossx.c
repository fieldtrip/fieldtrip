 /* CrossCorr
 * Cross-correlation function
 * MEX file
 * 
 * martin vinck 2010, adapted from Francesco Battaglia's crosscorr script
 *
 * input: t1, t2: two time series to cross correlate in seconds
 * Important: the timeseries should be sorted! 
 *        binsize: the binsize for the cross corr histogram in seconds
 *        nbins: the number of bins
 * output: C the cross correlation histogram
 *
 */

#include "mex.h"
#include <math.h>
#include <matrix.h>

void mexFunction(
  int nOutputs, mxArray *pointerOutputs[],
  int nInputs, const mxArray *pointerInputs[])
{
  
  double *tX;
  double *tY;
  double binSize;
  double *C;
  double lBound,uBound, minLag;
  int nX,nY;
  int nBins;
  int iBin,iX,jY,yStartIndx;
  double d,t1;
  int binNum;
  int c;
          
  nX = mxGetM(pointerInputs[0]) * mxGetN(pointerInputs[0]);
  nY = mxGetM(pointerInputs[1]) * mxGetN(pointerInputs[1]);

  tX = mxGetPr(pointerInputs[0]);
  tY = mxGetPr(pointerInputs[1]);
  
  binSize = mxGetScalar(pointerInputs[2]);
  nBins   = (int)mxGetScalar(pointerInputs[3]);
  double bins[nBins];
  
  /* create the pointer to the first output matrix */
  pointerOutputs[0] = mxCreateDoubleMatrix(nBins, 1, mxREAL);
  C                 = mxGetPr(pointerOutputs[0]);  
  minLag = -binSize*(nBins-1.0)/2.0;
              
  for(iBin = 0; iBin < nBins; iBin++)
  {  
    bins[iBin] = minLag + ((double)(iBin))*binSize;
  }
   
  /* compute the actual cross-correlations */
  minLag = -(nBins-1)*binSize/2.0;
  yStartIndx = 0;
  for(iX=0; iX<nX; iX++)
  {
      if (tX[iX]<(tY[yStartIndx]+minLag))
        continue;
      
      /* first determine where to start */
      t1 = tX[iX];
      lBound = t1 + minLag; 
      uBound = t1 - minLag;
      while(tY[yStartIndx] < lBound && yStartIndx < (nY-1))
        yStartIndx++;
      
      if (tY[yStartIndx]>uBound || tY[yStartIndx]<lBound )
      {
        continue;
      }
      for(jY = yStartIndx; jY < nY; jY++)
      {
        /* find the binnumber associated with the distance */
        d      = tY[jY] - t1;
        binNum =  (int)(floor((d - minLag)/binSize));              
        if (binNum<0)
        {
          binNum = 0;
        }                
        if (binNum>(nBins-2))
        {
          break;
        }
        C[binNum]++;
      }
  }     
}
  
  
		 
