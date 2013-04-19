function [xn, xmean, xstd] = normdata(x,xmean,xstd)
%NORMDATA  Normalize input to zero mean and unit variance
%
%  Description
%    [XN, XMEAN, XSTD] = NORMDATA(X) normalizes X to zero mean and
%    unit variance. Returns normalized XN, mean of X XMEAN and
%    standard deviance of X XSTD.
%  
%    XN = NORMDATA(X,XMEAN,XSTD) normalizes X using precomputed
%    XMEAN and XSTD.
%
%  See also DENORMDATA
%
  
% Copyright (c) 2010 Aki Vehtari

% This software is distributed under the GNU General Public 
% License (version 3 or later); please refer to the file 
% License.txt, included with the software, for details.

  if nargin<=1
    xmean=nanmean(x);
    xstd=nanstd(x);
  end
  xn=bsxfun(@rdivide,bsxfun(@minus,x,xmean),xstd);
