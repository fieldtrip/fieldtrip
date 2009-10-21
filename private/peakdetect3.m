function [pindx, pval] = peakdetect3(dat, threshold, mindist)

% PEAKDETECT3 detects peaks above a certain threshold in single-channel data
%
% Use as
%   [pindx, pval] = peakdetect3(dat, threshold, mindist)
%
% See also PEAKDETECT, PEAKDETECT2

% Copyright (C) 2000-2005, Robert Oostenveld
%
% $Log: peakdetect3.m,v $
% Revision 1.2  2005/11/22 16:12:28  roboos
% fixed bug in mindist, did not work if no peaks were found
%
% Revision 1.1  2005/11/22 16:01:10  roboos
% new version for cvs, based on old code
%

% threshold the data
tr = dat>threshold;

% the derivative of the data changes its sign at the peak
td = diff(dat);
td = td(1:(end-1))>0 & td(2:end)<0;
td = [0 td 0];

pindx = find(td & tr);

if nargin>2 && length(pindx)>0
  % find the peaks that are too close to each other
  pd = [inf diff(pindx)];
  pindx = pindx(pd>mindist);
end

if nargout>1
  pval = dat(pindx);
end

