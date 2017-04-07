function [pindx, pval] = peakdetect3(dat, threshold, mindist)

% PEAKDETECT3 detects peaks above a certain threshold in single-channel data
%
% Use as
%   [pindx, pval] = peakdetect3(dat, threshold, mindist)
%
% See also PEAKDETECT, PEAKDETECT2

% Copyright (C) 2000-2005, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

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

