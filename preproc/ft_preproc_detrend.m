function [dat,beta,x] = ft_preproc_detrend(dat, begsample, endsample)

% FT_PREPROC_DETREND removes mean and linear trend from the
% data using using a General Linear Modeling approach.
%
% Use as
%   [dat] = ft_preproc_detrend(dat, begin, end)
% where
%   dat        = data matrix (Nchans X Ntime)
%   begsample  = index of the begin sample for the trend estimate
%   endsample  = index of the end sample for the trend estimate
%
% If no begin and end sample are specified for the trend estimate, it
% will be estimated on the complete data.
%
% See also FT_PREPROC_BASELINECORRECT, FT_PREPROC_POLYREMOVAL

% Copyright (C) 2008-2014, Robert Oostenveld
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

% take the whole segment if begsample and endsample are not specified
if nargin<2 || isempty(begsample)
  begsample = 1;
end
if nargin<3|| isempty(endsample)
  endsample = size(dat,2);
end

[dat,beta,x] = ft_preproc_polyremoval(dat, 1, begsample, endsample);
