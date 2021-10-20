function [state, y] = ft_preproc_online_downsample_apply(state, x)

% FT_PREPROC_ONLINE_DOWNSAMPLE_APPLY passes a signal through the online downsampler
% and returns the downsampler state and the downsampled signal. The state keeps track
% of the number of samples to be skipped in the next call.
%
% Use as
%    [state, dat] = ft_preproc_online_downsample_apply(state, x)
% where
%   dat   = Nchan x Ntime
%   state = downsampler state, see FT_PREPROC_ONLINE_DOWNSAMPLE_INIT
%
% See also PREPROC

% Copyright (C) 2010, Stefan Klanke
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

N = size(x,2);

% to get number K of sample we can write out, subtract the skipped samples,
% and then add maximum possible number of skip samples for next time (=state.factor-1)
K  = floor((N - state.numSkip + state.factor-1)/state.factor);

startIdx = 1+state.numSkip;
endIdx   = 1+state.numSkip + (K-1)*state.factor;

y = x(:,startIdx:state.factor:endIdx);
state.numSkip = state.factor-1-(N-endIdx); % for next time

