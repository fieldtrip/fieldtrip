function [DM, xd] = ft_preproc_online_downsample_apply(DM, x)

% function [DM, xd] = ft_preproc_online_downsample_apply(DM, x)
%
% Passes signal x (channels times samples) through the downsampler.
% Returns updated downsample model (numSkip!) and downsampled signal.

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

[dimX, numX] = size(x);

N  = size(x,2);
% to get number K of sample we can write out, subtract the skipped samples, 
% and then add maximum possible number of skip samples for next time (=DM.factor-1)
K  = floor((N - DM.numSkip + DM.factor-1)/DM.factor);
		
startIdx = 1+DM.numSkip;
endIdx   = 1+DM.numSkip + (K-1)*DM.factor;

xd = x(:,startIdx:DM.factor:endIdx);
DM.numSkip = DM.factor-1-(N-endIdx); % for next time
	
