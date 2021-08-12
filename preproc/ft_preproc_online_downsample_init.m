function state = ft_preproc_online_downsample_init(factor)

% FT_PREPROC_ONLINE_DOWNSAMPLE_INIT initializes an downsampling model with the given factor.
% 
% Use as
%  state = ft_preproc_online_downsample_init(factor)
%
% See also PREPROC, FT_PREPROC_ONLINE_DOWNSAMPLE_APPLY

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

if factor < 1 || factor~=round(factor)
	error('Argument ''factor'' must be a positive integer number');
end

state.factor  = factor;
state.numSkip = 0;
