function M = ft_omri_volume_to_mosaic(V)

% function M = ft_omri_volume_to_mosaic(V)
% 
% Reshuffle [MxNxS] volume to [XxY] mosaic

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

[h,w,ns] = size(V);
nn = ceil(sqrt(ns));
mm = ceil(ns/nn);

M = zeros(w*mm, h*nn, class(V));
row = 1;
col = 1;
for n = 1:ns
	is = 1 + (row-1)*h;
	ie = row*h;
	js = 1 + (col-1)*w;
	je = col*w;
		
	M(is:ie, js:je) = rot90(V(:,:,n));
	col = col + 1;
	if col>nn
	   row = row + 1;
	   col = 1;
	end
end
