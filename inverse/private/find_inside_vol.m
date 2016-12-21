function [inside, outside] = find_inside_vol(pos, vol)

% FIND_INSIDE_VOL locates dipole locations inside/outside the source
% compartment of a volume conductor model.
% 
% [inside, outside] = find_inside_vol(pos, vol)
%
% This function is obsolete and its use in other functions should be replaced 
% by inside_vol

% Copyright (C) 2003-2007, Robert Oostenveld
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

warning('find_inside_vol is obsolete and will be removed, please use ft_inside_vol');
inside  = ft_inside_vol(pos, vol);
% replace boolean vector with indexing vectors
outside = find(~inside);
inside  = find(inside);
