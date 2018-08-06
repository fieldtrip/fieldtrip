function [H] = traditional(f)

% TRADITIONAL creates the homogenous spatial transformation matrix
% for a 9 parameter traditional "Talairach-model" transformation
%
% Use as
%   [H] = traditional(f)
%
% The transformation vector f should contain the 
%   x-shift
%   y-shift
%   z-shift
% followed by the
%   pitch (rotation around x-axis)
%   roll  (rotation around y-axis)
%   yaw   (rotation around z-axis)
% followed by the 
%   x-rescaling factor
%   y-rescaling factor
%   z-rescaling factor
%
% The order in which the transformations are done is exactly opposite as
% the list above, i.e. first z-rescale, ... and finally x-shift.

% Copyright (C) 2000-2005, Robert Oostenveld
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

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

% compute the homogenous transformation matrix for the translation
T = translate(f([1 2 3]));

% compute the homogenous transformation matrix for the rotation
R = rotate(f([4 5 6]));

% compute the homogenous transformation matrix for the scaling
S = scale(f([7 8 9]));

% compute the homogenous transformation matrix for the combination
H = T*R*S;
