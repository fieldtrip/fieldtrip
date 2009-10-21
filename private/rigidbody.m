function [H] = rigidbody(f);

% RIGIDBODY creates the homogenous spatial transformation matrix
% for a 6 parameter rigid-body transformation 
%
% Use as
%   [H] = rigidbody(f)
%
% The transformation vector f should contain the 
%   x-shift
%   y-shift
%   z-shift
% followed by the
%   pitch (rotation around x-axis)
%   roll  (rotation around y-axis)
%   yaw   (rotation around z-axis)

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

% $Log: rigidbody.m,v $
% Revision 1.4  2006/04/13 10:37:38  roboos
% added a ; to the end of a line
%
% Revision 1.3  2005/08/15 08:15:32  roboos
% reimplemented the rotate function, which contained an error (the error is in the AIR technical reference)
% changed all functions to be dependent on the rotate, translate and scale function
% all functions now behave consistenly, which also means that they are not compleetly backward compatible w.r.t. the order of the rotations
%
% Revision 1.2  2004/05/19 09:57:07  roberto
% added GPL copyright statement, added CVS log item
%

% compute the homogenous transformation matrix for the translation
T = translate(f([1 2 3]));

% compute the homogenous transformation matrix for the rotation
R = rotate(f([4 5 6]));

% compute the homogenous transformation matrix for the combination
H = T*R;

