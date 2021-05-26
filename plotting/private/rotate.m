function [H] = rotate(f)

% ROTATE returns the homogenous coordinate transformation matrix
% corresponding to a rotation around the x, y and z-axis. The direction of
% the rotation is according to the right-hand rule.
%
% Use as
%   [H] = rotate(R)
% where
%   R       [rx, ry, rz] in degrees
%   H       corresponding homogenous transformation matrix
%
% Note that the order in which the rotations are performs matters. The
% rotation is first done around the z-axis, then the y-axis and finally the
% x-axis.
%
% See also TRANSLATE, SCALE, QUATERNION

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

if numel(f)~=3
  ft_error('incorrect input vector');
end

% convert degrees to radians
f = f*pi/180;

% get the individual angles (in radians)
rx = f(1);
ry = f(2);
rz = f(3);

% precompute the sin/cos values of the angles
cX = cos(rx);
cY = cos(ry);
cZ = cos(rz);
sX = sin(rx);
sY = sin(ry);
sZ = sin(rz);

% according to Roger Woods' http://bishopw.loni.ucla.edu/AIR5/homogenous.html
% it should be this, but I cannot reproduce his rotation matrix
% H = eye(4,4);
% H(1,1) = cZ*cY + sZ*sX*sY;
% H(1,2) = sZ*cY - cZ*sX*sY;
% H(1,3) =            cX*sY;
% H(2,1) = -sZ*cX;
% H(2,2) =  cZ*cX;
% H(2,3) =     sX;
% H(3,1) =  sZ*sX*cY - cZ*sY;
% H(3,2) = -cZ*sX*cY - sZ*sY;
% H(3,3) =             cX*cY;

% instead, the following rotation matrix does work according my
% expectations. It rotates according to the right hand rule and first
% rotates around z, then y and then x axis
H = [
           cZ*cY,          -sZ*cY,              sY,               0
  cZ*sY*sX+sZ*cX, -sZ*sY*sX+cZ*cX,          -cY*sX,               0
 -cZ*sY*cX+sZ*sX,  sZ*sY*cX+cZ*sX,           cY*cX,               0
               0,               0,               0,               1
];

if 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code can be used to construct the combined rotation matrix
% for either xyz or zyx ordering (using the MATLAB symbolic math toolbox)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  syms sX sY sZ cX cY cZ
  % this is for only rotating around x
  Rx = [
    1   0    0   0
    0   cX  -sX  0
    0   sX   cX  0
    0   0    0   1
    ];
  % this is for only rotating around y
  Ry = [
    cY  0   sY  0
    0   1   0   0
    -sY  0   cY  0
    0   0   0   1
    ];
  % this is for only rotating around z
  Rz = [
    cZ -sZ  0   0
    sZ  cZ  0   0
    0   0   1   0
    0   0   0   1
    ];
  % combine them
  Rzyx = Rz * Ry * Rx  % rotate around x, y, then z
  Rxyz = Rx * Ry * Rz  % rotate around z, y, then x
end
