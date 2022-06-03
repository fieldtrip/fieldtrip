function f = homogenous2traditional(H)

% HOMOGENOUS2TRADITIONAL estimates the traditional translation, rotation
% and scaling parameters from a homogenous transformation matrix. It will
% give an error if the homogenous matrix also describes a perspective
% transformation.
%
% Use as
%   f = homogenous2traditional(H)
% where H is a 4x4 homogenous transformation matrix and f is a vector with
% nine elements describing
%   x-shift
%	  y-shift
%	  z-shift
% followed by the
%	  pitch (rotation around x-axis in degrees)
%	  roll  (rotation around y-axis in degrees)
%	  yaw   (rotation around z-axis in degrees)
% followed by the
%   x-rescaling factor
%   y-rescaling factor
%	  z-rescaling factor
%
% The order in which the transformations would be done is exactly opposite
% as the list above, i.e. first z-rescale ... and finally x-shift.
%
% Example use:
%   t0 = [1 2 3]; r0 = [10 20 30]; s0 = [1.1 1.2 1.3]
%   H0 = translate(t0) * rotate(r0) * scale(s0)
%   f = homogenous2traditional(H0)
%   t1 = f(1:3); r1 = f(4:6); s1 = f(7:9);
%   H1 = translate(t1) * rotate(r1) * scale(s1)
%
% See also TRANSLATE, ROTATE, SCALE, HOMOGENOUS2QUATERNION, QUATERNION

% Copyright (C) 2005-2022, Robert Oostenveld
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

tol = 100*eps;

% The homogenous transformation matrix is built up according to
%   H = T * R * S
% where
%   R = Rx * Ry * Rz

% make a copy of the input homogenous transformation matrix
TRS = H;

% estimate the translation
tx = H(1,4);
ty = H(2,4);
tz = H(3,4);
T = [
  1 0 0 tx
  0 1 0 ty
  0 0 1 tz
  0 0 0 1
  ];
% recompute the homogenous transformation excluding the translation
RS = T\TRS;

% estimate the scaling
sx = norm(H(1:3,1));
sy = norm(H(1:3,2));
sz = norm(H(1:3,3));
S = [
  sx 0 0 0
  0 sy 0 0
  0 0 sz 0
  0 0 0  1
  ];
% recompute the homogenous matrix excluding the scaling
R = RS/S;

% the difficult part is to determine the rotations, since the order of the rotations matters
%   R = Rx * Ry * Rz

if norm(R-eye(4))<tol
  % there is no rotation
  rx = 0;
  ry = 0;
  rz = 0;

else
  % perform the rotation of a probe point along the z-axis
  probez = [0 0 1 1]';
  probez = R * probez;

  % the rotation around the y-axis results in an offset along the positive x-direction
  ry = asin(probez(1));

  % the rotation around the x-axis can subsequently be estimated by the projection on the yz-plane
  rx = -atan(probez(2)/probez(3));

  % recompute the individual rotation matrices
  Rx = rotate(rad2deg([rx 0 0]));
  Ry = rotate(rad2deg([0 ry 0]));
  Rz = inv(Ry) * inv(Rx) * R; % R = Rx * Ry * Rz, use left side multiplication

  % compute the remaining rotation using a probe point on the x-axis
  probex = [1 0 0 1]';
  probex = Rz * probex;
  rz = asin(probex(2));
end

% these should be in degrees
rx = rad2deg(rx);
ry = rad2deg(ry);
rz = rad2deg(rz);

% the complete rotation matrix can be reconstructed like this
R = rotate([rx 0 0]) * rotate([0 ry 0]) * rotate([0 0 rz]);

% compare the original translation with the one that was estimated
if norm(H - T * R * S)>tol
  error('cannot split the homogenous transformation into scale, rotate and translate')
end

f = [tx ty tz rx ry rz sx sy sz];


