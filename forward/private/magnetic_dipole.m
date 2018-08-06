function [lf] = magnetic_dipole(R, pos, ori)

% MAGNETIC_DIPOLE leadfield for a magnetic dipole in an infinite medium
%
% [lf] = magnetic_dipole(R, pos, ori)
%
% with input arguments
%   R           position dipole
%   pos         position magnetometers
%   ori         orientation magnetometers
%
% See also CURRENT_DIPOLE

% Copyright (C) 2003, Robert Oostenveld
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

mu0   = 4*pi*1e-7;
nchan = size(pos,1);

% ensure that the dipole position is a row vector
R = reshape(R, [1 3]);

% shift the magnetometer coils so that the dipole is in the origin
pos(:,1) = pos(:,1) - R(1);
pos(:,2) = pos(:,2) - R(2);
pos(:,3) = pos(:,3) - R(3);

% change the variable names for convenience
% R position of magnetometer, relative to dipole
% r distance of magnetometer, relative to dipole
R = pos;
r = sqrt(sum(R.^2,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the original code, which follows the physical formulation closely
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% magnetic field on all magnetometer positions for an x-oriented dipole
% Bmx = mu0/(4*pi) * (3 * repmat(R(:,1), [1 3]) .* R - repmat([1 0 0], [nchan 1]) .* repmat(r.^2, [1 3])) ./ repmat(r.^5, [1 3]);
% magnetic field on all magnetometer positions for an y-oriented dipole
% Bmy = mu0/(4*pi) * (3 * repmat(R(:,2), [1 3]) .* R - repmat([0 1 0], [nchan 1]) .* repmat(r.^2, [1 3])) ./ repmat(r.^5, [1 3]);
% magnetic field on all magnetometer positions for an z-oriented dipole
% Bmz = mu0/(4*pi) * (3 * repmat(R(:,3), [1 3]) .* R - repmat([0 0 1], [nchan 1]) .* repmat(r.^2, [1 3])) ./ repmat(r.^5, [1 3]);
% compute the field along the orientation of each magnetometer for an x/y/z oriented dipole
% lf(:,1) = dot(Bmx, ori, 2);
% lf(:,2) = dot(Bmy, ori, 2);
% lf(:,3) = dot(Bmz, ori, 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the optimized code, which make the previous computation approximately
% three times more efficent by introducing temporary variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r2 = repmat(r.^2, [1 3]);
r5 = repmat(r.^5, [1 3]);
x = R(:,1); x = [x x x];
y = R(:,2); y = [y y y];
z = R(:,3); z = [z z z];
mx = zeros(nchan,3); mx(:,1) = 1;
my = zeros(nchan,3); my(:,2) = 1;
mz = zeros(nchan,3); mz(:,3) = 1;
Tx = (3 * x .* R - mx .* r2);
Ty = (3 * y .* R - my .* r2);
Tz = (3 * z .* R - mz .* r2);
lf(:,1) = dot(Tx, ori, 2);
lf(:,2) = dot(Ty, ori, 2);
lf(:,3) = dot(Tz, ori, 2);
lf = mu0/(4*pi) * lf ./ r5;

