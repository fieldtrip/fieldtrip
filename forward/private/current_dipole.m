function [lf] = current_dipole(R, pos, ori)

% CURRENT_DIPOLE leadfield for a current dipole in an infinite homogenous medium
%
% [lf] = current_dipole(R, pos, ori)
%
% with input arguments
%   R           position dipole
%   pos         position magnetometers
%   ori         orientation magnetometers
%
% This implements equation 9.3-1 from R.M. Gulrajani (1998) Bioelectricity and
% Biomagnetism, John Wiley and Sons, ISBN 04712485252.
%
% See also MAGNETIC_DIPOLE

% Copyright (C) 2013, Robert Oostenveld
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

lf = zeros(nchan,3);

% for i=1:nchan
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % this is the original code, which follows the physical formulation closely
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Bx  = dot((mu0/(4*pi)) * cross([1 0 0], pos(i,:))/norm(pos(i,:))^3, ori(i,:));
%   By = dot((mu0/(4*pi)) * cross([0 1 0], pos(i,:))/norm(pos(i,:))^3, ori(i,:));
%   Bz = dot((mu0/(4*pi)) * cross([0 0 1], pos(i,:))/norm(pos(i,:))^3, ori(i,:));
%   lf(i,:) = [Bx By Bz];
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is the slightly optimized implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nchan
  lf(i,:) = cross(pos(i,:), ori(i,:))/norm(pos(i,:))^3;
end
lf = lf * mu0/(4*pi);
