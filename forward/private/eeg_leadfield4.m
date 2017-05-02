function [lf, vol] = eeg_leadfield4(R, elc, vol)

% EEG_LEADFIELD4 electric leadfield for a dipole in 4 concentric spheres
% 
% [lf] = eeg_leadfield4(R, elc, vol)
%
% with input arguments
%   R          position of the dipole
%   elc        position of the electrodes
% and vol being a structure with the elements
%   vol.r      radius of the 4 spheres 
%   vol.cond   conductivity of the 4 spheres
%   vol.t      constant factors for series expansion (optional)
%
% The center of the spheres should be at the origin.
%
% This implementation is adapted from
%   Lutkenhoner, Habilschrift 1992.
% The original reference is
%  Cuffin BN, Cohen D. Comparison of the magnetoencephalogram and electroencephalogram. Electroencephalogr Clin Neurophysiol. 1979 Aug;47(2):132-46. 
%
% See also EEG_LEADFIELD4_PREPARE for precomputing the constant factors,
% which can save time when multiple leadfield computations are done.

% Copyright (C) 2002, Robert Oostenveld
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

% sort the spheres from the smallest to the largest
[vol.r, indx] = sort(vol.r);
vol.cond      = vol.cond(indx);

% use more convenient names for the radii and conductivities
r1 = vol.r(1); c1 = vol.cond(1);
r2 = vol.r(2); c2 = vol.cond(2);
r3 = vol.r(3); c3 = vol.cond(3);
r4 = vol.r(4); c4 = vol.cond(4);

% check whether the electrode ly on the sphere, allowing 0.5% tolerance
dist = sqrt(sum(elc.^2,2));
if any(abs(dist-r4)>r4*0.005)
  warning('electrodes do not ly on sphere surface -> using projection')
end
elc = r4 * elc ./ [dist dist dist];

% check whether the dipole is inside the brain [disabled for EEGLAB]
% if sqrt(sum(R.^2))>=r1
%   error('dipole is outside the brain compartment');
% end

% rotate everything so that the dipole is along the pos. z-axis
% only if the dipole is not in the origin or along the positive z-axis
if R(1)~=0 || R(2)~=0
  % compute the rotation matrix
  % the inverse rotation matrix is the transposed of this one
  val1 = norm(R);
  val2 = norm(R(1:2));
  rot(1,1) = R(1) * R(3) / (val1 * val2); 
  rot(1,2) = R(2) * R(3) / (val1 * val2);
  rot(1,3) = -1.0 * val2 / val1;
  rot(2,1) = -1.0 * R(2) / val2;
  rot(2,2) =        R(1) / val2;
  rot(2,3) =                  0; 
  rot(3,:) = R ./ val1;
  % rotate the electrodes
  elc = elc*rot';
elseif R(1)==0 && R(2)==0 && R(3)<0
  % dipole on negative z-axis, this case is very simple: reflect on xy-plane
  elc(:,3) = -elc(:,3);
else
  % dipole is on positive z-axis, nothing has to be done
end

% compute the constant factors for the sphere configuration if needed
if ~isfield(vol, 't')
  vol.t = eeg_leadfield4_prepare(vol);
end

Nchans = size(elc,1);
lf     = zeros(Nchans,3);
Nmax   = length(vol.t);
n      = 1:Nmax;
f      = norm(R)/r4;        % following cuffin1979
% c      = r2/r4;       % following cuffin1979
% d      = r3/r4;       % following cuffin1979

% this code is to cross-validate the lutkenhoner and cuffin implementations
% [lut_t, cuf_t] = eeg_leadfield4_prepare(vol);
% lut_c = (2*n+1).^4.*f.^(n-1) ./ (lut_t.*4*pi*c4*r4^2);
% cuf_c = (2*n+1).^4.*f.^(n-1) .*(c*d).^(2.*n+1) ./ (cuf_t.*4*pi*c4*r4^2);

% given a fixed volume conductor, these only need to be computed once for all electrodes
const = (2*n+1).^4.*f.^(n-1) ./ (vol.t.*4*pi*c4*r4^2);

for i=1:Nchans
  % convert the position of the electrodes to spherical coordinates
  [phi, el] = cart2sph(elc(i,1), elc(i,2), elc(i,3));

  % change from colatitude to latitude and compute the cosine 
  cos_theta = cos(pi/2-el);

  % the series summation starts at zero
  s_x = 0;
  s_z = 0;

  for n=1:Nmax
    P0  = plgndr(n,0,cos_theta);        % zero'th order Legendre
    P1  = plgndr(n,1,cos_theta);        % first order Legendre
    s_x = s_x + const(n)*P1/n;          % s_y is identical
    s_z = s_z + const(n)*P0;
  end

  lf(i,1) = -cos(phi) * s_x;
  lf(i,2) = -sin(phi) * s_x;            % s_y is identical to s_x
  lf(i,3) = 1         * s_z;
end

% apply the inverse rotation to the leadfield matrix
if R(1)~=0 || R(2)~=0
  lf = lf*rot;
elseif R(1)==0 && R(2)==0 && R(3)<0
  lf(:,3) = -lf(:,3);
end

