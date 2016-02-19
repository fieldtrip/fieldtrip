function lf = eeg_leadfield1(R, elc, vol)

% EEG_LEADFIELD1 electric leadfield for a dipole in a single sphere
%
% [lf] = eeg_leadfield1(R, elc, vol)
%
% with input arguments
%   R         position dipole (vector of length 3)
%   elc       position electrodes
% and vol being a structure with the elements
%   vol.r     radius of sphere
%   vol.cond  conductivity of sphere
%
% The center of the sphere should be at the origin.
%
% This implementation is adapted from
%   Luetkenhoener, Habilschrift '92
% The original reference is
%   R. Kavanagh, T. M. Darccey, D. Lehmann, and D. H. Fender. Evaluation of methods for three-dimensional localization of electric sources in the human brain. IEEE Trans Biomed Eng, 25:421-429, 1978.

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

Nchans = size(elc, 1);
lf = zeros(Nchans,3);

% always take the outermost sphere, this makes comparison with the 4-sphere computation easier
[vol.r, indx] = max(vol.r);
vol.cond = vol.cond(indx);

% check whether the electrode ly on the sphere, allowing 0.5% tolerance
dist = sqrt(sum(elc.^2,2));
if any(abs(dist-vol.r)>vol.r*0.005)
  warning('electrodes do not ly on sphere surface -> using projection')
end
elc = vol.r * elc ./ [dist dist dist];

% check whether the dipole is inside the brain [disabled for EEGLAB]
% if sqrt(sum(R.^2))>=vol.r
%   error('dipole is outside the brain compartment');
% end

c0 = norm(R);
c1 = vol.r;
c2 = 4*pi*c0^2*vol.cond;

if c0==0
  % the dipole is in the origin, this can and should be handeled as an exception
  [phi, el] = cart2sph(elc(:,1), elc(:,2), elc(:,3));
  theta = pi/2 - el;
  lf(:,1) = sin(theta).*cos(phi);
  lf(:,2) = sin(theta).*sin(phi);
  lf(:,3) = cos(theta);
  % the potential in a homogenous sphere is three times the infinite medium potential
  lf = 3/(c1^2*4*pi*vol.cond)*lf;

else
  for i=1:Nchans
    % use another name for the electrode, in accordance with lutkenhoner1992
    r = elc(i,:);

    c3 = r-R;
    c4 = norm(c3);
    c5 = c1^2 * c0^2 - dot(r,R)^2;   % lutkenhoner A.11
    c6 = c0^2*r - dot(r,R)*R;        % lutkenhoner, just after A.17

    % the original code reads (cf. lutkenhoner1992 equation A.17)
    % lf(i,:) = ((dot(R, r/norm(r) - (r-R)/norm(r-R))/(norm(cross(r,R))^2) + 2/(norm(r-R)^3)) * cross(R, cross(r, R)) + ((norm(r)^2-norm(R)^2)/(norm(r-R)^3) - 1/norm(r)) * R) / (4*pi*vol.cond(1)*norm(R)^2);

    % but more efficient execution of the code is achieved by some precomputations
    if c5<1000*eps
      % the dipole lies on a single line with the electrode
      lf(i,:) = (2/c4^3 * c6 + ((c1^2-c0^2)/c4^3 - 1/c1) * R) / c2;
    else
      % nothing wrong, do the complete computation
      lf(i,:) = ((dot(R, r/c1 - c3/c4)/c5 + 2/c4^3) * c6 + ((c1^2-c0^2)/c4^3 - 1/c1) * R) / c2;
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast cross product
function [c] = cross(a,b)
c = [a(2)*b(3)-a(3)*b(2) a(3)*b(1)-a(1)*b(3) a(1)*b(2)-a(2)*b(1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fast dot product
function [c] = dot(a,b)
c = sum(a.*b);

