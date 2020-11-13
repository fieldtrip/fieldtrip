function [lf] = eeg_halfspace_monopole(monpos, elc, vol)

% EEG_HALFSPACE_MONOPOLE calculate the leadfield on positions elc for a monopole at
% position monpos. The halfspace solution requires a plane dividing a conductive zone
% (cond > 0), from a non-coductive zone (cond = 0).
%
% Use as
%   [lf] = eeg_halfspace_monopole(monpos, elc, vol)
%
% See also EEG_INFINITE_DIPOLE, EEG_INFINITE_MONOPOLE, EEG_HALFSPACE_DIPOLE

% Copyright (C) 2011, Cristiano Micheli
% Copyright (C) 2019, Robert Oostenveld
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

if isfield(vol, 'pos')
  % this is for forward/backward compatibility
  vol.pnt = vol.pos;
  vol = rmfield(vol, 'pos');
end

siz = size(monpos);
if any(siz==1)
  % positions are specified as a single vector
  Npoles = prod(siz)/3;
  monpos = monpos(:)'; % ensure that it is a row vector
elseif siz(2)==3
  % positions are specified as a Nx3 matrix -> reformat to a single vector
  Npoles = siz(1);
  monpos = monpos';
  monpos = monpos(:)'; % ensure that it is a row vector
else
  ft_error('incorrect specification of monopole locations');
end

Nelc     = size(elc,1);
lf       = zeros(Nelc,Npoles);

for i=1:Npoles
  % this is the position of monopole "i"
  monopole1 = monpos((1:3) + 3*(i-1));
  
  % find the position of a mirror monopole symmetric to the plane
  monopole2 = get_mirror_pos(monopole1, vol);
  
  % compute the potential of the original and the mirror dipole
  lf1 = eeg_infinite_monopole(monopole1, elc, vol);
  lf2 = eeg_infinite_monopole(monopole2, elc, vol);
  
  % take the sum of the two monopoles
  lf(:,i) = lf1 + lf2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P2 = get_mirror_pos(P1,vol)
% calculates the position of a point symmetric to pnt with respect to a plane

% define the plane
pnt = vol.pnt;
ori = vol.ori; % already normalized

if abs(dot(P1-pnt,ori))<eps
  ft_warning(sprintf ('point %f %f %f lies in the symmetry plane',P1(1),P1(2),P1(3)))
  P2 = P1;
else
  % define the plane in parametric form
  % define a non colinear vector vc with respect to the plane normal
  vc = [1 0 0];
  if abs(cross(ori, vc, 2))<eps
    vc = [0 1 0];
  end
  % define plane's direction vectors
  v1 = cross(ori, vc, 2);  v1 = v1/norm(v1);
  v2 = cross(pnt, ori, 2); v2 = v2/norm(v2);
  plane = [pnt v1 v2];
  
  % distance plane-point P1
  d = abs(dot(ori, plane(:,1:3)-P1(:,1:3), 2));
  
  % symmetric point
  P2 = P1 + 2*d*ori;
end
