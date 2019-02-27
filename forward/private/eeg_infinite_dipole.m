function [lf] = eeg_infinite_dipole(dippos, elc, vol)

% EEG_INFINITE_DIPOLE calculate the infinite medium leadfield on electrode positions
% elc for a dipole at dippos and with the conductivity cond.
%
% Use as
%   [lf] = eeg_infinite_dipole(R, elc, vol)
%
% See also EEG_INFINITE_MONOPOLE, EEG_HALFSPACE_DIPOLE, EEG_HALFSPACE_MONOPOLE

% Copyright (C) 1998, Robert Oostenveld
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

if ~isstruct(vol)
  % it only represents the conductivity, make a structure out of it
  vol = struct('cond', vol);
end

siz = size(dippos);
if any(siz==1)
  % positions are specified as a single vector
  Ndipoles = prod(siz)/3;
  dippos = dippos(:)'; % ensure that it is a row vector
elseif siz(2)==3
  % positions are specified as a Nx3 matrix -> reformat to a single vector
  Ndipoles = siz(1);
  dippos = dippos';
  dippos = dippos(:)'; % ensure that it is a row vector
else
  ft_error('incorrect specification of dipole locations');
end

cond     = vol.cond;
Nelc     = size(elc,1);
lf       = zeros(Nelc,3*Ndipoles);
s1       = size(dippos);

if s1(1)>s1(2)
  % make sure that the dipole position is a row-vector
  dippos = dippos';
end

for i=1:Ndipoles
  r = elc - ones(Nelc,1) * dippos((1:3) + 3*(i-1));
  R = (4*pi*cond) * (sum(r' .^2 ) .^ 1.5)';
  if any(R==0)
    ft_warning('dipole collides with one of the electrodes');
  end
  lf(:,(1:3) + 3*(i-1)) = r ./ [R R R];
end
