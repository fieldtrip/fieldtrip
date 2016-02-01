function dip = fixdipole(dip)

% FIXDIPOLE ensures that the dipole position and moment are
% consistently represented throughout FieldTrip functions.

% Copyright (C) 2009, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
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

[m, n] = size(dip.pos);

if n==3
  % the input representation is Nx3, which is what we want
elseif m==3
  % it is possible to translate it into a Nx3 unambiguously
  warning('input dipole positions should be specified as Nx3 matrix');
  dip.pos = dip.pos';
elseif m==1
  % it is possible to translate it into a Nx3 unambiguously
  warning('input dipole positions should be specified as Nx3 matrix');
  dip.pos = reshape(dip.pos, 3, n/3)';
else
  % it is not clear how to convert to a Nx3 matrix
  error('input dipole positions should be specified as Nx3 matrix');
end

if isfield(dip, 'mom')
  ndip = size(dip.pos,1);
  if numel(dip.mom)==ndip*3
    ntime = 1;
  else
    ntime = numel(dip.mom)/(ndip*3);
  end
  dip.mom = reshape(dip.mom, ndip*3, ntime);
else
  dip.mom = [];
end

