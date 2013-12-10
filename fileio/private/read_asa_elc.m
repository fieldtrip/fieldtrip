function elec = read_asa_elc(fn)

% READ_ASA_ELC reads electrodes from an ASA electrode file
% converting the units to mm

% Copyright (C) 2002-2013, Robert Oostenveld
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

% the older *.elc files have an Nx3 matrix with positions and separate labels
% the newer *.elc files are formatted like this
%    Fp1:    94.9    30.7    14.0
% and also include Positions2D


Npnt = read_asa(fn, 'NumberPositions=', '%d');
Ntri = read_asa(fn, 'NumberPolygons=', '%d');
Unit = read_asa(fn, 'UnitPosition', '%s', 1);
pnt  = read_asa(fn, 'Positions', '%f', Npnt, ':');
prj  = read_asa(fn, 'Positions2D', '%f', Npnt, ':'); % only in newer files
tri  = read_asa(fn, 'Polygons', '%d', Ntri);
lab  = read_asa(fn, 'Labels', '%s', Npnt);
ref  = read_asa(fn, 'ReferenceChannel', '%s', 1); % only in newer files

if strcmpi(Unit,'mm')
  pnt = 1*pnt;
elseif strcmpi(Unit,'cm')
  pnt = 100*pnt;
elseif strcmpi(Unit,'m')
  pnt = 1000*pnt;
elseif ~isempty(Unit)
  error('Unknown unit of distance for electrodes (%s)', Unit);
end

tmp = tokenize(lab{1});
if length(tmp)==size(pnt,1)
  % the electrode labels were on a single line
  % reformat the electrode labels into an appropriately sized cell array
  lab = tmp;
end

elec.elecpos = pnt;
elec.label   = lab(:);
elec.unit    = Unit;
