function vol = read_asa_vol(fn);

% READ_ASA_VOL reads an ASA volume conductor file
%
% all data is converted to the following units
%   vertices        mm
%   conductivities  S/m

% Copyright (C) 2002, Robert Oostenveld
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

Nbnd  = read_asa(fn, 'NumberBoundaries=', '%d');
UnitC = read_asa(fn, 'UnitConduct', '%s');
UnitP = read_asa(fn, 'UnitPosition', '%s');
cond  = read_asa(fn, 'Conductivities', '%f');
radii = read_asa(fn, 'Radii', '%f');
pos   = read_asa(fn, 'Positions', '%f');
bnd1  = read_asa(fn, 'Boundary1', '%s');
bnd2  = read_asa(fn, 'Boundary2', '%s');
bnd3  = read_asa(fn, 'Boundary3', '%s');
bnd4  = read_asa(fn, 'Boundary4', '%s');

if ~isempty(radii) | ~isempty(pos)
  % this appears to be a spherical volume conductor
  if strcmpi(UnitP,'mm')
    radii = 1*radii;
    pos   = 1*pos;
  elseif strcmpi(UnitP,'cm')
    radii = 100*radii;
    pos   = 100*pos;
  elseif strcmpi(UnitP,'m')
    radii = 1000*radii;
    pos   = 1000*pos;
  else
    error(sprintf('Unknown unit of distance for volume (%s)', UnitP));
  end
end

if strcmpi(UnitC,'s/m')
  cond = cond/1;
elseif strcmpi(UnitC,'s/cm')
  cond = cond/100;
elseif strcmpi(UnitC,'s/mm')
  cond = cond/1000;
else
  error(sprintf('Unknown unit of conductivity for volume (%s)', UnitC));
end

if ~isempty(radii)
  % this appears to be a spherical volume conductor
  vol.radius = radii;
  vol.cond   = cond;
  vol.center = pos;
else
  % this appears to be a realistical volume conductor
  [path, name, ext] = fileparts(fn);
  if Nbnd>=1
    vol.bnd(1) = read_asa_bnd(fullfile(path, bnd1));
  end
  if Nbnd>=2
    vol.bnd(2) = read_asa_bnd(fullfile(path, bnd2));
  end
  if Nbnd>=3
    vol.bnd(3) = read_asa_bnd(fullfile(path, bnd3));
  end
  if Nbnd>=4
    vol.bnd(4) = read_asa_bnd(fullfile(path, bnd4));
  end
  if Nbnd>=5
    error('cannot read more than 4 boundaries');
  end

  % if there is a precomputed matrix, read it from an external file
  mat_file = read_asa(fn, 'Matrix', '%s');
  if ~isempty(mat_file)
    nr = read_asa(fullfile(path, mat_file), 'NumberRows=', '%d');
    nc = read_asa(fullfile(path, mat_file), 'NumberColumns=', '%d');
    mab_file = read_asa(fullfile(path, mat_file), 'Matrix', '%s');
    fid = fopen(fullfile(path, mab_file), 'rb', 'ieee-le');
    if fid==-1
      error(sprintf('could not open file %s', mab_file));
    else
      vol.mat = fread(fid, [nr nc], 'float32');
      fclose(fid);
      % remove the factor 2 that ASA assumes in the system matrix
      vol.mat = vol.mat/2;
      % scale the system matrix corresponding to vertex coordinates in mm
      vol.mat = vol.mat*100;
    end
  end

  vol.cond   = cond;
end

% remove all empty fields
field=fieldnames(vol);
for i=1:length(field)
  if isempty(getfield(vol, field{i}))
    vol = rmfield(vol, field{i});
  end
end

