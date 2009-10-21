function vol = read_asa_vol(fn);

% READ_ASA_VOL reads an ASA volume conductor file
%
% all data is converted to the following units
%   vertices		mm
%   conductivities	S/m

% Copyright (C) 2002, Robert Oostenveld
% 
% $Log: read_asa_vol.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.4  2008/11/14 07:36:24  roboos
% use strcmpi instead of strcmp(lower())
%
% Revision 1.3  2003/12/16 10:24:31  roberto
% added ieee-le to binary reading of mat file to fix problem on Mac
%
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

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

