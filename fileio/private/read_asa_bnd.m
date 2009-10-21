/function bnd = read_asa_bnd(fn);

% READ_ASA_BND reads an ASA boundary triangulation file
% converting the units of the vertices to mm

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: read_asa_bnd.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.5  2008/11/14 07:36:24  roboos
% use strcmpi instead of strcmp(lower())
%
% Revision 1.4  2004/03/29 15:16:08  roberto
% added semicolum (;) to end of line
%
% Revision 1.3  2003/06/05 08:34:23  roberto
% fixed bug in reading of vertices from binary file
%
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

Npnt = read_asa(fn, 'NumberPositions=', '%d');
Ndhk = read_asa(fn, 'NumberPolygons=', '%d');
Unit = read_asa(fn, 'UnitPosition', '%s');

pnt = read_asa(fn, 'Positions', '%f');
if any(size(pnt)~=[Npnt,3])
  pnt_file = read_asa(fn, 'Positions', '%s');
  [path, name, ext] = fileparts(fn);
  fid = fopen(fullfile(path, pnt_file), 'rb', 'ieee-le');
  pnt = fread(fid, [3,Npnt], 'float')';
  fclose(fid);
end

dhk = read_asa(fn, 'Polygons', '%f');
if any(size(dhk)~=[Ndhk,3])
  dhk_file = read_asa(fn, 'Polygons', '%s');
  [path, name, ext] = fileparts(fn);
  fid = fopen(fullfile(path, dhk_file), 'rb', 'ieee-le');
  dhk = fread(fid, [3,Ndhk], 'int32')';
  fclose(fid);
end

if strcmpi(Unit,'mm')
  pnt   = 1*pnt;
elseif strcmpi(Unit,'cm')
  pnt   = 100*pnt;
elseif strcmpi(Unit,'m')
  pnt   = 1000*pnt;
else
  error(sprintf('Unknown unit of distance for triangulated boundary (%s)', Unit));
end

bnd.pnt = pnt;
bnd.tri = dhk + 1;		% 1-offset instead of 0-offset
% bnd.tri = fliplr(bnd.tri);	% invert order of triangle vertices

