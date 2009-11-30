function [pnt, dhk] = read_off(fn);

% READ_OFF reads vertices and triangles from a OFF format triangulation file
%
% [pnt, dhk] = read_off(filename)
%
% See also READ_TRI, READ_BND

% Copyright (C) 1998, Robert Oostenveld
%
% Subversion does not use the Log keyword, use 'svn log <filename>' or 'svn -v log | less' to get detailled information

fid = fopen(fn, 'rt');
if fid~=-1

  % scan the file type
  [s, count] = fscanf(fid, '%s\n', 1);
  if ~strcmp(s,'OFF')
    msg = sprintf('wrong file type %s', s);
    error(msg);
  end

  % read the number of vertex points and triangles
  [val, count] = fscanf(fid, '%d', 3);
  Npnt = val(1)
  Ndhk = val(2)

  % read the vertex points
  pnt  = fscanf(fid, '%f', [3, Npnt]);
  pnt  = pnt(1:3,:)';

  % read the triangles
  dhk = fscanf(fid, '%d', [4, Ndhk]);
  dhk = (dhk(2:4,:)+1)';
  fclose(fid);

else
  error('unable to open file');
end


