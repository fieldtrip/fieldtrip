function [pnt, dhk] = read_vtk(fn)

% READ_VTK reads a triangulation from a VTK (Visualisation ToolKit) format file
%
% [pnt, dhk] = read_vtk(filename)
%
% See also SAVE_VTK, READ_TRI, SAVE_TRI

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: read_vtk.m,v $
% Revision 1.1  2009/01/14 09:33:10  roboos
% moved even more files from fileio to fileio/private, see previous log entry
%
% Revision 1.2  2003/03/11 15:24:52  roberto
% updated help and copyrights
%

fid = fopen(fn, 'rt');
if fid~=-1

  npnt = 0;
  while (~npnt)
    line = fgetl(fid);
    if ~isempty(findstr(line, 'POINTS'))
       npnt = sscanf(line, 'POINTS %d float');
    end
  end
  pnt = zeros(npnt, 3);
  for i=1:npnt
    pnt(i,:) = fscanf(fid, '%f', 3)';
  end

  ndhk = 0;
  while (~ndhk)
    line = fgetl(fid);
    if ~isempty(findstr(line, 'POLYGONS'))
       tmp = sscanf(line, 'POLYGONS %d %d');
       ndhk = tmp(1);
    end
  end
  dhk = zeros(ndhk, 4);
  for i=1:ndhk
    dhk(i,:) = fscanf(fid, '%d', 4)';
  end
  dhk = dhk(:,2:4) + 1;

  fclose(fid);

else
  error('unable to open file');
end


