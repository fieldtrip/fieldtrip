function [pnt, tri, nrm] = read_stl(filename);

% READ_STL reads a triangulation from an ascii *.stl file, which is a file
% format native to the stereolithography CAD software created by 3D Systems.
%
% Use as
%   [pnt, tri, nrm] = read_stl(filename)
%
% See also WRITE_STL

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: read_stl.m,v $
% Revision 1.1  2009/01/14 09:33:10  roboos
% moved even more files from fileio to fileio/private, see previous log entry
%
% Revision 1.3  2008/11/14 07:49:19  roboos
% use standard matlab strtrim function instead of deblank2
%
% Revision 1.2  2008/11/12 16:58:30  roboos
% open as text, fixed typo in solLid
%
% Revision 1.1  2006/04/19 12:54:47  roboos
% new implementation, to be used in construction of solid model of monkey brain
%

% solid   testsphere
%   facet normal -0.13 -0.13 -0.98
%     outer loop
%       vertex 1.50000 1.50000 0.00000
%       vertex 1.50000 1.11177 0.05111
%       vertex 1.11177 1.50000 0.05111
%     endloop
%   endfacet
%   ...

fid = fopen(filename, 'rt');

ntri = 0;
while ~feof(fid)
  line = fgetl(fid);
  ntri = ntri + ~isempty(findstr('facet normal', line));
end
fseek(fid, 0, 'bof');

tri = zeros(ntri,3);
pnt = zeros(ntri*3,3);

line = fgetl(fid);
name = sscanf(line, 'solid %s');

for i=1:ntri
  line1 = fgetl(fid);
  line2 = fgetl(fid);
  line3 = fgetl(fid);
  line4 = fgetl(fid);
  line5 = fgetl(fid);
  line6 = fgetl(fid);
  line7 = fgetl(fid);
  i1 = (i-1)*3+1;
  i2 = (i-1)*3+2;
  i3 = (i-1)*3+3;
  tri(i,:) = [i1 i2 i3];
  dum = sscanf(strtrim(line1), 'facet normal %f %f %f'); nrm(i,:) = dum(:)';
  dum = sscanf(strtrim(line3), 'vertex %f %f %f'); pnt(i1,:) = dum(:)';
  dum = sscanf(strtrim(line4), 'vertex %f %f %f'); pnt(i2,:) = dum(:)';
  dum = sscanf(strtrim(line5), 'vertex %f %f %f'); pnt(i3,:) = dum(:)';
end

fclose(fid);

