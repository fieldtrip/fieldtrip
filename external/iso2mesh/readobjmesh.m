function [node, face] = readobjmesh(fname)
%
% [node,face]=readobjmesh(fname)
%
% read Wavefront obj-formatted surface mesh files (.obj)
%
% author: Qianqian Fang, <q.fang at neu.edu>
%
% input:
%    fname: name of the .obj data file
%
% output:
%    node: node coordinates of the mesh
%    face: list of elements of the surface mesh
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)
%

str = fileread(fname);
nodestr = regexp(str, 'v\s+([0-9.\-e]+\s+[0-9.\-e]+\s+[0-9.\-e]+)', 'tokens');
nodestr = [nodestr{:}];
node = sscanf(sprintf('%s ', nodestr{:}), '%f %f %f', [3, inf])';
facestr = regexp(str, 'f\s+(\d+(/\d+)*\s+\d+(/\d+)*\s+\d+(/\d+)*)', 'tokens');
facestr = [facestr{:}];
facestr = sprintf('%s ', facestr{:});
facestr = regexprep(facestr, '/\d+', '');
face = sscanf(facestr, '%d %d %d', [3, inf])';
