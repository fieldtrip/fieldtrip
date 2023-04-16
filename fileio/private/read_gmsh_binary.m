function [nodes, elements, nodedata, elementdata] = read_gmsh_binary(filename)

% READ_GMSH_BINARY reads a gmsh .msh binary file. Current support is only
% for version 2. There are some ASCII-readers floating around on the net,
% but they do not seem to work with the primary use case in FieldTrip (and
% the test data that I have available), which is SimNibs generated data.
%
% See also MESH_LOAD_GMSH4

% Copyright (C) 2023, Jan-Mathijs Schoffelen
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

fid = fopen(filename, 'rb');
txt = fgetl(fid);

if startsWith(txt, '$MeshFormat')
  % the file starts, as expected with the file header
  txt    = fgetl(fid);
  format = sscanf(txt, '%f %d %d');
  version  = format(1);
  isbinary = format(2)==1;
  dtype    = format(3);
  
  % read until the end of this chunk
  while 1
    % in a binary file, there's an endianness check byte, which is
    % currently not handled here.
    txt = fgetl(fid);
    if startsWith(txt, '$EndMeshFormat')
      break;
    end
  end
else
  ft_error('File %s cannot be read, header information is missing', filename);
end

if floor(version)~=2 || ~isbinary
  if ~isbinary
    str = 'not';
  else
    str = '';
  end
  ft_error('Currently only binary version 2 files can be read, the version is: %1.1f, and the file is %s binary', version, str);
end
if dtype~=8
  ft_error('Currently only 8-byte precision data is supported');
end

nodes       = [];
elements    = [];
nodedata    = [];
elementdata = [];

% read the contents of the file
while 1
  txt = fgetl(fid);
  if txt==-1
    % end of file
    fclose(fid);
    break;
  end
  % Nodes and Elements are expected to be minimally present, needed to
  % describe the mesh topology
  if startsWith(txt, '$Nodes')
    nodes = getnodes(fid);
  elseif startsWith(txt, '$Elements')
    elements = getelements(fid);
  elseif startsWith(txt, '$PhysicalNames')
    keyboard
  elseif startsWith(txt, '$Entities')
    keyboard
  elseif startsWith(txt, '$PartitionedEntities')
    keyboard
  elseif startsWith(txt, '$Periodic')
    keyboard
  elseif startsWith(txt, '$GhostElements')
    keyboard
  elseif startsWith(txt, '$Parametrizations')
    keyboard
  elseif startsWith(txt, '$NodeData')
    keyboard
  elseif startsWith(txt, '$ElementData')
    elementdata = getelementdata(fid);
  elseif startsWith(txt, '$ElementNodeData')
    keyboard
  elseif startsWith(txt, '$InterpolationScheme')
    keyboard
  else
    keyboard

  end
end

%%%%%%%%%%%%%%
% subfunctions

% $Nodes
function [nodes] = getnodes(fid)

% the chunk starts with an ascii-line with the number of elements
%
% the format of the nodes is 4-byte indx followed by 3x a float for the
% coordinates

txt = fgetl(fid);
N   = sscanf(txt, '%f');
ptr = ftell(fid);

fprintf('Reading position information for %d nodes\n', N);
nodes.indx = fread(fid, N, 'uint32', 24);
fseek(fid, ptr+4, 'bof');
nodes.nodes = fread(fid, 3*N, '3*double', 4);
nodes.nodes = reshape(nodes.nodes, 3, [])';
fseek(fid, ptr+28*N, 'bof');
txt = fgetl(fid);
assert(isequal(txt, '$EndNodes'), 'Reading of the nodes unexpectedly failed');

% $Elements
function [elements] = getelements(fid)

% the chunk starts with an ascii-line mentioning the overall number of elements, then 3 integers
% per sub-chunk of elements, that indicate the type, the number of elements
% in the sub-chunk, and the number of tags

txt = fgetl(fid);
N   = sscanf(txt, '%f');

fprintf('Reading element information for %d elements\n', N);
Ntotal = 0;
while N>Ntotal

  hdr  = fread(fid, 3, 'uint32');
  type = hdr(1);
  Nsub = hdr(2);
  ntag = hdr(3);

  [nnode, type_str] = type2nnode(type);
  
  fprintf('Reading %d %s with %d tags\n', Nsub, type_str, ntag);
  tmp = fread(fid, Nsub.*(nnode+ntag+1), '*uint32');
  tmp = reshape(tmp, nnode+ntag+1, [])';

  elements.(type_str)                     = tmp(:, (ntag+2):end);
  elements.(sprintf('%s_indx', type_str)) = tmp(:, 1);
  elements.(sprintf('%s_tag',  type_str)) = tmp(:, 1 + (1:ntag));
  Ntotal = Ntotal + Nsub;
end
txt = fgetl(fid);
assert(isequal(txt, '$EndElements'), 'Reading of the elements unexpectedly failed');

% $ElementData
function [elementdata] = getelementdata(fid)

% stringTags
txt = fgetl(fid);
N   = sscanf(txt, '%d');
for k = 1:N
  stringTag{k,1} = fgetl(fid);
end

% realTags
txt = fgetl(fid);
N   = sscanf(txt, '%d');
for k = 1:N
  txt = fgetl(fid);
  realTag(k,1) = sscanf(txt, '%d');
end

% integerTags
txt = fgetl(fid);
N   = sscanf(txt, '%d'); % N is kind of expected to be 4
for k = 1:N
  txt = fgetl(fid);
  integerTag(k,1) = sscanf(txt, '%d');
end
nfc = integerTag(2); % number of field components
N   = integerTag(3); % number of elements

fprintf('Reading elementdata for %d elements\n', N);
ptr = ftell(fid);
elementdata.indx = fread(fid, N, 'uint32', nfc*8);
fseek(fid, ptr+4, 'bof');
elementdata.value = fread(fid, nfc*N, sprintf('%d*double', nfc), 4);
elementdata.value = reshape(elementdata.value, nfc, [])';
fseek(fid, ptr+(4+nfc*8)*N, 'bof');

txt = fgetl(fid);
assert(isequal(txt, '$EndElementData'), 'Reading of the elementdata unexpectedly failed');

% helper function to get the number of nodes per element + a human interpretable name
function [nnode, type_str] = type2nnode(type)

switch type
  case 1 
    % 2-node line
    nnode = 2; type_str = 'lines';
  case 2
    % 3-node triangle
    nnode = 3; type_str = 'triangles';
  case 3
    % 4-node quadrangle
    nnode = 4; type_str = 'quadrangles';
  case 4
    % 4-node tetrahedron
    nnode = 4; type_str = 'tetrahedra';
  case 5
    % 8-node hexahedron
    nnode = 8; type_str = 'hexahedra';
  case 6
    % 6-node prism
    nnode = 6; type_str = 'prisms';
  case 7
    % 5-node pyramid
    nnode = 5; type_str = 'pyramids';
  case 8 
    % 3-node seconde order line (2 nodes associated with the vertices and 1 with the edge)
    nnode = 3; type_str = 'line';
  case 9
    % 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)
    nnode = 6; type_str = 'triangles';
  case 10
    % 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face)
    nnode = 9; type_str = 'quadrangles';
  case 11
    % 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges)
    nnode = 10; type_str = 'tetrahedra';
  case 12
    % 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges,  6 with the faces and 1 with the volume)
    nnode = 27; type_str = 'hexahedra';
  case 13
    % 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces)
    nnode = 18; type_str = 'prisms';
  case 14
    % 14-node seconde order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face)
    nnode = 14; type_str = 'pyramids';
  case 15
    % 1-node point
    nnode = 1; type_str = 'points';
  otherwise
    %% FIXME add types 16-31, and 92 93
    keyboard
end

